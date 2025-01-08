#Preprcesamiento de datos y generacion del archivo .rds
#Monica Mitjans
#TFM-Shiny App

#Filtrado de los datos de expresion
library(edgeR)

#filtrar genes con más de 10 CPM en al menos el 30% de las muestras
keep_genes <- rowSums(cpm(counts) > 10) >= (ncol(counts) * 0.3)

#aplicar el filtro a los datos de expresion
filtered_counts <- counts[keep_genes, ]

#Normalizacion de elos datos de expresion
#crear un objeto DGEList y normalizar con TMM
dge<- DGEList(counts = filtered_counts)
dge<- calcNormFactors(dge, method = "TMM") #usar metodo TMM

#obtener los valores normalizados como CPM
normalized_counts<- cpm(dge, normalized.lib.sizes = TRUE)

#obtener los genes que fueron quedaron retenidos luego del filtrado
filtered_genes<- rownames(filtered_counts)  #genes que se mantuvieron

#filtrar el objeto SummarizedExperiment para que coincida con los genes filtrados
tcga_data_filtered<- tcga_data[filtered_genes, ]

# Verificar que el numero de filas entre el summarizedExperiment filtrado y filtered_counts sea igual
dim(assay(tcga_data_filtered, "unstranded"))
dim(filtered_counts)  #es igual en ambos: 11185  1231

#asignar los datos normalizados
rownames(normalized_counts)<- rownames(assay(tcga_data_filtered, "unstranded"))
colnames(normalized_counts)<- colnames(assay(tcga_data_filtered, "unstranded"))

#actualizar el assay "normalized" en el objeto SummarizedExperiment filtrado
assay(tcga_data_filtered, "normalized", withDimnames = FALSE) <- normalized_counts

#verificar los assays en el objeto actualizado
assayNames(tcga_data_filtered)  #aparece "normalized"


#Filtrar los datos clinicos por los subtipos que interesan (ER status, HER2 status, PR status - para TNBC)

clinical_filtered <- clinical_tab_all$clinical_patient_brca[
  (clinical_tab_all$clinical_patient_brca$her2_status_by_ihc %in% c("Positive", "Negative")) &
    (clinical_tab_all$clinical_patient_brca$er_status_by_ihc %in% c("Positive", "Negative")) &
    (clinical_tab_all$clinical_patient_brca$pr_status_by_ihc %in% c("Positive", "Negative")), 
]

#para filtrar el summarizedExperiment por paciente utilizar "patient"
#para filtrar los datos clinicos por paciente utlizar "bcr_patient_barcode"

head(tcga_data_filtered$patient)
head(clinical_filtered$bcr_patient_barcode)


#filtrar clinical_filtered para incluir solo pacientes que están en tcga_data_filtered.
library(SummarizedExperiment)
common_samples<-intersect(tcga_data_filtered$patient, clinical_filtered$bcr_patient_barcode)
common_samples

#filtrar datos de expresión usando solo los IDs comunes
expression_filtered<-tcga_data_filtered[, tcga_data_filtered$patient %in% common_samples]

#filtrar datos clínicos para incluir solo los IDs comunes y en el mismo orden que en los datos de expresión
clinical_filtered_common <- clinical_filtered[clinical_filtered$bcr_patient_barcode %in% common_samples, ]

#verificar el numero de pacientes en cada conjunto
length(unique(expression_filtered$patient))
length(expression_filtered$patient)  #hay pacientes duplicados
length(unique(clinical_filtered_common$bcr_patient_barcode))
length(clinical_filtered_common$bcr_patient_barcode) #no hay pacientes duplicados

#eliminar los duplicados
duplicated_samples<-duplicated(colData(expression_filtered)$patient)
#eliminar muestras duplicadas, manteniendo la primera ocurrencia
expression_filtered <- expression_filtered[, !duplicated_samples]  #quedan 720
#verificar que las duplicaciones se hayan eliminado
duplicated_patients<-expression_filtered$patient[duplicated(expression_filtered$patient)]
table(duplicated_patients)

#Modificar el clinical_filtered_common para incluir los pacientes que tenemos en los datos de expresion
#Obtener los IDs únicos de pacientes en expression_filtered
expression_patient_ids<-unique(colData(expression_filtered)$patient)
#filtrar el conjunto de datos clínicos para mantener solo los pacientes en expression_filtered
clinical_filtered_common<-clinical_filtered_common[clinical_filtered_common$bcr_patient_barcode %in% expression_patient_ids, ]
#verificar el número de pacientes después del filtrado
nrow(clinical_filtered_common) #hay 720

#verificar que hay entradas
clinical_filtered_common$er_status_by_ihc 

#verificar si todos los IDs de pacientes en expression_filtered están también en clinical_filtered_common
all(expression_filtered$patient %in% clinical_filtered_common$bcr_patient_barcode)

#verificar si todos los IDs de pacientes en clinical_filtered_common están también en expression_filtered
all(clinical_filtered_common$bcr_patient_barcode %in% expression_filtered$patient)

#ordenar expression_filtered por el identificador de paciente
expression_filtered<-expression_filtered[, order(colData(expression_filtered)$patient)]

#ordena clinical_filtered_common por el identificador de paciente
clinical_filtered_common <- clinical_filtered_common[order(clinical_filtered_common$bcr_patient_barcode), ]

#verificar
head(expression_filtered$patient)
head(clinical_filtered_common$bcr_patient_barcode)

#cerificar que ambos conjuntos estén ahora en el mismo orden de pacientes
if (identical(colData(expression_filtered)$patient, clinical_filtered_common$bcr_patient_barcode)) {
  # Asigna los datos clínicos al SummarizedExperiment si están alineados
  colData(expression_filtered)$HER2_status <- clinical_filtered_common$her2_status_by_ihc
  colData(expression_filtered)$ER_status <- clinical_filtered_common$er_status_by_ihc
  colData(expression_filtered)$PR_status <- clinical_filtered_common$pr_status_by_ihc
}

#verificar que los datos se han asignado correctamente
table(colData(expression_filtered)$HER2_status)
table(colData(expression_filtered)$ER_status)
table(colData(expression_filtered)$PR_status)

#filtrar pacientes TNBC (HER2 negativo, ER negativo, PR negativo)
tnbc_patient_ids<-colData(expression_filtered)$patient[
  colData(expression_filtered)$HER2_status == "Negative" &
    colData(expression_filtered)$ER_status == "Negative" &
    colData(expression_filtered)$PR_status == "Negative"
]
length(tnbc_patient_ids) #numero de pacientes tnbc

#A este punto, se tiene el SummarizedExperiment con los datos necesarios

#Curva de supervivencia inicial

OS_time <- ifelse(!is.na(as.numeric(clinical_filtered_common$death_days_to)), as.numeric(clinical_filtered_common$death_days_to), as.numeric(clinical_filtered_common$last_contact_days_to))
OS_status <- ifelse(clinical_filtered_common$vital_status == "Dead", 1, 0)


library(survival)
library(survminer)

# Curva de Overall Survival
surv_object_OS <- Surv(time = OS_time, event = OS_status)
fit_OS <- survfit(surv_object_OS ~ 1)

ggsurvplot(
  fit_OS,
  data = clinical_filtered_common,
  risk.table = FALSE,             # Desactivar la tabla de riesgo
  pval = FALSE,                   # Desactivar p-valor (sin grupos para comparar)
  conf.int = TRUE,                # Mostrar intervalo de confianza
  xlab = "Tiempo (días)",         # Etiqueta del eje X
  ylab = "Probabilidad de Supervivencia", # Etiqueta del eje Y
  title = "Curva de Supervivencia General (OS)",
  palette = "blue",               # Color de la curva
  legend.title = "Supervivencia",
  legend.labs = "Población General",
  ggtheme = theme_minimal(),      # Tema limpio
  surv.median.line = "hv"         # Mostrar línea de mediana de supervivencia
)



#################
#Debido a que el usuario porveera una lista con el gene symbol, se convierte los identificadores ensembl de los datos de expresion a gene_symbol usando biomart
# Convertir Ensembl IDs a nombres de genes usando biomaRt
library(biomaRt)
ensembl_ids<- gsub("\\..*", "", rownames(expression_filtered)) #obtener Ensembl IDs sin las versiones
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")#conectar con biomart
#ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "uswest") #usar mirror site en caso de que el sitio principal no funcione
conversion <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

#reemplazar rownames de expression_filtered con nombres de genes
converted_gene_names <- conversion$hgnc_symbol[match(ensembl_ids, conversion$ensembl_gene_id)]

#acceder al rowData del objeto expression_filtered 
row_data <- rowData(expression_filtered)
gene_names <- row_data$gene_name #acceder al gene_name

#crear un summarizedExperiment a partir de expression_filtered al que cambiaremos los rownames
expression_data<-expression_filtered #crear una copia de expression filtered
rownames(expression_data) <- gene_names #cambiar los identificadores ensembl por gene_names
#generar un vector que contiene los genes comunes entre converted_gene_names y gene_names - para disolver discrepancias
commongenes<-intersect(converted_gene_names,gene_names)

#filtrar expression_data para incluir solo los genes coincidentes
#obtener indices de los genes comunes en 'gene_names'
commongenes_indices <- which(gene_names %in% commongenes)

#filtrar el expression_data en base a los indices
expr_data <- assay(expression_data, "normalized") #obtener datos de expresión normalizada
rownames(expr_data)
filtered_expdata <- expr_data[commongenes_indices, ]

#obtener primeros 12 caracteres de los id de la matriz de expresion
clinical_filtered_common$bcr_patient_barcode #contiene 12 caracteres
expression_patient_ids<-substr(colnames(filtered_expdata), 1, 12)

#corroborar que los ids extraidos y los id de datos clinicos coincidan
#identificar los IDs comunes
common_patients_exp_clin <- intersect(expression_patient_ids, clinical_filtered_common$bcr_patient_barcode)
# Filtrar datos clínicos
filtered_clinical1 <- clinical_filtered_common[clinical_filtered_common$bcr_patient_barcode %in% common_patients_exp_clin, ]
nrow(filtered_clinical1) #720

# Filtrar datos de expresión

colnames(filtered_expdata)<-expression_patient_ids #12 caracteres

filtered_expression1 <- filtered_expdata[, expression_patient_ids %in% common_patients_exp_clin] #mantener los coincidentes
colnames(filtered_expression1)
ncol(filtered_expression1) #720

#asegurar de que ambos estén en el mismo orden
filtered_clinical1<- filtered_clinical1[match(common_patients_exp_clin, filtered_clinical1$bcr_patient_barcode), ]
filtered_expression1 <- filtered_expression1[, match(common_patients_exp_clin, expression_patient_ids)]
dim(filtered_expression1)
#verificar que son identicos
identical(colnames(filtered_expression1), filtered_clinical1$bcr_patient_barcode) #coinciden


#Generacion del dataset a ser utilizado por Shiny
library(SummarizedExperiment)

# Filtrar y crear un nuevo objeto SummarizedExperiment
#filtrar el rowdata con los indices de los genes
fil_row_data <- rowData(expression_data)[commongenes_indices, ]
#filtrar el coldata con los pacientes comunes
fil_col_data <- colData(expression_data)

#extraer los 12 primeros caracteres
rownames(fil_col_data)<-substr(rownames(fil_col_data), 1, 12)

# Verificar las coincidencias
identical(rownames(filtered_expdata),rownames(fil_row_data))
identical(colnames(filtered_expdata),rownames(fil_col_data))


#generar el summarizedExperiment
fil_exp_data_ <- SummarizedExperiment(
  assays = list(normalized = filtered_expdata),
  rowData = fil_row_data,
  colData = fil_col_data
)
dim(filtered_expdata)

#filtered_expression1 contiene datos de expresión
#filtered_clinical1 contiene datos clínicos emparejados con filtered_expression1
#integrar HER2 status, ER status, PR status, days to last contact, days to death y vital status

fil_exp_data_$vital_status<-filtered_clinical1$vital_status #asignar los vital status de los datos clinicos
fil_exp_data_$days_to_death<-filtered_clinical1$death_days_to #asignar los days_to death de los datos clinicos

#asignar last_contact_days_to de los datos clinicos a los pacientes con vital status alive
fil_exp_data_$days_to_death<-ifelse(
  fil_exp_data_$vital_status == "Alive", 
  filtered_clinical1$last_contact_days_to, 
  filtered_clinical1$death_days_to
)

#pasar days to death a numerico
fil_exp_data_$days_to_death<-as.numeric(fil_exp_data_$days_to_death)
# Convertir 'status' a binario (1 = fallecido, 0 = vivo)
fil_exp_data_$vital_status <- ifelse(fil_exp_data_$vital_status == "Dead", 1, 0)

#guardar el summarized.experiment
# Filtrar para incluir solo ciertas columnas en colData
summ_filtrado<- fil_exp_data_
colData(summ_filtrado) <- colData(fil_exp_data_)[, c("HER2_status", "ER_status", "PR_status", "vital_status", "days_to_death")]
colnames(colData(summ_filtrado)) <- c("HER2_status", "ER_status", "PR_status", "event", "time")
#guardar el summarizedExperiment
saveRDS(summ_filtrado, "datosdeexpresion.rds")






