#Validacion tecnica - Analisis de supervivencia
#Monica Mitjans
#TFM-Shiny App


#Analisis de suoervivencia en funcion a firmas genicas
#cargar librerias necesarias
library(survival)
library(survminer)
library(dplyr)

#Extraer datos clínicos y de expresión desde internal_data
datosclinicos<-as.data.frame(colData(internal_data))

datosdeexpresion<-t(scale(t(assay(internal_data, "normalized")))) #extraer datos de expresión  escalados 

#Comprobar que las columnas necesarias existen
required_cols<- c("HER2_status", "ER_status", "time", "event")
if (!all(required_cols %in% colnames(datosclinicos))) {
  stop("Las columnas necesarias no están presentes en los datos clínicos.")
}


#############################################
#--Analisis de supervivencia - p53 pathway--#
#############################################

#todos los pacientes
#Filtrar Genes de Interés (P53 SIGNALING PATHWAY)
library(msigdbr)

#descargar firmas génicas KEGG_P53_SIGNALING_PATHWAY
c2_gene_sets<- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg<- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest<- c2_kegg[c2_kegg$gs_name == "KEGG_P53_SIGNALING_PATHWAY", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol) #pasar a caracter

#filtrar genes de interés en los datos de expresión
selected_genes_indices <- which(rownames(datosdeexpresion) %in% genes_of_interest_vector)
selected_genes_expr <- datosdeexpresion[selected_genes_indices, ]

#calcular expresión promedio de la firma
expr_signature <- colMeans(selected_genes_expr, na.rm = TRUE)

#agrupar pacientes en base al nivel de expresion de la firma
group <- ifelse(
  expr_signature > median(expr_signature, na.rm = TRUE),
  "Alta", "Baja"
)

#agregar grupo al dataframe clínico
datosclinicos$group <- group

#verificar distribución de grupos
print(table(datosclinicos$group))

#asegurar que las longitudes coinciden
stopifnot(length(group) == nrow(datosclinicos))

#preparar el dataframe de supervivencia
dfsurv <- data.frame(
  time = as.numeric(datosclinicos$time),
  status = as.numeric(datosclinicos$event),
  group = datosclinicos$group
)

#eliminar filas con valores NA y duplicados
dfsurv <- na.omit(dfsurv)
dfsurv <- dfsurv[!duplicated(dfsurv[c("time", "status")]), ]

#comprobar valores NA restantes
print(paste("NA en time:", sum(is.na(dfsurv$time))))
print(paste("NA en status:", sum(is.na(dfsurv$status))))

# Cargar librerías necesarias
library(survival)
library(survminer)

#crear objeto de supervivencia
surv_obj <- Surv(time = dfsurv$time, event = dfsurv$status)

#ajustar modelo Kaplan-Meier
fit <- survfit(surv_obj ~ group, data = dfsurv)

#raficar la curva Kaplan-Meier

ggsurvplot(
  fit, data = dfsurv,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - Todos (P53 SIGNALING PATHWAY)",
  risk.table = FALSE
)

#filtrar muestras por subtipo de cancer

#HER2

her2_samples <- datosclinicos$HER2_status == "Positive" #selecciona Her2 positivos - cambiar de acuerdo al subtipo de interes
her2_expressiondata <- datosdeexpresion[, her2_samples]
her2_clinicaldata <- datosclinicos[her2_samples, ]

#Filtrar Genes de Interés (P53 SIGNALING PATHWAY)
library(msigdbr)

#descargar firmas génicas KEGG_P53_SIGNALING_PATHWAY
c2_gene_sets<- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg<- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest<- c2_kegg[c2_kegg$gs_name == "KEGG_P53_SIGNALING_PATHWAY", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol) #pasar a caracter

#filtrar genes de interés en los datos de expresión
selected_genes_indices <- which(rownames(her2_expressiondata) %in% genes_of_interest_vector)
selected_genes_expr <- her2_expressiondata[selected_genes_indices, ]

#verificar coherencia de las dimensiones
stopifnot(identical(colnames(her2_expressiondata), rownames(her2_clinicaldata)))

#calcular expresion promedio de la firma
expr_signature_her2 <- colMeans(selected_genes_expr, na.rm = TRUE)

#agrupar pacientes por nivel de expresión de la firma
group_her2 <- ifelse(
  expr_signature_her2 > median(expr_signature_her2, na.rm = TRUE),
  "Alta", "Baja"
)

#añadir el grupo al dataframe clínico
her2_clinicaldata$group <- group_her2

# Verificar distribución
table(her2_clinicaldata$group)

#asegurar que las longitudes coinciden
stopifnot(length(group_her2) == nrow(her2_clinicaldata))

#preparar el dataframe para el análisis de supervivencia
dfsurv_her2 <- data.frame(
  time = as.numeric(her2_clinicaldata$time),
  status = as.numeric(her2_clinicaldata$event),
  group = her2_clinicaldata$group
)
#eliminar NA y duplicados en columnas clave
dfsurv_her2 <- na.omit(dfsurv_her2)
dfsurv_her2 <- dfsurv_her2[!duplicated(dfsurv_her2[c("time", "status")]), ]

#verificar valores NA
sum(is.na(dfsurv_her2$time))
sum(is.na(dfsurv_her2$status))

#análisis de Supervivencia Kaplan-Meier
surv_obj_her2 <- Surv(time = dfsurv_her2$time, event = dfsurv_her2$status) #objeto de supervivencia
fit_her2 <- survfit(surv_obj_her2 ~ group, data = dfsurv_her2)

# Graficar la Curva Kaplan-Meier
ggsurvplot(
  fit_her2, data = dfsurv_her2,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - HER2+ (P53 SIGNALING PATHWAY)",
  risk.table = FALSE
)

#filtrar muestras: ER+
er_samples <- datosclinicos$ER_status == "Positive"
er_expressiondata <- datosdeexpresion[, er_samples]
er_clinicaldata <- datosclinicos[er_samples, ]

#filtrar genes de interes (P53 SIGNALING PATHWAY)
library(msigdbr)

#descargar firmas génicas KEGG_P53_SIGNALING_PATHWAY
c2_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg <- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest <- c2_kegg[c2_kegg$gs_name == "KEGG_P53_SIGNALING_PATHWAY", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol)

#filtrar genes de interés en los datos de expresión
selected_genes_indices<- which(rownames(er_expressiondata) %in% genes_of_interest_vector)
selected_genes_expr<- er_expressiondata[selected_genes_indices, ]

#verificar coherencia de las dimensiones
stopifnot(identical(colnames(er_expressiondata), rownames(er_clinicaldata)))

#calcular la firma de expresión promedio
expr_signature_er <- colMeans(selected_genes_expr, na.rm = TRUE)

#agrupar pacientes por expresión de firma génica
group_er <- ifelse(
  expr_signature_er > median(expr_signature_er, na.rm = TRUE),
  "Alta", "Baja"
)

#verificar distribución de grupos
table(group_er)

#asegurar que las longitudes coinciden
stopifnot(length(group_er) == nrow(er_clinicaldata))

#preparar el dataframe para el análisis de supervivencia
dfsurv_er<- data.frame(
  time = as.numeric(er_clinicaldata$time),
  status = as.numeric(er_clinicaldata$event),
  group = group_er
)

#verificar valores NA
sum(is.na(dfsurv_er$time))
sum(is.na(dfsurv_er$status))

#análisis de Supervivencia Kaplan-Meier
fit_er <- survfit(Surv(time = dfsurv_er$time, event = dfsurv_er$status) ~ group, data = dfsurv_er)

#graficar la Curva Kaplan-Meier
ggsurvplot(
  fit_er, data = dfsurv_er,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - ER+ (P53 SIGNALING PATHWAY)",
  risk.table = FALSE
)

#filtrar Muestras: TNBC (HER2-, ER-, PR-)
tnbc_samples <- datosclinicos$HER2_status == "Negative" & 
  datosclinicos$ER_status == "Negative" & 
  datosclinicos$PR_status == "Negative"
tnbc_expressiondata <- datosdeexpresion[, tnbc_samples]
tnbc_clinicaldata <- datosclinicos[tnbc_samples, ]

#descargar firmas génicas KEGG_P53_SIGNALING_PATHWAY
c2_gene_sets<- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg<- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest<- c2_kegg[c2_kegg$gs_name == "KEGG_P53_SIGNALING_PATHWAY", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol) #pasar a caracter

#filtrar Genes de Interés (P53 SIGNALING PATHWAY)
selected_genes_indices<- which(rownames(tnbc_expressiondata) %in% genes_of_interest_vector)
selected_genes_expr_tnbc <- tnbc_expressiondata[selected_genes_indices, ]

#verificar coherencia de las dimensiones
stopifnot(identical(colnames(tnbc_expressiondata), rownames(tnbc_clinicaldata)))

# Calcular la firma de expresión promedio
expr_signature_tnbc <- colMeans(selected_genes_expr_tnbc, na.rm = TRUE)

#agrupar Pacientes por Expresión de Firma Génica
group_tnbc <- ifelse(
  expr_signature_tnbc > median(expr_signature_tnbc, na.rm = TRUE),
  "Alta", "Baja"
)

#verificar distribución de grupos
table(group_tnbc)

#asegurar que las longitudes coinciden
stopifnot(length(group_tnbc) == nrow(tnbc_clinicaldata))

#preparar datos para el análisis de supervivencia
dfsurv_tnbc <- data.frame(
  time = as.numeric(tnbc_clinicaldata$time),
  status = as.numeric(tnbc_clinicaldata$event),
  group = group_tnbc
)

#verificar valores NA
sum(is.na(dfsurv_tnbc$time))
sum(is.na(dfsurv_tnbc$status))

#análisis de Supervivencia Kaplan-Meier
fit_tnbc <- survfit(Surv(time = dfsurv_tnbc$time, event = dfsurv_tnbc$status) ~ group, data = dfsurv_tnbc)

# Graficar la Curva Kaplan-Meier
ggsurvplot(
  fit_tnbc, data = dfsurv_tnbc,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - TNBC (P53 SIGNALING PATHWAY)",
  risk.table = FALSE
)

#guardar genes de interés para la App
genes_p53 <- data.frame(Gene_Symbol = genes_of_interest_vector)
write.csv(genes_p53, file = "genes_p53_er.csv", row.names = FALSE)


######################################################
#---Analisis de supervivencia - Cell cycle pathway---#
######################################################
# Cargar librerías necesarias
library(survival)
library(survminer)
library(dplyr)
library(msigdbr)

#Extraer datos clínicos y de expresión desde internal_data
datosclinicos<-as.data.frame(colData(internal_data))

datosdeexpresion<-t(scale(t(assay(internal_data, "normalized")))) #extraer datos de expresión  escalados 

#Comprobar que las columnas necesarias existen
required_cols<- c("HER2_status", "ER_status", "time", "event")
if (!all(required_cols %in% colnames(datosclinicos))) {
  stop("Las columnas necesarias no están presentes en los datos clínicos.")
}
#todos los pacientes
#Filtrar Genes de Interés (Cell cycle pathway)
library(msigdbr)

#descargar firmas génicas KEGG_CELL_CYCLE
c2_gene_sets<- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg<- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest<- c2_kegg[c2_kegg$gs_name == "KEGG_CELL_CYCLE", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol) #pasar a caracter

#filtrar genes de interés en los datos de expresión
selected_genes_indices <- which(rownames(datosdeexpresion) %in% genes_of_interest_vector)
selected_genes_expr <- datosdeexpresion[selected_genes_indices, ]

#calcular expresión promedio de la firma
expr_signature <- colMeans(selected_genes_expr, na.rm = TRUE)

#agrupar pacientes en base al nivel de expresion de la firma
group <- ifelse(
  expr_signature > median(expr_signature, na.rm = TRUE),
  "Alta", "Baja"
)

#agregar grupo al dataframe clínico
datosclinicos$group <- group

#verificar distribución de grupos
print(table(datosclinicos$group))

#asegurar que las longitudes coinciden
stopifnot(length(group) == nrow(datosclinicos))

#preparar el dataframe de supervivencia
dfsurv <- data.frame(
  time = as.numeric(datosclinicos$time),
  status = as.numeric(datosclinicos$event),
  group = datosclinicos$group
)

#eliminar filas con valores NA y duplicados
dfsurv <- na.omit(dfsurv)
dfsurv <- dfsurv[!duplicated(dfsurv[c("time", "status")]), ]

#comprobar valores NA restantes
print(paste("NA en time:", sum(is.na(dfsurv$time))))
print(paste("NA en status:", sum(is.na(dfsurv$status))))

# Cargar librerías necesarias
library(survival)
library(survminer)

#crear objeto de supervivencia
surv_obj <- Surv(time = dfsurv$time, event = dfsurv$status)

#ajustar modelo Kaplan-Meier
fit <- survfit(surv_obj ~ group, data = dfsurv)

#graficar la curva Kaplan-Meier

ggsurvplot(
  fit, data = dfsurv,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - Todos (CELL CYCLE)",
  risk.table = FALSE
)

#filtrar muestras por subtipo de cancer

#HER2

her2_samples <- datosclinicos$HER2_status == "Positive" #selecciona Her2 positivos - cambiar de acuerdo al subtipo de interes
her2_expressiondata <- datosdeexpresion[, her2_samples]
her2_clinicaldata <- datosclinicos[her2_samples, ]

#Filtrar Genes de Interés (KEGG_CELL_CYCLE)
library(msigdbr)

#descargar firmas génicas KEGG_CELL_CYCLE
c2_gene_sets<- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg<- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest<- c2_kegg[c2_kegg$gs_name == "KEGG_CELL_CYCLE", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol) #pasar a caracter

#filtrar genes de interés en los datos de expresión
selected_genes_indices <- which(rownames(her2_expressiondata) %in% genes_of_interest_vector)
selected_genes_expr <- her2_expressiondata[selected_genes_indices, ]

#verificar coherencia de las dimensiones
stopifnot(identical(colnames(her2_expressiondata), rownames(her2_clinicaldata)))

#calcular expresion promedio de la firma
expr_signature_her2 <- colMeans(selected_genes_expr, na.rm = TRUE)

#agrupar pacientes por nivel de expresión de la firma
group_her2 <- ifelse(
  expr_signature_her2 > median(expr_signature_her2, na.rm = TRUE),
  "Alta", "Baja"
)

#añadir el grupo al dataframe clínico
her2_clinicaldata$group <- group_her2

# Verificar distribución
table(her2_clinicaldata$group)

#asegurar que las longitudes coinciden
stopifnot(length(group_her2) == nrow(her2_clinicaldata))

#preparar el dataframe para el análisis de supervivencia
dfsurv_her2 <- data.frame(
  time = as.numeric(her2_clinicaldata$time),
  status = as.numeric(her2_clinicaldata$event),
  group = her2_clinicaldata$group
)
#eliminar NA y duplicados en columnas clave
dfsurv_her2 <- na.omit(dfsurv_her2)
dfsurv_her2 <- dfsurv_her2[!duplicated(dfsurv_her2[c("time", "status")]), ]

#verificar valores NA
sum(is.na(dfsurv_her2$time))
sum(is.na(dfsurv_her2$status))

#análisis de Supervivencia Kaplan-Meier
surv_obj_her2 <- Surv(time = dfsurv_her2$time, event = dfsurv_her2$status) #objeto de supervivencia
fit_her2 <- survfit(surv_obj_her2 ~ group, data = dfsurv_her2)

# Graficar la Curva Kaplan-Meier
ggsurvplot(
  fit_her2, data = dfsurv_her2,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - HER2+ (CELL CYCLE)",
  risk.table = FALSE
)

#filtrar muestras: ER+
er_samples <- datosclinicos$ER_status == "Positive"
er_expressiondata <- datosdeexpresion[, er_samples]
er_clinicaldata <- datosclinicos[er_samples, ]

#filtrar genes de Interés (cell cycle signaling pathway)
library(msigdbr)

#descargar firmas génicas KEGG_CELL_CYCLE
c2_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg <- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest <- c2_kegg[c2_kegg$gs_name == "KEGG_CELL_CYCLE", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol)

#filtrar genes de interés en los datos de expresión
selected_genes_indices<- which(rownames(er_expressiondata) %in% genes_of_interest_vector)
selected_genes_expr<- er_expressiondata[selected_genes_indices, ]

#verificar coherencia de las dimensiones
stopifnot(identical(colnames(er_expressiondata), rownames(er_clinicaldata)))

#calcular la firma de expresión promedio
expr_signature_er <- colMeans(selected_genes_expr, na.rm = TRUE)

#agrupar pacientes por expresión de firma génica
group_er <- ifelse(
  expr_signature_er > median(expr_signature_er, na.rm = TRUE),
  "Alta", "Baja"
)

#añadir el grupo al dataframe clínico
er_clinicaldata$group <- group_er

# Verificar distribución
table(er_clinicaldata$group)

#asegurar que las longitudes coinciden
stopifnot(length(group_er) == nrow(er_clinicaldata))

#preparar el dataframe para el análisis de supervivencia
dfsurv_er <- data.frame(
  time = as.numeric(er_clinicaldata$time),
  status = as.numeric(er_clinicaldata$event),
  group = er_clinicaldata$group
)
#eliminar NA y duplicados en columnas clave
dfsurv_er <- na.omit(dfsurv_er)
dfsurv_er <- dfsurv_er[!duplicated(dfsurv_er[c("time", "status")]), ]

#verificar valores NA
sum(is.na(dfsurv_er$time))
sum(is.na(dfsurv_er$status))


#análisis de Supervivencia Kaplan-Meier
fit_er <- survfit(Surv(time = dfsurv_er$time, event = dfsurv_er$status) ~ group, data = dfsurv_er)

# Graficar la Curva Kaplan-Meier
ggsurvplot(
  fit_er, data = dfsurv_er,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - ER+ (CELL CYCLE)",
  risk.table = FALSE
)

#filtrar Muestras: TNBC (HER2-, ER-, PR-)
tnbc_samples <- datosclinicos$HER2_status == "Negative" & 
  datosclinicos$ER_status == "Negative" & 
  datosclinicos$PR_status == "Negative"
tnbc_expressiondata <- datosdeexpresion[, tnbc_samples]
tnbc_clinicaldata <- datosclinicos[tnbc_samples, ]

#descargar firmas génicas CELL CYCLE
c2_gene_sets<- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg<- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest<- c2_kegg[c2_kegg$gs_name == "KEGG_CELL_CYCLE", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol) #pasar a caracter

#filtrar Genes de Interés (KEGG_CELL_CYCLE)
selected_genes_indices<- which(rownames(tnbc_expressiondata) %in% genes_of_interest_vector)
selected_genes_expr_tnbc <- tnbc_expressiondata[selected_genes_indices, ]

#verificar coherencia de las dimensiones
stopifnot(identical(colnames(tnbc_expressiondata), rownames(tnbc_clinicaldata)))

#calcular la expresion promedio de la firma
expr_signature_tnbc <- colMeans(selected_genes_expr_tnbc, na.rm = TRUE)

#agrupar Pacientes por Expresión de Firma Génica
group_tnbc <- ifelse(
  expr_signature_tnbc > median(expr_signature_tnbc, na.rm = TRUE),
  "Alta", "Baja"
)

#añadir el grupo al dataframe clínico
tnbc_clinicaldata$group <- group_tnbc

# Verificar distribución
table(tnbc_clinicaldata$group)

#asegurar que las longitudes coinciden
stopifnot(length(group_tnbc) == nrow(tnbc_clinicaldata))

#preparar el dataframe para el análisis de supervivencia
dfsurv_tnbc <- data.frame(
  time = as.numeric(tnbc_clinicaldata$time),
  status = as.numeric(tnbc_clinicaldata$event),
  group = tnbc_clinicaldata$group
)
#eliminar NA y duplicados en columnas clave
dfsurv_tnbc <- na.omit(dfsurv_tnbc)
dfsurv_tnbc <- dfsurv_tnbc[!duplicated(dfsurv_tnbc[c("time", "status")]), ]

#verificar valores NA
sum(is.na(dfsurv_tnbc$time))
sum(is.na(dfsurv_tnbc$status))

#análisis de Supervivencia Kaplan-Meier
fit_tnbc <- survfit(Surv(time = dfsurv_tnbc$time, event = dfsurv_tnbc$status) ~ group, data = dfsurv_tnbc)

#graficar la Curva Kaplan-Meier
ggsurvplot(
  fit_tnbc, data = dfsurv_tnbc,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - TNBC (CELL CYCLE)",
  risk.table = FALSE
)

#Guardar Genes de Interés para la App
genes_cell_cycle <- data.frame(Gene_Symbol = genes_of_interest_vector)
write.csv(genes_cell_cycle, file = "genes_cell_cycle.csv", row.names = FALSE)


#####################################################
#-Analisis de supervivencia-DNA replication pathway-#
#####################################################

# Cargar librerías necesarias
library(survival)
library(survminer)
library(dplyr)
library(msigdbr)

#Extraer datos clínicos y de expresión desde internal_data
datosclinicos<-as.data.frame(colData(internal_data))

datosdeexpresion<-t(scale(t(assay(internal_data, "normalized")))) #extraer datos de expresión  escalados 

#Comprobar que las columnas necesarias existen
required_cols<- c("HER2_status", "ER_status", "time", "event")
if (!all(required_cols %in% colnames(datosclinicos))) {
  stop("Las columnas necesarias no están presentes en los datos clínicos.")
}


#todos los pacientes
#Filtrar Genes de Interés (KEGG_DNA_REPLICATION)
library(msigdbr)

#descargar firmas génicas KEGG_DNA_REPLICATION
c2_gene_sets<- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg<- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest<- c2_kegg[c2_kegg$gs_name == "KEGG_DNA_REPLICATION", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol) #pasar a caracter

#filtrar genes de interés en los datos de expresión
selected_genes_indices <- which(rownames(datosdeexpresion) %in% genes_of_interest_vector)
selected_genes_expr <- datosdeexpresion[selected_genes_indices, ]

#calcular expresión promedio de la firma
expr_signature <- colMeans(selected_genes_expr, na.rm = TRUE)

#agrupar pacientes en base al nivel de expresion de la firma
group <- ifelse(
  expr_signature > median(expr_signature, na.rm = TRUE),
  "Alta", "Baja"
)

#agregar grupo al dataframe clínico
datosclinicos$group <- group

#verificar distribución de grupos
print(table(datosclinicos$group))

#asegurar que las longitudes coinciden
stopifnot(length(group) == nrow(datosclinicos))

#preparar el dataframe de supervivencia
dfsurv <- data.frame(
  time = as.numeric(datosclinicos$time),
  status = as.numeric(datosclinicos$event),
  group = datosclinicos$group
)

#eliminar filas con valores NA y duplicados
dfsurv <- na.omit(dfsurv)
dfsurv <- dfsurv[!duplicated(dfsurv[c("time", "status")]), ]

#comprobar valores NA restantes
print(paste("NA en time:", sum(is.na(dfsurv$time))))
print(paste("NA en status:", sum(is.na(dfsurv$status))))

# Cargar librerías necesarias
library(survival)
library(survminer)

#crear objeto de supervivencia
surv_obj <- Surv(time = dfsurv$time, event = dfsurv$status)

#ajustar modelo Kaplan-Meier
fit <- survfit(surv_obj ~ group, data = dfsurv)

#graficar la curva Kaplan-Meier
ggsurvplot(
  fit, data = dfsurv,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - Todos (DNA Replication)",
  risk.table = FALSE
)

#filtrar muestras por subtipo de cancer

#HER2

her2_samples <- datosclinicos$HER2_status == "Positive" #selecciona Her2 positivos - cambiar de acuerdo al subtipo de interes
her2_expressiondata <- datosdeexpresion[, her2_samples]
her2_clinicaldata <- datosclinicos[her2_samples, ]

#Filtrar Genes de Interés (KEGG_DNA_REPLICATION)
library(msigdbr)

#descargar firmas génicas KEGG_DNA_REPLICATION
c2_gene_sets<- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg<- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest<- c2_kegg[c2_kegg$gs_name == "KEGG_DNA_REPLICATION", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol) #pasar a caracter

#filtrar genes de interés en los datos de expresión
selected_genes_indices <- which(rownames(her2_expressiondata) %in% genes_of_interest_vector)
selected_genes_expr <- her2_expressiondata[selected_genes_indices, ]

#verificar coherencia de las dimensiones
stopifnot(identical(colnames(her2_expressiondata), rownames(her2_clinicaldata)))

#calcular expresion promedio de la firma
expr_signature_her2 <- colMeans(selected_genes_expr, na.rm = TRUE)

#agrupar pacientes por nivel de expresión de la firma
group_her2 <- ifelse(
  expr_signature_her2 > median(expr_signature_her2, na.rm = TRUE),
  "Alta", "Baja"
)

#añadir el grupo al dataframe clínico
her2_clinicaldata$group <- group_her2

# Verificar distribución
table(her2_clinicaldata$group)

#asegurar que las longitudes coinciden
stopifnot(length(group_her2) == nrow(her2_clinicaldata))

#preparar el dataframe para el análisis de supervivencia
dfsurv_her2 <- data.frame(
  time = as.numeric(her2_clinicaldata$time),
  status = as.numeric(her2_clinicaldata$event),
  group = her2_clinicaldata$group
)
#eliminar NA y duplicados en columnas clave
dfsurv_her2 <- na.omit(dfsurv_her2)
dfsurv_her2 <- dfsurv_her2[!duplicated(dfsurv_her2[c("time", "status")]), ]

#verificar valores NA
sum(is.na(dfsurv_her2$time))
sum(is.na(dfsurv_her2$status))

#análisis de Supervivencia Kaplan-Meier
surv_obj_her2 <- Surv(time = dfsurv_her2$time, event = dfsurv_her2$status) #objeto de supervivencia
fit_her2 <- survfit(surv_obj_her2 ~ group, data = dfsurv_her2)

#graficar la Curva Kaplan-Meier
ggsurvplot(
  fit_her2, data = dfsurv_her2,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - HER2+ (DNA Replication)",
  risk.table = FALSE
)

#filtrar muestras: ER+
er_samples <- datosclinicos$ER_status == "Positive"
er_expressiondata <- datosdeexpresion[, er_samples]
er_clinicaldata <- datosclinicos[er_samples, ]

#filtrar genes de Interés (KEGG_DNA_REPLICATION)
library(msigdbr)

#descargar firmas génicas KEGG_DNA_REPLICATION
c2_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg <- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest <- c2_kegg[c2_kegg$gs_name == "KEGG_DNA_REPLICATION", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol)

#filtrar genes de interés en los datos de expresión
selected_genes_indices<- which(rownames(er_expressiondata) %in% genes_of_interest_vector)
selected_genes_expr<- er_expressiondata[selected_genes_indices, ]

#verificar coherencia de las dimensiones
stopifnot(identical(colnames(er_expressiondata), rownames(er_clinicaldata)))

#calcular la expresion promerdio de la firma
expr_signature_er <- colMeans(selected_genes_expr, na.rm = TRUE)

#agrupar pacientes por expresión de firma génica
group_er <- ifelse(
  expr_signature_er > median(expr_signature_er, na.rm = TRUE),
  "Alta", "Baja"
)

#añadir el grupo al dataframe clínico
er_clinicaldata$group <- group_er

# Verificar distribución
table(er_clinicaldata$group)

#asegurar que las longitudes coinciden
stopifnot(length(group_er) == nrow(er_clinicaldata))

#preparar el dataframe para el análisis de supervivencia
dfsurv_er <- data.frame(
  time = as.numeric(er_clinicaldata$time),
  status = as.numeric(er_clinicaldata$event),
  group = er_clinicaldata$group
)
#eliminar NA y duplicados en columnas clave
dfsurv_er <- na.omit(dfsurv_er)
dfsurv_er <- dfsurv_er[!duplicated(dfsurv_er[c("time", "status")]), ]

#verificar valores NA
sum(is.na(dfsurv_er$time))
sum(is.na(dfsurv_er$status))

#análisis de Supervivencia Kaplan-Meier
fit_er <- survfit(Surv(time = dfsurv_er$time, event = dfsurv_er$status) ~ group, data = dfsurv_er)

# Graficar la Curva Kaplan-Meier
ggsurvplot(
  fit_er, data = dfsurv_er,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - ER+ (DNA Replication)",
  risk.table = FALSE
)

#filtrar Muestras: TNBC (HER2-, ER-, PR-)
tnbc_samples <- datosclinicos$HER2_status == "Negative" & 
  datosclinicos$ER_status == "Negative" & 
  datosclinicos$PR_status == "Negative"
tnbc_expressiondata <- datosdeexpresion[, tnbc_samples]
tnbc_clinicaldata <- datosclinicos[tnbc_samples, ]

#descargar firmas génicas DNA replication
c2_gene_sets<- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg<- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest<- c2_kegg[c2_kegg$gs_name == "KEGG_DNA_REPLICATION", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol) #pasar a caracter

#filtrar Genes de Interés (KEGG_DNA_REPLICATION)
selected_genes_indices<- which(rownames(tnbc_expressiondata) %in% genes_of_interest_vector)
selected_genes_expr_tnbc <- tnbc_expressiondata[selected_genes_indices, ]

#verificar coherencia de las dimensiones
stopifnot(identical(colnames(tnbc_expressiondata), rownames(tnbc_clinicaldata)))

#calcular la expresion promedio de la firma
expr_signature_tnbc <- colMeans(selected_genes_expr_tnbc, na.rm = TRUE)

#agrupar Pacientes por Expresión de Firma Génica
group_tnbc <- ifelse(
  expr_signature_tnbc > median(expr_signature_tnbc, na.rm = TRUE),
  "Alta", "Baja"
)

#añadir el grupo al dataframe clínico
tnbc_clinicaldata$group <- group_tnbc

# Verificar distribución
table(tnbc_clinicaldata$group)

#asegurar que las longitudes coinciden
stopifnot(length(group_tnbc) == nrow(tnbc_clinicaldata))

#preparar el dataframe para el análisis de supervivencia
dfsurv_tnbc <- data.frame(
  time = as.numeric(tnbc_clinicaldata$time),
  status = as.numeric(tnbc_clinicaldata$event),
  group = tnbc_clinicaldata$group
)
#eliminar NA y duplicados en columnas clave
dfsurv_tnbc <- na.omit(dfsurv_tnbc)
dfsurv_tnbc <- dfsurv_tnbc[!duplicated(dfsurv_tnbc[c("time", "status")]), ]

#verificar valores NA
sum(is.na(dfsurv_tnbc$time))
sum(is.na(dfsurv_tnbc$status))

#análisis de Supervivencia Kaplan-Meier
fit_tnbc <- survfit(Surv(time = dfsurv_tnbc$time, event = dfsurv_tnbc$status) ~ group, data = dfsurv_tnbc)

#graficar la Curva Kaplan-Meier
ggsurvplot(
  fit_tnbc, data = dfsurv_tnbc,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - TNBC (DNA Replication)",
  risk.table = FALSE
)

#Guardar genes de interés para la App
genes_dna_replication <- data.frame(Gene_Symbol = genes_of_interest_vector)
write.csv(genes_dna_replication, file = "genes_dna_replication.csv", row.names = FALSE)

########################################################
#-Analisis de supervivencia-TGF Beta signaling pathway-#
########################################################

# Cargar librerías necesarias
library(survival)
library(survminer)
library(dplyr)
library(msigdbr)

#Extraer datos clínicos y de expresión desde internal_data
datosclinicos<-as.data.frame(colData(internal_data))

datosdeexpresion<-t(scale(t(assay(internal_data, "normalized")))) #extraer datos de expresión  escalados 

#Comprobar que las columnas necesarias existen
required_cols<- c("HER2_status", "ER_status", "time", "event")
if (!all(required_cols %in% colnames(datosclinicos))) {
  stop("Las columnas necesarias no están presentes en los datos clínicos.")
}

#todos los pacientes
#Filtrar Genes de Interés (KEGG_TGF_BETA_SIGNALING_PATHWAY)
library(msigdbr)

#descargar firmas génicas KEGG_TGF_BETA_SIGNALING_PATHWAY
c2_gene_sets<- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg<- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest<- c2_kegg[c2_kegg$gs_name == "KEGG_TGF_BETA_SIGNALING_PATHWAY", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol) #pasar a caracter

#filtrar genes de interés en los datos de expresión
selected_genes_indices <- which(rownames(datosdeexpresion) %in% genes_of_interest_vector)
selected_genes_expr <- datosdeexpresion[selected_genes_indices, ]

#calcular expresión promedio de la firma
expr_signature <- colMeans(selected_genes_expr, na.rm = TRUE)

#agrupar pacientes en base al nivel de expresion de la firma
group <- ifelse(
  expr_signature > median(expr_signature, na.rm = TRUE),
  "Alta", "Baja"
)

#agregar grupo al dataframe clínico
datosclinicos$group <- group

#verificar distribución de grupos
print(table(datosclinicos$group))

#asegurar que las longitudes coinciden
stopifnot(length(group) == nrow(datosclinicos))

#preparar el dataframe de supervivencia
dfsurv <- data.frame(
  time = as.numeric(datosclinicos$time),
  status = as.numeric(datosclinicos$event),
  group = datosclinicos$group
)

#eliminar filas con valores NA y duplicados
dfsurv <- na.omit(dfsurv)
dfsurv <- dfsurv[!duplicated(dfsurv[c("time", "status")]), ]

#comprobar valores NA restantes
print(paste("NA en time:", sum(is.na(dfsurv$time))))
print(paste("NA en status:", sum(is.na(dfsurv$status))))

# Cargar librerías necesarias
library(survival)
library(survminer)

#crear objeto de supervivencia
surv_obj <- Surv(time = dfsurv$time, event = dfsurv$status)

#ajustar modelo Kaplan-Meier
fit <- survfit(surv_obj ~ group, data = dfsurv)

#graficar la curva Kaplan-Meier
ggsurvplot(
  fit, data = dfsurv,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - Todos (TGF-B signaling pathway)",
  risk.table = FALSE
)

#filtrar muestras por subtipo de cancer

#HER2

her2_samples <- datosclinicos$HER2_status == "Positive" #selecciona Her2 positivos - cambiar de acuerdo al subtipo de interes
her2_expressiondata <- datosdeexpresion[, her2_samples]
her2_clinicaldata <- datosclinicos[her2_samples, ]

#Filtrar Genes de Interés (TGF-B signaling pathway)
library(msigdbr)

#descargar firmas génicas TGF-B signaling pathway
c2_gene_sets<- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg<- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest<- c2_kegg[c2_kegg$gs_name == "KEGG_TGF_BETA_SIGNALING_PATHWAY", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol) #pasar a caracter

#filtrar genes de interés en los datos de expresión
selected_genes_indices <- which(rownames(her2_expressiondata) %in% genes_of_interest_vector)
selected_genes_expr <- her2_expressiondata[selected_genes_indices, ]

#verificar coherencia de las dimensiones
stopifnot(identical(colnames(her2_expressiondata), rownames(her2_clinicaldata)))

#calcular expresion promedio de la firma
expr_signature_her2 <- colMeans(selected_genes_expr, na.rm = TRUE)

#agrupar pacientes por nivel de expresión de la firma
group_her2 <- ifelse(
  expr_signature_her2 > median(expr_signature_her2, na.rm = TRUE),
  "Alta", "Baja"
)

#añadir el grupo al dataframe clínico
her2_clinicaldata$group <- group_her2

# Verificar distribución
table(her2_clinicaldata$group)

#asegurar que las longitudes coinciden
stopifnot(length(group_her2) == nrow(her2_clinicaldata))

#preparar el dataframe para el análisis de supervivencia
dfsurv_her2 <- data.frame(
  time = as.numeric(her2_clinicaldata$time),
  status = as.numeric(her2_clinicaldata$event),
  group = her2_clinicaldata$group
)
#eliminar NA y duplicados en columnas clave
dfsurv_her2 <- na.omit(dfsurv_her2)
dfsurv_her2 <- dfsurv_her2[!duplicated(dfsurv_her2[c("time", "status")]), ]

#verificar valores NA
sum(is.na(dfsurv_her2$time))
sum(is.na(dfsurv_her2$status))

#análisis de Supervivencia Kaplan-Meier
surv_obj_her2 <- Surv(time = dfsurv_her2$time, event = dfsurv_her2$status) #objeto de supervivencia
fit_her2 <- survfit(surv_obj_her2 ~ group, data = dfsurv_her2)

#graficar la Curva Kaplan-Meier
ggsurvplot(
  fit_her2, data = dfsurv_her2,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - HER2+ (TGF-B signaling pathway)",
  risk.table = FALSE
)


#filtrar muestras: ER+
er_samples <- datosclinicos$ER_status == "Positive"
er_expressiondata <- datosdeexpresion[, er_samples]
er_clinicaldata <- datosclinicos[er_samples, ]

#filtrar genes de Interés (TGFB signaling pathway)
library(msigdbr)

#descargar firmas génicas KEGG_TGF_BETA_SIGNALING_PATHWAY
c2_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg <- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest <- c2_kegg[c2_kegg$gs_name == "KEGG_TGF_BETA_SIGNALING_PATHWAY", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol)

#filtrar genes de interés en los datos de expresión
selected_genes_indices<- which(rownames(er_expressiondata) %in% genes_of_interest_vector)
selected_genes_expr<- er_expressiondata[selected_genes_indices, ]

#verificar coherencia de las dimensiones
stopifnot(identical(colnames(er_expressiondata), rownames(er_clinicaldata)))

#calcular la expresion promedio de la firma
expr_signature_er <- colMeans(selected_genes_expr, na.rm = TRUE)

#agrupar pacientes por expresión de firma génica
group_er <- ifelse(
  expr_signature_er > median(expr_signature_er, na.rm = TRUE),
  "Alta", "Baja"
)

#añadir el grupo al dataframe clínico
er_clinicaldata$group <- group_er

# Verificar distribución
table(er_clinicaldata$group)

#asegurar que las longitudes coinciden
stopifnot(length(group_er) == nrow(er_clinicaldata))

#preparar el dataframe para el análisis de supervivencia
dfsurv_er <- data.frame(
  time = as.numeric(er_clinicaldata$time),
  status = as.numeric(er_clinicaldata$event),
  group = er_clinicaldata$group
)
#eliminar NA y duplicados en columnas clave
dfsurv_er <- na.omit(dfsurv_er)
dfsurv_er <- dfsurv_er[!duplicated(dfsurv_er[c("time", "status")]), ]

#verificar valores NA
sum(is.na(dfsurv_er$time))
sum(is.na(dfsurv_er$status))

#análisis de Supervivencia Kaplan-Meier
fit_er <- survfit(Surv(time = dfsurv_er$time, event = dfsurv_er$status) ~ group, data = dfsurv_er)

# Graficar la Curva Kaplan-Meier
ggsurvplot(
  fit_er, data = dfsurv_er,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - ER+ (TGF-B signaling pathway)",
  risk.table = FALSE
)

#filtrar Muestras: TNBC (HER2-, ER-, PR-)
tnbc_samples <- datosclinicos$HER2_status == "Negative" & 
  datosclinicos$ER_status == "Negative" & 
  datosclinicos$PR_status == "Negative"
tnbc_expressiondata <- datosdeexpresion[, tnbc_samples]
tnbc_clinicaldata <- datosclinicos[tnbc_samples, ]

#descargar firmas génicas KEGG_TGF_BETA_SIGNALING_PATHWAY
c2_gene_sets<- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg<- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest<- c2_kegg[c2_kegg$gs_name == "KEGG_TGF_BETA_SIGNALING_PATHWAY", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol) #pasar a caracter

#filtrar Genes de Interés (KEGG_TGF_BETA_SIGNALING_PATHWAY)
selected_genes_indices<- which(rownames(tnbc_expressiondata) %in% genes_of_interest_vector)
selected_genes_expr_tnbc <- tnbc_expressiondata[selected_genes_indices, ]

#verificar coherencia de las dimensiones
stopifnot(identical(colnames(tnbc_expressiondata), rownames(tnbc_clinicaldata)))

# Calcular la expresion promedio de la firma
expr_signature_tnbc <- colMeans(selected_genes_expr_tnbc, na.rm = TRUE)

#agrupar pacientes por expresión de firma génica
group_tnbc <- ifelse(
  expr_signature_tnbc > median(expr_signature_tnbc, na.rm = TRUE),
  "Alta", "Baja"
)

#añadir el grupo al dataframe clínico
tnbc_clinicaldata$group <- group_tnbc

# Verificar distribución
table(tnbc_clinicaldata$group)

#asegurar que las longitudes coinciden
stopifnot(length(group_tnbc) == nrow(tnbc_clinicaldata))

#preparar el dataframe para el análisis de supervivencia
dfsurv_tnbc <- data.frame(
  time = as.numeric(tnbc_clinicaldata$time),
  status = as.numeric(tnbc_clinicaldata$event),
  group = tnbc_clinicaldata$group
)
#eliminar NA y duplicados en columnas clave
dfsurv_tnbc <- na.omit(dfsurv_tnbc)
dfsurv_tnbc <- dfsurv_tnbc[!duplicated(dfsurv_tnbc[c("time", "status")]), ]

#verificar valores NA
sum(is.na(dfsurv_tnbc$time))
sum(is.na(dfsurv_tnbc$status))

#análisis de Supervivencia Kaplan-Meier
fit_tnbc <- survfit(Surv(time = dfsurv_tnbc$time, event = dfsurv_tnbc$status) ~ group, data = dfsurv_tnbc)

#graficar la Curva Kaplan-Meier
ggsurvplot(
  fit_tnbc, data = dfsurv_tnbc,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - TNBC (TGF-B signaling pathway)",
  risk.table = FALSE
)

#guardar genes de interés para la App
genes_tgf_beta <- data.frame(Gene_Symbol = genes_of_interest_vector)
write.csv(genes_tgf_beta, file = "genes_tgf_beta_signaling.csv", row.names = FALSE)

##########################################################
#-Analisis de Supervivencia - JAK-STAT signaling pathway-#
##########################################################

# Cargar librerías necesarias
library(survival)
library(survminer)
library(dplyr)
library(msigdbr)

#Extraer datos clínicos y de expresión desde internal_data
datosclinicos<-as.data.frame(colData(internal_data))

datosdeexpresion<-t(scale(t(assay(internal_data, "normalized")))) #extraer datos de expresión  escalados 

#Comprobar que las columnas necesarias existen
required_cols<- c("HER2_status", "ER_status", "time", "event")
if (!all(required_cols %in% colnames(datosclinicos))) {
  stop("Las columnas necesarias no están presentes en los datos clínicos.")
}

#todos los pacientes
#Filtrar Genes de Interés (KEGG_JAK_STAT_SIGNALING_PATHWAY)
library(msigdbr)

#descargar firmas génicas KEGG_JAK_STAT_SIGNALING_PATHWAY
c2_gene_sets<- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg<- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest<- c2_kegg[c2_kegg$gs_name == "KEGG_JAK_STAT_SIGNALING_PATHWAY", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol) #pasar a caracter

#filtrar genes de interés en los datos de expresión
selected_genes_indices <- which(rownames(datosdeexpresion) %in% genes_of_interest_vector)
selected_genes_expr <- datosdeexpresion[selected_genes_indices, ]

#calcular expresión promedio de la firma
expr_signature <- colMeans(selected_genes_expr, na.rm = TRUE)

#agrupar pacientes en base al nivel de expresion de la firma
group <- ifelse(
  expr_signature > median(expr_signature, na.rm = TRUE),
  "Alta", "Baja"
)

#agregar grupo al dataframe clínico
datosclinicos$group <- group

#verificar distribución de grupos
print(table(datosclinicos$group))

#asegurar que las longitudes coinciden
stopifnot(length(group) == nrow(datosclinicos))

#preparar el dataframe de supervivencia
dfsurv <- data.frame(
  time = as.numeric(datosclinicos$time),
  status = as.numeric(datosclinicos$event),
  group = datosclinicos$group
)

#eliminar filas con valores NA y duplicados
dfsurv <- na.omit(dfsurv)
dfsurv <- dfsurv[!duplicated(dfsurv[c("time", "status")]), ]

#comprobar valores NA restantes
print(paste("NA en time:", sum(is.na(dfsurv$time))))
print(paste("NA en status:", sum(is.na(dfsurv$status))))

# Cargar librerías necesarias
library(survival)
library(survminer)

#crear objeto de supervivencia
surv_obj <- Surv(time = dfsurv$time, event = dfsurv$status)

#ajustar modelo Kaplan-Meier
fit <- survfit(surv_obj ~ group, data = dfsurv)

#graficar la curva Kaplan-Meier
ggsurvplot(
  fit, data = dfsurv,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - Todos (JAK-STAT signaling pathway)",
  risk.table = FALSE
)


#filtrar muestras por subtipo de cancer
#her2
her2_samples <- datosclinicos$HER2_status == "Positive" #selecciona Her2 positivos - cambiar de acuerdo al subtipo de interes
her2_expressiondata <- datosdeexpresion[, her2_samples]
her2_clinicaldata <- datosclinicos[her2_samples, ]

#filtrar Genes de Interés (KEGG_JAK_STAT_SIGNALING_PATHWAY)
library(msigdbr)

#descargar firma génica KEGG_JAK_STAT_SIGNALING_PATHWAY
c2_gene_sets<- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg<- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest<- c2_kegg[c2_kegg$gs_name == "KEGG_JAK_STAT_SIGNALING_PATHWAY", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol)

#filtrar genes de interés en los datos de expresión
selected_genes_indices <- which(rownames(her2_expressiondata) %in% genes_of_interest_vector)
selected_genes_expr <- her2_expressiondata[selected_genes_indices, ]

#verificar coherencia de las dimensiones
stopifnot(identical(colnames(her2_expressiondata), rownames(her2_clinicaldata)))

#calcular la expresion promedio de la firma
expr_signature_her2 <- colMeans(selected_genes_expr, na.rm = TRUE)

# agrupar pacientes por expresión de firma genica
group_her2 <- ifelse(
  expr_signature_her2 > median(expr_signature_her2, na.rm = TRUE),
  "Alta", "Baja"
)

#añadir el grupo al dataframe clínico
her2_clinicaldata$group <- group_her2

# Verificar distribución
table(her2_clinicaldata$group)

#asegurar que las longitudes coinciden
stopifnot(length(group_her2) == nrow(her2_clinicaldata))

#preparar el dataframe para el análisis de supervivencia
dfsurv_her2 <- data.frame(
  time = as.numeric(her2_clinicaldata$time),
  status = as.numeric(her2_clinicaldata$event),
  group = her2_clinicaldata$group
)
#eliminar NA y duplicados en columnas clave
dfsurv_her2 <- na.omit(dfsurv_her2)
dfsurv_her2 <- dfsurv_her2[!duplicated(dfsurv_her2[c("time", "status")]), ]

#verificar valores NA
sum(is.na(dfsurv_her2$time))
sum(is.na(dfsurv_her2$status))

#análisis de Supervivencia Kaplan-Meier
surv_obj_her2 <- Surv(time = dfsurv_her2$time, event = dfsurv_her2$status) #objeto de supervivencia
fit_her2 <- survfit(surv_obj_her2 ~ group, data = dfsurv_her2)

#graficar la Curva Kaplan-Meier
ggsurvplot(
  fit_her2, data = dfsurv_her2,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - HER2+ (JAK-STAT SIGNALING PATHWAY)",
  risk.table = FALSE
)

#filtrar muestras: ER+
er_samples <- datosclinicos$ER_status == "Positive"
er_expressiondata <- datosdeexpresion[, er_samples]
er_clinicaldata <- datosclinicos[er_samples, ]

#filtrar genes de Interés (JAK-STAT signaling pathway)
library(msigdbr)

#descargar firmas génicas KEGG_JAK_STAT_SIGNALING_PATHWAY
c2_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg <- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest <- c2_kegg[c2_kegg$gs_name == "KEGG_JAK_STAT_SIGNALING_PATHWAY", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol)

#filtrar genes de interés en los datos de expresión
selected_genes_indices<- which(rownames(er_expressiondata) %in% genes_of_interest_vector)
selected_genes_expr<- er_expressiondata[selected_genes_indices, ]

#verificar coherencia de las dimensiones
stopifnot(identical(colnames(er_expressiondata), rownames(er_clinicaldata)))

#calcular la expresion promedio de la firma
expr_signature_er <- colMeans(selected_genes_expr, na.rm = TRUE)

#agrupar pacientes por expresión de firma génica
group_er <- ifelse(
  expr_signature_er > median(expr_signature_er, na.rm = TRUE),
  "Alta", "Baja"
)

#añadir el grupo al dataframe clínico
er_clinicaldata$group <- group_er

# Verificar distribución
table(er_clinicaldata$group)

#asegurar que las longitudes coinciden
stopifnot(length(group_er) == nrow(er_clinicaldata))

#preparar el dataframe para el análisis de supervivencia
dfsurv_er <- data.frame(
  time = as.numeric(er_clinicaldata$time),
  status = as.numeric(er_clinicaldata$event),
  group = er_clinicaldata$group
)
#eliminar NA y duplicados en columnas clave
dfsurv_er <- na.omit(dfsurv_er)
dfsurv_er <- dfsurv_er[!duplicated(dfsurv_er[c("time", "status")]), ]

#verificar valores NA
sum(is.na(dfsurv_er$time))
sum(is.na(dfsurv_er$status))

#análisis de Supervivencia Kaplan-Meier
fit_er <- survfit(Surv(time = dfsurv_er$time, event = dfsurv_er$status) ~ group, data = dfsurv_er)

# Graficar la Curva Kaplan-Meier
ggsurvplot(
  fit_er, data = dfsurv_er,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - ER+ (JAK-STAT signaling pathway)",
  risk.table = FALSE
)

#filtrar Muestras: TNBC (HER2-, ER-, PR-)
tnbc_samples <- datosclinicos$HER2_status == "Negative" & 
  datosclinicos$ER_status == "Negative" & 
  datosclinicos$PR_status == "Negative"
tnbc_expressiondata <- datosdeexpresion[, tnbc_samples]
tnbc_clinicaldata <- datosclinicos[tnbc_samples, ]

#descargar firmas génicas KEGG_JAK_STAT_SIGNALING_PATHWAY
c2_gene_sets<- msigdbr(species = "Homo sapiens", category = "C2")
c2_kegg<- subset(c2_gene_sets, grepl("KEGG", gs_name))
genes_of_interest<- c2_kegg[c2_kegg$gs_name == "KEGG_JAK_STAT_SIGNALING_PATHWAY", "gene_symbol"]
genes_of_interest_vector <- as.character(genes_of_interest$gene_symbol) #pasar a caracter

#filtrar Genes de Interés (KEGG_JAK_STAT_SIGNALING_PATHWAY)
selected_genes_indices<- which(rownames(tnbc_expressiondata) %in% genes_of_interest_vector)
selected_genes_expr_tnbc <- tnbc_expressiondata[selected_genes_indices, ]

#verificar coherencia de las dimensiones
stopifnot(identical(colnames(tnbc_expressiondata), rownames(tnbc_clinicaldata)))

#calcular la expresion promedio de la firma
expr_signature_tnbc <- colMeans(selected_genes_expr_tnbc, na.rm = TRUE)

#agrupar pacientes por expresión de firma génica
group_tnbc <- ifelse(
  expr_signature_tnbc > median(expr_signature_tnbc, na.rm = TRUE),
  "Alta", "Baja"
)

#añadir el grupo al dataframe clínico
tnbc_clinicaldata$group <- group_tnbc

# Verificar distribución
table(tnbc_clinicaldata$group)

#asegurar que las longitudes coinciden
stopifnot(length(group_tnbc) == nrow(tnbc_clinicaldata))

#preparar el dataframe para el análisis de supervivencia
dfsurv_tnbc <- data.frame(
  time = as.numeric(tnbc_clinicaldata$time),
  status = as.numeric(tnbc_clinicaldata$event),
  group = tnbc_clinicaldata$group
)
#eliminar NA y duplicados en columnas clave
dfsurv_tnbc <- na.omit(dfsurv_tnbc)
dfsurv_tnbc <- dfsurv_tnbc[!duplicated(dfsurv_tnbc[c("time", "status")]), ]

#verificar valores NA
sum(is.na(dfsurv_tnbc$time))
sum(is.na(dfsurv_tnbc$status))

#análisis de Supervivencia Kaplan-Meier
fit_tnbc <- survfit(Surv(time = dfsurv_tnbc$time, event = dfsurv_tnbc$status) ~ group, data = dfsurv_tnbc)

#graficar la Curva Kaplan-Meier
ggsurvplot(
  fit_tnbc, data = dfsurv_tnbc,
  pval = TRUE, conf.int = FALSE,
  title = "Curva de Supervivencia - TNBC (JAK-STAT signaling pathway)",
  risk.table = FALSE
)
#Guardar genes de interes para la App
genes_jak_stat <- data.frame(Gene_Symbol = genes_of_interest_vector)
write.csv(genes_jak_stat, file = "genes_jak_stat_signaling.csv", row.names = FALSE)








