#pruebas de enriquecimiento funcional
setwd("C:/Users/mitja/OneDrive/master bioinformatica/TFM/R_files/set de datos finales")
getwd()
#cargar los datos a internal_data
internal_data <- readRDS("datosdeexpresion.rds")
identical(rownames(colData(internal_data)), colnames(internal_data)) #coinciden


# Cargar las librerías necesarias
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(limma)

# Cargar la lista de genes desde un archivo CSV

genes <- read.csv("genes_jak_stat_signaling.csv", header = FALSE, stringsAsFactors = FALSE)[, 1] #cambiar segun la firma genica a analizar

#extraer datos de expresion y clinicos
expr_matrix <- t(scale(t(assay(internal_data, "normalized")))) #extraer datos de expresión  escalados 
column_data <- as.data.frame(colData(internal_data))

#clasificar por subtipos clínicos
column_data$clinical_group <- with(column_data,
                                   ifelse(HER2_status == "Negative" &
                                            ER_status == "Negative" &
                                            PR_status == "Negative", "TNBC",
                                          ifelse(HER2_status == "Positive", "HER2 positivo",
                                                 ifelse(ER_status == "Positive", "ER positivo", "Otro"))))

#seleccionar el subtipo de cancer
selected_clinical_group <- "ER positivo" #cambiar segun el subtipo a analizar

#filtrar datos clunicos segun subtipo seleccionado
column_data_filtered <- column_data[column_data$clinical_group == selected_clinical_group, ]

#asegurar que haya datos disponibles
if (nrow(column_data_filtered) == 0) stop("No hay datos para el grupo clínico seleccionado.")

#filtrar la matriz de expresión por pacientes seleccionados
expr_filtered <- expr_matrix[, rownames(column_data_filtered), drop = FALSE]

#asegurarse de que las muestras en expr_filtered coincidan con column_data_filtered
common_Patients <- intersect(rownames(column_data_filtered), colnames(expr_filtered))

if (length(common_Patients) == 0) stop("No hay pacientes en común entre datos clínicos y datos de expresión.")

#filtrar genes válidos
valid_genes <- intersect(rownames(expr_filtered), genes)
if (length(valid_genes) == 0) stop("No hay genes válidos en la matriz de expresión para el análisis de enriquecimiento.")

#calcular la expresión promedio de los genes seleccionados
patient_avg_expression <- colMeans(expr_filtered[valid_genes, common_Patients, drop = FALSE], na.rm = TRUE)

#asignar la expresión promedio al dataset clínico
column_data_filtered$avg_expression <- NA
column_data_filtered[common_Patients, "avg_expression"] <- patient_avg_expression[common_Patients]

#eliminar filas con valores NA
column_data_filtered <- na.omit(column_data_filtered)

#verificar si hay suficientes pacientes para continuar
if (nrow(column_data_filtered) == 0) stop("No hay suficientes pacientes para continuar con el análisis.")

#crear los grupos de expresión
column_data_filtered$expression_group <- ifelse(
  column_data_filtered$avg_expression > median(column_data_filtered$avg_expression, na.rm = TRUE),
  "Alta", "Baja"
)

#crear un diseño para el análisis con limma
design <- model.matrix(~ expression_group, data = column_data_filtered)

#realizar el análisis de expresión diferencial
fit <- lmFit(expr_filtered[valid_genes, ], design)
fit <- eBayes(fit)
res <- topTable(fit, coef = 2, number = Inf, adjust = "fdr")  # coef=2 asume el segundo nivel

#filtrar los genes significativos
sig_genes <- rownames(res)[res$adj.P.Val < 0.05]
if (length(sig_genes) == 0) stop("No se encontraron genes significativamente diferentes entre los grupos de expresión.")

#conversión de identificadores de genes
conversion <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
if (nrow(conversion) == 0) stop("No se pudieron convertir los identificadores de genes a ENTREZID.")

# Análisis de enriquecimiento funcional GO
go_enrich <- enrichGO(
  gene = conversion$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",  # para procesos biológicos
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

#verificar si se encontraron términos GO enriquecidos
if (is.null(go_enrich) || nrow(as.data.frame(go_enrich)) == 0) {
  print("No se encontraron términos GO enriquecidos.")
} else {
  #mostrar tabla de resultados
  print(head(as.data.frame(go_enrich)[, c("ID", "Description", "p.adjust")], 10))
  
  #crear gráfico de barras
  barplot(go_enrich, showCategory = 10) +
    ggtitle(paste("Top 10 términos GO enriquecidos -", selected_clinical_group)) +
    theme_minimal()
  
  #crear gráfico de puntos
  dotplot(go_enrich, showCategory = 10) +
    ggtitle(paste("Gráfico de Puntos - Top 10 términos GO -", selected_clinical_group)) +
    theme_minimal()
}


#Enriquecimiento funcional para "Todos"
# Cargar las librerías necesarias
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(limma)

# Cargar la lista de genes desde un archivo CSV
genes <- read.csv("genes_jak_stat_signaling.csv", header = FALSE, stringsAsFactors = FALSE)[, 1] 

# Extraer datos de expresión y clínicos
expr_matrix <- t(scale(t(assay(internal_data, "normalized")))) # Datos de expresión escalados
column_data <- as.data.frame(colData(internal_data)) # Datos clínicos

# Asegurarse de que hay datos disponibles
if (nrow(column_data) == 0) stop("No hay datos clínicos disponibles.")

# Asegurarse de que las muestras coincidan entre datos clínicos y datos de expresión
common_Patients <- intersect(rownames(column_data), colnames(expr_matrix))
if (length(common_Patients) == 0) stop("No hay pacientes en común entre datos clínicos y datos de expresión.")

# Filtrar la matriz de expresión y datos clínicos por pacientes comunes
expr_filtered <- expr_matrix[, common_Patients, drop = FALSE]
column_data_filtered <- column_data[common_Patients, , drop = FALSE]

# Filtrar genes válidos
valid_genes <- intersect(rownames(expr_filtered), genes)
if (length(valid_genes) == 0) stop("No hay genes válidos en la matriz de expresión para el análisis de enriquecimiento.")

# Calcular la expresión promedio de los genes seleccionados
patient_avg_expression <- colMeans(expr_filtered[valid_genes, , drop = FALSE], na.rm = TRUE)

# Asignar la expresión promedio al dataset clínico
column_data_filtered$avg_expression <- NA
column_data_filtered$avg_expression <- patient_avg_expression

# Eliminar filas con valores NA
column_data_filtered <- na.omit(column_data_filtered)

# Verificar si hay suficientes pacientes para continuar
if (nrow(column_data_filtered) == 0) stop("No hay suficientes pacientes para continuar con el análisis.")

# Crear los grupos de expresión
column_data_filtered$expression_group <- ifelse(
  column_data_filtered$avg_expression > median(column_data_filtered$avg_expression, na.rm = TRUE),
  "Alta", "Baja"
)

# Crear un diseño para el análisis con limma
design <- model.matrix(~ expression_group, data = column_data_filtered)

# Realizar el análisis de expresión diferencial
fit <- lmFit(expr_filtered[valid_genes, ], design)
fit <- eBayes(fit)
res <- topTable(fit, coef = 2, number = Inf, adjust = "fdr")  # coef=2 asume el segundo nivel

# Filtrar los genes significativos
sig_genes <- rownames(res)[res$adj.P.Val < 0.05]
if (length(sig_genes) == 0) stop("No se encontraron genes significativamente diferentes entre los grupos de expresión.")

# Conversión de identificadores de genes
conversion <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
if (nrow(conversion) == 0) stop("No se pudieron convertir los identificadores de genes a ENTREZID.")

# Análisis de enriquecimiento funcional GO
go_enrich <- enrichGO(
  gene = conversion$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",  # "BP" para procesos biológicos
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# Verificar si se encontraron términos GO enriquecidos
if (is.null(go_enrich) || nrow(as.data.frame(go_enrich)) == 0) {
  print("No se encontraron términos GO enriquecidos.")
} else {
  # Mostrar tabla de resultados
  print(head(as.data.frame(go_enrich)[, c("ID", "Description", "p.adjust")], 10))
  
  # Crear gráfico de barras
  barplot(go_enrich, showCategory = 10) +
    ggtitle("Top 10 términos GO enriquecidos") +
    theme_minimal()
  
  # Crear gráfico de puntos
  dotplot(go_enrich, showCategory = 10) +
    ggtitle("Gráfico de Puntos - Top 10 términos GO") +
    theme_minimal()
}


#################################################

#asegurar que los genes son validos - gene symbol
library(org.Hs.eg.db)
#TGFB1 SEC14L3
# Lista de símbolos de genes
genesdepurueba <- c("MKI67", "SEC14L3")

# Verificar si son válidos en la base de datos
valid_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
genesvalidos <- genesdepurueba[genesdepurueba %in% valid_symbols]

# Mostrar genes válidos
print(genesvalidos)


library(org.Hs.eg.db)

# Buscar símbolos similares en la base de datos
pattern <- "PDCD1LG2"  # Reemplaza con el símbolo no válido
possible_matches <- grep(pattern, keys(org.Hs.eg.db, keytype = "SYMBOL"), value = TRUE)

# Mostrar resultados
print(possible_matches)

# Buscar usando alias
alias_matches <- select(org.Hs.eg.db, keys = "TGFB1", keytype = "ALIAS", columns = c("SYMBOL", "GENENAME"))

# Mostrar resultados
print(alias_matches)


