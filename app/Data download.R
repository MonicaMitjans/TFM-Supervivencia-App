
getwd()

library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)
library(dplyr)

#Descarga de datos de expresion
query <- GDCquery( project = "TCGA-BRCA", 
                   data.category = "Transcriptome Profiling", 
                   data.type = "Gene Expression Quantification", 
                   workflow.type = "STAR - Counts" )

GDCdownload(query)

tcga_data<- GDCprepare(query)

#guardar los datos descargados en un archivo .RData
save(tcga_data, file = "tcga_datos.RData")

#cargar los datos
load("tcga_datos.RData")

#descarga de datos clinicos- se descargan a parte
#BCR biotabs files

library(TCGAbiolinks)
query_clinical <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab"
)

GDCdownload(query_clinical)
clinical_tab_all <- GDCprepare(query_clinical)

#see all tables
names(clinical_tab_all)

# columns from clinical_patient
dplyr::glimpse(clinical_tab_all$clinical_patient_brca)

unique(clinical_tab_all$clinical_patient_brca$er_status_by_ihc)

unique(clinical_tab_all$clinical_patient_brca$her2_status_by_ihc)

nrow(rowData(tcga_data))
nrow(rowData(tcga_data))
