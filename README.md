# Herramienta Web para Análisis de Supervivencia y Enriquecimiento Funcional en Cáncer de Mama Basada en la Expresión Génica 

Esta es una aplicación web Shiny que permite realizar análisis de supervivencia y enriquecimiento funcional en cáncer de mama, con un enfoque en subtipos específicos como HER2-positivo, ER-positivo y Triple-Negativo (TNBC). 

## Funcionalidades

- **Cargar datos**: Permite cargar un archivo con una lista de genes de interés.
- **Análisis de supervivencia**: Generación de curvas de supervivencia Kaplan-Meier basadas en distintos subtipos de cáncer de mama.
- **Análisis de enriquecimiento funcional**: Realiza análisis de enriquecimiento utilizando términos GO para identificar posibles procesos biológicos asociados.

## Cómo usar la aplicación

1. **Instalación**:
   - Clona el repositorio: 
     ```bash
     git clone  https://github.com/MonicaMitjans/TFM-Supervivencia-App.git
     ```
   - Instala las dependencias necesarias:
     ```r
     install.packages("shiny")
     install.packages("survival")
     install.packages("survminer")
     install.packages("limma")
     install.packages("ClusterProfiler")
     if (!require("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
     BiocManager::install("org.Hs.eg.db")
     if (!require("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
     BiocManager::install("SummarizedExperiment")
     ```
    - Guarda el archivo datosdeexpresion.rds en el mismo directorio donde está el script de la aplicación (app.R).
    - Carga el archivo .rds
     ```r
     internal_data <- readRDS("datosdeexpresion.rds")
     ```

2. **Ejecutar la aplicación**:
   - Navega hasta el directorio donde se encuentra el archivo `app.R`.
   - Ejecuta el siguiente comando en R:
     ```r
     shiny::runApp("app.R")
     ```
3. **Sigue las intrucciones del manual de usuario**

## Contenido del repositorio

- `app.R`: Código principal de la aplicación Shiny.
- `data/`: Carpeta que contiene los archivos csv de datos de ejemplo.
- `docs/`: Manual de usuario.
- `README.md`: Este archivo.



