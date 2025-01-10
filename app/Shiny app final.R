#Shiny app para analsisi de supervivencia y enriquecimiento funcional
#Monica Mitjans
#TFM-Shiny App

library(shiny)
library(SummarizedExperiment)
library(survival)
library(survminer)


#cargar los datos a internal_data
internal_data <- readRDS("datosdeexpresion.rds")
identical(rownames(colData(internal_data)), colnames(internal_data)) #verificar cpincidencias
scaled_matrix <- t(scale(t(assay(internal_data, "normalized"))))

#Interfaz de usuario
ui <- fluidPage(
  titlePanel("Análisis de Supervivencia y Enriquecimiento Funcional en Cáncer de Mama"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("geneList", "Subir Lista de Genes (CSV o TXT)", 
                accept = c(".csv", ".txt")),
      
      selectInput("groupVar", "Seleccionar Grupo Clínico",
                  choices = c("Todos", "TNBC", "HER2 positivo", "ER positivo"),
                  selected = "Todos"),
      
      selectInput("enrichGroup", "Seleccionar Subtipo para Enriquecimiento",
                  choices = c("Todos", "TNBC", "HER2 positivo", "ER positivo"),
                  selected = "Todos"),
      
      hr(),
      actionButton("survButton", "Generar Curvas de Supervivencia", 
                   icon = icon("chart-line"), 
                   class = "btn-primary"),
      
      actionButton("enrichButton", "Análisis de Enriquecimiento Funcional", 
                   icon = icon("chart-bar"), 
                   class = "btn-success"),
      
      hr(),
      p("Seleccione un archivo CSV o TXT con una lista de genes, el grupo clínico de interés y realice un análisis específico por subtipo de cáncer."),
      # Botón para descargar resultados de enriquecimiento
      downloadButton("downloadEnrichResults", "Descargar Resultados de Enriquecimiento")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Curvas de Supervivencia", 
                 h4("Gráfico de Curvas de Supervivencia"),
                 plotOutput("survPlot"),
                 hr(),
                 p("El gráfico muestra la comparación de supervivencia entre los diferentes niveles de expresión genética.")),
        
        tabPanel("Resultados de Enriquecimiento", 
                 h4("Tabla de Resultados de Enriquecimiento"),
                 tableOutput("enrichResults"),
                 hr(),
                 
                 h4("Gráfico de Barras de Enriquecimiento Funcional"),
                 plotOutput("enrichPlot"),
                 hr(),
                 
                 h4("Gráfico de Puntos de Enriquecimiento Funcional"),
                 plotOutput("enrichDotPlot"),
                 
                 p("Los gráficos muestran los términos GO más enriquecidos para los genes seleccionados."))
      )
    )
  )
)

#server
server <- function(input, output, session) {
  
  # Cargar lista de genes del usuario
  user_genes <- reactive({
    req(input$geneList)
    file <- input$geneList
    if (grepl("\\.csv$", file$name)) {
      genes <- read.csv(file$datapath, header = FALSE, stringsAsFactors = FALSE)[, 1]
    } else {
      genes <- read.table(file$datapath, header = FALSE, stringsAsFactors = FALSE)[, 1]
    }
    # Limpiar espacios en blanco y valores vacíos
    genes <- trimws(genes)
    genes <- genes[genes != ""]
    genes <- as.character(genes)
    
    cat("Número de genes leídos desde el archivo:", length(genes), "\n")
    cat("Primeros genes:", head(genes), "\n")
    
    unique(genes)
  })
  
  observeEvent(input$survButton, {
    req(user_genes())
    genes <- user_genes()
    
    sum_exp <- internal_data
    scaled_matrix <- t(scale(t(assay(internal_data, "normalized")))) #datos escalados
    column_data <- as.data.frame(colData(sum_exp))
    
    column_data$clinical_group <- with(column_data,
                                       ifelse(HER2_status == "Negative" &
                                                ER_status == "Negative" &
                                                PR_status == "Negative", "TNBC",
                                              ifelse(HER2_status == "Positive", "HER2 positivo",
                                                     ifelse(ER_status == "Positive", "ER positivo",
                                                            "Otro"))))
    
    selected_clinical_group <- input$groupVar
    if (selected_clinical_group != "Todos") {
      column_data_filtered <- column_data[column_data$clinical_group == selected_clinical_group, ]
    } else {
      column_data_filtered <- column_data
    }
    req(nrow(column_data_filtered) > 0, "No hay datos para el grupo clínico seleccionado.")
    
    filtered_genes <- rownames(scaled_matrix)[
      rownames(scaled_matrix) %in% genes]
    
    req(length(filtered_genes) > 0, "No hay genes en común con el conjunto interno.")
    
    patient_avg_expression <- colMeans(scaled_matrix[filtered_genes, , drop = FALSE], na.rm = TRUE)
    common_Patients <- intersect(rownames(column_data_filtered), names(patient_avg_expression))  
    
    req(length(common_Patients) > 0, "No hay pacientes en común entre datos clínicos y datos de expresión.")
    
    column_data_filtered$avg_expression <- NA
    column_data_filtered[common_Patients, "avg_expression"] <- patient_avg_expression[common_Patients]
    column_data_filtered <- na.omit(column_data_filtered)
    
    # Eliminar duplicados en 'time' y 'event'
    column_data_filtered <- column_data_filtered[!duplicated(column_data_filtered[c("time", "event")]), ]
    
    
    column_data_filtered$expression_group <- ifelse(
      column_data_filtered$avg_expression > median(column_data_filtered$avg_expression, na.rm = TRUE),
      "Alta", "Baja"
    )
    cat("Número de pacientes por grupo de expresión:\n")
    print(table(column_data_filtered$expression_group))
    print(head(column_data_filtered[, c("time", "event", "expression_group")]))
    cat("Número de pacientes en el grupo seleccionado:", nrow(column_data_filtered), "\n")
    cat("Dimensiones de col_data_filtered antes de ggsurvplot:", nrow(column_data_filtered), "\n")
    
    
    library(survival)
    library(survminer)
    # Validar que no hay valores NA en las columnas necesarias
    req(all(!is.na(column_data_filtered$time)))
    req(all(!is.na(column_data_filtered$event)))
    
    output$survPlot <- renderPlot({
      req(column_data_filtered$time, column_data_filtered$event, column_data_filtered$expression_group)
      
      # Validar que no existan valores NA en las columnas necesarias
      req(all(!is.na(column_data_filtered$time)))
      req(all(!is.na(column_data_filtered$event)))
      req(all(!is.na(column_data_filtered$expression_group)))
      
      # Ajustar el modelo de supervivencia
      
      fit <- survfit(Surv(time, event) ~ expression_group, data = column_data_filtered)
      
      cat("Número de filas en el modelo survfit:", length(fit$time), "\n")
      
      # Generar el gráfico de supervivencia
      surv_plot <- ggsurvplot(
        fit,
        data = column_data_filtered,
        pval = TRUE,
        risk.table = FALSE,
        title = paste("Supervivencia para", selected_clinical_group)
      )
      
      # Asegurarse de imprimir el gráfico
      cat("Curva de supervivencia creada", "\n" )
      print(surv_plot$plot)
    })
  })
  
  ## --------------------------------------------
  ## 🔹 Análisis de Enriquecimiento Funcional
  ## --------------------------------------------
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(limma)
  
  observeEvent(input$enrichButton, {
    req(user_genes())
    genes <- user_genes()
    
    # Extraer matriz de expresión normalizada y datos clínicos
    scaled_matrix <- t(scale(t(assay(internal_data, "normalized"))))
    column_data <- as.data.frame(colData(internal_data))
    
    # Clasificar por subtipos clínicos
    column_data$clinical_group <- with(column_data,
                                       ifelse(HER2_status == "Negative" &
                                                ER_status == "Negative" &
                                                PR_status == "Negative", "TNBC",
                                              ifelse(HER2_status == "Positive", "HER2 positivo",
                                                     ifelse(ER_status == "Positive", "ER positivo","Otros"))))
    
    selected_clinical_group <- input$enrichGroup
    
    if (selected_clinical_group != "Todos") {
      column_data_filtered <- column_data[column_data$clinical_group == selected_clinical_group, ]
    } else {
      column_data_filtered <- column_data
    }
    
    req(nrow(column_data_filtered) > 0, "No hay datos para el grupo clínico seleccionado.")
    
    # Filtrar la matriz de expresión por pacientes seleccionados
    expr_filtered <- scaled_matrix[, rownames(column_data_filtered), drop = FALSE]
    
    # Asegurarse de que las muestras en expr_filtered coincidan con column_data_filtered
    common_Patients <- intersect(rownames(column_data_filtered), colnames(expr_filtered))
    
    req(length(common_Patients) > 0, "No hay pacientes en común entre datos clínicos y datos de expresión.")
    
    # Filtrar genes válidos
    valid_genes <- intersect(rownames(expr_filtered), genes)
    req(length(valid_genes) > 0, "No hay genes válidos en la matriz de expresión para el análisis de enriquecimiento.")
    
    # Calcular la expresión promedio de los genes seleccionados
    patient_avg_expression <- colMeans(expr_filtered[valid_genes, common_Patients, drop = FALSE], na.rm = TRUE)
    
    # Asignar la expresión promedio al dataset clínico
    column_data_filtered$avg_expression <- NA
    column_data_filtered[common_Patients, "avg_expression"] <- patient_avg_expression[common_Patients]
    
    # Eliminar filas con valores NA
    column_data_filtered <- na.omit(column_data_filtered)
    
    # Verificar si hay suficientes pacientes para continuar
    req(nrow(column_data_filtered) > 0, "No hay suficientes pacientes para continuar con el análisis.")
    
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
    if (length(sig_genes) == 0) {
      showNotification(
        "No se encontraron genes significativamente diferentes entre los grupos de expresión.",
        type = "error",
        duration = 5
      )
      return()
    }
    
    req(length(sig_genes) > 0, "No se encontraron genes significativamente diferentes entre los grupos de expresión.")
    
    # Conversión de identificadores de genes
    library(org.Hs.eg.db)
    conversion <- tryCatch({
      bitr(sig_genes, fromType = "ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb = org.Hs.eg.db)
    }, error = function(e) {
      showNotification("Fallo con ENSEMBL, intentando con SYMBOL...", type = "warning", duration = 5)
      tryCatch({
        bitr(sig_genes, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
      }, error = function(e2) {
        showNotification(
          paste("Fallo con SYMBOL:", e2$message),
          type = "error",
          duration = 5
        )
        return(NULL)
      })
    })
    
    # Validar el resultado de la conversión
    if (is.null(conversion) || nrow(conversion) == 0) {
      showNotification(
        "No se pudieron convertir los identificadores de genes a ENTREZID.",
        type = "error",
        duration = 5
      )
      stop("No se encontraron identificadores válidos después de la conversión.")
    }
    
    # Mostrar genes convertidos
    print("Genes convertidos exitosamente:")
    print(conversion)
    # Análisis de enriquecimiento funcional GO
    go_enrich <- enrichGO(
      gene = conversion$ENTREZID,
      OrgDb = org.Hs.eg.db,
      ont = "BP",  # "BP" para procesos biológicos, puedes cambiar a "CC" o "MF" según necesites
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
    
    if (is.null(go_enrich) || nrow(as.data.frame(go_enrich)) == 0) {
      output$enrichResults <- renderTable({ data.frame(Mensaje = "No se encontraron términos GO enriquecidos.") })
      output$enrichPlot <- renderPlot({ plot.new(); text(0.5, 0.5, "No hay términos GO enriquecidos para mostrar.") })
      output$enrichDotPlot <- renderPlot({ plot.new(); text(0.5, 0.5, "No hay términos GO enriquecidos para mostrar.") })
      return()
    }
    
    # Mostrar tabla de resultados
    output$enrichResults <- renderTable({
      head(as.data.frame(go_enrich)[, c("ID", "Description", "p.adjust")], 10)
    })
    
    # Mostrar gráfico de barras
    output$enrichPlot <- renderPlot({
      barplot(go_enrich, showCategory = 10) +
        ggtitle(paste("Top 10 términos GO enriquecidos -", selected_clinical_group)) +
        theme_minimal()
    })
    
    # Mostrar gráfico de puntos

    output$enrichDotPlot <- renderPlot({
      dotplot(go_enrich, showCategory = 10) +
        ggtitle(paste("Gráfico de Puntos - Top 10 términos GO -", selected_clinical_group)) +
        theme_minimal()
    })
    # Habilitar descarga de resultados de enriquecimiento
    output$downloadEnrichResults <- downloadHandler(
      filename = function() { paste("enrichment_results.csv") },
      content = function(file) {
        write.csv(as.data.frame(go_enrich), file)
      }
    )
  })
  
  
}


shinyApp(ui, server)
