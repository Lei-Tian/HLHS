library(shiny)
library(Seurat)

load("./all.Robj")
options(repos = BiocManager::repositories())

ui <- navbarPage("Human fetal HLHS single-cell RNA-seq",
  navbarMenu("iPSC-EC",
    tabPanel("Metadata",
      titlePanel(h2("Metadata")),

      sidebarLayout(
        sidebarPanel(
          helpText(h3("Create uMAP colored by metadata.")),

          radioButtons("radio",
                      h4("Set color:"),
                      choices = list("Sample" = 1,
                                      "#Gene" = 2,
                                      "#UMI" = 3,
                                      "%MT" = 4,
                                      "Cluster" = 5),
                      selected = 1)
        ),
        mainPanel(plotOutput("meta"))
      )
    ),

    tabPanel("Gene expression",
      titlePanel(h2("Gene expression")),

      sidebarLayout(
        sidebarPanel(
          helpText(h3("Create uMAP colored by the
          expression of a selected gene.")),

          textInput("text",
                    h4("Enter gene name below ..."),
                    value = "CDH5"),
          br(),
          br(),

          radioButtons("radio",
                      h4("Combine HLHS and Control?"),
                      choices = list("Yes" = 1,
                                      "No" = 2),
                      selected = 1)
          ),

        mainPanel(plotOutput("gene_umap"))
      )
    ),

    tabPanel("DEG",
      titlePanel(h2("Differential expression")),

      sidebarLayout(
        sidebarPanel(
          helpText(h3("Create uMAP colored by metadata.")),

          checkboxGroupInput("checkGroup1",
                      h4("Select cluster(s) A:"),
                      choices = list("Cluster 0" = 0,
                                     "Cluster 1" = 1,
                                     "Cluster 2" = 2,
                                     "Cluster 3" = 3,
                                     "Cluster 4" = 4,
                                     "Cluster 5" = 5,
                                     "Cluster 6" = 6
                                    ),
                      selected = 0),

          checkboxGroupInput("checkGroup2",
                       h4("Select cluster(s) B:"),
                       choices = list("Cluster 0" = 0,
                                      "Cluster 1" = 1,
                                      "Cluster 2" = 2,
                                      "Cluster 3" = 3,
                                      "Cluster 4" = 4,
                                      "Cluster 5" = 5,
                                      "Cluster 6" = 6,
                                      "All others" = 7
                                      ),
                        selected = 1)
        ),
        mainPanel(tableOutput("detect_DEGs"))
      )
    )
  ),

  navbarMenu("Heart tissue",
    tabPanel("Metadata"),
    tabPanel("Gene expression"),
    tabPanel("Differential expression")
  )
)

server <- function(input, output){
  getWidth = reactive({
    if(input$radio == 1){
      return(500)
    }else{
      return(900)
    }
  })
  show_legend = FALSE
  output$meta = renderPlot({
    if(input$radio == 1){
      p = DimPlot(all, reduction.use = "umap", do.label = F, do.return = T,
                  pt.size = 2, group.by = "orig.ident")
    }else if(input$radio == 2){
      p = FeaturePlot(all, features.plot = ("nGene"), reduction.use = "umap", no.legend = F)
    }else if(input$radio == 3){
      p = FeaturePlot(all, features.plot = ("nUMI"), reduction.use = "umap", no.legend = F)
    }else if(input$radio == 4){
      p = FeaturePlot(all, features.plot = ("percent.mito"), reduction.use = "umap", no.legend = F)
    }else{
      p = DimPlot(all, reduction.use = "umap", do.label = F, do.return = T, pt.size = 2)
    }

    p + theme(axis.title.x = element_text(vjust = 0.5, size = 24, face = "bold"),
              axis.text.x  = element_text(vjust = 0.5, size = 24, face = "bold"),
              axis.title.y = element_text(vjust = 0.5, size = 24, face = "bold"),
              axis.text.y  = element_text(vjust = 0.5, size = 24, face = "bold"),
              legend.text = element_text(size = 20, face = "bold"))
  }, height = 600, width =600)

  output$gene_umap = renderPlot({
    if(input$radio == 1){
      groupInfo = "combined"
    }else{
      groupInfo = "orig.ident"
    }

    FeatureHeatmap(object = all, features.plot = input$text,
                   group.by = groupInfo, reduction.use = "umap",
                   cols.use = c("lightgrey", "darkblue"), pt.size = 2,
                   plot.horiz = F, key.position = "top", do.return = T) +
    theme(axis.title.x = element_text(vjust = 0.5, size = 24),
          axis.text.x  = element_text(vjust = 0.5, size = 24),
          axis.title.y = element_text(vjust = 0.5, size = 24),
          axis.text.y  = element_text(vjust = 0.5, size = 24),
          legend.text  = element_text(size = 24),
          legend.title = element_text(size = 24),
          strip.text.x = element_text(size = 24, face = "bold"),
          strip.text.y = element_text(size = 24, face = "bold.italic"))


  }, height = 600, width = getWidth)

  output$detect_DEGs = renderTable({
    if(input$checkGroup2 == 7){
      group2 = NULL
    }else{
      group2 = input$checkGroup2
    }
    DEGs = FindMarkers(all, ident.1 = input$checkGroup1, ident.2 = group2, logfc.threshold = log(2))
    DEGs$gene = row.names(DEGs)
    print(DEGs)
  })
}

shinyApp(ui = ui, server = server)
