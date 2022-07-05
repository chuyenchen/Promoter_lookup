library(shiny)
library(dplyr)
library(tidyverse)
library(stringr)
library(readxl)

#BUG:1)promoterlib50k has oct-1. 2)download not in csv format
#UI
ui <-  fluidPage(tags$h1("Synthetic Promoter Library Lookup", align = "center"),
                 tabsetPanel(
                   tabPanel("Overview",
                            h5("We are a synthetic immunology lab at Dana-Farber Cancer Institute. Aim for developing novel therapies that can overcome key challenges in cancer immunotherapy by the integration of synthetic biology, computational biology and immunology. For more information, please visit:", a(href = 'https://www.syntheticimmunity.net/', "Lab website")),
                            br(),
                            hr(),
                            h5("1. This app allows you to search synthetic promoter sequence(s) of interest, look their counts and export in a csv file."),
                            h5("2. Database", strong("PromoterLibrary6100")," contains about 6100 sequenceses of human transcription factors from position weight matrix (PWM) deposit in the ENCODE and Motif databases (The MEME Suite)."),
                            h5("3. Database", strong("MYC_E2F1_SOX8_synpro"), "contains 110 top consensus sequenceses of human transcription factor MYC, E2F1 and SOX8."),
                            h5("4. Database", strong("PromoterLibrary50000"), "contains about 17K sequenceses which are derived from top 3 PWMs of human transcription factors."),
                            h5("5. References:"),
                            h6(em("a.Wu MR*, Nissim L*, Stupp D*, Pery E, Binder-Nissim A, Weisinger K, Enghuus C, Palacios SR, Humphrey M, Zhang Z, Novoa EM, Kellis M, Weiss R, Rabkin SD, Tabach Y, and Lu TK. A High-Throughput Screening and Computation Platform for Identifying Synthetic Promoters with Enhanced Cell-State Specificity (SPECS). Nature Communications 10(1): 2880, 2019. * denotes equal contributions."), 
                               p(em("b.Wu MR, Jusiak B, and Lu TK. Engineering advanced cancer therapies with synthetic biology. Nature Reviews Cancer 19(4): 187-195, 2019."))),
                            br(),
                            hr(),
                            hr(),
                            h6(p(code("This is a Shinny app")), "Developed by Chu-Yen Chen"),
                            h6("Last updates: 2020-03-07"),
                            tags$img(src="logo_v3.png", height = '20%', width = '20%'),
                            tags$img(src="DFCI_logo.png", height = '20%', width = '20%')),
                   tabPanel("Search",
                            titlePanel(""),
                            sidebarLayout(position = "left",
                                          sidebarPanel(selectInput('db',"Databases", c("PromoterLibrary6100", "MYC_E2F1_SOX8_synpro", "PromoterLibrary17000"), selected = "PromoterLibrary6100"),
                                                       textInput("goi", "Gene of Interest (use '|' to make multiple selections. E.g. E2F1|PAX8)",""),
                                                       actionButton("submit", ("Submit"))
                                          ),
                                          mainPanel(
                                            strong("Selected databse:"), 
                                            br(),
                                            textOutput("db_name"),
                                            hr(),
                                            strong("Total number of genes:"), 
                                            br(),
                                            textOutput("length"),
                                            hr(),
                                            strong("Number of gene of interest:"), textOutput("length_goi"),
                                            br(),
                                            hr(),
                                            strong("Summary of the database:"),
                                            tableOutput("db_all")
                                          )
                            )
                   ),
                   tabPanel("Plots",
                            titlePanel(""),
                            "Number of TFs",
                            numericInput("c", "Counts greater than", value = 20, min = 0, max = 300, step = 10), 
                            plotOutput("TF_count")
                            ),
                   tabPanel("Results",
                            titlePanel(""),
                            mainPanel(
                              downloadButton("download", ("Download")),
                              tableOutput("db_filtered")
                            )
                   )
                 ))


#SERVER
server <- function(input, output) {
  db1 <- eventReactive(input$submit,{
    if(input$db == "PromoterLibrary6100"){
      as.data.frame(read_excel("promoter sequences and names.xlsx"))
    } else if (input$db == "MYC_E2F1_SOX8_synpro"){
      as.data.frame(read_excel("Lib1_final.xlsx"))
    } else if (input$db == "PromoterLibrary50000"){
      as.data.frame(read_excel("allTF_synpro_top10_final.xlsx"))
    }
  })
  
  output$length <- renderText(length(db1()$promoter_name))
  output$db_all <- renderTable(summary(db1()))
  output$db_name <- renderText(input$db)
  output$db_filtered <- renderTable(db1() %>% filter(str_detect(promoter_name, input$goi)))
  output$length_goi <- renderText(sum(str_detect(db1()$promoter_name, input$goi)))
  
  db2 <- reactive(db1() %>% filter(str_detect(promoter_name, input$goi)))
  output$download <- downloadHandler(
    filename = function(){paste("Results",".csv")},
    content = function(fname){write.csv(db2(), fname)}
  )
  
  db3 <- reactive(filter(count(db1(), TF), n>input$c))
  output$TF_count <- renderPlot({
    ggplot(db3(), aes(x=TF, y=n)) + geom_histogram(stat = "identity") + theme(axis.text.x = element_text(angle =90, vjust = 0.5, hjust = 1), plot.title = element_text(hjust = 0.5)) + theme_bw()
  })
  }

shinyApp(ui = ui, server = server)