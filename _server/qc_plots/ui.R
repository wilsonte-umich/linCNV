
#----------------------------------------------------------------------
# ui.R defines the html page layout and defines user inputs
#----------------------------------------------------------------------

ui <- fluidPage(
    
    # some page styles
    style="padding: 10px; font-size: 12px;",
    tags$head(tags$style(HTML( # hack fixes stupid default Shiny checkbox group indenting
        ".checkbox-inline, .radio-inline  { 
            margin-left: 4px;
            margin-right: 4px;
        }
        .checkbox-inline+.checkbox-inline, .radio-inline+.radio-inline {
            margin-left: 4px;
            margin-right: 4px;
        }
        hr {
            margin-top: 0;
            margin-bottom: 10px;
            border: 0;
            border-top: 1px solid #666;
        }
        #downloadExcel {
            color: rgb(0,0,200);
        }
        h4 {
            text-align: center;
            color: rgb(0,0,200);
        }
        "
    ))),
    
    # standard Shiny main page layout
    sidebarLayout(
        
        # top level inputs that control the output
        sidebarPanel(
            
            # filters for selecting the plotted data
            selectInput(inputId = 'project',
                        label = 'Project',
                        choices = projects),
            selectInput(inputId = "sample",
                         label = "Sample",
                         choices = c('-')),       
            
            # some action buttons
            hr(),
            #actionButton(inputId = "resetFilters",  "Reset Page"),
            actionButton(inputId = "refreshServer", "Refresh Server"),

            # width of the option section
            width=2
        ),

        # the output elements top to bottom
        mainPanel(
            fluidRow(
                column(6, plotOutput(outputId="qcPlot1")),
                column(6, plotOutput(outputId="qcPlot2"))
            ),
            fluidRow(
                column(6, plotOutput(outputId="qcPlot3")),
                column(6, plotOutput(outputId="qcPlot4"))
            ),   
            fluidRow(
                column(6, plotOutput(outputId="qcPlot5")),
                column(6, plotOutput(outputId="qcPlot6"))
            ),
            fluidRow(
                column(6, plotOutput(outputId="qcPlot7")),
                column(6, plotOutput(outputId="qcPlot8"))
            ),
            fluidRow(
                column(6, plotOutput(outputId="qcPlot9")),
                column(6, plotOutput(outputId="qcPlot10"))
            ),   
            fluidRow(
                column(6, plotOutput(outputId="qcPlot11")),
                column(6, plotOutput(outputId="qcPlot12"))
            ) 
        )
    )
)
