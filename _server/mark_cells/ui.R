
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
        .buttons-label {
            font-size: 14px;
            font-weight: bold;
            text-align: center;
            margin: auto;
        }
        .discarded-button {
            border-color: rgb(200,100,100);
        }
        .accepted-button {
            border-color: rgb(100,160,100);
        }
        .navigation-button {
            border-color: #666;
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
            radioButtons(inputId = "cellTypeFilter",
                         label = "Cell Type",
                         inline = FALSE,
                         choiceNames  = c('All',cellTypes),
                         choiceValues = c(0,    1:nCellTypes)
                         ),

            # current cell #
            hr(),
            textInput("cellNumber", "Cell Number", value=1),         
            
            # some action buttons
            hr(),
            actionButton(inputId = "resetFilters",  "Reset Page"),
            actionButton(inputId = "refreshServer", "Refresh Server"),

            # width of the option section
            width=2
        ),

        # the output elements top to bottom
        mainPanel(
            plotOutput(outputId="binsPlot"), #, height='900px'
            br(),
            hr(),
            fluidRow(
                column(12, uiOutput("cellType"))
            ),
            br(),
            hr(),
            fluidRow(
                column(2, "Discarded", class="buttons-label"),
                lapply(rev(unacceptedButtons), function(i){
                    id <- paste0("cellType", i)
                    column(2, actionButton(id, cellTypes[[i]], class="discarded-button")) 
                }),
                column(2, ""),
                style="display: flex;"
            ),
            br(),
            fluidRow(
                column(2, "Accepted", class="buttons-label"),
                lapply(rev(acceptedButtons), function(i){
                    id <- paste0("cellType", i)
                    column(2, actionButton(id, cellTypes[[i]], class="accepted-button")) 
                }),
                column(2, ""),
                style="display: flex;"
            ),
            br(),
            fluidRow(
                column(3, ""),
                lapply(cellTypeCodes['Unmarked'], function(i){
                    id <- paste0("cellType", i)
                    column(2, actionButton(id, cellTypes[[i]], class="navigation-button")) 
                }),
                column(1, ""),
                column(1, actionButton("previousCell", "  <  ", class="navigation-button")),            
                column(1, actionButton("nextCell",     "  >  ", class="navigation-button")),
                column(4, "")
            ),
            br(),
            hr(),
            fluidRow(
                column(4, downloadButton('downloadExcel', label="Excel")),
                column(4, DT::dataTableOutput("summaryTable")),
            )
        )
    )
)
