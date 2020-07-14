
#----------------------------------------------------------------------
# ui.R defines the html page layout and defines user inputs
#----------------------------------------------------------------------

# page layout
ui <- fluidPage(
    
    # some page styles
    style="padding: 10px; font-size: 12px;",
    tags$head(
        tags$style(HTML( # hack fixes stupid default Shiny checkbox group indenting
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
        "
    ))),

    # standard Shiny main page layout
    sidebarLayout(
        
        # top level inputs that control the output
        sidebarPanel(
            
            # sample selection
            selectInput(inputId = 'project',
                        label = 'Project | Sample',
                        choices = projects),
            selectInput(inputId = "sample",
                         label = NULL,
                         choices = c('-')),
            hr(),
            
            # filters for selecting the plotted cells
            checkboxGroupInput(inputId = 'cellTypeFilter',
                               label = "Include Cells Marked As",
                               inline=TRUE,
                               choiceNames = cellTypes[acceptedCellTypes],
                               choiceValues = acceptedCellTypes,
                               selected = acceptedCellTypes),
            #radioButtons(inputId = "containsCNV",
            #             inline = TRUE,
            #             label = "Cell Has HMM CNV",
            #             choices = c('all', 'yes', 'no'),
            #             selected='all'),
            hr(),
            
            # plotting options
            radioButtons(inputId = "colorHeatMapBy",
                label = "Color Heat Map By",
                choices = c("CNC", "CN"), #, "HMM", "calls"
                selected = "CNC",
                inline = TRUE),

            # clustering options
            radioButtons(inputId = "viewportClusterValue",
                label = "Cluster Heat Map By",
                choices = c("CNC", "Z0", "coverage"),
                selected = "CNC",
                inline = TRUE),
            radioButtons(inputId = "distMethod",
                label = "Cluster Distance Method",
                choices = c("euclidean", "manhattan", "pearson"),
                selected = "euclidean",
                inline = TRUE),
            hr(),

            # navigation options
            span("Show All "),
            actionButton("resetBins", "Bins", inline=TRUE),
            actionButton("resetCells", "Cells", inline=TRUE),
            actionButton("resetBinsCells", "Both", inline=TRUE),
            br(),br(),
            hr(),
            
            # a tiny reads per allele plot
            plotOutput(outputId="readPerAllelePlot", height='250px'),
            hr(),
            
            # global action buttons
            actionButton(inputId = "reload", "Refresh Sample Names"),

            # width of the option section
            width=3
        ),
        
        # one row of sample selector inputs for each available plot color
        mainPanel(
            # viewport selectors            
            fluidRow( 
                style="font-weight: bold;",
                column(2, textOutput(outputId = "genome")),
                column(2, "Start"),
                column(2, "End"),
                column(1, "Left"),
                column(1, "Right"),
                column(1, "In"),
                column(1, "Out"),
                column(2, "JumpTo")
            ),
            fluidRow(
                column(2, selectInput(inputId = "chrom",
                                      label = NULL,
                                      choices = c("-") )), # -
                column(2, textInput(inputId = "start",
                                    label = NULL,
                                    value = "" )),
                column(2, textInput(inputId = "end",
                                    label = NULL,
                                    value = "" )),
                column(1, actionButton(inputId = "jumpLeft",  "<")),
                column(1, actionButton(inputId = "jumpRight", ">")),
                column(1, actionButton(inputId = "zoomIn",    "+")),
                column(1, actionButton(inputId = "zoomOut",   "-")),
                column(2, textInput(inputId = "jumpTo",
                                    label = NULL,
                                    value = "" ))
            ),

            # the output plots
            plotOutput(outputId="heatMap", height='auto',
                       brush=brushOpts("heatMapBrush", fill="#9cf", stroke="#036", opacity=0.25,
                                        delay=500, delayType="debounce", resetOnNew=TRUE)
            ),
            plotOutput(outputId="dendogram", click="dendogramClick"),
            plotOutput(outputId="aggregatePlot"),
            plotOutput(outputId="aggregatePlot_uncorrected")
        )   
    )
)

