
#----------------------------------------------------------------------
# server.R returns a function that defines reactive server actions
#----------------------------------------------------------------------

# handle user actions and update output
`%notin%` <- Negate(`%in%`)
server <- function(input, output, session){
    sessionEnv <- environment()
    verbose <<- FALSE

    # declare session-specific data resources; only visible within server block
    # see: https://shiny.rstudio.com/articles/scoping.html
    for(script in c(
        "data_sources.R"
    )) source(paste(serverEnv$ACTIONS_PATH, script, sep="/"), local=TRUE)
    
    # load all samples and supporting data again when requested (e.g. when new data are available)
    observeEvent(input$refreshServer, { initalizeApp() })
    
    # clear the current filter selections
    observeEvent(input$resetFilters, {
        #updateSelectInput(session, 'project', selected='-')
        updateTextInput(session, 'project', value='')
        updateSelectInput(session, 'sample',  selected='-')        
        updateRadioButtons(session, 'cellTypeFilter', selected=0)
        #updateCheckboxGroupInput(session, 'targetTypeFilter', selected=defaultTargetPairTypes)        
    })
    
    # handle sample selection
    initSample <- function(){
        reportProgress('initSample')
        loadSample_marks(input)
        marks_$marks <- marks
        ctf <- as.numeric(input$cellTypeFilter)
        newN <- if(ctf == 0){ 1 } else {
            cllN <- 1:length(marks)
            min(cllN[marks==ctf])
        }
        updateTextInput(session, 'cellNumber', value=newN)
    }
    observeEvent(input$project, {
        if(input$project == "-" || input$project == ""){
            updateSelectInput(session, 'sample', choices=c('-'))
        } else {
            crr <- input$sample
            xs  <- samples[[ input$project ]]
            if(!is.null(xs) & length(xs) > 0){
                selected <- if(crr %in% xs) crr else '-'
                updateSelectInput(session, 'sample', choices=c('-', xs), selected=selected)                     
            } else {
                updateSelectInput(session, 'sample', choices=c('-'))
            }
        }
        initSample()
    })
    observeEvent(input$sample, {
        initSample()
    })
    
    # handle a cell type filter change
    observeEvent(input$cellTypeFilter, {
        initSample()
    })

    # handle forward back actions
    getCurrentI <- function(){
        ctf <- as.numeric(input$cellTypeFilter)
        cllN <- 1:length(marks)
        filteredCellN <- if(ctf == 0) cllN else cllN[marks==ctf]            
        currentN <- as.numeric(input$cellNumber)
        currentI <- which(filteredCellN==currentN)
        list(I=currentI, filteredCellN=filteredCellN)
    }
    advanceCellNumber <- function(increment, forceI=NULL){
        reportProgress('advanceCellNumber')
        if(!is.null(marks)){
            current = getCurrentI()
            prevI <- if(is.null(forceI)) current$I else forceI
            newI <- prevI + increment
            if(newI >= 1 & newI <= length(current$filteredCellN)) {
                updateTextInput(session, "cellNumber", value=current$filteredCellN[newI])
            }
        } 
    }    
    observeEvent(input$previousCell, { advanceCellNumber(-1) })
    observeEvent(input$nextCell,     { advanceCellNumber(+1) })
    
    # handle cell mark actions
    lapply(1:nCellTypes, function(i){
        id <- paste0("cellType", i)
        observeEvent(input[[id]], {
            if(!is.null(marks)){
                currentN <- as.numeric(input$cellNumber)                
                prevMark <- marks[currentN]
                if(prevMark == i){ # same mark, nothing to do, just advance
                    advanceCellNumber(+1)
                } else { 
                    current = getCurrentI()
                    marks[currentN] <<- i
                    saveMarks(input)
                    marks_$marks <- marks    
                    if(input$cellTypeFilter == "0") { # unfiltered list, just a simple advance
                        advanceCellNumber(+1)
                    } else if(current$I < length(current$filteredCellN) - 1) { # mark action removed a cell from the filtered list
                        advanceCellNumber(0, current$I)
                    }                    
                }
            }
        })
    })

    # activate the data plot
    output$binsPlot <- renderImage({
        list(
            src = getCellPlot(input),
            height = '100%'
        )
    }, deleteFile = FALSE)
    
    # activate the current cell type display
    output$cellType <- renderUI({
        sampleKey <- crrSample
        current   <- as.numeric(input$cellNumber)
        h4(if(is.null(marks)){
            "Please load a sample"
        } else {
            if(marks[current] >= 1){
                cellTypes[marks[current]]
            } else {
                "Not marked yet"
            }
        })
    })

    # activate the summary table
    marks_ <- reactiveValues(marks=marks)
    getMarksTable <- function(){
        reportProgress('getMarksTable')
        if(is.null(marks_$marks)) NULL else {
            df <- aggregate(marks_$marks, list(marks_$marks), length)
            colnames(df) <- c('CellType', 'Count')
            df$Percent <- round(df$Count / sum(df$Count) * 100, 1)
            df$CellType <- cellTypes[df$CellType]
            df
        }
    }
    output$summaryTable <- DT::renderDataTable(
        getMarksTable(),
        options=list(
            pageLength=nCellTypes,
            paging=FALSE,
            searching=FALSE
        ),
        selection='none',
        rownames=FALSE            
    )
    output$downloadExcel <- downloadHandler(
        filename = function() {
            paste(input$project, input$sample, Sys.Date(), "xlsx", sep=".")
        },
        content = function(file) {
            write_xlsx(getMarksTable(), path=file, col_names=TRUE, format_headers=TRUE)  
        },
        contentType = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
    )
}

