
#----------------------------------------------------------------------
# server.R returns a function that defines reactive server actions
#----------------------------------------------------------------------

# handle top level user actions
server <- function(input, output, session){
    sessionEnv <- environment()
    verbose <<- FALSE

    # declare session-specific data resources; only visible within server block
    # see: https://shiny.rstudio.com/articles/scoping.html
    loadSessionScripts <- function(){
        for(script in c(
            "data_sources.R",
            "plot_data.R",
            "make_plots.R",
            "interactive.R",
            "../imager.R"
        )) source(paste(serverEnv$SLAVE_PATH, script, sep="/"), local=sessionEnv)
    }
    loadSessionScripts()
    
    # load all samples and supporting data again when requested (e.g. when new data are available)
    observeEvent(input$reload, {
        initalizeApp()
        loadSessionScripts()     
    })

    ## clear the current filter selections
    #observeEvent(input$resetFilters, {
    #    updateSelectInput(session, 'project', selected='-')        
    #    updateSelectInput(session, 'sample',  selected='-')        
    #    updateRadioButtons(session, 'cellTypeFilter', selected=0)
    #    #updateCheckboxGroupInput(session, 'targetTypeFilter', selected=defaultTargetPairTypes)        
    #})    
    
    # handle sample selection
    initSample <- function(){
        reportProgress('initSample')
        setSampleGenome(input, output)
        chroms_ <- if(is.null(genomeInfo$chroms)) c("-") else c('all', genomeInfo$chroms)
        #chroms_ <- if(is.null(genomeInfo$chroms)) c("-") else c('all', genomeInfo$chroms, 'special')
        crr <- input$chrom # keep the chromosome the same when changing samples, if possible
        selected <- if(crr %in% chroms_) crr else 'chr19' #'chr3' 'all' 
        updateSelectInput(session, 'chrom', choices=chroms_, selected=selected)        
        heatMapColors <- list()
        resetHeatMap()
    }
    observeEvent(input$project, {
        reportProgress('input$project')
        if(input$project == "-"){
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

    # handle chrom change
    observeEvent(input$chrom, { resetHeatMap() })
    
    # handle page resets
    observeEvent(input$resetBins, {
        reportProgress('input$resetBins')
        clearPos(session)        
    })
    observeEvent(input$resetCells, {
        reportProgress('input$resetCells')
        clearBrushSelection()
    })
    observeEvent(input$resetBinsCells, {
        reportProgress('input$resetBinsCells')
        clearBrushSelection()
        clearPos(session)
    })
    # handle viewport actions    
    observeEvent(input$jumpLeft,  { jumpWindow(input, session, -1) })
    observeEvent(input$jumpRight, { jumpWindow(input, session,  1) })
    observeEvent(input$zoomIn,    { zoomWindow(input, session,  1) })
    observeEvent(input$zoomOut,   { zoomWindow(input, session, -1) })
    observeEvent(input$jumpTo,    { jumpTo(input, session, genomeInfo$name) })
}

