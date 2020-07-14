
#----------------------------------------------------------------------
# server.R returns a function that defines reactive server actions
#----------------------------------------------------------------------

# handle user actions and update output
`%notin%` <- Negate(`%in%`)
server <- function(input, output, session){
    sessionEnv <- environment()
    verbose <<- FALSE

    # load all samples and supporting data again when requested (e.g. when new data are available)
    observeEvent(input$refreshServer, { initalizeApp() })
    
    # handle sample selection
    initSample <- function(){
        reportProgress('initSample')
        setSampleGenome()
    }
    observeEvent(input$project, {
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
    
    # enumerate the qc plots
    qc_plots <- c(
        'readCounts',
        'readsPerBin',     
        'binSize_Initial_vs_GC',        
        'binSize_GC-Adjusted_vs_Mappability',        
        'binSize_Mappability-Adjusted_vs_Mappability',        
        'binSize_CNV-Adjusted_vs_Mappability',
        'binSize_Initial',
        'binSize_GC-Adjusted',        
        'binSize_Mappability-Adjusted',
        'binSize_CNV-Adjusted',        
        'binAdjustments',
        'modalCopyNumbers'
    )
    nQCPlots <- length(qc_plots)

    # activate the qc plots
    lapply(1:nQCPlots, function(i){
        id <- paste0("qcPlot", i)
        output[[id]] <- renderImage({
            reportProgress('renderImage')
            #if(input$sample == "-" | genomeInfo$name == '-') return(list())
            reportProgress(getQCPlot(input, genomeInfo$name, qc_plots[i]))
            list(
                src = getQCPlot(input, genomeInfo$name, qc_plots[i]),
                height = '100%'
            )
        }, deleteFile = FALSE)
    })
    
    # get the sample's genome, for plot file names
    genomeInfo <- reactiveValues(name='-')
    setSampleGenome <- function(){
        reportProgress('setSampleGenome')
        genomeInfo$name <<- '-'
        if(input$sample != "-"){
            load(getChromsFile(input)) # has genome and chroms objects
            genomeInfo$name <<- genome
        }
    }    
}

