
#----------------------------------------------------------------------
# plot_data.R extracts the current data to display
#----------------------------------------------------------------------

# reactive variables used to control plot dependencies
viewport_ <- reactiveValues(region=NULL, cells=NULL, stack=NULL)
hclust_   <- reactiveValues(model=NULL)
brushSelection_ <- reactiveValues(bins=NULL, cells=NULL)
resetHeatMap <- function(){
    reportProgress('resetHeatMap')
    viewport_$region <- NULL
    viewport_$cells  <- NULL
    viewport_$stack  <- NULL
    hclust_$model    <- NULL
    clearPos(session)
}
clearBrushSelection <- function(){
    reportProgress('clearBrushSelection')
    brushSelection_$bins <- NULL
    brushSelection_$cells <- NULL
}

# return the vector of typical CN for the displayed genome or chrom
getCNInt <- function(){
    reportProgress('getCNInt')
    # whole genome uses condensed/merged bins
    if(input$chrom == 'all') {
        mergedBins$CNInt
    # a single selected chromosome (not yet filtered for position range)    
    } else { 
        bin$modalCNsInt[bin$accepted & bin$data[[chromI]] == input$chrom]
    }    
}

# return the n x p matrix of the data visible on screen; cells not ordered yet
getViewportData <- function(bins, allowedCells, 
                            level=input$viewportClusterValue, noLowerCase=FALSE){
    reportProgress('getViewportData')
    if(!noLowerCase) level <- tolower(level)
    if(level == "coverage") level <- "z0"
    zLayers[[input$chrom]][[level]][bins,allowedCells]
}
    
# cluster cells by the requested attribute
getCellStack <- function(bins, allowedCells, set=FALSE){ # return the ordered numbers of plotted cells
    reportProgress('getCellStack')
    vals <- getViewportData(bins, allowedCells)
    if(input$viewportClusterValue == "coverage"){ # coverage clusters cells by average CNC across window
        vals <- apply(vals, 2, function(v) rep(sum(v), length(v)) )
    }
    dist <- if(input$distMethod == "pearson"){
      pearson.dist(t(vals))
    } else {
      dist(t(vals), method=input$distMethod)
    }
    hc <- hclust(dist, method="ward.D2") # "complete"  
    if(set) hclust_$model <- hc # heat map passes cluster model to dendogram
    allowedCells[hc$order]
}

# allow user to show only a subset of all accepted cells
# NB: non-accepted cells are NOT viewable (they are gone by this point)
applyCellFilters <- function(bins){
    reportProgress('applyCellFilters')
    
    # the indices of _accepted_ cells; only these are present in zLayers
    allCells <- 1:cell$N_accepted 

    # filter for or against cells that have HMM CNVs in viewport
    cnvFilter <- if(input$containsCNV == 'all'){ 
        TRUE       
    } else {        
        hasCNV <- sapply(allCells, function(cellI){
            sum(zLayers[[input$chrom]]$hmm[bins,cellI] != 0) != 0 # only considers CNVs in the viewport
        })
        if(input$containsCNV == 'yes') hasCNV else !hasCNV
    }
    #cnvFilter <- TRUE
    
    # only show cells on a requested list of user mark types
    allowedCellTypes <- as.integer(input$cellTypeFilter)
    nAllowedTypes <- length(allowedCellTypes)  
    cellTypeFilter <- if(nAllowedTypes > 0 & nAllowedTypes < length(acceptedCellTypes)){
        marks_accepted %in% allowedCellTypes
    } else {
        TRUE    
    }
    
    # if brush selected, only show cells in the box
    #print(brushSelection_$cells)
    brushFilter <- if(!is.null(brushSelection_$cells)){
        allCells %in% brushSelection_$cells
    } else {
        TRUE
    }

    # return the result
    allCells[cnvFilter & cellTypeFilter & brushFilter]
}

# return information on displayed genome bins, the x-axis of the heat map plot
getChromBinData <- function(){
    chromBins <- bin$accepted & bin$data[[chromI]] == input$chrom
    bin$data[chromBins,]    
}
getPosBins <- function(){
    reportProgress('getPosBins')
    chromBinData <- getChromBinData()
    minPos <- as.integer(input$start) # chromosomes plot by coordinate, not bin
    maxPos <- as.integer(input$end)
    bins <- chromBinData$start <= maxPos & chromBinData$end >= minPos
    posBinData <- chromBinData[bins,] 
    list(chromBinData=chromBinData,
         minPos=minPos,
         maxPos=maxPos,
         bins=bins,
         posBinData=posBinData)
}
getGenomeBins <- function(){
    reportProgress('getGenomeBins')
    minBin <- as.integer(input$start) # genome view plots by bin, not coordinate
    maxBin <- as.integer(input$end)
    bins <- minBin:maxBin
    list(minBin=minBin,
         maxBin=maxBin,
         bins=bins)
}

