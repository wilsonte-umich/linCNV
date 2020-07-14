
#----------------------------------------------------------------------
# make_plots.R executes the data plotting actions
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# side plot = reads per allele histogram for active sample
#----------------------------------------------------------------------
output$readPerAllelePlot <- renderPlot({
    loadSampleChromData('readPerAllelePlot')
    if(!exists('cell')) return(NULL)
    reportProgress('readPerAllelePlot')
    hist(cell$readsPerAllele_working[cell$Is_accepted], breaks=25,
         xlab="Reads Per Allele", ylab="No. of Cells", main=input$sample)
})   

#----------------------------------------------------------------------
# plot level 1 = heat map of cells x bin depth
#----------------------------------------------------------------------
# establish the pixel pseudo-colors using baseline-corrected values
# store in user session until sample or chrom changes
heatMapColors <- list() # large object, lives in sessionEnv
fillHeatMapColors <- function(){
    if(!is.null(heatMapColors[[input$chrom]][[input$colorHeatMapBy]])) return()
    reportProgress('fillHeatMapColors')

    # generate pseudo-colors
    if(input$colorHeatMapBy == "CNC"){ # whether a cell's bin was different that a typical cell; a calculated value
        G <- normZ(zLayers[[input$chrom]]$z0)     # green = CNC0, i.e. neutral
        R <- normZ(zLayers[[input$chrom]]$zG,  1) # red = gain
        B <- normZ(zLayers[[input$chrom]]$zL, -1) # blue = loss
        
    } else if(input$colorHeatMapBy == "CN"){ # the actual predicted copy number, regardless of cell ploidy
        CNInt <- getCNInt()
        exp1 <- matrix(sapply(1:cell$N_accepted, function(cellI){
            zLayers[[input$chrom]]$exp0[,cellI] / CNInt
        }), ncol=cell$N_accepted)
        #exp0 <- exp1 * 0        
        exp2 <- exp1 * 2
        exp3 <- exp1 * 3
        zG <- (zLayers[[input$chrom]]$raw - exp2) / sqrt(exp2) # green is neutral = CN2 (assuming diploid cells)
        zR <- (zLayers[[input$chrom]]$raw - exp3) / sqrt(exp3) # red is gain, so CN3
        zB <- (zLayers[[input$chrom]]$raw - exp1) / sqrt(exp1) # blue is loss, so CN1
        G <- normZ(zG[,])
        R <- normZ(zR[,],  1)
        B <- normZ(zB[,], -1)
        
    } else if(input$colorHeatMapBy == "HMM"){ # the CNC states called for each cell; a classification
        gain <- zLayers[[input$chrom]]$hmm[,] ==  1 # coloring scheme same as CNC, but pure
        loss <- zLayers[[input$chrom]]$hmm[,] == -1
        G <- ifelse(!gain & !loss, 0.5, 0)
        R <- ifelse(gain, 1, 0)
        B <- ifelse(loss, 1, 0)
    }
    
    # return the result
    heatMapColors[[input$chrom]][[input$colorHeatMapBy]] <<- matrix(rgb(R, G, B), ncol=cell$N_accepted) # color matrix includes all _accepted_ bins and cells
}
getHeatMapColors <- function(bins){
    reportProgress('getHeatMapColors')
    
    # get the order set of cells to plot
    allowedCells <- applyCellFilters(bins) # a subset of the indices of _accepted_ cells; total set = 1:cell$N_accepted
    cellStack <- getCellStack(bins, allowedCells, set=TRUE) # a reordering of allowedCells
    viewport_$cells <- allowedCells
    viewport_$stack <- cellStack
    
    # generate pseudo-colors
    fillHeatMapColors()

    # send the final color matrix for current viewport plot, filtered by bin and sorted for cells
    heatMapColors[[input$chrom]][[input$colorHeatMapBy]][bins,cellStack] 
}
plotHeatMap_genome <- function(){
    reportProgress('plotHeatMap_genome')

    # initialize display coordinates as needed
    if(input$start == "" | input$end == "" | prevChrom != input$chrom){
        prevChrom <<- input$chrom
        updateTextInput(session, 'start', value=1)
        updateTextInput(session, 'end',   value=dim(zLayers[[input$chrom]]$raw)[[1]])
        return(0)
    }
    
    # get bin data for restricting genome plot
    region <- getGenomeBins()
    viewport_$region <- region
    xOffset <- region$minBin - 1
    
    # make the bin-based matrix plot; use rectangle to plot each pixel
    colors <- getHeatMapColors(region$bins) # is no position filtering in whole genome view
    nBins  <- nrow(colors)
    nCells <- ncol(colors)   
    par(mar=c(4.1,4.1,0.1,2.1))
    reportProgress('plotting heat map')
    plot(0, 0, typ="n",
         xlim=c(1,nBins+1)+xOffset, ylim=c(1,nCells+1),
         xlab="Variable Width Bin", ylab="Cell",
         xaxs="i", yaxs="i", yaxt="n")
    majorUnit <- if(nCells >= 200) 50 else { if(nCells >= 100) 25 else 10}
    ticks <- seq(0, nCells, by=majorUnit)
    axis(side=2, labels=ticks, at=ticks+0.5)
    axis(side=4, labels=ticks, at=ticks+0.5)
    xleft <- rep(1:nBins, nCells) + xOffset
    ybot  <- as.vector(sapply(1:nCells, function(y) rep(y, nBins)))
    rect(xleft, ybot, xleft + 1, ybot + 1, col=colors, border=NA) # the slowest step!
    
    # demarcate the chromosomes
    abline(v=chromInfo$starts, lwd=1, col="white") 
    mtext(sub("chr", "", chroms), side=3, line=0.5, at=chromInfo$mids)
    nCells
}

# coordinate-based heat map of a chromosome or smaller region
prevChrom <- "NULL"
plotHeatMap_region <- function(){
    reportProgress('plotHeatMap_region')

    # initialize display coordinates as needed
    chromBinData <- getChromBinData()
    if(input$start == "" | input$end == "" | prevChrom != input$chrom){
        prevChrom <<- input$chrom
        if(!jumpStatus$inProgress){
            updateTextInput(session, 'start', value=min(chromBinData$start, na.rm=TRUE))
            updateTextInput(session, 'end',   value=max(chromBinData$end,   na.rm=TRUE))
            return(0)
        }
    }
    jumpStatus$inProgress <<- FALSE
    
    # get bin data for converting bins to coordinates
    region <- getPosBins()
    viewport_$region <- region
    
    # get cell data and scale so gene bar stays proper size
    colors <- getHeatMapColors(region$bins)
    nCells <- ncol(colors)
    geneTrackFig <- 1 - getGeneTrackPixels() / heatMapHeight(nCells)

    # add a gene track to the plot output
    par(fig=c(0,1,geneTrackFig,1), new=TRUE, mar=c(0.1,4.1,0.1,0.1))
    plot(0, 0, typ="n",
         xlim=c(region$minPos,region$maxPos), ylim=c(0,4), ylab="",
         xaxt="n", yaxt="n", xaxs="i", yaxs="i", axes=FALSE)
    vpWidth <- vpWidth(region)
    genes <- genes[[genome]]
    if(!is.null(genes)){
        geneData <- genes[genes$chrom==input$chrom & genes$start<=region$maxPos & genes$end>=region$minPos,]
        if(vpWidth <= 50e6 & nrow(geneData > 0)){
            rect(geneData$start, geneData$ybot, geneData$end, geneData$ybot + 1,
                 col=geneData$color, border=NA)
            geneData$name <- ifelse((geneData$end-geneData$start)>=vpWidth*0.025, geneData$name, "")
            text(geneData$xname, geneData$yname, geneData$name, col=geneData$color) 
        }        
    }

    # make the coordinate-based matrix plot; use rectangle to plot each pixel
    par(fig=c(0,1,0,geneTrackFig), new=TRUE, mar=c(4.1,4.1,0.1,0.1))
    reportProgress('plotting heat map')
    xScalar <- 1e6
    
    plot(0, 0, typ="n",
         xlim=c(region$minPos,region$maxPos) / xScalar, ylim=c(1,nCells+1),
         xlab="Genome Coordinate (Mpb)", ylab="Cell",
         xaxs="i", yaxs="i")

    xleft  <- rep(region$posBinData$start + 1, nCells)
    xright <- rep(region$posBinData$end,       nCells)
    ybot   <- as.vector(sapply(1:nCells, function(y) rep(y, nrow(region$posBinData))))
    rect(region$minPos / xScalar, 1, region$maxPos / xScalar, nCells+1, col="grey", border=NA)
    rect(xleft / xScalar, ybot, xright / xScalar, ybot + 1, col=colors, border=NA)
    nCells
}
heatMap_ <- reactiveValues(nCells=100)
getGeneTrackPixels <- function(){
    if(input$chrom == "all") 100 else 50
}
getPixelsPerCell <- function(nCells){
    if(nCells <= 100){ 3
    } else if(nCells <= 500){ 2
    } else { 1.5 }
}
heatMapHeight <- function(nCells=NULL){ # enable dynamic heat map plot height based on plotted cell number
    if(is.null(nCells)) nCells <- heatMap_$nCells
    nCells * getPixelsPerCell(nCells) + getGeneTrackPixels() + 50
}
output$heatMap  <- renderPlot({
    loadSampleChromData('plotHeatMap')
    heatMap_$nCells <- if(!is.null(zLayers[[input$chrom]])){
        if(input$chrom == 'all') plotHeatMap_genome() else plotHeatMap_region()
    }  else { # dummy plot to give visual feedback during initial plot actions
        plot(0, 0, typ="n",
             xlim=c(1,1000), ylim=c(1,100),
             xlab="Bin", ylab="Cell",
             xaxs="i", yaxs="i")
        rect(1, 1, 1000, 100, col="grey", border=NA)
        100
    }
}, height = heatMapHeight)

#----------------------------------------------------------------------
# plot level 2 = interactive dendogram of the clustered viewport
#----------------------------------------------------------------------
output$dendogram <- renderPlot({
    if(is.null(hclust_$model)) return(NULL)
    reportProgress('dendogram')
    par(mar=c(4.1,4.1,0.1,0.1))
    #printObjectSizes(sessionEnv, minSizeMb=1)
    plot(hclust_$model, labels=FALSE, main="", xlab="Cell")
    axis(1)
})

#----------------------------------------------------------------------
# plot level 3 = aggregated CNC plot for a cell cluster selected on the dendogram
#----------------------------------------------------------------------
getDendogramClick <- function(level){
    if(is.null(input$dendogramClick)) return(NULL)
    if(is.null(hclust_$model)) return(NULL)
    stackI <- round(input$dendogramClick$x, 0) # the row on displayed heat map, from bottom
    groups <- cutree(hclust_$model, h=input$dendogramClick$y)
    group  <- groups[hclust_$model$order[stackI]] # the cluster group the index cell belongs to
    inGroup <- groups==group
    d <- getViewportData(viewport_$region$bins, viewport_$cells, level) # all cells in stack
    gd <- d[,inGroup,drop=FALSE] # only those cells in the selected group
    x <- 1:nrow(gd)    
    list(
        groups=groups,        
        group=group,
        inGroup=inGroup,
        d=d,
        gd=gd,
        x=x
    )
}
makeDendogramPlot <- function(x, y, yLab=""){
    plot(x, y, typ="n", xlab="Bin", ylab=paste("Mean CNC", yLab), ylim=c(-2.5,2.5))
    abline(h=-2:2, col="grey40")
    abline(h=0, lwd=2)
    points(x, y, pch=16, cex=0.5, col="blue")    
}
output$aggregatePlot <- renderPlot({
    dendogram <- getDendogramClick("CNC") # plotted as average CNC
    if(is.null(dendogram)) return(NULL)
    reportProgress('aggregatePlot')
    makeDendogramPlot(dendogram$x, rowMeans(dendogram$gd))
})
#----------------------------------------------------------------------
# plot level 4 = aggregated CNC plot, independent of all prior baseline corrections
#----------------------------------------------------------------------
output$aggregatePlot_uncorrected <- renderPlot({
    dendogram <- getDendogramClick("raw") 
    if(is.null(dendogram)) return(NULL)
    if(input$chrom == "all") return(NULL)
    reportProgress('aggregatePlot_uncorrected')
    
    #return(rawAggregatePlot(bin, cell, viewport_, dendogram, mergedBins))
    
    #binAdj <- viewport_$region$posBinData$sizes / bin$peakSize_autosome
    binAdj <- 1
    readsPerAllele <- cell$readsPerBin[cell$Is_accepted] / ploidy
    groupCells <- viewport_$cells[dendogram$inGroup] # values are indices of zLayers and thus RPA
    nGroupCells <- length(groupCells)
    y <- rowMeans(matrix(unlist(sapply(1:nGroupCells, function(i){
        CN <- dendogram$gd[,i] / readsPerAllele[groupCells[i]] / binAdj
        CN - viewport_$region$posBinData$CNInt
    })), ncol=nGroupCells))
    makeDendogramPlot(dendogram$x, y, "uncorrected")
})

