
#----------------------------------------------------------------------
# 'plot.R' has functions used to create visualization images of genome regions
# adapted from R Shiny web tool
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# file names
#----------------------------------------------------------------------
getPlotName <- function(type, bins){
    paste(env$PLOT_PREFIX, 
          chrom, sprintf("%04d", bins[1]), sprintf("%04d", bins[length(bins)]),
          type, "png", sep=".")
}
getPlotName2 <- function(type){
    paste(env$PLOT_PREFIX, 
          chrom, type, "png", sep=".")
}

#----------------------------------------------------------------------
# heat map of CNC per bin+cell combination
#----------------------------------------------------------------------
plotZLayer <- function(bins, cellStack, type, zL=NULL){
    if(is.null(zL)) zL <- zLayers
    saveHeatMap_three_color(
        zToIntensity_neg_decay(-zL$zL[bins,cellStack]),
        zToIntensity_symmetric( zL$z0[bins,cellStack]), 
        zToIntensity_neg_decay( zL$zG[bins,cellStack]),
        getPlotName(paste('heatmap', type, sep="_"), bins),
        mirror=FALSE, xScaleFactor=2, yScaleFactor=1
    )        
}

#----------------------------------------------------------------------
# square heat maps of metrics that define bin-bin relationships
#----------------------------------------------------------------------
addLocalMaxima <- function(cimg, scaleFactor, unused, Z_sd, Z_se){

    # make sure we only look in the lower triangle (where Z_se is plotted)
    Z_se[upper.tri(Z_sd)] <- NA
    Z_se[upper.tri(Z_se)] <- NA    
    diag(Z_sd) <- NA
    diag(Z_se) <- NA
    
    # find local and absolute maxima in Z_se as data frames of points
    lmax   <- findMatrixMaxima(Z_se, threshold=threshold_se)
    lmaxXY <- which(lmax, arr.ind=TRUE)
    max    <- max(Z_se, na.rm=TRUE)
    maxXY  <- which(Z_se == max, arr.ind=TRUE)

    # set all local maxima to grey pixels
    if(nrow(lmaxXY) > 0) {
        for(i in 1:nrow(lmaxXY)){ # this is slow simply because there are so many maxima
            x <- lmaxXY[i,1]
            y <- dim(cimg)[1] - lmaxXY[i,2] + 1            
            cimg[x,y,,] <- 0 # or color by z_se intensity, or rank?
        }
    }

    # set the maximum value to a black pixel
    #if(max >= threshold_se)
    #cimg <- setPixel(cimg, maxXY[1,1], maxXY[1,2], 0)

    # return the image
    cimg
}
addLocalMaxima2 <- function(cimg, scaleFactor, unused, Z_sd, Z_se){

    # set all local maxima to grey pixels
    if(nrow(maximaIs) > 0) {
        #x <- maximaIs[[1]]
        #y <- dim(cimg)[1] - maximaIs[[2]] + 1
        
        x <- maximaIs[[2]]
        y <- dim(cimg)[1] - maximaIs[[1]] + 1

        nx <- length(x)
        for(c in 1:3){
            dim <- matrix(c(x,y,rep(1,nx),rep(c,nx)), nrow=nx, ncol=4)
            cimg[dim] <- 0
        }
    }

    # return the image
    cimg
}


# composite square matrix of bin-to-bin relationships, all cells independently of cell clustering
plotSquareMatrix <- function(sType, bins, pairs, decorate=NULL){ # intensities

    # combine three representations into a single square matrix plot
    nPlotBins <- nrow(pairs$Z)
    upperTri <- pairs$Z
    lowerTri <- pairs$Z_sd
    Z_se     <- pairs$Z_se
    Z <- matrix(0, nPlotBins, nPlotBins)
    lowerTri <- mirrorPyramid(lowerTri) # make sure any lower.tri values are mirrored
    Z[upper.tri(Z)] <- upperTri[upper.tri(upperTri)]
    Z[lower.tri(Z)] <- lowerTri[lower.tri(lowerTri)]

    # parse the yellow highlighting in the lower triangle based on Z_se
    Z_se <- mirrorPyramid(Z_se)
    highlight <- abs(Z_se) 
    highlight[upper.tri(highlight, diag=TRUE)] <- 0

    # create the heat map
    saveHeatMap_two_color(zToIntensity_pos(-Z),
                          zToIntensity_pos( Z),
                          highlight / max(highlight[!is.nan(highlight)]) / 2,
                          getPlotName(sType, bins),
                          lowerTri, Z_se, mirror=FALSE,
                          xScaleFactor=1, yScaleFactor=1,
                          decorate=decorate) 
}

#----------------------------------------------------------------------
# find_cells MA style plot
#----------------------------------------------------------------------
plotFindCells_RLL <- function(bins, cnv){
    png(getPlotName('cell_RLL', bins),
        width = 2400, height = 2400, units = "px", pointsize = 8,
         bg = "white",  res = 600, type = "cairo")
    nCnvCells <- sum(cnv$cell$hasCNV == 1, na.rm=TRUE)
    main <- paste(
        'cnc', cnv$cnc,
        ', bins', cnv$size$bins,
        ', cells', nCnvCells, '/', cell$N_accepted
    )
    cnvColor <- if(cnv$cnc == 1) "red3" else "blue"
    plot(cnv$cell$exp0, cnv$cell$RLL, main=main,
         xlab='Cell Depth', ylab='RLL',
         pch = ifelse(is.na(cnv$cell$hasCNV), 16, # ambiguous cells
               ifelse(cnv$cell$hasCNV==-1,4, # masked cells (i.e. with wider CNV)
               16)), # unmasked cells
         col=ifelse(is.na(cnv$cell$hasCNV), "grey50", # ambiguous cells
             ifelse(cnv$cell$hasCNV!=0, cnvColor, # CNV cells
            'green3' )) # non-CNV cells
    )    
    abline(h=0)
    graphics.off()
}

#----------------------------------------------------------------------
# find_cells MA style plot
#----------------------------------------------------------------------
plotZByCnvSize <- function(type, bins, Zs, maxima, h=0.2){
    png(getPlotName(type, bins),
        width = 2400, height = 2400, units = "px", pointsize = 8,
         bg = "white",  res = 600, type = "cairo")
    #plot(pyramid$N_paired, Zs[[type]], 
    #     xlab='CNV Size (# of bins)', ylab=type,
    #     pch = 16, cex=0.1,
    #     col=ifelse(maxima, "red3", 'blue' )
    #)
    plot(
         pyramid$N_single[maxima],
         #Zs[['Z_se']][maxima], 
         Zs[[type]][maxima], 
         xlab='CNV Size (# of bins)',
         ylab=type, ylim=c(0,max(Zs[[type]][maxima], na.rm=TRUE)),
         pch = 16, cex=0.35,
         col="red3"
    )
    #abline(h=c(0,h))
    graphics.off()
}


#----------------------------------------------------------------------
# linear plot along bins axis
#----------------------------------------------------------------------
plotAlongBins <- function(type, x, y, h=integer()){
    png(getPlotName2(type),
        width = 4800, height = 2400, units = "px", pointsize = 8,
         bg = "white",  res = 600, type = "cairo")
    plot(x,
         y, 
         xlab='Bin',
         ylab=type,
         #ylim=c(0,max(y, na.rm=TRUE)),
         typ="l"
         #,
         #pch = 16, cex=0.35,
         #col="red3"
    )
    abline(h=c(0,h))
    graphics.off()
}

