
# some function used to make plots and such during code development

getDebugPngFile <- function(fileType){
  fileName <- paste(env$DATA_NAME, fileType, "png", sep=".")
  paste(env$PLOT_DIR, "/", fileName, sep="")  
}
plotDebug_correlation <- function(x, y, fileType, xlab, ylab, ylim){
  png(getDebugPngFile(fileType), width=4800, height=4800, units="px", pointsize=6, res=600, type="cairo")
  layout(matrix(1:4, 2, 2, byrow=TRUE))
  for(cellI in 1:4) plot(x, y[,cellI], pch=16, cex=0.5, xlab=xlab, ylab=ylab, ylim=ylim)
  dev.off()       
}
plotDebug_cellBins <- function(y, fileType, xlab, ylab, ylim, hmm=NULL){
  png(getDebugPngFile(fileType), width=3600, height=4800, units="px", pointsize=6, res=600, type="cairo")
  nBins  <- nrow(y)
  nCells <- ncol(y)
  layout(matrix(1:nCells, nCells, 1))
  for(cellI in 1:nCells) {
    plot(1:nBins, y[,cellI], pch=16, cex=0.5, xlab=xlab, ylab=ylab, ylim=ylim)
    if(!is.null(hmm)) points(1:nBins, hmm[,cellI], pch=16, cex=1, col="orange")
  }
  dev.off()       
}
plotDebug_PCA12 <- function(corrected=FALSE, col="black", suffix=NULL){
  fileType <- if(corrected) 'PC1_PC2_corrected' else 'PC1_PC2'
  if(!is.null(suffix)) fileType <- paste(fileType, suffix, sep="_")
  png(getDebugPngFile(fileType), width=2400, height=2400, units="px", pointsize=6, res=600, type="cairo")
  plot(pca$x[,1], pca$x[,2], xlab="PC1", ylab="PC2", col=col)
  #plot(pca$x[,1], pca$x[,2], xlab="PC1", ylab="PC2", typ="n")
  #text(pca$x[,1], pca$x[,2], col=col, pch=as.character(1:cell$N_accepted))
  if(length(col) > 1){
    legend('topright', cellTypes[acceptedCellTypes], col=acceptedCellTypes, pch=1)
  }
  dev.off()       
}
plotDebug_TSNE12 <- function(corrected=FALSE, col="black", suffix=NULL){
  fileType <- if(corrected) 'TSNE1_TSNE2_corrected' else 'TSNE1_TSNE2'
  if(!is.null(suffix)) fileType <- paste(fileType, suffix, sep="_")
  png(getDebugPngFile(fileType), width=2400, height=2400, units="px", pointsize=6, res=600, type="cairo")
  plot(rtsne_out$Y, xlab="tSNE-1", ylab="tSNE-2", col=col)
  #plot(rtsne_out$Y, xlab="tSNE-1", ylab="tSNE-2", typ="n")
  #text(jitter(rtsne_out$Y[,1], amount=1), jitter(rtsne_out$Y[,2], amount=1), col=col, pch=as.character(1:cell$N_accepted))
  if(length(col) > 1){
    legend('topright', cellTypes[acceptedCellTypes], col=acceptedCellTypes, pch=1)
  }
  dev.off()       
}


plotDebug_CN_2_times <- function(d){
  maxI <- 0
  vLines <- c(1)
  mids <- numeric()
  uniqueChroms <- unique(mergedBins$chrom)
  for (chrom in uniqueChroms){
      nChrom <- length(mergedBins$chrom[mergedBins$chrom==chrom])
      vLines <- c(vLines, nChrom + maxI)
      mids <- c(mids, maxI + nChrom / 2)
      maxI <- maxI + nChrom
  }
  #sink <- mclapply(1:cell$N_accepted, function(cellI){
  sink <- lapply(1:cell$N_accepted, function(cellI){
      fileName <- paste(env$DATA_NAME, 'CN_2_times', sprintf("%04d", cellI), "png", sep=".")
      pngFile <- paste(env$PLOT_DIR, 'CN_2_times', fileName, sep="/")  
      png(pngFile, width=2400, height=3600, units="px", pointsize=6, res=600, type="cairo")
      layout(matrix(1:2, 2, 1))
      for(i in 1:2){
        par(mar=c(4.1, 4.1, 1.1, 1.1))
        plot(0, 0, typ="n",
             ylim=c(-0.5,5.5), ylab="Copy Number", yaxs="i",
             xlim=c(1,bin$N_collapsed), xlab="Bin", xaxs="i")        
        # add demarcation lines and labels
        abline(v=vLines[2:(length(vLines)-1)], col="grey40", lwd=1)
        abline(h=0:5, col="grey40", lwd=1)
        text(mids, y = 5.25, labels = gsub("chr", "", uniqueChroms),
             cex = 0.65, col = "darkred")
        # add the bin data points for one cell
        points(1:bin$N_collapsed, d[[i]][,cellI,zLN$cn],  pch=16, cex=0.5, col="blue")        
        points(1:bin$N_collapsed, d[[i]][,cellI,zLN$hmm] + 2, pch=16, cex=1, col="orange")
      }
      dev.off()
      1
  #}, mc.cores = env$N_CPU)
  })
}

