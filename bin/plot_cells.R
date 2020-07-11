
# use parallel processing
library(parallel)

# get passed arguments
env <- as.list(Sys.getenv())

# load data (remember, this overwrites the environment...)
rDataFile <- paste(env$BIN_PREFIX, "layers", "RData", sep=".")
load(rDataFile)

# save object data for cell marking
rDataFile <- paste(env$BIN_PREFIX, "cell_marks", "RData", sep=".")
marks <- rep(-1, cell$N)
save(marks, file=rDataFile)

# parse the chromosome boundaries
maxI <- 0
vLines <- c(1)
mids <- numeric()
chroms <- bin$data[[1]]
uniqueChroms <- unique(chroms)
for (chrom in uniqueChroms){
    nChrom <- length(chroms[chroms==chrom])
    vLines <- c(vLines, nChrom + maxI)
    mids <- c(mids, maxI + nChrom / 2)
    maxI <- maxI + nChrom
}

# calculate the fractional predicted copy number per cell x bin
CNInt <- matrix(rep(bin$modalCNsInt, cell$N), ncol=cell$N)
rpa   <- layers$exp / CNInt # i.e. reads per allele
cn    <- layers$raw / rpa

# make one optimized plot for each cell
CELL_PLOT_DIR <- paste(env$PLOT_DIR, "cells", sep="/")
minBinN <- 100
sink <- mclapply(cell$Is, function(j){
    
    # condense the bins to give a statistically interpretable plot
    binAvg <- cell$sums[j] / bin$N  
    nCollapseBins <- ceiling(minBinN / binAvg) 
    collapsed <- collapseVector(cn[,j], nCollapseBins) / nCollapseBins
    nPoints <- length(collapsed)
    
    # initialize the plot    
    pngFile <- paste(CELL_PLOT_DIR, "/", env$DATA_NAME, ".cell.", sprintf("%05d", j), ".png", sep="") 
    png(pngFile, width=2400, height=1200, units="px", pointsize=6, res=600, type="cairo")
    par(mar=c(4.1, 4.1, 1.1, 1.1))
    plot(0, 0, typ="n",
         ylim=c(-0.5,5.5), ylab="Copy Number", yaxs="i",
         xlim=c(1,nPoints),xlab=paste("Bin at mean(read count)~", minBinN, sep=""), xaxs="i")

    # add demarcation lines and labels
    abline(v=vLines[2:(length(vLines)-1)]/nCollapseBins, col="grey40", lwd=1)
    abline(h=0:5, col="grey40", lwd=1)
    text(mids/nCollapseBins, y = 5.25, labels = gsub("chr", "", uniqueChroms),
         cex = 0.65, col = "darkred")
    
    # add the bin data points for one cell
    points(1:nPoints, collapsed, pch=16, cex=0.5, col="blue")
    dev.off()
    
}, mc.cores = env$N_CPU)

message("plot_cells.R done")

