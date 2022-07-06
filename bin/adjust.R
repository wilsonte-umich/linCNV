
#----------------------------------------------------------------
# construct the initial, unnormalized bin x cell matrix layers
# then adjust bin properties to account for CNVs and ploidy (but not waviness)
#----------------------------------------------------------------

# use parallel processing
library(parallel)

# get passed arguments
env <- as.list(Sys.getenv())
ploidy <- as.integer(env$PLOIDY)
if(is.null(env$INCLUDE_Y)) env$INCLUDE_Y = ""

# load common functions and parameters
command <- "merge"
source(paste(env$PIPELINE_DIR, "common", "utilities.R", sep="/"))

# load bins x cells as data frame
message("loading data into R")
matrixFile <- paste(env$BIN_PREFIX, "counts", "gz", sep=".")
d <- read.table(matrixFile, header=TRUE, sep="\t",
                stringsAsFactors=FALSE, comment.char="")
endI <- 3 # file format is chrom, start, end, cell1, cell2, ...

# load count of excluded bases in bins and their fraction mabbability and GC
exclFile <- paste(env$BIN_PREFIX, "gc_fraction", "gz", sep=".")
x <- read.table(exclFile, header=TRUE, sep="\t",
                stringsAsFactors=FALSE, comment.char="")

# calculate cell metrics
message("calculating cell and bin metrics")
cell        <- list()
cell$dIs    <- (endI+1):ncol(d)
cell$Is     <- cell$dIs - endI
cell$N      <- length(cell$Is)
cell$sums   <- colSums(d[,cell$dIs], na.rm=TRUE)  # total number of reads per cell (yields 0 for unaccepted cells)
cell$Is_accepted  <- which(!is.na(d[1,cell$dIs]))
cell$dIs_accepted <- cell$Is_accepted + endI
cell$N_accepted   <- length(cell$Is_accepted)

# calculate bin metrics
bin        <- list()
bin$N      <- nrow(d)
autosomeIs <- getAutosomeIs(d[[1]])
chrXI <- toupper(d[[1]]) == "CHRX"
chrYI <- toupper(d[[1]]) == "CHRY"
bin$N_autosome <- sum(autosomeIs)
bin$sizes  <- (d$end - d$start) - x$excluded # the number of usable bases that contributed to the bin
bin$medianSize <- median(bin$sizes)
bin$medianSize_autosome <- median(bin$sizes[autosomeIs])
bin$peakSize_autosome <- distributionPeak(bin$sizes[autosomeIs])$x
bin$data <- d[,1:endI]
bin$data$excluded <- x$excluded
bin$data$mappability <- x$mappability
bin$data$gc_fraction <- x$gc_fraction
bin$data$all_cells <- rowSums(d[,cell$dIs], na.rm=TRUE)

# and a few more cell metrics that require bin metrics
cell$sums_autosome <- colSums(d[autosomeIs,cell$dIs], na.rm=TRUE)
cell$allSum <- sum(bin$data$all_cells, na.rm=TRUE) # total number of reads over all cells
cell$avgSum <- cell$allSum / cell$N_accepted
rm(x)

# report various values
reportCount(bin$medianSize,          "medianBinSize",         "median bin size")
reportCount(bin$medianSize_autosome, "medianBinSizeAutosome", "median autosome bin size")
reportCount(bin$peakSize_autosome,   "peakBinSizeAutosome",   "peak autosome bin size")
reportCount(cell$allSum,             "cellAllSum",            "read pairs over all cells")
reportCount(cell$avgSum,             "cellAvgSum", paste("read pairs per each of", commifyInt(cell$N_accepted), "normalized cells", sep=" "))
reportCount(cell$avgSum/bin$N,       "cellBinAvg", paste("read pairs per each of", commifyInt(bin$N), "bins per cell", sep=" "))

# calculate chromosome coverage depths
message("calculating read depth by chromosome")
chromAgg <- merge(
    aggregate(bin$data$all_cells, by=list(chrom=bin$data[[1]]), sum, na.rm=TRUE),
    aggregate(bin$sizes,          by=list(chrom=bin$data[[1]]), sum, na.rm=TRUE),
    by=c("chrom")
)    
colnames(chromAgg) <- c('chrom','count','size')
chromAgg$depth <- round(chromAgg$count / chromAgg$size, 3) # outside of if below, in case we want this later
chromAggAutosomeIs <- getAutosomeIs(chromAgg$chrom)
chromAggXI <- toupper(chromAgg$chrom) == "CHRX"
chromAggYI <- toupper(chromAgg$chrom) == "CHRY"
print(chromAgg) # for log file

# assign expected copy number of autosomes and sex chromosomes
message("making ploidy and sex determinations")
reportCount(ploidy, "ploidy", "expected sample ploidy (autosome copy number per cell)")
ploidy_working <- if(ploidy == 1) 2 else ploidy # helps deal with sex chromosome variable behavior
ploidyRatio <- ploidy / ploidy_working
bin$data$modalCNInt_working <- ploidy_working
meanAutosomeDepth <- weighted.mean(chromAgg$depth[chromAggAutosomeIs],
                                   chromAgg$size [chromAggAutosomeIs])
XOverAutosome <- chromAgg$depth[chromAggXI] / meanAutosomeDepth
XCNInt_working <- max(1, round(XOverAutosome * ploidy_working, 0))
YCNInt_working <- if(env$INCLUDE_Y == ""){
    0 # will be never used
} else {
    YOverAutosome <- chromAgg$depth[chromAggYI] / meanAutosomeDepth
    round(YOverAutosome * ploidy_working, 0)
}
bin$data[chrXI,'modalCNInt_working'] <- XCNInt_working 
bin$data[chrYI,'modalCNInt_working'] <- YCNInt_working
bin$data$modalCNInt <- bin$data$modalCNInt_working * ploidyRatio # one value per chromosome, NOT sample adjusted
chrXWeight <- XCNInt_working * ploidyRatio
chrYWeight <- YCNInt_working * ploidyRatio
sampleSex <- paste0(paste(rep("X", as.integer(XCNInt_working)), collapse=""),
                    paste(rep("Y", as.integer(YCNInt_working)), collapse=""))
reportFraction(chrXWeight, "chrXWeight", "weight of chrX on same scale as autosome ploidy")
reportFraction(chrYWeight, "chrYWeight", "weight of chrY on same scale as autosome ploidy")
reportString(sampleSex, "sampleSex", "inferred sex of sample over all cells")

# apply GC and mappability corrections to bin sizes
# fit using autosomes, then predict for all bins in genome
message("regressing GC and mappability effects from bin sizes")
modalCNRatio <- ploidy / bin$data$modalCNInt # helps adjust for different CN expectations on sex chromosomes
sizeV <- ploidy_working / (ploidy_working + c(-1, -0.5, 0, 0.5, 1)) # for plotting, loss to gain marks
plotBinSizeRegression <- function(x, y, yPrefix, type, xlim, peak, fit=NULL){
    binSizeLabel <- paste(yPrefix, "Bin Size (bp)")
    histFile <- paste0('binSize_', yPrefix)
    plotHistogram(y[y < 3.5 * peak], histFile, binSizeLabel,
                  c(0,  3.5 * peak), sizeV * peak)
    plotCorrelation(x, y, paste0(histFile, '_vs_', type),
                    paste('Fraction', type), binSizeLabel, cex=0.75,
                    xlim=xlim, ylim=c(0, 5*bin$peakSize_autosome), 
                    hLine=sizeV * peak, fit=fit)    
}
binSizeRegression <- function(regressor, y, yPrefix, type, xlim, prob=0.001){
    x <- bin$data[[regressor]]
    g <- factor(round(x * 25, 0))   
    ya <- y[autosomeIs]
    xa <- x[autosomeIs]
    ga <- g[autosomeIs]
    peak <- distributionPeak(ya)$x    
    yqLow  <- aggregate(ya, list(ga), quantile, prob)
    yqHigh <- aggregate(ya, list(ga), quantile, 1-prob)
    trimOutliers <- unlist(sapply(1:length(ya), function(i){
        yqI <- which(yqLow[[1]] == ga[i])
        qLow  <- yqLow [yqI,2]
        qHigh <- yqHigh[yqI,2]
        ya[i] >= qLow & ya[i] <= qHigh & ya[i] < peak * 2
    }))
    #fit <- loess(y ~ x, # NB: in testing, the polynomial fit works better than loess
    #             data = data.frame(y=ya[trimOutliers], x=xa[trimOutliers]))
    fit <- lm(y ~ poly(x,3), data = data.frame(y=ya[trimOutliers], x=xa[trimOutliers]))
    plotBinSizeRegression(x, y, yPrefix, type, xlim, peak, fit)
    residual <- y / modalCNRatio - predict(fit, data.frame(x=x))
    yOut <- (peak + residual) * modalCNRatio
    yOut <- ifelse(is.na(yOut), y, yOut)
    yOut
}
bin$sizes_corrected <- binSizeRegression('gc_fraction', bin$sizes, 'Initial', 'GC', c(0.25,0.75))
bin$sizes_corrected <- binSizeRegression('mappability', bin$sizes_corrected, 'GC-Corrected', 'Mappability', c(0,1))
plotBinSizeRegression(bin$data$mappability, bin$sizes_corrected,
                      'Mappability-Corrected', 'Mappability', c(0,1),
                      distributionPeak(bin$sizes_corrected)$x)

# initialize sample data layers, each a matrix with: row=bin(x|i), column=cell(y|j)
# non-accepted cells persist as NA values
message("initializing data layers")
layers <- list(
    raw = as.matrix(d[,cell$dIs]) # actual read count for each cell-bin
)
rm(d)

# examine the count range for every cell to scale cells to a common basis
# assumes >50% of autosome bins in a cell have the expected modal CN based on ploidy
message("scaling each cell to quantal copy numbers")
minBinN  <- 200  # merge bins until they have at least this many average counts in a cell
maxAllowedCN <- 20 # account for cancer cells with high copy amplifications
searchCN <- 0:maxAllowedCN
optimizeReadN <- function(readsPerAllele){ # sum of a cell's bin distances to their nearest integer copy number
    sum(sapply(collapsed, function(readN) min(abs(readN / readsPerAllele - searchCN))))
}
cell$readsPerAllele_working <- unlist(mclapply(cell$Is, function(j){ # how many reads in a bin corresponds to one allele for a cell
    if(j %in% cell$Is_accepted){
        binAvg <- cell$sums_autosome[j] / bin$N_autosome
        nCollapseBins <- ceiling(minBinN / binAvg) # merge adjacent bins as needed to get the count up
        collapsed <<- collapseVector(layers$raw[autosomeIs,j], nCollapseBins)
        readsPerAlleleGuess <- binAvg * nCollapseBins / ploidy_working
        optimize(optimizeReadN, c(0.5, 1.5) * readsPerAlleleGuess)$minimum / nCollapseBins
    } else {
        NA
    }  
}, mc.cores = env$N_CPU))
cell$readsPerBin <- cell$readsPerAllele_working * ploidy_working

# calculate a 3-point moving average to give more precision to bin adjustments
message("calculating 3-point moving average")
ratio_ma <- matrix(unlist(mclapply(cell$Is_accepted, function(j){
    ratio <- layers$raw[,j] / cell$readsPerBin[j]
    mav <- ma(ratio, 3)
    mav[1] <- mean(ratio[1:2])
    mav[bin$N] <- mean(ratio[(bin$N-1):bin$N])
    mav
}, mc.cores = env$N_CPU)), ncol=cell$N_accepted)

# examine the count range for every bin to account for impact of CNVs in some cells
#   bin$adjustments are <1 when a minority of cells have a copy number gain (and vice versa)
#   because cells with higher CN make bins get smaller, so we expect fewer reads for other cells
message("adjusting bins for outlier cells") # this does not mean "removing bad/noisy cells"!
rpb <- cell$readsPerBin[cell$Is_accepted]
bin$adjustments <- unlist(mclapply(1:bin$N, function(i){
    js <- which(ratio_ma[i,] >= 0.1) # this line prevents modal CN from being zero (even if it is)
    distributionPeak(ratio_ma[i,js], weights=rpb[js])$x
}, mc.cores = env$N_CPU))
rm(rpb, ratio_ma)

# use bin adjustments to estimate bin modal copy numbers
message("assigning modal copy numbers to bins")
bin$sizes_adjusted <- bin$sizes_corrected / bin$adjustments # thus, bins with some high copy cells are effectively larger than their actual size
bin$peakSize_autosome_adjusted <- distributionPeak(bin$sizes_adjusted[autosomeIs])$x
plotBinSizeRegression(bin$data$mappability, bin$sizes_adjusted,
                      'CNV-Adjusted', 'Mappability', c(0,1), bin$peakSize_autosome_adjusted)
bin$modalCNs <- bin$peakSize_autosome_adjusted / bin$sizes_adjusted * ploidy
bin$modalCNsInt <- pmax(1, round(bin$modalCNs, 0)) # most likely integer copy number of the predominant cell state

# use bin adjustments to alter read count expectations per bin to achieve the modal copy number
message("calculating bin expected values per cell")
layers$exp <- sapply(cell$readsPerBin, '*', bin$adjustments)

# plot value histograms and correlations
message("plotting aggregate plots (histograms and correlations)")
plotHistogram(cell$sums, "readCounts", "Cell Total Read Count", breaks=25)
plotHistogram(cell$readsPerBin, "readsPerBin", "Cell Bin Read Count", breaks=25)
plotHistogram(bin$adjustments, "binAdjustments", "Bin Adjustment", vLine=1)
plotHistogram(bin$modalCNs, "modalCopyNumbers", "Bin Modal Copy Number", c(0, 5), vLine=ploidy)

# save object data for later steps
rDataFile <- paste(env$BIN_PREFIX, "layers", "RData", sep=".")
save.image(rDataFile)

