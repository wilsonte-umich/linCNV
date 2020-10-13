
# use parallel processing
library(parallel)

# initialize script
env <- as.list(Sys.getenv())
source(paste(env$PIPELINE_DIR, "common", "utilities.R", sep="/"))

# load cell labels
manifest <- read.csv(env$MANIFEST_FILE, header=TRUE, stringsAsFactors=FALSE)
manifest <- manifest[manifest$Lane==1,3:4]
rownames(manifest) <- sapply(manifest$Sample_ID, function(sampleId){
    sampleId <- gsub('-','.',sampleId)
    paste('Sample',sampleId,sep="_")
})
manifest$sampleId <- rownames(manifest)
manifest[,c('sampleName','cellNumber')] <- t(as.data.frame(strsplit(unname(unlist(manifest$Description)), '_#')))
manifest$cellNumber <- as.integer(gsub('_','',manifest$cellNumber)) # fixed a typo in one name
manifest <- manifest[order(manifest$sampleName,manifest$cellNumber),]

# load crosstab csv
csv <- paste('crosstab', env$BIN_SIZE, 'csv', sep=".")
csv <- paste(env$OUTPUT_DIR, csv, sep="/")
ct <- read.csv(csv, header=TRUE, stringsAsFactors=FALSE)
env$BIN_SIZE <- as.integer(env$BIN_SIZE)

# set cell indices
cellJs <- 7:(ncol(ct)-3)
cellKs <- cellJs - 6
cellIds <- colnames(ct)[cellJs]
cellNames <- manifest[cellIds,'Description']
nCells <- length(cellIds)
cols <- rainbow(nCells)[order(cellNames)]

# filter and re-order to canonical chromosomes
chroms <- paste0('chr', c(as.character(1:22), 'X')) #, 'Y'
chroms <- chroms[chroms %in% unique(ct$chrom)]
nChroms <- length(chroms)
ct <- do.call('rbind', lapply(chroms, function(chrom) ct[ct$chrom==chrom,]))

# enforce minimum mappability and effective bin size: TODO, expose as options
minMappability <- 0.8
minEffectiveSize <- env$BIN_SIZE / 10
ct$binSize <- env$BIN_SIZE - ct$excluded_bases
ct <- ct[ct$mappability >= minMappability & ct$binSize >= minEffectiveSize,]

# done filtering bins, set additional needed parameters
nBins <- nrow(ct)
autosomes <- ct$chrom %notin% c('chrX', 'chrY')
chromBoundaries <- c(0, cumsum(rle(ct$chrom)$lengths)) + 0.5
ct$binSizeCorrection <- env$BIN_SIZE / ct$binSize
obs <- apply(ct[,cellJs], 2, '*', ct$binSizeCorrection)

# clear prior images
for(sampleName in manifest$sampleName) {
    dir <- paste(env$PLOT_DIR, sampleName, sep="/")
    dir.create(dir, showWarnings=FALSE)
    unlink(paste(dir, '*', sep="/"))
}

# calculate cell statistics
cellSums <- colSums(ct[autosomes,cellIds])
#cellBinMedians <- apply(ct[autosomes,cellIds], 2, median)

# plot cell read counts
plotCellRanks <- function(name, y, ylab){
    file <- paste(name, 'png', sep=".")
    file <- paste(env$PLOT_DIR, file, sep="/")
    png(file, width = 4, height = 4, units = "in", pointsize = 8, res = 600,  type = "cairo")
    order <- order(y)
    plot(cellKs, y[order], xlab="Cell Rank", ylab=ylab, main=name, pch=16, col=cols[order])
    graphics.off()
}
plotCellRanks('cellSums',       cellSums,       'Total Read Pairs')
#plotCellRanks('cellBinMedians', cellBinMedians, 'Bin Median Count')

# make a composite lorenz plot
file <- paste('lorenz', 'png', sep=".")
file <- paste(env$PLOT_DIR, file, sep="/")
png(file, width = 4, height = 4, units = "in", pointsize = 8, res = 600,  type = "cairo")
plot(c(0,1), c(0,1), typ="l", xlab="Fraction of Genome", xlim=c(0,1), ylab="Fraction of Reads", ylim=c(0,1))
for(k in cellKs){
    y <- obs[autosomes,k]
    order <- order(y)    
    x <- ct[autosomes,'binSize'][order]    
    y <- y[order]
    x <- cumsum(as.numeric(x)) / sum(x)
    y <- cumsum(as.numeric(y)) / sum(y)
    lines(x, y, col=cols[k])
}
graphics.off()

# set fit and plot parameters
toFit <- data.frame(
    x  = ct[autosomes,'gc_fraction'],
    x2 = ct[autosomes,'gc_fraction'] ^ 2,
    w  = ct[autosomes,'binSize']
)
toPredict <- data.frame(
    x  = ct[,'gc_fraction'],
    x2 = ct[,'gc_fraction'] ^ 2
)
gcPch <- 16
gcCex <- 0.4

# fit and plot each cell
plotCell <- function(k){
    plotDir <- paste(env$PLOT_DIR, manifest$sampleName[k], sep="/")
    
    # fit a curve to account for GC bias
    obs_k <- obs[,k]
    toFit$y = obs_k[autosomes]
    #fit <- loess(y ~ x, toFit, weights=toFit$w) # more prone to outlier influence
    fit <- lm(y ~ x + x2, toFit, weights=toFit$w)
    
    # use fit to calculate inferred copy number
    toPredict$y <- obs_k
    exp0 <- predict(fit, toPredict)
    rpa  <- exp0 / 2     
    cn   <- obs_k / rpa       
    
    # plot the GC fit
    file <- paste('gc_fit', manifest$Description[k], k, 'png', sep=".")
    file <- paste(plotDir, file, sep="/")
    png(file, width = 4, height = 4, units = "in", pointsize = 8, res = 600,  type = "cairo") 
    plot(toPredict$x, toPredict$y, pch=gcPch, cex=gcCex, col="red3",
         xlab="Fraction GC", ylab="Bin Read Count") # sex chroms
    points(toFit$x, toFit$y, pch=gcPch, cex=gcCex, col="blue") # autosomes
    points(toPredict$x, exp0, pch=gcPch, cex=gcCex, col="green3")
    graphics.off()

    # plot the genome scatter plot
    file <- paste('bins', manifest$Description[k], k, 'png', sep=".")
    file <- paste(plotDir, file, sep="/")
    png(file, width = 8, height = 4, units = "in", pointsize = 8, res = 600,  type = "cairo")
    minY <- 0
    maxY <- 5
    cn <- ifelse(cn<minY, minY, cn)
    cn <- ifelse(cn>maxY, maxY, cn)
    plot(0, 0, typ="n",
         xlab="Bin Number", xlim=c(0,nBins+1),
         ylab="Inferred Copy Number", ylim=c(minY, maxY),
         main=manifest$Description[k])
    abline(h=minY:maxY, col="grey")
    abline(v=chromBoundaries, col="grey")   
    points(1:nBins, cn, pch=16, cex=0.7)
    graphics.off()
}
sink <- mclapply(cellKs, plotCell)   

q('no')

## fit a curve to GC data using loess based on autosomes only (sex chromosomes unreliable modal CN)
## mappability has not proven to be a significant factor
#ct$all_cells_corrected <- ct$all_cells * ct$binSizeCorrection # adjust values in bins tha overlap gap edges
#fit <- loess(all_cells_corrected ~ gc_fraction, ct[autosomes,])
#ct$fitted <- predict(fit, ct)

## plot normalization summaries
#plotCorrelation <- function(name, x, xlab){
#    file <- paste(name, 'png', sep=".")
#    file <- paste(env$PLOT_DIR, file, sep="/")
#    png(file, width = 4, height = 4, units = "in", pointsize = 8, res = 600,  type = "cairo") 
#    plot(x, ct$all_cells_corrected, xlab=xlab, ylab="Bin Count All Cells", main=name, pch=".")
#    if(name == "gc_fraction") points(x, ct$fitted, pch=".", col="blue")
#    graphics.off()
#}
#plotCorrelation('mappability', ct$mappability, 'Mappability')
#plotCorrelation('gc_fraction', ct$gc_fraction, 'Fraction GC')

## calculate and optimize cell fractions relative to all_cells used for GC correction
#allCellsSum <- sum(ct[autosomes,'all_cells']) # again, only use autosomes when estimating relative cell coverage
#cellSums <- colSums(ct[autosomes,cellJs])
#cellFractions <- cellSums / allCellsSum
#layers <- list()
#layers$obs <- apply(ct[,cellJs], 2, '*', ct$binSizeCorrection)

## plot each cell
#plotCell <- function(k){
#    file <- paste(manifest$Description[k], k, 'png', sep=".")
#    file <- paste(env$PLOT_DIR, manifest$sampleName[k], file, sep="/")
#    png(file, width = 8, height = 4, units = "in", pointsize = 8, res = 600,  type = "cairo")
#    minY <- 0
#    maxY <- 5
#    y <- layers$cn[,k]
#    y <- ifelse(y<minY, minY, y)
#    y <- ifelse(y>maxY, maxY, y)
#    plot(0, 0, typ="n",
#         xlab="Bin Number", xlim=c(0,nBins+1),
#         ylab="Inferred Copy Number", ylim=c(minY, maxY),
#         main=manifest$Description[k])
#    abline(h=minY:maxY, col="grey")
#    #abline(h=log2((1:3)/2), col="blue")
#    abline(v=chromBoundaries, col="grey")    
#    points(1:nBins, y, pch=16, cex=0.7)
#    graphics.off()
#}
#sink <- mclapply(cellKs, plotCell)    

#binDensity <- ct[,cellJs] / ct$binSize
#cellByChrom <- list(counts=matrix(sapply(1:nrow(comb), function(i){
#    chrom  <- comb[i,1]
#    cellId <- comb[i,2]
#    sum(ct[ct$chrom==chrom,cellId])
#}), nrow=nChroms))
#cellByChrom$density <- apply(cellByChrom$counts, 2, '/', chromSizes)            
#cellChromMedians <- apply(cellByChrom$density, 2, median)            
#cellByChrom$normalized <- t(apply(cellByChrom$density, 1, '/', cellChromMedians))

