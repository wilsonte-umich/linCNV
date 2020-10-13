
# get environment variables
env <- as.list(Sys.getenv())

# utility functions
`%notin%` <- Negate(`%in%`)

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

# filter and re-order to canonical chromosomes
chroms <- paste0('chr', c(as.character(1:22), 'X')) #, 'Y'
chroms <- chroms[chroms %in% unique(ct$chrom)]
nChroms <- length(chroms)
ct <- do.call('rbind', lapply(chroms, function(chrom) ct[ct$chrom==chrom,]))

# calculate effective bin sizes
ct$binSize <- env$BIN_SIZE - ct$excluded_bases

# enforce minimum mappability and effective bin size: TODO, expose as options
minMappability <- 0.8
minEffectiveSize <- env$BIN_SIZE / 10
ct <- ct[ct$mappability >= minMappability & ct$binSize >= minEffectiveSize,]

# done filtering bins, set various needed parameters
nBins <- nrow(ct)
autosomes <- ct$chrom %notin% c('chrX', 'chrY')
chromBoundaries <- c(0, cumsum(rle(ct$chrom)$lengths)) + 0.5

# fit a curve to GC data using loess based on autosomes only (sex chromosomes unreliable modal CN)
# mappability has not proven to be a significant factor
ct$all_cells_corrected <- ct$all_cells * (env$BIN_SIZE / ct$binSize) # adjust values in bins tha overlap gap edges
fit <- loess(all_cells_corrected ~ gc_fraction, ct[autosomes,])
ct$fitted <- predict(fit, ct)

# plot normalization summaries
plotCorrelation <- function(name, x, xlab){
    file <- paste(name, 'png', sep=".")
    file <- paste(env$PLOT_DIR, file, sep="/")
    png(file, width = 4, height = 4, units = "in", pointsize = 8, res = 600,  type = "cairo") 
    plot(x, ct$all_cells_corrected, xlab=xlab, ylab="Bin Count All Cells", main=name, pch=".")
    if(name == "gc_fraction") points(x, ct$fitted, pch=".", col="blue")
    graphics.off()
}
plotCorrelation('mappability', ct$mappability, 'Mappability')
plotCorrelation('gc_fraction', ct$gc_fraction, 'Fraction GC')

# apply GC correction to individual cell bin densities
allCellsSum <- sum(ct[autosomes,'all_cells']) # again, only use autosomes when estimating relative cell coverage
cellSums <- colSums(ct[autosomes,cellJs])
cellFractions <- cellSums / allCellsSum
layerNames <- c('corrected','expected','lrr')
layers <- array(0.0, dim=c(nBins,nCells,length(layerNames)), dimnames=list(NULL,NULL,layerNames))
layers[,,'corrected'] <- apply(ct[,cellJs], 2, '/', (env$BIN_SIZE / ct$binSize))
layers[,,'expected']  <- sapply(cellFractions, '*', ct$fitted)
layers[,,'lrr'] <- log2(layers[,,'corrected'] / layers[,,'expected'])

# plot each cell
plotCell <- function(k){
    file <- paste(manifest$Description[k], k, 'png', sep=".")
    file <- paste(env$PLOT_DIR, manifest$sampleName[k], file, sep="/")
    png(file, width = 4, height = 4, units = "in", pointsize = 8, res = 600,  type = "cairo")
    minY <- -3
    maxY <- 3
    y <- layers[,k,'lrr']
    y <- ifelse(y<minY, minY, y)
    y <- ifelse(y>maxY, maxY, y)
    plot(0, 0, typ="n",
         xlab="Bin Number", xlim=c(0,nBins+1),
         ylab="Log2 (observed / expected)", ylim=c(minY, maxY),
         main=manifest$Description[k])
    abline(h=minY:maxY, col="grey")
    abline(v=chromBoundaries, col="grey")    
    points(1:nBins, y, pch=16, cex=0.5)
    graphics.off()
}
for(sampleName in manifest$sampleName) dir.create(paste(env$PLOT_DIR,sampleName,sep="/"), showWarnings=FALSE)
sink <- sapply(cellKs, plotCell)

q('no')











as.matrix(sapply(cellIds, function(cellId){
    ct[,cellId] * (env$BIN_SIZE / ct$binSize)
    expected  <- ct[,'fitted'] * cellFractions[k]
    
    expectedDensity <- 
    binDensity[,k] - 
}), nrow=nBins)

binCorrecctions <- as.matrix(sapply(cellIds, function(cellId){
    corrected <- ct[,cellId] * (env$BIN_SIZE / ct$binSize)
    expected  <- ct[,'fitted'] * cellFractions[k]
    
    expectedDensity <- 
    binDensity[,k] - 
}), nrow=nBins)

q('no')


binDensity <- ct[,cellJs] / ct$binSize



binCorrecctions <- as.matrix(sapply(cellKs, function(k){
    expectedDensity <- ct[,'fitted'] * cellFractions[k]
    binDensity[,k] - 
}), nrow=nBins)


print(str(binResiduals))


# plot each cell
plotCell <- function(k){
    file <- paste(manifest$Description[k], k, 'png', sep=".")
    file <- paste(env$PLOT_DIR, manifest$sampleName[k], file, sep="/")
    png(file, width = 4, height = 4, units = "in", pointsize = 8, res = 600,  type = "cairo") 
    plot(1:nBins, binResiduals[,k], xlab="Bin Number", ylab="Residual Density", pch=16)
    graphics.off()
}
for(sampleName in manifest$sampleName) dir.create(paste(env$PLOT_DIR,sampleName,sep="/"), showWarnings=FALSE)
sink <- sapply(cellKs, plotCell)


#debugCol <- c('chrom','start','end','binN','all_cells','excluded_bases','mappability','gc_fraction',
#              'binSize','binDensity')
#
#print(ct[1:10,debugCol])
q('no')



# calculate cell statistics
cellSums <- colSums(ct[autosomes,cellIds])
cellBinMedians <- apply(ct[autosomes,cellIds], 2, median)
comb <- expand.grid(chroms,cellIds,stringsAsFactors=FALSE)


cellByChrom <- list(counts=matrix(sapply(1:nrow(comb), function(i){
    chrom  <- comb[i,1]
    cellId <- comb[i,2]
    sum(ct[ct$chrom==chrom,cellId])
}), nrow=nChroms))
cellByChrom$density <- apply(cellByChrom$counts, 2, '/', chromSizes)            
cellChromMedians <- apply(cellByChrom$density, 2, median)            
cellByChrom$normalized <- t(apply(cellByChrom$density, 1, '/', cellChromMedians))





#chromSizes <- sapply(chroms, function(chrom) max(ct[ct$chrom==chrom,'end']))

 #$ fitted   : num [1:2767] 0.241 0.233 0.245 0.256 0.254 ...
 #$ residuals: Named num [1:2767] -0.0437 0.0241 0.0603 0.0633 0.0255 ...
 # ..- attr(*, "names")= chr [1:2767] "2" "3" "4" "5" ...

#   chrom   start      end binN all_cells excluded_bases mappability gc_fraction
#1   chr1       0  1000000    1     55329         792500   0.7535743   0.5315096
#2   chr1 1000000  2000000    2    196862              0   0.8567237   0.5727892
#3   chr1 2000000  3000000    3    245694          43509   0.9084559   0.5836914
#4   chr1 3000000  4000000    4    305822              0   0.9703750   0.5646823
#5   chr1 4000000  5000000    5    319340              0   0.9619000   0.4790939
#6   chr1 5000000  6000000    6    279092              0   0.9664240   0.4723347
#7   chr1 6000000  7000000    7    249410              0   0.9412420   0.5163806
#8   chr1 7000000  8000000    8    272522              0   0.9529970   0.4778302
#9   chr1 8000000  9000000    9    187138              0   0.9304540   0.4454565
#10  chr1 9000000 10000000   10    197836              0   0.9231368   0.5014825









                

# per-cell plot function
plotCellRanks <- function(name, y, ylab){
    file <- paste(name, 'png', sep=".")
    file <- paste('plots', file, sep="/")
    png(file, width = 4, height = 4, units = "in", pointsize = 8, res = 600,  type = "cairo") 
    plot(cellKs, sort(y), xlab="Cell Rank", ylab=ylab, main=name, pch=16)
    graphics.off()
}
plotCellRanks('cellSums',       cellSums,       'Total Read Pairs')
plotCellRanks('cellBinMedians', cellBinMedians, 'Bin Median Count')

# per chrom plot function
plotChroms <- function(name, ylab, ylim=NULL, log=FALSE){
    file <- paste(name, 'png', sep=".")
    file <- paste('plots', file, sep="/")
    png(file, width = 4, height = 4, units = "in", pointsize = 8, res = 600,  type = "cairo") 
    #plot(cellKs, y, xlab="Cell Rank", ylab=ylab, main=name, pch=16)
    y <- cellByChrom[[name]]
    if(log) y <- log2(y)
    if(is.null(ylim)) ylim <- range(y)
    y <- ifelse(y < ylim[1], ylim[1], y)
    y <- ifelse(y > ylim[2], ylim[2], y)    
    plot(0, 0, xlim=c(1,nChroms), ylim=ylim, xlab="Chrom # (23=X)", ylab=ylab)
    for(k in cellKs) points(1:nChroms, y[,k], pch=16)
    graphics.off() 
}
plotChroms('counts', 'Counts')
plotChroms('normalized', 'log2 Normalized Density', c(-2,2), log=TRUE)

# plot each cell
plotCell <- function(k){
    file <- paste(manifest$Description[k], k, 'png', sep=".")
    file <- paste('plots', manifest$sampleName[k], file, sep="/")
    png(file, width = 4, height = 4, units = "in", pointsize = 8, res = 600,  type = "cairo") 
    plot(1:nChroms, cellByChrom$density[,k], xlab="Chrom # (23=X)", ylab="Read Density", pch=16)
    graphics.off()
}
for(sampleName in manifest$sampleName) dir.create(paste('plots',sampleName,sep="/"), showWarnings=FALSE)
sink <- sapply(cellKs, plotCell)

#print(cellByChrom$normalized)
q('no')
