
#----------------------------------------------------------------
# filter bins and cells more aggresively to remove unreliable data
# calculate Z scores by Poisson distribution for three-state CNC model
# normalize cell baseline for waviness by adjusting read count expectations
#----------------------------------------------------------------

# use parallel processing
library(parallel)

# initialize script
env <- as.list(Sys.getenv())
source(paste(env$PIPELINE_DIR, "common", "utilities.R", sep="/"))

# load data (remember, this overwrites the environment...)
message("loading data into R")
marksFile <- paste(env$REBIN_PREFIX, "cell_marks", "RData", sep=".")
load(marksFile)
rDataFile <- paste(env$REBIN_PREFIX, "layers", "RData", sep=".")
load(rDataFile)

# re-initialize script so env etc. are current
env <- as.list(Sys.getenv())
env$N_CPU <- as.integer(env$N_CPU)
source(paste(env$PIPELINE_DIR, "common", "utilities.R", sep="/"))
#source(paste(env$PIPELINE_DIR, "normalize", "debug.R", sep="/"))
command <- "analyze"

# disregard potentially untrustworthy bins
# good bin have reliable read counts and don't cross large excluded/unmappable regions
bin$accepted <- bin$modalCNs > as.numeric(env$MIN_MODAL_CN) &
                bin$data$mappability > as.numeric(env$MIN_MAPPABILITY) &
                bin$data$excluded < as.numeric(env$MAX_EXCLUDED_BASES) # purges bins that cross large gaps
bin$N_accepted  <- sum(bin$accepted)
bin$Is_accepted <- which(bin$accepted)
for(i in 1:length(layers)) layers[[i]][-bin$Is_accepted,] <- NA # clear unacceptable bins to NA
reportCount(bin$N,          'binN',         "total input bins")
reportCount(bin$N_accepted, 'binNAccepted', "accepted bins used for analysis")
reportCount(cell$N_accepted, 'cellNAccepted', "accepted cells used for analysis")

# get output chromosomes; save genome info for server
genome <- env$GENOME
chromI <- 1
chroms <- unique(bin$data[[chromI]]) # previously sorted in natural order
genomeInfoFile <- paste(env$ANALYZE_PREFIX, "genome_info", "RData", sep=".")
save(genome, chroms, file = genomeInfoFile)

# set some processing and output variables
minExp  <- 0.01 # cannot expect 0 for Poisson, so use a minimal non-zero value
minBinN <- 100  # merge bins until they have at least this many average counts in a typical cell in genome view
collapseFactor <- ceiling(minBinN / median(layers$raw, na.rm=TRUE)) # for whole genome view

# create a layers array (now a list) with bin counts and Z-scores
collapseByChrom <- function(layer, nMergedBins){ # when collapsing genome, don't merge bins on different chromosomes
    do.call(rbind, mclapply(chroms, function(chrom){
        chromBins <- bin$data[[chromI]] == chrom
        vs <- apply(layers[[layer]][bins_accepted & chromBins,cell$Is_accepted], 2, collapseVector, collapseFactor) 
        matrix(vs, ncol=cell$N_accepted)
    }, mc.cores = env$N_CPU))   
}
#fillZLayerValues <- function(zL, CNInt){ # requires only that raw and exp0 have been set
#  CNInt   <- matrix(rep(CNInt, cell$N_accepted), ncol=cell$N_accepted)
#  rpa     <- zL$exp0 / CNInt # i.e. reads per allele  
#  zL$cn   <- zL$raw / rpa  
#  zL$cnc  <- zL$cn - CNInt  
#  zL$expG <- zL$exp0 + rpa
#  zL$expL <- zL$exp0 - rpa
#  zL$expL <- pmax(zL$expL, minExp) # prevent divide by zero
#  zL$z0   <- (zL$raw - zL$exp0) / sqrt(zL$exp0) # Poisson variance = mean
#  zL$zG   <- (zL$raw - zL$expG) / sqrt(zL$expG)
#  zL$zL   <- (zL$raw - zL$expL) / sqrt(zL$expL) 
#  zL
#}
assembleArray <- function(bins, chrom=NULL){

  # output only has good bins and good cells, must account for this in future handling
  bins_accepted <<- bins & bin$accepted
  binIs <- which(bins_accepted)

  # set base counts for each bin x cell combination
  # collapse (i.e. merge) adjacent bins as requested
  if(!is.null(chrom)){ # the view on an individual chromosome, higher resolution
    zChroms <<- bin$data[binIs,chromI]
    nChromBins <- length(zChroms)
    mergedBins_chrom <- which(mergedBins$chrom == chrom) # already subjected to bins_accepted filter
    exp0_correction_chrom <<- matrix(sapply(1:cell$N_accepted, function(cellI){
      uncollapseVector(exp0_correction[mergedBins_chrom,cellI], collapseFactor, nChromBins) # already subjected to cells_accepted filter
    }), ncol=cell$N_accepted) # keep this value in output to use in CNV calling   
    zL <- list()
    zL$raw  <- layers$raw[binIs,cell$Is_accepted]  
    zL$exp0 <- layers$exp[binIs,cell$Is_accepted] * exp0_correction_chrom  # applies the waviness correction 
    CNInt <- bin$modalCNsInt[binIs]  
  } else { # the low resolution view over the whole genome, use to determine waviness corrections
    mergedBins <<- do.call(rbind, mclapply(chroms, function(chrom){
        chromBins <- bin$data[[chromI]] == chrom
        cv <- collapseVector(bin$modalCNsInt[bins_accepted & chromBins], collapseFactor) / collapseFactor
        cs <- rep(chrom, length(cv))
        data.frame(chrom=cs, CNInt=round(cv, 0), stringsAsFactors=FALSE)
    }, mc.cores = env$N_CPU))      
    nMergedBins <- nrow(mergedBins)    
    zL <- list()
    zL$raw  <- collapseByChrom('raw', nMergedBins)
    zL$exp0 <- collapseByChrom('exp', nMergedBins)
    CNInt <- mergedBins$CNInt
  }

  # calculate CNV counts and Z scores for each bin x cell combination
  fillZLayerValues(zL, CNInt)
}

# print an RData object for CNV finding and R Shiny visualization
saveZLayers <- function(chrom){
  objs <- ls(env=globalenv())
  objs <- objs[!(objs %in% "layers")] # remove prior layers object for memory
  rDataFile <- paste(env$ANALYZE_PREFIX, "layers", chrom, "RData", sep=".")
  save(list = objs, file = rDataFile)
}

# calculate Z scores for the condensed genome bins and call large-scale anomalies
message("constructing uncorrected zLayers object")
zLayers <- assembleArray(rep(TRUE, bin$N))
bin$N_collapsed <- dim(zLayers$raw)[1]

# perform uncorrected HMM by chrom to find aneuploidy and large CNVs, to aid in baseline correction
message("determining bin copy number by cell, genome view")
#HMM_set_persistence(1 - as.numeric(env$TRANSITION_PROB) / 10) # an especially resistant HMM for 1st pass ??
HMM_set_persistence(1 - as.numeric(env$TRANSITION_PROB))
zLayers$hmm <- sapply(1:cell$N_accepted, HMMByChrom, mergedBins$chrom)
hmm_cn <- zLayers$hmm + ploidy
    #plotDebug_cellBins(zLayers[,1:4,zLN$cn], "Bin_CN", "Bin", "Copy Number",
    #                   ylim=c(0,4), hmm=hmm_cn[,1:4])

# use the HMM to establish the reference copy number per cell for PCA baseline correction
# this step keeps us from normalizing away large CNVs!
message("performing baseline correction using z-score PCA and quadratic fit")
getPcaZLevels   <- Vectorize(function(hmm) switch(hmm+(nHMMStates-1), 'zL',   'z0',   'zG'))
getPcaExpLevels <- Vectorize(function(hmm) switch(hmm+(nHMMStates-1), 'expL', 'exp0', 'expG'))
pcaLevels_z   <- getPcaZLevels  (zLayers$hmm)
pcaLevels_exp <- getPcaExpLevels(zLayers$hmm)
ij <- expand.grid(b=1:bin$N_collapsed, c=1:cell$N_accepted)
pcaZs   <- matrix(mapply(function(b,c,l) zLayers[[l]][b,c], ij$b, ij$c, pcaLevels_z),   bin$N_collapsed)
pcaExps <- matrix(mapply(function(b,c,l) zLayers[[l]][b,c], ij$b, ij$c, pcaLevels_exp), bin$N_collapsed)
rm(ij)
    #plotDebug_cellBins(pcaZs[,1:4], "Bin_Z", "Bin", "Input Z Score", ylim=c(-4,4))

# temporarily remove rows (bins) where any cell has a called HMM CN=0
# cells with CN0 mess up the PCA baseline correction due to a lack of variance
pcaZs_na <- ifelse(hmm_cn == 0, NA, pcaZs)
pcaZs_na_omit <- na.omit(pcaZs_na)
pcaZs_omitted_bins <- na.action(pcaZs_na_omit)
rm(pcaZs_na)

# perform PCA-based baseline correction using the condensed genome bins
zLayers_corrected <- zLayers
pca <- prcomp(t(pcaZs_na_omit), center=FALSE, scale.=FALSE, rank.=2) # zRefs are already both scaled and centered
    #plotDebug_PCA12(col=marks[cell$Is_accepted])
    #plotDebug_correlation(pca$rotation[,1], pcaZs_na_omit[,1:4], "Rotation_ZScore", "Rotation", "Z Score", ylim=c(-5,5))
    #z_residual <- matrix(sapply(1:4, function(cellI) {
    #    lm(pcaZs_na_omit[,cellI] ~ pca$rotation[,1] + pca$rotation[,2])$residual
    #}), ncol=4)
    #plotDebug_cellBins(z_residual, "Bin_Z_PCA", "Bin", "PCA-Corrected Z Score", ylim=c(-4,4))

# perform an additional level of baseline correction by fitting a quadratic per cell, per chrom
z_residual <- matrix(unlist(mclapply(1:cell$N_accepted, function(cellI){
  if(is.null(pcaZs_omitted_bins)){
    z_working <- lm(pcaZs_na_omit[,cellI] ~ pca$rotation[,1] + pca$rotation[,2])$residual   
  } else {
    z_working <- pcaZs[,cellI] # the original values not adjusted by PCA; those adjustments added next line
    z_working[-pcaZs_omitted_bins] <- lm(pcaZs_na_omit[,cellI] ~ pca$rotation[,1] + pca$rotation[,2])$residual   
  }
  sapply(chroms, function(chrom){
    mergedBins_chrom <- which(mergedBins$chrom == chrom)
    nChromBins <- length(mergedBins_chrom)
    if(nChromBins > 2){
        lm(z_working[mergedBins_chrom] ~ poly(1:nChromBins, 2))$residual
    } else {
        z_working[mergedBins_chrom]
    } 
  })
}, mc.cores = env$N_CPU)), ncol=cell$N_accepted)
    #plotDebug_cellBins(z_residual[,1:4], "Bin_Z_Quadratic", "Bin", "Quadratic-Corrected Z Score", ylim=c(-4,4))
rm(pcaZs, pcaZs_na_omit, pcaZs_omitted_bins)

# back-calculate new exp0 values from corrected z0 values
# NB: we adjust the _expected_ counts for a bin+cell, not the raw/observed values
message("calculating revised read count expectations for cells x bins")
zLayers_corrected$exp0 <- unlist(mclapply(1:cell$N_accepted, function(cellI){
  z_corrected <- z_residual[,cellI]
  optimizeExp <- function(exp){ # exp could be at any HMM CNC
      z <- (rawIJ - exp) / sqrt(exp)
      abs(z - targetZ)
  }  
  sapply(1:bin$N_collapsed, function(binI){ # NB: we adjust state expectations, not raw data
      exp0 <- zLayers$exp0[binI,cellI] # read count we previously expected at CNC=0
      if(hmm_cn[binI,cellI] == 0){
        exp0 
      } else {
        rawIJ <<- zLayers$raw[binI,cellI]
        targetZ <<- z_corrected[binI] # the new/corrected z value at the HMM CN state  
        expX <- pcaExps[binI,cellI] # read count we previously expected at the HMM CN state
        optimize(optimizeExp, c(1/2, 2) * expX)$minimum * exp0 / expX        
      }
  })
}, mc.cores = env$N_CPU))
rm(hmm_cn, z_residual)

# remember the correction ratios to apply later to each chromosome
exp0_correction <- zLayers_corrected$exp0 / zLayers$exp0 # this is collapsed still, must uncollapse by chrom

# swap the zLayers arrays so that zLayers default is the corrected array
#zLayers_uncorrected <- zLayers # keep zLayers_uncorrected for future reference??
zLayers <- zLayers_corrected
rm(zLayers_corrected)

# fill the rest of the corrected zLayers array and save it
message("constructing corrected zLayers object")
zLayers <- fillZLayerValues(zLayers, mergedBins$CNInt)
zLayers$hmm <- NULL
saveZLayers('all') # all means all chromosome, i.e. the low resolution genome view

# apply the baseline normalization to each chromosome and save it individually for scanning
message("applying baseline correction to each chromosome's unmerged bins")
sink <- mclapply(chroms, function(chrom){ 
  zLayers <<- assembleArray(bin$data[[chromI]] == chrom, chrom=chrom)
  saveZLayers(chrom)
  1
}, mc.cores = env$N_CPU)

#============================================

# steps below (HMM and hierachical clustering) no longer performed up front
# they are deferred until after 1st-pass lineage assembly by bin correlation methods

## perform corrected HMM by chrom to find aneuploidy and large CNVs, as definitive output
#message("determining corrected bin copy number by cell, genome view")
#HMM_set_persistence(1 - as.numeric(env$TRANSITION_PROB))
#zLayers$hmm <- sapply(1:cell$N_accepted, HMMByChrom, mergedBins$chrom)
#    #plotDebug_cellBins(zLayers[,1:4,zLN$cn], "Bin_CN_corrected", "Bin", "Copy Number",
#    #                   ylim=c(0,4), hmm=zLayers[,,zLN$hmm] + ploidy)

## perform genome-level hierachical clustering
#message("clustering cells over entire genome (merged bins)")
#hClustGenome <- hclust(pearson.dist(t(zLayers[,,zLN$z0])))
#saveZLayers('all')

## cluster each chromosome; segment individual cells (will work on clustered cell HMM later)
#message("performing HMM and cell clustering on each chromosome, unmerged bins")
#sink <- mclapply(chroms, function(chrom){ 
#  zLayers <<- assembleArray(bin$data[[chromI]] == chrom, chrom=chrom)
#  hClustChrom <<- hclust(pearson.dist(t(zLayers[,,zLN$z0])))
#  zLayers[,,zLN$hmm] <<- sapply(1:cell$N_accepted, HMM_viterbi, TRUE)
#  saveZLayers(chrom)
#  1
#}, mc.cores = env$N_CPU)

#============================================
# insights from Jun Li:
#If you have “data” as an n-p matrix, for read depths over n samples and p bins along the genome.
#Run PCA with: pca<-prcomp(data)
#prcomp(data)$x and prcomp(data)$rotation will give you the eigen-score and eigen-vector, respectively.  If the row-column needs to be transposed, use t(data).
# You can look at the cluster of the n samples in PC1-PC2 by
#plot(pca$x[,1],pca$x[,2])
#Here prcomp(data)$x[,1:3] and prcomp(data)$rotation[,1:3] will give you the first three PCs.
# I did a quick check:
#> data<-matrix(rnorm(5000,0,1),50,100)
#> dim(data)
#[1]  50 100
#> pca<-prcomp(data)
#> dim(pca$x)
#[1] 50 50
# To regress the read depth for sample i, which is a p-element vector, against the 1st pc-vector:
#lm(data[i,]~pca$rotation[,1])
# The residual of data[i,], after correcting by the regression against PC1, is
#lm(data[i,]~pca$rotation[,1])$residual.
# To regress out the top 2 PCs:
#lm(data[i,]~pca$rotation[,1]+pca$rotation{,2})$residual.
# You can explore regressing by the eigen-score, for the i-th bin of the p bins, and see how the result looks.
#lm(data[,i]~pca$x[,1]+pca$x{,2})$residual.
#========================================

