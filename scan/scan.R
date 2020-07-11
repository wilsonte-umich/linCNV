
#----------------------------------------------------------------
# 'scan.R' analyzes a single chromosome for CNVs in a detectable fraction of cells
# using methods that are NOT sensitive to prior hierarchical cell clustering
# TODO: could move scanning to immediately follow normalize (while zLayers in memory)
#----------------------------------------------------------------

#----------------------------------------------------------------
# initialize script
#----------------------------------------------------------------
# use parallel processing
library(parallel)

# load environment
env <- as.list(Sys.getenv())
source(paste(env$PIPELINE_DIR, "common", "utilities.R", sep="/"))
source(paste(env$PIPELINE_DIR, "scan", "plot.R", sep="/"))
source(paste(env$PIPELINE_DIR, "scan", "imager.R", sep="/"))

# set CNV-calling thresholds (TODO: make settable?)
# rigorous to avoid false positive CNV calls, but with inevitable ambiguity in cell assignments
nullDistributionPairs <- 1e4 # use at least this many bin pairs to establish the null distribution
threshold_se_maxima <- 4 # initial Z_se threshold used when searching for maxima (more lenient)
threshold_se <- 5.5  # Z_se at maxima must be >threshold_se ...
threshold_sx <- 3.5
threshold_p  <- 1e-8 # ... OR p-value of T-test applied to bin correlation Z-scores must be <threshold_p
cell_p_threshold <- 1/100 # i.e. the fraction of mislabeled cells we will tolerate 
nSamples_ecdf <- 1 / cell_p_threshold * 10 # how hard we try to establish cell quantile values

# set variables that coordinate the process of scanning at increasing merge levels (i.e. small to large CNVs)
maxCnvBinsPerLevel <- 100
foldChangePerLevel <- maxCnvBinsPerLevel / 10

# set variable used to assess CNV status on the flanks of a candidate CNV
nFlankBins <- 15 # empirically determined to give good CNV/non-CNV separation without being TOO large
minFlankBins <- 10 # allow less than nFlankBins, but not too much less

# load genome information for looping chromosomes
genomeInfoFile <- paste(env$ANALYZE_PREFIX, "genome_info", "RData", sep=".")
load(genomeInfoFile)

#----------------------------------------------------------------------
# use called cells with the same CNV to update the CNV boundaries by HMM applied to summed cells
#----------------------------------------------------------------------
refineCnvBoundaries <- function(cnvBins, cnc, cnvCells){ # called at most once per window

    # sum each bin in chromosome over all cells believed to have the CNV
    zL <- list()
    for(layer in c('raw','expL','exp0','expG')){ # the values needed for HMM
        zL[[layer]] <- rowSums(zLayers[[layer]][,cnvCells])  
    }  
    
    # set the emission probabilities and run HMM on the merged cells
    maxCount <- round(2 * zL$expG, 0)
    ep <- sapply(c('expL','exp0','expG'), function(exp){
      log(dpois(pmin(maxCount, round(zL$raw,0)), zL[[exp]]))
    })
    hmm <- HMM_viterbi_ep(ep)
    
    # extract the runs of HMM CNC that match our target
    cncMatch <- hmm == cnc
    if(!any(cncMatch)) return(cnvBins) # unexpected, but just in case HMM didn't call our CNV
    runs <- rle(cncMatch)   
    ends <- cumsum(runs$lengths) # endpoints of the HMM segments    
    starts <- c(1,(ends + 1)[1:length(ends)-1])       
    runsMatch <- runs$values # whether or not the bin run has a matching CNV
    crrStart <- min(cnvBins) # end points of the CNV call that we are refining
    crrEnd   <- max(cnvBins)
    runsOverlap <- crrStart <= ends & starts <= crrEnd # the runs that overlap our target CNV call
    matchingRuns <- runsMatch & runsOverlap
    if(!any(matchingRuns)) return(cnvBins)

    # return the widest span of all matching HMM segments that overlap our called CNV
    # technically, HMM could drop out of our CNV in the middle, but we already know it is there!
    min(starts[matchingRuns]):max(ends[matchingRuns])
}

#----------------------------------------------------------------------
# functions that update the zLayers/zL objects as scanning proceeds
#----------------------------------------------------------------------
# set new read count expectations based on CNV calls
updateBinExpectations <- function(revisedBins, cnv){
    nullCells <- which(cnv$cell$hasCNV == 0)  # based on the read counts observed in NON-CNV cells    
    ratio <- t(apply(zLayers$raw[revisedBins,nullCells], 1, '/', cell_RPB[nullCells]))
    ratio <- ratio / exp0_correction_chrom[revisedBins,nullCells]
    binAdj <- apply(ratio, 1, median, na.rm=TRUE)     
    zLayers$exp0[revisedBins,] <<- sapply(cell_RPB, '*', binAdj) * # accounts for presence of CNV cells
                                   exp0_correction_chrom[revisedBins,] # plus prior baseline drift corrections
    zLayers <<- fillZLayerValues(zLayers, bin_modalCN)      
}
# mask CNVs we have already called so that they don't interfere with finding next CNVs
maskCalledCnv <- function(zL, cnv, cnvBinIs, bmcn){
    rawInc <- cnv$cnc * zL$exp0[cnvBinIs,cnv$cell$Is] / bmcn[cnvBinIs]    
    zL$raw[cnvBinIs,cnv$cell$Is] <- zL$raw[cnvBinIs,cnv$cell$Is] - rawInc
    fillZLayerValues(zL, bmcn)
}

#----------------------------------------------------------------------
# worker functions common to assignCellsToCnv and maskCellsWithWiderCnv
#----------------------------------------------------------------------
# establish random distributions for RLL values based on read depth
# TODO: could we calculate this from statistical principles (i.e. not use random sampling)?
getCellEcdf <- function(exp, expX, exp0){ # exp is the state we are modeling, expX is the CNV state we are testing against
    raw <- sapply(exp, function(x) rpois(nSamples_ecdf, x)) # yields cells in columns
    RLLX <- apply(raw, 1, function(r){ # yields cells in rows
        dpois(r, expX, log=TRUE) - dpois(r, exp0, log=TRUE) # 'exp|x' means 'expected' counts
    })                                                      # 'raw|r' mean 'raw', i.e. observed, counts
    apply(RLLX, 1, ecdf) # return one ecdf per cell
}
# this function actually makes the final cell assignments
getCellsWithCnv <- function(is0, isX, pXX, p0X, RLLX){ 
    ifelse(is0 & isX, NA, # ambiguous state consistent with null AND CNV = no call
        ifelse(is0, 0, # CNV is absent
        ifelse(isX, 1, # CNV is present (1 is a boolean, not a CNC value)
        ifelse(pXX / p0X > 5 & RLLX > 0,  1, # high cell depth, be more lenient for some outliers            
        ifelse(p0X / pXX > 5, 0,
        NA # high cell counts too far away from any one state; ambiguous
    )))))
}

#----------------------------------------------------------------------
# decide whether flanking bins in a cell have the same CNV state we are trying to call
# they shouldn't, we want bounded CNVs only
#----------------------------------------------------------------------
maskCellsWithWiderCnv <- function(cnvBins, cnc, hasCNV, side){
    
    # parse the query bins flanking the CNV span
    if(side == 'left'){
        cnvStart <- min(cnvBins) # relative to chrom, i.e. zLayers
        flankBins <- if(cnvStart > 1) max(1, cnvStart - nFlankBins):(cnvStart-1) else integer()  
    } else {
        cnvEnd <- max(cnvBins)
        nBins <- nrow(zLayers$raw)
        flankBins <- if(cnvEnd < nBins) (cnvEnd+1):min(nBins, cnvEnd + nFlankBins) else integer()
    }
    if(nFlankBins < minFlankBins) return(hasCNV) # not enough flanking bins to care about...    
    
    # calculate Relative Log Likelihood that each cell has in the flanks for the called CNV type
    cnvCells <- which(hasCNV == 1)
    nCnvCells <- length(cnvCells)
    if(nCnvCells < 2) return(hasCNV) # not likely, but don't mask our only CNV cell
    expX <- if(cnc==1) "expG" else "expL" # 'X' refers to the CNV state, either Gain or Loss here
    raw  <- round(colSums(zLayers$raw    [flankBins,cnvCells]), 0) # sum observed and expected over candidate region
    expX <-  pmax(colSums(zLayers[[expX]][flankBins,cnvCells]), minExp)
    exp0 <-  pmax(colSums(zLayers$exp0   [flankBins,cnvCells]), minExp)
    RLLX <- dpois(raw, expX, log=TRUE) - dpois(raw, exp0, log=TRUE)    
    
    # for each cell, establish RLL distributions for null and CNV states based on coverage depth
    ecdfXX <- getCellEcdf(expX, expX, exp0)
    ecdf0X <- getCellEcdf(exp0, expX, exp0)    
    
    # use ECDF to calculate the probabilities of the actual RLL values
    pXX <- sapply(1:nCnvCells, function(i)     ecdfXX[[i]](RLLX[i]) ) # one probability for every cell 
    p0X <- sapply(1:nCnvCells, function(i) 1 - ecdf0X[[i]](RLLX[i]) )

    # use probablities to classify cells
    isX <- pXX > cell_p_threshold
    is0 <- p0X > cell_p_threshold
    flankHasCNV <- getCellsWithCnv(is0, isX, pXX, p0X, RLLX)

    # if flank is consistent with the same CNV type as CNV span, mask cell to value -1
    hasCNV[cnvCells] <- ifelse(flankHasCNV == 1, -1, hasCNV[cnvCells])
    hasCNV
}

#----------------------------------------------------------------------
# assign the states of all cells relative to a given CNV span
# use RLL distributions to characterize each cell's CNV likelihood
# recognize three states: has CNV, lacks CNV, and ambiguous "in between" state
# preserve state probabilities to aid in later lineage conflict resolution
#----------------------------------------------------------------------
assignCellsToCnv <- function(cnvBins){ # potentially called twice per window

    # calculate Relative Log Likelihood that each cell has each type of candidate CNV
    raw  <- round(colSums(zLayers$raw [cnvBins,]), 0) # sum observed and expected over candidate region
    expG <-  pmax(colSums(zLayers$expG[cnvBins,]), minExp)
    exp0 <-  pmax(colSums(zLayers$exp0[cnvBins,]), minExp)
    expL <-  pmax(colSums(zLayers$expL[cnvBins,]), minExp)
    RLLG <- dpois(raw, expG, log=TRUE) - dpois(raw, exp0, log=TRUE) # positive if CNV is more likely than null
    RLLL <- dpois(raw, expL, log=TRUE) - dpois(raw, exp0, log=TRUE)

    # for each cell, establish RLL distributions for null and CNV states based on coverage depth
    ecdfGG <- getCellEcdf(expG, expG, exp0)
    ecdf0G <- getCellEcdf(exp0, expG, exp0) # model null hypothesis state (0), calculate RLL against gain CNV (G) expectations
    ecdfLL <- getCellEcdf(expL, expL, exp0) # etc.
    ecdf0L <- getCellEcdf(exp0, expL, exp0)
    
    # use ECDF to calculate the probabilities of the actual RLL values
    pGG <- sapply(1:cell$N_accepted, function(i)     ecdfGG[[i]](RLLG[i]) ) # one probability for every cell 
    p0G <- sapply(1:cell$N_accepted, function(i) 1 - ecdf0G[[i]](RLLG[i]) ) # '1 - ' yields upper tail
    pLL <- sapply(1:cell$N_accepted, function(i)     ecdfLL[[i]](RLLL[i]) ) # as is yields lower tail
    p0L <- sapply(1:cell$N_accepted, function(i) 1 - ecdf0L[[i]](RLLL[i]) ) 

    # use probablities to classify cells
    isG <- pGG > cell_p_threshold # i.e. cell is NOT significantly different than the hypothesis
    is0 <- p0G > cell_p_threshold & p0L > cell_p_threshold
    isL <- pLL > cell_p_threshold
    calls <- # make an initial cell CNC assignment based on comparison of null, gain AND loss at once
        ifelse(is0 & (isG | isL), NA, # low read depth, RLL consistent with both null a CNV state; ambiguous
        ifelse(is0,  0, # RLL only consistent with neutral state ... 
        ifelse(isG,  1, # ... or a CNV state
        ifelse(isL, -1, # values here are CNCs
        ifelse(pGG / p0G > 5 & RLLG > 0,  1, # high cell depth, be more lenient for some outliers
        ifelse(pLL / p0L > 5 & RLLL > 0, -1,               
        ifelse(p0G / pGG > 5 | p0L / pLL > 5, 0,
        NA # high cell counts too far away from any one state; ambiguous
    )))))))
    
    # determine the predominant CNV type called; presume that it drove the window CNV call
    nCellsGain <- sum(calls ==  1, na.rm=TRUE) # excludes ambiguous cells
    nCellsLoss <- sum(calls == -1, na.rm=TRUE)
    if(nCellsGain == 0 & nCellsLoss == 0) return(NULL) # unusual, but no confidence in even one cell

    # return the result
    fillHasCnv <- function(cnc, isX, pXX, p0X, RLLX){ # a boolean flag for presence/absence of output CNV
        hasCNV <- getCellsWithCnv(is0, isX, pXX, p0X, RLLX) # a vector of booleans
        hasCNV <- maskCellsWithWiderCnv(cnvBins, cnc, hasCNV, 'left')
        hasCNV <- maskCellsWithWiderCnv(cnvBins, cnc, hasCNV, 'right')
        # TODO: do we also need to maskCellsWithNarrowerCnv ?
        N <- sum(hasCNV == 1, na.rm=TRUE)
        list(
            cnc = cnc,
            cell = list(
                cnc = colMeans(zLayers$cnc[cnvBins,]),
                hasCNV = hasCNV,
                Is = which(hasCNV==1),
                N = N,
                fraction = N / ncol(zLayers$raw),
                exp0 = exp0,
                RLL = RLLX,
                p0 = p0X,
                pX = pXX            
            )
        )
    }
    if(nCellsGain > nCellsLoss){ # declare that the output CNV is a gain ...
        return(fillHasCnv( 1, isG, pGG, p0G, RLLG))
    } else { # ... or a loss
        return(fillHasCnv(-1, isL, pLL, p0L, RLLL))
    }
}

#----------------------------------------------------------------------
# call zero to one predominant CNVs in a scan window, at the maximum of Z_se
# i.e. determine if a CNV exists in the window and declare its boundaries
#----------------------------------------------------------------------
# assemble the final call based on the composite over both metrics
assembleCalledCnv <- function(cnvBins, stats){ # 'cnvBins' are zLayers indices (not zL)
    
    # parse initial CNV bin calls
    size_bins <- length(cnvBins)
    if(size_bins == 0) return(NULL) # i.e. if no CNV was called
    cnvStart <- min(cnvBins)
    cnvEnd   <- max(cnvBins)  
    
    # make initial cell assignments
    tmp <- assignCellsToCnv(cnvBins)
    
    # abort if unable to confidently assign any cells to the CNV
    # acts as an additional quality check on validity of the initial CNV call
    # TODO: might speed thing up to use higher threshold values at larger CNV calls? (where cell check is time consuming)
    if(is.null(tmp)) return(NULL)
    if(tmp$cell$N < 2 | # NB: single-cell CNVs will be called elsewhere
       tmp$cell$N < cell$N_accepted * 0.01) return(NULL) # scan seeks CNVs in >=1% of cells (to aid lineage assembly)
    
    # update the CNV boundaries based on the cells with the CNV using HMM
    revisedBins <- refineCnvBoundaries(cnvBins, tmp$cnc, which(tmp$cell$hasCNV == 1))    
    
    # revise CNV-based bin adjustments, and in turn our per-cell count expectations, in the CNV bins
    updateBinExpectations(revisedBins, tmp) # modifies zLayers (not zL)

    # repeat cell assignment using the newly refined CNV boundaries and bin adjustments
    #if(!identical(revisedBins, cnvBins)){ # would use this line if we skip the revised bin adjustments above
        cnvStart <- min(revisedBins)
        cnvEnd   <- max(revisedBins)
        size_bins <- cnvEnd - cnvStart + 1
        tmp <- assignCellsToCnv(revisedBins)          
    #}

    # find the outermost coordinates of the bins called as having the CNV
    cnvStart_bp <- min(bin_data[revisedBins,2])
    cnvEnd_bp   <- max(bin_data[revisedBins,3])
    size_bp     <- cnvEnd_bp - cnvStart_bp + 1

    # assemble and report the called CNV
    list(
            chrom = chrom,
            stats = stats,      
            bin = list(
                Is_initial = cnvBins,
                Is = revisedBins
            ), 
            endpoints = list(
                bins = c(cnvStart,cnvEnd),
                coordinates = c(cnvStart_bp,cnvEnd_bp) 
            ),
            size = list(
                bins = size_bins,
                coordinates = size_bp 
            ),
            region = paste0(chrom, ":", cnvStart_bp, "-", cnvEnd_bp),
            cnc  = tmp$cnc,
            cell = tmp$cell
    )
}
analyzeCandidateCNV <- function(cnvBins_zL, cnvBins_zLayers, Zs, Z_se){

    # calculate a p-value for the set of correlation Z scores
    # TODO: use a different mu than 0??
    corZs <- Zs$Z[cnvBins_zL,cnvBins_zL]
    corZs <- corZs[upper.tri(corZs)]
    p <- if(length(corZs) > 1) t.test(corZs, alternative="greater", mu=0)$p.value else 1    

    # assemble the called bins in a completed CNV call (or reject as a false CNV)
    # TODO: should thresholds be adjusted based on working bin merging??
    if(Z_se > threshold_se | # OR, not AND, is the correct intent
       p < threshold_p){
        assembleCalledCnv(cnvBins_zLayers, list(
            Z_se = Z_se,
            p = p 
        ))
    } else {
        NULL   
    }
}

#----------------------------------------------------------------------
# the callCNVs worker function finds zero to one CVNs on each pass
#   it only calls CNVs that are 100 working bins or less (actual CNV size determined by bin merging by caller)
#   if none, the search ends
#   if one, the CNV is committed and then masked prior to recalling this function
#----------------------------------------------------------------------
callCNVs <- function(zL, iter=1, ignore=NULL){ # zL is a local version of zLayers, used for scanning

    # calculate all bin-bin weighted correlations, with current bin merging and CNV masking
    binCor <- list(values = cor.wt(t(zL$cnc), cell_RPB))
    
    # establish the distribution of correlations according to the null hypothesis
    nullDist <- normalize(binCor$values[!is.na(binCor$values)], 0, 0.2, peakOnly=TRUE)
    
    # calculate mean correlation over all bin pairs within all possible CNV endpoint combinations
    nZLBins <- nrow(zL$raw)
    setPyramidStructure(nZLBins)
    pyramid$sum <- sumPyramid(binCor$values, diag=FALSE)    
    binCor$means <- pyramid$sum / pyramid$N_paired 

    # calculate Z scores for all possible CNV calls relative to the null distribution
    delta <- binCor$means - nullDist$mean # i.e. the excess local correlation
    Zs <- list(
        Z = (binCor$values - nullDist$mean) / nullDist$sd, # Z score for the paired bin values themselves
        Z_sd = delta / nullDist$sd, # Z score for the mean of all paired bin values within the CNV
        Z_se = delta / (nullDist$sd / sqrt(pyramid$N_paired)), # same as above, but calculated using standard error (not deviation)  
        Z_sx = delta / (nullDist$sd / sqrt(pyramid$N_single))  # and again, with standard error by CNB bin # (not paired bin #)
    )
    
    # initialize the matrix used to track failed CNV end pairs
    if(iter == 1) ignore <- matrix(FALSE, nZLBins, nZLBins)

    # find all local maxima of Z_se above a calling threshold, i.e. putative CNV endpoints
    maxima <- findMatrixMaxima(Zs$Z_se, threshold=threshold_se_maxima) # returns a logical matrix
    maxima <- maxima &
                !ignore & 
                pyramid$N_single >= 2 & # need 3 bins to calcuate p-value by t.test (i.e. at least 2 correlations)
                pyramid$N_single <= maxCnvBinsPerLevel & # our scan size max
                Zs$Z_sx >= threshold_sx # find maxima based on Z_se (more precise), but filter based on Z_sx (more stable)
    maximaIs <<- as.data.frame(which(maxima, arr.ind=TRUE)) # one row per local maximum, columns = row,col
    maximaIs$size <<- maximaIs[[2]] - maximaIs[[1]]
    
    # plot summary images for the entire original chromosome
    #if(iter == 1){
        cellStack <- hclust(pearson.dist(t(zL$z0)), method="ward.D2")$order
        plotZLayer(1:nrow(zL$cnc), cellStack, paste(zLFoldMerged,iter,sep="_"), zL=zL)
        plotSquareMatrix(paste('correlation',chrom,zLFoldMerged,iter,sep="_"),
                         1:1, Zs, decorate=addLocalMaxima2)        
    #}

    # abort if no maxima
    if(nrow(maximaIs) == 0) return(NULL)
    
    # start with the smallest candidate CNV until we find one that is significant and called as real
    for(i in order(maximaIs$size)){
        cnvStart <- maximaIs[i,1]
        cnvEnd   <- maximaIs[i,2] 

        message(paste(iter, i, cnvStart, cnvEnd)) 
    
        cnvStart_zLayers <- max(1, (cnvStart - 1) * zLFoldMerged + 1)
        cnvEnd_zLayers   <- min(nChromBins, cnvEnd * zLFoldMerged)
        cnv <- analyzeCandidateCNV(cnvStart:cnvEnd,
                                   cnvStart_zLayers:cnvEnd_zLayers,
                                   Zs, Zs$Z_se[cnvStart,cnvEnd])
        if(is.null(cnv)){
            ignore[cnvStart,cnvEnd] <- TRUE # ignore rejected candidate on the next pass
        } else {
 
            # report the found CNV
            message(paste(cnv$region,
                           commifyInt(cnv$size$coordinates),
                           sum(cnv$cell$hasCNV == 1, na.rm=TRUE),
                           cnv$cnc,
                           cnv$size$bins))
            
            # plot the found CNV
            cellStack <- rev(order(cnv$cnc * cnv$cell$cnc)) # always place CNV at bottom of image
            windowBins <- max(1, min(cnv$bin$Is) - 10):min(nChromBins, max(cnv$bin$Is) + 10)
            plotZLayer(windowBins, cellStack, 'window')
            plotFindCells_RLL(windowBins, cnv)        
            cellTypes <- c( 
                'noCNV',              
                'hasCNV',
                'largerCNV',            
                'ambiguous'
            )
            for(cellType in cellTypes){
                cells <- if(cellType == 'noCNV'){
                    which(cnv$cell$hasCNV == 0) # cells that do not match the CNV
                } else if(cellType == 'hasCNV') {
                    which(cnv$cell$hasCNV == 1) # cells that match the CNV, including sharing both boundaries
                } else if(cellType == 'largerCNV') {
                    which(cnv$cell$hasCNV == -1) # matched CNV, then masked since extends beyond _this_ CNV
                } else {
                    which(is.na(cnv$cell$hasCNV)) # amgibuous cells the algorithm could not call
                }
                if(length(cells) > 0) {
                    cellStack <- cells[rev(order(cnv$cnc * cnv$cell$cnc[cells]))]
                    plotZLayer(windowBins, cellStack, cellType)
                }
            }         

            # "undo"/mask the CNV we just found to prevent it from influencing further scanning
            zLayers <<- maskCalledCnv(zLayers, cnv, cnv$bin$Is, bin_modalCN) # use revised bins
            zL <- collapseZLayers(zLayers, bin_modalCN, zLFoldMerged)

            # then try to find a next CNV
            return(callCNVs(zL, iter + 1, ignore))
        }            
    }
}

#----------------------------------------------------------------------
# main loop: process one chromosome at a time, in parallel
#----------------------------------------------------------------------
#chroms <- 'chr3'
chroms <- c('chr6', 'chr3', 'chr9','chr1','chr18','chrX')
cnvs <- unlist(lapply(chroms, function(chrom){
    message(paste('scanning', chrom))
    chrom <<- chrom

    # load this chromosome's data
    rDataFile <- paste(env$ANALYZE_PREFIX, "layers", chrom, "RData", sep=".")
    load(rDataFile, envir=.GlobalEnv)
    
    # re-initialize script so env etc. are current
    env <- as.list(Sys.getenv())
    source(paste(env$PIPELINE_DIR, "common", "utilities.R", sep="/"))
    command <- "analyze"

    # collect bin information used for parsing CNVs on this chrom
    binIs <- which(bin$data[[1]] == chrom & bin$accepted)
    binSizes_corrected <<- bin$sizes_corrected[binIs] # size corrected for GC and mappability, but NOT CNVs   
    bin_modalCN <<- bin$modalCNsInt[binIs]
    bin_data <<- bin$data[binIs,]
    nChromBins <<- nrow(zLayers$raw)
    chromBinIs <<- 1:nChromBins
    
    # collect cell information for assisting with revised bin adjustments
    cell_RPB <<- cell$readsPerBin[cell$Is_accepted]

    # initiate the recursive CNV calling process
    # NB: zLayers is a global object that is repeatedly masked and merged
    # zL is a local object reflecting the current bin merge level (and masking)
    # zL is used to initially find candidate CNVs
    # zLayers is used to analyze them, in particular to make cell assignments
    zLFoldMerged <<- 1 # the fold factor by which zL is currently merged
    callCNVs(zLayers, bin_modalCN) # first find and mask the smallest CNVs with unmerged bins
    while(nChromBins / zLFoldMerged > maxCnvBinsPerLevel){
        zLFoldMerged <<- zLFoldMerged * foldChangePerLevel # then collapse (another) 10-fold          
        zL <- collapseZLayers(zLayers, bin_modalCN, zLFoldMerged) # inherits the masking applied to zLayers    
        bmcn <- round(collapseVector(bin_modalCN, zLFoldMerged) / zLFoldMerged, 0)
        callCNVs(zL, bmcn) # and try to find larger CNVs
    }
    
    NA
    
}), recursive=FALSE)
#, mc.cores = env$N_CPU


# very different approach would use a multi-state HMM
# hs = all combinations of nStackedCNVs,nCNVCells
# Z_se maxima set tp above 0, in an oriented manner (leading edge of CNV span allows nStackedCNVs to go up)
# Z_sd would somehow need to set ep
# the constant problem is how to set ep's based on total span of CNV (narrowest of all possible triangle peaks?)
#(base) [wilsonte@wilsonte-n1 lib]$ pwd
#/garage/wilsonte_lab/bin/wilson-ware/msvtools-0.0.1/lib
#(base) [wilsonte@wilsonte-n1 lib]$ ls -l run_HMM.R
#-rw-rw-r-- 1 wilsonte wilsonte_lab 14988 May 14  2019 run_HMM.R


    #nAdjBins <- 3
    #nAdjBinsOffset <- nAdjBins - 1
    #adjBinIs <- nAdjBins:(nChromBins-nAdjBins)
    #adjVals <- sapply(adjBinIs, function(i){
    #    bins_left  <- (i-nAdjBinsOffset):i
    #    bins_right <- (i+1):(i+nAdjBinsOffset)
    #    cnc_left   <- colMeans(zLayers$cnc[bins_left,])
    #    cnc_right  <- colMeans(zLayers$cnc[bins_right,])
    #    #raw_left   <- colSums(zLayers$raw[bins_left,])
    #    #raw_right  <- colSums(zLayers$raw[bins_right,])
    #    
    #    
    #    
    #    
    #    
    #    #raw_int_left  <- round(raw_left,  0)
    #    #raw_int_right <- round(raw_right, 0)
    #    #raw_pos_left  <- pmax(raw_left,  0.01)
    #    #raw_pos_right <- pmax(raw_right, 0.01) 
    #    #exp0_left  <- colSums(zLayers$exp0[bins_left,])
    #    #exp0_right <- colSums(zLayers$exp0[bins_right,])
    #    #depthAdj <- exp0_left / exp0_right # account for different cell read depths, including waviness 
    #    c(
    #        cov.wt(data.frame(cnc_left, cnc_right), cor=TRUE, wt=cell_RPB)$cor[1,2],
    #        sum(cnc_right - cnc_left)
    #        #sum(dpois(
    #        #    c(raw_int_left,             raw_int_right), 
    #        #    c(raw_pos_right * depthAdj, raw_pos_left / depthAdj), # note flipped 1/2 order!
    #        #    log=TRUE
    #        #), na.rm=TRUE)               
    #    )
    #})
    #
    #print(str(adjVals))
    #
    #plotAlongBins('adjCor', adjBinIs, adjVals[1,])
    #plotAlongBins('adjPL',  adjBinIs, adjVals[2,])


## set the scan window properties
#nScanBins <- as.integer(env$N_SCAN_BINS) # i.e. the window size in bins
#scanOverlap <- 10 # overlap the bins to prevent missing CNVs at edges
#scanStep <- nScanBins - scanOverlap

#----------------------------------------------------------------------
# define the paired bin scores and dependencies
# functions return one value for one pair of bins defined by the ij vector
#----------------------------------------------------------------------
## how well bins correlate over all cells (bins with no CNVs show no correlation, i.e. 0)
#correlationFN <- function(ij){  # weighted correlation
#    cov.wt(data.frame(cnc[ij[1],],cnc[ij[2],]), cor=TRUE, wt=cellWeights)$cor[1,2]
#}
## create global objects carrying the data needed for calculating paired bin scores
#initializeScoreObjects <- function(bins=TRUE, cells=TRUE){
#    #exp0    <<- zLayers$exp0[bins,cells]
#    #raw_int <<- round(zLayers$raw[bins,cells], 0)   # dpois counts must be integers    
#    #raw_pos <<- pmax(zLayers$raw[bins,cells], 0.01) # dpois lambdas should be >0, otherwise returns -Inf
#    cnc <<- zLayers$cnc[bins,cells]
#    cellWeights <<- cell$readsPerBin[cell$Is_accepted] # make sure low coverage cell don't skew correlations
#}

#----------------------------------------------------------------------
# get the properties of the distribution of paired bin scores throughout a chromosome
# presumption is that most bin pairs do NOT track CNVs, so main peak matches the null hypothesis
#----------------------------------------------------------------------
#getNullDistribution <- function(nBins){
#    nUniqPairs <- (nBins ** 2 - nBins) / 2 # how many possible correlations on the chromosome
#    if(nUniqPairs < nullDistributionPairs * 2){
#        setPyramidStructure(nBins) # more efficient to just calculate all pairs for small chromosomes
#        m <- applyPyramidPairs(correlationFN)
#    } else {
#        m <- matrix(NA, nBins, nBins) # work in chunks until we have enough samples
#        pairsPerIter <- min(nBins, 1000)
#        while(length(m[!is.na(m)]) < nullDistributionPairs){
#            mapply(function(i, j){
#                if(j > i) m[i,j] <<- correlationFN(c(i,j)) 
#            }, sample(1:nBins, pairsPerIter), sample(1:nBins, pairsPerIter) )     
#        }          
#    }
#    normalize(m[!is.na(m)], 0, 0.2, peakOnly=TRUE) # find the main peak of values
#}

#----------------------------------------------------------------------
# calculate paired scores, and associated Z-scores, for all bin pairs in a scan window
#----------------------------------------------------------------------
#processPairedBins <- function(){
#    
#    # list of square output matrices, i.e. bin x bin
#    # start by calculating the score values for each pair of bins in the window
#    x <- list( value = applyPyramidPairs(correlationFN) ) # diagonal is NA
#
#    # calculate mean values over all bin pairs within all possible CNV endpoint combinations
#    pyramid$sum <- sumPyramid(x$value, diag=FALSE)
#    x$mean <- pyramid$sum / pyramid$N_paired    
#
#    # calculate and return Z scores for CNV calls relative to the null distribution
#    delta <- x$mean - nullDist$mean
#    list(
#        Z = (x$value - nullDist$mean) / nullDist$sd, # Z score for the paired bin values themselves
#        Z_sd = delta / nullDist$sd, # Z score for the mean of all paired bin values within the CNV
#        Z_se = delta / (nullDist$sd / sqrt(pyramid$N_paired)) # same as above, but calculated using standard error (not deviation)  
#    )
#}

    ## calculate all bin-bin correlations in the window
    #initializeScoreObjects(bins=windowBins)
    #setPyramidStructure(nScanBins)
    #cors <- processPairedBins()   
    #
    ## find maximum in Z_se, which maps to putative CNV endpoints
    #max_se <- max(cors$Z_se, na.rm=TRUE) # one absolute maximum (to find only one CNV)
    #XY_at_max <- which(cors$Z_se == max_se, arr.ind=TRUE)
    #cnvBins <- XY_at_max[1,1]:XY_at_max[1,2] 

    ## map out the bin window strategy based on preferred window width
    #windowStarts <- seq(1, nBins - nScanBins + 1, scanStep)
    #windowEnds   <- windowStarts + nScanBins - 1
    #nWindows <- length(windowStarts)
    #windowStarts[nWindows] <- windowStarts[nWindows] - (windowEnds[nWindows] - nBins)
    #windowEnds[nWindows] <- nBins # shift last window left (i.e. more overlap) to avoid overstepping end of chrom

    ## use random pair sampling to define a population mean and std dev. for paired bin scores
    ## don't do this per window, as some windows will be entirely filled with a CNV!
    ##   of course, aneuploidy might fill an entire chromosome, which becomes the baseline correlation, etc.
    ##   in future, might want to eliminate whole chromosome aneuploidies prior to chromosome scanning
    #initializeScoreObjects()
    #nullDist <<- getNullDistribution(nBins)
    #
    ## scan for CNVs, one window at a time, using paired bin properties relative to null distributions
    #chromCnvs <- list()
    #for(i in 1:nWindows){
    #    windowBins <- windowStarts[i]:windowEnds[i] # as zLayer indices
    #    cnv <- callCNV(windowBins)
    #    if(!is.null(cnv)) chromCnvs <- c(chromCnvs, cnv)
    #}
    ##message(paste('done with', chrom))

#message(length(cnvs))


## how well bins predict each others cell values, expressed as a crossed likelihood
#pairedLikelihoodFN <- function(ij){ 
#    depthAdj <- exp0[ij[1],] / exp0[ij[2],] # account for different cell read depths, including waviness 
#    sum(dpois(
#        c(raw_int[ij[1],],            raw_int[ij[2],]), 
#        c(raw_pos[ij[2],] * depthAdj, raw_pos[ij[1],] / depthAdj), # note flipped 1/2 order!
#        log=TRUE
#    ), na.rm=TRUE)   
#}
# each score metric can nominate one CNV
#callCNVsByMetric <- function(sType, pairs){
#    
#    # invert Z sign on paired likelihoods
#    # since CNV-positive values are lower than null hypothesis (correlations are higher)
#    sign <- if(sType == 'correlation') 1 else -1
#    
#    # find maxima in Z_se as data frames of points
#    #lmax   <- findMatrixMaxima(sign * pairs$Z_se, threshold=threshold_se) # local maxima
#    #lmaxXY <- which(lmax, arr.ind=TRUE)
#    max_se  <- max(sign * pairs$Z_se, na.rm=TRUE) # one absolute maximum (to find only one CNV)
#    maxXY   <- which(pairs$Z_se == sign * max_se, arr.ind=TRUE)
#    max_sd  <- abs(pairs$Z_sd[maxXY[1,1],maxXY[1,2]])
#
#    # if passing threshold, report as a found CNV
#    boolean <- rep(0, nScanBins) # a span of bins, without or without a CNV call as 1's
#    size <- abs(maxXY[1,2] - maxXY[1,1] + 1)
#    range <- maxXY[1,1]:maxXY[1,2]
#    pairsZ <- if(sType == 'correlation'){
#        x <- pairs$Z[range,range]
#        x[upper.tri(x)] # each pairwise combination is unique information
#    } else {
#        if(size == nScanBins){
#            NA
#        } else {
#            x <- mirrorPyramid(pairs$Z)
#            x[range,range] <- NA 
#            rowMeans(x[range,,drop=FALSE], na.rm=TRUE) # otherwise each unusual bin is compared to many ==> artificially high N          
#        }
#    }  
#    p <- if(length(pairsZ) > 1) t.test(sign * pairsZ, alternative="greater", mu=0)$p.value else 1
#    if(max_se > threshold_se &
#       max_sd > threshold_sd[size] &
#       p      < threshold_p[size]){
#        boolean[maxXY[1,1]:maxXY[1,2]] <- 1
#    }
#    list(boolean=boolean, stats=list(Z_sd=max_sd, Z_se=max_se, p=p)) 
#}

