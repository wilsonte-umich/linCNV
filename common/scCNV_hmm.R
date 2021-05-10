
#----------------------------------------------------------------
# Hidden Markov Model: wholly generic functions
#----------------------------------------------------------------

# control transition probabilities
HMM_set_persistence <- function(persistence, nHMMStates){
    persistence <<- persistence # i.e. the probability of remaining in state
    tp <<- matrix(log((1-persistence)/(nHMMStates-1)), nHMMStates, nHMMStates)
    for(i in 1:nHMMStates) tp[i,i] <<- log(persistence)
}

# solve HMM by viterbi algorithm
HMM_viterbi_ep <- function(ep){ 
    ep[is.na(ep)] <- -Inf
    
    # 1. initialization (observation t=1)
    T         <- nrow(ep) # length of the sequence of observations
    N         <- ncol(ep) # number of states
    delta     <- log(matrix(0, nrow=T, ncol=N))
    delta[1,] <- sapply(1:N, function(i) log(1/N) + ep[1,i])
    phi       <- matrix(NA, nrow=T, ncol=N)
    
    # 2. recursion;
    # NB: these 'for' loops are faster than apply methods with array as implemented and given recursion restrictions
    for (t in 2:T){
        pt <- t - 1
        for (j in 1:N){   # j = this hs
            ep_j <- ep[t,j]
            for (i in 1:N){ # i = prev hs
                delta_ <- delta[pt,i] + tp[i,j] + ep_j
                if(delta[t,j] < delta_){
                    delta[t,j] <- delta_
                    phi[pt,j]  <- i
                }
            }
        }
    }
    
    # 3. termination
    prob <- -Inf
    hsi  <- rep(1, T)
    for (j in 1:N){
        if(prob < delta[T,j]){
            prob <- delta[T,j]
            hsi[T] <- j
        }
    }
    
    # 4. reconstruction and return the hidden state indices
    for (t in (T-1):1) hsi[t] <- phi[t,hsi[t+1]]
    hsi
}

#----------------------------------------------------------------
# run HMM on the bins from one chromosome (or a subset) in one cell
#----------------------------------------------------------------
 
#HMM_viterbi <- function(cellI, chromBins){ 
#  maxCount <- round(2 * zLayers$expG[chromBins,cellI], 0)
#  ep <- sapply(c('expL','exp0','expG'), function(exp){ # set emission probabilities
#    log(dpois(pmin(maxCount, round(zLayers$raw[chromBins,cellI],0)), zLayers[[exp]][chromBins,cellI]))
#  })
#  HMM_viterbi_ep(ep)
#}
#HMM_viterbi <- function(zL, bins, cell){
#    # set emission probabilities and solve CN HMM
#    maxCount <- 2 * zL$exp[[maxCN+1]][bins, cell]
#    raw <- round(pmin(zL$raw[bins, cell], maxCount), 0)  
#    ep <- sapply(0:maxCN, function(cn) log(dpois(raw, zL$exp[[cn+1]][bins, cell])) )
#    cn <- HMM_viterbi_ep(ep) - 1
#    # extract runs of common CN
#    x <- rle(cn)
#    d$run <- data.frame(length=x$lengths, cn=x$values)
#    d$run$end_   <- cumsum(d$run$length) # endpoints of the HMM segments
#    d$run$start_ <- c(1,(d$run$end_ + 1)[1:length(d$run$end_)-1])
#    d$run$start  <- d$run$start_ + bin1I_zL - 1 # put bin numbers back to zLayers indices
#    d$run$end    <- d$run$end_   + bin1I_zL - 1
#    d$run[,c('sign','sign_exp')] <- sapply(1:nrow(d$run), function(i){
#        range <- d$run[i,'start_']:d$run[i,'end_']
#        c( round(mean(d$sign    [range]), 0), # i.e. neutral (0), gain (+1) or loss (-1)
#           round(mean(d$sign_exp[range]), 0) )
#    }) 
#}
# run HMM on all chromosomes of a single cell
HMMByChrom <- function(cellI, zLayerChroms){ 
    unlist(mclapply(chroms, function(chrom){
      HMM_viterbi(cellI, which(zLayerChroms == chrom))
    }, mc.cores = env$N_CPU))
}

#----------------------------------------------------------------
# solve a CN HMM for counts summed by sumZLayerByBin or sumZLayerByCell
# expects global value for maxCN, and that HMM_set_persistence was called for states 0:maxCN
#----------------------------------------------------------------

# solveHMM for a single vector of ordered values (bins or cells)
solveHMM_sums <- function(sums, refCN=0){
    
    # set emission probs and run HMM
    maxCount <- 2 * sums$exp[[maxCN+1]] 
    raw <- round(pmin(sums$raw, maxCount), 0)    
    ep <- sapply(sums$exp, function(exp) log(dpois(raw, exp)) ) 
    d <- list( # one assigned CN value for each bin or cell
        cn = pmax(pmin(HMM_viterbi_ep(ep) - 1 + refCN, maxCN), 0)
    )
    d$cnc_modal  <- d$cn - sums$modalCN
    d$cnc_ploidy <- d$cn - sums$ploidy     
    d$sign <- sign(d$cnc_ploidy) # i.e. gain or loss relative to ploidy
    
    # extract runs of common CN
    x <- rle(d$cn)
    d$run <- data.frame(length=x$lengths, cn=x$values)
    d$run$end_   <- cumsum(d$run$length) # endpoints of the HMM segments as bin or cell indices in sums
    d$run$start_ <- c(1,(d$run$end_ + 1)[1:length(d$run$end_)-1])
    d$run[,c('sign','ploidy')] <- sapply(1:nrow(d$run), function(i){
        range <- d$run[i,'start_']:d$run[i,'end_']
        c( round(mean(d$sign[range]), 0), # i.e. neutral (0), gain (+1) or loss (-1)
           round(mean(sums$ploidy[range]), 0) )
    })
    cnvRuns <- d$run$cn != d$run$ploidy
    d$nCnvRuns  <- nrow(d$run[cnvRuns,])
    d$nCnvItems <- sum(d$run[cnvRuns,'length'])
    
    d
}

# HMM for a set of bins (i.e. one value per bin)
solveHMM_binSums <- function(sums, bin1I_zL=NULL, # bin1I_zL is the zLayers index of sums[1]
        queryRange=NULL, refCN=0, expectedCN=NA){ # query is an internal bin span we care about
    if(is.null(bin1I_zL)) bin1I_zL <- min(sums$bins)
    
    # solve the HMM
    d <- solveHMM_sums(sums, refCN)
    d$run$start  <- d$run$start_ + bin1I_zL - 1 # put bin numbers back to zLayers indices
    d$run$end    <- d$run$end_   + bin1I_zL - 1    
    
    # determine which runs match the query
    if(is.null(queryRange)) queryRange <- range(1:length(sums$raw) + bin1I_zL - 1)
    d$run$overlaps    <- queryRange[1] <= d$run$end   & d$run$start <= queryRange[2]
    d$run$internal    <- queryRange[1] <  d$run$start & d$run$end   <  queryRange[2]    
    #d$run$matchesCN   <- d$run$cn == expectedCN
    #d$run$matchesSign <- d$run$sign == d$run$sign_exp
    #d$run$isMatch     <- d$run$matchesCN   & d$run$overlaps
    #d$run$isSignMatch <- d$run$matchesSign & d$run$overlaps
    #d$hasGaps <- sum(!d$run$matchesSign & d$run$internal) > 0 # TRUE if HMM not continuous through query
    
    # establish runs of interest that overlap the query region
    d$ROI <- which(d$run$overlaps & d$run$cn != ploidy)
    #d$ROI <- which(d$run$overlaps) 
    d$ROI <- if(length(d$ROI) == 0) d$run[FALSE,] # all segments between the outermost overlapping CNV segments
             else d$run[min(d$ROI):max(d$ROI),]   # includes internal CNC=0 segments
    d$nROI <- nrow(d$ROI)
    d$ROI_range <- if(d$nROI > 0) as.integer(c(min(d$ROI$start), max(d$ROI$end))) else NULL
    
    ## establish CNV segments that could explain the ROI (segments and ROI are NOT synonymous)
    #transitions <- diff(c(sums$ploidy[min(d$ROI$start)], # include transitions into and out of the outermost runs
    #                      d$ROI$cn,
    #                      sums$ploidy[max(d$ROI$end)] ))
    #nTransitions <- length(transitions)    
    #d$segments <- cbind(expand.grid(1:nTransitions, 1:nTransitions),
    #                    expand.grid(transitions,    transitions))
    #colnames(d$segments) <- c('i','j','startSign','endSign')
    #d$segments <- d$segments[d$segments$j > d$segments$i &
    #                     d$segments$startSign != 0 &
    #                     d$segments$endSign   != 0 &
    #                     sign(d$segments$startSign) != sign(d$segments$endSign),]    
    #d$segments$start <- d$ROI[d$segments$i,  'start']
    #d$segments$end   <- d$ROI[d$segments$j-1,'end']
    
    ## mask the non-ROI CNVs from the cnc map
    #d$cn_roi <- d$cn # result only carries CN for runs that overlap the query region; is true CN
    #nonRoi_runs <- which(!d$run$overlaps & d$run$cn != ploidy)
    #for(i in nonRoi_runs) d$cn_roi[d$run$start_[i]:d$run$end_[i]] <- ploidy

    d  
}

# HMM for a set of cells (i.e. one value per cell)
solveHMM_cellSums <- function(sums, cellStack){

    # solve the HMM
    d <- solveHMM_sums(sums)

    # establish groups of cells with CNVs / aneuploidy
    d$ROI <- which(d$run$cn != ploidy)
    d$ROI <- if(length(d$ROI) == 0) d$run[FALSE,] # all runs of cells with a CNV
             else d$run[d$ROI,] 
    d$nROI <- nrow(d$ROI)
    
    d
}

## solve a CNC HMM of each individual cell vs. a CN state model
#solveRegionHMM_vs_model <- function(bins, model){ # bins are zLayers indices, model$cn same nBins as zLayers
#    refCN <- model$cn[bins]
#    cn_offset <- refCN - bin$ploidy[bins]
#    bins_ <- 1:length(cn_offset)
#    sapply(1:cell$N, function(j){
#        exp <- sapply(0:maxCN, function(cn_out){
#            expI <- cn_out + cn_offset + 1            
#            sapply(bins_, function(i){            
#                if(expI[i] < 1 | expI[i] > maxCN+1) NA else zLayers$exp[[expI[i]]][i,j]  
#            })
#        }, simplify=FALSE)
#        maxCount <- 2 * exp[[maxCN+1]] 
#        raw <- round(pmin(zLayers$raw[bins,j], maxCount), 0)        
#        ep <- sapply(exp, function(exp_) log(dpois(raw, exp_)) )
#        cn <- pmax(pmin(HMM_viterbi_ep(ep) - 1 + cn_offset, maxCN), 0)
#        identical(cn, refCN) # return a logical for each cell; TRUE means the cell matched model at all bins  
#    })
#}


#===========================================================

## solve the HMM for previously summed counts for a set of bins (i.e. one value per bin)
#solveRegionHMM <- function(binSums, bin1I, # bin1I_zL is the zLayers index of binSums[1]
#                           queryRange=NULL, refCNC=0, expectedCNC=NA){ # query is an internal bin span we care about
#    if(is.null(queryRange)) queryRange <- range(1:length(binSums$raw) + bin1I_zL - 1)
#
#    # set emission probs and run HMM on merged cells
#    maxCount <- round(2 * binSums$expG, 0) 
#    ep <- sapply(c('expL','exp0','expG'), function(exp){
#      log(dpois(pmin(maxCount, round(binSums$raw,0)), binSums[[exp]]))
#    })
#    d <- list( # one assigned CNC value for each padded bin  
#        cnc = pmax(-1, pmin(1, HMM_viterbi_ep(ep) + refCNC))
#    )
#
#    # extract runs of common CNC; abs(CNC) can be > 1 if refCNC is provided!
#    x <- rle(d$cnc)
#    d$run <- data.frame(length=x$lengths, cnc=x$values)
#    d$run$end_   <- cumsum(d$run$length) # endpoints of the HMM segments
#    d$run$start_ <- c(1,(d$run$end_ + 1)[1:length(d$run$end_)-1])
#    d$run$start  <- d$run$start_ + bin1I_zL - 1 # put bin numbers back to zLayers indices
#    d$run$end    <- d$run$end_   + bin1I_zL - 1
#
#    # determine which runs match the query
#    d$run$matchesCNC <- d$run$cnc == expectedCNC
#    d$run$overlaps   <- queryRange[1] <= d$run$end & d$run$start <= queryRange[2]
#    d$run$matches    <- d$run$matchesCNC & d$run$overlaps
#    d$ROI <- which(d$run$overlaps & d$run$cnc != 0) # runs of interest
#    d$ROI <- if(length(d$ROI) == 0){
#        d$run[FALSE,]
#    } else { # all segments between the outermost overlapping CNV segments
#        d$run[min(d$ROI):max(d$ROI),] # includes internal CNC=0 segments
#    }
#    d$nROI <- nrow(d$ROI)
#    
#    # mask the non-ROI CNVs from the cnc map
#    d$cnc_roi <- d$cnc # result only carries CNC for runs that overlap the query region
#    for(i in which(!d$run$overlaps & d$run$cnc != 0)) d$cnc_roi[d$run$start_[i]:d$run$end_[i]] <- 0
#    d  
#}

