
#----------------------------------------------------------------
# miscellaneous functions
`%notin%` <- Negate(`%in%`)
#----------------------------------------------------------------

#----------------------------------------------------------------
# reporting functions
#----------------------------------------------------------------
commify    <- function(x) format(x, big.mark=",", scientific=FALSE)
commifyInt <- function(x) commify(round(x, 0))
reportCount <- function(count, lbl, msg){ # pipeline log feedback
    message(paste(command, commifyInt(count), lbl, msg, sep="\t"))
}
reportFraction <- function(fraction, lbl, msg, round=1){ # pipeline log feedback
    message(paste(command, round(fraction, round), lbl, msg, sep="\t"))
}
reportString <- function(str, lbl, msg){ # pipeline log feedback
    message(paste(command, str, lbl, msg, sep="\t"))
}
head.table <- function(table, nrow=10, ncol=10){
    print(table[1:nrow,1:ncol])
}
setStartTime <- function() startTime <<- proc.time()
reportEndTime <- function(message){ # for timing actions during development
    message(message)
    print(proc.time() - startTime)
    setStartTime()
}

#----------------------------------------------------------------
# plotting functions
#----------------------------------------------------------------
plotDim   <- 1800
pointsize <- 7
plotRes   <- 600
plotHistogram <- function(v, type, xlab,
                          xlim=range(v), vLine=NULL, breaks=100){
    v <- v[!is.na(v)]
    p <- v[v >= xlim[1] & v <= xlim[2]]
    file <- paste(env$PLOT_PREFIX, env$PIPELINE_COMMAND, type, "png", sep=".")
    png(file, width=plotDim, height=plotDim,
        units="px", res=plotRes, pointsize=pointsize, type = c("cairo"))
    hist(p, breaks=breaks, main="", xlab=xlab, xlim=xlim, ylab="Frequency")
    if(!is.null(vLine)) abline(v=vLine, col="blue")
    graphics.off()
}
plotCorrelation <- function(x, y, type, xlab, ylab, xlim=range(x), ylim=range(y), cex=1,
                            vLine=NULL, hLine=NULL, fit=NULL){
    file <- paste(env$PLOT_PREFIX, env$PIPELINE_COMMAND, type, "png", sep=".")
    png(file, width=plotDim, height=plotDim,
        units="px", res=plotRes, pointsize=pointsize, type = c("cairo"))  
    plot(x, y, pch=16, cex=0.5, main="", xlab=xlab, xlim=xlim, ylab=ylab, ylim=ylim)
    if(!is.null(vLine)) abline(v=vLine, col="blue")
    if(!is.null(hLine)) abline(h=hLine, col="blue")
    if(!is.null(fit)) {
      x <- seq(xlim[1], xlim[2], (xlim[2] - xlim[1]) / 100 )
      y <- predict(fit, data.frame(x=x))
      lines(x, y, col="red", lwd=2)
    }
    graphics.off()
}

#----------------------------------------------------------------
# chromosome functions
#----------------------------------------------------------------
getOrderedChroms <- function(chroms, includeY=TRUE){
  chroms_ <- sub("chr", "",  unique(chroms))
  isAutosome <- !is.na(suppressWarnings(as.numeric(chroms_)))
  chroms <- paste("chr", c(sort(as.numeric(chroms_[isAutosome])), sort(chroms_[!isAutosome])), sep="")
  if(!includeY) chroms <- chroms[!(chroms %in% "chrY")]
  chroms
}
getAutosomeIs <- function(chroms){ # logical vector of chromosomes that are autosomes
    numChroms  <- sub("CHR", "", toupper(chroms)) # make sure chroms are 1,2...X,Y
    suppressWarnings(!is.na(as.numeric(numChroms)))    
}

#----------------------------------------------------------------
# vector functions
#----------------------------------------------------------------
collapseVector <- function(v, n) { # sum every n adjacent elements of a vector
    cv <- unname(tapply(v, (seq_along(v)-1) %/% n, sum))
    tailLength <- length(v) %% n # number of input elements summed into incomplete last output element    
    if(tailLength != 0){
        cvLength <- length(cv) # expand incomplete last element to same scale as all others
        cv[cvLength] <- cv[cvLength] * n / tailLength          
    }
    cv
}
uncollapseVector <- function(v, n, len) { # reverse the actions of collapseVector
    ucv <- as.vector(sapply(v, rep, n))
    extra <- length(ucv) - len # user must remember how long the original vector was    
    if(extra > 0) ucv <- ucv[1:len]    
    ucv
}
ma <- function(x, n = 5){ # moving average
    filter(x, rep(1 / n, n), sides = 2)
}

#----------------------------------------------------------------
# zLayer functions
#----------------------------------------------------------------
fillZLayerValues <- function(zL, CNInt){ # requires only that raw and exp0 have been set
  CNInt   <- matrix(rep(CNInt, cell$N_accepted), ncol=cell$N_accepted)
  rpa     <- zL$exp0 / CNInt # i.e. reads per allele  
  zL$cn   <- zL$raw / rpa
  zL$cnc  <- zL$cn - CNInt
  zL$expG <- zL$exp0 + rpa
  zL$expL <- zL$exp0 - rpa
  zL$expL <- pmax(zL$expL, minExp) # prevent divide by zero
  zL$z0   <- (zL$raw - zL$exp0) / sqrt(zL$exp0) # Poisson variance = mean
  zL$zG   <- (zL$raw - zL$expG) / sqrt(zL$expG)
  zL$zL   <- (zL$raw - zL$expL) / sqrt(zL$expL)
  zL
}
collapseZLayers <- function(zL, CNInt, n){
    if(n == 1) return(zL)
    zL$exp0 <- apply(zL$exp0, 2, collapseVector, n)   # sum of expected and observed counts
    zL$raw  <- apply(zL$raw , 2, collapseVector, n)
    CNInt   <- round(collapseVector(CNInt, n) / n, 0) # mean of CNInt     
    fillZLayerValues(zL, CNInt)
}

#----------------------------------------------------------------
# 3-state Hidden Markov Model based on Poisson emissions
#----------------------------------------------------------------
nHMMStates <- 3
HMM_set_persistence <- function(persistence){
    persistence <<- persistence
    tp <<- matrix(log((1-persistence)/(nHMMStates-1)), nHMMStates, nHMMStates)
    for(i in 1:nHMMStates) tp[i,i] <<- log(persistence)
}
HMM_viterbi <- function(cellI, chromBins){ # run HMM on the bins from one chromosome in one cell

  # 0. set emission probabilities
  maxCount <- round(2 * zLayers$expG[chromBins,cellI], 0)
  ep <- sapply(c('expL','exp0','expG'), function(exp){
    log(dpois(pmin(maxCount, round(zLayers$raw[chromBins,cellI],0)), zLayers[[exp]][chromBins,cellI]))
  })
  
  # run algorithm
  HMM_viterbi_ep(ep)
}
HMM_viterbi_ep <- function(ep){ # run HMM on the bins from one chromosome in one cell

  # 1. initialization (observation t=1)
  T         <- nrow(ep) # length of the sequence of observations
  N         <- ncol(ep)
  delta     <- log(matrix(0, nrow=T, ncol=N))
  delta[1,] <- sapply(1:N, function(i) log(1/N) + ep[1,i])
  phi       <-     matrix(NA, nrow=T, ncol=N)

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
  hsi  <- rep(-1, T)
  for (j in 1:N){
    if(prob < delta[T,j]){
      prob <- delta[T,j]
      hsi[T] <- j
    }
  }

  # 4. reconstruction
  for (t in (T-1):1) hsi[t] <- phi[t,hsi[t+1]]

  # return loss(-1), neutral(0), gain(+1)
  hsi - (nHMMStates - 1)
  # list(prob=prob, hsi=hsi) # overall likelihood and ordered set of most likely states
}
HMMByChrom <- function(cellI, zLayerChroms){ # run HMM on all chromosomes of a single cell
    unlist(mclapply(chroms, function(chrom){
      HMM_viterbi(cellI, which(zLayerChroms == chrom))
    }, mc.cores = env$N_CPU))
}

#----------------------------------------------------------------
# distribution functions
#----------------------------------------------------------------
# the position of the main (max) peak in a distribution
distributionPeak <- function(v, weights=NULL){ 
  weights <- if(!is.null(weights)) weights <- weights/sum(weights)
  dens <- density(v, weights=weights)
  i <- which.max(dens$y)
  list(dens=dens, x=dens$x[i], y=dens$y[i])
}
# return a distribution as individual Z scores relative to the main distribution peak
normalize <- function(d, peakCenter=NULL, peakWidth=NULL, peakOnly=FALSE){
    
    # use R 'density' function to examine the value distribution
    dens <- density(d[d > -Inf & d < Inf], na.rm=TRUE)
    
    # find the main peak and make initial mean+SD guesses    
    peak <- if(is.null(peakCenter)){
        iAtMaxDens <- which.max(dens$y)     
        maxDens <- dens$y[iAtMaxDens]
        valAtMaxDens <- dens$x[iAtMaxDens]
        if(is.null(peakWidth)) peakWidth <- abs(valAtMaxDens / 5) 
        dens$x[which(dens$y >= maxDens * 0.2 &
                     abs(dens$x - valAtMaxDens) <= peakWidth   )]        
    } else {
        dens$x[which(abs(dens$x - peakCenter) <= peakWidth   )]  
    }
    median <- median(peak, na.rm=TRUE)
    sd <- sd(peak, na.rm=TRUE)
    
    # fatten out around the initial guess and calculate refined mean+SD values
    delta <- d - median    
    peak <- d[delta > sd * -4 & delta < sd * 4] 
    median <- median(peak, na.rm=TRUE)
    sd <- sd(peak, na.rm=TRUE)
    
    # return all results
    N <- length(na.omit(as.vector(d)))
    if(peakOnly){
        list( 
            median = median,
            mean = median,
            sd = sd,
            N = N
        )    
    } else {
        list( 
            data = d,
            dens = dens,
            median = median,
            mean = median,
            sd = sd,
            N = N,
            Z = (d - median) / sd
        )         
    }
}

#----------------------------------------------------------------
# matrix functions
#----------------------------------------------------------------
# construct a set of tools useful for filling a matrix with stacked values (i.e. a pyramid)
setPyramidStructure <- function(N){ # N is the size of the expected vector of element values
    m_single <- diag(N) # identity matrix
    m_paired <- matrix(0,N,N)
    m_paired[upper.tri(m_paired)] <- 1
    pyramid <<- list(
        #m_paired = m_paired,
        N_single = sumPyramid(m_single, diag=TRUE), # number of diagonal elements summed in simple single-value stacks 
        N_paired = sumPyramid(m_paired, diag=FALSE), # number element pair values summed in complex paired-value stacks
        N = N
        #, # size of one axis of the sqaure matrix
        #pairs = combn(1:N,2), # all possible element index combinations, for apply functions
        #selfPairs = matrix(1:N,2,N,byrow=TRUE)
    )
    #pyramid$N_outerBlocks <<- sumOuterBlocks(pyramid$N_paired, pyramid$m_paired, diag=TRUE)
}
# sum a matrix of (paired) values into a new pyramid (i.e. triangle)
# assumes values are in the upper triangle +/- diagonal; lower triangle should be NA
sumPyramid <- function(m, diag=TRUE){
    m[is.na(m)] <- 0 # to allow cumsum to work properly     
    x <- m 
    x[1,] <- cumsum(colSums(m)) # initialize first row of output x
    for(i in 2:nrow(m)) x[i,] <- x[i-1,] - cumsum(m[i-1,]) # iteratively fill remaining rows
    if(!diag) diag(x) <- NA # ensure output again matches the upper triangle structure
    x[lower.tri(x)] <- NA
    x
}
# sum the rectangles extending from an internal triangle to the outer matrix edges
sumOuterBlocks <- function(pyramidSum, pyramidVal, diag=FALSE){
    diag(pyramidSum) <- 0
    all <- pyramidSum[1,pyramid$N]  # the biggest triangle sum (i.e. all pairs)
    applyPyramidPairs(function(ij){ # this leaves the rectangles where we expect bad paired likelihood
        thisTri  <- pyramidSum[ij[1],ij[2]] # this triangle's sum (where expect acceptable paired likelihood)         
        leftTri  <- if(ij[1] > 1) pyramidSum[1,ij[1]-1] else 0 # triangles away from us
        rightTri <- if(ij[2] < pyramid$N) pyramidSum[ij[2]+1,pyramid$N] else 0
        outerRect <- if(ij[1] > 1 & ij[2] < pyramid$N) { # the rectangle toward the apex of the big triangle
            sum(pyramidVal[1:(ij[1]-1),(ij[2]+1):pyramid$N], na.rm=TRUE)
        } else {
            0
        }
        all - thisTri - leftTri - rightTri - outerRect
    }, diag)   
}
# fill a pyramid with values calculated from paired elements
applyPyramidPairs <- function(FN, diag=FALSE){ # function must expect ij vector
    m <- matrix(NA, pyramid$N, pyramid$N)
    m[lower.tri(m)] <- apply(pyramid$pairs, 2, FN)
    if(diag) diag(m) <- apply(pyramid$selfPairs, 2, FN)
    t(m) # place results in upper triangle
}
# copy the upper into the lower triangle of a matrix and set the diagonal
mirrorPyramid <- function(m, diagVals=NULL){
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    if(!is.null(diagVals)) diag(m) <- diagVals
    m
}
twoWayCumSum <- function(m, diagVals=0){
    m <- mirrorPyramid(m, diagVals)
    m <- t(apply(m, 1, cumsum))
    apply(m, 2, cumsum)       
}

#----------------------------------------------------------------
# translate a matrix by one step in the left,right,up,down directions
# fill the now-empty spaces as requested
shiftMatrix <- function(m, direction="right", fill=NA){
    switch(direction,
        "right" = cbind(fill, matrix(as.vector(m)[1:(length(m)-nrow(m))], nrow(m))),
        "left"  = cbind(matrix(as.vector(m)[(nrow(m)+1):length(m)], nrow(m)), fill),
        "down"  = t(shiftMatrix(t(m), "right")),
        "up"    = t(shiftMatrix(t(m), "left"))
    )
}
# take a moving average of a matrix of a cell and all cells adjacent to it
# thus, returns a moving average of a 9-cell (3 x 3) square
ma33Matrix <- function(m){
    rgt <- shiftMatrix(m, "right")    
    lft <- shiftMatrix(m, "left")
    up_ <- shiftMatrix(m, "up")    
    dwn <- shiftMatrix(m, "down")
    rgt_up_ <- shiftMatrix(rgt, "up")
    rgt_dwn <- shiftMatrix(rgt, "down")
    lft_up_ <- shiftMatrix(lft, "up")
    lft_dwn <- shiftMatrix(lft, "down")
    apply(
        abind(m, rgt, lft, dwn, up_, rgt_up_, rgt_dwn, lft_up_, lft_dwn, along=3),
        c(1,2), mean, na.rm=TRUE
    )
}
# take the second derivative of either the rows (1) or columns (2) of a matrix
# NB: only the _sign_ of the value is correct!
deriv2Matrix <- function(m, axis=1:2){
    switch(axis,
        t(deriv2Matrix(t(m), 2)),           
        rbind(NA, diff(sign(diff(m))), NA) # diff acts by column
    ) 
}
# determine which cells in a matrix represent local maxima
findMatrixMaxima <- function(m, ma33=FALSE, threshold=NULL){
    if(ma33) m <- ma33Matrix(m) # offer a bit of smoothing for stability
    LOW <- -1e9
    M <- cbind(LOW, rbind(LOW, m, LOW), LOW) # padding allows finding of maxima at edges
    M[is.na(M)] <- LOW # supports use for triangles or other partial matrices
    fxx <- deriv2Matrix(M, 1)
    fyy <- deriv2Matrix(M, 2)
    boo <- fxx * fyy > 0 & # Hessian test to avoid saddle points
           fxx < 0 & # negative 2nd derivative identifies maxima
           if(is.null(threshold)) TRUE else M >= threshold # enforce minimum peak height 
    boo[(1:nrow(m))+1,(1:ncol(m))+1]
}

#----------------------------------------------------------------
# correlation-based distance of a set of already-centered Z scores
#----------------------------------------------------------------
# correlation-based distance of a set of already-centered Z scores
pearson.dist <- function (m) { # m is a matrix, computes distance between rows
  m <- m / sqrt (rowSums (m^2))
  m <- tcrossprod (m)
  m <- as.dist(m)
  0.5 - m / 2
}
pearson.matrix <- function (m) { # return the full matrix, not a dist formatted object
  m <- m / sqrt (rowSums (m^2))
  m <- tcrossprod (m)
  0.5 - m / 2 
}

#----------------------------------------------------------------
# weighted Pearson correlations of CNC values between columns (typically bins)
#----------------------------------------------------------------
cor.wt <- function(m, weights){ # weights are typically cell depths
    m <- cov.wt(m, cor=TRUE, wt=weights)$cor
    diag(m) <- NA # use upper triangle only
    m[lower.tri(m)] <- NA
    m
}

