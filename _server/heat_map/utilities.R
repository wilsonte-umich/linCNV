
# vector functions
collapseVector <- function(v, n) { # sum every n adjacent elements of a vector
    cv <- unname(tapply(v, (seq_along(v)-1) %/% n, sum))
    tailLength <- length(v) %% n # number of input elements summed into incomplete last output element    
    if(tailLength != 0){
        cvLength <- length(cv) # expand incomplete last element to same scale as all others
        cv[cvLength] <- cv[cvLength] * n / tailLength          
    }
    cv
}

# correlation-based distance of a set of already-centered Z scores
pearson.dist <- function (m) { # m is a matrix
  m <- m / sqrt (rowSums (m^2))
  m <-  tcrossprod (m)
  m <- as.dist(m)
  0.5 - m / 2
}
pearson.matrix <- function (m) { # m is a matrix
  m <- m / sqrt (rowSums (m^2))
  m <-  tcrossprod (m)
  0.5 - m / 2 # return the full matrix, not a dist formatted object
}

# function to create a normalized pseudo-color image
normZ.zLimit <- 3 
normZ <- function(z, edge=NULL){
  # scale Z from -3:3 into -1:1
  z <- ifelse(z < -normZ.zLimit, -normZ.zLimit, ifelse(z > normZ.zLimit, normZ.zLimit, z)) / normZ.zLimit
  # force extreme gain and loss value to z = 1 in final output
  if(!is.null(edge)){
    if(edge == 1) z <- apply(z, 2, pmin, 0) # gains
    else  z <- apply(z, 2, pmax, 0) # losses
  }
  # use abs values to reflect probabiliy of difference, and flip for proper plotting
  val <- 1 - abs(z)
  ifelse(is.na(val), 1, val) # unclear why this sometimes happens on "all" plot
}

# vertically expand a heatmap to multiple pixels per cell (NOT interpolation, just interleaved repetition)
expandImgV <- function(img, n){
  imgD <- dim(img)
  lrg <- as.cimg(array(0, dim=c(imgD[1], imgD[2]*n, 1, 3)))
  for(x in 0:(n-1)) lrg[,1:imgD[2] * n - x,,] <- img[,,,]    
  lrg
}

expandImg <- function(img, h=1, v=1){
  if(h > 1){
    imgD <- dim(img)
    lrg <- as.cimg(array(0, dim=c(imgD[1]*h, imgD[2], 1, 3)))   
    for(x in 0:(h-1)) lrg[1:imgD[1] * h - x,,,] <- img[,,,]
    img <- lrg
  }
  if(v > 1){
    imgD <- dim(img)
    lrg <- as.cimg(array(0, dim=c(imgD[1], imgD[2]*v, 1, 3)))
    for(x in 0:(v-1)) lrg[,1:imgD[2] * v - x,,] <- img[,,,]
    img <- lrg
  }
  img
}

# code development utility for exploring object sizes and environments
printObjectSizes_ <- function(envir, name, minSizeMb=1){
    message(name)
    OS <- round(sort( sapply(ls(envir=envir),function(x){
        object.size(get(x, envir=envir))
    })) / 1e6, 3)
    print(OS[OS>minSizeMb])     
}
printObjectSizes <- function(sessionEnv, minSizeMb=1){
    message('object sizes in MB')
    printObjectSizes_(.GlobalEnv, '.GlobalEnv', minSizeMb)    
    printObjectSizes_(sessionEnv, 'sessionEnv', minSizeMb)
}
# code development utility for finding slow steps
verbose <- FALSE
reportProgress <- function(message){
    if(verbose) message(message)
}



## Hidden Markov Model based on Poisson emissions
#persistence <- 1 - 1e-6
#tp <- matrix(log((1-persistence)/(3-1)), 3, 3)
#for(i in 1:3) tp[i,i] <- log(persistence)
#HMM_viterbi <- function(cellI){
#
#  # 0. emission probabilities
#  maxCount <- round(2 * cL[,cellI,cLN$expG], 0)
#  ep <- sapply(c(cLN$expL, cLN$exp0, cLN$expG), function(exp){
#    log(dpois(pmin(maxCount, round(cL[,cellI,cLN$raw],0)), cL[,cellI,exp]))
#  })
#
#  # 1. initialization (observation t=1)
#  T         <- nrow(ep) # length of the sequence of observations
#  N         <- ncol(ep)
#  delta     <- log(matrix(0, nrow=T, ncol=N))
#  delta[1,] <- sapply(1:N, function(i) log(1/N) + ep[1,i])  
#  phi       <-     matrix(NA, nrow=T, ncol=N)  
#  
#  # 2. recursion;
#  # NB: these 'for' loops are faster than apply methods with array as implemented and given recursion restrictions
#  for (t in 2:T){
#    pt <- t - 1
#    for (j in 1:N){   # j = this hs
#      ep_j <- ep[t,j] 
#      for (i in 1:N){ # i = prev hs
#        delta_ <- delta[pt,i] + tp[i,j] + ep_j      
#        if(delta[t,j] < delta_){
#          delta[t,j] <- delta_
#          phi[pt,j]  <- i
#        }
#      }
#    }
#  }   
#
#  # 3. termination
#  prob <- -Inf
#  hsi  <- rep(-1, T)
#  for (j in 1:N){
#    if(prob < delta[T,j]){
#      prob <- delta[T,j]
#      hsi[T] <- j
#    }
#  }  
#  
#  # 4. reconstruction
#  for (t in (T-1):1) hsi[t] <- phi[t,hsi[t+1]]  
#  
#  # return
#  hsi - 2
#  # list(prob=prob, hsi=hsi) # overall likelihood and ordered set of most likely states  
#}

