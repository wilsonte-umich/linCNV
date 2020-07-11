
# load dependencies
library(imager)
library(abind)

#-------------------------------------------------------------------
# define available colors heat map colors
#-------------------------------------------------------------------
# define available colors
hues <- list( 
    red   = list(
        max = 0.8, # max intensity prevents excessive image brightness
        value = rgb(0.8, 0, 0, maxColorValue = 1),
        index = 1
    ),  
    green = list(
        max = 0.65,
        value = rgb(0, 0.65, 0, maxColorValue = 1),
        index = 2
    ), 
    blue  = list(
        max = 0.9,
        value = rgb(0, 0, 0.9, maxColorValue = 1),
        index = 3
    )
)
lowColor  = "blue" # consistent with prior work, e.g. blue = deletion
highColor = "red"
highlightColor = "green" # yields yellow for red+gree, cyan for blue+green

#-------------------------------------------------------------------
# save intensity data as png heat map using ImageR
# much faster than heat map plotting in R (e.g. by rect or other)
#-------------------------------------------------------------------
# convert an intensity matrix to an array of RGB colors
parseIntensity <- function(intensity, triangle, mirror){ # step 1: apply triangle, mirror and inversion
    intensity[is.na(intensity)] <- 0 # ImageR requires actual values
    if(triangle) {
        intensity <- as.matrix(imrotate(as.cimg(intensity), -45)) # -45 is the angle that works!
        nc <- ncol(intensity)
        intensity <- intensity[,ceiling(nc/2):nc]
    } else if(mirror) {
        intensity <- mirrorPyramid(intensity) # make rectangular view be symmetric
    }    
    t(intensity) # transpose for ImageR coordinate consistency
}
parseHue <- function(intensity, hue){ # step 2: create the colors array
    intensity <- 1 - intensity # for color, 0 is max intensity
    max <- hues[[hue]]$max
    hue_ <- max + (1 - max) * intensity
    switch(hue,
        "red"   = abind(hue_, intensity, intensity, along=3), # abind to merge matrices
        "green" = abind(intensity, hue_, intensity, along=3),
        "blue"  = abind(intensity, intensity, hue_, along=3)
    )    
}
finishAndPrintImage <- function(img, file, xScaleFactor, yScaleFactor, transpose, decorate, ...){
    #reportProgress('saveImg_png')
    if(transpose){
        tmp <- xScaleFactor
        xScaleFactor <- yScaleFactor
        yScaleFactor <- tmp
        img <- aperm(img, c(2,1,3))
    }
    cimg <- imrotate(suppressWarnings(as.cimg(img)), -90) # rotate so 1,1 is bottom left of rectangle, or left of triangle
    if(xScaleFactor != 1 | yScaleFactor != 1) { # expand the image when instructed
        cimg <- resize(cimg, size_x=-100*xScaleFactor, size_y=-100*yScaleFactor)
    }
    if(!is.null(decorate)) cimg <- decorate(cimg, xScaleFactor, yScaleFactor, ...) # decorate applied AFTER expansion
    save.image(cimg, file)
}
# save intensity matrix [0,1=high] as a single-color heat map png
# triangle = TRUE yields a rotated triangle corresponding to the upper triangle of identity matrix
# mirror = TRUE makes a rectangular matrix be symmetrix from upper triangle
saveHeatMap_one_color <- function(intensity, hue, file, ..., # ... passed to decorateFN
                                  triangle=FALSE, mirror=TRUE,
                                    xScaleFactor=1, yScaleFactor=1, transpose=FALSE,
                                    decorate=NULL){
    #reportProgress('saveHeatMap_one_color')   
    intensity <- parseIntensity(intensity, triangle, mirror)
    cimg <- suppressWarnings(as.cimg(parseHue(intensity, hue)))
    finishAndPrintImage(cimg, file, xScaleFactor, yScaleFactor, transpose, decorate, ...)
}
# same thing for two color intensity on a diverging axis (i.e. low color, high color)
# when input low intensity is >0, high should be 0, and vice versa
saveHeatMap_two_color <- function(intensityLow, intensityHigh, highlight, file, ..., # ... passed to decorateFN
                                  triangle=FALSE, mirror=TRUE,
                                    xScaleFactor=1, yScaleFactor=1, transpose=FALSE,
                                    decorate=NULL){
    #reportProgress('saveHeatMap_two_color')
    intensityLow  <- parseIntensity(intensityLow,  triangle, mirror)
    intensityHigh <- parseIntensity(intensityHigh, triangle, mirror)
    highlight     <- parseIntensity(highlight,     triangle, mirror)
    colorsLow     <- parseHue(intensityLow,  lowColor) 
    colorsHigh    <- parseHue(intensityHigh, highColor)
    for(i in 1:3){ # execute the merge of the two, mutually exclusive, colors
        colorsHigh[,,i] <- ifelse(intensityLow > 0, colorsLow[,,i], colorsHigh[,,i])
    }
    colorsHigh[,,hues[[highlightColor]]$index] <- pmax(highlight[,],
                                                       colorsHigh[,,hues[[highlightColor]]$index])
    finishAndPrintImage(colorsHigh, file, xScaleFactor, yScaleFactor, transpose, decorate, ...)
}
# more comlicated pseudo-coloring where low distribution = blue, middle = green, high = red
saveHeatMap_three_color <- function(intensityLow, intensityMid, intensityHigh,
                                    file, ..., # ... passed to decorateFN
                                    triangle=FALSE, mirror=TRUE,
                                    xScaleFactor=1, yScaleFactor=1, transpose=FALSE,
                                    decorate=NULL){
    #reportProgress('saveHeatMap_three_color')
    finishAndPrintImage(abind(
            parseIntensity(intensityHigh, triangle, mirror), # here, do not invert intensities
            parseIntensity(intensityMid,  triangle, mirror),
            parseIntensity(intensityLow,  triangle, mirror),
            along=3
    ), file, xScaleFactor, yScaleFactor, transpose, decorate, ...)
}

#-------------------------------------------------------------------
# convert Z scores to color intensities in range [0,1]
#-------------------------------------------------------------------
# ramp up color signal only when Z >= 0
zToIntensity_pos <- function(Z, minZ=1, maxZ=4){
    Z[Z<0] <- 0
    deltaZ <- pmin(pmax(Z, minZ), maxZ)
    (deltaZ - minZ) / (maxZ - minZ)  
}
# decay signal for negative Z, full intensity for Z > 0
zToIntensity_neg_decay <- function(Z, maxZ=3){
    fullOn <- Z >= 0
    deltaZ <- pmin(pmax(-Z, 0), maxZ)
    ifelse(fullOn, 1, 1 - deltaZ / maxZ)
}
# decay signal on both sides of peak, full intensity for Z == 0
zToIntensity_symmetric <- function(Z, maxZ=3){
    deltaZ <- pmin(abs(Z), maxZ)
    1 - deltaZ / maxZ
}
# 0.5 intensity for Z == 0, 0 for -maxZ, 1 for +maxZ
zToIntensity_midpoint <- function(Z, maxZ=3){
    pmax(pmin(Z / maxZ / 2 + 0.5, 1), 0)
}
# ramp up signal on both sides of peak, no intensity for Z == 0
zToIntensity_inverted <- function(Z, minZ=1, maxZ=4){
    deltaZ <- pmin(pmax(abs(Z), minZ), maxZ)
    (deltaZ - minZ) / (maxZ - minZ)  
}

