
rawAggregatePlot <- function(bin, cell, viewport_, dendogram, mergedBins){
    #binAdj <- viewport_$region$posBinData$sizes / bin$medianSizeAutosome
    #RPA <- cell$readsPerAllele[cell$Is_accepted] # same number of cells as zLayers    
    #groupCells <- viewport_$cells[dendogram$inGroup] # values are indices of zLayers   
    #y <- rowMeans(matrix(unlist(sapply(1:length(groupCells), function(i){
    #    CN <- dendogram$gd[,i] / RPA[groupCells[i]] / binAdj
    #    CN - viewport_$region$posBinData$CNInt
    #})), ncol=ncol(dendogram$gd)))

    #binAdj <- bin$sizes / bin$medianSizeAutosome
    #exp0_bins <- CNInt_bins * RPA_cell *  binAdj_bins
    # should only correcct for binAdj_bins when mixed CN types
    # when bin is large due to mappability, all cells in fact expect CNInt_bins * RPA_cell

    CNInt <- viewport_$region$posBinData$CNInt

    
    #CNInt <- mergedBins$CNInt[viewport_$region$bins]
    
    plot(1:length(CNInt), CNInt)    
}

