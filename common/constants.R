
#----------------------------------------------------------------
# cell types
#----------------------------------------------------------------
cellTypes <- c(
    # accepted cell types
    "__Good__",
    "__Wavy__",    
    "Aneuploid",    
    "Anomalous",  
    # rejected cell types
    "Hyperseg",
    "Doublet",
    "Replicating",
    "Discard",
    "Unmarked"
)
nCellTypes <- length(cellTypes)
cellTypeCodes <- 1:nCellTypes
names(cellTypeCodes) <- cellTypes
acceptedCellTypes   <- cellTypeCodes['__Good__']:cellTypeCodes['Anomalous']
unacceptedCellTypes <- cellTypeCodes['Hyperseg']:cellTypeCodes['Unmarked']

