
#----------------------------------------------------------------------
# constant values
#----------------------------------------------------------------------

# define cell types and acceptance patterns
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
acceptedButtons   <- acceptedCellTypes
unacceptedButtons <- cellTypeCodes['Hyperseg']:cellTypeCodes['Discard']
cellTypeColors <- c(
    'grey50',
    'green2',
    'blue',    
    'red'
    #,    
    #'darkorange'     
)

