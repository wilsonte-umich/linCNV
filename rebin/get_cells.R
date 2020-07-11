
# simple script to extract the accepted cell list from the user-modified RData file

# get passed arguments
env <- as.list(Sys.getenv())
command <- 'rebin'

# load common functions and parameters
source(paste(env$PIPELINE_DIR, "common", "constants.R", sep="/"))

# load cell marks
rDataFile <- paste(env$BIN_PREFIX, "cell_marks", "RData", sep=".")
load(rDataFile)

# make a hard copy of the marks in force at time of rebinning
rDataFile <- sub('\\.bin\\.', '\\.rebin\\.', rDataFile)
save(marks, file=rDataFile)

# get list of accepted cells based on user marks
isAccepted_get_cells <- marks %in% acceptedCellTypes
N_accepted <- sum(isAccepted_get_cells)
if(N_accepted == 0) stop('No cells marked as accepted! Did you manually mark cells in web interface?')

# apply one last filter to disregard untrustworthy cells with low read depth
rDataFile <- paste(env$BIN_PREFIX, "layers", "RData", sep=".")
load(rDataFile)
env <- as.list(Sys.getenv())
highDepthCellIs <- which(cell$readsPerAllele_working >= as.numeric(env$MIN_ALLELE_DEPTH))
isAccepted_get_cells <- isAccepted_get_cells & cell$Is %in% highDepthCellIs
N_accepted <- sum(isAccepted_get_cells)
if(N_accepted == 0) stop('All marked cells were rejected based on read depth; please lower MIN_ALLELE_DEPTH.')

# return list to caller
# output values are Cell Ranger style 0-reference cell_ids
writeLines(paste(which(isAccepted_get_cells) - 1, collapse=" ")) 

