
# simple script to extract the accepted cell list from the user-modified RData file

# get passed arguments
env <- as.list(Sys.getenv())

# load common functions and parameters
source(paste(env$PIPELINE_DIR, "common", "constants.R", sep="/"))

# load cell marks
rDataFile <- paste(env$BIN_PREFIX, "cell_marks", "RData", sep=".")
load(rDataFile)

# get list of accepted cells based on user marks
isAccepted_get_cells <- marks %in% normalCellTypes
N_accepted <- sum(isAccepted_get_cells)
if(N_accepted == 0) stop('No cells marked as accepted! Did you manually mark cells in web interface?')

# return list to caller
# output values are Cell Ranger style 0-reference cell_ids
writeLines(paste(which(isAccepted_get_cells) - 1, collapse=" ")) 

