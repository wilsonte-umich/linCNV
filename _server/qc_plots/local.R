
#----------------------------------------------------------------------
# local.R launches the web tool on a desktop or laptop via localhost
#----------------------------------------------------------------------
# you just need to edit the proper file paths for your computer, below
#----------------------------------------------------------------------

# get the script path
actionsPath <- dirname(sys.frame(1)$ofile)

# set required environment variables
Sys.setenv(
    # you should not need to change these values
    MODE            = "local",
    PIPELINE_NAME   = "linCNV",
    LOCAL_PORT      = "10000",
    ACTIONS_PATH    = actionsPath,
    
    # edit the following specific paths for local mode
    DATA_PATH       = "Z:\\data\\scCNV",  # required
    GENOMES_DIR     = "Z:\\data\\genomes" # can be "", but you won't see gene annotations    
)

# call app.R
source(paste(actionsPath, 'app.R', sep="/"))

