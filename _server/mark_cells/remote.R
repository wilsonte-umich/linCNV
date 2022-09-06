
#----------------------------------------------------------------------
# remote.R launches the web tool on a remote server accessed via SSH
#----------------------------------------------------------------------
# you just need to edit the proper file paths for your server, below
#----------------------------------------------------------------------

# get the script path
actionsPath <- file.path(Sys.getenv("SERVER_DIR"), Sys.getenv("APP"))

# set required environment variables
Sys.setenv(
    # you should not need to change these values
    MODE            = "remote",
    PIPELINE_NAME   = "linCNV",
    # REMOTE_PORT     = "11000",
    ACTIONS_PATH    = actionsPath
    # ,
    
    # # edit the following specific paths for local mode
    # DATA_PATH       = "Z:\\data\\scCNV",  # required
    # GENOMES_DIR     = "" # can be "", but you won't see gene annotations    
)

# call app.R
source(paste(actionsPath, 'app.R', sep="/"))

