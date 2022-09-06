
#----------------------------------------------------------------------
# app.R assembles and launches the Shiny server application
#----------------------------------------------------------------------

# load dependencies
library(shiny)
#library(DT) # data table add-on to shiny
#library(writexl) # too allow download of data tables
#library(Rtsne)

# get or set required variables
serverEnv <- as.list(Sys.getenv()) # uncomment if script is called by a wrapper script

# collect the server mode
if(is.null(serverEnv$MODE)) serverEnv$MODE = 'dev' # default to developer mode (inline execution)
isDev <- serverEnv$MODE == 'dev'
urlSuffix <- if(isDev) 'dev' else ''
port <- as.integer(switch(serverEnv$MODE, # get and check the port
    dev   = serverEnv$DEV_PORT,
    prod  = serverEnv$PROD_PORT,
    local = serverEnv$LOCAL_PORT,
    remote = serverEnv$REMOTE_PORT
))
if(is.null(port)) stop(paste("port missing for mode:", serverEnv$MODE))
host <- switch(serverEnv$MODE, # get and check the host
    dev   = "0.0.0.0",
    prod  = "0.0.0.0",
    local = "127.0.0.1",
    remote = "127.0.0.1"
)
if(is.null(host)) stop(paste("host missing for mode:", serverEnv$MODE))
for (var in c('PIPELINE_NAME','ACTIONS_PATH','DATA_PATH')){ # check required paths
    if(is.null(serverEnv[[var]])) stop(paste("missing variable:", var))
}

# load data resources and server scripts on server initialization
#sampleData <- list() # existing sample data does NOT reload, but new samples will
isInitialized <- FALSE
initalizeApp <- function(){
    if(isDev | !isInitialized){ # developer mode always reloads new code
        for(script in c(
            "../constants.R",   # ../<scipts> pattern >> shared by multiple server apps used by a pipeline
            "../data_sources.R",
            "utilities.R",      # non-prefixed scripts >> specific to this server app
            "genome.R",  
            "viewport.R",
            "temp.R"
        )) source(paste(serverEnv$ACTIONS_PATH, script, sep="/"))
    }
    setProjects()
    setSamples()
    if(isDev | !isInitialized){ # order of sourced scripts is important
        for(script in c(
            "ui.R",
            "server.R"            
        )) source(paste(serverEnv$ACTIONS_PATH, script, sep="/"))
    }
    isInitialized <<- TRUE
}
initalizeApp()

# how to load app in a browser (configured externally)
if(!is.null(serverEnv$URL)) message(paste("\n", serverEnv$URL, urlSuffix, "/", sep=""))
if(serverEnv$MODE == 'local') message(paste("\n", host, ":", port, sep=""))

# start the server with public access
runApp(
    shinyApp(ui=ui, server=server),
    host=host,
    port=port
)

