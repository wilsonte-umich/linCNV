
#----------------------------------------------------------------
# code development utility for reporting progress to find bugs and slow steps
#----------------------------------------------------------------
verbose <- FALSE
reportProgress <- function(message){
    if(verbose) message(message)
}

#----------------------------------------------------------------
# code development utility for exploring object sizes and environments
#----------------------------------------------------------------
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

