
#----------------------------------------------------------------------
# data_sources.R defines data paths and functions to load data
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# file paths
#----------------------------------------------------------------------
getDataSubPath   <- function(subFolder) {
    dataPath <- if(exists('env')){
        if(!is.null(env$DATA_PATH)) env$DATA_PATH else serverEnv$DATA_PATH
    } else {
        serverEnv$DATA_PATH
    }
    paste(dataPath, subFolder, sep="/")
}
getProjectFolder <- function(project) paste(getDataSubPath("projects"), project, sep="/")
getSampleFolder  <- function(project, sample) paste(getProjectFolder(project), sample, sep="/")
getDataFile <- function(project, sample, type, subType, extension){ # only return a file if it exists
    folder <- getSampleFolder(project, sample)
    file   <- paste(sample, "*", type, subType, extension, sep=".") # of whatever genome a sample used
    list.files(folder, file, full.names=TRUE)[1]
}
getDataFileName <- function(project, sample, type, subType, extension){ # returns the proper name of any file, no genome
    folder <- getSampleFolder(project, sample)
    file   <- paste(sample, type, subType, extension, sep=".")
    paste(folder, file, sep="/")
}
getMarksFile <- function(input, type='bin') {
    getDataFile(input$project, input$sample, type, 'cell_marks', 'RData')
}
getCellPlot <- function(input) {
    plotFolder <- paste(getSampleFolder(input$project, input$sample), "plots", "cells", sep="/")
    paste(plotFolder, "/", input$sample, ".cell.", sprintf("%05d", as.integer(input$cellNumber)), ".png", sep="")
}
#----------------------------------------------------------------------
getChromsFile <- function(input){
    getDataFile(input$project, input$sample, 'analyze', 'genome_info', 'RData')
}
getAnalyzeFile <- function(input){
    subType <- paste('layers', input$chrom, sep=".")
    getDataFile(input$project, input$sample, 'analyze', subType, 'RData')
}
getFindFile <- function(input){
    subType <- paste('find', input$chrom, sep=".")
    getDataFile(input$project, input$sample, 'analyze', subType, 'RData')
}
#----------------------------------------------------------------------
getQCPlot <- function(input, genome, type) {
    plotFolder <- paste(getSampleFolder(input$project, input$sample), "plots", sep="/")
    pngName <- paste(input$sample, genome, 'analyze', type, 'png', sep=".")
    paste(plotFolder, pngName, sep="/")
}

#----------------------------------------------------------------------
# initialize the data sources for the overall app, e.g. all samples
# stays in memory, shared by all user sessions
#----------------------------------------------------------------------
# all available projects
projects <- character()
setProjects <- function(){
    projects <<- c('-', list.dirs(getDataSubPath("projects"), full.names=FALSE, recursive=FALSE))
}
# all available samples
samples <- list()
setSamples <- function(){
    samples <<- lapply(projects, function(project){
        list.dirs(getProjectFolder(project), full.names=FALSE, recursive=FALSE)
    })
    names(samples) <<- projects
}

