
#----------------------------------------------------------------------
# load and save data for cell marking
#----------------------------------------------------------------------

# working variables
getSampleKey <- function(input) paste(input$project, input$sample, sep=", ")
nllSample <- "-, -"
crrSample <- nllSample
marks <- NULL

# load a sample's marks, called by plot_data function
loadSample_marks <- function(input){
    if(input$project != '-' & input$sample != '-'){
        if(crrSample != getSampleKey(input)){
            reportProgress(paste('loadSample_marks', input$sample))
            crrSample <<- nllSample         
            marks <<- NULL
            rDataFile <- getMarksFile(input)
            if(length(rDataFile) > 0) {
                load(rDataFile)
                marks[marks<=0] <- cellTypeCodes['Unmarked']
                crrSample <<- getSampleKey(input)                
                marks <<- marks
            }              
        }
    } else {
        crrSample <<- nllSample         
        marks <<- NULL
    }
}

# save the cell marks
saveMarks <- function(input){
    reportProgress('saveMarks')
    if(input$project != '-' & input$sample != '-' & !is.null(marks) & crrSample != nllSample){
        rDataFile <- getMarksFile(input)
        save(marks, file = rDataFile)
    }
}

