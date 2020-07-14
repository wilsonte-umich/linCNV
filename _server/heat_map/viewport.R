
#----------------------------------------------------------------------
# viewport.R controls the x-axis coordinates of the display window
#----------------------------------------------------------------------

# current viewport x-axis range
getRegion <- function(input){
    list(
        minPos = as.integer(input$start),
        maxPos = as.integer(input$end)        
    )
}
vpWidth <- function(region=NULL, input=NULL) {
    if(is.null(region)) region <- getRegion(input)
    region$maxPos - region$minPos + 1
}

# clear position values when changing chromosomes
clearPos <- function(session){
    if(!jumpStatus$inProgress){
        updateTextInput(session, 'start', value = "")
        updateTextInput(session, 'end',   value = "")        
    }
}

# move window left and right
jumpWindow <- function(input, session, multiplier){
    region <- getRegion(input)
    vpWidth <- vpWidth(region)
    min <- region$minPos + multiplier * vpWidth
    max <- region$maxPos + multiplier * vpWidth
    if(min >= 1) {
        updateTextInput(session, 'start', value = min)
        updateTextInput(session, 'end',   value = max)
    } 
}

# zoom window in and out
zoomFactor <- 3
getDw <- function(inWidth, zoomFactor=3, exponent=-1){
    newWidth <- round(inWidth / zoomFactor ** exponent, 0)
    (newWidth - inWidth) / 2  
}
zoomWindow <- function(input, session, exponent){
    region <- getRegion(input)
    vpWidth <- vpWidth(region)
    dw <- getDw(vpWidth, zoomFactor, exponent)
    min <- region$minPos - dw
    max <- region$maxPos + dw
    if(min >= 1) {
        updateTextInput(session, 'start', value = min)
        updateTextInput(session, 'end',   value = max)
    } 
}

# jump window to coordinates requested by user chr3:50000000-75000000   chr3 50000000 75000000 FHIT
jumpStatus <- reactiveValues(inProgress=FALSE) # prevents reset to whole chromosome on chrom change
jumpTo <- function(input, session, genome){
    if(input$jumpTo != ""){
        if(grepl(":", input$jumpTo)){ # samtools style coordinates, chr:start-end
            tmp <- strsplit(input$jumpTo, ":")[[1]]
            jumpTo_(session, c(tmp[1], strsplit(tmp[2], "-")[[1]]))
        } else if(grepl("\\s+", input$jumpTo)){ # space-delimited, chr start end
            jumpTo_(session, strsplit(input$jumpTo, "\\s+")[[1]])  
        } else if(toupper(input$jumpTo) %in% genes[[genome]]$NAME) { # a gene name
            gene <- genes[[genome]][genes[[genome]]$NAME==toupper(input$jumpTo),][1,]
            dw <- getDw(gene$end-gene$start)
            jumpTo_(session, c(gene$chrom, gene$start-dw, gene$end+dw))
        }
    }
}
jumpTo_ <- function(session, jumpTo){
    if(length(jumpTo) == 3){
        jumpStatus$inProgress <<- TRUE
        updateSelectInput(session, 'chrom', selected=jumpTo[1])
        updateTextInput(session,   'start', value = jumpTo[2])
        updateTextInput(session,   'end',   value = jumpTo[3])
        updateTextInput(session,   'jumpTo',value = "")
    }    
}

