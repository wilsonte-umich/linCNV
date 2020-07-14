
#----------------------------------------------------------------------
# data_sources.R has functions to load scCNV sample data
# sourced locally by 'server', so variables are scoped to each user session
#----------------------------------------------------------------------

# load genome information for a plottable sample
genomeInfo <- list(name="")
chromI <- 1 # since has an X in front of the column name...
setSampleGenome <- function(input, output){
    genomeInfo <<- list(name = "")
    if(input$sample != "-"){
        load(getChromsFile(input)) # has genome and chroms objects
        loadAnnotation(genome)
        genomeInfo <<- list(
            name   = genome,
            chroms = chroms
        )   
    }
    output$genome <- renderText({
        genomeInfo$name })    
}

# load the bin x cell data for a sample + chrom (sc)
crrProject <- "NULL"
crrSample  <- "NULL"
zLayers <- list()   # incoming data object that changes with the sample and chromosome
sampleObjects <- c( # incoming data objects that change with the sample, but not the chromosome
    #'zLN',
    'cell',
    'bin',
    'mergedBins',
    'genome',
    'chroms',
    'ploidy'
)
clearSampleData <- function(){
    crrProject <<- "NULL"
    crrSample  <<- "NULL"
    zLayers <<- list()    
    suppressWarnings(rm(list=sampleObjects, envir=sessionEnv))
    heatMapColors <<- list()
}
sampleFileError <- function(file){
    message(paste('file load failed:', file))
    clearSampleData()
} 
loadSampleChromData <- function(caller=NULL){
    if(input$sample == "-" | input$chrom == "-") return( clearSampleData() ) # nothing to do
    if(crrProject != input$project | crrSample != input$sample) clearSampleData()
    crrProject <<- input$project
    crrSample  <<- input$sample
    if(!is.null(zLayers[[input$chrom]])) return() # already loaded this sample chrom
    
    #if(!is.null(caller)) message(caller)
    reportProgress(paste("loading data:", input$sample, input$chrom))

    # load just the sample cell marks
    rDataFile <- getMarksFile(input, 'rebin')
    if(length(rDataFile) == 0) return(sampleFileError('rebin marks'))
    load(rDataFile)
    assign('marks', marks, envir=sessionEnv)    

    # then load the entire set of analysis RData objects
    rDataFile <- getAnalyzeFile(input)
    if(length(rDataFile) == 0) return(sampleFileError(paste('analyze', input$sample, input$chrom)))
    load(rDataFile)
    bin$data$sizes <- bin$sizes # simplify downstream code
    bin$data$CNInt <- bin$modalCNsInt

    # keep required sample-specific objects on first sample load
    if(!exists('cell', envir=sessionEnv)){
        for(x in sampleObjects) assign(x, get(x), envir=sessionEnv)
        assign('marks_accepted', marks[cell$Is_accepted], envir=sessionEnv)
        assign('chromInfo',  list(
            starts = sapply(chroms, function(chrom){
                min(which(mergedBins$chrom==chrom))
            }),
            mids = sapply(chroms, function(chrom){
                mean(which(mergedBins$chrom==chrom))
            })            
        ), envir=sessionEnv)          
    }    
 
    # keep chomosome-specific objects
    zLayers[[input$chrom]] <<- zLayers
}


#-rw-rw-r-- 1 wilsonte wilsonte_lab         153 Aug 30 11:42 alarms_summary.txt
#-rw-rw-r-- 1 wilsonte wilsonte_lab   772817810 Aug 30 11:42 cnv_data.h5
#-rw-rw-r-- 1 wilsonte wilsonte_lab   132200535 Aug 30 11:43 dloupe.dloupe
#-rw-rw-r-- 1 wilsonte wilsonte_lab       71421 Aug 30 10:55 mappable_regions.bed
#-rw-rw-r-- 1 wilsonte wilsonte_lab    12111730 Aug 30 11:40 node_cnv_calls.bed
#-rw-rw-r-- 1 wilsonte wilsonte_lab   184532396 Aug 30 11:40 node_unmerged_cnv_calls.bed
#-rw-rw-r-- 1 wilsonte wilsonte_lab      103483 Aug 30 11:41 per_cell_summary_metrics.csv
#-rw-rw-r-- 1 wilsonte wilsonte_lab 72153790037 Aug 30 11:52 possorted_bam.bam
#-rw-rw-r-- 1 wilsonte wilsonte_lab     8176416 Aug 30 11:58 possorted_bam.bam.bai
#-rw-rw-r-- 1 wilsonte wilsonte_lab        1134 Aug 30 11:41 summary.csv
#-rw-rw-r-- 1 wilsonte wilsonte_lab     1705876 Aug 30 11:42 web_summary.html

#per_cell_summary_metrics.csv
#             barcode cell_id total_num_reads num_unmapped_reads
#1 AAACCTGCACCACACG-1       0         1401732               6661
#2 AAACCTGCAGCGATGA-1       1          960204               4016
#3 AAACCTGTCAGTGTGT-1       2         2041926               7735
#4 AAACGGGTCTCGAGCG-1       3         1212198               5405
#5 AAAGCAACACGCACGT-1       4         1426874               5681
#6 AAAGCAAGTCGGAAGT-1       5         1836188               6359
#  num_lowmapq_reads num_duplicate_reads num_mapped_dedup_reads
#1            222057              163321                1009693
#2            133983              113308                 708897
#3            288588              242671                1502932
#4            198413              135741                 872639
#5            240543              167098                1013552
#6            286227              215738                1327864
#  frac_mapped_duplicates effective_depth_of_coverage effective_reads_per_1Mbp
#1              0.1165137                  0.07004812                      370
#2              0.1180041                  0.04789969                      260
#3              0.1188442                  0.10246717                      551
#4              0.1119792                  0.05971592                      320
#5              0.1171077                  0.07141811                      372
#6              0.1174923                  0.09193981                      487
#    raw_mapd normalized_mapd raw_dimapd normalized_dimapd mean_ploidy
#1 0.11524713      0.11524713  1.0562024         1.0562024    1.949663
#2 0.12158889      0.12158889  0.9575596         0.9575596    1.944958
#3 0.09572153      0.09572153  1.0071066         1.0071066    1.944778
#4 0.11828912      0.11828912  1.0239070         1.0239070    1.944844
#5 0.11145977      0.11145977  1.0172727         1.0172727    1.946156
#6 0.09389352      0.09389352  0.9349878         0.9349878    1.942019
#  ploidy_confidence is_high_dimapd is_noisy est_cnv_resolution_mb
#1                -2              0        0             0.9989746
#2                -2              0        0             1.0991699
#3                14              0        0             0.7084961
#4                -2              0        0             1.0469238
#5                -2              0        0             0.9382813
#6                -2              0        0             0.6821777

