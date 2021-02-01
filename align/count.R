
# load resources
library(jsonlite)
env <- as.list(Sys.getenv())

# load counts data
d <- read.table(env$COUNTS_FILE, sep="\t", header=FALSE, stringsAsFactors=FALSE)
names(d) <- c('sampleName','filtered','grouped')
d$libraryName <- gsub('Sample_', '', d$sampleName)
d$sampleN <- as.integer(sapply(strsplit(d$libraryName, '-'), function(v) v[length(v)]))

# collect read pair counts from fastp
d$readPairs <- sapply(d$sampleName, function(sampleName){
    json <- paste(env$LOGS_PREFIX, sampleName, 'fastp', 'json', sep=".")
    json <- read_json(json)
    json$summary$before_filtering$total_reads / 2
})

# add a few read pair rates
d$alignRate <- d$filtered/d$readPairs
d$dupRate   <- 1 - d$grouped/d$filtered

# load and add source cell labels
manifest <- read.csv(env$MANIFEST_FILE, header=TRUE, stringsAsFactors=FALSE)
manifest <- manifest[manifest$Lane==1,3:4]
names(manifest) <- c('libraryName','cellName')
d <- merge(d, manifest, by='libraryName', all.x=TRUE, all.y=FALSE)
d <- d[order(d$sampleN),
       c('sampleN','sampleName','libraryName','cellName','readPairs','filtered','grouped','alignRate','dupRate')]

# write out results
write.table(d, env$RATES_FILE, quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

