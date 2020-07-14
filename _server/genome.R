
#----------------------------------------------------------------------
# genome.R handles genome-level metadata and annotations
#----------------------------------------------------------------------

# set the standard annotations
annotations <- list(
    hg38 = "gencode_26",
    mm10 = "gencode_m12"
)

# set genome strand colors
topStrandColor <- rgb(0,200,0,maxColorValue=255)
botStrandColor <- rgb(200,0,0,maxColorValue=255)

# set genome file paths
getGenomeDir <- function(genome){
    paste(serverEnv$GENOMES_DIR, genome, sep="/")
}
getAnnotationFile <- function(genome){
    paste(getGenomeDir(genome), "/", genome, ".",  annotations[[genome]], ".map.genes.bed.bgz", sep="")
}

# load the reference genome annotation set
# stored in shared memory in global environment for all user sessions
genes <- list()
loadAnnotation <- function(genome){    
    if(!is.null(genes[[genome]])) return()
    genomeFile <- getAnnotationFile(genome)
    if(!file.exists(genomeFile)) return()
    g <- read.table(genomeFile, header=FALSE, sep="\t", stringsAsFactors=FALSE)    
    colnames(g) <- c("chrom", "start", "end", "name", "score", "strand")    
    g$name  <- sapply(strsplit(g$name, "/"), function(x) x[2])
    g$NAME  <- toupper(g$name)
    g$xname <- g$start + (g$end - g$start)/2
    g$yname <- ifelse(g$strand=="+", 3, 0) + 0.5
    g$ybot  <- ifelse(g$strand=="+", 2, 1)
    g$color <- ifelse(g$strand=="+", topStrandColor, botStrandColor)
    genes[[genome]] <<- g
}


#chr1    11868   14409   ENSG00000223972.5/DDX11L1       0       +
#chr1    14403   29570   ENSG00000227232.5/WASH7P        0       -
#chr1    17368   17436   ENSG00000278267.1/MIR6859-1     0       -
#chr1    29553   31109   ENSG00000243485.5/MIR1302-2HG   0       +
#chr1    30365   30503   ENSG00000284332.1/MIR1302-2     0       +
#chr1    34553   36081   ENSG00000237613.2/FAM138A       0       -
#chr1    52472   53312   ENSG00000268020.3/OR4G4P        0       +
#chr1    62947   63887   ENSG00000240361.1/OR4G11P       0       +
#chr1    69090   70008   ENSG00000186092.4/OR4F5 0       +
#chr1    89294   133723  ENSG00000238009.6/RP11-34P13.7  0       -

#  chrom   start     end          name score strand   xname yname ybot   color
#1  chr1 3073252 3074322 4933401J01Rik     0      + 3073787     3    2 #00C800
#2  chr1 3102015 3102125       Gm26206     0      + 3102070     3    2 #00C800
#3  chr1 3205900 3671498          Xkr4     0      - 3438699     0    1 #C80000
#4  chr1 3252756 3253236       Gm18956     0      + 3252996     3    2 #00C800
#5  chr1 3365730 3368549       Gm37180     0      - 3367140     0    1 #C80000
#6  chr1 3375555 3377788       Gm37363     0      - 3376672     0    1 #C80000

#-rw-rw-r-- 1 bmagnuso wilsonte_lab       4500 Dec  9  2016 chromInfo.txt
#-rw-rw-r-- 1 wilsonte wilsonte_lab      42999 Mar  6  2014 gap.txt
#-rw-rw-r-- 1 wilsonte wilsonte_lab   37607127 Aug 21  2017 gencode.v26.annotation.gtf.gz
#-rw-rw-r-- 1 wilsonte wilsonte_lab   23333604 Aug  6  2017 gencode.v26.basic.annotation.gtf.gz
#-rw-rw-r-- 1 wilsonte wilsonte_lab      38196 Mar 29  2017 hg38.bands.bed
#-rw-r--r-- 1 wilsonte wilsonte_lab        389 Mar 18  2017 hg38.chromSizes.standard.txt
#-rw-r--r-- 1 wilsonte wilsonte_lab        428 Mar 18  2017 hg38.chromSizes.txt
#-rw-rw-r-- 1 wilsonte wilsonte_lab 3161924569 May  5  2014 hg38.fa
#-rw-rw-r-- 1 wilsonte wilsonte_lab      21294 Mar 29  2017 hg38.gap.bed
#-rw-rw-r-- 1 wilsonte wilsonte_lab    7542890 Aug  6  2017 hg38.gencode_26.blocks.stranded.bed.bgz
#-rw-rw-r-- 1 wilsonte wilsonte_lab     140445 Aug  6  2017 hg38.gencode_26.blocks.stranded.bed.bgz.tbi
#-rw-rw-r-- 1 wilsonte wilsonte_lab    7569686 Aug  6  2017 hg38.gencode_26.blocks.unstranded.bed.bgz
#-rw-rw-r-- 1 wilsonte wilsonte_lab     247817 Aug  6  2017 hg38.gencode_26.blocks.unstranded.bed.bgz.tbi
#-r--r--r-- 1 saurabha wilsonte_lab   37607127 Aug 19  2017 hg38.gencode_26.gtf.gz
#-rw-rw-r-- 1 wilsonte wilsonte_lab    4079247 Aug  6  2017 hg38.gencode_26.map.exons.bed.bgz
#-rw-rw-r-- 1 wilsonte wilsonte_lab     599074 Aug  6  2017 hg38.gencode_26.map.exons.bed.bgz.tbi
#-rw-rw-r-- 1 wilsonte wilsonte_lab    1081757 Aug  6  2017 hg38.gencode_26.map.genes.bed.bgz
#-rw-rw-r-- 1 wilsonte wilsonte_lab     172411 Aug  6  2017 hg38.gencode_26.map.genes.bed.bgz.tbi
#-rw-rw-r-- 1 wilsonte wilsonte_lab    9809405 Aug  6  2017 hg38.gencode_26.map.stranded.bed.bgz
#-rw-rw-r-- 1 wilsonte wilsonte_lab     308433 Aug  6  2017 hg38.gencode_26.map.stranded.bed.bgz.tbi
#-rw-rw-r-- 1 wilsonte wilsonte_lab    7675607 Aug  6  2017 hg38.gencode_26.map.transcripts.bed.bgz
#-rw-rw-r-- 1 wilsonte wilsonte_lab     189475 Aug  6  2017 hg38.gencode_26.map.transcripts.bed.bgz.tbi
#-rw-rw-r-- 1 wilsonte wilsonte_lab    6995198 Aug  6  2017 hg38.gencode_26.map.unstranded.bed.bgz
#-rw-rw-r-- 1 wilsonte wilsonte_lab     299981 Aug  6  2017 hg38.gencode_26.map.unstranded.bed.bgz.tbi
#-rw-rw-r-- 1 wilsonte wilsonte_lab      45544 Apr  2  2017 hg38.ideogram.bed
#-rw-rw-r-- 1 wilsonte wilsonte_lab       9696 Apr  2  2017 hg38.ideogram.bed.bgz
#-rw-rw-r-- 1 wilsonte wilsonte_lab      15643 Apr  2  2017 hg38.ideogram.bed.bgz.tbi


#-rw-rw-r-- 26 wilsonte wilsonte_lab       2791 Jan 18  2014 chromInfo.txt
#-rw-rw-r-- 26 wilsonte wilsonte_lab        598 Mar  7  2012 chromInfo.txt.gz
#-rw-rw-r-- 26 wilsonte wilsonte_lab      33467 Jan 18  2014 gap.txt
#-rw-rw-r-- 26 wilsonte wilsonte_lab   37896812 May  4  2014 genomicSuperDups.txt.gz
#-rw-rw-r-- 26 wilsonte wilsonte_lab        331 Jan 18  2014 mm10.chromSizes.chrM.txt
#-rw-r--r-- 26 wilsonte wilsonte_lab        341 Mar 18  2017 mm10.chromSizes.standard.txt
#-rw-r--r-- 26 wilsonte wilsonte_lab        368 Mar 18  2017 mm10.chromSizes.txt
#-rw-rw-r-- 26 wilsonte wilsonte_lab 2780048573 Jan 18  2014 mm10.fa
#-rw-rw-r-- 26 wilsonte wilsonte_lab      22639 Jan 12  2017 mm10.gap.bed
#-rw-r--r-- 26 wilsonte wilsonte_lab       7122 Jan  3  2017 mm10.gap.bed.bgz
#-rw-r--r-- 26 wilsonte wilsonte_lab       9131 Jan  3  2017 mm10.gap.bed.bgz.tbi
#-rw-rw-r-- 26 wilsonte wilsonte_lab    4505109 Mar 23  2017 mm10.gencode_m12.blocks.stranded.bed.bgz
#-rw-rw-r-- 26 wilsonte wilsonte_lab     115298 Mar 23  2017 mm10.gencode_m12.blocks.stranded.bed.bgz.tbi
#-rw-rw-r-- 26 wilsonte wilsonte_lab    4495438 Mar 23  2017 mm10.gencode_m12.blocks.unstranded.bed.bgz
#-rw-rw-r-- 26 wilsonte wilsonte_lab     199863 Mar 23  2017 mm10.gencode_m12.blocks.unstranded.bed.bgz.tbi
#-rw-rw-r-- 26 wilsonte wilsonte_lab   24160787 Mar 14  2017 mm10.gencode_m12.gtf.gz
#-rw-rw-r-- 26 wilsonte wilsonte_lab    3487376 Mar 26  2017 mm10.gencode_m12.map.exons.bed.bgz
#-rw-rw-r-- 26 wilsonte wilsonte_lab     495029 Mar 26  2017 mm10.gencode_m12.map.exons.bed.bgz.tbi
#-rw-rw-r-- 26 wilsonte wilsonte_lab     911391 Mar 26  2017 mm10.gencode_m12.map.genes.bed.bgz
#-rw-rw-r-- 26 wilsonte wilsonte_lab     160556 Mar 26  2017 mm10.gencode_m12.map.genes.bed.bgz.tbi
#-rw-rw-r-- 26 wilsonte wilsonte_lab    8528970 Mar 26  2017 mm10.gencode_m12.map.stranded.bed.bgz
#-rw-rw-r-- 26 wilsonte wilsonte_lab     275195 Mar 26  2017 mm10.gencode_m12.map.stranded.bed.bgz.tbi
#-rw-rw-r-- 26 wilsonte wilsonte_lab    5278103 Mar 26  2017 mm10.gencode_m12.map.transcripts.bed.bgz
#-rw-rw-r-- 26 wilsonte wilsonte_lab     180855 Mar 26  2017 mm10.gencode_m12.map.transcripts.bed.bgz.tbi
#-rw-rw-r-- 26 wilsonte wilsonte_lab    6061369 Mar 26  2017 mm10.gencode_m12.map.unstranded.bed.bgz
#-rw-rw-r-- 26 wilsonte wilsonte_lab     265817 Mar 26  2017 mm10.gencode_m12.map.unstranded.bed.bgz.tbi
#-r--r--r-- 26 saurabha wilsonte_lab   24953748 Aug 19  2017 mm10.gencode_m14.gtf.gz
#-rw-rw-r-- 26 wilsonte wilsonte_lab   15214014 Jan  7  2017 mm10.genomicSuperDups.bed
#-rw-r--r-- 26 wilsonte wilsonte_lab    1026771 Jan  7  2017 mm10.genomicSuperDups.bed.bgz
#-rw-r--r-- 26 wilsonte wilsonte_lab      98778 Jan  7  2017 mm10.genomicSuperDups.bed.bgz.tbi
#-rw-rw-r-- 26 wilsonte wilsonte_lab      15276 Apr  2  2017 mm10.ideogram.bed
#-rw-rw-r-- 26 wilsonte wilsonte_lab       4595 Apr  2  2017 mm10.ideogram.bed.bgz
#-rw-rw-r-- 26 wilsonte wilsonte_lab       6986 Apr  2  2017 mm10.ideogram.bed.bgz.tbi
#-rw-rw-r-- 26 wilsonte wilsonte_lab       5911 Jan 18  2014 mm10.stats.csv
#-rw-rw-r-- 26 wilsonte wilsonte_lab       4305 Nov 25  2016 mm10.svtools.exclude.CN_10.bed
#-rw-rw-r-- 26 wilsonte wilsonte_lab       8085 Nov 25  2016 mm10.svtools.exclude.CN_5.bed
#-rw-rw-r-- 26 wilsonte wilsonte_lab    1120463 Jan 18  2014 mm10.UCSC.genes.bed
#-rw-rw-r-- 26 wilsonte wilsonte_lab    1156905 Jan 18  2014 mm10.UCSC.map.genes.bed
#-rw-rw-r-- 26 wilsonte wilsonte_lab   48542799 Jan 18  2014 mm10.UCSC.map.stranded.bed
#-rw-rw-r-- 26 wilsonte wilsonte_lab   21481487 Jan 18  2014 mm10.UCSC.map.unstranded.bed
#-rw-rw-r-- 26 wilsonte wilsonte_lab   19658277 Jan 18  2014 mm10.UCSC.transcripts.bed
#-rw-rw-r-- 26    10171 wilsonte_lab  455496876 Apr  1  2014 Mus_musculus.GRCm38.75.gtf

