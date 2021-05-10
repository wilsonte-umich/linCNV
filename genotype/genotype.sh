#!/bin/bash

echo "find heterozygous SNPs from good and wavy (not aneuploid, anomalous or discarded) cells"

# collect information from 'bin' analysis and user's manual cell calling
export ACCEPTED_CELLS=`Rscript $PIPELINE_DIR/genotype/get_cells.R`
checkPipe

# determine the species and thus analysis chromosomes from the genome
if [[ "$GENOME" = hg* ]]; then
    export SPECIES="human"
    export MAX_CHROM=22
else
    export SPECIES="mouse"
    export MAX_CHROM=19
fi

# determine if input bam reheader action is needed to add "chr"
# script assumes that resource, exclusion and other files are "chr1" format
TEST_CHROM=`samtools view $BAM_FILE | head -n 1 | cut -f 3`
if [[ "$TEST_CHROM" = chr* ]]; then
    export CHROMS=`perl -e 'print join("\n", map { "chr$_" } 1..$ENV{MAX_CHROM})'`    
    export REHEADER="cat"    
else
    export CHROMS=`perl -e 'print join("\n", 1..$ENV{MAX_CHROM})'`    
    REHEADER_FILE=$GENOTYPE_PREFIX.reheader.sam
    samtools view -H $BAM_FILE | perl -ne 's/(\@SQ\tSN:)/$1chr/g; print' > $REHEADER_FILE
    checkPipe
    export REHEADER="samtools reheader $REHEADER_FILE -"      
fi

# determine if we need to force VCF output to "1" instead of "chr1" for imputation server
if [[ "$GENOME" = "hg19" && "$WILL_IMPUTE" == "TRUE" ]]; then
    export STRIP_CHR="s/chr//"
else
    export STRIP_CHR=""
fi

# code to run bcftools mpileup + call per chromosome
function runChrom {

    # filter against low quality read pairs and bad genomic regions
    slurp samtools view -b -q $MIN_MAPQ -f 3 -F 3852 $BAM_FILE $1 |
    $REHEADER | 
    bedtools intersect -v -a stdin -b $EXCLUSIONS_FILE | # EXCLUSIONS_FILE was created by bin.sh
    
    # filter for reads from good+wavy cells only
    samtools view -h - |
    perl $PIPELINE_DIR/genotype/filter_reads.pl |   

    # execute variant calling; keep homozgyous variants to assist with phasing
    bcftools mpileup -Ou -f $GENOME_DIR/$GENOME.fa -q 0 -Q 0 - | # -Q,--min-BQ; at least some CellRanger bam files have only 'F:,' values in QUAL
    bcftools call -Ou -mv | 
    bcftools filter -s LowQual -e '%QUAL<20 || DP<10' | #  || DP>100
        
    # ensure that hg19 chromosomes are converted without "chr" prefix for submission to imputation server
    # all other genomes retain the "chr" prefix (even if the input bam didn't have the prefix)
    perl -ne "$STRIP_CHR; print" | 

    # compress and finish; inplicitly already sorted since bam was sorted
    bgzip -c |
    slurp -s 1M -o $GENOTYPE_PREFIX.$1.vcf.gz
    # NB: cannot do checkPipe in this exported function
    
    # requery the vcf files to establish base values for each read covering SNPs from cells with large deletions

    # YES, do filter this for dups!
    slurp samtools view -b $BAM_FILE $1 |
    perl $PIPELINE_DIR/genotype/filter_cells.pl |   
    
    $ST mpileup --positions <(filter and query vcf) --output-extra XC - |
    cut -f 1,2,5,7 |
    sed 's/\^.//g'
    
#/garage/wilsonte_lab/bin/samtools/samtools-1.12/bin/samtools    
#--positions
#BED or position list file containing a list of regions or sites where pileup or BCF should be generated. Position list files contain two columns (chromosome and position) and start counting from 1. BED files contain at least 3 columns (chromosome, start and end position) and are 0-based half-open.
#While it is possible to mix both position-list and BED coordinates in the same file, this is strongly ill advised due to the differing coordinate systems. [null]
#1       9994    tttt    CACATTTGTGTCCGAC-1,CACATTTGTGTCCGAC-1,CACATTTGTGTCCGAC-1,CACATTTGTGTCCGAC-1
#1       9995    ggggGg  CACATTTGTGTCCGAC-1,CACATTTGTGTCCGAC-1,CACATTTGTGTCCGAC-1,CACATTTGTGTCCGAC-1,CCGTACTAGTTTCCTT-1,TCTGAGATCCCAAAGT-1
}
export -f runChrom

# run bcftools in parallel per chromosome
echo "$CHROMS" | parallel -j $N_CPU runChrom {}

# finish up 
echo "done"
echo




#(base) [wilsonte@wilsonte-n1 hg19]$ zcat HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz | head
##CHROM  POS     ID      REF     ALT     AC      AN      AF      AC_EXCLUDING_1000G      AN_EXCLUDING_1000G      AF_EXCLUDING_1000G      AA
#1       13380   rs571093408     C       G       5       64940   7.69941e-05     0       59950   0       .
#1       16071   rs541172944     G       A       8       64940   0.000123191     0       59950   0       G
#1       16141   rs529651976     C       T       9       64940   0.000138589     4       59950   6.67223e-05     C
#1       16280   .       T       C       43      64940   0.00066215      2       59950   3.33611e-05     T
#1       49298   rs200943160     T       C       41571   64940   0.640145        38155   59950   0.636447        T
#1       54353   rs140052487     C       A       59      64940   0.000908531     22      59950   0.000366972     C
#1       54564   rs558796213     G       T       17      64940   0.00026178      3       59950   5.00417e-05     G
#1       54591   rs561234294     A       G       14      64940   0.000215584     0       59950   0       A
#1       54676   rs2462492       C       T       24175   64940   0.372267        23218   59950   0.387289        T



