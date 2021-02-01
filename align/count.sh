#!/bin/bash

echo "counting reads before and after filtering and grouping"

# set tmp dir
TMP_DIR_WRK=$TMP_DIR_LARGE/$PIPELINE_NAME"_count_"$DATA_NAME
mkdir -p $TMP_DIR_WRK

# search all subdirectories of INPUT_DIR for cell FASTQ files
cd $INPUT_DIR
export CELL_DIRS=`ls -1 .` # e.g. Sample_2251-LK-1

# set counts file
export COUNTS_FILE=$DATA_GENOME_PREFIX.counts.txt
export RATES_FILE=$DATA_GENOME_PREFIX.rates.txt

# collect the DUP_RATE information per cell
slurp -s 100M samtools view $BAM_FILE |
cut -f 12,13 |
sed -e 's/XC:i://' -e 's/CB:Z://' |
sort --parallel=$N_CPU -T $TMP_DIR_WRK -S $MAX_SORT_RAM -k2,2 |
bedtools groupby -g 2 -c 1,1 -o sum,count > $COUNTS_FILE
checkPipe

# use R to add the initial read pair counts (before filtering and grouping)
Rscript $PIPELINE_DIR/align/count.R
checkPipe

# clean up
rm -rf $TMP_DIR_WRK

echo "done"
echo


