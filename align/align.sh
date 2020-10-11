#!/bin/bash

echo "aligning cells to genome and sort per cell"

# temp directory for sorting
export TMP_DIR_WRK=$TMP_DIR_LARGE/$PIPELINE_NAME"_align_"$DATA_NAME
mkdir -p $TMP_DIR_WRK

# search all subdirectories of INPUT_DIR for cell FASTQ files
cd $INPUT_DIR
export CELL_DIRS=`ls -1 .`

# align one cell at a time, with parallelization on each cell
perl $PIPELINE_DIR/align/align.pl
checkPipe

# merge the final output; one file for all cells, sorted
echo "merging to single final bam file"
samtools merge -b $LOGS_PREFIX.align.bam.list --threads $N_CPU - | 
slurp -s 100M -o $BAM_FILE
checkPipe

# index the final bam file
echo "indexing the final bam file"
samtools index -@ $N_CPU $BAM_FILE
checkPipe

# clean up
rm -r $TMP_DIR_WRK

echo "done"
echo

