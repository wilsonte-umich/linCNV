#!/bin/bash

echo "aligning cells to genome and merging to sorted bam"

# temp directory for sorting
export TMP_DIR_WRK=$TMP_DIR_LARGE/$PIPELINE_NAME"_align_"$DATA_NAME
mkdir -p $TMP_DIR_WRK

# help coordinate read merging via named pipe
export NAMED_PIPE=$TMP_DIR_WRK/align_$DATA_NAME
rm -f $NAMED_PIPE
mkfifo $NAMED_PIPE

# search all subdirectories of INPUT_DIR for cell FASTQ files
cd $INPUT_DIR
export CELL_DIRS=`ls -1 .`

# merge and sort the final output
cat $NAMED_PIPE |
samtools view -b - | 
samtools sort --threads $N_CPU -m $RAM_PER_CPU -T $TMP_DIR_WRK/tmp - | 
slurp -s 100M -o $BAM_FILE &

# align one cell at a time, with parallelization on each cell
perl $PIPELINE_DIR/align/align.pl

# index the final bam file
echo "indexing the sorted bam file"
samtools index -@ $N_CPU $BAM_FILE

# clean up
rm $NAMED_PIPE
rm -r $TMP_DIR_WRK

echo "done"
echo

