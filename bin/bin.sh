#!/bin/bash

echo "setting bin endpoints and counting reads per bin per cell"

# create a single exclusion file from the gap and bad-region files
cat \
<(cut -f 1-3 $GAP_FILE       | awk 'BEGIN{OFS="\t"}{print $0, "gap"}') \
<(cut -f 1-3 $BAD_REGIONS_FILE | awk 'BEGIN{OFS="\t"}{print $0, "bad_region"}') |
sort -k1,1 -k2,2n |
bedtools merge -i stdin -c 4 -o distinct > $EXCLUSIONS_FILE
checkPipe

# execute the binning process
source $PIPELINE_DIR/bin/run_binning.sh

# plot cells for marking
export CELL_PLOT_DIR=$PLOT_DIR/cells
mkdir -p $CELL_PLOT_DIR
Rscript $PIPELINE_DIR/bin/plot_cells.R
checkPipe

# clean up
rm $BIN_PREFIX.counts.*.gz
rm $BIN_PREFIX.counts.gz
rm $BIN_PREFIX.gc_fraction.gz

echo "done"
echo

