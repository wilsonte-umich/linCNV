#!/bin/bash

echo "setting bin endpoints and counting reads per bin per cell"

# collect information from 'bin' analysis and user's manual cell calling
# also apply a cell depth filter here to purge v. low coverage cells
export ACCEPTED_CELLS=`Rscript $PIPELINE_DIR/rebin/get_cells.R`
checkPipe

# force log and file naming to rebin (allows code re-use)
export FORCE_SCRIPT=rebin 
export BIN_PREFIX=$REBIN_PREFIX

# execute the re-binning process
source $PIPELINE_DIR/bin/run_binning.sh

# clean up
rm $BIN_PREFIX.counts.*.gz
rm $BIN_PREFIX.counts.gz
rm $BIN_PREFIX.gc_fraction.gz

echo "done"
echo

