#!/bin/bash

echo "filtering cells and bins and normalizing data baselines"

# NB: BIN_PREFIX still carries the value $REBIN_PREFIX

# process the bin x cell matrix in R
Rscript $PIPELINE_DIR/normalize/normalize.R
checkPipe

echo "done"
echo

