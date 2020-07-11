#!/bin/bash

echo "scanning genome for CNVs independently of cell clustering"

# NB: BIN_PREFIX still carries the value $REBIN_PREFIX

# process the bin x cell matrix in R
Rscript $PIPELINE_DIR/scan/scan.R
checkPipe

echo "done"
echo

