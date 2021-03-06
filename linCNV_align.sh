#!/bin/bash

# application-specific data_script called by workflow launcher
# organizes a set of serial pipeline steps, with staged execution

# always call this support script first
source $PIPELINE_DIR/_workflow/workflow.sh

# align reads to genome and create output files that mimic CellRanger DNA
runWorkflowStep 1 align/align.sh

# make a table with cells counts to help summarize relative cell quality
runWorkflowStep 2 align/count.sh

#!!!!!!!!!!!!!!!!!!!
#exit 666
#!!!!!!!!!!!!!!!!!!!

# when done, provide the directory
#   $OUTPUT_DIR/outs
# to commands 'bin' and 'analyze' as variable
#   $CELL_RANGER_DIR

