#!/bin/bash

# application-specific master shell script called by workflow launcher
# organizes a set of serial pipeline steps, with staged execution

# always call this support script first
source $PIPELINE_DIR/_workflow/workflow.sh

# parse CellRanger bam data into variable-width bins in preparation for cell marking
runWorkflowStep 1 bin/bin.sh            

