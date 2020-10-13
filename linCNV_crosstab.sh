#!/bin/bash

# application-specific data_script called by workflow launcher
# organizes a set of serial pipeline steps, with staged execution

# always call this support script first
source $PIPELINE_DIR/_workflow/workflow.sh

# make a crosstab table of fixed-width bin x cell, plus bin metrics
runWorkflowStep 1 crosstab/make_crosstab.sh

# apply normalization and make summary plots and reports
runWorkflowStep 2 crosstab/plot_crosstab.sh

