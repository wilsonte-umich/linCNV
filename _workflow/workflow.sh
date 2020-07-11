
# workflow.sh has utility functions for managing a bash-based workflow script

#--------------------------------------------------------------------
# get and set the last successfully completed step in a multi-step serial workflow
#--------------------------------------------------------------------
function setStatusFile { # construct a rule-based status file name
    if [ "$OUTPUT_DIR" = "" ]; then
        echo "missing variable: OUTPUT_DIR"
        exit 1
    fi
    if [ ! -d "$OUTPUT_DIR" ]; then
        echo "directory does not exist: $OUTPUT_DIR"
        exit 1
    fi
    if [ "$PIPELINE_NAME" = "" ]; then
        echo "missing variable: PIPELINE_NAME"
        exit 1
    fi
    if [ "$PIPELINE_COMMAND" = "" ]; then
        echo "missing variable: PIPELINE_COMMAND"
        exit 1
    fi 
    if [ "$DATA_NAME" = "" ]; then
        echo "missing variable: DATA_NAME"
        exit 1
    fi
    LOGS_DIR=$OUTPUT_DIR/$PIPELINE_NAME"_logs"
    if [ ! -d "$LOGS_DIR" ]; then mkdir $LOGS_DIR; fi    
    STATUS_FILE=$LOGS_DIR/$DATA_NAME.$PIPELINE_NAME.status
    if [ ! -e "$STATUS_FILE" ]; then touch $STATUS_FILE; fi        
}
function getWorkflowStatus {
    setStatusFile
    LAST_SUCCESSFUL_STEP=`awk '$1=="'$PIPELINE_COMMAND'"' $STATUS_FILE | tail -n1 | cut -f2`
    if [ "$LAST_SUCCESSFUL_STEP" = "" ]; then LAST_SUCCESSFUL_STEP=0; fi
}
function setWorkflowStatus {
    setStatusFile
    STEP_NUMBER=$1
    STEP_SCRIPT=$2
    DATE=`date`
    echo -e "$PIPELINE_COMMAND\t$STEP_NUMBER\t$STEP_SCRIPT\t$DATE" >> $STATUS_FILE
}
function resetWorkflowStatus { # override any prior job outcomes and force a new status
    setStatusFile
    if [ "$LAST_SUCCESSFUL_STEP" != "" ]; then
        awk '$1!="'$PIPELINE_COMMAND'"||$2<='$LAST_SUCCESSFUL_STEP $STATUS_FILE >> $STATUS_FILE.tmp
        mv -f $STATUS_FILE.tmp $STATUS_FILE
    fi
}
function showWorkflowStatus {
    setStatusFile
    STATUS_LINE_LENGTH=`cat $STATUS_FILE | wc -l`
    if [ "$STATUS_LINE_LENGTH" -gt "0" ]; then
        echo -e "COMMAND\tSTEP\tSCRIPT\tDATE"
        cat $STATUS_FILE        
    fi
}
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# execute a pipeline step by sourcing its script, if not already successfully completed
#--------------------------------------------------------------------
function runWorkflowStep {   
    STEP_NUMBER=$1
    STEP_SCRIPT=$2
    getWorkflowStatus    
    if [ "$LAST_SUCCESSFUL_STEP" -lt "$STEP_NUMBER" ]; then    
        source $PIPELINE_DIR/$STEP_SCRIPT # NB: script is responsible for calling checkPipe        
        setWorkflowStatus $STEP_NUMBER $STEP_SCRIPT
    else
        echo "already succeeded: $PIPELINE_COMMAND"," step $STEP_NUMBER, $STEP_SCRIPT"
    fi    
}
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# ensure that all commands in a pipe had exit_status=0
#--------------------------------------------------------------------
function checkPipe {  
    local PSS=${PIPESTATUS[*]}
    for PS in $PSS; do
       if [[ ( $PS > 0 ) ]]; then
           echo "pipe error: [$PSS]"
           exit 99
       fi;
    done
}
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# ensure a job will have data to work on
#--------------------------------------------------------------------
function checkForData { # ensure that a data stream will have at least one line of data
    local COMMAND="$1"
    if [ "$COMMAND" = "" ]; then
        echo "checkForData error: system command not provided"
        exit 100
    fi
    local LINE_1="`$COMMAND | head -n1`"
    if [ "$LINE_1" = "" ]; then
        echo "no data; exiting quietly"
        exit 0
    fi
}
function waitForFile {  # wait for a file to appear on the file system; default timeout=60 seconds
    local FILE="$1"
    local TIME_OUT="$2"
    if [ "$FILE" = "" ]; then
        echo "waitForFile error: file not provided"
        exit 100
    fi
    if [ "$TIME_OUT" = "" ]; then
        local TIME_OUT=60
    fi
    local ELAPSED=0
    while [ ! -s $FILE ]
    do
        sleep 2;
        let "ELAPSED += 2"
        if [ "$ELAPSED" -gt "$TIME_OUT" ]; then
            echo "waitForFile error: $FILE not found after $TIME_OUT seconds"
            exit 100
        fi
    done
}
function checkFileExists {  # verify non-empty file, or first of glob if called as checkFileExists $GLOB
    local FILE="$1"
    if [ "$FILE" = "" ]; then
        echo "checkFileExists error: file not provided"
        exit 100
    fi
    if [ ! -s "$FILE" ]; then
        echo "file empty or not found on node "`hostname`
        echo $FILE
        exit 100
    fi
}
#--------------------------------------------------------------------

