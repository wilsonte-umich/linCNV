#!/bin/bash

export REMOTE_PORT=$1
export SERVER_DIR=$2
export DATA_PATH=$3
export R_DIRECTORY=$4
export APP=$5
export GENOMES_DIR=""

# check for server currently running on REMOTE_PORT
PID_FILE=linCnv-remote-pid-$REMOTE_PORT.txt
MDI_PID=""
if [ -e $PID_FILE ]; then MDI_PID=`cat $PID_FILE`; fi
EXISTS=""
if [ "$MDI_PID" != "" ]; then EXISTS=`ps -p $MDI_PID | grep -v PID`; fi 

# launch Shiny if not already running
SEPARATOR="---------------------------------------------------------------------"
WAIT_SECONDS=15
if [ "$EXISTS" = "" ]; then
    echo $SEPARATOR 
    echo "Please wait $WAIT_SECONDS seconds for the web server to start"
    echo $SEPARATOR
    $R_DIRECTORY/Rscript $SERVER_DIR/$APP/remote.R &
    MDI_PID=$!
    echo "$MDI_PID" > $PID_FILE
    sleep $WAIT_SECONDS # give Shiny time to start up before showing further prompts   
fi
trap "kill -9 $MDI_PID; exit" SIGINT SIGQUIT SIGHUP

# report the PID to the user
echo $SEPARATOR
echo "Web server process running on remote port $REMOTE_PORT as PID $MDI_PID"

# report on browser usage within the command shell on user's local computer
echo $SEPARATOR
echo "Please point any web browser to:"
echo
echo "http://127.0.0.1:$REMOTE_PORT"
echo

# prompt for exit action, with or without killing of the R web server process
USER_ACTION=""
while [ "$USER_ACTION" != "1" ]; do
    echo $SEPARATOR
    echo "To close the remote server connection:"
    echo
    echo "  1 - close the connection AND stop the web server"
    echo
    echo "Select an action (type '1' and hit Enter):"
    read USER_ACTION
done

# kill the web server process if requested
if [ "$USER_ACTION" = "1" ]; then
    kill -9 $MDI_PID
fi 

# send a final helpful message
# note: the ssh process on client will NOT exit when this script exits since it is port forwarding still
echo
echo "You may now safely close this command window."
