#!/bin/bash

function usage {

 echo 'This script takes the log file name, timestamp, dst directory and user id as inputs'
 echo 'cpfiles.sh <log file name> <YYMMDD_HHMMDD> <dst directory> <user id>'


}

dir='st'
exeDir='.'
list=$dir/*.txt


while getopts "h" arg; do
  case $arg in
    h)
      usage
      exit
      ;;
  esac
done

logFile=$1
timeStamp=$2
dst=$3
user=$4

echo $timeStamp\_$user.txt

cp $logFile.txt $dst/$timeStamp\_$user.txt
cp $logFile.cfg $dst/$timeStamp\_$user.cfg
cp $logFile.BIN $dst/$timeStamp\_$user.BIN