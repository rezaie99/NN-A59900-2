#!/bin/bash

function usage {

 echo 'This script takes an input directory and runs convert.py and process.py to run the time domain processing'
 echo ' and plot the raw data'

}

dir='st'
exeDir='.'
csv=''



while getopts "hd:e:f:c" arg; do
  case $arg in
    h)
      usage
      exit
      ;;
    d)
      dir=$OPTARG
      echo $dir
      ;;
    e)
       execDir=$OPTARG
       echo $execDir
       ;;
    f)
       list=$OPTARG
       echo $list
       ;;
    c)
        csv=".csv"
  esac
done

list=$dir/*.txt

for file1 in $list
do
    echo $file1
    filepart=`basename -s .txt $file1`
    echo $filepart
    truth="$dir/$filepart.interpolated.truth.csv"
    echo $truth


    $exeDir/process.py $dir/$filepart.txt --truthFile $truth &
    processPid=$!
    $exeDir/sensorLog.py $dir/$filepart
    kill $processPid



done
