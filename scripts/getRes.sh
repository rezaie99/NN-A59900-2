#!/bin/bash

function usage {

 echo 'This script takes an input directory and runs convert.py and process.py to run the time domain processing'
 echo ' and plot the raw data'

}

exeDir='.'
dir='st/'
saveDir='save/'
resType='ios'


while getopts "hd:e:f:r:t:s:" arg; do
  case $arg in
    h)
      usage
      exit
      ;;
    d)
      dir=$OPTARG
      #echo $dir
      ;;
    e)
       execDir=$OPTARG
       #echo $execDir
       ;;
    f)
       list=$OPTARG
       #echo $list
       ;;
    r)
       resDir=$OPTARG
       #echo $resDir
       ;;

    t)
       resType=$OPTARG
       #echo $resDir
       ;;
    o)
        options=$OPTARG
        #echo $options
        ;;

    s)
        saveDir=$OPTARG
        #echo $options
        ;;

  esac
done

list=$dir/*.txt

mkdir -p $saveDir

for file1 in $list
do
    echo $file1
    filepart=`basename -s .txt $file1`
    echo $filepart

    $exeDir/plotTruth.py  --truthFile $dir/$filepart".truth" --resultFile $resDir/$filepart$resType --save $saveDir/$filepart --noPlots

    $exeDir/sensorLog.py  $dir/$filepart --save $saveDir/$filepart --noPlots



done
