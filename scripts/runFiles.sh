#!/bin/bash

function usage {

 echo 'This script takes an input directory and runs through the desired script (python,c,plotTruth)

        Required inputs:
                    -d input directory
                    -s save directory

        Requires at least one of:
                    -y runs it through python : trackStreaming.py
                    -c runs it through C : sigprocDataTest
                    -p runs plotTruth.py with C results and saves it to the (-s) save directory

        Optional Arguments:
                    -q runs the script in parallel
                    -l runs specific files from the input directory
                    -o additional options to be appended to trackStreaming.py
                    -e execution directory
                    
        '
}

exeDir='.'
dir='st/'
saveDir=save
channelNum=2

runAll=false
runIos=false
runPlots=false
quickRun=false
runPython=false

while getopts "hd:e:s:o:l:a:yqpc" arg; do
  case $arg in

    h)
      usage
      exit
      ;;

    d)
      dir=$OPTARG
      echo sourceDir = $dir
      ;;

    e)
       execDir=$OPTARG
       echo $execDir
       ;;

    s)
       saveDir=$OPTARG
       echo saveDir = $saveDir
       ;;

    o)
        options=$OPTARG
        echo $options
        ;;

    c)
        runIos=true
        echo runIos = $runIos
        ;;

    p)
        runPlots=true
        echo runPlots = $runPlots
        ;;

    l)
        list=$OPTARG
        echo "Selected file(s) = $list"
        ;;

    a)
        runAll=$OPTARG
        echo runAll =  $runAll
        ;;

    q)
        quickRun=true
        echo "Quick Run In the Background"
        ;;

    y)
        runPython=true
        echo runPython = $runPython
        ;;

  esac

done

if [ x$list = x ] ; then
    list=$dir/*.opt
fi


mkdir -p $saveDir


for file1 in $list
do
    echo $file1
    filepart=`basename -s .opt $file1`
    echo $filepart
    truth="$dir/$filepart.truth"


    if [ "$runPython" = true ] ; then

            if [ "$quickRun" = true ] ; then

                $exeDir/trackStreaming.py $dir/$filepart --save $saveDir/$filepart --noPlots $options &

            else

                $exeDir/trackStreaming.py $dir/$filepart --save $saveDir/$filepart $options

            fi

    fi


    if [ "$runIos" = true ] ; then


            if [ "$quickRun" = true ] ; then

                ../targets/iOs/Build/Products/Debug/sigprocDataTest -i $dir/$filepart -o $saveDir/$filepart &

            else

                ../targets/iOs/Build/Products/Debug/sigprocDataTest -i $dir/$filepart -o $saveDir/$filepart -d


            fi

    fi



    if [ "$runPlots" = true ] ; then

        ./plotTruth.py --truthFile $truth --resultFile $saveDir/$filepart.c.hr.txt --save $saveDir/$filepart

    fi



done
