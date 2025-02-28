#!/bin/bash -l

## vars
proj=b2013155
mail="nicolas.delhomme@umu.se"

## define a function
usage(){
    echo "This function take one argument as parameter; one of 'raw','trimmomatic','sortmerna'"
    echo "The UPSCb env. var. needs to be set to your Git UPSCb checkout directory."
    exit 1
}

## args number
if [ $# != 1 ]; then
    usage
fi

## process the argument
dir=
case "$1" in
    raw)
	dir="/proj/$proj/somaticEmbryogenesis/FastQC/$1"
	;;
    trimmomatic) 
	dir="/proj/$proj/nobackup/somaticEmbryogenesis/FastQC/$1"
	;;
    sortmerna)
	dir="/proj/$proj/nobackup/somaticEmbryogenesis/FastQC/$1"
	;;
esac

## stop if no dir
if [ -z $dir ]; then
    usage
fi

## check that the dir exists
if [ ! -d $dir ]; then
    usage
fi

## check that the UPSCb env var exists
if [ -z $UPSCb ]; then
    usage
fi

## submit
bash $UPSCb/pipeline/runFastQCMultiviewer.sh $dir
