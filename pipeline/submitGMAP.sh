#!/bin/bash -l

## error and verbose
set -ex

## module
module load myPipe

## define a function
usage () {
echo >&2 \
"This function take one argument as parameter:
     the location; one of 'umea','uppmax'
Note:You need to set the UPSCb env. variable to your UPSCb git checkout directory."
exit 1
}

## args number
if [ $# != 1 ]; then
    usage
    exit 1
fi

## check for the UPSCb env. var.                                                                                                  
if [ -z $UPSCb ]; then
    echo "You need to set the UPSCb environment variable"
    usage
fi

## process the argument
runLoc=
cfgNum=1
case "$1" in
    umea)
	runLoc=$1
	in=/mnt/picea/projects/spruce/22_Somatic_Embryogenesis_Project/fasta/somEmb.fa
	out=/mnt/picea/projects/spruce/22_Somatic_Embryogenesis_Project/gmap
	inxDir=/mnt/picea/storage/reference/Picea-abies/v1.0/indices/GMAP
	inx=Pabies1.0
	cfg=/mnt/picea/projects/spruce/22_Somatic_Embryogenesis_Project/cfg
	;;
    uppmax) 
	echo Not implemented yet.
	exit 1
	;;
    *)
	usage;;
esac

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## clean
if [ -f $cfg/$cfgNum.in ]; then
 rm $cfg/$cfgNum.in
fi

## process
if [ ! -d $out ]; then
    mkdir $out
fi

if [ $runLoc == "umea" ]; then
    echo bash $UPSCb/pipeline/runGmapl.sh -i 69000 $in $inxDir $inx $out >> $cfg/$cfgNum.in
else
    echo some sbatch
fi

if [ $runLoc == "umea" ]; then        
    cd $cfg
    runParallel.pl -t 1 -l $cfgNum -m $cfgNum > log-$cfgNum.out 2> log-$cfgNum.err
fi
