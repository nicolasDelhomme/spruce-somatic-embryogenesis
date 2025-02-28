#!/bin/bash

## define a function
usage () {
    echo "This function take one argument as parameter; one of 'umea','uppmax'"
    echo "The UPSCb env. var. needs to be defined and point to your Git UPSCb checkout directory"
    exit 1
}

## args number
if [ $# != 1 ]; then
    usage    
fi

## env var
if [ -z $UPSCb ]; then
    echo "Set your UPSCb env. var."
    usage
fi

## vars
mail="nicolas.delhomme@umu.se"
proj=b2013155
export TRINITY_RNASEQ_ROOT=~/opt/trinityrnaseq_r20131110

## process the argument
runLoc=
case "$1" in
    umea)
	echo Not implemented yet.
	exit 1
	;;
    uppmax)	
	runLoc=$1
	in=/proj/$proj/nobackup/somaticEmbryogenesis/trimmomatic
	out=/proj/$proj/nobackup/somaticEmbryogenesis/diginorm
	;;
esac

## stop if no dir
if [ -z $runLoc ]; then
    usage
    exit 1
fi

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out/pair
    mkdir -p $out/left
    mkdir -p $out/right
fi

## clean
## umea only
# if [ -f $cfg/$cfgNum.in ]; then
#     rm -f $cfg/$cfgNum.in
# fi

## Combine the files
if [ ! -f $out/left.fq ]; then
    zcat $in/*_trimmomatic_1.fq.gz > $out/left.fq
fi

if [ ! -f $out/right.fq ]; then
    zcat $in/*_trimmomatic_2.fq.gz > $out/right.fq
fi

if [ ! -f $out/left-single.fq ]; then
    zcat $in/*_unpaired_1.fq.gz > $out/left-single.fq
fi

if [ ! -f $out/right-single.fq ]; then
    zcat $in/*_unpaired_2.fq.gz > $out/right-single.fq
fi

## process
## for f in `find $in -name $pattern -type f`; do echo `basename ${f//_[1,2].f*q.gz/}` ; done | sort | uniq | while read line;
## do
if [ $runLoc == "umea" ]; then
    zcat $in/*sortmerna_trimmomatic_1.fq.gz > $out/left.fq
    zcat $in/*sortmerna_trimmomatic_2.fq.gz > $out/right.fq
    echo bash `pwd`/../../../pipeline/runDigitalNormalization.sh -p 32 $out/left.fq $out/right.fq $out >> $cfg/$cfgNum.in
    rm $out/*.success
    echo bash `pwd`/../../../pipeline/runDigitalNormalization.sh -p 40 -s $out/left-single.fq $out >> $cfg/$cfgNum.in
    rm $out/*.success
    echo bash `pwd`/../../../pipeline/runDigitalNormalization.sh -p 40 -s $out/right-single.fq $out >> $cfg/$cfgNum.in
    cd $cfg
    runParallel.pl -t 1 -l $cfgNum -m $cfgNum > log-$cfgNum.out 2> log-$cfgNum.err
    ## more cleanup here!
    rm jellyfish* single* both* left-single.fq right-single.fq left.fq right.fq
else
    sbatch -A $proj --mail-user $mail -C "usage_mail" -e $out/diginorm-pair.err -o $out/diginorm-pair.out -J diginorm-pair-$proj $UPSCb/pipeline/runDigitalNormalization.sh -p 16 -m 120G $out/left.fq $out/right.fq $out/pair
    sbatch -A $proj --mail-user $mail -C "usage_mail" -e $out/diginorm-left.err -o $out/diginorm-left.out -J diginorm-left-$proj $UPSCb/pipeline/runDigitalNormalization.sh -p 16 -m 120G -s $out/left-single.fq $out/left
    sbatch -A $proj --mail-user $mail -C "usage_mail" -e $out/diginorm-right.err -o $out/diginorm-right.out -J diginorm-right-$proj $UPSCb/pipeline/runDigitalNormalization.sh -p 16 -m 120G -s $out/right-single.fq $out/right
fi

