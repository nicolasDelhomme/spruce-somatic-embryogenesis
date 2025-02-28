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
export TRINITY_RNASEQ_ROOT=~/opt/trinityrnaseq_r20140413

## process the argument
runLoc=
case "$1" in
    umea)
	echo Not implemented yet.
	exit 1
	;;
    uppmax)	
	runLoc=$1
	in=/proj/$proj/somaticEmbryogenesis/STAR
	out=/proj/$proj/nobackup/somaticEmbryogenesis/trinity
	tmp=/glob/delhomme/tmp
	;;
    *) usage;;
esac

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## create the forward list of files
fwd=
for f in `find $in -type f -name "*mate1.gz"`; do
    if [ -z $fwd ]; then
	fwd=$f
    else
	fwd=$fwd","$f
    fi
done

## create the reverse list of files
rev=
for f in `find $in -type f -name "*mate2.gz"`; do
    if [ -z $rev ]; then
	rev=$f
    else
	rev=$rev","$f
    fi
done

if [ $runLoc == "umea" ]; then
    if [ -f $cfg/$cfgNum.in ]; then
	rm $cfg/$cfgNum.in
    fi
    echo bash $UPSCb/pipeline/runTrinity.sh -p 32 -m 50G $in/left.fq.normalized_K25_C2_pctSD200.fq.gz $in/right.fq.normalized_K25_C2_pctSD200.fq.gz $in/left-single.fq.normalized_K25_C2_pctSD200.fq.gz $in/right-single.fq.normalized_K25_C2_pctSD200.fq.gz $out $tmp >> $cfg/$cfgNum.in
    cd $cfg
    runParallel.pl -t 1 -l $cfgNum -m $cfgNum > log-$cfgNum.out 2> log-$cfgNum.err
else
    module load bioinfo-tools
    module load bowtie
    module load samtools
    sbatch -t 4-00:00:00 -A $proj --mail-user $mail -C "fat&usage_mail" -e $out/trinity.err -o $out/trinity.out -J trinity-$proj $UPSCb/pipeline/runTrinity.sh -m 240G -p 16 -n $out $tmp $fwd $rev
fi


