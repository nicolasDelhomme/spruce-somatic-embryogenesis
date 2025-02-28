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
export TRINITY_RNASEQ_ROOT=~/opt/trinity_beta_Jan28_2014_SGEplus

## process the argument
runLoc=
case "$1" in
    umea)
	echo Not implemented yet.
	exit 1
	;;
    uppmax)	
	runLoc=$1
	fasta=/proj/$proj/nobackup/somaticEmbryogenesis/trinity/others/Trinity.fasta
	fastq=/proj/$proj/nobackup/somaticEmbryogenesis/trimmomatic
	out=/proj/$proj/nobackup/somaticEmbryogenesis/RSEM/others
	tmp=/glob/delhomme/tmp
	;;
    *) usage;;
esac

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## umea
for f in `find $fastq -name "*_sortmerna_trimmomatic_[1,2].fq.gz" -type f`; do echo `basename ${f//_[1,2].f*q.gz/}` ; done | sort | uniq | while read line;
    do
    if [ $runLoc == "umea" ]; then
	if [ -f $cfg/$cfgNum.in ]; then
	    rm $cfg/$cfgNum.in
	fi
	echo bash $UPSCb/pipeline/runAlignReads.sh $out $tmp $fasta $fastq/${line}_1.fq.gz $fastq/${line}_2.fq.gz >> $cfg/$cfgNum.in
    else
	sbatch -A $proj --mail-user $mail -C "usage_mail" -e $out/$line.err -o $out/$line.out -J realign-$line $UPSCb/pipeline/runAlignReads.sh $out $tmp $fasta $fastq/${line}_1.fq.gz $fastq/${line}_2.fq.gz
    fi
done

if [ $runLoc == "umea" ]; then
	cd $cfg
	runParallel.pl -t 4 -l $cfgNum -m $cfgNum > log-$cfgNum.out 2> log-$cfgNum.err
fi


