#!/bin/bash

## args
mail="nicolas.delhomme@umu.se"
proj=b2013155
in=/proj/$proj/nobackup/somaticEmbryogenesis/sortmerna
out=/proj/$proj/nobackup/somaticEmbryogenesis/trimmomatic

if [ ! -d $out ]; then
    mkdir -p $out
fi

if [ -z $UPSCb ]; then
    echo "The UPSCb environment variable needs to be set."
    exit 1
fi

if [ ! -f $UPSCb/data/adapter.fa ]; then
    echo "Either your UPSC env. var. is not set correctly or your checkout is too old."
    exit 1
fi

for f in `find $in -type f -name "*_sortmerna_[1,2].fq.gz"`; do echo `basename ${f//_[1,2].fq.gz*/}` ; done | sort | uniq | while read line;
do sbatch -A $proj --mail-user $mail -C "usage_mail" -e $out/$line.err -o $out/$line.out -J Trim-$line $UPSCb/pipeline/runTrimmomatic.sh -c "ILLUMINACLIP:$UPSCb/data/adapter.fa:1:30:9" $in/${line}_1.fq.gz $in/${line}_2.fq.gz $out 33 SLIDINGWINDOW:7:20 MINLEN:50
done

