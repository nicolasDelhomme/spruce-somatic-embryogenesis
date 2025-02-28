#!/bin/bash -l

proj=b2013155
tmp=/glob/delhomme/tmp
mail="nicolas.delhomme@umu.se"
out=/proj/$proj/nobackup/somaticEmbryogenesis/sortmerna
in=/proj/$proj/somaticEmbryogenesis/raw

if [ ! -d $out ]; then
    mkdir -p $out
fi

export SORTMERNADIR=/home/delhomme/sortmerna

if [ -z $UPSCb ]; then
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir."
fi

for f in `find $in -name "*_[1,2].fastq.gz" -type l`; do echo "${f//_[1,2].fastq.gz/}" ; done | sort | uniq | while read line;
do 
fnam=`basename $line`
sbatch -A $proj --mail-user $mail -C "usage_mail" -e $out/$fnam.err -o $out/$fnam.out -J smr-$fnam $UPSCb/pipeline/runSortmerna.sh -m $out $tmp ${line}_1.fastq.gz ${line}_2.fastq.gz
exit
done

