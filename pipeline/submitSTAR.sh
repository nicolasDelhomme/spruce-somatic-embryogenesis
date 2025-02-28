#!/bin/bash -l

## global vars
proj=b2013155
mail="nicolas.delhomme@umu.se"
genome=/proj/b2011227/indices/STAR/Pabies01-genome

## create the out dir
in=/proj/$proj/nobackup/somaticEmbryogenesis/trimmomatic
out=/proj/$proj/somaticEmbryogenesis/STAR

if [ ! -d $out ]; then
    mkdir -p $out
fi

## run
for f in `find $in -name "*_sortmerna_trimmomatic_[1,2].fq.gz" -type f`; do echo `basename ${f//_[1,2].f*q.gz/}` ; done | sort | uniq | while read line;
do
    sbatch -A $proj -C "fat&usage_mail" --mail-user $mail -e $out/$line.err -o $out/$line.out -J STAR-$line /proj/b2011227/nobackup/spruceRoots/STAR-Alignments/runSTAR.sh $in/${line}_1.fq.gz $in/${line}_2.fq.gz $genome $out/$line
done
