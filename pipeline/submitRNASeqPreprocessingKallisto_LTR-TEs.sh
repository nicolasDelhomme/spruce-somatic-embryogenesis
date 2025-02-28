#!/bin/bash -l

set -ex # 'set -e' arrête le script quand il y a une erreur à une étape pour éviter de continuer les étapes suivantes pour rien avec des erreurs
		# 'set -x' print et interprete les variables pour se rendre compte là où on a mal tapé un nom de variable

# environment variables
proj=u2015037

mail="camillacanovi@yahoo.it"

in=/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/trimmomatic

out=/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/kallisto_LTR-TEs

kallistoFasta=/mnt/picea/storage/reference/Picea-abies/v1.1/fasta/Pabies1.0-all-phase.gff3.CDSandLTR-TE.fa
			  
kallistoIndex=/mnt/picea/storage/reference/Picea-abies/v1.1/indices/kallisto/Pabies1.0-all-phase.gff3.CDSandLTR-TE.inx

if [ ! -d $out ]; then
    mkdir -p $out
  fi

#module load bioinfo-tools kallisto

if [ -z $UPSCb ]; then # -z teste si on a la variable UPSCb dans notre environnement
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir." 
    exit 1
fi

# loop to get every PE fastq file and apply the script runRNASeqPreprocessing.sh to it
for f in `find $in -name "*_[1,2].fq.gz"`; do echo "${f//_[1,2].fq.gz/}" ; done | sort | uniq | while read line; 
# uniq - report or omit repeated lines (With no options, matching lines are merged to the first occurrence.)
do 
fnam=$(basename $line)
sbatch -A $proj --mail-user=$mail -e $out/$fnam.err -o $out/$fnam.out ~/Git/UPSCb/pipeline/runKallisto.sh -r \
${line}_1.fq.gz ${line}_2.fq.gz $kallistoIndex $kallistoFasta $out
# options of runRNASeqPreprocessing.sh, run it for every PE fastq file called f in the loop
done
