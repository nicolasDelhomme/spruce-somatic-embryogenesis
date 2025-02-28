#!/bin/bash

# Load the tools
module load bioinfo-tools seidr-devel

# arguments
mail=nicolas.delhomme@umu.se
network=/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/seidr/aggregate/aggregated.sf
out=/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/seidr/threshold

if [ ! -d $out ]; then
  mkdir -p $out
fi

sbatch -o $out/threshold.out -e $out/threshold.err -J SSE-threshold --mail-user=$mail \
  $UPSCb/pipeline/runSeidrThreshold.sh $network $out/threshold.txt
