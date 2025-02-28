#!/bin/bash

# define the thresholds
thresholds=( 2.33 2.05 1.88 1.75 1.64 1.55 1.48 1.41 1.34 1.28 )

# Load the tools
module load bioinfo-tools seidr-devel

# process the argument
network=/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/seidr/aggregate/aggregated.sf
out=/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/seidr/backbone

if [ ! -d $out ]; then
  mkdir -p $out
fi

# submit
for i in {0..9}; do
  j=$(expr $i + 1)
  
  sbatch -o $out/backbone-${j}-percent.out \
  -e $out/backbone-${j}-percent.err \
  $UPSCb/pipeline/runSeidrBackbone.sh $network ${thresholds[$i]} \
  $out/backbone-${j}-percent.sf
done
