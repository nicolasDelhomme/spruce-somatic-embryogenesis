#!/bin/bash

# vars
account=facility
mail=nicolas.delhomme@umu.se
in=/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/seidr/sf
out=/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/seidr/aggregate

# check
if [ -z $UPSCb ]; then
  echo "Your UPSCb env. var. needs to be set to your Git UPSCb checkout directory path" && exit 1
fi

# modules
module load bioinfo-tools seidr-devel

# dir
if [ ! -d $out ]; then
  mkdir -p $out
fi

# submit
sbatch -A $account --mail-user=$mail -e $out/aggregate.err -o $out/aggregate.err \
-J SSE-aggregate $UPSCb/pipeline/runSeidrAggregate.sh $out $in/*.sf 
