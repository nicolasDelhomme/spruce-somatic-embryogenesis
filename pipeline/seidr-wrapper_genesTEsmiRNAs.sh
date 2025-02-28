#!/bin/bash -l

# safeguards
set -ex

# project vars
account=u2015040
mail=katja.stojkovic@umu.se

# check
if [ -z $UPSCb ]; then
  echo "The UPSCb environment variable needs to be set to your UPSCb Git checkout path"
  exit 1
fi

# source
source $UPSCb/src/bash/functions.sh

# modules
module load bioinfo-tools seidr-devel
#source /pfs/nobackup/home/b/bastian/seidr/build/sourcefile

# Variables
gset=/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/seidr_geneTEmiRNA/miRNAs.tsv
inference=(aracne clr genie3 llr-ensemble narromi pcor pearson plsnet spearman tigress el-ensemble)

# 14 workers on 2 nodes (kk has 28 per node)
# --tasks-per-node is max number of node / c, so for c 7, it's 4
# narromi is not thread safe, hence -c 1
default="-n 2 -c 14 -t 1-00:00:00 --tasks-per-node 2"
arguments=([0]=$default [1]=$default [2]="-n 2 -c 28 -t 2-00:00:00" [3]=$default [4]="-n 28 -t 1-00:00:00" [5]="-n 7 -t 12:00:00" [6]="-n 7 -t 12:00:00" [7]=$default [8]="-n 7 -t 12:00:00" [9]="-n 1 -c 28 -t 2-00:00:00" [10]="-n 1 -c 28 -t 2-00:00:00")
parallel="-O "'$SLURM_CPUS_PER_TASK'
command=([0]="mi -m aracne "$parallel [1]="mi -m CLR "$parallel [2]="genie3 "$parallel [3]="llr-ensemble "$parallel [4]="narromi "$parallel [5]="pcor" [6]="correlation -m pearson" [7]="plsnet "$parallel [8]="correlation -m spearman" [9]="tigress "$parallel [10]="el-ensemble"$parallel)

# usage
USAGETXT=\
"
$0 <expression-matrix.tsv> <genes.tsv>

Methods to be run: ${inference[@]}

"

# number of arguments
if [ $# -ne 2 ]; then
  abort "This script expects 2 arguments"
fi

# input: gene and data
if [ ! -f $1 ]; then
  abort "The first argument needs to be the expression matrix tab delimited file"
fi

if [ ! -f $2 ]; then
  abort "The second argument needs to be the gene names tab delimited file"
fi

# Create template script and submit
len=${#inference[@]}
for ((i=0;i<len;i++)); do
  inf=${inference[$i]}
  if [ ! -f results/$inf/$inf.tsv ]; then
    mkdir -p results/$inf
    echo "#!/bin/bash" > results/$inf/$inf.sh
    echo "unset OMP_NUM_THREADS" >> results/$inf/$inf.sh
    echo "srun ${command[$i]} -i $1 -g $2 -t $gset -o results/$inf/$inf.tsv" >> results/$inf/$inf.sh
    sbatch --mail-type=ALL --mail-user=$mail -A $account -J $inf -e results/$inf/$inf.err -o results/$inf/$inf.out ${arguments[$i]} results/$inf/$inf.sh
  fi
done

