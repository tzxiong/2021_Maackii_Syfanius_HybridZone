#!/bin/bash
#SBATCH -J elai
#SBATCH -n 5             # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 7-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p unrestricted,shared           # Partition to submit to
#SBATCH --mem 25gb           # Memory pool for all cores (see also --me$
#SBATCH -o elai_%j.out  # File to which STDOUT will be written, %j i$
#SBATCH -e elai_%j.err  # File to which STDERR will be written, %j i$

module load parallel

seed=$1
gen=$2
sampleSize=$3
trial=$4
chr=$5
dir=$6

HOME_DIR=/home_dir

cd ${dir}/slim_seed${seed}

entries=($(shuf -i 0-100000000 -n 150)) # 150 random numbers, need to change if different number of combinations of generations and repeats are launched

parallel "echo {1} {2}" ::: 5000 500000 ::: {1..50} | parallel --colsep ' ' -j5 "$HOME_DIR/Software/elai/elai-lin -g slim_seed${seed}_gen${gen}.sampleSize.${sampleSize}.trial.${trial}.chr.${chr}.bimbam.pop1.txt -p 10 -g slim_seed${seed}_gen${gen}.sampleSize.${sampleSize}.trial.${trial}.chr.${chr}.bimbam.pop2.txt -p 11 -g slim_seed${seed}_gen${gen}.sampleSize.${sampleSize}.trial.${trial}.chr.${chr}.bimbam.admx.txt -p 1 -pos slim_seed${seed}_gen${gen}.sampleSize.${sampleSize}.trial.${trial}.chr.${chr}.bimbam.position.txt -s 30 -o slim_seed${seed}_gen${gen}.sampleSize.${sampleSize}.trial.${trial}.chr.${chr}.bimbam.admx.output.elaiGen.{1}.run.{2} -mg {1} -C 2 -c 10 -R {3} --exclude-nopos" :::: - :::+ "${entries[@]}"
