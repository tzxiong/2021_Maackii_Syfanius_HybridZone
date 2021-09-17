#!/bin/bash

seed=$1
gen=$2
trial=$3
sampleSize=$4
dir=$5

module load parallel

cd ${dir}/slim_seed${seed}


## split target VCF into per-chromosome files

echo compressing..

bgzip -c slim_seed${seed}_gen${gen}.sampleSize.${sampleSize}.trial.${trial}.vcf > slim_seed${seed}_gen${gen}.sampleSize.${sampleSize}.trial.${trial}.vcf.gz

bcftools index slim_seed${seed}_gen${gen}.sampleSize.${sampleSize}.trial.${trial}.vcf.gz

echo separating VCFs..

parallel -j1 "bcftools view -r {1} slim_seed${seed}_gen${gen}.sampleSize.${sampleSize}.trial.${trial}.vcf.gz -o slim_seed${seed}_gen${gen}.sampleSize.${sampleSize}.trial.${trial}.chr.{2}.vcf -Ov" :::: ../../Intervals.list :::+ {1..10}

sleep 1s

echo converting to bimbam..

## split each per-chromosome VCF into bimbam files for each population
parallel -j1 "/n/home00/txiong/miniconda3/envs/EnvPy3.9/bin/python3 /n/home00/txiong/Research/2020_BarriersToGeneFlow/SLiM_simulation/Entropy/ELAI/vcf2bimbam.py --seed $seed --gen $gen --sampleSize $sampleSize --trial $trial --chrom {1} --dir $dir" ::: {1..10}

## run ELAI on the bimbam output

echo running ELAI

parallel -j1 "sbatch /n/home00/txiong/Research/2020_BarriersToGeneFlow/SLiM_simulation/Entropy/ELAI/bimbam2elai.sh $seed $gen $sampleSize $trial {} $dir" ::: {1..10}