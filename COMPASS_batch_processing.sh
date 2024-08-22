#!/bin/bash

#SBATCH --job-name=hello_world
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=interactive
#SBATCH --account=default
#SBATCH --time=12:00:00

# This allows batch processing for a WT and mutant strain(s).
# The input (or raw) fastq files should have the format:
# $SAMPLE"_R1.fastq.gz" and $SAMPLE"_R2.fastq.gz"
# and be placed in the directory fastq/ which is underneath
# the main COMPASS_DIR such that RAW_FASTQ_DIR=$COMPASS_DIR"fastq/"
# This script is set up for sbatch job submissions so the user
# needs to be familiar with Slurm.

biological_samples="wt mutant"
replicates="1 2 3"
for biological_sample in $biological_samples;
do for replicate in $replicates;
do SAMPLE=$biological_sample"_"$replicate;
sbatch COMPASS.sh $SAMPLE
done
done
