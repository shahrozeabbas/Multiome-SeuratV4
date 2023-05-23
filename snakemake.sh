#!/bin/sh

#SBATCH --mem=400g
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=20
#SBATCH --partition=largemem
#SBATCH --mail-type=ALL,TIME_LIMIT_50

source /data/abbass2/Apps/conda/bin/activate snakes

snakemake --use-conda --cores 16