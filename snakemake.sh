#!/bin/sh

#SBATCH --mem=800g
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=36
#SBATCH --partition=largemem
#SBATCH --mail-type=ALL,TIME_LIMIT_50

module load snakemake

snakemake --use-conda --cores 36
