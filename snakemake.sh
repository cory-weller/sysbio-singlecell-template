#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 8G
#SBATCH --time 24:00:00

module purge
module load snakemake/9.23.1

snakemake ${@} --nolock --verbose -p --profile settings/profile.v9+.yaml  --configfile 'settings/config.yaml'  --conda-prefix envs --snakefile scripts/snakefile
