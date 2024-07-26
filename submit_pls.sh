#!/bin/bash

#SBATCH -D /home/pd/sstuenzi/CG_CryoThermokarst
#SBATCH -o slurm_output_%A_%a.out
#SBATCH -e slurm_error_%A_%a.err
#SBATCH -p PermaRisk
#SBATCH -J Simone
#SBATCH --exclude=pls11
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=simone.stuenzi@awi.de

matlab -nodisplay -nosplash -nodesktop -r "run('run_CG.m');exit;"