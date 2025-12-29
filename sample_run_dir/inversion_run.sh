#!/bin/bash
#SBATCH --job-name=IS_2016
#SBATCH --partition=parallel
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=12000
#SBATCH --time=72:00:00
#SBATCH --mail-user=kmorga46@jhu.edu
#SBATCH --dependency=singleton

module load matlab/R2025b

matlab -nosplash -nodesktop -nodisplay < /scratch4/smill191/kmorgan/OCO2_MIP_v11/IS/IS_2016/inversion_run.m > GIM.out
