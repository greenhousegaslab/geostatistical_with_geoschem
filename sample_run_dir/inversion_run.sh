#!/bin/bash
#SBATCH -J GIM_LNLG_2021
#SBATCH --partition defq
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=20000
#SBATCH --time=72:00:00
#SBATCH --output=GIM.out
#SBATCH --dependency=singleton

module load matlab

matlab -nosplash -nodesktop -nodisplay < /scratch4/smill191/smiller/geoschem_adjoint_co2/LNLG_2021/inversion_run.m
