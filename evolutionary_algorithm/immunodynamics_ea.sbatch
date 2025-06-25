#!/bin/bash -l
#SBATCH --job-name=immunodynamics_evolutionary_algorithm
#SBATCH --account=def-nathoo # adjust this to match the accounting group you are using to submit jobs
#SBATCH --time=0-08:00         # adjust this to match the walltime of your job
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12      # adjust this if you are using parallel commands
#SBATCH --mem=212G             
#SBATCH --mail-user=pvelidi@uvic.ca
#SBATCH --mail-type=ALL

# Choose a version of MATLAB by loading a module:
module load matlab/2023b.2
matlab -batch "EvolutionaryAlgorithm"
