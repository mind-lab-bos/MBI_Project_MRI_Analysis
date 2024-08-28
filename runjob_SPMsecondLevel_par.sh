#!/bin/bash

#SBATCH --job-name=SPMsecondLevel
#SBATCH --partition=short
#SBATCH --output='logs/slurm_SPMsecondLevel-%j.out'
#SBATCH --error='logs/slurm_SPMsecondLevel-%j.err'
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G        # memory per cpu-core (4G is default) #1
#SBATCH --time=05:00:00          # total run time limit (HH:MM:SS) 
#SBATCH --mail-type=all          # send email on job start, end and fault. 
#SBATCH --mail-user=xi.wu@northeastern.edu
#SBATCH --array=6

# Choose a version of MATLAB by loading a module:
module load matlab/R2023a

# Remove -singleCompThread below if you are using parallel commands:
matlab -batch "runPipeline_SPMsecondLevel($SLURM_ARRAY_TASK_ID)"

