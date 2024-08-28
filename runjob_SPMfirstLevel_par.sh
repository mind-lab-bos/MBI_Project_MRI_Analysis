#!/bin/bash

#SBATCH --job-name=SPMfirstLevel
#SBATCH --partition=short
#SBATCH --output='logs/slurm_SPMfirstLevel-%j.out'
#SBATCH --error='logs/slurm_SPMfirstLevel-%j.err'
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G        # memory per cpu-core (4G is default) #1
#SBATCH --time=00:45:00          # total run time limit (HH:MM:SS) 12 indTT 26cTT
#SBATCH --mail-type=all          # send email on job start, end and fault. 
#SBATCH --mail-user=xi.wu@northeastern.edu
#SBATCH --array=1-24            # change according to num of subs (=isub)

# Choose a version of MATLAB by loading a module:
module load matlab/R2023a

# Remove -singleCompThread below if you are using parallel commands:
matlab -batch "runPipeline_SPMfirstLevel($SLURM_ARRAY_TASK_ID)"

