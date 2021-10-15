#!/bin/bash
#SBATCH --time=23:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p short

matlab -nodisplay -batch "MBI_GAMMA_CONN_Submit"