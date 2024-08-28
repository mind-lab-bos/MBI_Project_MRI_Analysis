change directory to the /code folder

go to compute node first on terminal

`srun --partition=short --nodes=1 --cpus-per-task=1 --pty /bin/bash`

edit the runjob_SPMfirstLevel_par.sh file - change username to your own, change number of subjects; edit params file

then `sbatch runjob_SPMfirstLevel_par.sh` in terminal

check .err and .out files in the logs folder. if the size of .err is 0 it means no error