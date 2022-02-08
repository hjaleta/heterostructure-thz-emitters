#!/bin/bash -l
#SBATCH -J W3-2-L30
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 160:00:00
#SBATCH -A snic2021-1-22

srun /proj/oppeneer/users/x_hjali/src-V2.6/tetralith/main.e < input.dat> oput.dat