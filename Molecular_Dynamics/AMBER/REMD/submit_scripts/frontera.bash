#!/bin/bash
#SBATCH -J Mremd-OL3
#SBATCH -p rtx
#SBATCH -N 16
#SBATCH -n 192
#SBATCH -t 24:00:00
#SBATCH -A MCB20008
#SBATCH --mail-type=all
#SBATCH --mail-user=lauren.winkler@utah.edu
set -x
module load mvapich2-x/2.3
module load cuda/11.0

DIR="OUT"
ibrun -np 192 /home1/07450/lgw_19/amber22/amber22/bin/pmemd.cuda.MPI -ng 192  -groupfile groupfile -remd-file remd.dim -remlog remlog
