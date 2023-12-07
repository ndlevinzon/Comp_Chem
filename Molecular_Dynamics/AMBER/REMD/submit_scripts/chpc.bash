#!/bin/bash
#SBATCH --job-name=hremd-dihedral
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --tasks=8
#SBATCH --mem=0
#SBATCH --gres=gpu:2080ti:1
#SBATCH --account=cheatham-gpu-np
#SBATCH --partition=cheatham-gpu-np
set -x
module load gcc/8.5.0
module load mvapich2/2.3.7
module load intel-oneapi-mpi/2021.4.0
module load cuda/11.0
module load amber/20.20-gpu
mpirun -np 8 pmemd.cuda.MPI -rem 3 -ng 8 -groupfile groupfile
cd ../run.102/
sbatch prod.bash
