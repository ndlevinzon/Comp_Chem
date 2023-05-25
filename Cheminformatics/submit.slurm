#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --job-name=NDL_ChEMBL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=75G
#SBATCH --output=output.log
###################################################################

# ACTIVATE ANACONDA
eval "$(conda shell.bash hook)"
conda activate venv

python3 projB.py