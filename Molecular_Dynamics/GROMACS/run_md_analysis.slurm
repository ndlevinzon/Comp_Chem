#!/bin/bash
#SBATCH --job-name=gmx_analysis
#SBATCH --output=slurm/md_analysis_%j.out
#SBATCH --error=slurm/md_analysis_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=36   # One per replica
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --exclusive 
#SBATCH --mem=0

# Load GROMACS module
module purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 gromacs/2024.4


# ----------- STEP 1: LOOP THROUGH REPLICAS ---------------

echo "Processing MD run..."

OUTDIR="md_analysis_results"
mkdir -p "$OUTDIR"

# 2. Concatonate Independent Replicas
gmx trjcat -f step5_*.xtc -o "$OUTDIR/combined.xtc"

# 3. RMSD and RMSD Matrix
echo 4 4 | gmx_mpi rms -s step5_1.tpr -f "$OUTDIR/combined.xtc" -o "$OUTDIR/rmsd.xvg"
echo 4 4 | gmx_mpi rms -s step5_1.tpr -f "$OUTDIR/combined.xtc" -m "$OUTDIR/rmsd-matrix.xpm"
echo "Converting to .txt" | gmx_mpi xpm2mat -f "$OUTDIR/rmsd-matrix.xpm" -o "$OUTDIR/rmsd-matrix.txt"


# 4. RMSF
echo 4 | gmx_mpi rmsf -s step5_1.tpr -f "$OUTDIR/combined.xtc" -o "$OUTDIR/rmsf.xvg" -res

# 5. Radius of Gyration
echo 1 | gmx_mpi gyrate -s step5_1.tpr -f "$OUTDIR/combined.xtc" -o "$OUTDIR/gyrate.xvg"

# 6. PCA
echo 4 4 | gmx_mpi covar -s step5_1.tpr -f "$OUTDIR/combined.xtc" -o "$OUTDIR/eigenval.xvg" -v "$OUTDIR/eigenvec.trr"
echo 4 4 | gmx_mpi anaeig -v "$OUTDIR/eigenvec.trr" -s step5_1.tpr -f "$OUTDIR/combined.xtc" -first 1 -last 10 -proj "$OUTDIR/proj.xvg"

# 7. Energy
gmx_mpi eneconv -f step5_*.edr -o "$OUTDIR/combined.edr"

echo "Potential" | gmx_mpi energy -f "$OUTDIR/combined.edr" -o "$OUTDIR/kbT_scalar.xvg"

# Remove Cheackpoint Files
# find . -type f -name '*#*' -exec rm -i {} \;

echo "All MD replicas analyzed. Results in md_analysis_results/"
