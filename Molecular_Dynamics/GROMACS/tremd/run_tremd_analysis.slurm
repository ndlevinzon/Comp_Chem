#!/bin/bash
#SBATCH --job-name=remd_analysis
#SBATCH --output=slurm/remd_analysis_%j.out
#SBATCH --error=slurm/remd_analysis_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00

# A Decent Tutorial on TREMD with GROMACS: https://www.researchgate.net/publication/335453547_Tutorial_on_setting_up_the_Replica_Exchange_Molecular_Dynamics_REMD_simulation_using_alanine_dipeptide_as_a_toy_model #

# Load GROMACS module
module purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 gromacs/2024.4

# ----------- SETUP -------------------

REPLICAS=36
DEFFNM="remd"  # Trajectory/log base name

# ----------- STEP 0: DEMUX INDEX ---------------

echo "Generating replica index..."

mkdir -p remd_analysis_results

# Create a temporary directory for logs
LOG_TMP="remd_logs"
mkdir -p "$LOG_TMP"

# Copy remd.log from each replica directory and rename it numerically
for i in $(seq 0 $((REPLICAS - 1))); do
    cp $i/remd.log "$LOG_TMP/remd_$i.log"
done

# Concatenate all logs into a single file
cat "$LOG_TMP"/remd_*.log > "$LOG_TMP/REMD.log"

# De-multiplexing the logs, move them to analysis directory
demux.pl "$LOG_TMP/REMD.log" > remd_analysis_results/remd_replica_index.log
mv replica_index.xvg remd_analysis_results/
mv replica_temp.xvg remd_analysis_results/

# De-multiplexing the trajectories
echo "Demultiplexing trajectories..."

# Build list of input .xtc files
XTCLIST=()
for i in $(seq 0 $((REPLICAS - 1))); do
    XTCLIST+=("$i/remd.xtc")
done

# Uncomment if downsampling is required
#awk '($1 % 2000 == 0)' remd_analysis_results/replica_index.xvg > remd_analysis_results/replica_index_filtered.xvg

# Make demuxed .xtc files
gmx_mpi trjcat -f "${XTCLIST[@]}" -demux remd_analysis_results/replica_index.xvg

# ----------- STEP 1: LOOP THROUGH REPLICAS ---------------

for TARGET_REPLICA in $(seq 0 $((REPLICAS - 1))); do
    echo "Processing replica $TARGET_REPLICA..."

    OUTDIR="remd_analysis_results/replica_${TARGET_REPLICA}"
    mkdir -p "$OUTDIR"

    # Move .xtc to the correct analysis directory
    mv ${TARGET_REPLICA}_trajout.xtc ${OUTDIR}

    DEMUX_TRAJ="$OUTDIR/${TARGET_REPLICA}_trajout.xtc"
    DEMUX_TPR="$OUTDIR/demuxed_${TARGET_REPLICA}.tpr"

    # 2. Copy matching .tpr
    cp ${TARGET_REPLICA}/${DEFFNM}.tpr "$DEMUX_TPR"

    # 3. RMSD
    echo 4 4 | gmx_mpi rms -s "$DEMUX_TPR" -f "$DEMUX_TRAJ" -o "$OUTDIR/rmsd.xvg"

    # 4. RMSF
    echo 4 | gmx_mpi rmsf -s "$DEMUX_TPR" -f "$DEMUX_TRAJ" -o "$OUTDIR/rmsf.xvg" -res

    # 5. Radius of Gyration
    echo 1 | gmx_mpi gyrate -s "$DEMUX_TPR" -f "$DEMUX_TRAJ" -o "$OUTDIR/gyrate.xvg"

done

# ----------- STEP 2: ANALYZE COMBINED TRAJECTORY ---------------

# 6. PCA
echo 4 4 | gmx_mpi covar -s step4.1_equilibration.tpr -f remd_analysis_results/demuxed_0.xtc -n index.ndx -o remd_analysis_results/eigenval.xvg -v remd_analysis_results/eigenvec.trr
echo 4 4 | gmx_mpi anaeig -v remd_analysis_results/eigenvec.trr -s step4.1_equilibration.tpr -f remd_analysis_results/demuxed_0.xtc -first 1 -last 10 -proj remd_analysis_results/proj.xvg

# 7. Energy
# Build list of input .edr files
XTCLIST=()
for i in $(seq 0 $((REPLICAS - 1))); do
    XTCLIST+=("$i/remd.edr")
done
gmx_mpi eneconv -f "${XTCLIST[@]}" -o remd_analysis_results/combined.edr

echo "Potential" | gmx_mpi energy -f remd_analysis_results/combined.edr -o remd_analysis_results/kbT_scalar.xvg


# Remove Cheackpoint Files
# find . -type f -name '*#*' -exec rm -i {} \;

echo "All REMD replicas analyzed. Results in remd_analysis_results/"
