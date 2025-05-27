#!/bin/bash
#SBATCH --job-name=gmx_md
#SBATCH --output=remd_job_%j.out
#SBATCH --error=remd_job_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=36   # One per replica
#SBATCH --cpus-per-task=1

# Load GROMACS module
module purge
module load openmpi
source ~/gromacs_no_simd/bin/GMXRC

# Variables
init="step3_input"
mini_prefix="step4.0_minimization"
equi_prefix="step4.1_equilibration"
prod_prefix="step5_production"
prod_step="step5"

# ------------------------
# Step 1: Minimization
# ------------------------
echo "Starting energy minimization..."
gmx grompp -f ${mini_prefix}.mdp -o ${mini_prefix}.tpr -c ${init}.gro -r ${init}.gro -p topol.top -n index.ndx -maxwarn -1
gmx_d mdrun -v -deffnm ${mini_prefix}

# ------------------------
# Step 2: Equilibration
# ------------------------
echo "Starting equilibration..."
gmx grompp -f ${equi_prefix}.mdp -o ${equi_prefix}.tpr -c ${mini_prefix}.gro -r ${init}.gro -p topol.top -n index.ndx
gmx mdrun -v -deffnm ${equi_prefix}

# ------------------------
# Step 3: Production (10 segments)
# ------------------------
echo "Starting production MD..."
cnt=1
cntmax=10

while [ $cnt -le $cntmax ]; do
    pcnt=$((cnt - 1))
    istep="${prod_step}_${cnt}"
    
    if [ $cnt -eq 1 ]; then
        # First chunk starts from equilibration output
        pstep=$equi_prefix
        gmx grompp -f ${prod_prefix}.mdp \
                   -o ${istep}.tpr \
                   -c ${pstep}.gro \
                   -p topol.top \
                   -n index.ndx \
                   -maxwarn 1
    else
        # Subsequent chunks use checkpoint to continue
        pstep="${prod_step}_${pcnt}"
        gmx grompp -f ${prod_prefix}.mdp \
                   -o ${istep}.tpr \
                   -c ${pstep}.gro \
                   -t ${pstep}.cpt \
                   -p topol.top \
                   -n index.ndx \
                   -maxwarn 1
    fi

    gmx mdrun -v -deffnm ${istep}
    cnt=$((cnt + 1))
done

echo "All GROMACS steps completed."
