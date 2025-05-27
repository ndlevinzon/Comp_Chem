#!/bin/bash
#SBATCH --job-name=gmx_remd
#SBATCH --output=slurm/remd_job_%j.out
#SBATCH --error=slurm/remd_job_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=36   # One per replica
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --exclusive 
#SBATCH --mem=0

# Load GROMACS module
module purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 gromacs/2024.4

# Variables
init="step3_input"
mini_prefix="step4.0_minimization"
equi_prefix="step4.1_equilibration"
topol="topol.top"
index="index.ndx"
temperatures=(300.00 301.37 302.74 304.12 305.51 306.89 308.29 309.69 311.10 312.51 313.92 315.34 316.77 318.20 319.63 321.07 322.52 323.97 325.43 326.89 328.36 329.83 331.31 332.79 334.28 335.77 337.27 338.77 340.28 341.80 343.32 344.84 346.39 347.93 349.47 350.00)  # 36 replicas

# ------------------------
# Step 1: Minimization
# ------------------------
# echo "Starting energy minimization..."
# gmx grompp -f ${mini_prefix}.mdp -o ${mini_prefix}.tpr -c ${init}.gro -r ${init}.gro -p ${topol} -n ${index} -maxwarn 1
# gmx mdrun -v -deffnm ${mini_prefix}

# # ------------------------
# # Step 2: Equilibration
# # ------------------------
# echo "Starting equilibration..."
# gmx grompp -f ${equi_prefix}.mdp -o ${equi_prefix}.tpr -c ${mini_prefix}.gro -r ${init}.gro -p ${topol} -n ${index}
# gmx mdrun -v -deffnm ${equi_prefix}

# # ------------------------
# # Step 3: Prepare REMD Directories
# # ------------------------
# echo "Preparing REMD directories..."
# for i in ${!temperatures[@]}; do
#     temp=${temperatures[$i]}
#     dir=$i
#     mkdir -p $dir
#     cp ${topol} ${index} ${equi_prefix}.gro $dir/

#     # Generate MDP for this temperature
#     sed "s/XXX/${temp}/g" remd_template.mdp > $dir/remd.mdp
    
#     # grompp for each replica
#     pushd $dir
#     gmx grompp -f remd.mdp -c ${equi_prefix}.gro -p ${topol} -n ${index} -o remd.tpr -maxwarn 1
#     popd

# done

# ------------------------
# Step 4: Run T-REMD with checkpointing
# ------------------------
echo "Launching T-REMD with ${#temperatures[@]} replicas..."
timeout 3300 \
mpirun -np ${#temperatures[@]} gmx_mpi mdrun -multidir $(seq 0 $((${#temperatures[@]} - 1))) \
    -replex 1000 \
    -deffnm remd \
    -ntomp 1 \
    -cpt 5 \
    -cpi remd.cpt

# ------------------------
# Step 5: Resubmit Job
# ------------------------
echo "Job ran for 55 minutes. Resubmitting self..."
sbatch run_gmx_tremd_beluga.slurm
