#!/bin/bash

#SBATCH --job-name=vina_docking
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=22000
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=cheatham
#SBATCH --partition=lonepeak

# Directory paths
receptor_file="6IRT_clean.pdbqt"
ligands_dir="ligands/"
output_dir="output/"
job_scripts_dir="job_scripts/"

# Ensure output and job scripts directories exist
mkdir -p $output_dir
mkdir -p $job_scripts_dir

# Center and size of the grid box
center_x=144.498
center_y=142.652
center_z=134.552
size_x=40
size_y=40
size_z=40

# Initialize the dependency variable
prev_job_id=""

# Loop through each ligand file
for ligand_file in $ligands_dir/*.pdbqt; do
    # Extract the base name of the ligand file
    ligand_base=$(basename $ligand_file)
    output_file="${output_dir}/out_${ligand_base}"
    config_file="${job_scripts_dir}/config_${ligand_base}.txt"
    job_file="${job_scripts_dir}/job_${ligand_base}.sh"

    # Create a specific config file for each ligand
    cat <<EOL > $config_file
receptor = $receptor_file
ligand = $ligand_file
center_x = $center_x
center_y = $center_y
center_z = $center_z
size_x = $size_x
size_y = $size_y
size_z = $size_z
out = $output_file
EOL

    # Create a SLURM job script for each ligand
    cat <<EOL > $job_file
#!/bin/bash
#SBATCH --job-name=vina_${ligand_base}
#SBATCH --output=${output_dir}/slurm-%j.out
#SBATCH --error=${output_dir}/slurm-%j.err
#SBATCH --partition=lonepeak
#SBATCH --account=cheatham
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=4G

module load autodock-vina/1.2.5
vina --config $config_file
EOL

# Submit the job script to SLURM with dependency if prev_job_id is set
if [ -z "$prev_job_id" ]; then
    # Submit the first job without dependency
    job_id=$(sbatch $job_file | awk '{print $4}')
else
    # Submit subsequent jobs with dependency
    job_id=$(sbatch --dependency=afterok:$prev_job_id $job_file | awk '{print $4}')
fi

# Update the prev_job_id for the next iteration
prev_job_id=$job_id
done
