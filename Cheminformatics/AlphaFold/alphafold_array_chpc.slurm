#!/bin/bash
#SBATCH --job-name=alphafold_array
#SBATCH -p notchpeak-gpu
#SBATCH -A notchpeak-gpu
#SBATCH --gres=gpu:A100:1          # Request 1 GPU per job
#SBATCH --cpus-per-task=8          # Allocate 8 CPU cores per job
#SBATCH --mem=32G                  # Allocate 32GB RAM per job
#SBATCH --time=6:00:00             # Max execution time per job
#SBATCH --array=0-15999%50         # Submit 16,000 jobs, running max 50 at a time
#SBATCH --output=logs/af_job_%A_%a.out
#SBATCH --error=logs/af_job_%A_%a.err

# Debugging: Print SLURM task ID
echo "Running SLURM Array Job ID: $SLURM_ARRAY_TASK_ID"

# Load required modules
ml purge
ml alphafold/3.0.0
ml cuda/11.7  # Reload CUDA to ensure compatibility

# Define directories
INPUT_JSONS_DIR="/uufs/chpc.utah.edu/common/home/cheatham-group1/nate/diehl/biosensors/json"
OUTPUT_BASE_DIR="/uufs/chpc.utah.edu/common/home/cheatham-group1/nate/diehl/biosensors/struct"

# Ensure the JSON directory exists and is not empty
if [ ! -d "$INPUT_JSONS_DIR" ] || [ -z "$(ls -A $INPUT_JSONS_DIR/*.json 2>/dev/null)" ]; then
    echo "Error: JSON input directory is empty or does not exist."
    exit 1
fi

# Get the list of JSON files
mapfile -t JSON_FILES < <(find "$INPUT_JSONS_DIR" -maxdepth 1 -name "*.json" | sort)
NUM_FILES=${#JSON_FILES[@]}

# Check if SLURM task ID is within range
if [[ $SLURM_ARRAY_TASK_ID -ge $NUM_FILES ]]; then
    echo "Error: SLURM task ID ($SLURM_ARRAY_TASK_ID) exceeds available JSON files ($NUM_FILES)."
    exit 1
fi

# Get the corresponding JSON file
JSON_FILE=${JSON_FILES[$SLURM_ARRAY_TASK_ID]}

# Extract filename without extension for unique output directory
BASENAME=$(basename "$JSON_FILE" .json)
OUTPUT_DIR="$OUTPUT_BASE_DIR/$BASENAME"

# Define the expected output files that indicate completion
EXPECTED_FILES=("ranked_0.pdb" "features.pkl" "result_model_1.pkl")

# Check if job is already completed
COMPLETED=1
for file in "${EXPECTED_FILES[@]}"; do
    if [[ ! -f "$OUTPUT_DIR/$file" ]]; then
        COMPLETED=0
        break
    fi
done

if [[ $COMPLETED -eq 1 ]]; then
    echo "Job for $BASENAME already completed. Skipping."
    exit 0
fi

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Debugging: Print job-specific information
echo "Processing JSON File: $JSON_FILE"
echo "Output Directory: $OUTPUT_DIR"

# Ensure correct GPU allocation
export CUDA_VISIBLE_DEVICES=0

# Run AlphaFold inside Apptainer container
srun apptainer exec --nv /uufs/chpc.utah.edu/sys/installdir/alphafold/3.0.0/alphafold-3.0.0.sif \
    python /app/alphafold/run_alphafold.py \
    --db_dir=/scratch/general/vast/app-repo/alphafold3 \
    --model_dir=/scratch/general/vast/app-repo/alphafold3 \
    --jax_compilation_cache_dir=$HOME/.cache \
    --input_json "$JSON_FILE" --output_dir "$OUTPUT_DIR"
