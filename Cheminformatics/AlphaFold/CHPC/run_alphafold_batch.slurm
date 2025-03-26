#!/bin/bash
#SBATCH --job-name=alphafold_batch
#SBATCH -p notchpeak-gpu
#SBATCH -A notchpeak-gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=logs/af_batch_%A_%a.out
#SBATCH --error=logs/af_batch_%A_%a.err
#SBATCH --array=1-100

# Debugging: Print SLURM task ID
echo "Running SLURM Array Job ID: $SLURM_ARRAY_TASK_ID"

# Load required modules
ml purge
ml alphafold/3.0.0
ml cuda/11.7  # Reload CUDA to ensure compatibility

# Define directories
INPUT_JSONS_DIR="/uufs/chpc.utah.edu/common/home/cheatham-group1/nate/diehl/biosensors/json"
OUTPUT_BASE_DIR="/uufs/chpc.utah.edu/common/home/cheatham-group1/nate/diehl/biosensors/struct"

# Get the list of JSON files
mapfile -t JSON_FILES < <(find "$INPUT_JSONS_DIR" -maxdepth 1 -name "*.json" | sort)
NUM_FILES=${#JSON_FILES[@]}

# Convert SLURM 1-based indexing to 0-based array indexing
ARRAY_INDEX=$((SLURM_ARRAY_TASK_ID - 1))

# Check if SLURM task ID is within range
if [[ $ARRAY_INDEX -ge $NUM_FILES ]]; then
    echo "Error: SLURM task ID ($SLURM_ARRAY_TASK_ID) exceeds available JSON files ($NUM_FILES)."
    exit 1
fi

# Get the corresponding JSON file
JSON_FILE=${JSON_FILES[$ARRAY_INDEX]}

# Extract filename for output directory
BASENAME=$(basename "$JSON_FILE" .json)
OUTPUT_DIR="$OUTPUT_BASE_DIR/$BASENAME"

# Skip if already completed
EXPECTED_FILES=("ranked_0.pdb" "features.pkl" "result_model_1.pkl")
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

mkdir -p "$OUTPUT_DIR"

echo "Array Index: $ARRAY_INDEX"
echo "Total JSON Files: $NUM_FILES"
echo "Selected JSON File: $JSON_FILE"
echo "Output will be written to: $OUTPUT_DIR"

# Print a few lines of the JSON file to confirm contents
echo "----- JSON Preview -----"
head -n 20 "$JSON_FILE"
echo "------------------------"

# Ensure container exists
if [[ ! -f /uufs/chpc.utah.edu/sys/installdir/alphafold/3.0.0/alphafold-3.0.0.sif ]]; then
    echo "Error: Apptainer container not found!"
    exit 1
fi

# Run AlphaFold inside Apptainer container
srun apptainer exec --nv /uufs/chpc.utah.edu/sys/installdir/alphafold/3.0.0/alphafold-3.0.0.sif \
    python /app/alphafold/run_alphafold.py \
    --db_dir=/scratch/general/vast/app-repo/alphafold3 \
    --model_dir=/scratch/general/vast/app-repo/alphafold3 \
    --jax_compilation_cache_dir=$HOME/.cache \
    --json_path "$JSON_FILE" --output_dir "$OUTPUT_DIR"
