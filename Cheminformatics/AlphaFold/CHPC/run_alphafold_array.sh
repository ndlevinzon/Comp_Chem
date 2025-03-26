#!/bin/bash

# SLURM's max array size
MAX_ARRAY_SIZE=100
TOTAL_JOBS=16000
JOBS_AT_A_TIME=1

# Load required modules
ml purge
ml alphafold/3.0.0
ml cuda/11.7

# Define directories
INPUT_JSONS_DIR="/uufs/chpc.utah.edu/common/home/cheatham-group1/nate/diehl/biosensors/json"
OUTPUT_BASE_DIR="/uufs/chpc.utah.edu/common/home/cheatham-group1/nate/diehl/biosensors/struct"

# Debugging: Check if JSON directory exists and print its path
if ! find "/uufs/chpc.utah.edu/common/home/cheatham-group1/nate/diehl/biosensors" -maxdepth 1 -type d -name "json" | grep -q .; then
    echo "Error: JSON input directory does not exist!"
    echo "Checked under: /uufs/chpc.utah.edu/common/home/cheatham-group1/nate/diehl/biosensors"
    exit 1
fi

# Detect JSON directory dynamically
INPUT_JSONS_DIR=$(find "/uufs/chpc.utah.edu/common/home/cheatham-group1/nate/diehl/biosensors" -maxdepth 1 -type d -name "json" | head -n 1)

# Debugging: Print detected directory
echo "Detected JSON directory: $INPUT_JSONS_DIR"

# Ensure the directory contains JSON files
if ! find "$INPUT_JSONS_DIR" -maxdepth 1 -type f -name "*.json" | grep -q .; then
    echo "Error: No JSON files found in detected directory: $INPUT_JSONS_DIR"
    ls -lah "$INPUT_JSONS_DIR"
    exit 1
fi

# Debugging: List a few JSON files
echo "Found JSON files:"
find "$INPUT_JSONS_DIR" -maxdepth 1 -type f -name "*.json" | head -n 10

# Get the list of JSON files
mapfile -t JSON_FILES < <(find "$INPUT_JSONS_DIR" -maxdepth 1 -name "*.json" | sort)
NUM_FILES=${#JSON_FILES[@]}


for (( START=1; START<=TOTAL_JOBS; START+=MAX_ARRAY_SIZE )); do
    END=$((START + MAX_ARRAY_SIZE - 1))
    if [[ $END -gt $TOTAL_JOBS ]]; then
        END=$TOTAL_JOBS
    fi

    echo "Submitting job array: $START-$END"

    sbatch --array=$START-$END%$JOBS_AT_A_TIME run_alphafold_batch.slurm
done

