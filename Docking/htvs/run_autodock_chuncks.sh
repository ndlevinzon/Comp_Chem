#!/usr/bin/env bash
set -euo pipefail

LIGLIST="/scratch/rai/vast1/u1116818/kash/dock/decoys/ligands.list"
JOBSCRIPT="run_autodock.slurm"
RESULTS_DIR="/scratch/rai/vast1/u1116818/kash/dock/results/pore"

MAX_CONCURRENT=50
OUTPUT_LAYOUT="perdir"   # "perdir" or "flat"

MAX_ARRAY=999  # submit arrays as 1..N where N<=999

[[ -f "$LIGLIST" ]] || { echo "ERROR: ligands.list not found: $LIGLIST" >&2; exit 1; }
[[ -d "$RESULTS_DIR" ]] || { echo "ERROR: RESULTS_DIR not found: $RESULTS_DIR" >&2; exit 1; }
[[ -f "$JOBSCRIPT" ]] || { echo "ERROR: JOBSCRIPT not found in current dir: $JOBSCRIPT" >&2; exit 1; }

N=$(wc -l < "$LIGLIST" | tr -d '[:space:]\r')
[[ -n "$N" && "$N" -gt 0 ]] || { echo "ERROR: ligands.list appears empty (N=$N)" >&2; exit 1; }

echo "Total ligands in list: $N"
echo "Results dir: $RESULTS_DIR"
echo "Output layout: $OUTPUT_LAYOUT"

TMP_IDS="$(mktemp)"
trap 'rm -f "$TMP_IDS"' EXIT

done_count=0
todo_count=0

for i in $(seq 1 "$N"); do
  LIG=$(sed -n "${i}p" "$LIGLIST" | tr -d '\r')
  [[ -n "${LIG:-}" ]] || continue

  BASE=$(basename "$LIG" .pdbqt)

  is_done=0
  if [[ "$OUTPUT_LAYOUT" == "perdir" ]]; then
    if compgen -G "${RESULTS_DIR}/${BASE}/*.dlg" > /dev/null; then
      is_done=1
    fi
  else
    if compgen -G "${RESULTS_DIR}/${BASE}*.dlg" > /dev/null; then
      is_done=1
    fi
  fi

  if (( is_done )); then
    done_count=$((done_count + 1))
  else
    echo "$i" >> "$TMP_IDS"
    todo_count=$((todo_count + 1))
  fi
done

echo "Already completed: $done_count"
echo "To run: $todo_count"

if (( todo_count == 0 )); then
  echo "Nothing to submit. All ligands appear completed."
  exit 0
fi

# Stable master TODO file (for record)
MASTER_TODO="${RESULTS_DIR}/todo_ids.$(date +%Y%m%d_%H%M%S).txt"
cp "$TMP_IDS" "$MASTER_TODO"
echo "Wrote master todo ids: $MASTER_TODO"

# Split into chunk TODO files, chain arrays with dependencies
prev_jobid=""

start=1
chunk_idx=0
while (( start <= todo_count )); do
  end=$(( start + MAX_ARRAY - 1 ))
  (( end > todo_count )) && end=$todo_count
  chunk_len=$(( end - start + 1 ))
  chunk_idx=$(( chunk_idx + 1 ))

  CHUNK_TODO="${RESULTS_DIR}/todo_chunk_${chunk_idx}.$(basename "$MASTER_TODO")"
  sed -n "${start},${end}p" "$MASTER_TODO" > "$CHUNK_TODO"

  echo "Chunk ${chunk_idx}: TODO lines ${start}-${end} -> ${CHUNK_TODO} (len=${chunk_len})"

  sbatch_args=(
    --array="1-${chunk_len}%${MAX_CONCURRENT}"
    --export=ALL,TODO_FILE="$CHUNK_TODO",LIGLIST="$LIGLIST",OUTPUT_LAYOUT="$OUTPUT_LAYOUT",RESULTS_DIR="$RESULTS_DIR"
  )

  # Chain submission so only one array is queued/running at a time
  if [[ -n "$prev_jobid" ]]; then
    # Use afterany so the chain continues even if some tasks fail
    sbatch_args+=( --dependency="afterany:${prev_jobid}" )
  fi

  # Submit and capture jobid
  submit_out=$(sbatch "${sbatch_args[@]}" "$JOBSCRIPT")
  jobid=$(echo "$submit_out" | awk '{print $4}')
  echo "Submitted chunk ${chunk_idx} as job ${jobid}"
  prev_jobid="$jobid"

  start=$(( end + 1 ))
done

echo "All chunks chained. Final job in chain: ${prev_jobid}"
