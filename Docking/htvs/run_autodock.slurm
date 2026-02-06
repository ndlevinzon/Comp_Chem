#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=4
#SBATCH --job-name=kash_htvs
#SBATCH --gres=gpu:1
#SBATCH --output=slurm/SLURM-%A_%a.out-%N
#SBATCH --error=slurm/SLURM-%A_%a.err-%N
#SBATCH --account=cheatham
#SBATCH --partition=granite-gpu
#SBATCH --qos=granite-gpu

set -eo pipefail
echo "... Job started at $(date)"

# Ensure slurm/ exists BEFORE Slurm tries to write into it
mkdir -p slurm

module purge
module load miniforge3/25.11.0

# Robust conda activation in batch shells
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /uufs/chpc.utah.edu/common/home/u1116818/.conda/envs/htvs

# ---- USER SETTINGS (can be overridden by --export in sbatch) ----
LIGLIST="${LIGLIST:-/scratch/rai/vast1/u1116818/kash/dock/decoys/ligands.list}"
TODO_FILE="${TODO_FILE:-}"   # required when using TODO mapping
GRIDFLD="${GRIDFLD:-/scratch/rai/vast1/u1116818/kash/dock/receptor/pore/pore_adt.maps.fld}"
RESULTS_DIR="${RESULTS_DIR:-/scratch/rai/vast1/u1116818/kash/dock/results/pore}"
NRUN="${NRUN:-50}"
ADGPU_BIN="${ADGPU_BIN:-autodock_gpu_128wi}"
# ---------------------------------------------------------------

mkdir -p "$RESULTS_DIR"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID is unset. Submit with: sbatch --array=1-N $0" >&2
  exit 1
fi

# --- NEW: map array task -> ligands.list line number via TODO_FILE ---
if [[ -z "${TODO_FILE:-}" ]]; then
  echo "ERROR: TODO_FILE is not set. Submit with: sbatch --export=ALL,TODO_FILE=/path/to/todo.txt ..." >&2
  exit 1
fi
[[ -f "$TODO_FILE" ]] || { echo "ERROR: TODO_FILE not found: $TODO_FILE" >&2; exit 1; }
[[ -f "$LIGLIST"   ]] || { echo "ERROR: LIGLIST not found: $LIGLIST" >&2; exit 1; }

LINE_NUM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TODO_FILE" | tr -d '[:space:]\r')
[[ -n "${LINE_NUM:-}" ]] || { echo "ERROR: No LINE_NUM found in TODO_FILE for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
[[ "$LINE_NUM" =~ ^[0-9]+$ ]] || { echo "ERROR: Invalid LINE_NUM '${LINE_NUM}' (task ${SLURM_ARRAY_TASK_ID})" >&2; exit 1; }

LIG=$(sed -n "${LINE_NUM}p" "$LIGLIST" | tr -d '\r')
[[ -n "${LIG:-}" ]] || { echo "ERROR: No ligand found in LIGLIST at line ${LINE_NUM} (task ${SLURM_ARRAY_TASK_ID})" >&2; exit 1; }
# -------------------------------------------------------------------

[[ -f "$LIG" ]] || { echo "ERROR: Ligand file not found: $LIG" >&2; exit 1; }
[[ -f "$GRIDFLD" ]] || { echo "ERROR: Grid .fld not found: $GRIDFLD" >&2; exit 1; }

BASE=$(basename "$LIG" .pdbqt)

echo "Task: ${SLURM_ARRAY_TASK_ID} -> ligands.list line: ${LINE_NUM}"
echo "Docking ligand: $LIG"
echo "Results base: $BASE"
echo "Results dir: $RESULTS_DIR"

LIG_OUTDIR="${RESULTS_DIR}/${BASE}"
mkdir -p "$LIG_OUTDIR"
cd "$LIG_OUTDIR"

# Optional: skip if already has dlg(s) (helps if requeued/retried)
if ls -1 *.dlg >/dev/null 2>&1; then
  echo "SKIP: .dlg already exists in $LIG_OUTDIR"
  exit 0
fi

$ADGPU_BIN \
  -ffile "$GRIDFLD" \
  -lfile "$LIG" \
  -nrun "$NRUN" \
  --resnam "$BASE" \
  --dlgoutput 1 \
  --xmloutput 0 \
  > "${BASE}.stdout.txt" \
  2> "${BASE}.stderr.txt"

# quick check
ls -1 *.dlg 2>/dev/null || echo "WARNING: No .dlg found in $LIG_OUTDIR (check stdout/stderr)"
echo "Done: ${BASE}"
