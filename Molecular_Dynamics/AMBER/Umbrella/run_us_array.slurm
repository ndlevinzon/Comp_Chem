#!/bin/bash
#SBATCH --job-name=horvath_um
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --gres=gpu:2080ti:1
#SBATCH --account=cheatham-gpu-np
#SBATCH --partition=cheatham-gpu-np
#SBATCH --output=slurm/MD.%A_%a.out
#SBATCH --error=slurm/MD.%A_%a.err

set -euo pipefail
set -x

module purge
module load cuda/12.8.1
module load gcc/11.2.0
module load openmpi/5.0.3

export AMBERHOME="/uufs/chpc.utah.edu/sys/installdir/r8/amber/24"
ENGINE="$AMBERHOME/bin/pmemd.cuda"

cd "$SLURM_SUBMIT_DIR"
mkdir -p slurm

PRMTOP_REL="../../../build/repartitioned.topo"
START_RST_REL="../../min/step9.rst7"
REF_RST_REL="../../min/step9.rst7"

WINLIST="windows.list"

# ---------- Parent: build list + submit array ----------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "Parent job -> regenerating ${WINLIST} from chi_* dirs in: $(pwd)"

  shopt -s nullglob
  dirs=(chi_*)
  shopt -u nullglob

  if (( ${#dirs[@]} == 0 )); then
    echo "ERROR: No chi_* directories found in $(pwd)"
    exit 1
  fi

  IFS=$'\n' dirs_sorted=($(printf "%s\n" "${dirs[@]}" | sort))
  unset IFS

  : > "$WINLIST"
  idx=0
  for d in "${dirs_sorted[@]}"; do
    [[ -d "$d" ]] || continue
    tag="${d#chi_}"
    printf "%d\t%s\t%s\n" "$idx" "$d" "$tag" >> "$WINLIST"
    ((++idx))
  done

  NWIN=$(wc -l < "$WINLIST" | awk '{print $1}')
  echo "Wrote ${WINLIST} with ${NWIN} windows."
  sbatch --array=0-$((NWIN-1)) "$0"
  exit 0
fi

# ---------- Array task: run one window ----------
[[ -f "$WINLIST" ]] || { echo "ERROR: Missing $WINLIST"; exit 1; }

NWIN=$(wc -l < "$WINLIST" | awk '{print $1}')
if (( SLURM_ARRAY_TASK_ID >= NWIN )); then
  echo "Task ${SLURM_ARRAY_TASK_ID} >= NWIN ${NWIN}. Exiting."
  exit 0
fi

line=$(awk -v id="${SLURM_ARRAY_TASK_ID}" 'NR==id+1{print}' "$WINLIST")
wdir=$(echo "$line" | awk '{print $2}')
tag=$(echo "$line"  | awk '{print $3}')

echo "Running window: task=${SLURM_ARRAY_TASK_ID}, dir=${wdir}, tag=${tag}"
cd "$wdir"

[[ -f "$PRMTOP_REL"    ]] || { echo "ERROR: Missing $PRMTOP_REL (from $(pwd))"; exit 1; }
[[ -f "$START_RST_REL" ]] || { echo "ERROR: Missing $START_RST_REL (from $(pwd))"; exit 1; }
[[ -f "$REF_RST_REL"   ]] || { echo "ERROR: Missing $REF_RST_REL (from $(pwd))"; exit 1; }

# If production already exists, skip
if [[ -f "prod_${tag}.nc" ]]; then
  echo "prod_${tag}.nc exists -> skipping window ${tag}"
  exit 0
fi

# 1) Seed minimization (creates seed_min_TAG.rst7)
if [[ ! -f "seed_min_${tag}.rst7" ]]; then
  echo "Seeding MIN..."
  "$ENGINE" -O \
    -p "$PRMTOP_REL" \
    -c "$START_RST_REL" \
    -ref "$REF_RST_REL" \
    -i "seed_min_${tag}.mdin" \
    -o "seed_min_${tag}.out" \
    -r "seed_min_${tag}.rst7"
else
  echo "seed_min_${tag}.rst7 exists -> skipping seed MIN"
fi

# 2) Seed MD (creates seed_md_TAG.rst7 and seed_md_TAG.nc)
SEED_MD_IN="seed_min_${tag}.rst7"
if [[ -f "seed_md_${tag}.rst7" ]]; then
  echo "seed_md_${tag}.rst7 exists -> skipping seed MD"
else
  echo "Seeding MD..."
  "$ENGINE" -O \
    -p "$PRMTOP_REL" \
    -c "$SEED_MD_IN" \
    -ref "$REF_RST_REL" \
    -i "seed_md_${tag}.mdin" \
    -o "seed_md_${tag}.out" \
    -r "seed_md_${tag}.rst7" \
    -x "seed_md_${tag}.nc"
fi

# 3) Production umbrella
PROD_IN="seed_md_${tag}.rst7"
if [[ ! -f "$PROD_IN" ]]; then
  echo "ERROR: Missing $PROD_IN (seed MD did not produce restart?)" >&2
  exit 1
fi

echo "Production umbrella..."
"$ENGINE" -O \
  -p "$PRMTOP_REL" \
  -c "$PROD_IN" \
  -ref "$REF_RST_REL" \
  -i "prod_${tag}.mdin" \
  -o "prod_${tag}.out" \
  -r "prod_${tag}.rst7" \
  -x "prod_${tag}.nc"

echo "Done window ${tag}"
