#!/usr/bin/env bash
set -euo pipefail

# ---------------- USER SETTINGS ----------------
ROOT="/scratch/rai/vast1/u1116818/kash/dock/decoys"
MANIFEST="/scratch/rai/vast1/u1116818/kash/dock/decoys/pdbqt_out/pdbqt_manifest.csv"
OUTLIST="/scratch/rai/vast1/u1116818/kash/dock/decoys/ligands.list"
TYPELIST="/scratch/rai/vast1/u1116818/kash/dock/decoys/ligand_types.unique.txt"
# ------------------------------------------------

if [[ ! -f "$MANIFEST" ]]; then
  echo "ERROR: Manifest file not found: $MANIFEST" >&2
  exit 1
fi
if [[ ! -d "$ROOT" ]]; then
  echo "ERROR: ROOT directory not found: $ROOT" >&2
  exit 1
fi

# ---- build ligand list from manifest ----
awk -F',' '
NR==1 {
  for (i=1; i<=NF; i++) {
    tmp=$i
    gsub(/\r|"/, "", tmp)
    if (tmp == "pdbqt") { col = i; break }
  }
  if (!col) {
    print "ERROR: pdbqt column not found in manifest" > "/dev/stderr"
    exit 1
  }
  next
}
$col != "" { print $col }
' "$MANIFEST" \
| sed 's/\r$//' \
| awk -v root="$ROOT" '{
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1)
    if ($1 ~ /^\//) print $1
    else {
      r = root
      sub(/\/+$/, "", r)
      print r "/" $1
    }
  }' \
> "$OUTLIST"

N=$(wc -l < "$OUTLIST" | awk '{print $1}')
echo "Wrote $N ligand paths to $OUTLIST"

# ---- existence check ----
MISSING=0
while read -r f; do
  [[ -z "${f:-}" ]] && continue
  if [[ ! -f "$f" ]]; then
    echo "WARNING: Missing ligand file: $f" >&2
    ((MISSING++))
  fi
done < "$OUTLIST"

if (( MISSING > 0 )); then
  echo "WARNING: $MISSING ligand files listed but not found" >&2
else
  echo "All ligand files exist ✔"
fi

# ---- extract unique AutoDock atom types ----
echo "Extracting unique AutoDock atom types..."

awk '
$1=="ATOM" || $1=="HETATM" { print $NF }
' $(cat "$OUTLIST") \
| sort -u \
> "$TYPELIST"

NTYPES=$(wc -l < "$TYPELIST" | awk '{print $1}')
echo "Wrote $NTYPES unique ligand atom types to $TYPELIST"

echo "Suggested GPF ligand_types line:"
echo -n "ligand_types "
tr '\n' ' ' < "$TYPELIST"
echo
