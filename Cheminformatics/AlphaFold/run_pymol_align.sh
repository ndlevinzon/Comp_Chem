#!/bin/bash

# Load PyMol for "align" utility
module load pymol/2.5.7

# Set root directory containing all subdirectories and reference mobile PDB
ROOT_DIR="/uufs/chpc.utah.edu/common/home/cheatham-group1/nate/diehl/alpha_test/4"
MOBILE_PDB="${ROOT_DIR}/5ls7_clean.pdb"

# Sanity check
if [[ ! -f "$MOBILE_PDB" ]]; then
  echo "Error: Mobile PDB (5ls7_clean.pdb) not found at $MOBILE_PDB"
  exit 1
fi

# Traverse all subdirectories
for SUBDIR in "$ROOT_DIR"/*/; do
  # Strip trailing slash and get folder name only
  SUBDIR_NAME=$(basename "$SUBDIR")
  TARGET_CIF="$SUBDIR/modified_protein_data/modified_protein_data_model.cif"

  # Check if target .cif exists
  if [[ ! -f "$TARGET_CIF" ]]; then
    echo "Skipping $SUBDIR_NAME: no modified_protein_data_model.cif found"
    continue
  fi

  echo "Processing: $SUBDIR_NAME"

  # Create apo and holo directories in subdirectory
  APO_DIR="${SUBDIR}/apo"
  HOLO_DIR="${SUBDIR}/holo"
  mkdir -p "$APO_DIR" "$HOLO_DIR"

  # Construct output file paths
  APO_OUT="${APO_DIR}/${SUBDIR_NAME}_apo.pdb"
  HOLO_OUT="${HOLO_DIR}/${SUBDIR_NAME}_holo.pdb"
  PML_SCRIPT="${SUBDIR}/align_${SUBDIR_NAME}.pml"

  # Create PyMOL alignment script
  cat > "$PML_SCRIPT" <<EOF
load $TARGET_CIF, target
load $MOBILE_PDB, mobile

align mobile, target

remove mobile and not (resn ACO or resn MG)
save $HOLO_OUT

remove mobile and not (resn MG)
save $APO_OUT

quit
EOF

  # Run PyMOL alignment headlessly
  pymol -cq "$PML_SCRIPT"
done
