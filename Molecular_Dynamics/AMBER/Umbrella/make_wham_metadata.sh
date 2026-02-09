#!/usr/bin/env bash
set -euo pipefail

ROOT="."
OUT="metadata.wham"
OVERLAP="overlap_samples.tsv"

# Your generator: rk2 is kcal/mol/rad^2
RK2_KCAL_PER_RAD2=100.0

# Convert AMBER rk2 (kcal/mol/rad^2) -> WHAM k (kcal/mol/deg^2),
# including the factor-of-2 convention (kwham = 2*rk2*(pi/180)^2).
# Pure awk (no python): pi via atan2(0,-1)
KWHAM=$(awk -v rk2="$RK2_KCAL_PER_RAD2" 'BEGIN{
  pi = atan2(0,-1);
  kwham = 2.0 * rk2 * (pi/180.0)^2;
  printf "%.12g\n", kwham
}')

: > "$OUT"
: > "$OVERLAP"

for rst in "$ROOT"/chi_*/disang_prod*.RST; do
  wdir=$(dirname "$rst")                   # chi_090
  tag=${wdir##*/}                          # chi_090
  X=${tag#chi_}                            # 090

  dat="$wdir/chi_${X}.dat"                 # chi_090/chi_090.dat
  [[ -f "$dat" ]] || { echo "Missing $dat"; exit 1; }

  # Center from r2 (degrees)
  center=$(awk '
    BEGIN{IGNORECASE=1}
    match($0,/r2[[:space:]]*=[[:space:]]*([-+0-9.]+)/,m){print m[1]; exit}
  ' "$rst")
  [[ -n "${center:-}" ]] || { echo "Could not parse r2 from $rst"; exit 1; }

  # 1) WHAM metadata
  printf "%s\t%s\t%s\n" "$dat" "$center" "$KWHAM" >> "$OUT"

  # 2) Long-form samples for overlap checking
  #    - skips blank lines and comment lines starting with # or @
  #    - uses the *last column* as the chi value (works for "time chi" or "chi" formats)
  #    - wraps chi into [0,360) so 370 -> 10 and -5 -> 355 (important for circular dihedrals)
  awk -v win="$tag" -v ctr="$center" '
    function wrap360(x){
      x = x % 360.0
      if (x < 0) x += 360.0
      return x
    }
    NF==0 {next}
    $1 ~ /^[#@]/ {next}
    {
      v = $(NF) + 0.0
      v = wrap360(v)
      printf "%s\t%.6f\t%.6f\n", win, ctr, v
    }
  ' "$dat" >> "$OVERLAP"

done

echo "Wrote $OUT"
echo "Wrote $OVERLAP (window  center_deg  chi_deg)"
head -n 10 "$OUT"
echo "----"
head -n 10 "$OVERLAP"
