#!/usr/bin/env bash
# submit_sweep.sh
# Generate & submit Slurm jobs sweeping Ma, Kn, tmax.

set -euo pipefail

# === Edit these lists ===
# MA_LIST=(0.0)
# KN_LIST=(1.0)
# TMAX_LIST=(0.1)

# MA_LIST=(0.0 3.0)
MA_LIST=(0.0)
KN_LIST=(1.0 0.1 0.01 0.001)
TMAX_LIST=(0.04)

# Sweep mode:
#   cross -> cartesian product of all three lists
#   zip   -> pairwise (i-th of each list together); lists must be equal length
MODE="cross"

# Base Slurm script
BASE="hyqmom_base.sbatch"

# Optional: set DRY_RUN=1 to preview without submitting
DRY_RUN="${DRY_RUN:-0}"

# Sanitize floats so they’re safe in filenames/job names (e.g., 1.0 -> 1p0, 1e-3 -> 1e m3)
sanitize_float() {
  local s="$1"
  s="${s//./p}"   # dot -> p
  s="${s//-/m}"   # minus -> m
  s="${s//+/p}"   # plus -> p (rare in user input, but safe)
  echo "$s"
}

# Create a modified copy and submit
submit_combo() {
  local ma="$1" kn="$2" tmax="$3"
  local J_MA J_KN J_TMAX
  J_MA="$(sanitize_float "$ma")"
  J_KN="$(sanitize_float "$kn")"
  J_TMAX="$(sanitize_float "$tmax")"

  local job="hyqmom_Ma${J_MA}_Kn${J_KN}_t${J_TMAX}"
  local out="${job}.sbatch"

  # Generate the job file by editing the template in-place:
  # - replace the job name
  # - replace the CLI flags --Ma, --Kn, --tmax (numeric patterns)
  awk -v ma="$ma" -v kn="$kn" -v tmax="$tmax" -v job="$job" '
    {
      # Update job name line
      if ($0 ~ /^#SBATCH[[:space:]]+-J[[:space:]]+/) {
        print "#SBATCH -J " job
        next
      }
      # Replace numeric args in the run line(s)
      gsub(/--Ma[[:space:]]+[0-9.eE+-]+/,  "--Ma " ma)
      gsub(/--Kn[[:space:]]+[0-9.eE+-]+/,  "--Kn " kn)
      gsub(/--tmax[[:space:]]+[0-9.eE+-]+/, "--tmax " tmax)
      print
    }
  ' "$BASE" > "$out"

  echo "Prepared ${out}  (job: ${job}; Ma=${ma}, Kn=${kn}, tmax=${tmax})"
  if [[ "$DRY_RUN" != "1" ]]; then
    sbatch "$out"
  fi
}

# --- Main ---
if [[ ! -f "$BASE" ]]; then
  echo "ERROR: Template '$BASE' not found. Save your original batch as ${BASE}." >&2
  exit 1
fi

if [[ "$MODE" == "zip" ]]; then
  if (( ${#MA_LIST[@]} != ${#KN_LIST[@]} || ${#MA_LIST[@]} != ${#TMAX_LIST[@]} )); then
    echo "zip mode requires equal-length MA_LIST, KN_LIST, TMAX_LIST." >&2
    exit 1
  fi
  for i in "${!MA_LIST[@]}"; do
    submit_combo "${MA_LIST[$i]}" "${KN_LIST[$i]}" "${TMAX_LIST[$i]}"
  done
else
  # Cartesian product
  for ma in "${MA_LIST[@]}"; do
    for kn in "${KN_LIST[@]}"; do
      for tmax in "${TMAX_LIST[@]}"; do
        submit_combo "$ma" "$kn" "$tmax"
      done
    done
  done
fi