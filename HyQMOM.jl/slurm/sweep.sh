#!/usr/bin/env bash
# submit_sweep.sh
# Generate & submit Slurm jobs sweeping Ma, Kn, tmax, and grid resolution.
#
# Usage (run from HyQMOM.jl directory):
#   cd HyQMOM.jl
#   bash slurm/sweep.sh

set -euo pipefail

# === Edit these lists ===
# Grid resolution (will be used as Nx, Ny, Nz - cubic grids)
NX_LIST=(80)
NY_LIST=(80)
NZ_LIST=(80)

# Physical parameters
MA_LIST=(0.0)
KN_LIST=(1.0 0.1 0.01 0.001)
TMAX_LIST=(0.04)

# Sweep mode:
#   cross -> cartesian product of all lists
#   zip   -> pairwise (i-th of each list together); lists must be equal length
MODE="cross"

# Base Slurm script (template)
BASE="slurm/hyqmom_base.sbatch"

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
  local nx="$1" ny="$2" nz="$3" ma="$4" kn="$5" tmax="$6"
  local J_MA J_KN J_TMAX
  J_MA="$(sanitize_float "$ma")"
  J_KN="$(sanitize_float "$kn")"
  J_TMAX="$(sanitize_float "$tmax")"

  # Job name matches output file naming: Nx{nx}_Ny{ny}_Nz{nz}_Ma{ma}_Kn{kn}_t{tmax}
  local job="hyqmom_Nx${nx}_Ny${ny}_Nz${nz}_Ma${J_MA}_Kn${J_KN}_t${J_TMAX}"
  local out="slurm/${job}.sbatch"

  # Generate the job file by editing the template:
  # - replace the job name
  # - replace the CLI flags --Nx, --Ny, --Nz, --Ma, --Kn, --tmax
  awk -v nx="$nx" -v ny="$ny" -v nz="$nz" -v ma="$ma" -v kn="$kn" -v tmax="$tmax" -v job="$job" '
    {
      # Update job name line
      if ($0 ~ /^#SBATCH[[:space:]]+-J[[:space:]]+/) {
        print "#SBATCH -J " job
        next
      }
      # Replace numeric args in the run line(s)
      gsub(/--Nx[[:space:]]+[0-9]+/,       "--Nx " nx)
      gsub(/--Ny[[:space:]]+[0-9]+/,       "--Ny " ny)
      gsub(/--Nz[[:space:]]+[0-9]+/,       "--Nz " nz)
      gsub(/--Ma[[:space:]]+[0-9.eE+-]+/,  "--Ma " ma)
      gsub(/--Kn[[:space:]]+[0-9.eE+-]+/,  "--Kn " kn)
      gsub(/--tmax[[:space:]]+[0-9.eE+-]+/, "--tmax " tmax)
      print
    }
  ' "$BASE" > "$out"

  echo "Prepared ${out}"
  echo "  Grid: ${nx}×${ny}×${nz}, Ma=${ma}, Kn=${kn}, tmax=${tmax}"
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
  if (( ${#NX_LIST[@]} != ${#NY_LIST[@]} || ${#NX_LIST[@]} != ${#NZ_LIST[@]} || 
        ${#NX_LIST[@]} != ${#MA_LIST[@]} || ${#NX_LIST[@]} != ${#KN_LIST[@]} || 
        ${#NX_LIST[@]} != ${#TMAX_LIST[@]} )); then
    echo "zip mode requires equal-length NX_LIST, NY_LIST, NZ_LIST, MA_LIST, KN_LIST, TMAX_LIST." >&2
    exit 1
  fi
  for i in "${!MA_LIST[@]}"; do
    submit_combo "${NX_LIST[$i]}" "${NY_LIST[$i]}" "${NZ_LIST[$i]}" \
                 "${MA_LIST[$i]}" "${KN_LIST[$i]}" "${TMAX_LIST[$i]}"
  done
else
  # Cartesian product
  for nx in "${NX_LIST[@]}"; do
    for ny in "${NY_LIST[@]}"; do
      for nz in "${NZ_LIST[@]}"; do
        for ma in "${MA_LIST[@]}"; do
          for kn in "${KN_LIST[@]}"; do
            for tmax in "${TMAX_LIST[@]}"; do
              submit_combo "$nx" "$ny" "$nz" "$ma" "$kn" "$tmax"
            done
          done
        done
      done
    done
  done
fi