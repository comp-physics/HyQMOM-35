#!/usr/bin/env bash
# submit_sweep.sh
# Generate & submit Slurm jobs sweeping Ma, Kn, tmax, and grid resolution.
#
# Usage (run from HyQMOM.jl directory):
#   cd HyQMOM.jl
#   bash slurm/sweep.sh

set -euo pipefail

# Check we're in the correct directory
if [ ! -f "Project.toml" ]; then
    echo "ERROR: Must run from HyQMOM.jl directory (Project.toml not found)" >&2
    echo "Usage:" >&2
    echo "  cd HyQMOM.jl" >&2
    echo "  bash slurm/sweep.sh" >&2
    exit 1
fi

# === Edit these lists ===
# Grid resolution (cubic grids: N used for Nx, Ny, Nz)
N_LIST=(80)

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
  local n="$1" ma="$2" kn="$3" tmax="$4"
  local J_MA J_KN J_TMAX
  J_MA="$(sanitize_float "$ma")"
  J_KN="$(sanitize_float "$kn")"
  J_TMAX="$(sanitize_float "$tmax")"

  # Job name matches output file naming: N{n}_Ma{ma}_Kn{kn}_t{tmax}
  # (Output files use Nx/Ny/Nz but all equal to N for cubic grids)
  local job="hyqmom_N${n}_Ma${J_MA}_Kn${J_KN}_t${J_TMAX}"
  local out="slurm/${job}.sbatch"

  # Generate the job file by editing the template:
  # - replace the job name
  # - replace the CLI flags --Nx, --Ny, --Nz (all set to n), --Ma, --Kn, --tmax
  awk -v n="$n" -v ma="$ma" -v kn="$kn" -v tmax="$tmax" -v job="$job" '
    {
      # Update job name line
      if ($0 ~ /^#SBATCH[[:space:]]+-J[[:space:]]+/) {
        print "#SBATCH -J " job
        next
      }
      # Replace numeric args in the run line(s)
      gsub(/--Nx[[:space:]]+[0-9]+/,       "--Nx " n)
      gsub(/--Ny[[:space:]]+[0-9]+/,       "--Ny " n)
      gsub(/--Nz[[:space:]]+[0-9]+/,       "--Nz " n)
      gsub(/--Ma[[:space:]]+[0-9.eE+-]+/,  "--Ma " ma)
      gsub(/--Kn[[:space:]]+[0-9.eE+-]+/,  "--Kn " kn)
      gsub(/--tmax[[:space:]]+[0-9.eE+-]+/, "--tmax " tmax)
      print
    }
  ' "$BASE" > "$out"

  echo "Prepared ${out}"
  echo "  Grid: ${n}³, Ma=${ma}, Kn=${kn}, tmax=${tmax}"
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
  if (( ${#N_LIST[@]} != ${#MA_LIST[@]} || ${#N_LIST[@]} != ${#KN_LIST[@]} || 
        ${#N_LIST[@]} != ${#TMAX_LIST[@]} )); then
    echo "zip mode requires equal-length N_LIST, MA_LIST, KN_LIST, TMAX_LIST." >&2
    exit 1
  fi
  for i in "${!N_LIST[@]}"; do
    submit_combo "${N_LIST[$i]}" "${MA_LIST[$i]}" "${KN_LIST[$i]}" "${TMAX_LIST[$i]}"
  done
else
  # Cartesian product
  for n in "${N_LIST[@]}"; do
    for ma in "${MA_LIST[@]}"; do
      for kn in "${KN_LIST[@]}"; do
        for tmax in "${TMAX_LIST[@]}"; do
          submit_combo "$n" "$ma" "$kn" "$tmax"
        done
      done
    done
  done
fi