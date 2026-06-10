#!/usr/bin/env bash
#
# Incompressible mesh-mapping convergence driver.
# Drives Exec/RegTests/PipeFlow with several mesh-mapping configurations.
#
# Physical domain (identical across all runs):
#   x in [0,     0.04]  (periodic)
#   y in [-0.01, 0.01]  (no-slip walls)
#   z in [-0.01, 0.01]  (no-slip walls)
#
# Reference run (fac = 1,1,1):
#   AMReX domain = physical domain, n_cell = (N_x, N_y, N_z)
#
# Mapped-i run (fac = s in direction i, 1 elsewhere):
#   AMReX domain in direction i is halved (so L_AMReX_i = L_phys_i / s)
#   n_cell unchanged -- physical spacing dx_AMReX * s = original dx
#
# Fixed physical time; same number of cells per direction; same physical
# dx in each direction across configs.
#
# Note: If compile settings other gnu += OMP were used to build the
# executables to be run, pass a full path as BIN=... in the call to this 
# script.  Also note, a number of other settings can be modified, including
# NS, MAX_STEP, STOP_TIME, RESULTS_DIR, etc

set -euo pipefail

# -------- Configuration ---------------------------------------------------
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$HERE/../../.."
PIPEFLOW_DIR="$PROJECT_ROOT/Exec/RegTests/PipeFlow"
# Prefer the OpenMP binary if it was built; fall back to serial.  Can
# also be overridden via the BIN env variable.
if [[ -z "${BIN:-}" ]]; then
  if [[ -x "${PIPEFLOW_DIR}/PeleLMeX3d.gnu.OMP.ex" ]]; then
    BIN="${PIPEFLOW_DIR}/PeleLMeX3d.gnu.OMP.ex"
  else
    BIN="${PIPEFLOW_DIR}/PeleLMeX3d.gnu.ex"
  fi
fi
INP="${PIPEFLOW_DIR}/input.3d-Poiseuille"

# Resolutions to sweep for the convergence study.  Keep small enough
# to run quickly in CI; can be overridden via env.
: "${NS:=32 64}"          # space-separated list (add 128 for a full study)
: "${MAX_STEP:=20}"       # number of PeleLMeX steps per run
: "${STOP_TIME:=0.003}"   # fallback physical stop time (s)
: "${RESULTS_DIR:=${HERE}/results/incompressible}"

mkdir -p "$RESULTS_DIR"

if [[ ! -x "$BIN" ]]; then
  echo "ERROR: $BIN not found.  Build PipeFlow first:"
  echo "  cd $PIPEFLOW_DIR"
  echo "  make TPL && make -j"
  exit 1
fi

# -------- Helpers ---------------------------------------------------------
# Run a single configuration.
#   $1: tag (e.g. 'ref' / 'mapped_x')
#   $2: n_cell (space-separated "Nx Ny Nz")
#   $3: prob_lo "x y z"
#   $4: prob_hi "x y z"
#   $5: mesh_mapping extras (e.g. 'geometry.mesh_mapping=ConstantMap ConstantMap.scaling_factor=2.0 1.0 1.0')
run_case() {
  local tag="$1"
  local ncell="$2"
  local plo="$3"
  local phi="$4"
  local mm_args="$5"
  local outdir="$RESULTS_DIR/${tag}"
  mkdir -p "$outdir"
  cd "$outdir"

  echo "=== [$(date +%T)] $tag  n_cell=${ncell}  mapping='${mm_args}'"
  # -- run ---
  # Redirect all plt/chk output into the case-local dir.
  #shellcheck disable=SC2086
  "$BIN" "$INP" \
      amr.n_cell="$ncell" \
      "geometry.prob_lo=$plo" "geometry.prob_hi=$phi" \
      amr.max_level=0 \
      amr.max_step="$MAX_STEP" \
      amr.stop_time="$STOP_TIME" \
      amr.plot_file=plt_ \
      amr.check_file=chk_ \
      amr.plot_int="$MAX_STEP" \
      amr.check_int=-1 \
      amrex.fpe_trap_invalid=0 \
      amrex.fpe_trap_zero=0 \
      amrex.fpe_trap_overflow=0 \
      $mm_args \
      > run.log 2>&1
  cd "$HERE"
  echo "   ... done"
}

# -------- Per-resolution sweep --------------------------------------------
# Physical domain is 0.04 x 0.02 x 0.02, aspect ratio 2:1:1.  To keep
# AMReX dx isotropic in the reference (unmapped) run we use
# n_cell = 2*N, N, N (proportional to the domain lengths).  In the
# mapped runs the AMReX domain in the stretched direction is halved,
# so dx_AMReX becomes anisotropic -- but the *physical* spacing
# fac_i * dx_AMReX stays isotropic and matches the reference.
for N in $NS; do
  tag_suffix="N${N}"
  N2="$((2 * N))"

  # --- Reference, no mapping ---
  run_case "ref_${tag_suffix}" \
           "$N2 $N $N" \
           "0.0 -0.01 -0.01" \
           "0.04 0.01 0.01" \
           ""

  # --- Identity mapping (sanity) ---
  run_case "ident_${tag_suffix}" \
           "$N2 $N $N" \
           "0.0 -0.01 -0.01" \
           "0.04 0.01 0.01" \
           "geometry.mesh_mapping=ConstantMap ConstantMap.scaling_factor=1.0 1.0 1.0"

  # --- Mapped x  (fac_x = 2) : halve AMReX x ---
  run_case "mapped_x_${tag_suffix}" \
           "$N2 $N $N" \
           "0.0 -0.01 -0.01" \
           "0.02 0.01 0.01" \
           "geometry.mesh_mapping=ConstantMap ConstantMap.scaling_factor=2.0 1.0 1.0"

  # --- Mapped y  (fac_y = 2) : halve AMReX y ---
  run_case "mapped_y_${tag_suffix}" \
           "$N2 $N $N" \
           "0.0 -0.01 -0.01" \
           "0.04 0.0 0.01" \
           "geometry.mesh_mapping=ConstantMap ConstantMap.scaling_factor=1.0 2.0 1.0"

  # --- Mapped z  (fac_z = 2) : halve AMReX z ---
  run_case "mapped_z_${tag_suffix}" \
           "$N2 $N $N" \
           "0.0 -0.01 -0.01" \
           "0.04 0.01 0.0" \
           "geometry.mesh_mapping=ConstantMap ConstantMap.scaling_factor=1.0 1.0 2.0"
done

echo "=== incompressible sweep done.  Results in $RESULTS_DIR"
