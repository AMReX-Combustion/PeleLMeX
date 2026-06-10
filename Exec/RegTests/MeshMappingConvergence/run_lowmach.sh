#!/usr/bin/env bash
#
# Low-Mach mesh-mapping convergence driver.
# Drives Exec/RegTests/HotBubble 2D with several mesh-mapping configurations.
#
# The HotBubble case uses a 2D Cartesian domain.  Physical domain:
#   x in [0,   0.016]  (periodic)
#   y in [0,   0.032]  (symmetry lo, outflow hi)
#
# Reference run (fac = 1,1):   AMReX = physical, n_cell = (N_x, N_y)
# Mapped-i run (fac_i = 2):     AMReX domain in dir i halved; n_cell unchanged
#
# The same IC-vs-AMReX-coordinate caveat applies as in the incompressible
# runner -- results are directly comparable only after the initial
# transient has dissipated.
#
# Note: If compile settings other gnu += OMP were used to build the
# executables to be run, pass a full path as BIN=... in the call to this 
# script.  Also note, a number of other settings can be modified, including
# NS, MAX_STEP, STOP_TIME, RESULTS_DIR, etc


set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$HERE/../../.."
HOTBUBBLE_DIR="$PROJECT_ROOT/Exec/RegTests/HotBubble"
if [[ -z "${BIN:-}" ]]; then
  if [[ -x "${HOTBUBBLE_DIR}/PeleLMeX2d.gnu.OMP.ex" ]]; then
    BIN="${HOTBUBBLE_DIR}/PeleLMeX2d.gnu.OMP.ex"
  else
    BIN="${HOTBUBBLE_DIR}/PeleLMeX2d.gnu.ex"
  fi
fi
INP="${HOTBUBBLE_DIR}/input.2d-regt"

: "${NS:=32 64}"          # cells per direction sweep
: "${MAX_STEP:=20}"       # number of time steps
: "${STOP_TIME:=0.05}"    # physical stop time (s)
: "${RESULTS_DIR:=${HERE}/results/lowmach}"

mkdir -p "$RESULTS_DIR"

if [[ ! -x "$BIN" ]]; then
  echo "ERROR: $BIN not found.  Build HotBubble first:"
  echo "  cd $HOTBUBBLE_DIR"
  echo "  make TPL && make -j"
  exit 1
fi

run_case() {
  local tag="$1"
  local ncell="$2"     # "Nx Ny"
  local plo="$3"       # "xlo ylo"
  local phi="$4"       # "xhi yhi"
  local mm_args="$5"
  local outdir="$RESULTS_DIR/${tag}"
  mkdir -p "$outdir"
  cd "$outdir"

  echo "=== [$(date +%T)] $tag  n_cell=${ncell}  mapping='${mm_args}'"
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

# Physical domain: 0.016 x 0.032.  n_cell = N x 2N (preserve aspect ratio).
for N in $NS; do
  Ny=$((2 * N))
  tag_suffix="N${N}"

  # --- Reference, no mapping ---
  run_case "ref_${tag_suffix}" \
           "$N $Ny" \
           "0.0 0.0" \
           "0.016 0.032" \
           ""

  # --- Identity mapping (sanity) ---
  run_case "ident_${tag_suffix}" \
           "$N $Ny" \
           "0.0 0.0" \
           "0.016 0.032" \
           "geometry.mesh_mapping=ConstantMap ConstantMap.scaling_factor=1.0 1.0"

  # --- Mapped x (fac_x = 2) : halve AMReX x ---
  run_case "mapped_x_${tag_suffix}" \
           "$N $Ny" \
           "0.0 0.0" \
           "0.008 0.032" \
           "geometry.mesh_mapping=ConstantMap ConstantMap.scaling_factor=2.0 1.0"

  # --- Mapped y (fac_y = 2) : halve AMReX y ---
  run_case "mapped_y_${tag_suffix}" \
           "$N $Ny" \
           "0.0 0.0" \
           "0.016 0.016" \
           "geometry.mesh_mapping=ConstantMap ConstantMap.scaling_factor=1.0 2.0"
done

echo "=== low-Mach sweep done.  Results in $RESULTS_DIR"
