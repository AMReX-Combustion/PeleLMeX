#!/usr/bin/env bash
#
# Buoyancy-free low-Mach mesh-mapping convergence driver.
#
# Drives Exec/RegTests/HotBubble 2D with gravity OFF and transport
# coefficients turned on.  The hot Gaussian T bump diffuses into the
# surrounding fluid; the resulting thermal expansion drives divU != 0
# and forces the projection to do real work on mapped grids -- but
# without the buoyancy-feedback amplification that dominates the
# gravity-on HotBubble sweep.  Under this setup, mapped vs ref should
# agree at 2nd order as N grows.
#
# Physical domain: 0.016 x 0.032  (same as the standard HotBubble run).
# BC: periodic in x, symmetry lo / outflow hi in y (inherited from
# input.2d-regt).  With gravity = 0 and a Gaussian T bump centred in
# the domain, early-time dynamics are near-symmetric about x = 0.008.
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
: "${MAX_STEP:=40}"       # steps (short, to keep things manifestly in the
                          # thermal-diffusion-dominated regime)
: "${STOP_TIME:=0.05}"
: "${CFL:=0.5}"
: "${VISC:=1.0e-3}"        # Pa.s (~100x air, to speed viscous damping)
: "${COND:=2.5}"           # W/m/K  (~100x air, to compress the thermal
                            #         diffusion timescale of the Gaussian
                            #         T bump down to ~4 ms, well within
                            #         a 60-step convergence run)
: "${RESULTS_DIR:=${HERE}/results/lowmach_nograv}"

mkdir -p "$RESULTS_DIR"

if [[ ! -x "$BIN" ]]; then
  echo "ERROR: $BIN not found.  Build HotBubble first:"
  echo "  cd $HOTBUBBLE_DIR"
  echo "  make TPL && make -j"
  exit 1
fi

# Common input overrides (kill gravity, turn on transport coefficients).
# IMPORTANT: the HotBubble IC places the bubble at y=bubble_y0.  With
# gravity=0 and symmetry-lo at y=0 we want the bubble safely away from
# y=0 so its expansion isn't fighting the symmetry plane.  Default
# bubble_y0 = 0.01 sits fine within y in [0,0.032].
COMMON_OVERRIDES=(
  "peleLM.gravity=0.0 0.0 0.0"
  "transport.const_viscosity=${VISC}"
  "transport.const_bulk_viscosity=0.0"
  "transport.const_conductivity=${COND}"
  "transport.const_diffusivity=0.0"
)

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
      amr.cfl="$CFL" \
      amr.dt_shrink=1.0 \
      amr.init_dt=1.0e-5 \
      amr.max_dt=5.0e-4 \
      amr.dt_change_max=1.3 \
      amrex.fpe_trap_invalid=0 \
      amrex.fpe_trap_zero=0 \
      amrex.fpe_trap_overflow=0 \
      "${COMMON_OVERRIDES[@]}" \
      $mm_args \
      > run.log 2>&1
  cd "$HERE"
  echo "   ... done"
}

# n_cell = N x 2N (preserve aspect ratio).
for N in $NS; do
  Ny=$((2 * N))
  tag_suffix="N${N}"

  # Reference: no mapping
  run_case "ref_${tag_suffix}" \
           "$N $Ny" \
           "0.0 0.0" \
           "0.016 0.032" \
           ""

  # Identity mapping (sanity)
  run_case "ident_${tag_suffix}" \
           "$N $Ny" \
           "0.0 0.0" \
           "0.016 0.032" \
           "geometry.mesh_mapping=ConstantMap ConstantMap.scaling_factor=1.0 1.0"

  # Mapped x (fac_x = 2): halve AMReX x
  run_case "mapped_x_${tag_suffix}" \
           "$N $Ny" \
           "0.0 0.0" \
           "0.008 0.032" \
           "geometry.mesh_mapping=ConstantMap ConstantMap.scaling_factor=2.0 1.0"

  # Mapped y (fac_y = 2): halve AMReX y
  run_case "mapped_y_${tag_suffix}" \
           "$N $Ny" \
           "0.0 0.0" \
           "0.016 0.016" \
           "geometry.mesh_mapping=ConstantMap ConstantMap.scaling_factor=1.0 2.0"
done

echo "=== lowmach_nograv sweep done.  Results in $RESULTS_DIR"
