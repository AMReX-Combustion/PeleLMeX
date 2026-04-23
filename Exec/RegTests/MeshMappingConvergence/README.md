# Mesh-mapping convergence tests

This harness drives the existing `PipeFlow` (incompressible) and
`HotBubble` (low-Mach) regression cases with several mesh-mapping
configurations and reports velocity-field norms so we can check that
mapping produces numerically consistent results.

## Protocol

For each physical test case, three comparisons are run:

1. **Equivalence at matched physical resolution.**  The reference run
   uses `fac = (1,1,1)` on an AMReX grid that spans the physical
   domain `L_x x L_y x L_z`.  The three *mapped* runs use
   `fac = (2,1,1)`, `(1,2,1)`, and `(1,1,2)` respectively, with the
   AMReX domain in the stretched direction halved so
   `L_AMReX_i = L_phys_i / fac_i`.  The physical cell spacing
   `dx_AMReX * fac_i` is the same in every direction, so cells at
   the same `(i,j,k)` index occupy the same physical location in all
   four runs.

2. **Convergence rate.**  The same four configurations are re-run at
   `N = 32`, `64` and `128` cells-per-direction.  Global norms should
   decrease at second order under refinement.

3. **Identity-mapping check** (sanity).  `fac = (1,1,1)` with mesh
   mapping turned on must match the no-mapping reference to
   discretization order (often bit-for-bit, depending on which MG
   kernel path the nodal projector takes internally).

### Important caveat: initial conditions

The `PipeFlow` and `HotBubble` `pelelmex_prob.H` files evaluate the
initial condition from **AMReX-grid coordinates**, not physical
coordinates.  Under mesh mapping, the AMReX-grid extent is shrunk by
`1/fac_i` in the stretched direction, so the IC has a different
spatial structure in the mapped runs relative to the reference run.
Consequences:

- Early in the run, mapped and reference solutions differ because
  their ICs differ.
- Once the initial transient decays (viscous dissipation,
  projection, etc.) and the flow reaches a steady (or
  statistically-stationary) state, the QoIs should be comparable.
- We therefore run to a fixed physical final time `T_final` that
  is well past the transient decay time for the chosen physical
  case, and compare at that end state.

For a more rigorous test, the problem's IC code can be made
mapping-aware (evaluate from physical coords instead of AMReX coords)
-- a future refinement; not required for a first pass.

## Running

```
cd Exec/RegTests/MeshMappingConvergence
./run.sh               # runs both incompressible + low-Mach
./run_incompressible.sh  # just PipeFlow
./run_lowmach.sh         # just HotBubble
python3 analyze.py results/   # reads all plotfiles, prints table
```

Requires a python environment with `yt` and `numpy` installed:

```
pip install --user yt numpy
```

## Output structure

```
results/
  incompressible/
    ref_N32/          # fac=(1,1,1), N=32
    mapped_x_N32/     # fac=(2,1,1), N=32
    mapped_y_N32/     # fac=(1,2,1), N=32
    mapped_z_N32/     # fac=(1,1,2), N=32
    ref_N64/  ...
  lowmach/
    ... (same layout, HotBubble)
```

Each leaf directory contains PeleLMeX plotfiles + log.

## What to look for

The `analyze.py` table shows, per run:

- `L2(|u|)`: volume-averaged L2 norm of velocity magnitude
- `Linf(|u|)`: maximum velocity magnitude in the domain
- Ratio vs reference at same N

Expectation:

- Identity mapping matches reference closely (ideally bit-identical).
- Non-identity mapping differs from reference because the ICs differ;
  magnitude of difference should be bounded and decrease with N.
- At fixed `fac`, the norms should converge as N doubles.

## Expected outcomes

If the mesh-mapping implementation is consistent, the following should
hold at the chosen resolutions:

- `ref_Nx / mapped_*_Nx` ratios approach 1 as Nx grows.
- Log-log slope of `|QoI(N) - QoI(2N)|` vs `dx` has slope approx 2.

A failure mode (bug in scaling) would show either (a) blowup or
NaN under non-identity mapping, or (b) failure to converge / wrong
convergence rate.
