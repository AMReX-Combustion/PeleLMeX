# LidDrivenCavity

Classical 2D lid-driven cavity benchmark, used here to probe the
accuracy and convergence characteristics of single-level stretched
meshes (`TanhStretchMap`) on a problem with corner singularities and
boundary layers along all four walls.

## Geometry and BCs

Square domain `[0,1] x [0,1]`.  Three stationary walls with no-slip
BCs (`NoSlipWallAdiab`).  Moving lid at `y = 1` with tangential
velocity `u = U_lid = 1`, imposed by setting the high-y face to
`Inflow` and overriding the tangential velocity in `bcnormal`:

```cpp
if (idir == 1 && sgn == -1) {
  s_ext[VELX] = prob_parm.U_lid;
}
```

## Mesh mapping

`TanhStretchMap` with independent `beta_x` and `beta_y` (default
`2.0 2.0`) clusters cells against all four walls.  The local stretch
factor is
```
fac(eta) = beta * sech^2(beta * (2*eta - 1)) / tanh(beta),    eta in [0,1]
```
which is symmetric about the axis midpoint, smallest at the walls,
largest at the centerline.  The ratio of largest to smallest cell
along a stretched axis equals `cosh^2(beta)`:

| `beta` | max / min cell ratio | dy_wall (N=64) | dy_center (N=64) |
|-------:|---------------------:|--------------:|----------------:|
| 0      | 1.0                  | 0.0156        | 0.0156          |
| 1      | 2.38                 | 0.0117        | 0.0202          |
| 2      | 14.15                | 0.00229       | 0.0324          |
| 3      | 101.4                | 0.000459      | 0.0466          |
| 4      | 745.                 | 0.000060      | 0.0593          |

(Same formula and ratios apply along x, independently.)

## Reynolds-number sweep

Three different Reynolds numbers were tested in this folder, all at `N = 64x64`, `beta = 2 2`:

| file | Re  | `peleLM.mu` | `stop_time` |
|---|---:|---:|---:|
| `input.2d-Re100`  | 100  | 0.01   | 100 |
| `input.2d-Re400`  | 400  | 0.0025 | 100 |
| `input.2d-Re1000` | 1000 | 0.001  | 150 |

Each uses `peleLM.incompressible = 1`, `peleLM.rho = 1`, and the BC
combination `lo_bc = NoSlipWallAdiab NoSlipWallAdiab`,
`hi_bc = NoSlipWallAdiab Inflow`.  Re is controlled by changing `peleLM.mu`

## Reference data

Ghia, Ghia & Shin, "High-Re solutions for incompressible flow using
the Navier-Stokes equations and a multigrid method," *Journal of
Computational Physics* 48 (1982) 387-411.  They tabulate `u(x=0.5, y)`
and `v(x, y=0.5)` on the cavity centerlines at Re = 100, 400, 1000,
3200, 5000, 7500, and 10000.

Centerline profiles can be extracted from a plotfile with `yt` using
the provided script `postprocess_cavity.py`

## Convergence study

A useful study sweeps two axes independently:

  - `N` (resolution): pair `N = 32, 64, 128` at fixed `beta = 2`.
    Compares solution at each resolution to either Ghia tabulated
    values (preferred) or to a Richardson-extrapolated reference.
  - `beta` (clustering): pair `beta = 0, 1, 2, 3` at fixed `N = 64`.
    Quantifies the accuracy improvement at the wall and the loss of
    resolution in the cavity core as cells migrate outward.

At high `beta` the MAC linear-solver convergence degrades (see the
`MeshMappingConvergence/run_stretch_sweep.sh` analysis in the
`PipeFlow` regtest).  Beyond about `beta = 3-4` at N=64 you will need
to switch to `mac_proj.use_mlhypre = 1` with HYPRE / BoomerAMG and
`nodal_proj.bottom_solver = hypre` -- requires building PeleLMeX with
`USE_HYPRE = TRUE` (see the AMReX documentation for building an
appropriate version of hypre and setting `HYPER_DIR`)

## Mesh-mapping-aware initial data and bcnormal (opt-in API)

The driver passes a device-callable `MeshMapEvaluator mmap` to
`initdata_mapped` and uses it to compute the `x[]` array that the
driver hands to `bcnormal`.  Both behaviours are zero-cost on
unmapped runs (the evaluator falls back to `Kind::Identity` and the
helpers return the unmapped xi-space coordinate).

### `initdata` (legacy) vs `initdata_mapped` (opt-in)

The cavity here uses the legacy `initdata` (its IC is `u = 0`
everywhere, so no coordinate is needed).  For a problem that needs
the physical cell-centre coordinate, define `initdata_mapped`
instead — the driver detects this at compile time via the
`has_initdata_mapped_v<>` trait and dispatches accordingly.  Sketch:

```cpp
static void initdata_mapped(
  int i, int j, int k, int is_incomp,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& aux,
  amrex::GeometryData const& geomdata,
  MeshMapEvaluator const& mmap,        // <-- new arg
  MyProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const* pmf_data)
{
  // Physical cell-centre coordinates under whatever map is active
  // (Identity when none is active, so this is safe in both cases).
  const amrex::Real x = mmap.x_phys_cc(0, i, geomdata);
  const amrex::Real y = mmap.x_phys_cc(1, j, geomdata);
  // ... use (x, y) to set state(i, j, k, ...) ...
}
```

Defining both `initdata` and `initdata_mapped` is legal but
`initdata_mapped` wins; the legacy `initdata` becomes dead code in
that case.

### `bcnormal` (unchanged signature, smarter `x`)

The `x[]` array that the driver hands to `bcnormal` is now the
**physical face-centre position** for `x[idir]` (the boundary
direction being filled) and the physical cell-centre for the
tangential axes.  Under no mesh mapping these reduce to the standard
uniform-grid formulas, with one small correction: `x[idir]` is now
at the boundary face rather than the ghost-cell centre (a half-cell
shift in the boundary direction), matching the documented intent of
the `bcnormal` interface.  User code that branches on `idir` /
`sgn` (as this cavity does) is bit-exact unaffected; code that
queries `x[idir]` for the boundary location now sees the actual
wall position.

## Known limitations

  - **Single-level only.**  Multi-level AMR with mesh mapping is an
    inherited AmrWind limitation; the input files pin
    `amr.max_level = 0`.
  - **Time-step floor.**  At `beta = 2`, the smallest near-wall cell
    is ~7x smaller than the uniform-grid equivalent, so CFL-limited
    `dt` is ~7x smaller.  Each Re case takes several thousand steps
    to reach the Ghia-comparable steady state.
