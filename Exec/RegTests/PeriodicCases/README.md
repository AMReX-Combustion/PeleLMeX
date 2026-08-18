## PeriodicCases
A set of 2D/3D periodic cases: convected vortex in different directions, convected temperature/mixture Gaussian bump, and
pure diffusion of species/temperature. This is the base case to test PeleLMeX accuracy and order of convergence of
the spatial and temporal schemes.

### ShearLayer: a temporal mixing layer for interior mesh stretching

`prob.type = ShearLayer` adds a temporally-developing planar mixing
layer, and exists mainly as the worked example for
`InteriorStretchMap` (see `Source/Mesh/PeleLMeX_InteriorStretchMap.H`).
Despite living in `PeriodicCases`, it is periodic in **x only** — `y` is
closed with `SlipWallAdiab` far from the layer.

The initial condition is

```
u(y) = Uc + (dU/2) tanh(2 (y - y0) / dw)
```

with `dw` the vorticity thickness, so `du/dy` peaks at `dU/dw` on the
layer, plus a Kelvin–Helmholtz seed taken as the curl of a
streamfunction

```
psi' = A sin(kx x) exp(-((y-y0)/dw)^2),   kx = 2 pi nwave / Lx
```

so that `u' = -d psi'/dy`, `v' = +d psi'/dx` is analytically
divergence-free and `A` is scaled to give `max|v'| = shear_pert * dU`.

| input key | meaning | default |
|---|---|---|
| `prob.shear_y0` | physical y of the layer | 0.0 |
| `prob.shear_dw` | vorticity thickness | 0.004 |
| `prob.shear_dU` | velocity difference across the layer | 1.0 |
| `prob.shear_Uc` | mean (convective) velocity | 0.0 |
| `prob.shear_pert` | `max\|v'\| / dU` | 0.02 |
| `prob.shear_nwave` | KH wavelengths across the periodic x extent | 1 |
| `prob.shear_ampY` | optional scalar step amplitude (low-Mach only) | 0.0 |

#### Why the stretched axis must be the non-periodic one

`InteriorStretchMap` clusters at an interior point and coarsens towards
both ends of the axis, so `fac` is *largest* at both ends:
`fac(0) = cosh^2(beta_lo)/N` and `fac(1) = cosh^2(beta_hi)/N`.  Across a
periodic wrap those two are the same location, so the metric jumps there
unless `beta_lo == beta_hi`, and even then `fac` has a C1 kink at the
wrap — exactly the defect the `cosh^2` construction is designed to avoid
at the interior join.  Hence `x` is left uniform (`beta = 0`) and only
`y` is stretched.  A fully periodic double shear layer is doubly
unsuitable: periodicity forces an even number of layers, and one
interior cluster can only serve one of them.

#### What the stretching buys

Two inputs are provided, both at `64 x 128` (`dy_uniform = 0.625 mm`),
`dw = 2.8 mm`, layer at `y = 0` in `y in [-0.03, 0.05]`:

| input | `beta_lo, beta_hi` (y) | dy at layer | cells across `dw` | gain | `cosh^2` |
|---|---|---|---|---|---|
| — (uniform) | 0, 0 | 0.625 mm | 4.5 | 1.00 | 1.00 |
| `input.2d_ShearLayer` | 0.6, 1.0 | 0.485 mm | 6 | 1.29 | 2.38 |
| `input.2d_ShearLayer_beta2` | 1.6, 2.4 | 0.154 mm | 18 | 4.05 | 30.88 |

The resolution gain at the fine point is exactly the normalisation
`N = G(beta)` — `1.41x` at `beta = 1`, `3.91x` at `beta = 2`, `17x` at
`beta = 3`.  Staying inside the default MLMG comfort zone
(`cosh^2 ~ 2.4`) therefore caps the gain near `1.3x`.

`input.2d_ShearLayer_beta2` is the one that shows the map off, and it
**requires** an executable built with `USE_HYPRE = TRUE`: on the default
AMReX MLMG path it fails immediately.  With `mac_proj.use_mlhypre = 1`
and `nodal_proj.bottom_solver = hypre` it runs and the rollup is clean.

That measurement is worth keeping, because it contradicts the
explanation in `LidDrivenCavity/input.2d`, which blames that case's
`beta = 2` MLMG stall on the cavity's corner-singularity divergence
coupling through stretched cells.  This case is periodic in x with slip
walls in y and has no corner singularity, and MLMG still fails — so
operator anisotropy on its own is enough.  Treat a hypre bottom solver
as a requirement for any strongly stretched mesh, independent of
geometry.  Both inputs here have been run: `cosh^2 = 2.38` works on the
default MLMG path and `cosh^2 = 30.9` fails immediately on it, so the
threshold is bracketed between them but not located.

Setting both betas to `0 0` gives an exact identity map (`detJ = 1`,
`fac = 1`) and must reproduce the run with `geometry.mesh_mapping`
omitted — a cheap regression check on the mapping machinery.

#### Mapping-aware initial data

This case defines both `initdata` and `initdata_mapped`; the driver
prefers the latter.  Both delegate to a shared `set_state()`, so they
cannot disagree — the only difference is that the mapped path takes the
cell-centre coordinates from `MeshMapEvaluator::x_phys_cc()` and the
physical x-extent from `domain_size_phys()`.  That matters here: on a
stretched mesh the uniform-grid formula would place the shear layer at
the wrong `y`.  On an unmapped run the evaluator falls back to
`Kind::Identity` and the two paths agree bit-for-bit, so the existing
problem types in this directory are unaffected.
