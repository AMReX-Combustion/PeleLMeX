# Mesh-mapping convergence tests

This harness drives PeleLMeX with several mesh-mapping configurations
and reports velocity-field norms so we can check that mapping produces
numerically consistent results.

Three sweeps are available:

| driver                       | physical case | regime       | primary purpose |
|------------------------------|---------------|--------------|-----------------|
| `run_incompressible.sh`      | `PipeFlow`    | incompressible, inviscid | clean convergence check on the mechanical scaling (MAC proj, nodal proj, advection, CFL) |
| `run_lowmach.sh`             | `HotBubble`   | low-Mach, inert, gravity ON | end-to-end low-Mach path, *including* buoyancy-feedback amplification |
| `run_lowmach_nograv.sh`      | `HotBubble`   | low-Mach, inert, gravity OFF, conductivity + viscosity ON | low-Mach path WITHOUT buoyancy amplification: deviation plateaus instead of growing |
| `run_stretch_sweep.sh`       | `PipeFlow`    | incompressible, inviscid | Sweep through mesh stretching factors (requires executable built with hypre) |

## Protocol (all three sweeps)

For each of the first three cases, five runs are executed at each resolution `N`:

1. **ref** — `fac = (1,1,…)` on an AMReX grid matching the physical domain.
2. **ident** — `fac = (1,1,…)` with mesh mapping *enabled* (sanity / bit-identity check).
3. **mapped_x** — `fac_x = 2`, AMReX x-extent halved, physical extent preserved.
4. **mapped_y** — `fac_y = 2`, AMReX y-extent halved, physical extent preserved.
5. **mapped_z** — `fac_z = 2`, only in 3D.

Every run has the same physical domain, the same `n_cell` per direction,
and therefore the same physical cell spacing.  Cells at index `(i,j,k)`
occupy the same physical location in every configuration, so
cell-by-cell comparison is meaningful.

Initial conditions in both `PipeFlow/pelelmex_prob.H` and
`HotBubble/pelelmex_prob.{H,cpp}` are *mapping-aware*: they evaluate the
IC from physical coordinates (`x_phys = prob_lo + (i + 0.5) * dx *
fac`), so the starting state agrees bit-for-bit across ref and mapped
runs when the AMReX grid spans the same physical region.

Multi-level AMR + mesh mapping is an inherited AmrWind limitation and
is not exercised here; the driver scripts pin `amr.max_level = 0`.

For the stretch sweep, 32 and 64 cell cases are tried with the mesh stretching
beta = 0,1,2,3,4,5,6.  If AMReX GMG is used (when built with USE_HYPRE=FALSE)
most of these cases will fail. The set of knobs currently included here
will enable running with beta <=~ 4.

## Physical-space plotfile rendering

When mesh mapping is active, PeleLMeX's `WritePlotFile` emits a per-level
nodal `MultiFab` of node-displacement (`x_phys − x_xi`) alongside the
standard cell-centered data and appends the AMReX ParaView/VisIt
plugin's `amrexvec` trailer/block to the plotfile `Header`, with the
vector components written as `nu_x`, `nu_y`, and `nu_z`.  A
ParaView/VisIt session loading these plotfiles via the AMReX reader
automatically renders the solution on the curvilinear physical grid —
no user-side state file, calculator, or warp filter required.  This is
the same on-disk protocol used by ERF for its terrain-following output
(`WriteGenericPlotfileHeaderWithTerrain`).

Kill-switch: set `peleLM.plot_mesh_mapping = 0` to fall back to the
standard `WriteMultiLevelPlotfile` call, in which case the plotfile
is byte-compatible with the pre-mapping format and renders in
&xi;-space.

`yt` and other tools that only look at the standard cell-centered
variables (including `analyze.py` in this directory) are unaffected by
the extra metadata; they continue to read the plotfile in its native
&xi;-space coordinates.

## Running

```bash
cd Exec/RegTests/MeshMappingConvergence

# Build dependencies first:
#   cd ../PipeFlow && make TPL && make -j
#   cd ../HotBubble && make TPL && make -j

./run_incompressible.sh       # PipeFlow sweep
./run_lowmach.sh              # HotBubble sweep (buoyancy ON)
./run_lowmach_nograv.sh       # HotBubble sweep (buoyancy OFF)
./run_stretch_sweep.sh        # PipeFlow sweep

python3 analyze.py results/   # reads plotfiles, prints tables
```

Requires a Python environment with `yt` and `numpy`.

The `NS` environment variable overrides the resolution sweep
(default `"32 64"`); use `NS="32 64 128"` for a full study.  Likewise
`MAX_STEP`, `STOP_TIME`, `CFL`, `VISC`, `COND`, and `RESULTS_DIR` can
be overridden per-invocation.  See the stretch sweep file for additional
controls.

## Expected results

Indicative numbers at `N = 32` (from `analyze.py`, as of the last run
of this harness):

| sweep              | `ident / ref` | `mapped_x / ref` L2 | `mapped_y / ref` L2 | notes |
|--------------------|---------------|---------------------|---------------------|-------|
| incompressible     | bit-identical | 0.99985             | 0.99952             | excellent agreement, no amplification loop |
| lowmach (g ON)     | bit-identical | 0.974               | 0.975               | 2–3 % deviation, grows ~×2/step from buoyancy feedback |
| lowmach_nograv     | bit-identical | 1.048               | 1.048               | ~5 % deviation, **plateaus**; no amplification |

For the stretch sweep, even with hypre installed, cases with beta>4 will
fail in one of the projection steps.  Leading up to this failure in beta
one should observe dramatically increasing numbers of solver iterations,
particularly for the mac.  Larger runs at the smaller beta values will be
required to verify convergence order.

## Why the three sweeps give different magnitudes

Mathematically, ref and mapped are two consistent 2nd-order
discretisations of the same PDE on the same physical domain.  Their
per-step truncation errors agree to *higher* than leading order (empirically
O(dt³) per step at step 1), so the two solutions stay close on their
own.

When a physical mechanism *amplifies* small differences, those
higher-order truncation differences get magnified into visible
end-state deviations.  Here the mechanism is buoyancy feedback:

  δρ  →  δg · δρ  →  δu  →  advection of δρ  →  …

Around a hot low-density bubble in gravity, this is roughly a ×2
per-step multiplication at CFL 0.9.  Twenty steps at ×2 is a ×10⁶
amplification of an initial ~1e-8 truncation difference into the ~1e-2
end-state deviation we see.

Evidence that the mapping itself is numerically correct (not buggy):

- **`ident` == `ref` bit-identical** in all three sweeps.
- **Gravity off, no transport** (uniform density, at-rest) → mapped ==
  ref bit-identical (the mapping arithmetic is exact in that limit).
- **Fixed very small dt** → the step-1 mapped-vs-ref deviation drops to
  float round-off (~1e-16), confirming no dt-independent coding error.
- **Per-step error scales ~O(dt⁴) at step 1** — much steeper than
  2nd-order truncation; consistent with both schemes sharing their
  leading-order stencil and differing only in higher-order terms.
- **Incompressible sweep (no buoyancy loop)** agrees to ~0.05 %.
- **`lowmach_nograv` sweep (low-Mach, no buoyancy loop)** agrees to
  ~5 % *and plateaus*, vs. the gravity-on `lowmach` sweep that
  amplifies toward 2.5 % over 20 steps from a smaller per-step
  starting point.

The `lowmach_nograv` sweep is the cleanest low-Mach convergence
check: it exercises the divU-aware projection and scalar-advection
paths through thermal expansion of a diffusing hot bubble, without
being dominated by the buoyancy amplification of the standard
HotBubble case.

## Historical caveat (now resolved)

Earlier versions of this harness noted a suspected σ_x-only bug in
`MLNodeLaplacian::updateVelocity`/`getFluxes` in AMReX, which was
hypothesised to cause the 2–3 % low-Mach deviation.  The fix landed in
AMReX (`mlndlap_mknewu_ha`, feature macro
`AMREX_MLNODELAP_HAS_MKNEWU_HA`) and the corresponding PeleLMeX
post-hoc σ_x correction (in `velocityProjection` /
`initialProjection`) is now gated on that macro.  The fix made the
code cleaner but did **not** change the observed deviation — confirming
that the 2–3 % is buoyancy-amplification of higher-order truncation,
not a coding bug.

## Output structure

```
results/
  incompressible/
    ref_N32/        ident_N32/    mapped_{x,y,z}_N32/
    ref_N64/ ...
  lowmach/
    (same layout, HotBubble with gravity ON)
  lowmach_nograv/
    (same layout, HotBubble with gravity OFF + k, μ ON)
  stretch/
    stretch_N32_b0/
    stretch_N32_b1/
	...
    stretch_N64_b0/
    stretch_N64_b1/
	...
```

Each leaf directory contains PeleLMeX plotfiles and the run log.
