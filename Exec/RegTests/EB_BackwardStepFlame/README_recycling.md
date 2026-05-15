# Recycling-plane inflow smoke tests

This folder contains six short input files for the recycling-plane inflow
capability (see `Docs/sphinx/manual/LMeXControls.rst`
§ "Recycling-plane inflow"): one baseline/no-op case
(`eb_bfs_recycle_off.inp`) plus five variants that each touch a specific
code path of the new feature. Each file is a 10-step, 128x32, 4-rank
variant of the base `eb_bfs.inp` case. They are registered as CTest
entries in `Tests/CMakeLists.txt`.

## Tests

| Input file | What it exercises | Notes |
|---|---|---|
 | `eb_bfs_recycle_off.inp` | Baseline/no-op companion case — `peleLM.use_inlet_from_plane = 0`. | Reference plt; this file exists to confirm the new code path is a strict no-op when disabled. |
| `eb_bfs_recycle_basic.inp` | Happy path: plane at `x = 0.03`, no warmup, cumulative mean. | Snapshot, mean update, and `fillFromRecyclingPlane` ParallelAdd onto the lo face. |
| `eb_bfs_recycle_warmup.inp` | Warmup gate: `inlet_plane_warmup_steps = 15`, `max_step = 10`. | Snapshot path runs every step, but no fluctuation is ever injected. |
| `eb_bfs_recycle_window.inp` | EMA branch: `inlet_plane_avg_window = 1.0e-4`. | Forces `α = min(1, dt/window)` rather than the cumulative `1/n` fallback. |
| `eb_bfs_recycle_amr.inp` | AMR regrid: `max_level = 1` with a static refinement box over the recycling plane and `num_init_iter = 1`. | Exercises `buildRecyclingPlaneStorage` across regrids, the coarse → fine `InterpFromCoarseLevel` slab interpolation, and the `m_recyclingNeedsRebuild` hooks in `MakeNewLevelFromCoarse` / `RemakeLevel` / `ClearLevel`. |
| `eb_bfs_recycle_restart.inp` | Checkpoint round-trip via a two-phase CTest command (writes `chk00005`, restarts from it through step 10). | Exercises the `RecyclingPlane:` block in the checkpoint header and the per-level `recycle_mean` MultiFab read/write. |

## Invariants verified by `fcompare`

The smoke tests produce three relationships that act as functional checks
beyond "does it run":

1. **`off` ≡ `warmup` (bitwise identical).** With the warmup gate
   suppressing injection, the snapshot path still runs but never touches
   the state; `plt00010` from the two cases must agree to the last bit.
2. **`off` ≠ `basic` (recycling does inject).** Adding the cached
   fluctuation on top of the `ext_dir` inflow produces a measurable
   difference in the inlet-region velocity (~1e-4 m/s after 10 steps in
   our local check).
3. **`basic` ≈ `restart` (checkpoint round-trip).** The two-phase
   restart variant continues from `chk00005` and reaches `plt00010`
   matching the single-shot `basic` run to within plotfile I/O roundoff
   (~1e-7 relative).

## Running

Build the case in this directory the usual way (`make TPL && make`),
then invoke any of the inputs directly, e.g.

```
mpiexec -n 4 ./PeleLMeX2d.gnu.MPI.ex eb_bfs_recycle_basic.inp
```

The two-phase restart test cannot be reproduced by a single binary call;
the registered CTest entry runs the executable twice with phase-specific
overrides. To replicate manually:

```
mpiexec -n 4 ./PeleLMeX2d.gnu.MPI.ex eb_bfs_recycle_restart.inp \
    amr.max_step=5 amr.check_int=5 amr.plot_int=-1
mpiexec -n 4 ./PeleLMeX2d.gnu.MPI.ex eb_bfs_recycle_restart.inp \
    amr.restart=chk00005 amr.max_step=10 amr.check_int=-1
```

## Note on the AMR variant

The base `eb_bfs.inp` uses `peleLM.num_init_iter = 3`. At `max_level = 1`
on this geometry, three successive initial pressure iterations drive the
MAC projection divergent before the first real step (independent of the
recycling code — the same failure reproduces with
`use_inlet_from_plane = 0`). The AMR test therefore overrides
`num_init_iter = 1` and uses a static refinement box clear of the EB,
which is enough to exercise every recycling-plane regrid hook.
