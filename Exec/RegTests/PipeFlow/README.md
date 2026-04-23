# PipeFlow (non-EB channel flow)

A rectangular-channel reformulation of `EB_PipeFlow` that does **not** use
embedded boundaries.  The domain is a 3D box, periodic in the flow
direction (`x`) with no-slip walls on the `y` and `z` faces.

This case mirrors `EB_PipeFlow` as closely as possible so results can be
compared between the two; it exists primarily as the target for
radial/wall-normal mesh stretching (mesh mapping) in the incompressible
low-Mach solver, without the added complication of EB.

## Running

```
make TPL
make -j
./PeleLMeX3d.gnu.ex input.3d-Poiseuille
```

## Parameters

- `prob.meanFlowMag` — bulk flow speed used in the uniform-flow IC
- `prob.meanFlowDir` — Cartesian index of the mean flow direction (1 = x)
- `prob.problem_type` — `1` for uniform-mean + divergence-free
  perturbations (default), `2` for a turbulent-channel-like mean profile

## Differences from `EB_PipeFlow`

- `USE_EB = FALSE` in the `GNUmakefile` / `CMakeLists.txt`
- No `eb2.*` entries in the input file
- Walls implemented as Cartesian `NoSlipWallAdiab` BCs on both `y` and
  `z` faces rather than an EB cylinder
