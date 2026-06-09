# LidDrivenCavityUni

Classical 2D lid-driven cavity benchmark, used here to verify the
proper handling of Dirichlet boundary conditions on the velocity.

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

