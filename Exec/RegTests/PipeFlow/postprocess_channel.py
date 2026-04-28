#!/usr/bin/env python3
"""
Post-processor for channel-flow plotfiles.

Reproduces the U+ vs y+ and velocity-rms vs y+ plots from
``Docs/sphinx/manual/Validation.rst`` (the "Channel Flow using EB"
section, Re_tau = 180 case).  Designed to handle plotfiles produced
by the stretched-channel input ``input.3d-channel-Re180-stretched``,
the EB validation case in ``Exec/RegTests/EB_PipeFlow``, and uniform-
grid runs (no mesh mapping).

Workflow
========

1. Load each plotfile via yt's covering_grid and pull (u, v, w) as
   shape ``(Nx, Ny, Nz, 3)`` arrays.
2. Spatially average in the periodic directions (x, z by default) to
   get U_mean(y), V_mean(y), W_mean(y) plus the variances of each
   component.
3. Time-average those per-y profiles across all input plotfiles.
4. Compute the friction velocity ``u_tau`` either analytically from
   the imposed pressure gradient (``u_tau = sqrt(|gradP| * delta /
   rho)``, the cleanest choice for a periodic-channel run) or
   numerically from the wall stress mu*dU/dy (``--utau-source
   wall_stress``).
5. Wall-units rescale and plot.

Mesh-mapping awareness
======================

To convert AMReX index space (Xi-space) to physical y, supply mapping
parameters via CLI flags.  Without flags the script assumes uniform
grid (Xi == physical).  ConstantMap and ExpStretchMap are supported.

Examples
========

Stretched-channel run (matches input.3d-channel-Re180-stretched)::

    python3 postprocess_channel.py \\
        plt_*/Header \\
        --map exp --beta 4 --wall-end hi \\
        --delta 0.005 --mu 3.578e-5 --rho 0.4688 \\
        --gradP 709.79 --half-channel --output ./profiles_Re180

Uniform-grid run::

    python3 postprocess_channel.py plt_*/Header --delta 0.005 \\
        --mu 3.578e-5 --rho 0.4688 --gradP 709.79

Outputs (in --output dir, default ``./postprocess_channel_out``):

  - profiles.csv   : tabular y_phys, y+, U+, urms+, vrms+, wrms+
  - Uplus.png      : U+ vs y+ on a semilog x axis with viscous
                     sublayer (u+ = y+) and log law (u+ = 2.5*ln(y+)+5.5)
                     reference lines
  - VelRMSplus.png : urms+, vrms+, wrms+ vs y+
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------
# Mesh-mapping helpers
# --------------------------------------------------------------------


@dataclass
class MeshMapInfo:
    kind: str = "none"          # "none" | "constant" | "exp"
    direction: int = 1          # 0/1/2; default y
    wall_end: str = "hi"        # for "exp": "lo" or "hi"
    beta: float = 0.0           # for "exp"
    fac: tuple = (1.0, 1.0, 1.0)  # for "constant"


def _eta_of_xi(xi: np.ndarray, prob_lo: float, prob_hi: float,
               wall_end: str) -> np.ndarray:
    L_xi = prob_hi - prob_lo
    if wall_end == "lo":
        return (xi - prob_lo) / L_xi
    return (prob_hi - xi) / L_xi


def y_physical_from_xi(xi: np.ndarray, prob_lo: float, prob_hi: float,
                       mapping: MeshMapInfo) -> np.ndarray:
    """Map index-space coord to physical coord along the stretched axis."""
    if mapping.kind == "none":
        return xi
    if mapping.kind == "constant":
        return prob_lo + (xi - prob_lo) * mapping.fac[mapping.direction]
    if mapping.kind == "exp":
        L_xi = prob_hi - prob_lo
        eta = _eta_of_xi(xi, prob_lo, prob_hi, mapping.wall_end)
        beta = mapping.beta
        if abs(beta) < 1.0e-10:
            offset = eta
        else:
            offset = (np.exp(beta * eta) - 1.0) / (math.expm1(beta))
        if mapping.wall_end == "lo":
            return prob_lo + L_xi * offset
        return prob_hi - L_xi * offset
    raise ValueError(f"unknown mapping kind '{mapping.kind}'")


# --------------------------------------------------------------------
# Plotfile loading
# --------------------------------------------------------------------


def load_plotfile_velocity(plt_path: Path):
    """
    Load (u, v, w) on a covering_grid from a plotfile, plus geometry
    metadata.  Returns:
      u, v, w : np.ndarray, shape (Nx, Ny, Nz)
      prob_lo, prob_hi : tuples of length 3
      time : float
    """
    import yt  # type: ignore

    yt.set_log_level("error")
    # Allow user to pass the plotfile dir or the Header inside it.
    p = Path(plt_path)
    if p.is_file() and p.name == "Header":
        p = p.parent
    ds = yt.load(str(p))
    cg = ds.covering_grid(level=0,
                          left_edge=ds.domain_left_edge,
                          dims=ds.domain_dimensions)
    u = np.asarray(cg["boxlib", "x_velocity"])
    v = np.asarray(cg["boxlib", "y_velocity"])
    try:
        w = np.asarray(cg["boxlib", "z_velocity"])
    except Exception:
        w = np.zeros_like(u)
    prob_lo = tuple(float(x) for x in ds.domain_left_edge)
    prob_hi = tuple(float(x) for x in ds.domain_right_edge)
    return u, v, w, prob_lo, prob_hi, float(ds.current_time)


# --------------------------------------------------------------------
# Per-plotfile y-binned statistics
# --------------------------------------------------------------------


def per_y_stats(u, v, w, periodic_axes=(0, 2)):
    """
    Spatially average over periodic axes.  Returns dicts of arrays
    indexed by the wall-normal index j (length Ny):
      mean: U_bar, V_bar, W_bar
      var : U_var, V_var, W_var       (variance over the (x,z) plane)
    """
    U_bar = u.mean(axis=periodic_axes)
    V_bar = v.mean(axis=periodic_axes)
    W_bar = w.mean(axis=periodic_axes)
    U_var = u.var(axis=periodic_axes)
    V_var = v.var(axis=periodic_axes)
    W_var = w.var(axis=periodic_axes)
    return {"U_bar": U_bar, "V_bar": V_bar, "W_bar": W_bar,
            "U_var": U_var, "V_var": V_var, "W_var": W_var}


# --------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------


def plot_Uplus(yplus, Uplus, outpath):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.semilogx(yplus, Uplus, "o-", label="PeleLMeX", markersize=4)
    yp_visc = np.logspace(-1, 1, 50)
    ax.semilogx(yp_visc, yp_visc, "k--", alpha=0.6,
                label=r"$u^+ = y^+$ (viscous sublayer)")
    yp_log = np.logspace(0.7, 2.7, 50)
    ax.semilogx(yp_log, 2.5 * np.log(yp_log) + 5.5, "k:",
                alpha=0.6,
                label=r"$u^+ = 2.5\,\ln y^+ + 5.5$ (log law)")
    ax.set_xlabel(r"$y^+$")
    ax.set_ylabel(r"$U^+$")
    ax.set_xlim(0.1, max(200.0, yplus.max()))
    ax.set_ylim(0, max(20.0, Uplus.max() * 1.05))
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower right", frameon=False)
    ax.set_title(r"Mean streamwise velocity in wall units")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def plot_VelRMSplus(yplus, urms_p, vrms_p, wrms_p, outpath):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(yplus, urms_p, "o-", label=r"$u'_{rms}/u_\tau$", markersize=4)
    ax.plot(yplus, vrms_p, "s-", label=r"$v'_{rms}/u_\tau$", markersize=4)
    ax.plot(yplus, wrms_p, "^-", label=r"$w'_{rms}/u_\tau$", markersize=4)
    ax.set_xlabel(r"$y^+$")
    ax.set_ylabel(r"$\sqrt{\overline{u_i'\,u_i'}}/u_\tau$")
    ax.set_xlim(0, max(200.0, yplus.max()))
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", frameon=False)
    ax.set_title(r"Velocity fluctuation RMS in wall units")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------
# Main
# --------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotfiles", nargs="+",
                    help="one or more PeleLMeX plotfile directories or "
                         "their Header files; statistics are averaged "
                         "over the set")
    ap.add_argument("--delta", type=float, default=0.005,
                    help="channel half-width [m] (default: 0.005)")
    ap.add_argument("--mu", type=float, default=3.57816e-5,
                    help="dynamic viscosity [Pa s] (default: docs Re180)")
    ap.add_argument("--rho", type=float, default=0.4688,
                    help="density [kg/m^3] (default: docs Re180)")
    ap.add_argument("--gradP", type=float, default=709.79,
                    help="magnitude of imposed |dP/dx| [Pa/m] (default: "
                         "docs Re180; only used when --utau-source=gradP)")
    ap.add_argument("--utau-source", choices=("gradP", "wall_stress"),
                    default="gradP",
                    help="how to compute u_tau (default: gradP, the "
                         "analytic balance for a periodic channel)")
    ap.add_argument("--map", choices=("none", "constant", "exp"),
                    default="none",
                    help="mesh mapping kind (default: none, uniform grid)")
    ap.add_argument("--map-direction", type=int, default=1,
                    help="stretched axis (0=x, 1=y, 2=z; default 1=y)")
    ap.add_argument("--wall-end", choices=("lo", "hi"), default="hi",
                    help="for exp map: which end of the stretched axis "
                         "is the wall (default 'hi' to match the half-"
                         "channel input)")
    ap.add_argument("--beta", type=float, default=0.0,
                    help="for exp map: stretching strength (default 0 "
                         "= identity)")
    ap.add_argument("--fac", type=float, nargs=3, default=(1.0, 1.0, 1.0),
                    metavar=("FACX", "FACY", "FACZ"),
                    help="for constant map: per-direction scaling "
                         "(default 1 1 1)")
    ap.add_argument("--half-channel", action="store_true",
                    help="treat the y domain as a half-channel (y from "
                         "wall to symmetry plane).  Statistics are "
                         "reported for the full y range present in the "
                         "data; this flag only affects whether the "
                         "wall-distance for y+ is measured from y_lo or "
                         "y_hi.  Pair with --wall-end so they agree.")
    ap.add_argument("--periodic-axes", type=int, nargs="+", default=(0, 2),
                    help="indices of axes to spatially average over "
                         "(default: 0 2, i.e. x and z)")
    ap.add_argument("--output", default="./postprocess_channel_out",
                    help="output directory")
    args = ap.parse_args()

    mapping = MeshMapInfo(
        kind=args.map,
        direction=args.map_direction,
        wall_end=args.wall_end,
        beta=args.beta,
        fac=tuple(args.fac),
    )

    # --- Load and accumulate ----------------------------------------------
    accum = None
    accum_count = 0
    geom = None
    for pf in args.plotfiles:
        u, v, w, plo, phi, t = load_plotfile_velocity(Path(pf))
        if geom is None:
            geom = (u.shape, plo, phi)
        else:
            if u.shape != geom[0]:
                print(f"shape mismatch for {pf}: {u.shape} vs {geom[0]}",
                      file=sys.stderr)
                return 2
        stats = per_y_stats(u, v, w,
                            periodic_axes=tuple(args.periodic_axes))
        if accum is None:
            accum = {k: v.copy() for k, v in stats.items()}
        else:
            for k in accum:
                accum[k] += stats[k]
        accum_count += 1
        print(f"  loaded {pf}  t={t:.4e}  shape={u.shape}")

    if accum is None or accum_count == 0:
        print("no plotfiles loaded", file=sys.stderr)
        return 2

    for k in accum:
        accum[k] /= accum_count

    Ny = accum["U_bar"].size
    shape, plo, phi = geom

    # --- Cell-center physical y -----------------------------------------
    j = np.arange(Ny)
    dxi_y = (phi[1] - plo[1]) / shape[1]
    xi_cc = plo[1] + (j + 0.5) * dxi_y
    y_phys = y_physical_from_xi(xi_cc, plo[1], phi[1], mapping)

    # Wall distance.  For the half-channel input, the wall is at
    # y = phi[1] (high y).  Without --half-channel we treat the
    # smaller |y| from the centerline as the wall distance, which is
    # natural for full-channel runs centered at y=0.
    if args.half_channel:
        if args.wall_end == "hi":
            yw = phi[1] - y_phys
        else:
            yw = y_phys - plo[1]
    else:
        # Full channel: wall is at the high or low y boundary.  Use
        # the closer boundary cell-by-cell (i.e. min of distance to
        # either wall).
        yw = np.minimum(y_phys - plo[1], phi[1] - y_phys)

    # --- u_tau -----------------------------------------------------------
    if args.utau_source == "gradP":
        u_tau = math.sqrt(abs(args.gradP) * args.delta / args.rho)
    else:
        # use the velocity gradient at the wall:  tau_w = mu * dU/dy
        # Find first off-wall cell (yw smallest) and approximate dU/dy
        # by U_bar(j)/yw(j).
        order = np.argsort(yw)
        j0 = order[0]
        # Use linear extrapolation U(0) = 0 (no-slip wall)
        dU_dy_wall = accum["U_bar"][j0] / yw[j0]
        tau_w = args.mu * abs(dU_dy_wall)
        u_tau = math.sqrt(tau_w / args.rho)
    nu = args.mu / args.rho
    print(f"\n  u_tau = {u_tau:.5g}  (source: {args.utau_source})")
    print(f"  nu    = {nu:.5g}")
    print(f"  delta = {args.delta:.5g}")
    print(f"  Re_tau = {u_tau * args.delta / nu:.3f}")

    # --- Wall-units rescale ---------------------------------------------
    yplus = yw * u_tau / nu
    Uplus = accum["U_bar"] / u_tau
    urms_p = np.sqrt(np.maximum(accum["U_var"], 0.0)) / u_tau
    vrms_p = np.sqrt(np.maximum(accum["V_var"], 0.0)) / u_tau
    wrms_p = np.sqrt(np.maximum(accum["W_var"], 0.0)) / u_tau

    # Sort by y+ so the plot lines run wall-to-center cleanly
    order = np.argsort(yplus)
    y_phys_o = y_phys[order]
    yplus_o = yplus[order]
    Uplus_o = Uplus[order]
    urms_o = urms_p[order]
    vrms_o = vrms_p[order]
    wrms_o = wrms_p[order]

    # --- Output ----------------------------------------------------------
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    csv_path = outdir / "profiles.csv"
    with csv_path.open("w") as fh:
        fh.write("# y_phys, y_plus, U_plus, urms_plus, vrms_plus, wrms_plus\n")
        fh.write(f"# u_tau = {u_tau:.6g}  nu = {nu:.6g}  "
                 f"Re_tau = {u_tau * args.delta / nu:.3f}\n")
        for k in range(Ny):
            fh.write(f"{y_phys_o[k]:.6e}, {yplus_o[k]:.6e}, "
                     f"{Uplus_o[k]:.6e}, {urms_o[k]:.6e}, "
                     f"{vrms_o[k]:.6e}, {wrms_o[k]:.6e}\n")
    print(f"\n  wrote {csv_path}")

    plot_Uplus(yplus_o, Uplus_o, outdir / "Uplus.png")
    plot_VelRMSplus(yplus_o, urms_o, vrms_o, wrms_o,
                    outdir / "VelRMSplus.png")
    print(f"  wrote {outdir / 'Uplus.png'}")
    print(f"  wrote {outdir / 'VelRMSplus.png'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
