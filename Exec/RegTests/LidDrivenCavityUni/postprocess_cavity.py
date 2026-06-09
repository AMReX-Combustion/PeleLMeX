#!/usr/bin/env python3
"""
Post-processor for the 2D lid-driven cavity.

Reproduces the centerline-profile validation plots from Ghia, Ghia &
Shin, *J. Comp. Phys.* 48 (1982) 387-411, for the steady cavity flow
at Re = 100, 400, 1000, 3200, 5000, 7500, and 10000.  Designed to
handle plotfiles produced by the ``input.2d`` inputs file.

Workflow
========

1. Load the (one) plotfile via yt and pull u and v on a level-0
   covering grid.
2. Linearly interpolate u along the vertical centerline x = 0.5 and v
   along the horizontal centerline y = 0.5.
3. Overlay against Ghia tabulated points and save the comparison
   plots + a CSV.

Usage
=====

::

    python3 postprocess_cavity.py path/to/plt_NNNNN \\
        --Re 100 --output ./cavity_Re100_profiles

The CLI ``--Re`` selects which Ghia reference column is overlaid.
Outputs (in ``--output`` dir):

  - ``profiles.csv``           : u(y) and v(x) interpolated to the
                                 centerlines, plus the Ghia reference
                                 points for the chosen Re.
  - ``Uvertical.png``          : u along x=0.5 vs Ghia.
  - ``Vhorizontal.png``        : v along y=0.5 vs Ghia.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np


# --------------------------------------------------------------------
# Ghia et al. (1982) tabulated centerline data.
# Tables I and II for the three steady Re values produced by the
# bundled input files.  Columns: y or x, then u or v normalized by
# U_lid.  Values are typed in from the original publication.
# --------------------------------------------------------------------

GHIA_U_VERT = {
    # Table I: u(x = L/2, y), normalized by U_lid.
    100: [
        (0.0000, 0.00000), (0.0547, -0.03717), (0.0625, -0.04192),
        (0.0703, -0.04775), (0.1016, -0.06434), (0.1719, -0.10150),
        (0.2813, -0.15662), (0.4531, -0.21090), (0.5000, -0.20581),
        (0.6172, -0.13641), (0.7344, 0.00332),  (0.8516, 0.23151),
        (0.9531, 0.68717),  (0.9609, 0.73722),  (0.9688, 0.78871),
        (0.9766, 0.84123),  (1.0000, 1.00000),
    ],
    400: [
        (0.0000, 0.00000), (0.0547, -0.08186), (0.0625, -0.09266),
        (0.0703, -0.10338), (0.1016, -0.14612), (0.1719, -0.24299),
        (0.2813, -0.32726), (0.4531, -0.17119), (0.5000, -0.11477),
        (0.6172, 0.02135),  (0.7344, 0.16256),  (0.8516, 0.29093),
        (0.9531, 0.55892),  (0.9609, 0.61756),  (0.9688, 0.68439),
        (0.9766, 0.75837),  (1.0000, 1.00000),
    ],
    1000: [
        (0.0000, 0.00000), (0.0547, -0.18109), (0.0625, -0.20196),
        (0.0703, -0.22220), (0.1016, -0.29730), (0.1719, -0.38289),
        (0.2813, -0.27805), (0.4531, -0.10648), (0.5000, -0.06080),
        (0.6172, 0.05702),  (0.7344, 0.18719),  (0.8516, 0.33304),
        (0.9531, 0.46604),  (0.9609, 0.51117),  (0.9688, 0.57492),
        (0.9766, 0.65928),  (1.0000, 1.00000),
    ],
}

GHIA_V_HORIZ = {
    # Table II: v(x, y = L/2), normalized by U_lid.
    100: [
        (0.0000, 0.00000), (0.0625, 0.09233),  (0.0703, 0.10091),
        (0.0781, 0.10890), (0.0938, 0.12317),  (0.1563, 0.16077),
        (0.2266, 0.17507), (0.2344, 0.17527),  (0.5000, 0.05454),
        (0.8047, -0.24533),(0.8594, -0.22445), (0.9063, -0.16914),
        (0.9453, -0.10313),(0.9531, -0.08864), (0.9609, -0.07391),
        (0.9688, -0.05906),(1.0000, 0.00000),
    ],
    400: [
        (0.0000, 0.00000), (0.0625, 0.18360),  (0.0703, 0.19713),
        (0.0781, 0.20920), (0.0938, 0.22965),  (0.1563, 0.28124),
        (0.2266, 0.30203), (0.2344, 0.30174),  (0.5000, 0.05186),
        (0.8047, -0.38598),(0.8594, -0.44993), (0.9063, -0.33827),
        (0.9453, -0.22847),(0.9531, -0.19254), (0.9609, -0.15663),
        (0.9688, -0.12146),(1.0000, 0.00000),
    ],
    1000: [
        (0.0000, 0.00000), (0.0625, 0.27485),  (0.0703, 0.29012),
        (0.0781, 0.30353), (0.0938, 0.32627),  (0.1563, 0.37095),
        (0.2266, 0.33075), (0.2344, 0.32235),  (0.5000, 0.02526),
        (0.8047, -0.31966),(0.8594, -0.42665), (0.9063, -0.51550),
        (0.9453, -0.39188),(0.9531, -0.33714), (0.9609, -0.27669),
        (0.9688, -0.21388),(1.0000, 0.00000),
    ],
}


def axis_phys_positions(prob_lo: float, prob_hi: float, N: int,
                        beta: float, at: str) -> np.ndarray:
    """
    Compute physical positions along a stretched axis.

    Parameters
    ----------
    at : {"cc", "lo_node", "hi_node"}
        "cc"      : N cell-centers (indices 0..N-1, offset 0.5)
        "lo_node" : N+1 nodes; useful for face-aligned positions if
                    you ever need them
    """
    L = prob_hi - prob_lo
    if at == "cc":
        idx = np.arange(N) + 0.5
        n_out = N
    elif at == "lo_node":
        idx = np.arange(N + 1).astype(float)
        n_out = N + 1
    else:
        raise ValueError(f"unknown 'at': {at!r}")
    eta = idx / N
    return prob_lo + L * eta


# --------------------------------------------------------------------
# Plotfile loader
# --------------------------------------------------------------------


def load_uv_grid(plt_path: Path):
    """
    Load (u, v) from a 2D plotfile.  Returns:
      u, v : np.ndarray, shape (Nx, Ny)
      prob_lo, prob_hi : tuples of length 2
    """
    import yt  # type: ignore

    yt.set_log_level("error")
    p = Path(plt_path)
    if p.is_file() and p.name == "Header":
        p = p.parent
    ds = yt.load(str(p))
    cg = ds.covering_grid(
        level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
    u = np.asarray(cg["boxlib", "x_velocity"])
    v = np.asarray(cg["boxlib", "y_velocity"])
    # In 2D yt's arrays are (Nx, Ny, 1); squeeze the trailing dim.
    if u.ndim == 3 and u.shape[2] == 1:
        u = u[..., 0]
        v = v[..., 0]
    plo = (float(ds.domain_left_edge[0]),  float(ds.domain_left_edge[1]))
    phi = (float(ds.domain_right_edge[0]), float(ds.domain_right_edge[1]))
    return u, v, plo, phi


# --------------------------------------------------------------------
# Centerline interpolation
# --------------------------------------------------------------------


def interp_along_axis(values_2d: np.ndarray, x_axis: np.ndarray,
                      y_axis: np.ndarray, axis: str,
                      target: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Linearly interpolate values_2d along `axis` to the line where the
    perpendicular coordinate equals `target`.  Returns (coords_along,
    values_along).

    axis = "y" : extract along y (return y_axis), interpolating in x
                 between two columns straddling x = target.
    axis = "x" : extract along x (return x_axis), interpolating in y
                 between two rows straddling y = target.
    """
    Nx, Ny = values_2d.shape
    if axis == "y":
        # find columns i0, i1 with x_axis[i0] <= target <= x_axis[i1]
        if not (x_axis[0] <= target <= x_axis[-1]):
            raise ValueError(
                f"target x={target} outside [{x_axis[0]}, {x_axis[-1]}]")
        i1 = int(np.searchsorted(x_axis, target, side="left"))
        if i1 == 0:
            i0, i1 = 0, 1
        else:
            i0 = i1 - 1
        # if exactly on a cell center
        if x_axis[i1] == target:
            return y_axis.copy(), values_2d[i1, :].copy()
        w = (target - x_axis[i0]) / (x_axis[i1] - x_axis[i0])
        return y_axis.copy(), (1.0 - w) * values_2d[i0, :] + w * values_2d[i1, :]
    elif axis == "x":
        if not (y_axis[0] <= target <= y_axis[-1]):
            raise ValueError(
                f"target y={target} outside [{y_axis[0]}, {y_axis[-1]}]")
        j1 = int(np.searchsorted(y_axis, target, side="left"))
        if j1 == 0:
            j0, j1 = 0, 1
        else:
            j0 = j1 - 1
        if y_axis[j1] == target:
            return x_axis.copy(), values_2d[:, j1].copy()
        w = (target - y_axis[j0]) / (y_axis[j1] - y_axis[j0])
        return x_axis.copy(), (1.0 - w) * values_2d[:, j0] + w * values_2d[:, j1]
    else:
        raise ValueError(f"axis must be 'x' or 'y' (got {axis!r})")


# --------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------


def plot_u_vert(y_sim, u_sim, ghia_pts, outpath, Re: int):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot(u_sim, y_sim, "-", label=f"PeleLMeX", linewidth=2)
    if ghia_pts is not None:
        gy, gu = zip(*ghia_pts)
        ax.plot(gu, gy, "ko", label=f"Ghia et al. (1982), Re={Re}",
                markersize=5, markerfacecolor="none")
    ax.set_xlabel(r"$u / U_{\rm lid}$")
    ax.set_ylabel(r"$y / L$")
    ax.set_xlim(-0.6, 1.05)
    ax.set_ylim(0.0, 1.0)
    ax.axvline(0.0, color="k", linewidth=0.5, alpha=0.4)
    ax.grid(True, alpha=0.3)
    ax.set_title(f"u along x = 0.5 (Re = {Re})")
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def plot_v_horiz(x_sim, v_sim, ghia_pts, outpath, Re: int):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot(x_sim, v_sim, "-", label="PeleLMeX", linewidth=2)
    if ghia_pts is not None:
        gx, gv = zip(*ghia_pts)
        ax.plot(gx, gv, "ko", label=f"Ghia et al. (1982), Re={Re}",
                markersize=5, markerfacecolor="none")
    ax.set_xlabel(r"$x / L$")
    ax.set_ylabel(r"$v / U_{\rm lid}$")
    ax.set_xlim(0.0, 1.0)
    ax.axhline(0.0, color="k", linewidth=0.5, alpha=0.4)
    ax.grid(True, alpha=0.3)
    ax.set_title(f"v along y = 0.5 (Re = {Re})")
    ax.legend(loc="lower left", frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------
# Main
# --------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotfile",
                    help="path to a PeleLMeX plotfile directory or its "
                         "Header file (use the final-time plotfile)")
    ap.add_argument("--Re", type=int, default=100,
                    help="Reynolds number (selects Ghia reference column; "
                         "supported: 100, 400, 1000; default 100)")
    ap.add_argument("--U-lid", type=float, default=1.0,
                    help="lid velocity used to normalize the profiles "
                         "(default 1.0; matches the shipped inputs)")
    ap.add_argument("--output", default="./postprocess_cavity_out",
                    help="output directory")
    args = ap.parse_args()

    # ---- Load ----
    u, v, plo, phi = load_uv_grid(Path(args.plotfile))
    Nx, Ny = u.shape

    # ---- Cell-center physical coordinates ----
    x_phys = plo[0] + (np.arange(Nx) + 0.5) * (phi[0] - plo[0]) / Nx
    y_phys = plo[1] + (np.arange(Ny) + 0.5) * (phi[1] - plo[1]) / Ny

    # ---- Centerlines ----
    L_x = phi[0] - plo[0]
    L_y = phi[1] - plo[1]
    x_mid = 0.5 * (plo[0] + phi[0])
    y_mid = 0.5 * (plo[1] + phi[1])
    # Vertical centerline:  u along x = x_mid
    y_vert, u_vert = interp_along_axis(u, x_phys, y_phys, "y", x_mid)
    # Horizontal centerline:  v along y = y_mid
    x_horiz, v_horiz = interp_along_axis(v, x_phys, y_phys, "x", y_mid)

    # Normalize for plotting / comparison
    Ulid = args.U_lid
    u_vert_norm = u_vert / Ulid
    v_horiz_norm = v_horiz / Ulid
    y_vert_norm = (y_vert - plo[1]) / L_y
    x_horiz_norm = (x_horiz - plo[0]) / L_x

    # ---- Reference ----
    ghia_u = GHIA_U_VERT.get(args.Re)
    ghia_v = GHIA_V_HORIZ.get(args.Re)
    if ghia_u is None:
        print(f"  no bundled Ghia data for Re={args.Re} "
              "(supported: 100, 400, 1000); plots will skip the overlay")

    # ---- Output ----
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    csv_path = outdir / "profiles.csv"
    with csv_path.open("w") as fh:
        fh.write("# Lid-driven cavity centerline profiles "
                 f"(Re = {args.Re})\n")
        fh.write(f"# U_lid = {Ulid}, "
                 f"domain = [{plo[0]},{phi[0]}] x [{plo[1]},{phi[1]}]\n")
        fh.write("#\n# Section 1: u along x = 0.5 (vertical centerline)\n")
        fh.write("# y_norm, u_PeleLMeX, u_Ghia (interp at y_norm if avail)\n")
        for yi, ui in zip(y_vert_norm, u_vert_norm):
            ug = ""
            if ghia_u is not None:
                gy = np.array([p[0] for p in ghia_u])
                gu = np.array([p[1] for p in ghia_u])
                if gy[0] <= yi <= gy[-1]:
                    ug = f"{np.interp(yi, gy, gu):.6f}"
            fh.write(f"{yi:.6f}, {ui:.6f}, {ug}\n")
        fh.write("#\n# Section 2: v along y = 0.5 (horizontal centerline)\n")
        fh.write("# x_norm, v_PeleLMeX, v_Ghia (interp at x_norm if avail)\n")
        for xi, vi in zip(x_horiz_norm, v_horiz_norm):
            vg = ""
            if ghia_v is not None:
                gx = np.array([p[0] for p in ghia_v])
                gv = np.array([p[1] for p in ghia_v])
                if gx[0] <= xi <= gx[-1]:
                    vg = f"{np.interp(xi, gx, gv):.6f}"
            fh.write(f"{xi:.6f}, {vi:.6f}, {vg}\n")
    print(f"  wrote {csv_path}")

    plot_u_vert(y_vert_norm, u_vert_norm, ghia_u,
                outdir / "Uvertical.png", args.Re)
    plot_v_horiz(x_horiz_norm, v_horiz_norm, ghia_v,
                 outdir / "Vhorizontal.png", args.Re)
    print(f"  wrote {outdir / 'Uvertical.png'}")
    print(f"  wrote {outdir / 'Vhorizontal.png'}")

    # ---- L2 / Linf vs Ghia, where available ----
    if ghia_u is not None:
        gy = np.array([p[0] for p in ghia_u])
        gu = np.array([p[1] for p in ghia_u])
        # Interpolate sim onto Ghia y-points
        u_sim_at_gy = np.interp(gy, y_vert_norm, u_vert_norm)
        err = u_sim_at_gy - gu
        print(f"\n  u-profile vs Ghia (Re={args.Re}):  "
              f"L2 = {np.sqrt(np.mean(err**2)):.4e}, "
              f"Linf = {np.max(np.abs(err)):.4e}")
    if ghia_v is not None:
        gx = np.array([p[0] for p in ghia_v])
        gv = np.array([p[1] for p in ghia_v])
        v_sim_at_gx = np.interp(gx, x_horiz_norm, v_horiz_norm)
        err = v_sim_at_gx - gv
        print(f"  v-profile vs Ghia (Re={args.Re}):  "
              f"L2 = {np.sqrt(np.mean(err**2)):.4e}, "
              f"Linf = {np.max(np.abs(err)):.4e}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
