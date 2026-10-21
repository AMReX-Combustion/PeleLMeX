#!/usr/bin/env python3
"""Compare a DiagFramePlane slice from a mesh-mapped run against the same
slice from an unmapped run, to check that turbulent inflow is sampled at
PHYSICAL transverse positions.

usage:
  compare_inflow_plane.py <uniform_plane_dir> <mapped_plane_dir> \
      [--map tanh|constant] [--beta 2.0] [--field x_velocity] \
      [--require 0.99] [--plot out.png]

Both plane directories are 2D AMReX plotfiles written by PelePhysics'
DiagFramePlane (normal = 2, so plane axes are x, y).  The mapped run's plane
carries Xi (computational) coordinates; the script maps them to physical x
(TanhStretchMap or ConstantMap in x, parameter --beta) and interpolates the
uniform run onto those positions.  Each run's x extent comes from its own
Header, so the two domains may differ in size (ConstantMap).  Two
hypotheses are scored:

  H_phys : mapped(x_phys_i) == uniform(x_phys_i)   <- correct sampling
  H_xi   : mapped(x_phys_i) == uniform(xi_i)       <- the pre-fix bug

Correlation near 1 for H_phys and clearly lower for H_xi is a pass.  With
--require C the script exits non-zero unless the mean H_phys correlation is
at least C and exceeds the H_xi one by at least 0.2, so it can gate CI.
"""
import argparse
import os
import re
import sys

import numpy as np


# --------------------------------------------------------------------------
# Minimal reader for a single-level 2D AMReX plotfile (VisMF v1, native double)
# --------------------------------------------------------------------------
def read_plotfile_header(pltdir):
    with open(os.path.join(pltdir, "Header")) as f:
        lines = [ln.rstrip("\n") for ln in f]
    nvar = int(lines[1])
    varnames = lines[2 : 2 + nvar]
    p = 2 + nvar
    spacedim = int(lines[p]); p += 1
    time = float(lines[p]); p += 1
    finest = int(lines[p]); p += 1
    problo = [float(v) for v in lines[p].split()]; p += 1
    probhi = [float(v) for v in lines[p].split()]; p += 1
    p += 1  # ref ratios
    dom = re.findall(r"\(\((-?\d+),(-?\d+)(?:,-?\d+)?\) \((-?\d+),(-?\d+)(?:,-?\d+)?\)", lines[p])[0]
    domain = tuple(int(v) for v in dom)
    return dict(varnames=varnames, spacedim=spacedim, time=time, finest=finest,
                problo=problo, probhi=probhi, domain=domain)


def read_level0(pltdir, ncomp_expected=None):
    lev = os.path.join(pltdir, "Level_0")
    with open(os.path.join(lev, "Cell_H")) as f:
        txt = f.read()
    ncomp = int(txt.splitlines()[2])
    boxes = re.findall(r"\(\((-?\d+),(-?\d+)(?:,-?\d+)?\) \((-?\d+),(-?\d+)(?:,-?\d+)?\)", txt)
    fabs = re.findall(r"FabOnDisk:\s+(\S+)\s+(\d+)", txt)
    if len(fabs) == 0:
        sys.exit("no FabOnDisk entries in " + os.path.join(lev, "Cell_H"))
    boxes = boxes[: len(fabs)]  # box list precedes the FabOnDisk list
    return ncomp, [tuple(int(v) for v in b) for b in boxes], fabs


def read_fab(path, offset, ncomp):
    with open(path, "rb") as f:
        f.seek(offset)
        hdr = b""
        while not hdr.endswith(b"\n"):
            c = f.read(1)
            if not c:
                sys.exit("unexpected EOF reading FAB header in " + path)
            hdr += c
        h = hdr.decode()
        m = re.search(r"\(\((-?\d+),(-?\d+)(?:,-?\d+)?\) \((-?\d+),(-?\d+)(?:,-?\d+)?\)", h)
        ilo, jlo, ihi, jhi = (int(v) for v in m.groups())
        nc = int(h.strip().split()[-1])
        if ncomp_mismatch(nc, ncomp):
            sys.exit(f"FAB ncomp {nc} != Cell_H ncomp {ncomp}")
        nx, ny = ihi - ilo + 1, jhi - jlo + 1
        # Data descriptor "(8, (64 11 52 0 1 12 0 4))" is IEEE double; assume
        # native little-endian, which is what a Mac/Linux build writes.
        data = np.frombuffer(f.read(8 * nx * ny * nc), dtype="<f8")
        arr = data.reshape((nc, ny, nx)).transpose(0, 2, 1)  # -> [c, i, j]
        return (ilo, jlo, ihi, jhi), arr


def ncomp_mismatch(a, b):
    return b is not None and a != b


def load_plane(pltdir):
    hdr = read_plotfile_header(pltdir)
    ilo, jlo, ihi, jhi = hdr["domain"]
    nx, ny = ihi - ilo + 1, jhi - jlo + 1
    ncomp, boxes, fabs = read_level0(pltdir)
    full = np.full((ncomp, nx, ny), np.nan)
    for (bilo, bjlo, bihi, bjhi), (fname, off) in zip(boxes, fabs):
        fbox, arr = read_fab(os.path.join(pltdir, "Level_0", fname), int(off), ncomp)
        a, b, c, d = fbox
        full[:, a - ilo : c - ilo + 1, b - jlo : d - jlo + 1] = arr
    if np.isnan(full).any():
        sys.exit("plane not fully covered by FABs in " + pltdir)
    return hdr, full


# --------------------------------------------------------------------------
def tanh_map(xi, plo, L, beta):
    if abs(beta) < 1e-8:
        return xi
    eta = (xi - plo) / L
    return plo + L * 0.5 * (1.0 + np.tanh(beta * (2.0 * eta - 1.0)) / np.tanh(beta))


def constant_map(xi, plo, L, fac):
    return plo + fac * (xi - plo)


def periodic_interp(xq, x, f, plo, L):
    """Linear interpolation of f(x) (uniform, periodic on [plo, plo+L)) at xq."""
    xx = np.concatenate([x - L, x, x + L])
    ff = np.concatenate([f, f, f])
    return np.interp(xq, xx, ff)


def score(a, b):
    a = a - a.mean(); b = b - b.mean()
    corr = (a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum())
    rel_rms = np.sqrt(((a - b) ** 2).mean()) / np.sqrt((b * b).mean())
    return corr, rel_rms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("uniform")
    ap.add_argument("stretched")
    ap.add_argument("--map", choices=["tanh", "constant"], default="tanh",
                    help="mapping of the stretched run in x")
    ap.add_argument("--beta", type=float, default=2.0,
                    help="TanhStretchMap beta in x, or ConstantMap scaling factor in x")
    ap.add_argument("--field", default="x_velocity")
    ap.add_argument("--plot", default=None, help="write a PNG of one row (needs matplotlib)")
    ap.add_argument("--require", type=float, default=None,
                    help="exit 1 unless mean H_phys corr >= this and beats H_xi by 0.2")
    args = ap.parse_args()

    hu, du = load_plane(args.uniform)
    hs, ds = load_plane(args.stretched)
    if hu["varnames"] != hs["varnames"]:
        sys.exit("variable lists differ between the two planes")
    if du.shape != ds.shape:
        sys.exit(f"plane shapes differ: {du.shape} vs {ds.shape}")
    ic = hu["varnames"].index(args.field)
    nx, ny = du.shape[1:]
    # stretched run: Xi grid from its own header, mapped to physical x
    plo_s = hs["problo"][0]; L_s = hs["probhi"][0] - plo_s
    xi = plo_s + (np.arange(nx) + 0.5) * (L_s / nx)
    xph = tanh_map(xi, plo_s, L_s, args.beta) if args.map == "tanh" \
        else constant_map(xi, plo_s, L_s, args.beta)
    # uniform run: its own (physical) grid, periodic in x
    plo_u = hu["problo"][0]; L_u = hu["probhi"][0] - plo_u
    xu = plo_u + (np.arange(nx) + 0.5) * (L_u / nx)

    print(f"field {args.field}: {nx} x {ny} plane, time uniform={hu['time']:.6g} stretched={hs['time']:.6g}")
    print(f"{args.map} map, param={args.beta}: x_phys - xi displacement max = {np.abs(xph - xi).max():.4g} "
          f"({np.abs(xph - xi).max() / (L_u / nx):.2f} uniform cells); physical extent "
          f"stretched [{xph.min():.4g},{xph.max():.4g}] vs uniform [{xu.min():.4g},{xu.max():.4g}]")

    cp, rp, cx, rx = [], [], [], []
    for j in range(ny):
        us = ds[ic, :, j]
        uu = du[ic, :, j]
        u_at_xph = periodic_interp(xph, xu, uu, plo_u, L_u)  # H_phys
        u_at_xi = periodic_interp(xi, xu, uu, plo_u, L_u)    # H_xi
        c1, r1 = score(us, u_at_xph); c2, r2 = score(us, u_at_xi)
        cp.append(c1); rp.append(r1); cx.append(c2); rx.append(r2)
    cp, rp, cx, rx = map(np.array, (cp, rp, cx, rx))
    print("\n  hypothesis                corr (mean, min)      rel. rms diff (mean)")
    print(f"  H_phys (correct sampling) {cp.mean():8.4f} {cp.min():8.4f}      {rp.mean():8.4f}")
    print(f"  H_xi   (pre-fix bug)      {cx.mean():8.4f} {cx.min():8.4f}      {rx.mean():8.4f}")
    thresh = args.require if args.require is not None else 0.9
    ok = cp.mean() >= thresh and (cp.mean() - cx.mean()) > 0.2
    if ok:
        print(f"\n  => PASS: mapped-run plane matches the uniform run at physical x")
    else:
        print(f"\n  => FAIL: H_phys corr {cp.mean():.4f} (need >= {thresh}) and "
              f"separation from H_xi {cp.mean() - cx.mean():.4f} (need > 0.2)")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        j = ny // 2
        fig, ax = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
        ax[0].plot(xu, du[ic, :, j], "k.-", label="uniform run vs x")
        ax[0].plot(xph, ds[ic, :, j], "r.-", label="stretched run vs x_phys")
        ax[0].set_title(f"{args.field}, row j={j}: H_phys (should overlay)")
        ax[0].legend()
        ax[1].plot(xu, du[ic, :, j], "k.-", label="uniform run vs x")
        ax[1].plot(xi, ds[ic, :, j], "b.-", label="stretched run vs xi")
        ax[1].set_title("H_xi (should NOT overlay)")
        ax[1].legend(); ax[1].set_xlabel("x")
        fig.tight_layout(); fig.savefig(args.plot, dpi=120)
        print(f"  plot written to {args.plot}")

    if args.require is not None and not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
