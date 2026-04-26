#!/usr/bin/env python3
"""
Analyzer for the exponential-stretch mesh-mapping sweep.

Walks the results/ tree produced by run_stretch_sweep.sh, parses each
run.log for MLMG iteration counts (MAC + nodal projections), loads
the final plotfile for per-run state, and reports:

  - Whether each (N, beta) run converged within budget.
  - Total MLMG iterations per run, separately for MAC and nodal solves.
  - Self-convergence slope  log2( ||u_N - u_{N/2}||_L2 ) at a fixed
    physical time, per beta, to characterize how the formal order of
    accuracy degrades as stretching intensifies.

Usage:
    python3 analyze_stretch.py [results_dir]

results_dir defaults to results/stretch/ alongside this script.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------
# Plotfile loader (yt).  We use covering_grid to get the velocity field
# as a structured 3D array indexed by (i, j, k), so that fine-resolution
# data can be averaged onto coarse-resolution cells for proper
# cell-aligned comparison.
# ---------------------------------------------------------------------


def load_vel_grid(path: Path):
    """
    Load (x,y,z)-velocity from an AMReX plotfile as 4-D arrays.

    Returns
    -------
    u : np.ndarray, shape (Nx, Ny, Nz, 3)
        Velocity field on a structured grid.  In 2D, Nz = 1.
    dims : tuple of int
        (Nx, Ny, Nz) domain dimensions at level 0.
    time : float
    """
    import yt  # type: ignore

    yt.set_log_level("error")
    ds = yt.load(str(path))
    dims = tuple(int(x) for x in ds.domain_dimensions)
    cg = ds.covering_grid(
        level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
    vx = np.asarray(cg["boxlib", "x_velocity"])
    vy = np.asarray(cg["boxlib", "y_velocity"])
    try:
        vz = np.asarray(cg["boxlib", "z_velocity"])
    except Exception:
        vz = np.zeros_like(vx)
    u = np.stack([vx, vy, vz], axis=-1)
    return u, dims, float(ds.current_time)


def downsample_2x(arr: np.ndarray) -> np.ndarray:
    """
    Average-pool a (Nx, Ny, Nz, C) array by factor 2 in each spatial
    direction.  Each output cell is the mean of the 8 (in 3D) or 4 (in
    2D, with Nz==1) input cells.  Requires Nx, Ny, Nz to be even (Nz
    may be 1 to indicate 2D, in which case it is preserved).
    """
    Nx, Ny, Nz, nc = arr.shape
    assert Nx % 2 == 0 and Ny % 2 == 0
    if Nz == 1:
        return arr.reshape(Nx // 2, 2, Ny // 2, 2, 1, nc).mean(axis=(1, 3))
    assert Nz % 2 == 0
    return arr.reshape(
        Nx // 2, 2, Ny // 2, 2, Nz // 2, 2, nc).mean(axis=(1, 3, 5))


def cell_diff_norms(u_coarse: np.ndarray, u_fine: np.ndarray) -> dict:
    """
    Average-pool u_fine to the resolution of u_coarse, then return the
    cell-by-cell L2 and Linf norms of |u_coarse - u_fine_pooled|, plus
    a normalized form (relative to ||u_fine||_L2).
    """
    while u_fine.shape[:3] != u_coarse.shape[:3]:
        u_fine = downsample_2x(u_fine)
    diff = u_coarse - u_fine
    diff_mag = np.sqrt(np.sum(diff ** 2, axis=-1))
    fine_mag = np.sqrt(np.sum(u_fine ** 2, axis=-1))
    L2 = float(np.sqrt(np.mean(diff_mag ** 2)))
    Linf = float(np.max(diff_mag))
    fine_L2 = float(np.sqrt(np.mean(fine_mag ** 2)))
    rel = L2 / fine_L2 if fine_L2 > 0 else float("nan")
    return {"L2": L2, "Linf": Linf, "fine_L2": fine_L2, "rel": rel}


# ---------------------------------------------------------------------
# Log parsing: sum MLMG iteration counts per solve type.  PeleLMeX logs
# each MLMG solve with a header "MLMG: MAC Projection" or "Initial..." +
# a series of "MLMG: Iteration N ..." lines followed by either
# "Final Iter. N ..." or "Failed to converge ..." depending on outcome.
# ---------------------------------------------------------------------

MLMG_HEADER_MAC = re.compile(r"MLMG: MAC Projection")
MLMG_FINAL = re.compile(r"MLMG: Final Iter\.\s+(\d+)")
MLMG_FAILED = re.compile(r"MLMG: Fail(?:ed|ing) to converge after\s+(\d+)")


def parse_log(log_path: Path) -> dict:
    """
    Walk a run.log, return a dict with:
      - mac_total_iters, mac_solve_count, mac_max_iters, any_mac_failed
      - nodal_total_iters, nodal_solve_count, nodal_max_iters,
        any_nodal_failed
      - aborted (bool: run hit amrex::Abort)
      - finalized (bool: "AMReX ... finalized" seen)

    Log layout: each MLMG solve ends with either
      "MLMG: Final Iter. N ..."   (converged)
    or
      "MLMG: Fail(ed|ing) to converge after N ..."
    MAC-projection solves are preceded by the header "MLMG: MAC
    Projection".  Any non-MAC solve is counted as nodal.
    """
    out = dict(
        mac_total_iters=0,
        mac_solve_count=0,
        mac_max_iters=0,
        any_mac_failed=False,
        nodal_total_iters=0,
        nodal_solve_count=0,
        nodal_max_iters=0,
        any_nodal_failed=False,
        aborted=False,
        finalized=False,
    )
    if not log_path.is_file():
        return out

    text = log_path.read_text(errors="replace")
    out["aborted"] = ("amrex::Abort" in text) or ("SIGABRT" in text)
    out["finalized"] = "AMReX" in text and "finalized" in text

    pending_is_mac = False
    for line in text.splitlines():
        if MLMG_HEADER_MAC.search(line):
            pending_is_mac = True
            continue
        m = MLMG_FINAL.search(line)
        if m:
            n = int(m.group(1))
            bucket = "mac" if pending_is_mac else "nodal"
            out[f"{bucket}_total_iters"] += n
            out[f"{bucket}_solve_count"] += 1
            out[f"{bucket}_max_iters"] = max(out[f"{bucket}_max_iters"], n)
            pending_is_mac = False
            continue
        m = MLMG_FAILED.search(line)
        if m:
            n = int(m.group(1))
            bucket = "mac" if pending_is_mac else "nodal"
            out[f"{bucket}_total_iters"] += n
            out[f"{bucket}_solve_count"] += 1
            out[f"{bucket}_max_iters"] = max(out[f"{bucket}_max_iters"], n)
            out[f"any_{bucket}_failed"] = True
            pending_is_mac = False
            continue

    return out


# ---------------------------------------------------------------------
# Results gathering
# ---------------------------------------------------------------------


CASE_RE = re.compile(r"stretch_N(\d+)_b([0-9.]+)")


def latest_plotfile(case_dir: Path) -> Path | None:
    candidates = sorted(
        [p for p in case_dir.iterdir() if p.is_dir() and p.name.startswith("plt_")]
    )
    return candidates[-1] if candidates else None


def gather(root: Path) -> dict:
    """
    Returns { (N, beta) : case_info }.
    """
    cases = {}
    if not root.is_dir():
        return cases
    for child in sorted(root.iterdir()):
        if not child.is_dir():
            continue
        m = CASE_RE.match(child.name)
        if not m:
            continue
        N = int(m.group(1))
        beta = float(m.group(2))
        log_path = child / "run.log"
        plt = latest_plotfile(child)
        info = parse_log(log_path)
        info["plt"] = plt
        info["case_dir"] = child
        cases[(N, beta)] = info
    return cases


# ---------------------------------------------------------------------
# Self-convergence (across N at fixed beta)
# ---------------------------------------------------------------------


def self_converge(cases: dict) -> dict:
    """
    For each beta with at least two consecutive N (Nc, 2*Nc) present and
    finalized, compute  ||u_Nc - downsample_2x(u_{2*Nc})||_L2  using
    cell-aligned 2x2x2 (or 2x2 in 2D) averaging.  Returns

      { beta : [ (Nc, {'L2':..., 'Linf':..., 'rel':...}) ... ] }

    With three or more resolutions the printer derives the observed
    order from the ratio of consecutive L2 errors.
    """
    by_beta: dict = {}
    betas = sorted({b for (_, b) in cases.keys()})
    Ns = sorted({n for (n, _) in cases.keys()})
    for beta in betas:
        series = []
        for i in range(len(Ns) - 1):
            Nc = Ns[i]
            Nf = Ns[i + 1]
            if Nf != 2 * Nc:
                continue
            ck = (Nc, beta)
            fk = (Nf, beta)
            if ck not in cases or fk not in cases:
                continue
            c = cases[ck]
            f = cases[fk]
            if not c["finalized"] or not f["finalized"]:
                continue
            if c["plt"] is None or f["plt"] is None:
                continue
            try:
                uc, dims_c, _ = load_vel_grid(c["plt"])
                uf, dims_f, _ = load_vel_grid(f["plt"])
            except Exception as exc:
                series.append((Nc, None, f"load error: {exc}"))
                continue
            # Sanity: both grids should have factor-2 ratio.
            if any(2 * dc != df for dc, df in zip(dims_c[:3], dims_f[:3])
                   if df != 1 and dc != 1):
                series.append((Nc, None,
                               f"unexpected dims Nc={dims_c} Nf={dims_f}"))
                continue
            try:
                norms = cell_diff_norms(uc, uf)
            except Exception as exc:
                series.append((Nc, None, f"diff error: {exc}"))
                continue
            series.append((Nc, norms, None))
        if series:
            by_beta[beta] = series
    return by_beta


# ---------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------


def print_iter_table(cases: dict) -> None:
    print("\n=== MLMG iteration counts per (N, beta) ===\n")
    header = (
        f"  {'beta':>6} {'N':>5} {'status':>10}"
        f"  {'mac_tot':>8} {'mac_solves':>11} {'mac_max':>8} {'mac_fail':>8}"
        f"  {'nod_tot':>8} {'nod_solves':>11} {'nod_max':>8} {'nod_fail':>8}"
    )
    print(header)
    print("  " + "-" * (len(header) - 2))
    for (N, beta) in sorted(cases.keys(), key=lambda p: (p[1], p[0])):
        c = cases[(N, beta)]
        status = "ok" if c["finalized"] else ("abort" if c["aborted"] else "?")
        print(
            f"  {beta:>6.2f} {N:>5} {status:>10}"
            f"  {c['mac_total_iters']:>8} {c['mac_solve_count']:>11}"
            f" {c['mac_max_iters']:>8}"
            f" {'yes' if c['any_mac_failed'] else 'no':>8}"
            f"  {c['nodal_total_iters']:>8} {c['nodal_solve_count']:>11}"
            f" {c['nodal_max_iters']:>8}"
            f" {'yes' if c['any_nodal_failed'] else 'no':>8}"
        )


def print_self_convergence(cases: dict) -> None:
    series = self_converge(cases)
    if not series:
        print("\n=== Self-convergence: (no adjacent N pairs to compare) ===")
        return
    print(
        "\n=== Cell-aligned self-convergence per beta ===\n"
        "  e_Nc(beta) = || u_Nc - downsample_2x(u_{2*Nc}) ||_L2,\n"
        "  where downsample_2x averages each 2x2(x2) block of the fine\n"
        "  velocity field down to coarse resolution.  Observed order is\n"
        "  log2( e_Nc / e_{2*Nc} ) when at least three consecutive N are\n"
        "  available; with two N values only the absolute error magnitude\n"
        "  and its ratio to the fine ||u||_L2 are reported.\n")
    header = (
        f"  {'beta':>6} {'Nc':>4} -> {'2Nc':<4}"
        f"  {'e_Nc L2':>12} {'e_Nc Linf':>12}"
        f"  {'rel(L2)':>10} {'order':>8}")
    print(header)
    print("  " + "-" * (len(header) - 2))
    for beta, rows in sorted(series.items()):
        prev_L2 = None
        for (Nc, norms, err) in rows:
            tag = f"  {beta:>6.2f} {Nc:>4} -> {2*Nc:<4}"
            if err:
                print(f"{tag}  {err}")
                prev_L2 = None
                continue
            order_str = ""
            if prev_L2 is not None and norms["L2"] > 0 and prev_L2 > 0:
                order_str = f"{np.log2(prev_L2 / norms['L2']):.3f}"
            print(
                f"{tag}"
                f"  {norms['L2']:>12.4e} {norms['Linf']:>12.4e}"
                f"  {norms['rel']:>10.4e} {order_str:>8}")
            prev_L2 = norms["L2"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "results_dir",
        nargs="?",
        default=str(Path(__file__).with_name("results") / "stretch"),
    )
    args = ap.parse_args()
    root = Path(args.results_dir).resolve()
    if not root.is_dir():
        print(f"results dir not found: {root}", file=sys.stderr)
        return 2
    print(f"Reading results from: {root}")
    cases = gather(root)
    if not cases:
        print("no cases found")
        return 1
    print_iter_table(cases)
    print_self_convergence(cases)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
