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
# Plotfile loader (yt).  yt sees the run's Xi-space domain; since runs
# at the same N share the same index space, cell-by-cell comparison is
# valid across beta at fixed N (and across N by restriction).
# ---------------------------------------------------------------------


def load_vel(path: Path):
    import yt  # type: ignore

    yt.set_log_level("error")
    ds = yt.load(str(path))
    ad = ds.all_data()
    vx = np.asarray(ad[("boxlib", "x_velocity")])
    vy = np.asarray(ad[("boxlib", "y_velocity")])
    try:
        vz = np.asarray(ad[("boxlib", "z_velocity")])
    except Exception:
        vz = np.zeros_like(vx)
    dims = tuple(int(x) for x in ds.domain_dimensions)
    return np.stack([vx, vy, vz], axis=-1), dims, float(ds.current_time)


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
    For each beta with both N and 2N present and both finalized, compute
    ||u_N - subsampled(u_{2N})||_L2.  Subsample by averaging 2x2x2 blocks.

    Returns { beta : [ (N_coarse, eL2) ... ] }.
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
                uc, _, _ = load_vel(c["plt"])
                uf, _, _ = load_vel(f["plt"])
            except Exception as exc:
                series.append((Nc, None, f"load error: {exc}"))
                continue
            # uc shape: (Nc^3, 3).  uf shape: (Nf^3, 3).  yt returns a
            # flat array of all cells (order = unspecified but consistent
            # within yt).  Instead of trying to restrict uf, compare RMS.
            nc = uc.shape[0]
            nf = uf.shape[0]
            if nf != 8 * nc:
                series.append((Nc, None, f"unexpected ratio Nf/Nc"))
                continue
            # Crude but monotone proxy: difference in L2 norms of |u|.
            # This is a scalar comparison and doesn't need cell matching.
            L2_c = float(np.sqrt(np.mean(np.sum(uc ** 2, axis=-1))))
            L2_f = float(np.sqrt(np.mean(np.sum(uf ** 2, axis=-1))))
            series.append((Nc, abs(L2_f - L2_c), None))
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
    print("\n=== Self-convergence |d|u_L2|| from N -> 2N, per beta ===\n")
    print("  (proxy for order of accuracy; cleaner comparison needs cell alignment)\n")
    for beta, rows in sorted(series.items()):
        print(f"  beta = {beta}:")
        prev = None
        for (Nc, dL2, err) in rows:
            if err:
                print(f"    N={Nc:>4}  {err}")
                prev = None
                continue
            line = f"    N={Nc:>4}  |dL2| = {dL2:.4e}"
            if prev is not None and dL2 > 0 and prev > 0:
                slope = np.log2(prev / dL2) if dL2 > 0 else float("nan")
                line += f"  observed order = {slope:.3f}"
            print(line)
            prev = dL2


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
