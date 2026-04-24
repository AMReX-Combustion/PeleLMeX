#!/usr/bin/env python3
"""
Mesh-mapping convergence analyzer.

Scans the results/ tree produced by run.sh, loads the final plotfile
for each run via yt, and reports L2 / Linf norms of the velocity
magnitude plus ratios vs the reference (fac = identity) run at each
resolution.

Usage:
    python3 analyze.py [results_root]

results_root defaults to "./results" (relative to this script).

Expected directory layout (produced by run.sh):
    results/
      incompressible/
        ref_N32/      plt_00020/   run.log
        ident_N32/    plt_00020/   run.log
        mapped_x_N32/ plt_00020/   run.log
        mapped_y_N32/ plt_00020/   run.log
        mapped_z_N32/ plt_00020/   run.log
        ref_N64/ ...
      lowmach/
        ... (same convention)

If yt cannot read a plotfile (e.g. the run crashed before it wrote
anything), that row is flagged in the table.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np


# -----------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------

# Which norms to compute.  Each key is the label, value is a function
# (arr -> scalar) where arr is a 1D flattened array of |u|^2 cell values.
# We compute |u| once and then derive norms from it.


def compute_norms(vel_magnitude: np.ndarray, cell_volume: float) -> dict:
    """
    Returns L2 (volume-averaged sqrt) and Linf norms of |u| given a
    flat array of |u| values over a uniform grid.
    """
    n = vel_magnitude.size
    return {
        "L2": float(np.sqrt(np.mean(vel_magnitude ** 2))),
        "Linf": float(np.max(np.abs(vel_magnitude))),
    }


# -----------------------------------------------------------------------
# Plotfile loading
# -----------------------------------------------------------------------


def latest_plotfile(case_dir: Path) -> Path | None:
    """
    Return the highest-numbered plt_NNNNN directory in case_dir, or None.
    """
    candidates = sorted(
        [p for p in case_dir.iterdir() if p.is_dir() and p.name.startswith("plt_")]
    )
    return candidates[-1] if candidates else None


def load_velocity(plotfile: Path) -> tuple[np.ndarray, float, float, dict]:
    """
    Load the velocity magnitude from a PeleLMeX plotfile.

    Returns (|u|_values_flattened, final_physical_time, domain_volume, meta).
    Uses yt.  Covers all (refined) cells at the finest covering grid level
    so norms are taken over the physical domain uniformly.
    """
    import yt

    ds = yt.load(str(plotfile))
    ad = ds.all_data()

    # PeleLMeX plotfiles store velocity components as 'x_velocity',
    # 'y_velocity', 'z_velocity' (AMReX derived fields map
    # these to per-axis names under yt).
    # Fall back to alternate field names if needed.
    vx = np.asarray(ad[("boxlib", "x_velocity")])
    vy = np.asarray(ad[("boxlib", "y_velocity")])
    try:
        vz = np.asarray(ad[("boxlib", "z_velocity")])
    except Exception:
        vz = np.zeros_like(vx)

    vmag = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
    t = float(ds.current_time)
    # domain volume in AMReX (Xi) space -- we want norm weighting only
    # (not the absolute value of the norm); an even weighting by cell
    # count is fine since all cells have the same Xi-volume on a single
    # level.
    vol = 1.0
    meta = {
        "ncells": int(vmag.size),
        "time": t,
        "dims": tuple(int(x) for x in ds.domain_dimensions),
        "domain_lo": tuple(float(x) for x in ds.domain_left_edge),
        "domain_hi": tuple(float(x) for x in ds.domain_right_edge),
    }
    return vmag, t, vol, meta


# -----------------------------------------------------------------------
# Presentation
# -----------------------------------------------------------------------


RES_PATTERN = re.compile(r"_N(\d+)$")


def parse_case_name(name: str) -> tuple[str, int | None]:
    """
    Parse 'ref_N32' -> ('ref', 32), 'mapped_y_N64' -> ('mapped_y', 64).
    """
    m = RES_PATTERN.search(name)
    if not m:
        return name, None
    return name[: m.start()], int(m.group(1))


def format_float(x, width=11):
    if x is None or not np.isfinite(x):
        return " " * (width - 3) + "---"
    return f"{x:>{width}.4e}"


def summarize_suite(suite_dir: Path) -> dict:
    """
    Walk a suite directory (incompressible/ or lowmach/), load each
    case's latest plotfile, and return a dict keyed by (fac_tag, N)
    holding (norms, meta, plot_path).
    """
    results = {}
    if not suite_dir.is_dir():
        return results
    for case in sorted(suite_dir.iterdir()):
        if not case.is_dir():
            continue
        tag, N = parse_case_name(case.name)
        if N is None:
            continue
        plt = latest_plotfile(case)
        if plt is None:
            results[(tag, N)] = {
                "status": "no_plotfile",
                "case_dir": case,
            }
            continue
        try:
            vmag, tfinal, vol, meta = load_velocity(plt)
        except Exception as exc:  # noqa: BLE001
            results[(tag, N)] = {
                "status": f"load_error: {exc}",
                "case_dir": case,
                "plotfile": plt,
            }
            continue
        norms = compute_norms(vmag, vol)
        results[(tag, N)] = {
            "status": "ok",
            "case_dir": case,
            "plotfile": plt,
            "norms": norms,
            "meta": meta,
            "tfinal": tfinal,
        }
    return results


def print_suite_table(name: str, results: dict) -> None:
    """
    Print a per-suite table: rows = (fac_tag, N), columns = norms + ratios.
    """
    print(f"\n=== {name} ===")

    if not results:
        print("  (no results found)")
        return

    # Gather the set of N's and tags.
    Ns = sorted({N for (_, N) in results.keys()})
    tags = sorted({tag for (tag, _) in results.keys()})
    # Put 'ref' first; preserve order otherwise.
    tags = (["ref"] if "ref" in tags else []) + [t for t in tags if t != "ref"]

    # Header
    header = f"  {'tag':<12} {'N':>4} {'t_final':>11} {'L2(|u|)':>12} {'Linf(|u|)':>12} {'L2/ref':>8} {'Linf/ref':>9} {'status'}"
    print(header)
    print("  " + "-" * (len(header) - 2))

    for N in Ns:
        ref_key = ("ref", N)
        ref_norms = results.get(ref_key, {}).get("norms")
        for tag in tags:
            r = results.get((tag, N))
            if r is None:
                continue
            if r["status"] != "ok":
                print(
                    f"  {tag:<12} {N:>4} {'':>11} {'':>12} {'':>12} {'':>8} {'':>9} "
                    f"{r['status']}"
                )
                continue
            n = r["norms"]
            L2 = n["L2"]
            Linf = n["Linf"]
            if ref_norms is not None and ref_norms["L2"] > 0:
                L2_ratio = L2 / ref_norms["L2"]
                Linf_ratio = Linf / ref_norms["Linf"] if ref_norms["Linf"] > 0 else float("nan")
            else:
                L2_ratio = float("nan")
                Linf_ratio = float("nan")
            tfinal = r["tfinal"]
            print(
                f"  {tag:<12} {N:>4}"
                f" {format_float(tfinal, 11)}"
                f" {format_float(L2, 12)}"
                f" {format_float(Linf, 12)}"
                f" {format_float(L2_ratio, 8)}"
                f" {format_float(Linf_ratio, 9)}"
                f" ok"
            )


def print_convergence_rates(name: str, results: dict) -> None:
    """
    For each tag, look at successive N values and report the apparent
    self-convergence rate:
        slope = log2( |QoI(N) - QoI(2N)| / |QoI(2N) - QoI(4N)| )
    If only two N values are present we report a simple ratio and
    report the slope as "n/a".  If three or more are present we report
    the slope for the finest available triplet.
    """
    print(f"\n=== {name}: self-convergence check ===")
    tags = sorted({tag for (tag, _) in results.keys()})
    Ns = sorted({N for (_, N) in results.keys()})
    if len(Ns) < 2:
        print("  (need at least two resolutions to report convergence)")
        return

    header = f"  {'tag':<12} {'coarse->fine':>18} {'|dQ L2|':>12} {'|dQ Linf|':>12} {'slope L2':>10} {'slope Linf':>12}"
    print(header)
    print("  " + "-" * (len(header) - 2))

    for tag in tags:
        ns_here = [N for N in Ns if (tag, N) in results and results[(tag, N)]["status"] == "ok"]
        if len(ns_here) < 2:
            continue
        for i in range(len(ns_here) - 1):
            N_c, N_f = ns_here[i], ns_here[i + 1]
            nC = results[(tag, N_c)]["norms"]
            nF = results[(tag, N_f)]["norms"]
            dL2 = nF["L2"] - nC["L2"]
            dLinf = nF["Linf"] - nC["Linf"]
            # Slope -- needs 3 levels; here we report |.| + a two-level "ratio"
            if i + 2 < len(ns_here):
                N_ff = ns_here[i + 2]
                nFF = results[(tag, N_ff)]["norms"]
                dL2_next = nFF["L2"] - nF["L2"]
                dLinf_next = nFF["Linf"] - nF["Linf"]
                with np.errstate(divide="ignore", invalid="ignore"):
                    slope_L2 = (
                        np.log2(abs(dL2) / abs(dL2_next))
                        if abs(dL2_next) > 0
                        else float("nan")
                    )
                    slope_Linf = (
                        np.log2(abs(dLinf) / abs(dLinf_next))
                        if abs(dLinf_next) > 0
                        else float("nan")
                    )
            else:
                slope_L2 = float("nan")
                slope_Linf = float("nan")
            print(
                f"  {tag:<12} {N_c:>4} -> {N_f:<10}"
                f" {format_float(abs(dL2), 12)}"
                f" {format_float(abs(dLinf), 12)}"
                f" {format_float(slope_L2, 10)}"
                f" {format_float(slope_Linf, 12)}"
            )


# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "results_root",
        nargs="?",
        default=str(Path(__file__).with_name("results")),
        help="root directory containing 'incompressible/' and 'lowmach/'",
    )
    args = parser.parse_args()
    root = Path(args.results_root).resolve()
    if not root.is_dir():
        print(f"results root not found: {root}", file=sys.stderr)
        return 2

    print(f"Reading plotfiles from: {root}")

    # Silence yt's chatty logger
    try:
        import yt

        yt.set_log_level("error")
    except Exception:
        pass

    for suite in ("incompressible", "lowmach", "lowmach_nograv"):
        suite_dir = root / suite
        if not suite_dir.is_dir():
            continue
        results = summarize_suite(suite_dir)
        print_suite_table(suite, results)
        print_convergence_rates(suite, results)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
