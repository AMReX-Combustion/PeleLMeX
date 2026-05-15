#!/usr/bin/env python3
"""Analyse tempState CSV from the 3D BFS recycling run.

We call the flow "statistically stationary" at the first time at which the
running mean of kinetic energy over a window of W flow-throughs agrees with
the running mean over the *next* W flow-throughs to within a relative
tolerance of TOL. That is: the process is stationary once the forward-looking
and backward-looking KE averages have converged.

Reports stationarity time in absolute seconds, in flow-through units, and in
time steps. If no stationarity is reached in the available data, reports the
trend of the last window instead.
"""

import csv
import sys
from pathlib import Path

# Flow-through time for this case:
#   L_x = 0.08 m, U_bulk = 10 m/s  ->  tau = 8 ms
TAU = 0.008
# Stationarity threshold. Instantaneous domain-integrated KE in this box
# has a natural fluctuation envelope of ~30% (only a handful of
# independent large eddies fit in the post-step region), so comparing
# short windows will always look "non-stationary". We smooth over
# WINDOW_TAU flow-throughs and accept TOL agreement between adjacent
# windows -- loose enough not to trip on eddy-scale statistics but tight
# enough that a true transient will not cross it.
TOL = 0.10           # 10% relative tolerance
WINDOW_TAU = 4.0     # 4 flow-through averaging window


def read_tempstate(path: Path):
    with path.open() as f:
        rows = list(csv.DictReader(f))
    return [
        (int(r["iter"]), float(r["time"]), float(r["dt"]), float(r["kinEnergy"]),
         float(r["enstrophy"]))
        for r in rows
    ]


def running_mean(values, times, t0, t1):
    """Time-weighted mean of `values` over [t0, t1]."""
    num = 0.0
    den = 0.0
    for i in range(1, len(times)):
        dt_i = times[i] - times[i - 1]
        t_mid = 0.5 * (times[i] + times[i - 1])
        if t0 <= t_mid <= t1:
            num += values[i] * dt_i
            den += dt_i
    return num / den if den > 0 else None


def main(path: Path):
    rows = read_tempstate(path)
    if not rows:
        print("No data in tempState yet.")
        return

    steps = [r[0] for r in rows]
    times = [r[1] for r in rows]
    ke = [r[3] for r in rows]
    ens = [r[4] for r in rows]

    t_end = times[-1]
    print(f"Run progress: step {steps[-1]}, t = {t_end*1e3:.3f} ms  "
          f"({t_end/TAU:.2f} tau), {len(rows)} tempState rows.")
    print(f"KE: current {ke[-1]:.6g}, enstrophy: current {ens[-1]:.6g}")

    window = WINDOW_TAU * TAU
    if t_end < 2 * window:
        print(f"Need at least {2*window*1e3:.1f} ms "
              f"({2*WINDOW_TAU:.1f} tau) of data; have {t_end*1e3:.2f} ms.")
        return

    # Step through candidate stationarity times in 0.5-tau increments; flag
    # the first t* at which the trailing and leading W-tau means agree.
    step_dt = 0.5 * TAU
    t_candidate = window
    stationary_t = None
    while t_candidate + window <= t_end:
        trailing = running_mean(ke, times, t_candidate - window, t_candidate)
        leading = running_mean(ke, times, t_candidate, t_candidate + window)
        if trailing is None or leading is None:
            t_candidate += step_dt
            continue
        rel = abs(leading - trailing) / max(abs(trailing), 1e-30)
        if rel < TOL and stationary_t is None:
            stationary_t = t_candidate
        t_candidate += step_dt

    if stationary_t is not None:
        # Find the step index closest to stationary_t
        closest_step = min(rows, key=lambda r: abs(r[1] - stationary_t))[0]
        print()
        print(f"Stationary onset: t* = {stationary_t*1e3:.2f} ms  "
              f"= {stationary_t/TAU:.2f} tau  (step ~{closest_step})")
        print(f"(trailing vs leading KE means within {TOL*100:.1f}%, "
              f"window = {WINDOW_TAU:.1f} tau)")

        # Robustness check: sweep across all candidate t* >= onset and
        # report the fraction that pass the same tolerance. If the fit is
        # real, most/all later candidates should also pass.
        passes = 0
        checks = 0
        t_c = stationary_t
        while t_c + window <= t_end:
            trl = running_mean(ke, times, t_c - window, t_c)
            lead = running_mean(ke, times, t_c, t_c + window)
            if trl is not None and lead is not None:
                if abs(lead - trl) / max(abs(trl), 1e-30) < TOL:
                    passes += 1
                checks += 1
            t_c += step_dt
        if checks:
            print(f"Post-onset robustness: {passes}/{checks} "
                  f"window positions pass the same test.")

        # Statistics over the claimed stationary regime (onset -> end)
        def window_stats(values, times, t0, t1):
            num = 0.0
            den = 0.0
            num2 = 0.0
            for i in range(1, len(times)):
                dt_i = times[i] - times[i - 1]
                t_mid = 0.5 * (times[i] + times[i - 1])
                if t0 <= t_mid <= t1:
                    num += values[i] * dt_i
                    num2 += values[i] * values[i] * dt_i
                    den += dt_i
            if den <= 0:
                return None, None
            mean = num / den
            var = num2 / den - mean * mean
            return mean, max(var, 0.0) ** 0.5

        ke_m, ke_s = window_stats(ke, times, stationary_t, t_end)
        ens_m, ens_s = window_stats(ens, times, stationary_t, t_end)
        print()
        print(f"Stationary regime statistics ({stationary_t/TAU:.1f} -> "
              f"{t_end/TAU:.1f} tau):")
        if ke_m is not None:
            print(f"  KE:        mean {ke_m:.3e}, std {ke_s:.3e}  "
                  f"(CoV {100*ke_s/abs(ke_m):.1f}%)")
        if ens_m is not None:
            print(f"  Enstrophy: mean {ens_m:.3e}, std {ens_s:.3e}  "
                  f"(CoV {100*ens_s/abs(ens_m):.1f}%)")
    else:
        # Still developing -- report most recent window trend
        last_trail = running_mean(ke, times, t_end - 2 * window, t_end - window)
        last_lead = running_mean(ke, times, t_end - window, t_end)
        if last_trail is not None and last_lead is not None:
            rel = (last_lead - last_trail) / max(abs(last_trail), 1e-30)
            print()
            print(f"Not yet stationary. Most recent KE shift: {rel*100:+.2f}% "
                  f"over the last {WINDOW_TAU:.1f}-tau window.")


if __name__ == "__main__":
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
        "run01/temporals/tempState")
    main(path)
