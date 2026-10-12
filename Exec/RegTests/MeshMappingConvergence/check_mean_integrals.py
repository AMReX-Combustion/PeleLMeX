#!/usr/bin/env python3
"""Check the physical-volume means produced by run_mean_integrals.sh.

channel : the periodic-channel background pressure gradient printed by
          PeleLMeX ("Adding background pressure gradient in dir 0: X") on
          the uniform and the TanhStretchMap grids must agree with each
          other and with the exact value -16 mu U / Ly^2 of the quartic
          profile, to discretisation error.  Before the MFSum fix the
          stretched value was a Xi-weighted mean, off by O(1).
chamber : the closed-chamber dp0/dt printed at peleLM.v >= 3 must agree
          between the two grids, and the projection compatibility residual
          printed at peleLM.v >= 2 must be at round-off on both.

Exit status 1 on any failure; the per-case verdicts are printed.
"""
import re
import sys
from pathlib import Path

results = Path(sys.argv[1] if len(sys.argv) > 1 else "results/mean_integrals")
# The cell-centred Dirichlet wall treatment evaluates the wall gradient from
# the first cell only, which for the quartic profile is a first-order error
# of ~0.75 h_wall/(Ly/2): ~5% at N=32 on both grids (the stretched one has
# a smaller wall cell but its metric is sampled at the wall face, which
# biases the other way).  A Xi-weighted mean on the beta=2 tanh grid is off
# by ~100%, so 10% separates fixed from broken with margin, and the N=64
# stretched run must have a smaller error than the N=32 one.
TOL_CHANNEL = 1.0e-1   # relative to exact
TOL_DP0DT = 2.0e-2     # relative, uniform vs tanh
TOL_COMPAT = 1.0e-9    # relative compatibility residual

RE_GP0 = re.compile(r"Adding background pressure gradient in dir 0: (\S+)")
RE_DP0 = re.compile(r"Closed chamber pOld: (\S+), pNew: (\S+), dp0dt: (\S+)")
RE_COMPAT = re.compile(r"\.\.\. (\S+ projection) compatibility: .* = (\S+) \(")

ok = True


def log(tag):
    p = results / tag / "run.log"
    if not p.exists():
        print(f"MISSING {p}")
        return None
    return p.read_text()


def first(rx, text, grp=1):
    m = rx.search(text)
    return float(m.group(grp)) if m else None


def verdict(name, val, ref, tol):
    global ok
    if val is None or ref is None:
        print(f"  {name:38s} MISSING")
        ok = False
        return
    rel = abs(val - ref) / max(abs(ref), 1e-300)
    good = rel <= tol
    ok &= good
    print(f"  {name:38s} {val:+.6e} vs {ref:+.6e}  rel={rel:.3e}  "
          f"{'OK' if good else 'FAIL'}")


# ------------------------------------------------------------- channel
pf = results / "channel.params"
if pf.exists():
    prm = dict(kv.split("=") for kv in pf.read_text().split()[1:])
    mu, U, Ly = (float(prm[k]) for k in ("mu", "U", "Ly"))
    exact = -16.0 * mu * U / Ly**2
    print("channel: periodic-channel background pressure gradient (first step)")
    lu, lt = log("channel_uniform"), log("channel_tanh")
    gu = first(RE_GP0, lu) if lu else None
    gt = first(RE_GP0, lt) if lt else None
    l64 = log("channel_tanh_N64")
    g64 = first(RE_GP0, l64) if l64 else None
    verdict("uniform  N32 vs exact -16 mu U/Ly^2", gu, exact, TOL_CHANNEL)
    verdict("tanh     N32 vs exact -16 mu U/Ly^2", gt, exact, TOL_CHANNEL)
    verdict("tanh     N64 vs exact -16 mu U/Ly^2", g64, exact, TOL_CHANNEL)
    if None not in (gt, g64):
        e32 = abs(gt - exact)
        e64 = abs(g64 - exact)
        good = e64 < e32
        ok &= good
        print(f"  {'tanh error N64 < N32 (converging)':38s} {e64:.3e} < {e32:.3e}  "
              f"{'OK' if good else 'FAIL'}")

# ------------------------------------------------------------- chamber
lu, lt = log("chamber_uniform"), log("chamber_tanh")
if lu is not None or lt is not None:
    print("chamber: closed-chamber dp0/dt and projection compatibility")
    du = first(RE_DP0, lu, 3) if lu else None
    dt = first(RE_DP0, lt, 3) if lt else None
    verdict("tanh dp0dt vs uniform (first step)", dt, du, TOL_DP0DT)
    for tag, text in (("uniform", lu), ("tanh", lt)):
        if text is None:
            ok = False
            continue
        worst = {}
        for m in RE_COMPAT.finditer(text):
            who, rel = m.group(1), float(m.group(2))
            worst[who] = max(worst.get(who, 0.0), rel)
        if not worst:
            print(f"  {tag}: no compatibility lines found (peleLM.v >= 2?)  FAIL")
            ok = False
        for who, rel in sorted(worst.items()):
            good = rel <= TOL_COMPAT
            ok &= good
            print(f"  {tag:8s} {who:22s} max residual {rel:.3e}  "
                  f"{'OK' if good else 'FAIL'}")
        if "incompatible" in text:
            print(f"  {tag}: run aborted on the compatibility gate  FAIL")
            ok = False

print("ALL OK" if ok else "FAILURES")
sys.exit(0 if ok else 1)
