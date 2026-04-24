#!/usr/bin/env bash
#
# Run both incompressible and low-Mach mesh-mapping convergence sweeps.
# Accepts the same NS / MAX_STEP / STOP_TIME environment-variable overrides
# as the individual runners.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "### Incompressible (PipeFlow) ###"
"$HERE/run_incompressible.sh"

echo
echo "### Low-Mach (HotBubble, gravity ON) ###"
"$HERE/run_lowmach.sh"

echo
echo "### Low-Mach (HotBubble, gravity OFF + thermal diffusion) ###"
"$HERE/run_lowmach_nograv.sh"

echo
echo "### All runs complete.  Analyze with:"
echo "    python3 $HERE/analyze.py $HERE/results"
