#!/bin/bash
set -uo pipefail
cd /workspace/nucleus/work/home/sage-image
export PEAKRSS_LABEL=stageD-rate
export PEAKRSS_LIMIT_MIB=6144
export STACK_EVERY=120
export RATE_BLOCK=300
PORT=/workspace/nucleus/work/projects/p-165a6922b36b40789ddc90db34bb23c7/repo
PROV=/workspace/nucleus/work/projects/p-8c69fc7199b5456b9cfab0b695f65a6c/repo
timeout -s KILL 1560 ./run.sh \
    -w /work \
    -b "$PORT:/work" \
    -b "$PROV:/repos/sage-extensions" \
    -b "$PORT/.tmp/probe-stageD-conversion-rate.py:/probe.py" \
    -- /bin/sh -c 'export PYTHONPATH=/ci/tool-artifacts/pytest_plugins:/work/src:/work:/repos/sage-extensions/src; exec python -u -s /probe.py' \
  > "$PORT/.tmp/stageD-rate.log" 2>&1
