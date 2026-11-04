#!/bin/bash
set -uo pipefail
cd /workspace/nucleus/work/home/sage-image
export PEAKRSS_LABEL=stageD-atomB
export PEAKRSS_LIMIT_MIB=6144
export STACK_EVERY=120
PORT=/workspace/nucleus/work/projects/p-165a6922b36b40789ddc90db34bb23c7/repo
PROV=/workspace/nucleus/work/projects/p-8c69fc7199b5456b9cfab0b695f65a6c/repo
timeout -s KILL 2260 python3 peakrss.py -- \
  timeout -s KILL 2160 ./run.sh \
    -w /work \
    -b "$PORT:/work" \
    -b /workspace/nucleus/work/home/ai-review-ci:/ci \
    -b "$PROV:/repos/sage-extensions" \
    -b "$PORT/.tmp/probe-stageD-atomB.py:/probe.py" \
    -- /bin/sh -c 'export PYTHONPATH=/ci/tool-artifacts/pytest_plugins:/work/src:/work:/repos/sage-extensions/src; exec python -u -s /probe.py' \
  > "$PORT/.tmp/stageD-atomB.log" 2>&1
