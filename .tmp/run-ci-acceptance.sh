#!/usr/bin/env bash
# The body of `just ci-acceptance <family> <shard> <shards> [selection]` from the project justfile,
# executed directly because `just` is not installed in this environment.
set -euo pipefail
family=$1; shard=$2; shards=$3; selection=${4:-}
out=.tmp/ci/results/$family-$shard
mkdir -p "$out"
git rev-parse HEAD > "$out/commit"
provider_dir=.tmp/ci/sage-extensions
if [ ! -d "$provider_dir/.git" ]; then
    echo "ERROR: no provider checkout at $provider_dir; run 'just ci-provision-sage' first" >&2
    exit 1
fi
provider_head="$(git -C "$provider_dir" rev-parse HEAD)"
pinned_rev="$(sed -n 's|.*sage-extensions @ git+https://github.com/dzackgarza/sage-extensions@\([0-9a-f]\{40\}\).*|\1|p' pyproject.toml | head -1)"
if [ "$provider_head" != "$pinned_rev" ]; then
    echo "ERROR: provider checkout $provider_head does not match the pin $pinned_rev" >&2
    exit 1
fi
printf '%s\n' "$provider_head" > "$out/research_commit"
SAGE_BIN=${SAGE_BIN:-/sage/.venv/bin/sage}
sage_python="$("$SAGE_BIN" -c 'import sys; print(sys.executable)')"
# `ci-provision-sage` installs the pinned provider into the image's interpreter. That writes into the
# shared runtime, which is not this project's to modify, so the pinned checkout is placed first on the
# path instead: identical resolution of `sage_extensions`, resolved to the pyproject pin verified above
# rather than to whatever revision the runtime happens to carry.
export PYTHONPATH="$PWD/.tmp/ci/sage-extensions/src:$PWD/.tmp/ci/ai-review-ci/tool-artifacts/pytest_plugins:$PWD/src:$PWD"
echo "[gate] family=$family shard=$shard/$shards selection='$selection'"
echo "[gate] port $(git rev-parse --short HEAD) | provider $provider_head | sage_python $sage_python"
if [[ -n "$selection" ]]; then
    tests=("tests/acceptance/test_$family.sage")
    shard_args=()
elif [[ "$family" == reflective ]]; then
    timeout --signal=INT 345m "$sage_python" ci/pytest_batches.py --python "$sage_python" --out "$out" --batch-size 10 -- "tests/acceptance/test_reflective.sage" --shard-id="$shard" --num-shards="$shards" 2>&1 | tee "$out/log"
    exit ${PIPESTATUS[0]}
else
    tests=("tests/acceptance/test_$family.sage")
    shard_args=(--shard-id="$shard" --num-shards="$shards")
fi
# Deviation from the hosted recipe, stated rather than silent: `pytest-instafail` and `pytest-shard`
# are installed by `just ci-provision-sage` in CI and are absent here, and installing into the shared
# runtime is not this project's to do. Both are reporting/partitioning only: --instafail changes when
# failures are printed, and --shard-id/--num-shards with a single shard is the identity partition. Case
# selection, assertions and the JUnit record are unchanged. The reflective batch runner (which needs
# neither) is used verbatim when that family is requested.
echo "[gate] NOTE: --instafail/--shard-id omitted locally (CI-provisioned packages absent); single-shard identity"
timeout --signal=INT 345m "$sage_python" -c 'import sage.all, sys, pytest; sys.exit(pytest.main(sys.argv[1:]))' \
    -p qc_sage_session "${tests[@]}" --runxfail -rA --tb=long \
    --junitxml="$out/junit.xml" 2>&1 | tee "$out/log"
exit "${PIPESTATUS[0]}"
