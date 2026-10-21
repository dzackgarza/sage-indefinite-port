# ai-review-ci SageMath QC delegation justfile.
# The central implementation lives in ~/ai-review-ci/justfiles/sage.just.
# Public recipes delegate to that central justfile while preserving this repo as the caller root.

# ai-review-ci contract variables consumed by doctor and workflow installers.
ai_review_ci_schema_version := "1"
ai_review_ci_profile := "sage"
ai_review_ci_ref := "main"
ai_review_ci_release_channel := "main"
ai_review_ci_workflow_template_version := "1"
ai_review_ci_local_delegation := "global-justfile"
ai_review_ci_default_branch := "main"
# List available recipes.
default:
    @just --list

# Run ambient SageMath QC through the central implementation.
ambient:
    @just -f ~/ai-review-ci/justfiles/sage.just -d . ambient

# Run commit-tier SageMath QC through the central implementation.
test-commit:
    @just -f ~/ai-review-ci/justfiles/sage.just -d . test-commit

# Run the full SageMath test suite before pushing.
test-push:
    @just -f ~/ai-review-ci/justfiles/sage.just -d . test-push

# Run CI acceptance QC through the central implementation.
test-ci:
    @just -f ~/ai-review-ci/justfiles/sage.just -d . test-ci

# The acceptance workflow's jobs run inside the Sage image named in
# .github/workflows/acceptance.yml; these two recipes are its only steps.
ci_dir := ".tmp/ci"

# CI only: install the acceptance dependencies into the Sage image's interpreter.
ci-provision-sage:
    #!/usr/bin/env bash
    set -euo pipefail
    "$SAGE_BIN" -pip install --quiet "pytest>=9.1,<10" pyyaml "hypothesis>=6.165,<7" pytest-shard
    # dzack-research from its main branch, cloned without submodules: pip's
    # git install recurses into them, and one is an SSH-only workspace.
    mkdir -p {{ci_dir}}
    git clone --quiet --depth 1 https://github.com/dzackgarza/research {{ci_dir}}/research
    # uv, not pip: research declares tree-sitter-sage and sage-stubs only in
    # [tool.uv.sources], which pip ignores.
    sage_python="$("$SAGE_BIN" -python -c 'import sys; print(sys.executable)')"
    uv pip install --quiet --python "$sage_python" -e {{ci_dir}}/research
    uv pip install --quiet --python "$sage_python" --no-deps -e .
    # VinAl's Sage module is src/sage/vinal.py in its checkout (see .envrc), at
    # the commit uv.lock pins.
    vinal_rev="$(sed -n 's|^source = { git = "https://github.com/dzackgarza/vinal#\([0-9a-f]*\)" }$|\1|p' uv.lock)"
    test -n "$vinal_rev"
    git clone --quiet https://github.com/dzackgarza/vinal {{ci_dir}}/vinal
    git -C {{ci_dir}}/vinal -c advice.detachedHead=false checkout --quiet "$vinal_rev"
    git clone --quiet --depth 1 https://github.com/dzackgarza/ai-review-ci {{ci_dir}}/ai-review-ci

# CI only: run one shard of one frozen acceptance family; results land in .tmp/ci/results/<family>-<shard>.
ci-acceptance family shard shards:
    #!/usr/bin/env bash
    set -euo pipefail
    out={{ci_dir}}/results/{{family}}-{{shard}}
    mkdir -p "$out"
    git rev-parse HEAD > "$out/commit"
    git -C {{ci_dir}}/research rev-parse HEAD > "$out/research_commit"
    export PYTHONPATH="$PWD/{{ci_dir}}/ai-review-ci/tool-artifacts/pytest_plugins:$PWD/src:$PWD:$PWD/{{ci_dir}}/vinal/src/sage"
    # SIGINT before the job's hard limit, so pytest still writes the JUnit
    # report for every case that finished.
    timeout --signal=INT 345m "$SAGE_BIN" -python -c 'import sage.all, sys, pytest; sys.exit(pytest.main(sys.argv[1:]))' \
        -p qc_sage_session tests/acceptance/test_{{family}}.sage --runxfail -rA --tb=long \
        --shard-id={{shard}} --num-shards={{shards}} --junitxml="$out/junit.xml" 2>&1 | tee "$out/log"
