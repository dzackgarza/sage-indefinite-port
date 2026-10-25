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
    # Sage is a meson-python editable install that rebuilds on import; print
    # that rebuild, and fail here if Sage does not import.
    export MESONPY_EDITABLE_VERBOSE=1
    "$SAGE_BIN" -python -c 'import sage.all'
    "$SAGE_BIN" -pip install --quiet "pytest>=9.1,<10" pyyaml "hypothesis>=6.165,<7" pytest-shard
    # The acceptance contract is against this exact research API revision.
    # Tracking main is not reproducible and currently imports a Sage category
    # absent from the pinned CI Sage image (sage.categories.ordered_fields).
    # Keep the clone shallow and fetch the exact provider commit explicitly.
    mkdir -p {{ci_dir}}
    git clone --quiet --depth 1 https://github.com/dzackgarza/research {{ci_dir}}/research
    research_rev="f44e76f18e2e9cec43b09a846e2e74bca5fd50f6"
    git -C {{ci_dir}}/research fetch --quiet --depth 1 origin "$research_rev"
    git -C {{ci_dir}}/research -c advice.detachedHead=false checkout --quiet "$research_rev"
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
    # The preamble's Julia bridge (OSCAR comparisons) runs this Julia, the
    # version the local juliaup provides.
    mkdir -p {{ci_dir}}/julia
    # The image has no curl; Sage's Python fetches the tarball.
    "$sage_python" -c 'import sys, urllib.request; urllib.request.urlretrieve(sys.argv[1], sys.argv[2])' \
        https://julialang-s3.julialang.org/bin/linux/x64/1.12/julia-1.12.5-linux-x86_64.tar.gz {{ci_dir}}/julia.tar.gz
    tar -xzf {{ci_dir}}/julia.tar.gz -C {{ci_dir}}/julia --strip-components=1

# CI only: the bridge's Julia environment, instantiated and precompiled into the default depot.
ci-provision-julia:
    #!/usr/bin/env bash
    set -euo pipefail
    # sage_julia_bridge resolves its Julia at import.
    export SAGE_JULIA_COMMAND="$PWD/{{ci_dir}}/julia/bin/julia"
    sage_python="$("$SAGE_BIN" -python -c 'import sys; print(sys.executable)')"
    env_dir="$("$sage_python" -c 'import pathlib, sage_julia_bridge; print(pathlib.Path(sage_julia_bridge.__file__).parent / "julia_env")')"
    {{ci_dir}}/julia/bin/julia --project="$env_dir" -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# CI only: the cache key of the bridge's Julia environment (its Project; the
# installed package ships no Manifest, so the julia job resolves one and the
# cache carries it).
ci-julia-key:
    #!/usr/bin/env bash
    set -euo pipefail
    # sage_julia_bridge resolves its Julia at import.
    export SAGE_JULIA_COMMAND="$PWD/{{ci_dir}}/julia/bin/julia"
    sage_python="$("$SAGE_BIN" -python -c 'import sys; print(sys.executable)')"
    "$sage_python" -c 'import hashlib, pathlib, sage_julia_bridge; print("julia-depot-" + hashlib.sha256((pathlib.Path(sage_julia_bridge.__file__).parent / "julia_env" / "Project.toml").read_bytes()).hexdigest()[:16])'

# CI only: run one shard of one frozen acceptance family; results land in .tmp/ci/results/<family>-<shard>.
ci-acceptance family shard shards:
    #!/usr/bin/env bash
    set -euo pipefail
    out={{ci_dir}}/results/{{family}}-{{shard}}
    mkdir -p "$out"
    git rev-parse HEAD > "$out/commit"
    git -C {{ci_dir}}/research rev-parse HEAD > "$out/research_commit"
    export PYTHONPATH="$PWD/{{ci_dir}}/ai-review-ci/tool-artifacts/pytest_plugins:$PWD/src:$PWD:$PWD/{{ci_dir}}/vinal/src/sage"
    export SAGE_JULIA_COMMAND="$PWD/{{ci_dir}}/julia/bin/julia"
    # SIGINT before the job's hard limit, so pytest still writes the JUnit
    # report for every case that finished.
    timeout --signal=INT 345m "$SAGE_BIN" -python -c 'import sage.all, sys, pytest; sys.exit(pytest.main(sys.argv[1:]))' \
        -p qc_sage_session tests/acceptance/test_{{family}}.sage --runxfail -rA --tb=long \
        --shard-id={{shard}} --num-shards={{shards}} --junitxml="$out/junit.xml" 2>&1 | tee "$out/log"
