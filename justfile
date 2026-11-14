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

# Run the README example as a doctest; it is the smallest acceptance.
readme-example:
    #!/usr/bin/env bash
    set -euo pipefail
    if [ -n "${SAGE_BIN:-}" ]; then
        sage_python="$("$SAGE_BIN" -c 'import sys; print(sys.executable)')"
    else
        sage_python="$(just -f ~/ai-review-ci/justfiles/sage.just -d . _sage-python)"
    fi
    PYTHONPATH=src "$sage_python" -m doctest README.md

# The acceptance workflow's jobs run inside the Sage image named in
# .github/workflows/acceptance.yml; these recipes are its only steps.
ci_dir := ".tmp/ci"

# CI only: install the acceptance dependencies into the Sage image's interpreter.
ci-provision-sage:
    #!/usr/bin/env bash
    set -euo pipefail
    # The image's MAXIMA_PREFIX names the version root, but its installed
    # Maxima packages live one level below that root in share/.  maxima_lib
    # appends package names directly to MAXIMA_PREFIX, so correct the path
    # before Sage initializes ECL and preserve it for later workflow steps.
    if [[ -n "${MAXIMA_PREFIX:-}" && -d "${MAXIMA_PREFIX}/share/to_poly_solve" && ! -d "${MAXIMA_PREFIX}/to_poly_solve" ]]; then
        export MAXIMA_PREFIX="${MAXIMA_PREFIX}/share"
    fi
    if [[ -n "${MAXIMA_PREFIX:-}" && -n "${GITHUB_ENV:-}" ]]; then
        printf 'MAXIMA_PREFIX=%s\n' "$MAXIMA_PREFIX" >> "$GITHUB_ENV"
    fi
    # Sage is a meson-python editable install that rebuilds on import; print
    # that rebuild, and fail here if Sage does not import.
    export MESONPY_EDITABLE_VERBOSE=1
    "$SAGE_BIN" -c 'import sage.all; from sage.features.bliss import Bliss; Bliss().require(); from sage.interfaces.maxima_lib import maxima_lib; assert str(maxima_lib.eval("1+1")) == "2"'
    sage_python="$("$SAGE_BIN" -c 'import sys; print(sys.executable)')"
    uv pip install --quiet --python "$sage_python" "pytest>=9.1,<10" pyyaml "hypothesis>=6.165,<7" pytest-shard pytest-instafail "pytest-timeout>=2.4" "pytest-reportlog>=1.0" PyNormaliz "pycddlib>=3.0.2"
    "$sage_python" -c 'import cdd.gmp'
    # The acceptance contract is against the shared-library revision declared in
    # pyproject.toml; read it from there so the packaging pin and the CI pin
    # cannot drift apart. The library is sage-extensions: the pre-move
    # dzack_research.preamble.categories modules no longer exist.
    # Keep the clone shallow and fetch the exact provider commit explicitly.
    mkdir -p {{ci_dir}}
    git clone --quiet --depth 1 https://github.com/dzackgarza/sage-extensions {{ci_dir}}/sage-extensions
    provider_rev="$(sed -n 's|.*sage-extensions @ git+https://github.com/dzackgarza/sage-extensions@\([0-9a-f]\{40\}\).*|\1|p' pyproject.toml | head -1)"
    if [ -z "$provider_rev" ]; then
        echo "ERROR: pyproject.toml declares no sage-extensions revision to provision" >&2
        exit 1
    fi
    git -C {{ci_dir}}/sage-extensions fetch --quiet --depth 1 origin "$provider_rev"
    git -C {{ci_dir}}/sage-extensions -c advice.detachedHead=false checkout --quiet "$provider_rev"
    # uv, not pip: sage-extensions declares tree-sitter-sage and sage-stubs only
    # in [tool.uv.sources], which pip ignores.
    sage_python="$("$SAGE_BIN" -c 'import sys; print(sys.executable)')"
    # The image owns sage-julia-bridge and the other Sage projects, except that
    # the bridge itself is not installed in this image. Install it separately:
    # its only dependency is pydantic, and it does not re-resolve the image's
    # conflicting tree-sitter-sage sources.
    uv pip install --quiet --python "$sage_python" "sage-julia-bridge @ git+https://github.com/dzackgarza/sage-julia-bridge"
    # Installing the library with dependency resolution re-resolves the image's
    # already-installed editable projects and rejects their different
    # tree-sitter-sage source URLs. Install the exact provider checkout without
    # replacing the image's dependency graph.
    uv pip install --quiet --python "$sage_python" --no-deps -e {{ci_dir}}/sage-extensions
    uv pip install --quiet --python "$sage_python" --no-deps -e .
    git clone --quiet --depth 1 https://github.com/dzackgarza/ai-review-ci {{ci_dir}}/ai-review-ci

# CI only: the bridge's Julia environment, instantiated and precompiled into the default depot.
ci-provision-julia:
    #!/usr/bin/env bash
    set -euo pipefail
    sage_python="$("$SAGE_BIN" -c 'import sys; print(sys.executable)')"
    env_dir="$("$sage_python" -c 'import pathlib, sage_julia_bridge; print(pathlib.Path(sage_julia_bridge.__file__).parent / "julia_env")')"
    # The worker sees both its JSON project and the shared Oscar project.
    # Instantiate each selected manifest before exercising the same imports.
    JULIA_PKG_PRECOMPILE_AUTO=0 julia --project="$env_dir" -e 'using Pkg; bridge_project = Base.active_project(); Pkg.activate("v$(VERSION.major).$(VERSION.minor)"; shared=true); Pkg.instantiate(); Pkg.activate(dirname(bridge_project)); Pkg.instantiate(); using Oscar, JSON'

# CI only: run one long true witness from witnesses/ against the pinned provider.
ci-witness file test:
    #!/usr/bin/env bash
    set -euo pipefail
    out={{ci_dir}}/results/{{test}}
    mkdir -p "$out"
    git rev-parse HEAD > "$out/commit"
    git -C {{ci_dir}}/sage-extensions rev-parse HEAD > "$out/provider_commit"
    sage_python="$("$SAGE_BIN" -c 'import sys; print(sys.executable)')"
    export PYTHONPATH="$PWD/{{ci_dir}}/ai-review-ci/tool-artifacts/pytest_plugins:$PWD/src:$PWD"
    "$sage_python" -c 'import sage.all, sys, pytest; sys.exit(pytest.main(sys.argv[1:]))' \
        -p qc_sage_session "witnesses/{{file}}::{{test}}" -vv -rA --tb=long --instafail \
        -p ci.test_budgets --timeout-method=signal --report-log="$out/report.jsonl" --junitxml="$out/junit.xml" 2>&1 | tee "$out/log"
    exit "${PIPESTATUS[0]}"

# CI only: key Julia packages by the bridge project and shared image lock.
ci-julia-key:
    #!/usr/bin/env bash
    set -euo pipefail
    sage_python="$("$SAGE_BIN" -c 'import sys; print(sys.executable)')"
    env_dir="$("$sage_python" -c 'import pathlib, sage_julia_bridge; print(pathlib.Path(sage_julia_bridge.__file__).parent / "julia_env")')"
    julia --startup-file=no --project="$env_dir" -e 'using Pkg, SHA; bridge_project = Base.active_project(); Pkg.activate("v$(VERSION.major).$(VERSION.minor)"; shared=true); shared_project = Base.active_project(); inputs = vcat(read(bridge_project), read(shared_project), read(joinpath(dirname(shared_project), "Manifest.toml"))); println("julia-depot-" * bytes2hex(sha256(inputs))[1:16])'

# CI only: run job {{shard}} of the Acceptance plan for corpus {{selection}} (sample or full);
# results land in .tmp/ci/results/<family>-<shard>.
ci-acceptance family shard shards selection:
    #!/usr/bin/env bash
    set -euo pipefail
    out={{ci_dir}}/results/{{family}}-{{shard}}
    mkdir -p "$out"
    git rev-parse HEAD > "$out/commit"
    # The provider checkout is created by `ci-provision-sage` at
    # {{ci_dir}}/sage-extensions, at the revision pyproject.toml pins. The artifact keeps its
    # historical name `research_commit` because ci/acceptance_summary.py and its tests read that
    # file; only the path below is the live one. Record the revision actually checked out and
    # refuse a checkout that drifted from the pin, so the provenance the summary reports is the
    # provenance the tests ran against.
    provider_dir={{ci_dir}}/sage-extensions
    if [ ! -d "$provider_dir/.git" ]; then
        echo "ERROR: no provider checkout at $provider_dir; run 'just ci-provision-sage' first" >&2
        exit 1
    fi
    provider_head="$(git -C "$provider_dir" rev-parse HEAD)"
    pinned_rev="$(sed -n 's|.*sage-extensions @ git+https://github.com/dzackgarza/sage-extensions@\([0-9a-f]\{40\}\).*|\1|p' pyproject.toml | head -1)"
    if [ -z "$pinned_rev" ]; then
        echo "ERROR: pyproject.toml declares no sage-extensions revision" >&2
        exit 1
    fi
    if [ "$provider_head" != "$pinned_rev" ]; then
        echo "ERROR: provider checkout $provider_head does not match the pyproject.toml pin $pinned_rev" >&2
        exit 1
    fi
    printf '%s\n' "$provider_head" > "$out/research_commit"
    sage_python="$("$SAGE_BIN" -c 'import sys; print(sys.executable)')"
    export PYTHONPATH="$PWD/{{ci_dir}}/ai-review-ci/tool-artifacts/pytest_plugins:$PWD/src:$PWD"
    # ci/plan_acceptance.py packs the corpus's nodes into jobs by their derived budgets.
    mapfile -t tests < <(python3 ci/plan_acceptance.py nodes {{selection}} {{shard}})
    "$sage_python" ci/pytest_batches.py --python "$sage_python" --out "$out" --batch-size 10 -- "${tests[@]}" 2>&1 | tee "$out/log"
    exit "${PIPESTATUS[0]}"

# Run one acceptance node with a fixed cutoff and no budget table, for timing probes
# (.github/workflows/probe.yml); the report log records its outcome and duration.
ci-probe node cutoff:
    #!/usr/bin/env bash
    set -euo pipefail
    out={{ci_dir}}/probe
    mkdir -p "$out"
    git rev-parse HEAD > "$out/commit"
    git -C {{ci_dir}}/sage-extensions rev-parse HEAD > "$out/sage_extensions_commit"
    printf '%s\n' '{{node}}' > "$out/node"
    sage_python="$("$SAGE_BIN" -c 'import sys; print(sys.executable)')"
    export PYTHONPATH="$PWD/{{ci_dir}}/ai-review-ci/tool-artifacts/pytest_plugins:$PWD/src:$PWD"
    "$sage_python" -c 'import sage.all, sys, pytest; sys.exit(pytest.main(sys.argv[1:]))' \
        -p qc_sage_session -p no:cacheprovider -q --timeout={{cutoff}} --timeout-method=signal \
        --durations=0 --durations-min=0 --report-log="$out/report.jsonl" '{{node}}' 2>&1 | tee "$out/log"
    exit "${PIPESTATUS[0]}"
