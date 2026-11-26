ai_review_ci_schema_version := "1"
ai_review_ci_profile := "sage"
ai_review_ci_ref := "main"
ai_review_ci_release_channel := "main"
ai_review_ci_workflow_template_version := "1"
ai_review_ci_local_delegation := "global-justfile"
ai_review_ci_default_branch := "main"
[doc("List available recipes.")]
default:
    @just --list

[doc("Run ambient SageMath QC through the central implementation.")]
ambient:
    @just -f ~/ai-review-ci/justfiles/sage.just -d . ambient

[doc("Run commit-tier SageMath QC through the central implementation.")]
test-commit:
    @just -f ~/ai-review-ci/justfiles/sage.just -d . test-commit

[doc("Syntax checks on outgoing commits (not a full SageMath test run).")]
test-push:
    @just -f ~/ai-review-ci/justfiles/sage.just -d . test-push

[doc("Run CI acceptance QC through the central implementation.")]
test-ci:
    @just -f ~/ai-review-ci/justfiles/sage.just -d . test-ci

[doc("Inner mathematical feedback loop: run exactly one public regression node.")]
test-node node:
    #!/usr/bin/env bash
    set -euo pipefail
    node={{quote(node)}}
    case "$node" in
        tests/*::* ) ;;
        *) echo "Expected one explicit test node: tests/path.py::test_claim[parameter]" >&2; exit 2 ;;
    esac
    node_literal="$(python3 -c 'import json, sys; print(json.dumps(sys.argv[1]))' "$node")"
    export PYTHONPATH="$PWD/src:$PWD${PYTHONPATH:+:$PYTHONPATH}"
    "${SAGE_BIN:-sage}" -c "import sage.all, sys; from ci.one_node import run; sys.exit(run($node_literal))"

[doc("Run the README example as a doctest; it is the smallest acceptance.")]
readme-example:
    #!/usr/bin/env bash
    set -euo pipefail
    if [ -n "${SAGE_BIN:-}" ]; then
        sage_python="$("$SAGE_BIN" -c 'import sys; print(sys.executable)')"
    else
        sage_python="$(just -f ~/ai-review-ci/justfiles/sage.just -d . _sage-python)"
    fi
    PYTHONPATH=src "$sage_python" -m doctest README.md

ci_dir := ".tmp/ci"

[doc("CI only: install the acceptance dependencies into the Sage image's interpreter.")]
ci-provision-sage:
    #!/usr/bin/env bash
    set -euo pipefail
    if [[ -n "${MAXIMA_PREFIX:-}" && -d "${MAXIMA_PREFIX}/share/to_poly_solve" && ! -d "${MAXIMA_PREFIX}/to_poly_solve" ]]; then
        export MAXIMA_PREFIX="${MAXIMA_PREFIX}/share"
    fi
    if [[ -n "${MAXIMA_PREFIX:-}" && -n "${GITHUB_ENV:-}" ]]; then
        printf 'MAXIMA_PREFIX=%s\n' "$MAXIMA_PREFIX" >> "$GITHUB_ENV"
    fi
    export MESONPY_EDITABLE_VERBOSE=1
    "$SAGE_BIN" -c 'import sage.all; from sage.features.bliss import Bliss; Bliss().require(); from sage.interfaces.maxima_lib import maxima_lib; assert str(maxima_lib.eval("1+1")) == "2"'
    sage_python="$("$SAGE_BIN" -c 'import sys; print(sys.executable)')"
    uv pip install --quiet --python "$sage_python" "pytest>=9.1,<10" pyyaml "hypothesis>=6.165,<7" pytest-shard pytest-instafail "pytest-timeout>=2.4" "pytest-reportlog>=1.0" PyNormaliz "pycddlib>=3.0.2"
    "$sage_python" -c 'import cdd.gmp'
    mkdir -p {{ci_dir}}
    git clone --quiet --depth 1 https://github.com/dzackgarza/sage-extensions {{ci_dir}}/sage-extensions
    sage_python="$("$SAGE_BIN" -c 'import sys; print(sys.executable)')"
    uv pip install --quiet --python "$sage_python" "sage-julia-bridge @ git+https://github.com/dzackgarza/sage-julia-bridge"
    uv pip install --quiet --python "$sage_python" --no-deps -e {{ci_dir}}/sage-extensions
    uv pip install --quiet --python "$sage_python" --no-deps -e .
    git clone --quiet --depth 1 https://github.com/dzackgarza/ai-review-ci {{ci_dir}}/ai-review-ci

[doc("CI only: the bridge's Julia environment, instantiated and precompiled into the default depot.")]
ci-provision-julia:
    #!/usr/bin/env bash
    set -euo pipefail
    sage_python="$("$SAGE_BIN" -c 'import sys; print(sys.executable)')"
    env_dir="$("$sage_python" -c 'import pathlib, sage_julia_bridge; print(pathlib.Path(sage_julia_bridge.__file__).parent / "julia_env")')"
    JULIA_PKG_PRECOMPILE_AUTO=0 julia --project="$env_dir" -e 'using Pkg; bridge_project = Base.active_project(); Pkg.activate("v$(VERSION.major).$(VERSION.minor)"; shared=true); Pkg.instantiate(); Pkg.activate(dirname(bridge_project)); Pkg.instantiate(); using Oscar, JSON'

[doc("CI only: run one long true witness from witnesses/ against the provider at the head of its main branch.")]
ci-witness file test:
    #!/usr/bin/env bash
    set -euo pipefail
    out={{ci_dir}}/results/{{test}}
    mkdir -p "$out"
    git rev-parse HEAD > "$out/commit"
    git -C {{ci_dir}}/sage-extensions rev-parse HEAD > "$out/provider_commit"
    sage_python="$("$SAGE_BIN" -c 'import sys; print(sys.executable)')"
    export PYTHONPATH="$PWD/{{ci_dir}}/ai-review-ci/tool-artifacts/pytest_plugins:$PWD/src:$PWD"
    export ACCEPTANCE_CORPUS=witness
    "$sage_python" -c 'import sage.all, sys, pytest; sys.exit(pytest.main(sys.argv[1:]))' \
        -p qc_sage_session "witnesses/{{file}}::{{test}}" -vv -rA --tb=long --instafail \
        -p ci.test_budgets --timeout-method=signal --report-log="$out/report.jsonl" --junitxml="$out/junit.xml" 2>&1 | tee "$out/log"
    exit "${PIPESTATUS[0]}"

[doc("CI only: key Julia packages by the bridge project and shared image lock.")]
ci-julia-key:
    #!/usr/bin/env bash
    set -euo pipefail
    sage_python="$("$SAGE_BIN" -c 'import sys; print(sys.executable)')"
    env_dir="$("$sage_python" -c 'import pathlib, sage_julia_bridge; print(pathlib.Path(sage_julia_bridge.__file__).parent / "julia_env")')"
    julia --startup-file=no --project="$env_dir" -e 'using Pkg, SHA; bridge_project = Base.active_project(); Pkg.activate("v$(VERSION.major).$(VERSION.minor)"; shared=true); shared_project = Base.active_project(); inputs = vcat(read(bridge_project), read(shared_project), read(joinpath(dirname(shared_project), "Manifest.toml"))); println("julia-depot-" * bytes2hex(sha256(inputs))[1:16])'

[doc("CI only: run job {{shard}} of the Acceptance plan for corpus {{selection}} (tier1, tier2, full or probe); results land in .tmp/ci/results/<family>-<shard>.")]
ci-acceptance family shard shards selection:
    #!/usr/bin/env bash
    set -euo pipefail
    out={{ci_dir}}/results/{{family}}-{{shard}}
    mkdir -p "$out"
    git rev-parse HEAD > "$out/commit"
    provider_dir={{ci_dir}}/sage-extensions
    if [ ! -d "$provider_dir/.git" ]; then
        echo "ERROR: no provider checkout at $provider_dir; run 'just ci-provision-sage' first" >&2
        exit 1
    fi
    provider_head="$(git -C "$provider_dir" rev-parse HEAD)"
    printf '%s\n' "$provider_head" > "$out/research_commit"
    sage_python="$("$SAGE_BIN" -c 'import sys; print(sys.executable)')"
    export PYTHONPATH="$PWD/{{ci_dir}}/ai-review-ci/tool-artifacts/pytest_plugins:$PWD/src:$PWD"
    export ACCEPTANCE_CORPUS={{selection}}
    mapfile -t tests < <(python3 ci/plan_acceptance.py nodes {{selection}} {{shard}})
    "$sage_python" ci/pytest_batches.py --python "$sage_python" --out "$out" --batch-size 10 -- "${tests[@]}" 2>&1 | tee "$out/log"
    exit "${PIPESTATUS[0]}"
