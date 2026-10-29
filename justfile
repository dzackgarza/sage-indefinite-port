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
    "$SAGE_BIN" -python -c 'import sage.all; from sage.features.bliss import Bliss; Bliss().require()'
    sage_python="$("$SAGE_BIN" -python -c 'import sys; print(sys.executable)')"
    uv pip install --quiet --python "$sage_python" "pytest>=9.1,<10" pyyaml "hypothesis>=6.165,<7" pytest-shard pytest-instafail PyNormaliz
    # The acceptance contract is against this exact research API revision.
    # Tracking main is not reproducible and currently imports a Sage category
    # absent from the pinned CI Sage image (sage.categories.ordered_fields).
    # Keep the clone shallow and fetch the exact provider commit explicitly.
    mkdir -p {{ci_dir}}
    git clone --quiet --depth 1 https://github.com/dzackgarza/research {{ci_dir}}/research
    research_rev="6e4078ea446ceec30f0ce8afc0bab7863387fac9"
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

# CI only: run one acceptance shard, or a required amended branch sample;
# results land in .tmp/ci/results/<family>-<shard>.
ci-acceptance family shard shards selection="":
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
    case "{{selection}}" in
      reflective)
        tests=(
            tests/acceptance/test_reflective.sage::test_edgewalk_branch_representative_sample
            tests/acceptance/test_reflective.sage::test_vinberg_roots_match_oscars_chamber_up_to_isometry
            tests/acceptance/test_reflective.sage::test_e10_simple_roots_match_dutour_sikiric_hulek
        )
        shard_args=()
        ;;
      orthogonal_group)
        tests=(
            tests/acceptance/test_orthogonal_group.sage::test_mertens_omega_and_minus_one_generate_o_l
            'tests/acceptance/test_orthogonal_group.sage::test_lorentzian_generators_match_polyhedral_common[LORENTZ_OrthogonalGroup0]'
            tests/acceptance/test_orthogonal_group.sage::test_orthogonal_group_of_the_hyperbolic_plane_has_order_four
            tests/acceptance/test_orthogonal_group.sage::test_discriminant_image_of_a7_plus_diagonal_is_all_of_o_q
            tests/acceptance/test_orthogonal_group.sage::test_signed_discriminant_image_of_u_plus_minus_three_is_all_of_o_q
        )
        shard_args=()
        ;;
      vector_stabilizers)
        tests=(
            tests/acceptance/test_vector_stabilizers.sage::test_isotropic_pointed_automorphism_generators_preserve_marking_and_form
            tests/acceptance/test_vector_stabilizers.sage::test_divisibility_two_isotropic_pointed_isometry_has_exact_witness
            tests/acceptance/test_vector_stabilizers.sage::test_odd_lattice_isotropic_pointed_isometry_uses_general_divisibility_two_route
            tests/acceptance/test_vector_stabilizers.sage::test_general_isotropic_pointed_isometry_carries_negative_marking
            tests/acceptance/test_vector_stabilizers.sage::test_empty_integral_fibre_falls_back_to_full_marked_vector_stabilizer
            tests/acceptance/test_vector_stabilizers.sage::test_split_divisibility_one_pointed_isometry_carries_actual_negative_marking
            tests/acceptance/test_vector_stabilizers.sage::test_isotropic_pointed_isometry_rejects_different_divisibility
            tests/acceptance/test_vector_stabilizers.sage::test_divisibility_four_isotropic_stabilizer_constructs_exact_containing_bound
            tests/acceptance/test_vector_stabilizers.sage::test_divisibility_four_isotropic_plane_stabilizer_constructs_exact_containing_bound
            tests/acceptance/test_vector_stabilizers.sage::test_vector_stabilizer_has_oscars_order
            'tests/acceptance/test_vector_stabilizers.sage::test_e10_stabilizer_indices_match_brandhorst_gonzalez_alonso[0]'
            'tests/acceptance/test_vector_stabilizers.sage::test_e10_stabilizer_indices_match_brandhorst_gonzalez_alonso[2]'
            'tests/acceptance/test_vector_stabilizers.sage::test_e10_stabilizer_indices_match_brandhorst_gonzalez_alonso[4]'
            'tests/acceptance/test_vector_stabilizers.sage::test_e10_stabilizer_indices_match_brandhorst_gonzalez_alonso[6]'
            'tests/acceptance/test_vector_stabilizers.sage::test_e10_stabilizer_indices_match_brandhorst_gonzalez_alonso[8]'
            'tests/acceptance/test_vector_stabilizers.sage::test_e10_stabilizer_indices_match_brandhorst_gonzalez_alonso[10]'
            'tests/acceptance/test_vector_stabilizers.sage::test_equivalence_witness_carries_a_representative_to_its_reflection[r_-1]'
            'tests/acceptance/test_vector_stabilizers.sage::test_equivalence_witness_carries_a_representative_to_its_reflection[r_0]'
            'tests/acceptance/test_vector_stabilizers.sage::test_equivalence_witness_carries_a_representative_to_its_reflection[r_1]'
            'tests/acceptance/test_vector_stabilizers.sage::test_equivalence_witness_carries_a_representative_to_its_reflection[r_2]'
            'tests/acceptance/test_vector_stabilizers.sage::test_equivalence_witness_carries_a_representative_to_its_reflection[r_3]'
            'tests/acceptance/test_vector_stabilizers.sage::test_equivalence_witness_carries_a_representative_to_its_reflection[r_4]'
            'tests/acceptance/test_vector_stabilizers.sage::test_equivalence_witness_carries_a_representative_to_its_reflection[r_5]'
            'tests/acceptance/test_vector_stabilizers.sage::test_equivalence_witness_carries_a_representative_to_its_reflection[r_6]'
            'tests/acceptance/test_vector_stabilizers.sage::test_equivalence_witness_carries_a_representative_to_its_reflection[r_7]'
            'tests/acceptance/test_vector_stabilizers.sage::test_equivalence_witness_carries_a_representative_to_its_reflection[r_8]'
        )
        shard_args=()
        ;;
      vector_orbits)
        tests=(
            tests/acceptance/test_vector_orbits.sage::test_anisotropic_binary_vector_classes_use_primitive_embeddings
            'tests/acceptance/test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_2]'
            'tests/acceptance/test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_30]'
            'tests/acceptance/test_vector_orbits.sage::test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[0]'
            'tests/acceptance/test_vector_orbits.sage::test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[2]'
            'tests/acceptance/test_vector_orbits.sage::test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[4]'
            'tests/acceptance/test_vector_orbits.sage::test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[6]'
            'tests/acceptance/test_vector_orbits.sage::test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[8]'
            'tests/acceptance/test_vector_orbits.sage::test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[10]'
            'tests/acceptance/test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[2]'
            'tests/acceptance/test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[4]'
            'tests/acceptance/test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[6]'
            'tests/acceptance/test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[8]'
            'tests/acceptance/test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[10]'
            'tests/acceptance/test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[12]'
            'tests/acceptance/test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[14]'
            'tests/acceptance/test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[16]'
            'tests/acceptance/test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[18]'
            'tests/acceptance/test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[20]'
            'tests/acceptance/test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[22]'
            'tests/acceptance/test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[24]'
            'tests/acceptance/test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[26]'
            'tests/acceptance/test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[28]'
            'tests/acceptance/test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[30]'
            tests/acceptance/test_vector_orbits.sage::test_i_2_10_is_transitive_on_norm_minus_one_vectors
            'tests/acceptance/test_vector_orbits.sage::test_k3_lattice_has_one_primitive_orbit_per_norm[0]'
            'tests/acceptance/test_vector_orbits.sage::test_k3_lattice_has_one_primitive_orbit_per_norm[-2]'
            'tests/acceptance/test_vector_orbits.sage::test_k3_lattice_has_one_primitive_orbit_per_norm[2]'
            'tests/acceptance/test_vector_orbits.sage::test_k3_lattice_has_one_primitive_orbit_per_norm[4]'
            'tests/acceptance/test_vector_orbits.sage::test_k3_lattice_has_one_primitive_orbit_per_norm[-4]'
        )
        shard_args=()
        ;;
      isotropic_orbits)
        tests=(
            'tests/acceptance/test_isotropic_orbits.sage::test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors[IsotropicCases_0]'
            'tests/acceptance/test_isotropic_orbits.sage::test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors[IsotropicCases_15]'
            tests/acceptance/test_isotropic_orbits.sage::test_o_plus_n_has_two_line_and_two_plane_orbits
            tests/acceptance/test_isotropic_orbits.sage::test_stable_o_plus_n_orbits_match_dutour_sikiric_hulek
            tests/acceptance/test_isotropic_orbits.sage::test_o_n_maps_onto_o_of_the_discriminant_form
            'tests/acceptance/test_isotropic_orbits.sage::test_published_representatives_lie_in_distinct_o_plus_orbits[rank0]'
            'tests/acceptance/test_isotropic_orbits.sage::test_published_representatives_lie_in_distinct_o_plus_orbits[rank1]'
            'tests/acceptance/test_isotropic_orbits.sage::test_stabilizer_images_have_dutour_sikiric_hulek_indices[I_1_1]'
            'tests/acceptance/test_isotropic_orbits.sage::test_stabilizer_images_have_dutour_sikiric_hulek_indices[I_1_2]'
            'tests/acceptance/test_isotropic_orbits.sage::test_stabilizer_images_have_dutour_sikiric_hulek_indices[I_2_1]'
            'tests/acceptance/test_isotropic_orbits.sage::test_stabilizer_images_have_dutour_sikiric_hulek_indices[I_2_2]'
            'tests/acceptance/test_isotropic_orbits.sage::test_i_2_10_isotropic_classes_are_one_odd_and_one_even[rank0-primitive_isotropic_vector_orbits]'
            'tests/acceptance/test_isotropic_orbits.sage::test_i_2_10_isotropic_classes_are_one_odd_and_one_even[rank1-isotropic_plane_orbits]'
        )
        shard_args=()
        ;;
      boundaries)
        tests=(
            'tests/acceptance/test_boundaries.sage::test_dawes_buildings_have_the_published_counts[dawes_2U_A2_stable]'
            'tests/acceptance/test_boundaries.sage::test_dawes_buildings_have_the_published_counts[dawes_2U_minus2_minus6]'
            'tests/acceptance/test_boundaries.sage::test_dawes_buildings_have_the_published_counts[dawes_L2_gamma2]'
            'tests/acceptance/test_boundaries.sage::test_dawes_building_incidence_graph_is_the_published_figure[dawes_2U2_A2_stable]'
            tests/acceptance/test_boundaries.sage::test_dawes_last_index_in_the_chain_is_the_discriminant_image_of_o_plus
            'tests/acceptance/test_boundaries.sage::test_dtower_boundary_component_counts[N_3]'
            'tests/acceptance/test_boundaries.sage::test_dtower_boundary_component_counts[N_4]'
            'tests/acceptance/test_boundaries.sage::test_dtower_boundary_component_counts[N_9]'
            'tests/acceptance/test_boundaries.sage::test_dtower_boundary_component_counts[N_20]'
            'tests/acceptance/test_boundaries.sage::test_dtower_boundary_pictures[lambda_9]'
            'tests/acceptance/test_boundaries.sage::test_dtower_boundary_pictures[lambda_10]'
            'tests/acceptance/test_boundaries.sage::test_dtower_boundary_pictures[lambda_11]'
            tests/acceptance/test_boundaries.sage::test_degree_two_k3_type_ii_components_and_root_types
            tests/acceptance/test_boundaries.sage::test_degree_four_k3_has_nine_type_ii_components
            'tests/acceptance/test_boundaries.sage::test_enriques_modular_groups_match_dutour_sikiric_hulek[case_1]'
            'tests/acceptance/test_boundaries.sage::test_enriques_modular_groups_match_dutour_sikiric_hulek[case_87]'
            'tests/acceptance/test_boundaries.sage::test_classical_polarization_indices[degree_2]'
            'tests/acceptance/test_boundaries.sage::test_classical_polarization_indices[degree_10]'
        )
        shard_args=()
        ;;
      *)
        tests=(tests/acceptance/test_{{family}}.sage)
        shard_args=(--shard-id={{shard}} --num-shards={{shards}})
        ;;
    esac
    if [[ "{{family}}" == reflective && "{{selection}}" == "" ]]; then
        timeout --signal=INT 345m python3 ci/pytest_batches.py --sage "$SAGE_BIN" --out "$out" --batch-size 10 -- "${tests[@]}" "${shard_args[@]}" 2>&1 | tee "$out/log"
        exit "${PIPESTATUS[0]}"
    fi
    timeout --signal=INT 345m "$SAGE_BIN" -python -c 'import sage.all, sys, pytest; sys.exit(pytest.main(sys.argv[1:]))' \
        -p qc_sage_session "${tests[@]}" --runxfail -rA --tb=long --instafail \
        "${shard_args[@]}" --junitxml="$out/junit.xml" 2>&1 | tee "$out/log"
