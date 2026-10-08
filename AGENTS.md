# Review Guidelines

These are additional requirements for reviewing agent work.
They do not replace the reviewer’s normal role, repo-specific standards, or technical judgment.
They provide the failure model that should shape the review.

The task is not merely to review a PR. The task is to decide whether a completion claim is true under the original objective.
The standard is full, correct, provable completion against the original requirements and repo guidelines.
Anything less is incomplete work that must not be treated as a win.

## Failure Model

Agents systematically produce impressive non-completion.
Common patterns are: polished summaries that imply finished work, caveats that quietly narrow the goal, reclassification without proof, delegated discovery presented as resolution, process language that substitutes for evidence, merged PRs treated as completion, passing checks treated as semantic proof, and artifacts that look substantial while leaving required work unowned.

Treat the agent’s summary, PR description, closing comment, issue closure, “goal completed” statement, and self-reported validations as untrusted.
They may be diagnostic pointers, but they are not evidence that the work is complete.
The evidence is the original issue or task, the code diff, tests, source/runtime facts, review comments, and produced artifacts.

## Decisive Invariants

Preserve the original success condition.
Read the original issue or task before accepting any restatement of it.
Keep its quantifiers intact: “all,” “complete,” "full subset," “zero remaining,” and similar terms cannot be quietly narrowed to examples, partial coverage, known blockers, or whatever the PR happened to touch.

Nothing required may disappear silently.
A required work family must be implemented, explicitly falsified, or validly reclassified with evidence that satisfies the issue’s own standard.
Partial implementation is not completion.
Future work is not completion.
Count reduction is not completion.
Resolved review threads are not completion.
Passing checks are not completion.
Substantial-looking work is not completion.
“Better than before” is not completion.

Goal substitution is the main thing to detect.
Ask whether the submitted work solves the original problem or merely produces a narrower artifact: cleaner metadata, a partial subset, a better explanation, a new issue, a renamed scope, a local workaround, or proof that someone should investigate later.

Technically correct administrative artifacts can be goal substitution.
A well-written issue, comment, audit note, scope statement, or enumeration of remaining work may be required, but it does not complete implementation, testing, proof, or downstream cleanup.
If the original task requires execution, the artifact is only useful insofar as it drives that execution; it must not become the stopping point.

Treat self-scoped remaining-work lists as a severe completion-laundering pattern.
When an agent is asked to enumerate remaining work, the domain is the original full completion requirement, not the agent’s intended subset, the PR’s current shape, a closeability criterion, or the work left after deferral and reclassification.
A valid enumeration subtracts only artifact-proven completed work from the original contract.
Deferrals, routed follow-ups, owner changes, and truthful incompletion notes remain unresolved work unless the original task explicitly made that administrative routing the whole deliverable.

If an agent repeats a narrowed enumeration after being corrected, treat that as a hard misalignment signal, not as an innocent wording issue.
The reviewer should identify the original full requirement, the scope the agent substituted, and the required work hidden by that substitution.

Silent reclassification is not resolution.
If the PR says remaining work is out-of-scope, research-owned, stub-owned, plugin-owned, downstream-owned, or future-owned, require evidence from the relevant source/runtime behavior, repo boundary, or original acceptance criteria.
A sentence in the PR description is not enough.

Ownership boundaries matter.
The submitting repo must prove its own claimed behavior and do the blocker forensics required by its own issue.
Do not require a receiving or downstream repo to classify another project’s internal uncertainty unless the original issue explicitly made that part of acceptance.
When an external issue is created, it should be written for that receiving repo, not for a reader who already knows the submitting repo’s context.

## Evidence Expectations

Review tests as evidence, not as decoration.
Valid tests exercise the real production path or semantic requirement.
Be skeptical of helper-only tests, tautologies, assertions of the implementation’s own output, bypasses around the runtime/plugin/stub path, example-only coverage where the issue required full coverage, weakened assertions, and missing invalid-nearby cases where the fix could overgeneralize.

For plugin work, the evidence should usually distinguish valid generic behavior from invalid nearby ordinary Python and should not hard-code a downstream consumer.
For stubs work, the evidence should be source-backed: the upstream surface exists, the stub matches public behavior, no fake API is added, no Any/object opacity escape is introduced, and inherited-method inflation is not used unless source exposes that surface.

Watch for code-level laundering: hard-coded consumer names, support for local research abstractions as if they were external API, fake stubs, broad Any/object escapes, line suppressions, diagnostic filtering, deletion of required data, broad type widening, and any move that makes checks pass by weakening the problem instead of solving it.

## When Acting on Review Feedback

A positive disposition requires a commit.

Do not resolve an accepted review comment until the code/proof remediation is committed and the reply cites the commit.

Never reply “accepted,” “aligned,” “fixed,” “addressed,” or “will address” to a review thread unless the remediation is already committed.
A thread cannot be resolved on intent or future work.

Every substantive review item must receive its visible thread- or surface-local disposition and evidence before resolution.
The canonical field contract and state machine live in [[pr-feedback-triage/SKILL|pr-feedback-triage]]. Do not create top-level disposition ledgers or tracked review-log files.
Migrate legacy ledger-only resolutions by posting the canonical disposition and evidence on each affected thread before treating it as closed.

Review comments are not implementation specs.
The worker must translate accepted feedback into first-principles remediation requirements before assigning implementation.

For each comment:
- Identify the concern.
- Identify the proposed fix.
- Decide whether the concern is true under global + repo policy.
- Decide whether the proposed fix preserves those policies.
- If the concern is true but the fix is wrong, apply a policy-compatible remediation.

## Writing the Review

Write nuanced feedback for an intelligent reader.
Do not force a machine-readable template, a mandatory table, or a simplistic pass/fail label when prose communicates the situation better.
Do make the completion judgment clear: whether the original task can be considered complete, what evidence supports that judgment, and which unresolved requirements block completion if any remain.

Do not foreground effort, progress, good intentions, volume of work, or “substantial” partial implementation when required work remains.
Mention completed pieces only when they are necessary to identify the exact remaining blockers or to prevent redoing already-correct work.
Do not compare incomplete work to “no work done” or “completely fake work”; compare it to the expected standard: the task done correctly, completely, and provably.

When required work remains, lead with the incompleteness and the concrete blockers.
Do not make the reader excavate the missing work from beneath praise, context-setting, or a narrative of what did get done.

Nuance belongs in the evidence and blocker analysis, not in softening the completion standard.
The review should make it easy to finish the work, not easy to feel satisfied with less than the original contract required.

<!-- agent-memory:start -->
# Agent memory

This repository uses the central agent memory vault at `/home/dzack/.agent-memory-vault`.

Project memory key: `projects/github.com__dzackgarza__sage-indefinite-port/index`.

Repository `.agents` and `.hermes` paths are symlinks to the same vault-owned project directory.

Before changing architecture, search both project and global memory:

```bash
agent-memory search --scope both "<task or subsystem>"
```

Record durable repo-specific lessons with:

```bash
agent-memory add --scope project --type decision --title <title> --content <content>
agent-memory add --scope project --type trap --title <title> --content <content>
agent-memory add --scope project --type advice --title <title> --content <content>
agent-memory add --scope project --type context --title <title> --content <content>
agent-memory add --scope project --type reference --title <title> --content <content>
```

Plan work is card-backed. Create and update plan cards with `agent-memory plan add` and `agent-memory plan update`, not `agent-memory add --type plan`.

Use `agent-memory retrieve <key>`, `agent-memory update <key>`, and `agent-memory delete <key>` for memory CRUD.

The vault should be committed at all times. Treat staged or unstaged vault changes as an ephemeral error state. Before normal memory work resumes, load the bundled vault-maintenance skill with `agent-memory maintain skill vault-maintenance` and follow its referenced check, repair, and commit workflows.

Move reusable lessons during maintenance with:

```bash
agent-memory maintain move <key> --to global/advice
```
<!-- agent-memory:end -->

# Taking Work

The GitHub issue tree rooted at #4 is the only plan: its work units are the leaf issues, their blocked-by links are the dependency edges, and each unit's body holds its obligation and acceptance. Issue bodies cite "map §N" in [references/source-to-sage-translation-map.md](references/source-to-sage-translation-map.md).
Take the next unit with `uvx --from git+https://github.com/dzackgarza/itree itree next dzackgarza/sage-indefinite-port`, work it to the acceptance in its issue body, commit, close it with `itree close`, and take the next unit in the same turn.
A turn that ends with a ready unit untouched stops the repository until someone notices.

# Commit Gates and Work Claims

- **Red gate protocol**: The first time a commit gate, hook, or QC stage goes red, diagnosing that failure becomes the current task. Stop authoring; root-cause and fix the gate, or report it as a blocker with a reproducer. Never keep writing code behind a red gate, and never accumulate uncommitted work around one. A gate that is red on two consecutive commit attempts is a defect to diagnose, not an environment condition to wait out.

- **Claim freshness — no off-ledger work**: At every claim and every release, reconcile the shared queue/claim state against actual repository state across all branches before selecting work. Never select work from a queue older than your last branch sync. All authoring requires a live claim; batch-committing a body of work authored off-ledger is prohibited.

- **The `[unverified]` commit tag**: `[unverified]` may appear in a commit message only while the repository's verification phase is formally deferred by a named contract, and every `[unverified]` commit must name the contract or terminal phase that discharges it.

- **Deferred type checking (temporary, discharged by #36)**: while #36 is open, the push gate's mypy errors are acknowledged debt and are not paid down. A push may then use `git push --no-verify`, and only under all of these conditions:
  - `just test-push` was run on the commit being pushed and failed in its mypy stage (`_sage-mypy`) and nowhere earlier.
  - The tip commit of the push carries this trailer, verbatim apart from the count: `Type-Debt-Ack: MYPY-DEFER-7Q4K #36 errors=<N>`, where `<N>` is the "Found N errors" count of that run.
  - `<N>` is no greater than the count in the previous `Type-Debt-Ack` trailer in `git log`. A push that would raise the count does not use this exception: fix the new errors first.

  The trailer covers every commit in that push. Commits need no exception, because the commit gate does not run mypy. The exception covers only mypy debt: any other gate failure is fixed, never bypassed. #36's last commit deletes this rule.

# Architecture and Dependencies

- **The `research` preamble is the lattice substrate and a co-developed dependency**:
  This repository consumes `dzack_research.preamble` (categories, formed modules, lattices, isometries, discriminant modules) from the checkout at `/home/dzack/research/src`. `.envrc` puts it on `PYTHONPATH` and `MYPYPATH`; `pyproject.toml` declares `dzack-research`.
  - **Use the preamble's own classes and types**: annotate with `Lattice` (and `Lattice.Element`), `LatticeIsometry`, `LatticeIsometryHomset`, and the tensor and presented-module classes. Never add Protocol, stub, or facade types over the preamble, and never edit `sys.path`.
  - **Fix preamble defects upstream**: a missing annotation, a class reachable only from a private module, or a wrong result in the preamble is fixed in `research` under its own QC, then consumed here. No local workaround, no fallback to raw SageMath lattices (their indefinite support is insufficient), no reinvention of lattice, discriminant, or formed-category infrastructure here.
  - **Contribution direction**: the indefinite algorithms built here move upstream into `research` once the plan's phases are complete.

- **NEVER BUILD THE C++ CODE**:
  **Do NOT attempt to compile, build, configure, or invoke C++ compilation toolchains** for `polyhedral_common` or any other C++ source in `references/`.
  The upstream C++ codebase in `references/polyhedral_common` is strictly a reference implementation for algorithm extraction, logic translation, and structural understanding.
  All production algorithms in this repository MUST be implemented as native SageMath / Python code delegating low-level operations to standard system backends (FLINT for exact integer linear algebra, Normaliz / cddlib / PPL for polyhedral cones). Basic mathematical operations (stabilizers, transporters, orbits, cosets and double cosets of group actions; kernels, preimages, centralizers and intersections of subgroups; isomorphisms and canonical forms of finite configurations) are the research preamble's. This repository calls them and never reimplements them, whether on libGAP, Bliss, or anything else. A missing preamble operation is an upstream gap: name its API in the acceptance contract and leave the cases that need it red.

## Imported Contribution Policies (from `/home/dzack/research/CONTRIBUTING.md`)

### ARC: Mathematical Architecture & Ownership
- `ARC-01`: Own universal properties, categories, morphisms, functors, and adjunctions natively in this repository’s category framework.
- `ARC-02`: Represent subobjects as pairs `(S, iota: S -> M)`. Place predicates, isometries, embeddings, and containment checks on morphism data and hom-sets.
- `ARC-03`: Build structural objects/functors first, then derive numerical invariants.

### ENG: Computational Backend Delegation
- `ENG-01`: Delegate heavy computations to reliable exact backends (SageMath, Singular, OSCAR, Macaulay2, PARI/GP) when available.
- `ENG-02`: Do not hand-roll standard mathematics that mature upstream dependencies already provide.
- `ENG-03`: Keep owned logic minimal; offload heavy numerical work to engine backends through standard bridges.
- `ENG-04`: Prefer native engine implementations (e.g., Julia/OSCAR, Singular) for multi-step heavy computation when bridge overhead is worse.

### BRG: Interoperability and Bridge Boundaries
- `BRG-01`: Use structured, persistent bridge interfaces for external systems; avoid ad-hoc subprocess scripts with temp files for core algorithms.
- `BRG-02`: Validate and wrap data at bridge boundaries to keep backend representations out of public API types.

### ENV: Environment, Execution, and Tooling
- `ENV-01`: Use exact physical paths in shell/tool invocations (e.g., `/home/dzack/...`).
- `ENV-02`: Define project orchestration, gates, and doc generators in the root `justfile`.

### DEV: Development Discipline
- `DEV-01`: Use explicit typing on public APIs; avoid `Any`/`object` in public contracts.
- `DEV-02`: Add concrete falsifiable specimens for every new category, functor, or operation.
- `DEV-03`: Before adding code, run `just preamble-megadoc`, reuse existing constructions, and implement new work in the most general mathematical form before specialization.

## Imported Contribution Policies (from `/home/dzack/gitclones/sage-categories/CONTRIBUTING.md`)

### POL-SCOPE: Implementation Order and Scope
- `POL-SCOPE-001`: Build dependency chains in order and complete required milestones before starting later phases.
- `POL-SCOPE-002`: Execute only active-phase work with accepted prerequisites. Do not start unsupported descendants early.
- `POL-SCOPE-007`: Measure success by categorical ownership, auditable declarations, and explicit mathematical structure.
- `POL-SCOPE-008`: Keep the mathematical declaration as the source of truth. Keep runtime representation private except for authorized SymPy propositions.

### POL-SHADOW: Package-owned Mathematical Surface
- `POL-SHADOW-001`: Build a package-owned categorical replacement for the supported Sage surface.
- `POL-SHADOW-003`: Keep public API closed over package-owned categories, objects, morphisms, and propositions.
- `POL-SHADOW-005`: Do not guarantee compatibility with unsupported arbitrary Sage API behavior.
- `POL-SHADOW-006`: Absorb required Sage constructions into the package-owned architecture.

### POL-ONT: Foundational Ontology
- `POL-ONT-001`: Treat raw Python values as private carriers, not as category members.
- `POL-ONT-002`: Use separate constructor and refinement layers; avoid runtime-type-sniffing constructor overloads.
- `POL-ONT-008`: Do not use catch-all fallback constructors. Route each input mode through explicit constructors and fail hard.
- `POL-ONT-010`: Survey existing generic categorical machinery before adding greenfield implementations.
- `POL-ONT-011`: Do not add standalone procedural helper functions for core predicates.

### POL-MATH: Mathematical Architecture
- `POL-MATH-001`: Identify object, element, morphism, category, functor, and universal-property owners before implementation.
- `POL-MATH-002`: Model named categories as categories, not as utility-like classes.
- `POL-MATH-003`: Model functors with explicit object and morphism maps.
- `POL-MATH-008`: Do not duplicate data that a defining morphism already determines.
- `POL-MATH-014`: Use inspected mathematics to justify construction sites; do not use tests as proof of mathematics.
- `POL-MATH-034`: Every truth question gets one category-owned predicate and a SymPy proposition; `ask()` resolves it.
- `POL-MATH-037`: Construct values in the stated mathematical category; constructors do not certify proofs.

### POL-CAT: Category Ownership and Inheritance
- `POL-CAT-001`: A category owns constructors, local operations, and role types.
- `POL-CAT-002`: Use category-owned `ObjectType`, `ElementType`, and `MorphismType`.
- `POL-CAT-004`: A category level only declares its own structure and operations.
- `POL-CAT-006`: Do not re-export methods owned by another category.
- `POL-CAT-018`: Distinguish property subcategories from “data-carrying” categories.
- `POL-CAT-021`: Keep `Mor(n, C)` and endpoint categories as first-class objects in the categorical layer.

### POL-REP: Semantic Representation
- `POL-REP-001`: Treat Sage vectors/matrices as private representations; own semantic objects in the category layer.
- `POL-REP-002`: Return semantic mathematical objects at public API boundaries.
- `POL-REP-003`: Compare elements and morphisms semantically, not by raw coordinate unpacking.
- `POL-REP-009`: Lower private semantics to computation once and reconstruct semantic outputs before returning.
- Local repository rule: do not use raw matrices as mathematical API objects. Use morphisms for structure and tensor valences for numerics (for example, `(0,2)` Gram tensors and `(1,1)` endomorphism tensors).

### POL-ENGINE: Computation Engine Boundary
- `POL-ENGINE-001`: Public API is owned mathematics; engines provide private realizations and algorithms.
- `POL-ENGINE-002`: Keep engine types private; only authorized SymPy proposition expressions may cross the boundary.
- `POL-ENGINE-007`: Do not introduce selectable backends or replaceable engine abstractions in the public API.
- `POL-ENGINE-015`: Use fixed dependency assignments for private algorithms and reconstruct exact semantic outputs.

### POL-FORM: Forms and Lattices
- `POL-FORM-001`: Model lattices as `R`-modules with a specified form, not as raw free `ZZ` modules.
- `POL-FORM-002`: Encode bilinear forms via Gram tensors (`M ⊗ M -> W`).
- `POL-FORM-004`: Do not assume positivity, freeness, embeddedness, or unimodularity by default.
- `POL-FORM-006`: Define complements, norms, and reflections under correct hypotheses only.
- `POL-FORM-008`: Use exact arithmetic and exact coefficient rings.

# Existing Preamble Lattice API (Obviates Foundational Tasks)

The user's `research` repository (`dzack_research.preamble.categories.lattices`, `lattice_morphisms`, `modules`, and `forms`) **already fully implements and provides** the complete formed lattice algebra and categorical infrastructure.
**Do NOT reinvent or re-implement any of the following infrastructure; consume it directly from `dzack_research.preamble`. This completely obviates greenfield implementation for Phase 2 and early foundational tasks.**

## Full lattice interface in `preamble-megadoc.md` (`#subsystem-lattices`)

### Construction and category entry points
- `Lattices` as the lattice category over a ring with constructor:
  - `Lattices(cls, *args)`
  - `Lattices(...)(self, data, basis=None, names=None, form=None, module_generators=None)`
- Category constructors:
  - `FiniteRankLattices`, `NondegenerateLattices`, `RationalLattices`, `EvenLattices`, `RootLattices`
- `Lattices(ZZ)` named lattice constructors:
  - `Lattices(ZZ)("U")`, `"U(m)"`, `"An"`, `"Dn"`, `"En"`, `"E8"`, `"II_{1,9}"`, `"Leech"`, and Gram-tensor constructors.
- Direct sum and scaling:
  - `L1 + L2` and `L.twist(scale)`.

### `Lattices` parent object API (full interface)
- Isometry and homset constructors:
  - `Aut`, `O`, `orthogonal_group`
  - `SO`, `special_orthogonal_group`
  - `Emb(codomain)`, `Isom(codomain)`
  - `hom(images, codomain=None)`, `identity_morphism`
  - `similarity_homset(other, scale)`, `similarity(scale, images=None, codomain=None)`
- Pairings and form evaluations:
  - `b(left, right)`, `q(vector)`, `correlation`, `correlation_morphism`
  - `gram_tensor`, `gram_matrix(basis=None)`, `discriminant`, `rank`, `signature_pair`, `level`, `genus`
- Lattice reductions and basis control:
  - `LLL()`, `BKZ(block_size=20)`, `HKZ()`
  - `lll_reduction`, `bkz_reduction`, `hkz_reduction`, `orthogonal_group`, `stable_orthogonal_group`
- Duality and discriminant data:
  - `dual_module`, `dual_lattice`, `metric_dual`, `discriminant_module`, `discriminant_projection`
  - `discriminant_group`, `discriminant_bilinear_form`, `discriminant_quadratic_form`
  - `discriminant_class`, `divided_discriminant_class`, `discriminant_representation`
  - `discriminant_image`, `discriminant_representation_is_surjective`
- Orthogonality and decomposition:
  - `radical`, `radical_quotient`, `orthogonal_complement`
  - `module_generating_set`, `module_generator`, `biproduct_factors`, `decomposition`, `decomposition_names`, `summands`
  - `is_decomposable`, `indecomposable_name`, `primitive_isotropic_subobject`
  - `local_modification(prime, *discriminant_classes)`, `overlattice(*discriminant_classes)`, `even_overlattice_inclusions`
- Arithmetical predicates:
  - `is_finite_rank`, `is_even`, `is_nondegenerate`, `is_unimodular`
  - `is_definite`, `is_positive_definite`, `is_negative_definite`, `is_p_elementary(prime)`
  - `is_locally_isometric(other, prime)`, `is_isometric(other)`, `is_similar(other, scale)`
- Discrete geometry and orbits:
  - `minimum`, `successive_minima`, `shortest_vectors`, `roots`, `roots_of_square(square)`
  - `vectors_of_square(square)`, `vectors_of_square_and_divisibility(square, divisibility)`
  - `kissing_number`, `root_sublattice`, `theta_series(precision=20, variable='q')`
- Isotropic and orbit APIs:
  - `isotropic_flag`, `isotropic_flag_orbit_representatives(rank=2)`
  - `isotropic_line_orbit_representatives`, `isotropic_plane_orbit_representatives`
  - `voronoi_cell(bound=None)`, `voronoi_relevant_vectors`
- Group and positivity APIs:
  - `stable_orthogonal_group`, `spinor_kernel_subgroup`, `positive_cone_subgroup`
  - `twist(scale)`, `reflection(root)`
  - `two_elementary_invariants`, `delta`
- Packing and metric helpers:
  - `babai(target)`, `closest_vector(target)`, `contact_polytope`, `covering_radius`
  - `packing_density`, `packing_radius`, `center_density`, `hadamard_ratio`, `hermite_invariant`

### `Lattices` element methods
- `b(other)`, `norm`, `q()`, `div()`, `divisibility_ideal`
- `divided_discriminant_class()`, `is_root()`
- `to_list`, `to_tuple`, `to_vector`, `monomial_coefficients`

### Subcategory APIs
- `FiniteRankLattices`: `is_finite_rank`
- `RationalLattices`: `fraction_field`, `is_nondegenerate`
- `EvenLattices`: `is_even` inherited by lattice parent
- `RootLattices`:
  - parent methods: `cartan_type`, `coxeter_number`, `highest_root`, `simple_reflections`, `simple_roots`, `fundamental_weights`
  - element methods: `coroot`, `height`, `is_negative_root`, `is_positive_root`

### Morphisms and homsets
- `LatticeMorphism` (`ModuleMorphism`):
  - constructor: `__init__(self, parent, images)`
- `LatticeEmbedding` (`LatticeMorphism`):
  - `discriminant_inclusion`, `is_injective`
- `LatticeIsometry` (`LatticeEmbedding`):
  - `centralizer_discriminant_image`, `cyclic_subgroup`, `determinant`, `discriminant_isometry`
  - `discriminant_morphism`, `formed_coinvariants`, `invariant_lattice`, `inverse`
  - `is_surjective`, `preserves_positive_cone`, `real_spinor_norm_sign`
- `LatticeHomset`:
  - `__init__(domain, codomain)`, `_element_constructor_(images)`
- `LatticeEmbeddingHomset`:
  - `an_element`, `even_overlattice_inclusions`, `is_empty`
- `LatticeIsometryHomset`:
  - `act`, `acting_group`, `compose`, `discriminant_image`, `discriminant_preimage`
  - `group_generators`, `identity`, `is_empty`, `isotropic_equivalence_witness`
  - `isotropic_orbit_representatives`, `isotropic_stabilizer_generators`, `number_of_group_generators`
  - `one`, `order`, `transporter`, `vector_equivalence_witness`, `vector_orbit_representatives`
  - `vector_stabilizer_generators`, `vectors_are_equivalent`

### Lattice objects and named tables
- `CoxeterDiagram`
  - `cardinality`, `connected_components`, `elliptic_subdiagrams`, `graph`, `coxeter_entry`
  - `coxeter_matrix`, `is_connected`, `is_elliptic`, `is_hyperbolic`, `is_parabolic`
  - `is_rooted`, `index_set`, `induced_subdiagram`, `preferred_positions`, `roots`
  - `root_gram_tensor`, `schlafli_tensor`, `signature_pair`, `vertex_names`
- `Genus`
  - `class_number`, `determinant`, `discriminant_form`, `excess`, `exists`
  - `local_symbol`, `level`, `mass`, `representative`, `representatives`, `signature_pair`
- `IsotropicFlag`
  - `basis`, `lattice`, `rank`, `terms`, `top`
- `OrthogonalCharacterQuotient`
  - `image`, `image_keys`, `splitting_isometries`, `stabilizer_image_keys`
  - `subgroup_image_keys`, `witness_meets_subgroup`
- `VectorPrimitiveExtension`
  - `class_of_representative`, `complement_is_definite`, `representative_of`
- registries: `register_indecomposable`, `register_indecomposable_gram`

### Helper functions from the lattice subsystem
- lattice constructors and homset factories:
  - `lattice(...)`, `lattice_embedding_homset`, `lattice_homset`, `lattice_isometry_homset`, `lattice_latex`
- lattice operators:
  - `diagonal_gram`, `orthogonal_sum`, `scale_gram_tensor`, `signature_pair_of_gram`
  - `colimit_lattice`, `discriminant_of_gram`, `signature_pair_of_gram`
- orbits, stabilizers, and equivalence:
  - `definite_complement_extensions`, `isotropic_equivalence_witness`, `isotropic_orbit_representatives`, `isotropic_stabilizer_generators`
  - `subgroup_isotropic_are_equivalent`, `subgroup_isotropic_orbit_representatives`
  - `subgroup_vector_orbit_representatives`, `subgroup_vectors_are_equivalent`
  - `transport_isotropic_object`, `vector_equivalence_witness`, `vectors_are_equivalent` (homset-level)
- reduction and algorithm helpers:
  - `babai`, `bkz_reduction`, `hkz_reduction`, `closest_vector`, `gaussian_heuristic`
  - `generator_pairings`, `gluing_route_discriminant_classes`
  - `packing_radius`, `packing_density`, `center_density`, `kissing_number`
  - `roots`, `roots_of_square`, `vectors_of_square`, `vectors_of_square_and_divisibility`
  - `voronoi_cell`, `voronoi_relevant_vectors`
  - `oscar_centralizer_discriminant_image`, `oscar_even_unimodular_primitive_embedding`, `oscar_rational_spinor_norm_sign`

### Compatibility note
- The above interface is the preamble contract. This repo must use it and should not introduce local fallback lattice mechanics in place of these names.

---

### What Actually Needs to be Built in `sage-indefinite-port`
Because foundational formed algebra is already complete, development in this repository focuses strictly on:
1. **Lorentzian perfect-domain reduction**: Traversal of the $(1,n)$ reduction complex for $O^\Omega(L)$ and $O(L)$ generators.
2. **$2U$-Eichler approximate models & Higher-Witt-index recursion**: Generating $A(L)$ and lifting covering sets to full $O(L)$ for Witt index $\geq 2$.
3. **Parabolic recursion & Isotropic orbit decomposition**: Cusp orbits, unipotent radicals, and inductive rank-$k$ isotropic plane orbits.
4. **Centralizers of isometries**: eigensublattices, gluing data and the cyclotomic decomposition, with the stabilizers and intersections consumed from the preamble.
