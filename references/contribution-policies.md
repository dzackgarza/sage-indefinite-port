# Imported Contribution Policies

Moved verbatim from AGENTS.md. These bind every port edit.

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
