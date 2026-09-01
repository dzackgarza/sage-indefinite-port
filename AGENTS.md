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

# Architecture and Dependencies

- **Pure Dependency on `research` Preamble (DO NOT EDIT PREAMBLE CODE)**:
  This repository implicitly consumes the preamble code from the user's `research` repository (`dzack_research.preamble`, categories, formed modules, and semantic lattice interfaces) as an external upstream dependency.
  - **Do NOT modify preamble code**: Never edit, refactor, or touch files in `research/` or `dzack_research/preamble` during work on this repository.
  - **Do NOT probe internals**: Treat the preamble as a standard library dependency. Do not inspect its internals, audit its mechanisms, or probe its implementation unless a concrete, fatal runtime blocker occurs.
  - **Purpose of Dependency**:
    1. Prevent falling back to raw SageMath lattices, where indefinite lattice and isometry support is broken/insufficient.
    2. Avoid reinventing foundational lattice, discriminant module, and formed category infrastructure in this repository.
  - **Contribution Direction**: This repository consumes the preamble's basic lattice infrastructure now to implement indefinite algorithms, and will contribute completed indefinite capabilities back upstream to `research` at a later milestone.

- **NEVER BUILD THE C++ CODE**:
  **Do NOT attempt to compile, build, configure, or invoke C++ compilation toolchains** for `polyhedral_common` or any other C++ source in `references/`.
  The upstream C++ codebase in `references/polyhedral_common` is strictly a reference implementation for algorithm extraction, logic translation, and structural understanding.
  All production algorithms in this repository MUST be implemented as native SageMath / Python code delegating low-level operations to standard system backends (FLINT for exact integer linear algebra, Normaliz / cddlib / PPL for polyhedral cones, Bliss / Nauty for graph canonization, and `libgap` for finite quotient group actions).

# Existing Preamble Lattice API (Obviates Foundational Tasks)

The user's `research` repository (`dzack_research.preamble.categories.lattices`, `lattice_morphisms`, `modules`, and `forms`) **already fully implements and provides** the complete formed lattice algebra and categorical infrastructure.
**Do NOT reinvent or re-implement any of the following infrastructure; consume it directly from `dzack_research.preamble`. This completely obviates greenfield implementation for Phase 2 and early foundational tasks.**

### 1. Lattice Constructors & Objects
- `Lattices(ZZ)`: Category of integral formed quadratic lattices over $\mathbb{Z}$.
- Named constructors: `Lattices(ZZ)("U")` (hyperbolic plane), `"U(m)"`, `"An"`, `"Dn"`, `"En"`, `"E8"`, `"II_{1,9}"`, `"Leech"`, and arbitrary Gram matrices.
- Direct sums and scaling: `L1 + L2` (represented biproduct $\oplus$), `L.twist(scale)` (scaled lattice $L(a)$).

### 2. Bilinear & Quadratic Form Operations
- `L.gram_tensor()`: Symmetric $(0,2)$-tensor representing the bilinear form $b_L$.
- `L.gram_matrix(basis=None)`: Presentation of $B$ in a specified basis.
- `L.b(v, w)`: Bilinear pairing $b_L(v,w) \in \mathbb{Z}$.
- `L.q(v)`: Quadratic evaluation $q(v) = b_L(v,v) \in \mathbb{Z}$ (no $1/2$ factor).
- `L.rank()`: Free module rank $\mathrm{rk}(L)$.
- `L.signature_pair()`: Real signature $(p,q)$ over $\mathbb{R}$.
- `L.discriminant()`: Signed invariant determinant $d_\pm(b) = (-1)^{n(n-1)/2}\det G$.
- `L.is_even()`, `L.is_nondegenerate()`, `L.is_unimodular()`, `L.is_finite_rank()`.
- `L.level()`: Least $N > 0$ annihilating the discriminant form.
- `L.genus()`: Genus invariant object from signature and discriminant quadratic form.
- `L.is_locally_isometric(other, prime)`: Local $p$-adic isometry check over $\mathbb{Z}_p$.

### 3. Duals, Discriminant Modules, and Finite Quadratic Forms
- `L.dual_module()`: Algebraic dual $\mathrm{Hom}_{\mathbb{Z}}(L, \mathbb{Z})$.
- `L.dual_lattice()`, `L.metric_dual()`: Metric dual $L^\# \subset L \otimes \mathbb{Q}$.
- `L.correlation_morphism()`, `L.correlation()`, `L.metric_map()`: Canonical correlation $L \to L^\#$, $v \mapsto b(v,-)$.
- `L.discriminant_module()`: Finite formed quotient module $A_L = L^\#/L$ equipped with $q_{A_L}: A_L \to \mathbb{Q}/2\mathbb{Z}$ (or $b_{A_L}: A_L \times A_L \to \mathbb{Q}/\mathbb{Z}$).
- `L.discriminant_projection()`: Quotient morphism $\pi: L^\# \twoheadrightarrow A_L$.
- `L.discriminant_class(w)`: Projection of $w \in L^\#$ to its class in $A_L$.
- `L.divided_discriminant_class(v)`: Class $[v/\mathrm{div}(v)] \in A_L$ for $v \in L$.

### 4. Lattice Elements & Roots
- `v.b(w)`: Bilinear pairing $b_L(v,w)$.
- `v.q()`, `v.norm()`: Quadratic form value $q(v) = b_L(v,v)$.
- `v.div()`: Integer divisibility $\mathrm{div}(v) = \gcd(b(v,L))$.
- `v.divided_discriminant_class()`: Class $[v/\mathrm{div}(v)] \in A_L$.
- `v.is_root()`: Whether reflection $s_v(x) = x - \frac{2b(x,v)}{q(v)}v$ is an integral lattice isometry.
- `v.to_list()`, `v.to_tuple()`, `v.to_vector()`: Coordinate representations.

### 5. Subobjects, Saturations, Orthogonal Complements, and Reductions
- `L.subobject_on(vectors)`: Sublattice $S \subseteq L$ returned as a first-class subobject pair $(S, \iota: S \hookrightarrow L)$.
- `S.inclusion()`: Exact inclusion morphism $\iota: S \hookrightarrow L$.
- `S.saturation()`: Primitive saturation $S^{\text{sat}} \subseteq L$ via Smith normal form of $L/S$.
- `S.is_primitive()`: Boolean check whether $S$ is primitive.
- `S.orthogonal_complement()`, `iota.orthogonal_complement()`: Orthogonal complement $S^\perp \subseteq L$.
- `L.radical()`: Radical $\mathrm{rad}(L) = L^\perp \subseteq L$.
- `L.radical_quotient()`: Non-degenerate quotient $L/\mathrm{rad}(L)$.
- `S.isotropic_reduction()`: Non-degenerate quotient formed lattice $S^\perp/S$ for isotropic $S$.

### 6. Overlattices, Primitive Embeddings, and Glue Maps
- `L.overlattice(*classes)`: Overlattice $L \hookrightarrow L'$ generated by isotropic elements of $A_L$.
- `L.even_overlattice_inclusions()`: All even overlattice inclusions from isotropic subgroups of $A_L$.
- `L.local_modification(p, *classes)`: $p$-primary isotropic modification.
- `L.glue_map(S, R)`: Nikulin glue anti-isometry $H_S \to H_R(-1)$ for primitive orthogonal decompositions $S \oplus R \subseteq L$.
- `L.embeds_in_even_unimodular(p, q)`: Nikulin primitive embeddability decision into $\mathrm{II}_{p,q}$.
- `L.embed_in_even_unimodular(p, q)`: Witness primitive embedding into $\mathrm{II}_{p,q}$.

### 7. Orthogonal Groups and Homset Invariants
- `L.Emb(M)`: Form-preserving embeddings $L \hookrightarrow M$.
- `L.Isom(M)`: Isometries $L \xrightarrow{\sim} M$.
- `L.Aut()`, `L.O()`, `L.orthogonal_group()`: Full orthogonal group $O(L)$.
- `L.SO()`, `L.special_orthogonal_group()`: Special orthogonal group $\ker(\det: O(L) \to \{\pm 1\})$.
- `L.stable_orthogonal_group()`: $\ker(\rho_L: O(L) \to O(A_L))$.
- `L.spinor_kernel_subgroup()`: Kernel of the real spinor norm sign.
- `L.positive_cone_subgroup()`: For signature $(1,n)$, subgroup preserving the future cone component.
- `L.discriminant_representation()`: Functorial homomorphism $\rho_L: O(L) \to O(A_L)$.
- `L.discriminant_image()`: Image $\rho_L(O(L)) \subseteq O(A_L)$.
- `L.discriminant_representation_is_surjective()`: Boolean check $\rho_L(O(L)) = O(A_L)$.

---

### What Actually Needs to be Built in `sage-indefinite-port`
Because foundational formed algebra is already complete, development in this repository focuses strictly on:
1. **Lorentzian perfect-domain reduction**: Traversal of the $(1,n)$ reduction complex for $O^\Omega(L)$ and $O(L)$ generators.
2. **$2U$-Eichler approximate models & Higher-Witt-index recursion**: Generating $A(L)$ and lifting covering sets to full $O(L)$ for Witt index $\geq 2$.
3. **Parabolic recursion & Isotropic orbit decomposition**: Cusp orbits, unipotent radicals, and inductive rank-$k$ isotropic plane orbits.
4. **Rational group integralization & Bliss/Nauty configuration canonization**: Double-coset transporters via libGAP and graph canonization.
