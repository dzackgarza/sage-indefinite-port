# Goal: every required computation works and returns the correct mathematics

## Task

Make sage-extensions and this port compute every capability listed below, through the public Sage objects, on any input in the stated domain. The result is a working mathematical library. Correct output comes first, at any cost in time and memory.

Read `AGENTS.md`, `.tmp/HANDOFF.md` and the audit comment of the row in hand. The audits give the exact routes, gaps and prerequisites:
[#44](https://github.com/dzackgarza/sage-indefinite-port/issues/44#issuecomment-6727208334),
[#15](https://github.com/dzackgarza/sage-indefinite-port/issues/15#issuecomment-6727208509),
[#18](https://github.com/dzackgarza/sage-indefinite-port/issues/18#issuecomment-6727208683),
[#69 Tiers 0–1](https://github.com/dzackgarza/sage-indefinite-port/issues/69#issuecomment-6727207899),
[#69 Tiers 2–4](https://github.com/dzackgarza/sage-indefinite-port/issues/69#issuecomment-6727208132).
The specification of each #69 family is its section in #69 or the body of its closed issue (`gh issue view N`).

## The working rule

Before each step, name the row below that it advances. A step that advances no row is not part of this goal. Each step ends with a public operation that computes something it could not compute before, or that now returns a correct result where it returned a wrong one.

- **Any cost.** A correct slow algorithm is accepted. Set no timeout or memory limit on a computation. Run long computations in the background and continue with the next row. A timing gate, a slow call or a large stack sample carries no information for this goal.
- **Correct means** the specified mathematical object: complete integral generators, exact witnesses with inverses, correct negative decisions, the actual finite image, the original lattice parents, and every branch of the algorithm.
- **One example per branch.** Each new branch gets one small example whose answer has an independent mathematical derivation, committed with the implementation. The row's example column gives the first one. The example finishes the branch. A wrong answer or a crash is a defect to fix now, at the owning interface.
- **No stand-ins.** Bounded searches presented as complete, fabricated empty results, stubs, and an unsupported-input error on an input in the domain leave the row open. A lazy interface counts only where the family's specification calls for one.
- **Blocking repair.** Repair a callee only when a small example of the row in hand cannot finish at all or returns a wrong result. Then return to the row.
- **Missing algorithm.** Where no pinned source covers a row, use the literature route in its specification. If the specification names none, search the literature, record the route on the owning issue and implement it. Stop only when the mathematics itself is open, and state the exact open statement.

Commit and push each coherent functional change. Provider work goes to sage-extensions at its owning interface; publish it to provider `main` and update the port pin. Change the Status column in the same commit as the implementation. A row is `done` only when its public operation covers the row's whole domain and its examples return the expected answers.

## Work that waits until every row is done

Issues #64, #70 and #36. The frozen acceptance suite and the Acceptance workflow. A7, O_plus_N, Dawes and every other large case. Profiling, optimization, benchmarks, timing and call-time gates. Evidence files, JUnit records, revision-pair records, proof worktrees and receipts. Repeated audits, inventories and coverage tables. Typing, formatting and refactoring. CI state selects nothing.

## Work list, ordered by real dependency

Rows that do not name a dependency proceed independently. Work on several at once when one is running.

### Provider work in hand

| Capability | First example | Status |
|---|---|---|
| Decide finiteness of GL over a finite ring; the provider branch `fix/finite-orthogonal-cardinality` is red on `is_finite` of GL(F5²). Publish the branch and update the pin. | GL(F5²) has order 480 | done |
| Requested-subgroup Levi action for a general finite preimage (#44). Run its example once, with no bound, in the background. | 2U, I=Ze, preimage of {Id,−Id} modulo 3: pointwise kernel has parameter lattice 3 Hom(U,Z), index 9 | done |

### Missing routes of the pinned sources

| Capability | First example | Status |
|---|---|---|
| Right cosets of k-plane and k-flag stabilizers (#44; `CombinedAlgorithms.h:1565-1628,2461-2471`) | U+⟨2⟩, I=Ze: \|G/H\|=2 | done |
| Induced action on successive flag quotients (#44) | 3U complete flag: image (±1)³ of order 8 | done |
| Oriented and based marking variants (#44) | 2U oriented maximal null plane: one class; stabilizer of index 2 in the parabolic | done |
| Full edgewalk on non-reflective lattices: all vertex orbits, polyhedron isometries (#15) | U+⟨46⟩: every edge neighbour of a listed vertex is equivalent to a listed vertex | open |
| Edgewalk `OptionNorms="K3"`, W₂ walls (#15) | U+⟨4⟩: (0,0,1) is a wall in "all" mode and not in K3 mode | open |
| Edgewalk `OptionInitialVertex`, `reason_non_reflective`, isomorphism route (#15) | U+⟨2⟩ from a prescribed vertex; U+⟨46⟩ reason; U+⟨2⟩ and a basis change return a witness | open |
| Public invariant, indefinite reduction and approximate canonical form (#18) | Basis change of U+⟨2⟩ preserves invariant and canonical form | open |
| Routes with no finished example: spacelike Lorentzian vector orbits; rank ≥3 flags with proper Levi image; rank-3 discriminant image; odd-glued negative lift; SO splitting of (1,3) flags. Run each small example once, with no bound, in the background. | A1 in U+A1; U+U+U(3); 3U+⟨-6⟩; the #44 odd-glued case; (1,3) flags in 3U | open |

### Foundations that open several families

| Capability | Opens | First example | Status |
|---|---|---|---|
| Public representation L → Aut(L/mL) with actual image; congruence, level and residue subgroups, jointly with determinant and spinor (#42) | #43 and every congruence consumer | Γ(3) on 2U: actual image and index | open |
| Orbit splitting, stabilizer generators and transporters for general finite preimages, all markings (#43) | Requires the row above | Γ(3) on 2U: split a null-line orbit, with transporters | open |
| Stabilizers and integral transporters of non-isotropic and nonprimitive sublattice markings in indefinite L, through complement and glue (#47) | #46, #56, #61 | A1 in U+A1: setwise and pointwise stabilizers | open |
| Public Hermitian lattice object: order action, Hermitian form, trace form, isometries | #51, #58, #59, #60 | Eisenstein rank-two form with its rank-four trace lattice | open |

### Retained families with available inputs

| Capability | First missing piece | Status |
|---|---|---|
| Low rational Witt index witnesses and lazy O(L) (#40) | Fair GL_r(Z) witness search; exact membership; diag(2,2,-6,-6) and U+⟨2,-4⟩ changed-basis witnesses | open |
| Content, saturation and odd-form doubling (#48) | Public content parameter through the primitive part; doubling with original-discriminant semantics | open |
| Canonical representatives and duplicate-free iterators (#62) | Γ parameter, HNF encoding for submodules, provider consumer | open |
| Degenerate ambient consumers (#41) | Provider consumers on degenerate L; congruence-kernel generators | open |
| Commuting-involution centralizers (#49) | Joint eigensublattices and simultaneous gluing | open |
| Rational isotropic flag complex (#45) | Faces of every `FlagType`, face identities, stabilizer maps, incidence multiplicities | open |
| Gram-constrained embeddings (#46) | Rank-r fixed-source Isom, Aut and classes; requires #47 | open |
| Infinite-index orbit sets (#53) | Γ-parametrized orbit-set interface; pairing slices through the #47 stabilizer | open |
| Degenerate and nonprimitive equivariant markings (#61) | Nonprimitive through #48; prescribed Gram through #47 | open |
| Normalizers (#50) | Twisted transporter on cyclotomic summands; generators of N(H) | open |
| Unitary centralizers, definite components (#51) | Definite at every place, rank ≥2, through O(Tr(αh)) | open |
| Eisenstein E⁴,¹ and cubic-surface consumers (#52) | Rank-10 trace lattice, reflection generators, cusp and root orbits | open |
| Hermitian oceans and Picard groups (#59) | Gaussian four generators in trace coordinates, then the remaining source routes | open |
| Integral Hermitian parabolics and flags (#60) | n(Z,T) law, liftable Levi image, flag stabilizers | open |

### Arithmetic backends

| Capability | First missing piece | Status |
|---|---|---|
| Li–Margulis generators and isometry decisions (#55) | Explicit Theorem 2 constant; box searches with exact verification | open |
| Li–Margulis embedding orbits and isotropic reduction (#56) | Box, primitive filter and deduplication; requires #47 | open |
| Clifford arithmetic in ranks 4, 5, 6 and order units (#57) | (-1,3) order and its isometry action on diag(2,2,-6,-6); order-unit generators; transfer to O(L) with missing cosets | open |
| Low rational Witt index full outputs (#40) | Full O(L) generators, vector classes, Witt-one lines; recursive and boundary consumers | open |
| Hermitian groups over nonmaximal orders and commutants (#58) | Commutant linear algebra; full unitary groups | open |
| Unitary centralizers, indefinite components (#51) | Generators of indefinite Hermitian unitary groups | open |
| Grunewald–Segal groups, tensor orbit–stabilizer, arbitrary-isometry centralizers (#54) | Explicit Q-algebraic groups, then generation and orbit decisions | open |
| Full domains of #50, #53, #61 that consume #54 | After the row above | open |

## Completion

Every Status cell reads `done`. Then, and only then, the work listed under "Work that waits" begins.
