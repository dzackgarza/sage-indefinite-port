# Goal: finish v1 of the indefinite lattice port

## Task

Complete the open v1 issues of the port: #44, #15, #18, #68, #70 and #72. The ledger is #4. Its "Done when" is #64 together with #44.
v2 issues (#40–#61, #69) are retained post-v1 work. They do not gate v1. Do not work on them under this goal.

The status below was derived on 2026-11-11 by reading each issue and its comments against the source at port `9290d8d` and provider `1e58a2c`.
An issue comment, a handoff or an earlier goal is a claim. Before you trust a row, check it in the source.

## Working rule

- Name the row that each step advances. A step that advances no row is outside this goal.
- Scope comes from the issue's acceptance text and the first example in the row. A `NotImplementedError` outside that scope is not a work item.
- Missing functional code comes first. Verification rows come after it, and the slow frozen consumers last of all (#64).
- Subagents write code and do not run tests. The coordinator runs every test.
- Each fix of an observed defect starts from a red test that reproduces it.

## Missing code

| Issue | Item | First example or acceptance | Status |
|---|---|---|---|
| #72 | Seed every random source on the generator path, or record the seed, so that two processes return the same family | Seed or determinize `ind/eichler.py:106-119,194-201`, `ind/indefinite_reduction.py:108-119,187`, `ind/lorentzian_cells.py:4946-4966,5997-6005`. Two runs of `automorphism_group_generators` on T_Co, or a small indefinite lattice, give the same family | done |
| #18 | Isometry U⊕⟨−1⟩² ≅ [2 1; 1 0]⊕⟨−1⟩² | Witness carries the one Gram matrix to the other | done |
| #18 | Non-isometry of the Allcock rank-3 lattices with equal invariants | Isom is empty | done |
| #18 | Dawes §5 stabilizer and index data for 2U⊕⟨−2⟩⊕⟨−6⟩ | Dawes arXiv:2205.10601 §5 values | done |
| #18 | Orbit counts of norm-0, norm-2, plane and flag classes are invariant under a change of basis | 19_IndefiniteComp cases with a unimodular basis change | done |
| #18 | Public consumer and assertion for the flag invariant (`isotropic_induction.py:526`) | Equal on equivalent flags, separates two inequivalent flags | done |
| #18 | Public assertion for spacelike squares at positive index 1 (`marked_lattices.py:540`) | A spacelike vector class on a signature-(1,n) lattice | done |
| #68 | Negative isometry decisions through each port branch (Lorentzian, Eichler/approximate model), independent of the provider's spinor-genus prefilter | A same-genus non-isometric pair at index 1 and one at index ≥ 2, decided by the port route | index 1 done; index ≥ 2 needs a same-genus non-isometric pair with Witt index ≥ 2 |
| #44 | Marked and flag transporters and orbit splitting under a nontrivial congruence subgroup at rank ≥ 2 | Preimage of {±Id} in O(2U mod 3): orbit splitting of null planes or flags, with Γ-transporters | open |
| #44 | Graded flag-quotient action honours the requested Γ; partial flag | `isotropic_orbits.py:459` hard-codes `O()` | open |
| #44 | Oriented multi-block flags; oriented objects reachable through the marked-category API | 3U oriented complete flag | open |
| #44 | Precise unsupported-algorithm error where the object exists and no backend does | One input that reaches it | done |
| #15 | Tests for edgewalk options: full walk, `K3` norms, prescribed vertex, non-reflectivity reason, `edgewalk_isometry` | U+⟨46⟩ full walk; U+⟨4⟩ `K3` wall; U+⟨2⟩ from a vertex; U+⟨46⟩ reason; U+⟨2⟩ under a basis change | done |
| #15 | `edge_directions` at multi-component cusps (`edgewalk.py:1358-1372`, takes `omitted[0]` without proof); tied candidates; unverified corner vertices; raw `int(entry)` (`edgewalk.py:1073`, `edgewalk_rank2.py:48`) | Red reproducer for each observed defect, then the source-faithful fix | done |
| #15 | Loop defects: the `scaled_max ≤ 0` guard and the `_upper_bound` doubling (`lorentzian_cells.py:5789`) | Reproducer or termination proof | done |
| #70 | Regression: line-stabilizer generators on the Enriques N lattice against the exact stabilizer | Generated group equals the stabilizer | done |
| #4 | Generators of O⁺(L) and O^Ω(L) for the consumer, including O^Ω on T_Co; the O⁺ (2,2) reflection check (`prov/lattices.py:4131`) | Check the owning issue; then a small example of each | done |

## Verification

| Issue | Item | Status |
|---|---|---|
| #44 | Frozen rank-1 and rank-2 consumers at the final port/provider pair | open |
| #70 | I_1_1 and I_1_2 discriminant indices, O_plus_N two lines and two planes, sage-extensions #294, at the final pair | open |
| #68 | Conway–Sloane and positive witnesses at the final pair | open |
| #15 | `IsotropicCases_55` within its gate | open |
| #64 | Remaining families at one port/provider pair | open |
