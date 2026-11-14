# Complexity theory and test budgets

Each file in this directory states the cost of one ported algorithm: the pinned
source routine and its translation, the branches and the mathematical conditions
that select them, the work of each branch with the theorem or bound behind each
count, worst cases, and the unit costs of the backend primitives.

`ci/test_budget_models/` applies this theory to every collected test node, and
`ci/build_test_budgets.py` writes `ci/test_budgets.json`. The pytest plugin
`ci/test_budgets.py` sets each node's pytest-timeout from that table, refuses a
collected node without an entry, and logs one `BUDGET` line per test.
`ci/plan_acceptance.py` packs Acceptance jobs by the same budgets.

## Rules for a budget

- A budget states how long a correct implementation **should** take. It is
  derived by inspecting the mathematics and the code path. It is never read off
  a measured runtime: a measured runtime absorbs the slowness a budget exists to
  expose.
- Every parametrized case has its own budget. Two lattices of one family can
  differ in cost by orders of magnitude, so no per-function maximum and no
  sibling's value.
- Cost is not a function of input size. Classify the case into the branch the
  code takes, then price that branch from the invariants that drive it:
  reflective or not, isotropic or anisotropic and the Witt index, the
  discriminant group A_L and |O(q_L)| (a small-rank lattice with a large
  discriminant group can be the most expensive case), signature conditions,
  cusp, orbit, coset and cell counts, class numbers.
- Compute an invariant of the input exactly when it can be computed without
  running the computation under test; otherwise cite the literature, or state
  the assumption by name.
- No multiplicative safety factor. It scales the slack with the estimate and
  hides regressions in long tests. The budget is the derived maximum itself.
- No minimum and no default.
- Unit costs are what the primitive costs at the size used, in one shared Sage
  session. Anchor them to the reference implementation's own CI bounds where
  those exist: polyhedral_common's `ci_20_reflective.yml` runs all 8821
  ListReflect lattices within 120 minutes. A suite total in the thousands of
  hours means a model is wrong, not that the suite is large.
- Budgets catch order-of-magnitude regressions. Differences of seconds are noise.
- A budget above 20 minutes for a test or witness this repository writes needs a
  smaller input reaching the same branch, or an optimization issue under #64.
- Work shared by several cases (the same lattice's classes, generators or
  permutation model) is charged once per session, not to every case.
- No port factor and no "interpreted Python is slow" pricing. Sage's arithmetic
  is Cython over FLINT, PARI and GAP, and a correct implementation runs its hot
  loops compiled or vectorized. Price every step at the compiled cost of the
  same work the reference implementation does; a budget above the reference
  implementation's bound means the model is wrong. Per-element interpreted work
  in a hot loop is a defect the budget exposes, not a cost it includes.
- Every unit cost is a measured primitive from `ci/test_budget_models/unit_costs.json`
  (written by `unit_costs.py`, which times primitives and never a test), or an
  operation count from the theory times measured primitives. A unit cost typed in
  by hand is a guess and is not allowed. A primitive the models need and the
  benchmark lacks is added to the benchmark first.

## How to price a step

A price is derived, not stated. Do the analysis in this order.

1. **Name what the step does mathematically.** Constructing an object,
   checking a given map, and computing over a field are different problems
   from finding a map or enumerating a set. Start from these defaults, then
   confirm them against the code:
   - Constructing any object is negligible, because construction is lazy.
   - Linear algebra over a field is negligible. So is checking that a given
     map is an isometry, and so are rational eigenspaces and invariant
     subspaces.
   - The same constructions over ZZ are module theory. Their price is the
     price of the integral algorithm used (Hermite or Smith form,
     saturation).
   - Finding a map, and every search or enumeration, is hard.
   - A procedure that is not known to terminate needs a termination bound
     for the case at hand.
   - A routine taken from the pinned source is priced by reading it.
2. **Find the loops and what bounds them.** For a search, name the set it
   searches and the pruning that cuts it down. For an upstream routine, read
   its loops in the pinned source.
3. **Express each bound through invariants of the input.** For example: the
   rank, signature and definiteness; the determinant and its primes; A_L and
   the size of O(q_L); the genus size and the theorems that apply; the Witt
   index; the counts of vectors of fixed norm in a bounded region; orbit and
   coset counts. Compute the invariant when that costs less than the step.
   Otherwise bound it with a theorem.
4. **Check that the code computes the same thing the analysis prices.** A
   different decomposition, a repeated computation or an eager construction
   changes the invariants that the code actually pays for. Such a
   difference is a defect in the code, not a fact about the step's cost.
5. **Write down the function, the invariants used, and the measured unit
   costs.** A bare time with no function behind it skips the analysis and is
   not a price.

### Worked example: isometry of E10 vector complements

The case is test_equivalence_witness_carries_a_representative_to_its_reflection.

- **Step 1.** Finding an isometry of two definite rank-9 lattices is hard.
- **Step 2.** The search is PARI qfisom, the Plesken-Souvignier backtrack.
  Its candidate set is the S vectors of norm at most the largest diagonal
  entry of the LLL-reduced Gram matrix. The search depth is the rank, and
  the stabilizer chain is bounded by log2 of the size of Aut.
- **Step 3.** The complement of v in E10 with v^2 = 2d is E8 + <2d>. On the
  undecomposed Gram matrix the reduced diagonal contains 2d. The E8 theta
  series then gives S(d) = 1 + sum over k <= d of 240 sigma_3(k), so
  S(7) = 199,921, and S grows as d^4. On the orthogonal decomposition into
  indecomposables, which is unique by Eichler, S = 242 for every d.
- **Step 4.** The port decomposes the lattice before qfauto but not before
  qfisom. The price of the code is therefore the d^4 search, while the price
  of the algorithm is S = 242. The gap is a defect in the code.

Details are in isometry-testing.md, "Pointed isometries".

### Worked example: one vertex step of the rank-3 edgewalk

- **Step 1.** Each step of the walk is integral linear algebra on rank 2
  and rank 3 lattices: edge planes, cusp equations, divisible sublattices.
  The exception is the root enumeration at each vertex, which is a search.
- **Step 2.** The pinned loop counts are in vinberg-and-edgewalk.md. They
  are one edge plane per edge and norm, and one cusp equation per parabolic
  extension.
- **Step 3.** The walk's counts are set by the number of vertices and of
  simple roots of the lattice.
- **Step 4.** A profile measured about 20 s per lattice, while the
  operations in the step 2 counts total milliseconds. Most of the profile is
  object construction and category dispatch, which is negligible by step 1.
  That is a defect in the code, and the price stays the derived one.

Details are in vinberg-and-edgewalk.md and triage.md.
