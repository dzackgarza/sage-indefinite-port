# Acceptance repair handoff

## Task

Repair the failed mathematical computations across the frozen acceptance
families, in the port or in research according to operation ownership.
The scope is isometry, boundaries, centralizers, isotropic orbits,
orthogonal groups, reflective lattices, vector orbits and vector stabilizers.
Preserve perfect-domain correctness and every frozen fixture and assertion.

## Completion condition

Every reported failing case must compute its required result. Re-run cases
whose runners were lost and cases lacking reports. Final evidence must identify
the same published port and research revisions across all families. Repairs to
shared research operations require independent mathematical checks at their
owners. A passing helper or selected example leaves the other cases open.

The reference failure set is Acceptance
[40205836899](https://github.com/dzackgarza/sage-indefinite-port/actions/runs/40205836899),
port `5c990fd35f20f194afdcd95e0bb3647118298645`, research
`f44e76f18e2e9cec43b09a846e2e74bca5fd50f6`.
The full repair condition remains open.

## Now

Resolve the memory failure of unchanged `IsotropicCases_8716`. Acceptance
[40339136782](https://github.com/dzackgarza/sage-indefinite-port/actions/runs/40339136782)
at port `b3f01e016440b405081a61a8517be536475f8ec1` and provider
`07b24a2d3c02921d2a6098c32be0a1a171ce6406` killed isotropic shards 2, 7
and 10 without JUnit. Their progress logs place the active cases at 8716,
8806 and 8714 respectively. Later killed shards place isotropic 0 at 8666,
reflective 4 at 7799, reflective 5 at 7528 and reflective 13 at 7551.
These are progress-log locators, not complete case results. An isolated
8716 run with a 3 GiB address-space limit reaches `MemoryError` in
`Tensor.pullback`, materializing the family Gram tensor for the pairing graph.

Provider `daf20e327f0c9e61b72d4e98587d0ff143b6f4aa` is now pinned. It retains the owned lattice Gram tensor and family map at
the graph adapter, shares exact edge labels, and releases temporary graphs
and family realizations. Small positive/negative transporter assertions and
unchanged ListReflect_841/964 pass. The complete 8716 assertion is still open
in `.tmp/repair-isotropic-8716-certificates.log`, under a 6 GiB address-space
limit. The prior graph-lifetime probe failed and was stopped during expensive
pytest argument rendering. The current provider releases each input graph
and canonical edge list before the next graph, then certifies the composed
index permutation against the exact colours and pairing/relation tables.
A virtual-memory limit is not a hosted resident-memory failure certificate.

This provider includes the exact GAP coset realization and its dispatch
repair at `97a43580e3`. Its nonnormal S3 action, C4 preimage, anisotropic F3
plane and split four-dimensional F2 form pass their exact action, section
and Schreier assertions. These repairs are carried to sage-extensions.
Verify the complete isotropic and Enriques consumers before declaring their
failure paths repaired.

Next work, in order:

1. Obtain complete frozen reports with research
   `07b24a2d3c02921d2a6098c32be0a1a171ce6406` and the rebuilt image.
   Its independent split quadratic forms in dimensions 2, 4 and 10,
   exponent-four hyperbolic form, and bilinear U(2) group assertions pass.
   The unchanged `ci_indef_U_U_E7` isometry assertion also passes locally
   with Bliss. These checks leave the remaining consumers open.
2. Finish the unchanged Enriques case_1 boundary consumer with provider
   97a43580e3 in `.tmp/repair-transversal-enriques.log`. The earlier run's
   I_1_2 assertion passes, but its Enriques computation was superseded by
   the repaired finite-coset realization. Both require final-revision hosted confirmation. The
   unchanged D-tower N_6 boundary counts pass locally with provider 07b24a2;
   their hosted confirmation remains open.
3. Finish the A7 + diag(1,1,-1) discriminant-image assertion. Its initial
   perfect configuration has 49,861 vectors; the old pairing block retained
   51,057,664 exact scalars, and the hosted shard exited 137. Pairings now
   stream by row, and the hash counts the pinned upper-triangular multiset
   directly. Streamed columns use FLINT integer-matrix multiplication.
   The current A7 computation passed the initial hash and is computing the
   exact cone with PPL. The installed exact Normaliz backend is being probed
   on the same initial rays in `.tmp/repair-a7-normaliz.log`.
   The exact x=1 face of x^2-y^2 and a large integral basis change
   verify the invariant. The full A7 consumer is still open.
4. Read all completed artifacts from the newest Acceptance run and repair
   the remaining reflective, boundary, group, orbit and centralizer cases.
   Check JUnit and revision metadata, not job status alone. The Vinberg
   integer-quotient correction passes the unchanged ListReflect_549, 231,
   516, 241 and 880 cases. The exact pairing-label comparison repair also
   passes ListReflect_964 and ListReflect_841. Sage graph equality requires
   weighted graphs to compare labels; otherwise a canonical graph match can
   admit a nonisometry. The independent A2 pairing obstructions pass at the
   provider owner. Other reflective cases remain open. Immediate
   failure reporting now preserves tracebacks before a killed shard exits.
5. Obtain full frozen Acceptance evidence at one published port/provider
   pair, and carry any further provider repairs to the current library owner.

Every isometry, perfect-domain and vector-stabilizer shard has a complete
passing JUnit report at port `b3f01e0` / research `07b24a2` in
[40339136782](https://github.com/dzackgarza/sage-indefinite-port/actions/runs/40339136782).
That does not validate subsequent provider changes. Isotropic shards 8 and
9 and vector-orbit shard 0 also have passing reports at this pair.

## Provider and scope boundaries

The isolated provider checkout is `.tmp/research-acceptance-repair`, branch
`fix/acceptance-quotient-cardinality`. Repairs through 2e204df5 were merged
into research main at `fdfed1d61c4b2530abb9822c93e457c1d26bc825`.
The newer finite-form and Bliss repairs are published on the provider branch
and carried to sage-extensions main. The replacement library cannot yet
replace the pinned provider. Its selected rational value module is repaired
at 17a8390. The hyperbolic-plane invariants pass at 180c1e5.
Presentation-square and lattice-isometry arrow construction still fail, tracked in
[sage-extensions #1](https://github.com/dzackgarza/sage-extensions/issues/1).
The independent checks are committed at that owner in
`tests/test_finite_form_groups.py`. Research main has unrelated merge
conflicts from the library move. Preserve that work.

Polarized extension groups must restrict by O(L^f,h). Generic arithmetic,
finite configurations and finite-form groups belong to the library; the
indefinite searches and specialized group assembly belong to this port.
Preserve unrelated port work, including positive-vector validation and CVP.
Post-v1 requirements remain in TODO.md and issues #40/#41.
