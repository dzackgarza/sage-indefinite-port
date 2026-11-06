# General isotropic source parity

## Task

Complete the pinned indefinite algorithms over arbitrary isotropic sublattice ranks and flag dimension sequences. Repair computations that prevent required answers, including stabilizer-family conversion and the existing A7 facet-orbit consumer. Use the provider revision in pyproject.toml. The additional mathematical capabilities in issue #69 remain retained requirements; they are not established by isotropic source parity.

## Completion condition

| Required capability | Current evidence and gap |
| --- | --- |
| Arbitrary-rank sublattice Isom and Aut | General dispatch and radical GL actions are implemented. Maximal rank-three 3U transporters and the expected maximal-parabolic image over F3 pass. |
| Partial and complete finite flags | General chain dispatch, adapted bases, and radical flag actions are implemented. Transporters and single-orbit checks for ranks (1,3) and (1,2,3) in 3U pass. The public (2,3) marking-class check needs its result. |
| General orbit induction | The SeqDims loop, common invariant-frame lifts, common integral kernel, double cosets, and exact flag-equivalence deduplication are implemented. Maximal rank-three 3U enumeration passes. |
| Proper quotient-action image | The glued rank-four lattice in tests/test_isotropic_quotient_image.py has image order two versus full quotient order four modulo three. Its double-coset step produces two inequivalent flags. Both checks pass. |
| Impossible ranks | Rank-four sublattices in 3U and impossible partial flags return empty marking sets in the small checks. |
| Stabilizer computation | _integral_generator_set removes exact identities and duplicate full basis actions before constructing integral lattice isometries at all three conversion sites. Small checks retain sign and unipotent kernel actions invisible modulo two. Large-family runtime adequacy remains unverified; #70 is open. |
| A7 consumer | The existing facet walk now uses point_permutations_from_rows and supplies its certified action directly to facet enumeration. Small A2 checks preserve facet orbits, transport words, adjacent cells, and traversal generators. The unchanged A7 discriminant-image assertion needs a result at ed17907 or a later source revision. |
| Publication | Coherent source changes are on implement/general-isotropic-parity. Completion requires the remaining small mathematical checks and the actual computation-blocking repairs, not exhaustive family reruns. |

## Next work

1. Resolve the public (2,3) marking-class check, then check public marking-class construction for arbitrary sublattice ranks and longer finite-chain shapes. Complete the flag-stabilizer mathematical examples. Keep the source domain intact when a check fails.
2. Obtain the A7 discriminant-image result at the repaired revision. Diagnose its actual active stage before changing it. The old element-by-element ray-action path is superseded by the existing basis-action implementation; do not implement a new facet walk.
3. Measure the stabilizer-family reduction on the required large consumer. Full-basis equality preserves the integral kernel; selection solely by finite image does not. Exact duplicate removal is implemented, but runtime sufficiency for #70 is still unproven. Generic provider defects belong to sage-extensions with exact inputs and provider revision.
4. Reverify the chosen small assertions at the final published source/provider pair. Preserve the broader retained backlog and avoid claiming it complete from these specimens.

## Source and runtime

- Source: src/sage_indefinite_port/indefinite/isotropic_induction.py, marked_lattices.py, recursive.py, facet_orbits.py, lorentzian_cells.py.
- Authority: pinned CombinedAlgorithms.h:594-717 (SeqDims and flag actions), :1975-2190 (extension stabilizers and common-lift double cosets), :2192-2394 (orbit induction). The exact-map set stage is MatrixGroupSimplification.h:798-812 in polyhedral_common at 1592b246.
- Provider: cdafca84065405cb810bba047f8a965299b03762. Put sage-extensions-port-pinned/src first on PYTHONPATH; the default container import is a newer checkout.
- Run the selected mathematical checks with sage-docker-python and this worktree's src and root on PYTHONPATH. The tests are tests/test_general_isotropic_markings.py, test_general_isotropic_stabilizers.py, test_isotropic_radical_actions.py, test_common_isotropic_lift.py, test_isotropic_quotient_image.py, and test_integral_stabilizer_generating_set.py. Select the relevant individual assertions rather than repeatedly running these files wholesale.
- The shared container's Bliss extension is repaired and passes fresh-process canonical-form and known automorphism-order checks. Sage's libbliss-cxx dependency carries BLISS_USE_GMP; a library-only lookup corrupts bliss::Stats. The repair is owned by the Sage fork's docker/SAGE-ENVIRONMENT.md and results/sage-image-bliss-fix.md. Provider issue #277 has the standalone reproducer and repair evidence. Repeating the former crashing call before checking the installed extension is not useful.
