# General isotropic source parity

## Task

Complete the pinned indefinite algorithms over arbitrary isotropic sublattice ranks and flag dimension sequences. Repair computations that prevent required answers, including stabilizer-family conversion and the existing A7 facet-orbit consumer. Use the provider revision in pyproject.toml. The additional mathematical capabilities in issue #69 remain retained requirements; they are not established by isotropic source parity.

Select individual assertions with independently known answers: transporters in 3U, coprime unipotent generators, the glued rank-four quotient image, and A2/A3 facet walks. Each check exercises the changed operation directly on these bounded inputs.

## Completion condition

| Required capability | Current evidence and gap |
| --- | --- |
| Arbitrary-rank sublattice Isom and Aut | General dispatch and radical GL actions are implemented. Maximal rank-three 3U transporters and the expected maximal-parabolic image over F3 pass. |
| Partial and complete finite flags | General chain dispatch, adapted bases, and radical flag actions are implemented. Transporters and single-orbit checks for ranks (1,3) and (1,2,3) in 3U pass. Public (2,3) and complete (1,2,3) marking-class checks pass. Partial and complete flag stabilizers preserve their markings and have the expected parabolic images over F3. |
| General orbit induction | The SeqDims loop, common invariant-frame lifts, common integral kernel, double cosets, and exact flag-equivalence deduplication are implemented. Maximal rank-three 3U enumeration passes. |
| Proper quotient-action image | The glued rank-four lattice in tests/test_isotropic_quotient_image.py has image order two versus full quotient order four modulo three. Its double-coset step produces two inequivalent flags. Both checks pass. |
| Impossible ranks | Rank-four sublattices in 3U and impossible partial flags return empty marking sets in the small checks. |
| Fixed embeddings | Coslice Isom and Aut use every generator of a primitive isotropic source. In 2U the fixed-plane stabilizer has image order three over F3, and a transporter realizes the prescribed basis shear exactly. Higher-rank class dispatch now combines general subspace induction with pointwise-group double-coset splitting in one integral frame. The public 2U fixed-basis class assertion passes. The U plus U(2) proper-image specimen fails in the pinned provider section; sage-extensions #287 owns the lost word-membership witness. This branch remains unfinished. |
| Zero-rank markings | Moving and fixed zero sources have the unique empty embedding. Their Isom and Aut dispatch to the ambient lattice. The U class and transporter assertion passes. |
| Subgroup transporters | The public coslice Core Mor on U finds an O(U) transporter between its two null basis vectors and returns no transporter in SO(U). The fixed-vector SO(U) automorphism has determinant one. This verifies the small finite-index subgroup distinction and retains the original lattice parents. |
| Split-genus isometry | The public same-genus, non-isometric Conway-Sloane ternary pair from Chapter 15 Section 11 passes its unchanged individual assertion with the pinned provider. The native isometry route handles this specimen from #68. |
| Stabilizer computation | All three conversion sites use the supplied finite lattice action directly. _integral_generator_set applies exact product/inverse L1 reduction before constructing owned integral morphisms. Small checks retain sign and unipotent kernel actions invisible modulo two, reduce coprime unipotent powers to their primitive generator, and preserve the maximal 3U parabolic image and the proper glued quotient image. Large-family runtime adequacy remains unverified; #70 is open. |
| A7 consumer | The existing facet walk uses point_permutations_from_rows and supplies its certified action directly to facet enumeration. Facet seeding uses exact cddlib GMP lexicographic LPs; rational span/kernel calculations avoid unnecessary integral saturation and each cell caches its exact affine supporting row. Sage/GAP OnSets computes the finite facet set stabilizer, whose generators are restricted to the facet; ambient transport words retain the supplied ray generators. Small A2/A3 orbit and neighbour checks and exact orthant facets pass. Full A7 runtime remains unverified. |
| Publication | Coherent source changes are on implement/general-isotropic-parity. Completion requires the remaining small mathematical checks and the actual computation-blocking repairs, not exhaustive family reruns. |

## Next work

1. Resolve the generic section dependency tracked at https://github.com/dzackgarza/sage-extensions/issues/287, preserving the configured provider revision until its disposition is known. The small failing input and independent answer are tests/test_based_isotropic_orbits.py::test_bases_of_a_maximal_plane_in_u_plus_two_u_have_three_classes. A maximal plane in U plus U(2) has three based classes, distinguished by its divisor-two line modulo two. The provider rejects evaluation of a nontrivial double-coset section representative as undecided membership in its own generated rational group. Keep the full kernel and actual quotient action in the port caller.
2. Public maximal-rank and complete finite-chain marking-class construction have passing small examples. Reverify affected assertions only when subsequent source changes require it. Keep the source domain intact when a check fails.
3. Keep #70's large-input runtime gap and #69's additional mathematical implementations with their owning backlog. Full-basis Nielsen transformations preserve the integral kernel; selection solely by finite image does not. Generic provider defects belong to sage-extensions with exact inputs and provider revision.

## Source and runtime

- Source: src/sage_indefinite_port/indefinite/isotropic_induction.py, marked_lattices.py, recursive.py, facet_orbits.py, lorentzian_cells.py.
- Published code revision: 3466828 on implement/general-isotropic-parity. The based-class proper-image branch remains unfinished under provider #287.
- Authority: pinned CombinedAlgorithms.h:594-717 (SeqDims and flag actions), :1975-2190 (extension stabilizers and common-lift double cosets), :2192-2394 (orbit induction). MatrixGroupSimplification.h:650-751 and :1386-1487 in polyhedral_common at 1592b246 define the product/inverse complexity reduction.
- Provider: cdafca84065405cb810bba047f8a965299b03762. Put sage-extensions-port-pinned/src first on PYTHONPATH; the default container import is a newer checkout.
- pycddlib is a declared dependency. The shared Sage environment currently uses a temporary uv overlay for this dependency: inside the container, set VIRTUAL_ENV=/sage/.venv and run `uv run --active --no-project --python /sage/.venv/bin/python --with pycddlib python ...`, retaining the pinned-provider and port PYTHONPATH. This preserves the shared runtime and supplies cdd.gmp. Plain sage-docker-python does not yet supply that overlay.
- Run the selected mathematical checks with sage-docker-python and this worktree's src and root on PYTHONPATH. The tests are tests/test_general_isotropic_markings.py, test_general_isotropic_stabilizers.py, test_isotropic_radical_actions.py, test_common_isotropic_lift.py, test_isotropic_quotient_image.py, and test_integral_stabilizer_generating_set.py. Select the relevant individual assertions rather than repeatedly running these files wholesale.
- The shared container's Bliss extension is repaired and passes fresh-process canonical-form and known automorphism-order checks. Sage's libbliss-cxx dependency carries BLISS_USE_GMP; a library-only lookup corrupts bliss::Stats. The repair is owned by the Sage fork's docker/SAGE-ENVIRONMENT.md and results/sage-image-bliss-fix.md. Provider issue #277 has the standalone reproducer and repair evidence. Repeating the former crashing call before checking the installed extension is not useful.
