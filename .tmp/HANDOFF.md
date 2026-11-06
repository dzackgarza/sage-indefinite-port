# General isotropic induction

## Task

Complete the pinned isotropic sublattice and flag algorithms over their full rank and dimension-sequence domain. Then repair computations that prevent the required answers, including the stabilizer conversion and the existing A7 facet-orbit path. Work from current main with the provider revision in pyproject.toml. The retained extensions in issue #69 remain additional scope; source parity alone does not establish their completion.

## Completion condition

- Arbitrary-rank primitive isotropic sublattices: isometry dispatch is implemented; feasible higher-rank orbit enumeration remains unfinished.
- Partial and complete finite isotropic flags: isometry and stabilizer dispatch use the provider adapted basis and radical flag action. Small transporter checks exercise ranks (1,3) and (1,2,3) in 3U. General flag orbit enumeration remains unfinished. Stabilizer generators still need independent mathematical checks.
- Orbit induction: preserve the actual quotient action and the source double-coset splitting, including a proper-image example. This is unfinished.
- Impossible-rank requests: the small rank-four request on 3U returns the empty marking set.
- Computation: reduce the stabilizer generating family before constructing integral lattice isometries; diagnose the existing A7 implementation. Both remain unfinished.
- Verification: small mathematical examples on changed branches, exact provider pin, syntax/import checks, and publication of coherent changes. The local tests do not certify the full source domain.

## Now

Translate the common invariant-lattice lift for the two groups in CombinedAlgorithms.h:2076-2190 f_double_cosets. The construction is at :455-539: ComputeInvariantSublattice_method2, LiftToFullAutomorphism (:314-332), and MapOrthogonalSublatticeGroupUsingSublattice (:496-534). Both lifted groups include ComputeRelevantKernel of the same invariant lattice. For radical rank greater than one, extension from the perpendicular is not unique; independent arbitrary lifts do not prove subgroup inclusion.

The existing isotropic_stabilizer_action supplies a rational group, its stable containing bound, and IntegralStructureAction.invariant_module(). invariant_module_frame constructs its rational formed frame. lift_invariant_frame invokes the provider witt_extension in that frame and checks preservation in both directions. invariant_frame_kernel obtains the integral pointwise-perpendicular kernel on the same module by clearing the rational form denominators, invoking the provider kernel, and transporting its maps back. Next assemble the two lifted groups with this same kernel and use their common invariant module as the stable containing bound for the integral-structure double-coset action. Verify preservation of the common module for every lifted generator before forming double cosets.

The provider tensor engine _engine_coisotropic_witt_extension is explicitly the pinned Lorentzian.g extension and denominator-elimination computation. Invoke its public witt_extension owner; preserve the module framing through morphisms. Generic module/form/kernel defects belong to sage-extensions.

subspace_extensions assembles the common-lift double-coset step; its rank-two step in 3U produces one maximal extension. isotropic_flag_bases translates the iK loop at CombinedAlgorithms.h:2192-2394 with reduced cumulative ranks, quotient isotropic classes, extension stabilizers, double-coset sections, and exact flag-equivalence deduplication. primitive_marking_classes dispatches higher-rank arrows and finite-chain shapes to this loop and constructs the returned marking diagrams.

The end-to-end 3U orbit check exits 139 in the pinned provider coloured_vector_families.py:270 canonical_form call, reached through recursive quotient orthogonal-group generation. Provider issue https://github.com/dzackgarza/sage-extensions/issues/277 records the reproducer and stack. The native cause is unproven. Capture the failing graph and repair the owning provider/backend path under its own contribution rules; do not bypass the recursive group computation. Reverify the port against its required pin or an explicitly adopted repaired pin. The general orbit obligation remains unverified. Include a proper quotient-action image specimen before claiming it complete.

The captured graph reduces to a standalone Sage failure: `from sage.all import graphs; from sage.graphs.bliss import canonical_form; canonical_form(graphs.CompleteGraph(3), certificate=True)` exits 139 without provider or port imports. Labels and loops are unnecessary. The running container has libbliss-dev and libbliss2 0.77-3+b3, with the extension resolving /lib/x86_64-linux-gnu/libbliss.so.2. Diagnose and repair the installed Sage/Bliss path; repeating the full port computation cannot narrow this failure further. The captured original input is /home/dzack/.cache/sage-port-last-bliss-input.sobj, but the standalone triangle needs no saved artifact.

## Source and runtime entry points

- Branch: implement/general-isotropic-parity. Remote source revision before this handoff: 7a2ea68.
- Source: src/sage_indefinite_port/indefinite/isotropic_induction.py, marked_lattices.py, recursive.py.
- extension_stabilizer_images translates f_stab_plane_v on the perpendicular at arbitrary radical rank. Its rank-one caller remains the current plane-orbit path.
- flag_action_data supplies adapted radical bases and allowed radical generator images. _flag_action reads a finite chain into those inputs.
- Provider pin: cdafca84065405cb810bba047f8a965299b03762. The default container import points to a newer checkout. Use the isolated sage-extensions-port-pinned worktree when checking this branch.
- Targeted command: PYTHONPATH=/home/dzack/gitclones/sage-extensions-port-pinned/src:/home/dzack/gitclones/sage-indefinite-general/src sage-docker-python -m pytest tests/test_general_isotropic_markings.py tests/test_isotropic_radical_actions.py -q
- Issue #70 input evidence: .tmp/issue70-schreier-conversion-blocker.md. Check the current provider interface before relying on its older interface inventory. Selecting generators solely by finite image can discard integral kernel or marking actions.
