# Port design: mathematical operations and implementation routes

Research owns mathematical objects and interfaces, and independent, reasonably primitive computations.
This package implements the specialized indefinite-lattice algorithms that research imports.
A public preamble method and the algorithm answering it can have different owners.

Research implements independent basic operations on quotients, linear-map fibres, exact solving, supplied group actions and partial isometries under their stated hypotheses.
The port implements particular indefinite reduction, recursion and traversal algorithms that compose these operations.
A primitive remains general mathematics when it occurs in the reference source or is applied to an indefinite lattice.
Conversely, expressing an entire specialized orbit algorithm as a preamble method does not make research responsible for implementing that algorithm.
For each boundary, identify the general operation, the reusable computation, and any specific port implementation route separately.

## Assigning a computation

A request to research must satisfy both conditions:

1. Its implementation can be specified and carried out without understanding the specialized algorithms being extracted from this port's upstream.

2. It is a reasonably primitive mathematical construction, with explicit inputs and hypotheses, that can be composed and used independently in other problems.

A method name, a standard name for its result, or a possible generalization does not satisfy these conditions.
Record the given data, requested output and independent computational route.
Research implements the primitive; this port supplies the algorithm-specific choice of data, search, completeness argument and assembly.
These are computational assignments, not new semantic definitions.

In particular, a generic group image or finite-image section consumes supplied generators and their images.
It does not discover generators of an arithmetic parabolic from its defining predicate.
A general linear solver consumes supplied equations; forming a specialized extension problem and establishing its parameterization remain with the algorithm unless an independently justified primitive construction already does so.

Existing public preamble names may dispatch to the port.
Retaining a frozen name does not require research to implement the algorithm behind it.
Additional callback signatures require the existing #33 contract procedure; this document changes no frozen consumer signature or assertion.

For every group computation, distinguish the given data from the requested output.
A subgroup predicate does not supply generators.
Schreier computation on a supplied finite action does not construct that action, prove its orbit finite, find an arithmetic stabilizer, or find generators of its ambient orthogonal group.

Ownership labels below describe responsibility, not current runtime availability.
Research issues track primitives and dispatch; port issues track the specialized computations.
The acceptance corpus in #33 is unchanged.

Notation: `L` is an integral lattice, `G = L.O()`, `A_L` its discriminant form, and `q(v) = b(v,v)`.

Marked lattices.
Every orbit, stabilizer and equivalence question here is about isomorphisms in a category of marked lattices, built from the preamble's diagram categories.
A marking is a diagram of lattice morphisms ending at L, and an isomorphism of marked lattices is a natural isomorphism of diagrams: a commutative diagram whose components are isometries.
Which components may vary decides the category:

- a vector v with q(v) = n is a morphism ι: ⟨n⟩ → L. Isometries (L, ι) → (M, κ) satisfy g ∘ ι = κ exactly, so these are objects of the coslice category ⟨n⟩/Lattices (`CosliceCategory(Lattices(ZZ), ⟨n⟩)`), and Aut((L, ι)) is the group of isometries of L fixing v;

- a sublattice I is an embedding ι: I ↪ L, an object of the arrow category Ar(Lattices) = Fun([1], Lattices) (`Cat().Mor(FiniteOrdinalCategory(2), Lattices(ZZ))`) restricted to monomorphisms.
  Its isomorphisms are squares (α, g) with g ∘ ι = κ ∘ α and α ∈ Isom(I, J), so (L, I) ≅ (M, J) exactly when some isometry carries I onto J; Aut((L, ι)) is the setwise stabilizer of I;

- a finite vector configuration S is a morphism ι from the lattice on Z^S with Gram matrix (b(s, s')), with α restricted to automorphisms of that lattice preserving the colour classes, so g(S) = T as coloured sets;

- a flag I_1 ⊂ I_2 ⊂ L is a chain of embeddings, an object of Fun([2], Lattices) (`Cat().Mor(FiniteOrdinalCategory(3), Lattices(ZZ))`); its isomorphisms are natural isomorphisms of chains.

Deciding whether two marked lattices are isometric means deciding whether the core homset between them (`Core().Mor`) is empty, and producing an element when it is not.
"Orbits of vectors of norm n under O(L)" are the isomorphism classes of pointed lattices (L, ⟨n⟩ → L) with L fixed; "orbits of isotropic planes" are those of marked sublattices.
Research #403 owns these objects and dispatches their indefinite computations to this port.
A homset definition is not an isometry, orbit or generator algorithm.

## What the port owns

- **#15:** the Lorentzian perfect-cell search, flips, side pairings and quotient traversal; Allcock edgewalk, cusp handling, rank-two edge cycles, root selection and termination.

- **#18:** Eichler approximate groups, envelopes and covering representatives; recursive orthogonal-group and marked-isometry algorithms; the reference algorithm's composition for integral parabolic and flag stabilizers; isotropic orbit induction.
  Supply indefinite group/orbit computations to research's arithmetic subgroup and centralizer operations and retain their complete consumer acceptance.
  Independent extension and kernel primitives and finite-action computations from supplied data remain research operations.

- **#36:** type correctness of the resulting implementation and its dependency contracts.

The source algorithms are recorded in [the translation map](source-to-sage-translation-map.md).
Their use of linear algebra and finite groups does not transfer the surrounding algorithm to research.
Definite orthogonal groups and existing generic arithmetic operations remain dependency calls.

## Implementation assignments

Issue references #396–#404 in this section denote research issues; #15, #18 and #33 denote port issues.

Research paths below are relative to `src/dzack_research/preamble/categories`; port paths are relative to `src/sage_indefinite_port`. These rows assign concrete implementation work.
Existing method names do not certify that the computation works.

| Computation | Implement in research | Implement in this port |
| --- | --- | --- |
| Indefinite plain/pointed Isom and Aut | `lattice_morphisms.py` and marked-diagram homsets: construct the marking and call `indefinite.marked_lattices.isometry(category, source, target)` or `automorphism_group_generators(category, marked)` (#403) | `indefinite/recursive.py`: `EichlerRecursion.orthogonal_group_generators`, `isometry`, `pointed_isometry`, `pointed_automorphism_generators`; `indefinite/marked_lattices.py` selects the algorithm (#18) |
| Isotropic line/plane/flag stabilizer generators | Construct the marked diagram and dispatch through the same marked Aut entry (#403); compute the section/kernel primitives below | `indefinite/isotropic_induction.py`: `line_automorphism_generators`, `plane_automorphism_generators`, `flag_automorphism_generators` implement the reference recursion (#18) |
| Primitive vector, isotropic line/plane/flag classes and witnesses | `isotropic_orbits.py` and marked homsets call `primitive_marking_classes(category, source, lattice)` and `isometry(...)`; replace imports of the absent `indefinite.isotropic_flags` module (#403) | `marked_lattices.py`, `EichlerRecursion` and `IsotropicInduction` perform the actual classification/witness algorithms (#18) |
| Isometry-extension fibres and their integral members | `lattices.py`: `R.rational_lifts(target, isometry)` and the returned torsor's `integral_members()`; preserve and complete the reusable computation (#401) | Call those operations from `EichlerRecursion.pointed_isometry`; supply the recursively found reduced isometry (#18) |
| Isometries fixing I-perp pointwise | `lattices.py`: `R.pointwise_perpendicular_kernel()` computes ker(r) for restriction r: Stab_O(L)(I) → Aut(I-perp). Replace its import of the absent port `indefinite.isotropic_lifts` module with the research-owned computation (#398/#401) | Consume this kernel when the recursive assembly needs ker(r). Do not create a replacement port kernel module (#18) |
| Joint isotropic/quotient action and its kernel/image/sections | Retain the maps; compute images and finite-image sections from supplied generating families, and independently justified equation fibres/kernel realizations (#398/#401) | Construct the arithmetic compatibility action and liftable image from recursive quotient groups; find the required kernel and lift generators using those primitives and assemble the stabilizer (#18). Research must not discover this image from an unanswered parabolic predicate |
| Involution centralizers | Invariant decomposition; gluing compatibility for supplied factor generators/actions; finite-index extension (#397) | Supply indefinite component O(L) generators through `EichlerRecursion` (#18) |
| Non-scalar finite-order component centralizers | Public centralizer objects, polynomial-kernel decomposition and gluing of supplied component groups; an independently established backend route only under its documented hypotheses | #18 retains the full consumer obligation and must identify the missing computational route. Neither ordinary O(L) generators nor cyclotomic decomposition computes a non-scalar centralizer; research has no blanket obligation to invent that algorithm |
| Finite-action group computations and arithmetic orbit splitting | `group/groups.py`, `group/g_sets.py`, subgroup owners and `orthogonal_quotients.py`: compute from supplied generators/actions/ambient orbit data (#396) | Supply the indefinite group, orbit and stabilizer computations these operations consume (#18) |
| Integral structures for supplied rational groups | `rational_integral_stabilizers.py:IntegralStructureAction` computes the invariant module, residue action, stabilizer, witnesses and cosets from an owned generated group and lattice inclusion under the required boundedness hypotheses (#396/#401). Remove `_ported_integral_structure` and the port-specific `RationalMatrixGroup` admission requirement (#403) | Supply the rational group generated by the approximate model or recursive step, its lattice inclusion and the application-specific boundedness argument; consume the research action (#18) |
| Eichler approximate groups and covers | Individual transvections, module/gluing operations and exact solving remain here. Move `TwoUEichlerModel.approximate_generating_family`, `covering_vector_representatives`, `covering_orbit_data`, `recursive_stabilizer_data`, `full_generating_data` and `isometry_to` out of `eichler_criterion.py` (#404). Planned `eichler_subgroup_generators` and `primitive_vector_with_invariants` names do not transfer their specialized algorithms to research | `indefinite/eichler.py`: `EichlerEnvelope`, `OrbitCoverModel.approximate_generators`, `primitive_covering_vectors`, `covering_vectors`, `splitting_vector` (#18) |
| Lorentzian perfect domains | `reduction_complexes.py` calls `perfect_domain_traversal(lattice, mode)` with a semantic lattice and consumes its semantic result; replace the current Gram-row call and dictionary decoding (#403) | `indefinite/lorentzian_cells.py`: starting-cell search, flips, cell identifications, side pairings and traversal (#15) |
| Allcock edgewalk | `hyperbolic_lattices.py` calls `edgewalk_fundamental_domain(lattice)`; replace the current Gram argument and dictionary decoding (#403) | `indefinite/edgewalk.py` and `edgewalk_extensions.py`: corners, cusp/edge steps, root selection and traversal (#15) |
| Binary edgewalk cycle | Standard quadratic-form arithmetic and semantic period isometries/vector subsets | Implement `canonical_companion`, `promised_step`, `shorter_pair`, `reduced_start_pair`, `anisotropic_cycle` here (#15); move the copied body of research `lattice_engines._binary_reduction_cycle` here. Public `L.reduction_cycle` may dispatch |
| First-success scaled shell/sphere | `close_vectors(target,bound)` and fixed-boundary enumeration | `_PositiveVectorEnumerator.first_shell` owns the multiplier loop, target scaling and bound scaling (#15). Move research `_pari_first_affine_cvp_shell`'s search here; public `first_close_vector_shell`/`sphere` names may dispatch |
| GL flags, radical extensions, definite isometry, finite configurations, exact module algebra and polyhedra | Implement and reuse these general computations through their preamble owners and established backends (#398–#401) | Invoke them; no local implementations or carrier classes |

In the current source, the line and plane generator methods delegate straight to `parabolic_subgroup().group_generators()`, and the isotropic pointed-generator branch delegates to `pointwise_parabolic_subgroup().group_generators()`. Those are unfinished port algorithm bodies (#18), not requests for research to implement the whole recursive stabilizer.
The general section, fibre and kernel calls in the adjacent rows are research work.
Research #403 repairs the opposite direction: old imports and old argument/result shapes must reach the actual port entry points.

`IntegralStructureAction(G, i)` is the research interface, with G an owned group with supplied rational generators and i: L ↪ Res(V). Its general computation forms the orbit-generated module M under a justified finite stable containing bound, then the finite residue action on M/dM with dM contained in L. Research owns this computation, including the finite stabilizer and preimage.
The port supplies G, i and the hypotheses established by its particular algorithm.
The source `GroupAction.g:813–900` performs module closure and a finite-action computation without an indefinite-form-specific step.
Therefore research #403 removes `_ported_integral_structure`'s import of `groups.integral_structures`; no replacement port callback or port-owned action class is required.
The formed `G.invariant_overlattice(L)` operation has stronger integral-pairing hypotheses than a stable underlying module; do not silently substitute one for the other.

## Operations and implementation owners

This replaces the former blanket “Planned preamble API” migration table.
Research owns the public operations in every row; the last column identifies the port's particular implementation contribution, not exclusive ownership of the operation.
Exact internal callback signatures changed by implementation must be recorded through #33 before use.

| Operation / public surface | Research supplies | Port supplies / live obligation |
| --- | --- | --- |
| `D.object(diagram)`, `Core().Mor`, chain `.arrows()` | Diagram objects, commuting squares, homsets and dispatch (#403) | Indefinite decisions, witnesses and automorphism generators (#18); no call back into the same unanswered query |
| `D.primitive_marking_classes(shape,L)` | Result objects and dispatch (#403) | Primitive vector, isotropic sublattice and flag classification (#18) |
| `Isom(...).an_element_in(H)` | Homsets, given group actions and finite-image word lifts (#396, #403) | Marked-lattice witness search and assembly with the supplied subgroup condition (#18) |
| `G.subgroup(generators)` | Generated subgroup, exact morphisms and retained supplied generators (#400) | Find the generators of indefinite groups (#15, #18). Ambient isometry membership is not membership in an arbitrary generated subgroup |
| Finite G-sets, point stabilizers, transporters, kernels, preimages, intersections, cosets and double cosets | Schreier and finite-image computations from supplied generators and explicit finite action/transversal data (#396) | Construct the relevant arithmetic action, establish its finite-index hypotheses, and lift/recombine its results (#18) |
| `G.split_orbit(marked,H)` | Generic finite-action orbit splitting when the ambient orbit and point-stabilizer data are given (#396) | Compute the indefinite ambient orbits and their stabilizers; reconstruct marked-lattice representatives and stabilizers (#18) |
| Rational-group integral structures and `IntegralStructureAction` | Invariant module closure, finite residue action, lattice stabilizer, isometry witnesses and cosets, under explicit hypotheses (#396/#401) | Supply the algorithm's group, lattice inclusion and justification that the general construction applies (#18) |
| Products, restriction, finite-index extension, induced subquotient maps | General morphism constructions; transport supplied generators and words (#397) | Choose and assemble the factors and glue-compatible subgroup for the indefinite algorithm (#18) |
| `R.parabolic_subgroup()`, pointwise subgroup, restriction and quotient actions | Subgroup objects/maps; images and word sections from supplied generators; independent formed-module equation/kernel primitives (#398, #401) | Discover the arithmetic liftable image, construct its compatibility action, and assemble the recursive stabilizer (#18). A postcomputed image of supplied parabolic generators is not a prerequisite for finding those same generators |
| `GL(M)` and a saturated module flag | General linear groups and block triangular groups from a chosen free-module flag (#398) | Select the radical/flag data and combine them with the lattice form, quotient recursion and gluing (#18). A GL flag group is not the orthogonal flag stabilizer |
| Finite coloured spanning configurations | General finite-configuration isomorphism via pairing graphs, relations and integral checks (#399) | Select the perfect-cell configurations, facets, adjacencies and marking data; traverse and identify cells (#15) |
| Isotropic and signed witnesses, content, Bezout partner | General formed-module/module operations, e.g. a retained `qfsolve` witness and image fibre (#401) | Choose witnesses and splitting vectors for recursion and envelope construction (#18) |
| Hyperbolic partners and splitting | Partner locus with existence hypotheses; splitting from supplied valid data; independently implemented general two-plane splitting (#401) | The envelope/splitting search and coverage argument of the reference algorithm (#18). Primitive does not imply an integral partner exists |
| Perpendicular extension, Witt extension, integral descent, denominator clearing, twist | General morphism operations with their exact hypotheses (#401) | Compose these into recursive isometry and stabilizer computations (#18) |
| Sections/multisections of restriction and quotient-action maps | Independent partial-isometry extension and equation solving; finite-image word lifts from supplied generators (#398, #401) | Construct the algorithm's image conditions and compatibility problem, obtain reduced generators, and assemble the stabilizer (#18) |
| `family.integral_members()` | Affine integral inverse images of supplied linear maps; independently constructed rank-one Witt/transvection lift torsors with their complete integral loci (#401) | Supply reduced isometries and any specialized higher-rank extension equations/parameterization (#18). A nonlinear torsor cannot be assigned to the affine solver by name |
| Invariant orbit span | Generated submodule and finite/bounded-action realization (#401) | Establish a finite stable containing module for the rational groups used here (#18); commensuration alone does not imply one exists |
| Definite `close_vectors`, norms, affine-line points, divisibility sublattices | Exact primitive computations (#401) | Choose and advance the multiplier, shell/sphere and flip bounds (#15) |
| Indefinite binary arithmetic | Existing exact quadratic-form backend operations, fixed-norm fibres and rational null lines (#401) | Edgewalk-oriented reduction cycle, starting pair, stopping rule and candidate enumeration (#15); indefinite reduction strategy (#18) |
| Cone facets, rays, bounded integral points | General polyhedral operations | Establish bounds and enumerate fixed-norm marks in Lorentzian cells (#15, #18). Do not request the whole specialized enumeration as a cone primitive |
| Eichler transvection, direct sum and discriminant class | Individual formed morphisms and invariants (#401, #404) | Assemble the approximate generating family, choose prescribed-invariant covering vectors and prove coverage (#18); renaming a covering routine is not a migration |
| Order of a supplied finite quotient action, finite orbit closure | Generic finite-action computation with justified finiteness (#396) | Construct the action and prove the finite orbit needed by the edgewalk (#15) |
| `G.centralizer(I)` and equivariant lattices | Centralizer objects; independent finite-group/backend computations; invariant decomposition and gluing from supplied component groups (#397, #403) | Indefinite component O(L) computations and the complete former #23 consumer obligation (#18). Non-scalar component algorithms require an identified route, not a blanket request to research |

### Sections and kernels: specify the morphism

For an isotropic inclusion I ↪ L, put P = Stab_O(L)(I) and K = I-perp/I. Distinguish restriction r: P → Aut(I-perp) from the joint action ρ: P → GL(I) × O(K). The pointwise-perpendicular subgroup is ker(r). It is not, in general, ker(ρ), nor the whole unipotent radical.
Research #398 owns these maps and independently derived formed-module kernel computations; #401 owns independent isometry-extension fibres and their rational/integral sections.
General image and word-section algorithms require supplied generators.
Discovery of the arithmetic image needed to generate P belongs to #18. The port's recursion supplies the marked lattice and reduced group computations and uses those operations.

A selected preimage of h lies in r⁻¹(h) or ρ⁻¹(h). A section on the image selects such a preimage for each h; a multisection must specify the family retained in each fibre.
State the category: a set-theoretic section need not be a group-homomorphic splitting, and a rational section need not take integral values.
For isometries between distinct lattices, use the restriction map between the relevant Isom sets and its fibre over the specified partial isometry.
Its nonempty fibre is a torsor under the corresponding restriction kernel, not automatically an affine family of matrices.
Research computes the independent constructions, not merely their objects.
The request must give either supplied effective group data or an independent equation-based construction; the word "section" does not assign the whole arithmetic lifting algorithm to research.

The inspected `IndefiniteForms.g:1060–1128` composes the separate extension routine `Lorentzian.g:115–248` with a reduced-group computation.
The extension routine solves B = XA + AᵀXᵀ and computes homogeneous integral parameters under its isotropic-complement hypotheses.
Those are reusable extension/kernel computations (#398/#401); the reference's recursive choice and assembly remain #18.

## Algorithm composition

### Lorentzian perfect-domain traversal (#15)

Find an initial perfect cell, enumerate its facets, compute its finite configuration automorphisms, flip across facet-orbit representatives, identify neighbouring configurations, and record side pairings.
The cell automorphisms and side pairings generate the component-preserving group; adjoin minus identity for the full group.
Generic configuration and finite-action computations are inputs; the search and traversal are this port's algorithm.

For isotropic markings, glue local vertex orbits along side pairings.
For timelike vector classes, enumerate primitive vectors of the prescribed square in each cell using a proved finite bound and glue their local orbits (#18). Nonzero null rays have infinitely many multiples, so an unrestricted norm-zero fibre must not be advertised as a finite cone primitive.
Spacelike classes require their own valid algorithm; a missing route remains incomplete acceptance, not a research primitive request.

### Allcock edgewalk (#15)

Construct the starting corner, traverse the edges, solve the cusp and rank-two candidate problems, select the next wall, and establish the termination/reflectivity result.
Consume general roots, pairings, cone faces, exact arithmetic, definite close-vector and binary-form operations.
The first-success multiplier search and oriented binary cycle belong to this traversal.
Vinberg is an independent computation where applicable, not an implementation substitute for the required edgewalk.

### Eichler recursion and integral structures (#18)

Construct the envelope and approximate group, using individual transvections and ordinary module operations.
Supply the rational group and lattice inclusion to research's `IntegralStructureAction`; establish the boundedness hypotheses required by its invariant-module and finite-action computations.
Select a splitting vector v, recursively compute Aut((L,v)), and find isometries to the covering representatives.
Assemble O(L) from those generators and witnesses.

For nonisotropic v, recurse on its perpendicular and consume research's integral, glue-compatible extension computation.
For isotropic v, use research's restriction and quotient-action maps: recurse for the reduced group, supply the algorithm's compatibility data, compute the admissible image with the general action operations, evaluate a section on its generators, and adjoin the corresponding kernel generators.
Research supplies the maps and independent equation/kernel primitives.
#18 constructs the liftable arithmetic image and the generator assembly.
Postcomputed Levi images and lift tables may consume the resulting generators; they cannot supply those generators' own prerequisite.
Defining `Aut((L,v))` or calling its `group_generators()` is not an implementation of this step.

### Isotropic induction, rank at most two (#18)

For each line representative I, form K = I-perp/I with research's quotient construction.
Compute P_I = Aut((L,I)) here, then its image H_I in O(K). Recurse for isotropic classes u of K and their stabilizers Q_u. Construct the required finite coset action for Q_u \ O(K) / H_I, consume its generic double-coset computation, and lift and saturate the resulting lines to planes.
Identify planes and flags by the port's marked-isometry algorithm.
Compute their stabilizers by composing the research-owned section and kernel computations above with the port's recursive reduced-group computation.
Rank greater than two remains outside the frozen acceptance scope.

### Arithmetic subgroups and centralizers (#18)

Supply the full indefinite group generators, orbit representatives, point stabilizers and character/discriminant actions to generic finite-image preimage and orbit-splitting operations.
Reconstruct the lattice markings and their subgroup stabilizers.
For involution centralizers, research #397 computes the compatible subgroup of O(L_+) × O(L_-) preserving the gluing and extends it to L. This works through general invariant-subobject, group-action and finite-index extension operations; the indefinite O(L_±) calls use #15/#18 where applicable.
Map §15 explicitly introduces centralizers as new functionality, not a port of the source recursion.
No separate centralizer implementation is assigned here on that evidence.
For further finite-order cases, research constructs the decomposition and glues supplied component centralizers.
A non-scalar component-centralizer algorithm requires an identified independent backend or specialized route with explicit inputs and hypotheses.
Its name alone assigns no computation to research.
#18 retains responsibility for identifying and integrating the route required by its consumer cases; this is not evidence that an upstream centralizer algorithm exists.
Retain every finite-order case already present in the frozen suite.
#18 remains accountable for full consumer acceptance and its indefinite dependency calls, while #397 tracks the reusable computation; neither issue is complete merely because work is routed to the other.

## Boundary and acceptance

Research imports this package at computational call sites.
The port imports research's mathematical objects and primitive operations.
These imports are valid; a runtime cycle that asks the same unresolved orbit/stabilizer question again is not.

The existing marked entry points are `isometry(D,source,target)`, `automorphism_group_generators(D,marked)` and `primitive_marking_classes(D,shape,L)` in `indefinite.marked_lattices`. The category distinguishes fixed markings from markings movable by domain automorphisms.
The Lorentzian entries are `perfect_domain_traversal(L,mode)` and `edgewalk_fundamental_domain(L)`. Research #403 owns dispatch and removal of obsolete Gram/dictionary conversions.
A research-owned public operation may select a port implementation route.
The specialized body stays port-owned; an indefinite argument alone does not transfer an independent primitive to the port.

All frozen acceptance cases remain required.
#15 owns its traversal and reflective cases; #18 owns the group, isometry, vector, isotropic, boundary and centralizer cases, including the former #19–#23 obligations.
#36 closes type debt.
An issue transferred between repositories remains unfinished until its computation and acceptance are implemented; a renamed method, removed implementation or passing xfail does not complete it.
