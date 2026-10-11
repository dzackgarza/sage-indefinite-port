# Port design: computational ownership

Research owns the mathematical objects and reusable primitives.
This package implements the specialized indefinite-lattice algorithms that research imports.
A public preamble method and the algorithm answering it can have different owners.

An operation belongs in research when it is simple, composable mathematics with independent applications and a general implementation: module kernels and quotients, exact linear solving, finite group actions, or a standard backend operation.
A specialized construction extracted from the reference algorithms belongs here.
Giving it a standard name, a category-owned method, or a predicate-defined result does not turn its implementation into a preamble prerequisite.

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

- **#18:** Eichler approximate groups, envelopes and covering representatives; the recursive full orthogonal group, lattice and pointed isometries; exact integral parabolic generators, Levi-image construction, lifts and kernels; isotropic line, plane and flag orbits and stabilizers; arithmetic orbit splitting and indefinite centralizer assembly.

- **#36:** type correctness of the resulting implementation and its dependency contracts.

The source algorithms are recorded in [the translation map](source-to-sage-translation-map.md).
Their use of linear algebra and finite groups does not transfer the surrounding algorithm to research.
Definite orthogonal groups and existing generic arithmetic operations remain dependency calls.

## Operations and implementation owners

This replaces the former blanket “Planned preamble API” migration table.
Consumer-facing names remain mathematical; the computation column determines where implementation work belongs.
Exact internal callback signatures changed by implementation must be recorded through #33 before use.

| Operation / public surface | Research supplies | Port supplies / live obligation |
| --- | --- | --- |
| `D.object(diagram)`, `Core().Mor`, chain `.arrows()` | Diagram objects, commuting squares, homsets and dispatch (#403) | Indefinite decisions, witnesses and automorphism generators (#18); no call back into the same unanswered query |
| `D.primitive_marking_classes(shape,L)` | Result objects and dispatch (#403) | Primitive vector, isotropic sublattice and flag classification (#18) |
| `Isom(...).an_element_in(H)` | Homsets, given group actions and finite-image word lifts (#396, #403) | Marked-lattice witness search and assembly with the supplied subgroup condition (#18) |
| `G.subgroup(generators)` | Generated subgroup, exact morphisms and retained supplied generators (#400) | Find the generators of indefinite groups (#15, #18). Ambient isometry membership is not membership in an arbitrary generated subgroup |
| Finite G-sets, point stabilizers, transporters, kernels, preimages, intersections, cosets and double cosets | Schreier and finite-image computations from supplied generators and explicit finite action/transversal data (#396) | Construct the relevant arithmetic action, establish its finite-index hypotheses, and lift/recombine its results (#18) |
| `G.split_orbit(marked,H)` | Generic finite-action orbit splitting when the ambient orbit and point-stabilizer data are given (#396) | Compute the indefinite ambient orbits and their stabilizers; reconstruct marked-lattice representatives and stabilizers (#18) |
| Rational-group integral structures and `IntegralStructureAction` | Spans, inclusions, quotients, induced maps; finite-action solver (#396, #401) | The `MatrixIntegral_*` construction of a suitable integral model, finite action, lattice stabilizer, isometry witness and cosets (#18); prove the needed boundedness |
| Products, restriction, finite-index extension, induced subquotient maps | General morphism constructions; transport supplied generators and words (#397) | Choose and assemble the factors and glue-compatible subgroup for the indefinite algorithm (#18) |
| `R.parabolic_subgroup()`, pointwise subgroup, Levi and isotropic actions | Represent the subgroups and maps, and evaluate maps on supplied elements (#398, #403) | Compute full integral generators, actual Levi image, liftable subgroup, integral lifts and unipotent/pointwise kernels (#18) |
| `GL(M)` and a saturated module flag | General linear groups and block triangular groups from a chosen free-module flag (#398) | Select the radical/flag data and combine them with the lattice form, quotient recursion and gluing (#18). A GL flag group is not the orthogonal flag stabilizer |
| Finite coloured spanning configurations | General finite-configuration isomorphism via pairing graphs, relations and integral checks (#399) | Select the perfect-cell configurations, facets, adjacencies and marking data; traverse and identify cells (#15) |
| Isotropic and signed witnesses, content, Bezout partner | General formed-module/module operations, e.g. a retained `qfsolve` witness and image fibre (#401) | Choose witnesses and splitting vectors for recursion and envelope construction (#18) |
| Hyperbolic partners and splitting | Partner locus with existence hypotheses; splitting from supplied valid data; independently implemented general two-plane splitting (#401) | The envelope/splitting search and coverage argument of the reference algorithm (#18). Primitive does not imply an integral partner exists |
| Perpendicular extension, Witt extension, integral descent, denominator clearing, twist | General morphism operations with their exact hypotheses (#401) | Compose these into recursive isometry and stabilizer computations (#18) |
| Isotropic reduction lifts | Reduction objects, torsors, general map fibres and exact linear/congruence solving (#401) | Specialized parametrization, integral lift selection, kernel generation and parabolic assembly required by the reference algorithm (#18). Reuse independently valid existing research primitives |
| `family.integral_members()` | Integral inverse images of a specified affine linear map, retaining parameter inclusion and full fibre (#401) | Supply the actual extension equations (#18). Isotropic isometry lifts form a torsor; evaluation is generally nonlinear |
| Invariant orbit span | Generated submodule and finite/bounded-action realization (#401) | Establish a finite stable containing module for the rational groups used here (#18); commensuration alone does not imply one exists |
| Definite `close_vectors`, norms, affine-line points, divisibility sublattices | Exact primitive computations (#401) | Choose and advance the multiplier, shell/sphere and flip bounds (#15) |
| Indefinite binary arithmetic | Existing exact quadratic-form backend operations, fixed-norm fibres and rational null lines (#401) | Edgewalk-oriented reduction cycle, starting pair, stopping rule and candidate enumeration (#15); indefinite reduction strategy (#18) |
| Cone facets, rays, bounded integral points | General polyhedral operations | Establish bounds and enumerate fixed-norm marks in Lorentzian cells (#15, #18). Do not request the whole specialized enumeration as a cone primitive |
| Eichler transvection, direct sum and discriminant class | Individual formed morphisms and invariants (#401, #404) | Assemble the approximate generating family, choose prescribed-invariant covering vectors and prove coverage (#18); renaming a covering routine is not a migration |
| Order of a supplied finite quotient action, finite orbit closure | Generic finite-action computation with justified finiteness (#396) | Construct the action and prove the finite orbit needed by the edgewalk (#15) |
| Centralizers and equivariant lattices | Centralizer objects, eigensublattices, cyclotomic decomposition, gluing objects, products and finite-action primitives (#397, #403) | Indefinite component recursion, compatible-generator assembly and extension to L (#18, carrying former #23 acceptance) |

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
Construct the finite integral-structure action needed to restrict this rational group to L; pass that supplied finite action to the general Schreier solver.
Select a splitting vector v, recursively compute Aut((L,v)), and find isometries to the covering representatives.
Assemble O(L) from those generators and witnesses.

For nonisotropic v, recurse on its perpendicular and select the integral, glue-compatible extensions.
For isotropic v, compute the integral parabolic: construct the compatibility data, recurse for the quotient group, compute the admissible Levi image through a finite action, solve its lift problems, and adjoin the complete integral kernel.
Defining `Aut((L,v))` or calling its `group_generators()` is not an implementation of this step.

### Isotropic induction, rank at most two (#18)

For each line representative I, form K = I-perp/I with research's quotient construction.
Compute P_I = Aut((L,I)) here, then its image H_I in O(K). Recurse for isotropic classes u of K and their stabilizers Q_u. Construct the required finite coset action for Q_u \ O(K) / H_I, consume its generic double-coset computation, and lift and saturate the resulting lines to planes.
Identify planes and flags by the port's marked-isometry algorithm.
Compute their stabilizers using the integral lift and kernel constructions above.
Rank greater than two remains outside the frozen acceptance scope.

### Arithmetic subgroups and centralizers (#18)

Supply the full indefinite group generators, orbit representatives, point stabilizers and character/discriminant actions to generic finite-image preimage and orbit-splitting operations.
Reconstruct the lattice markings and their subgroup stabilizers.
For involution centralizers, form eigensublattices using research primitives, recursively compute their orthogonal groups, construct and solve the finite gluing condition, and extend the resulting generators.
Retain every finite-order case already present in the frozen suite; do not replace a component centralizer computation by a predicate.

## Boundary and acceptance

Research imports this package at computational call sites.
The port imports research's mathematical objects and primitive operations.
These imports are valid; a runtime cycle that asks the same unresolved orbit/stabilizer question again is not.

The existing marked entry points are `isometry(D,source,target)`, `automorphism_group_generators(D,marked)` and `primitive_marking_classes(D,shape,L)` in `indefinite.marked_lattices`. The category distinguishes fixed markings from markings movable by domain automorphisms.
The Lorentzian entries are `perfect_domain_traversal(L,mode)` and `edgewalk_fundamental_domain(L)`. Research #403 owns dispatch and removal of obsolete Gram/dictionary conversions.
Internal specialized operations remain port code even when exposed through research methods.

All frozen acceptance cases remain required.
#15 owns its traversal and reflective cases; #18 owns the group, isometry, vector, isotropic, boundary and centralizer cases, including the former #19–#23 obligations.
#36 closes type debt.
An issue transferred between repositories remains unfinished until its computation and acceptance are implemented; a renamed method, removed implementation or passing xfail does not complete it.
