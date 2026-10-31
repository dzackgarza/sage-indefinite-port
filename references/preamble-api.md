# Existing Preamble Lattice API (Obviates Foundational Tasks)

The user's `research` repository (`dzack_research.preamble.categories.lattices`, `lattice_morphisms`, `modules`, and `forms`) owns and provides the formed lattice algebra and categorical infrastructure below.
This surface inventory does not certify every listed computational method as implemented or assign its specialized algorithm to research.
**Do NOT reinvent or re-implement any of the following infrastructure; consume it directly from `dzack_research.preamble`. This completely obviates greenfield implementation for Phase 2 and early foundational tasks.**

## Full lattice interface in `preamble-megadoc.md` (`#subsystem-lattices`)

### Construction and category entry points

- `Lattices` as the lattice category over a ring with constructor:

  - `Lattices(cls, *args)`

  - `Lattices(...)(self, data, basis=None, names=None, form=None, module_generators=None)`

- Category constructors:

  - `FiniteRankLattices`, `NondegenerateLattices`, `RationalLattices`, `EvenLattices`, `RootLattices`

- `Lattices(ZZ)` named lattice constructors:

  - `Lattices(ZZ)("U")`, `"U(m)"`, `"An"`, `"Dn"`, `"En"`, `"E8"`, `"II_{1,9}"`, `"Leech"`, and Gram-tensor constructors.

- Direct sum and scaling:

  - `L1 + L2` and `L.twist(scale)`.

### `Lattices` parent object API (full interface)

- Isometry and homset constructors:

  - `Aut`, `O`, `orthogonal_group`

  - `SO`, `special_orthogonal_group`

  - `Emb(codomain)`, `Isom(codomain)`

  - `hom(images, codomain=None)`, `identity_morphism`

  - `similarity_homset(other, scale)`, `similarity(scale, images=None, codomain=None)`

- Pairings and form evaluations:

  - `b(left, right)`, `q(vector)`, `correlation`, `correlation_morphism`

  - `gram_tensor`, `gram_matrix(basis=None)`, `discriminant`, `rank`, `signature_pair`, `level`, `genus`

- Lattice reductions and basis control:

  - `LLL()`, `BKZ(block_size=20)`, `HKZ()`

  - `lll_reduction`, `bkz_reduction`, `hkz_reduction`, `orthogonal_group`, `stable_orthogonal_group`

- Duality and discriminant data:

  - `dual_module`, `dual_lattice`, `metric_dual`, `discriminant_module`, `discriminant_projection`

  - `discriminant_group`, `discriminant_bilinear_form`, `discriminant_quadratic_form`

  - `discriminant_class`, `divided_discriminant_class`, `discriminant_representation`

  - `discriminant_image`, `discriminant_representation_is_surjective`

- Orthogonality and decomposition:

  - `radical`, `radical_quotient`, `orthogonal_complement`

  - `module_generating_set`, `module_generator`, `biproduct_factors`, `decomposition`, `decomposition_names`, `summands`

  - `is_decomposable`, `indecomposable_name`, `primitive_isotropic_subobject`

  - `local_modification(prime, *discriminant_classes)`, `overlattice(*discriminant_classes)`, `even_overlattice_inclusions`

- Arithmetical predicates:

  - `is_finite_rank`, `is_even`, `is_nondegenerate`, `is_unimodular`

  - `is_definite`, `is_positive_definite`, `is_negative_definite`, `is_p_elementary(prime)`

  - `is_locally_isometric(other, prime)`, `is_isometric(other)`, `is_similar(other, scale)`

- Discrete geometry and orbits:

  - `minimum`, `successive_minima`, `shortest_vectors`, `roots`, `roots_of_square(square)`

  - `vectors_of_square(square)`, `vectors_of_square_and_divisibility(square, divisibility)`

  - `kissing_number`, `root_sublattice`, `theta_series(precision=20, variable='q')`

- Isotropic and orbit APIs:

  - `isotropic_flag`, `isotropic_flag_orbit_representatives(rank=2)`

  - `isotropic_line_orbit_representatives`, `isotropic_plane_orbit_representatives`

  - `voronoi_cell(bound=None)`, `voronoi_relevant_vectors`

- Group and positivity APIs:

  - `stable_orthogonal_group`, `spinor_kernel_subgroup`, `positive_cone_subgroup`

  - `twist(scale)`, `reflection(root)`

  - `two_elementary_invariants`, `delta`

- Packing and metric helpers:

  - `babai(target)`, `closest_vector(target)`, `contact_polytope`, `covering_radius`

  - `packing_density`, `packing_radius`, `center_density`, `hadamard_ratio`, `hermite_invariant`

### `Lattices` element methods

- `b(other)`, `norm`, `q()`, `div()`, `divisibility_ideal`

- `divided_discriminant_class()`, `is_root()`

- `to_list`, `to_tuple`, `to_vector`, `monomial_coefficients`

### Subcategory APIs

- `FiniteRankLattices`: `is_finite_rank`

- `RationalLattices`: `fraction_field`, `is_nondegenerate`

- `EvenLattices`: `is_even` inherited by lattice parent

- `RootLattices`:

  - parent methods: `cartan_type`, `coxeter_number`, `highest_root`, `simple_reflections`, `simple_roots`, `fundamental_weights`

  - element methods: `coroot`, `height`, `is_negative_root`, `is_positive_root`

### Morphisms and homsets

- `LatticeMorphism` (`ModuleMorphism`):

  - constructor: `__init__(self, parent, images)`

- `LatticeEmbedding` (`LatticeMorphism`):

  - `discriminant_inclusion`, `is_injective`

- `LatticeIsometry` (`LatticeEmbedding`):

  - `centralizer_discriminant_image`, `cyclic_subgroup`, `determinant`, `discriminant_isometry`

  - `discriminant_morphism`, `formed_coinvariants`, `invariant_lattice`, `inverse`

  - `is_surjective`, `preserves_positive_cone`, `real_spinor_norm_sign`

- `LatticeHomset`:

  - `__init__(domain, codomain)`, `_element_constructor_(images)`

- `LatticeEmbeddingHomset`:

  - `an_element`, `even_overlattice_inclusions`, `is_empty`

- `LatticeIsometryHomset`:

  - `act`, `acting_group`, `compose`, `discriminant_image`, `discriminant_preimage`

  - `group_generators`, `identity`, `is_empty`, `isotropic_equivalence_witness`

  - `isotropic_orbit_representatives`, `isotropic_stabilizer_generators`, `number_of_group_generators`

  - `one`, `order`, `transporter`, `vector_equivalence_witness`, `vector_orbit_representatives`

  - `vector_stabilizer_generators`, `vectors_are_equivalent`

### Lattice objects and named tables

- `CoxeterDiagram`

  - `cardinality`, `connected_components`, `elliptic_subdiagrams`, `graph`, `coxeter_entry`

  - `coxeter_matrix`, `is_connected`, `is_elliptic`, `is_hyperbolic`, `is_parabolic`

  - `is_rooted`, `index_set`, `induced_subdiagram`, `preferred_positions`, `roots`

  - `root_gram_tensor`, `schlafli_tensor`, `signature_pair`, `vertex_names`

- `Genus`

  - `class_number`, `determinant`, `discriminant_form`, `excess`, `exists`

  - `local_symbol`, `level`, `mass`, `representative`, `representatives`, `signature_pair`

- `IsotropicFlag`

  - `basis`, `lattice`, `rank`, `terms`, `top`

- `OrthogonalCharacterQuotient`

  - `image`, `image_keys`, `splitting_isometries`, `stabilizer_image_keys`

  - `subgroup_image_keys`, `witness_meets_subgroup`

- `VectorPrimitiveExtension`

  - `class_of_representative`, `complement_is_definite`, `representative_of`

- registries: `register_indecomposable`, `register_indecomposable_gram`

### Helper functions from the lattice subsystem

- lattice constructors and homset factories:

  - `lattice(...)`, `lattice_embedding_homset`, `lattice_homset`, `lattice_isometry_homset`, `lattice_latex`

- lattice operators:

  - `diagonal_gram`, `orthogonal_sum`, `scale_gram_tensor`, `signature_pair_of_gram`

  - `colimit_lattice`, `discriminant_of_gram`, `signature_pair_of_gram`

- orbits, stabilizers, and equivalence:

  - `definite_complement_extensions`, `isotropic_equivalence_witness`, `isotropic_orbit_representatives`, `isotropic_stabilizer_generators`

  - `subgroup_isotropic_are_equivalent`, `subgroup_isotropic_orbit_representatives`

  - `subgroup_vector_orbit_representatives`, `subgroup_vectors_are_equivalent`

  - `transport_isotropic_object`, `vector_equivalence_witness`, `vectors_are_equivalent` (homset-level)

- reduction and algorithm helpers:

  - `babai`, `bkz_reduction`, `hkz_reduction`, `closest_vector`, `gaussian_heuristic`

  - `generator_pairings`, `gluing_route_discriminant_classes`

  - `packing_radius`, `packing_density`, `center_density`, `kissing_number`

  - `roots`, `roots_of_square`, `vectors_of_square`, `vectors_of_square_and_divisibility`

  - `voronoi_cell`, `voronoi_relevant_vectors`

  - `oscar_centralizer_discriminant_image`, `oscar_even_unimodular_primitive_embedding`, `oscar_rational_spinor_norm_sign`

### Compatibility note

- The above interface is the preamble contract.
  This repo must use it and should not introduce local fallback lattice mechanics in place of these names.

* * *

### What Actually Needs to be Built in `sage-indefinite-port`

Consume the existing foundational formed algebra.
Development in this repository focuses on:

1. **Lorentzian perfect-domain reduction**: Traversal of the $(1,n)$ reduction complex for $O^\Omega(L)$ and $O(L)$ generators.

2. **$2U$-Eichler approximate models & Higher-Witt-index recursion**: Generating $A(L)$ and lifting covering sets to full $O(L)$ for Witt index $\geq 2$.

3. **Parabolic recursion & Isotropic orbit decomposition**: Compose research's group actions, sections and kernels with recursive reduced-group computations to obtain cusp and flag stabilizers and the frozen rank-at-most-two isotropic orbits.

4. **Centralizer and arithmetic-subgroup consumers**: Supply indefinite orthogonal-group and orbit computations to research's general centralizer, gluing and subgroup operations; verify all frozen consumer cases.
