# Implementation ownership inventory

The current split is specified in [port design](port-design.md).
This inventory identifies the source families affected by the earlier migration and their live owners.
Deleting a port implementation in #38 did not implement its replacement.
The source at `709f81a` is recovery material, not an instruction to restore raw matrices, duplicate primitives, or retain known defects.

| Source family | General operation and reusable computation in research | Particular port implementation contribution | Issue |
| --- | --- | --- | --- |
| `groups/integral_structures.py`: `IntegralStructureAction`, `FiniteIntegralRepresentation`, integral stabilizer/equivalence/coset assembly | Module spans and finite quotients; finite actions, Schreier generators and word lifting on supplied data | Construct bounded integral models and the finite action; reconstruct lattice isometries and arithmetic cosets | #18; research #396 retains only the primitive part |
| `RationalMatrixGroup`, `ArithmeticSubgroup`, generated/predicate subgroup carriers | Group objects, exact isometry elements, supplied generators and explicit membership capabilities | Compute indefinite generators and arithmetic membership data without a duplicate group class | #18; research #400 |
| `RightCosetDecomposition`, `DoubleCosetDecomposition`, kernel/preimage/intersection helpers | General finite-action coset, kernel, preimage and intersection computations | Establish finite-index/action hypotheses and apply them to the port's lattice data | #18; research #396 |
| `canonization.CellConfiguration`, Bliss pairing graph and configuration isomorphisms | Finite coloured spanning configurations, relation and integrality checks | Perfect-cell construction, choice of markings, adjacency and traversal | #15; research #399 |
| `polyhedral` facets, `definite` wrappers | Cone facets/rays and existing definite Aut/Isom computations | None: use existing primitives | #15/#18 |
| `invariants`, `presentations`, cache helpers | Content, norm, divisor, discriminant objects, exact change of basis | Dispatch, prefilter use, reduced-presentation strategy, caching and exact deduplication | #18 |
| `lorentzian_cells` | Definite close-vector queries, signed witnesses, finite configurations/actions | Initial cell, multiplier searches, flip bounds, adjacency, side pairings and termination | #15 |
| `recursive` | Complements, quotient objects, general isometry extension and integral descent | Recursive O(L), Isom, pointed Aut/Isom, liftable subgroup selection and coverage | #18 |
| `eichler`, `OrbitCoverModel`, `EichlerEnvelope`, `TwoUEichlerModel` algorithm routines in research | Individual transvections, direct sums, ordinary Witt/partner constructions and discriminant classes | Approximate generating sets, envelope search, prescribed-invariant covering representatives and full-group assembly | #18; research #404 moves these computations here |
| `IsotropicVectorSection`, `IsometryExtensionTorsor` | Restriction maps on isometry sets; fibres, sections/multisections and integral solution loci, with their reusable computations | Choose the marked reduction and compose these operations with recursive group/isometry computations | #18; research #401 owns extension computations, not just torsor objects |
| `MatrixEquationSolution`, `solve_isotropic_extension_equation`, `CodimensionOneIsotropicExtension` | Generic linear map fibres and valid Witt extension | Form the algorithm's equations and compose chosen solutions with the retained embeddings | #18 |
| `PointwisePerpendicularKernel`, `INDEF_FORM_GetRec_IsotropicKplane` | Kernel of r: Stab_O(L)(I) → Aut(I-perp), its reusable generator computation, quotient-action kernels and sections under their own hypotheses | Recursive reduced-group computation and the reference's composition of section values and kernel generators into flag stabilizers | #18; research #398/#401 own the general computations |
| `edgewalk_rank2`: `canonical_companion`, `promised_step`, `shorter_pair`, `reduced_start_pair`, `anisotropic_cycle` | Existing binary-form arithmetic and exact root computations | Oriented cycle, candidate enumeration and stopping rule for the edgewalk | #15 |
| `edgewalk`: cusp, corner, extension-root and next-edge routines | Affine-line lattice points, roots, Coxeter diagrams, finite quotient orders | Cusp equations' geometric setup, root selection, traversal and reflectivity certificate | #15 |
| Fixed-square vectors in Lorentzian cell cones | Cone and lattice objects, bounded integer-point enumeration | Bound the norm fibre and enumerate primitive cell marks | #15/#18 |
| Arithmetic subgroup orbit splitting (former #22) | Generic finite-image preimages and splitting from supplied orbit/stabilizer data | Compute the indefinite data and assemble the marked-lattice answers | #18 |
| Centralizers (former #23; map §15 adds functionality beyond the source port) | C_G(I) and applicable computations; invariant decomposition, compatible subgroup of component centralizers, and finite-index extension | Supply indefinite component orthogonal-group computations; preserve and verify the complete consumer result | research #397/#403; #18 retains former #23 acceptance |

Acceptance remains the complete #33 corpus and suite.
None of these rows is marked implemented by this ownership audit.
The defect observations in the historical record remain investigation inputs for #15/#18 until source and runtime evidence resolve them; their old destination labels do not assign ownership.

## Source evidence

- `references/vendor/Indefinite.jl@374a5ebb/indef/lib/IndefiniteForms.g:683`: `INDEF_FORM_GetApproximateModel` constructs the envelope, integral stabilizer/cosets and transformed covers.

- The same file, line 911: `GetAutomorphismOfFlag` is the elementary block-triangular module primitive.
  Lines 1060–1128: `INDEF_FORM_GetRec_IsotropicKplane` composes a reduced-group computation with calls to `Lorentzian.g:115–248` for extension fibres and kernel generators.
  Research owns the reusable section/kernel computation, including the linear equation B = XA + AᵀXᵀ under the source's hypotheses; the port owns the recursive composition.
  Source residency does not assign both to the port.

- `references/vendor/Indefinite.jl@374a5ebb/indef/lib/GroupAction.g:813`: `MatrixIntegral_GetInvariantSpace`; lines 869 onward: integral stabilizer, right-coset and equivalence constructions.
  Their finite-action engine is reusable; construction of the arithmetic problem is owned here.

- Translation map §§7–10 and §§12–15 retain the remaining source correspondences and acceptance obligations.
  Section 15 explicitly introduces centralizers as new functionality.
  Its invariant-decomposition/gluing composition belongs with research's centralizer operation; indefinite component O(L) computations call this port.
  No separate port centralizer route is established by that source record.

## Historical deletion record

The following is the inventory used by #37/#38 at port `709f81a` and research `9bc51845d1`. Its MIGRATE/REPLACE/DROP labels describe the old decision and the cited deletion commits, not current ownership or completed implementations.
Use the table above for all new work.

<details> <summary>Original inventory and defect observations, retained for source recovery</summary>

# Migration inventory

Every function and class in `src/sage_indefinite_port`, private ones included, classified against references/port-design.md.
The last column is the commit that applied the row (#38); — marks a KEEP row the cleanup left unchanged.
Source revision: 709f81a. Preamble revision: research@9bc51845d1. This is the work list of #37.

Dispositions:

- **MIGRATE**: general mathematics the preamble lacks; moves to research under the named issue, adapted to the preamble's categories and public operations.

- **REPLACE**: duplicates an existing preamble operation, named with its location; the port calls it instead.

- **KEEP-I**: indefinite-specific logic owned here (#15, #18).

- **DROP**: dead, unused, or a wrapper with no content (usage checked by grep over the port's src and tests and research's src).

Paths: port paths are relative to `src/sage_indefinite_port/`; preamble paths are relative to `dzack_research/preamble/categories/`. Line numbers are those of the source revisions above.

## Cross-cutting decisions and findings

- **Ownership boundary.** Research holds mathematical operations of general use, stated in standard terms.
  The port holds the algorithm glue that composes them as Dutour Sikirić's code and papers specify: approximate models, covering lists, the generation step O(L) = ⟨A, Stab(v), t_i⟩ (if A-translates of finitely many r_i = t_i(v) cover G·v, then G = ⟨A, Stab_G(v), t_i⟩), envelopes, prefilters and traversal bookkeeping.
  Glue that crossed into research moves to the port: `TwoUFullGeneratingDatum`, `full_generating_data`, `covering_orbit_data` and `recursive_stabilizer_data` (`eichler_criterion.py:197-261, 870-1062`), and the further items listed in dzackgarza/research#404 (the approximate-model framing, base-changed approximate actions, the Lorentzian traversal record classes and the decoding of the port's records).

- **The circular imports end with the migration.** Research imports the port at `rational_integral_stabilizers.py:122` (integral structures, research#396) and `lattices.py:5135` (pointwise kernel, research#398). Both bodies move to research.

- **The preamble's subgroup carriers stay; the port's capability moves into them.** The port's `RationalMatrixGroup`, `ArithmeticSubgroup`, `GeneratedSubgroup` and the shells of its `*Subgroup` classes duplicate `Groups.subgroup` / `GeneratedSubgroups` (`group/groups.py:3105, 3922`) and `PredicateSubgroups` (`group/predicate_subgroups.py`). The new capability migrates into the preamble's classes as generator resolutions (research#396): Schreier generators through a finite permutation image, word lifting, finite-index `index`, orbit with transversal, and cosets and double cosets in infinite groups.
  Prerequisite (research#400): the generated-subgroup engine must carry rational, non-integral isometries; `_from_backend_row_action` coerces coefficients to ZZ (`lattice_morphisms.py:1772`).

- **Configurations generalise an existing category.** `VectorConfigurations` (`vector_configurations.py`) already canonizes the Bliss pairing graph for basis framings.
  The port's `CellConfiguration` extends that to overcomplete spanning families with colours, plus integral correction, and migrates into `VectorConfigurations` rather than beside it (research#399).

- **Python version.** `except A, B:` (PEP 758) at `integral_structures.py:256` and `recursive.py:460` parses only on Python 3.14, while pyproject declares `>=3.12,<3.15`. Both sites are error-swallowing fallbacks that are removed in this migration anyway.

- **Reference source missing.** `references/README.md` pins a `references/polyhedral_common` mirror at a55fcb7b, but the mirror is not checked out.
  Code that cites the C++ cannot be checked against it until it is.

- **New preamble gap.** The Lorentzian positive-vector enumeration needs a public "first multiplier with a nonempty affine close-vector shell, with its vectors".
  Only the private `_first_close_vector_scale` (`lattices.py:3846`) and the engine-private `first_close_vector_scale_coordinates` (`definite_lattices.py:493`) exist.
  Added to research#401.

## Correctness defects found (fixed when the code migrates or is rewritten; each needs a red reproducer first)

| Location | Defect | Owner | Basis |
| --- | --- | --- | --- |
| `vector_sections.py:333-438` `IsotropicVectorSection.rational_lift` | The lifts of ψ: v^⊥/v → v'^⊥/v' form an affine family of dimension n−2 (transvections E_{v',x}), but the torsor carries no directions, so integral transporters can be missed. | research#401, #18 | algebra checked |
| `isotropic_lifts.py:316-451` `PointwisePerpendicularKernel` | Trivial for rank-1 I, its only use (`recursive.py:496`): an isometry fixing e^⊥ sends f to f + c·e with q(f) = q(f) + 2c·b(e, f), so c = 0. The odd-lattice branch adds no unipotent generators. | research#398, #18 | algebra checked |
| `recursive.py:481-499` `_rational_stabilizer_lifts` | For isotropic v, all of O(K) is lifted rather than the gluing-preserving Levi image; the generated rational group may miss Stab(v). | #18 (via research#398, #396) | inference |
| `eichler.py:112-170` `OrbitCoverModel.covering_representatives` | Norm 0 through the envelope asks only for vectors primitive in E, so vectors primitive in L but not in E are missed. | #18 | inference |
| `invariants.py:108-124` `VectorPrefilter.discriminant_class_orbit_key` | A function of (norm, div, content) only; it does not use the class [v/div]. | research#401 | read |
| `presentations.py` `ReducedLatticePresentation` | Calls `lll_reduction`, which raises on every indefinite lattice. | drop | read |
| `lorentzian_cells.py:921-929, 891-918, 958-1000, 1017-1069` | Unbounded `while True` loops: `_search_initial_vectors` doubles forever on an anisotropic lattice in isotropic mode; `_find_positive_vectors` has no stop when `scaled_max ≤ 0`; `_upper_bound` doubles unboundedly. | #15 | inference |
| `lorentzian_cells.py:77-87` | The rank-2 U special case depends on the basis: it fires only when both basis vectors are isotropic. | #15 | inference |
| `lorentzian_cells.py:543` `perfect_domain_traversal` | Takes Gram rows and returns records; the frozen test passes a lattice and reads `.cells()`; research passes negated rows. | #15, research#403 | read |
| `edgewalk.py` `_edge_directions` (517), `_isotropic_next_vertex` (486), `_edge_step` (504), `_corner_from_roots` (130), `_class_action_order` (266) | Unhandled cases: non-simplicial cones at multi-component cusps; anisotropic edge planes should return non-reflective but raise; tied candidates dropped; unverified corner vertices; an arbitrary bound `max(12, 12N²)`. | #15 | inference |
| `edgewalk.py:598-606` | The rank-2 branch reports reflective with one wall, against the preamble's treatment of signature (1,1). | #15 | read |
| `edgewalk.py:91,125,594` | `int(entry)` truncates a rational entry silently. | #15 | inference |
| `recursive.py:409-428`, `eichler.py:181-211, 376-396` | The splitting vector is chosen only among e_i and e_i ± e_j, and the code raises when none is positive, although positive vectors exist; three copies. | #18 | read |
| `recursive.py:861` | Vector-orbit deduplication discards its transporters, against the design. | #18 | read |
| `integral_structures.py:256`, `recursive.py:458-461`, `integral_structures.py:1404-1443` | Error-swallowing fallbacks (`except … : return False / () / None`). | migration | read |
| `canonization.py:263-264`, `polyhedral.py:63` | `libgap.Group([])` when a stabilizer has no generators (observed in #15's red cases). | research#399 | observed |
| research `eichler_criterion.py:314-316` `from_represented_biproduct` | Takes `factors[-1]` as the complement and drops the middle factors when there are four or more. | research#402 | inference |
| research `eichler_criterion.py:1064-1143` `TwoUEichlerModel.isometry_to` | False non-isometry for 2U + E8(-1)^2 against 2U + D16^+(-1). | research#402 | mathematics |
| research `lattices.py:3269` against `isotropic_parabolics.py:243` | Two Eichler transvection conventions with opposite parameter sign. | research#402 | read |

## groups/integral_structures.py (1515 lines)

| Symbol | Disposition | Destination or replacement | Adaptations and defects | Commit |
| --- | --- | --- | --- | --- |
| `RowLatticeKey`, `IntegralRows`, `RationalMatrixKey` (43-45) | MIGRATE (private) | research#396 module | — | 36e286e |
| `RationalMatrixGroup` (48) and `__init__`, `rational_lattice`, `generators`, `generators_and_inverses`, `_repr_` | REPLACE | `V.Aut().subgroup(generators)` (`group/groups.py:3105`) → `GeneratedSubgroups` (`:3922`); `group_generators()`, `supergroup().domain()` | blocked until the engine carries QQ isometries (research#400); private `_isometry_from_column_matrix`, `_row_action_matrix` | 36e286e |
| `RationalMatrixGroup._from_column_matrices` (72), `._generator_column_matrices` (107) | MIGRATE (private) | research#400 engine seam of `LatticeIsometryMor` | private `_engine_component_matrix`; used by `eichler.py:83,105,479` | 36e286e |
| `_embedded_basis` (126) | DROP | private duplicate at `rational_integral_stabilizers.py:407` | — | 36e286e |
| `_underlying_coordinates` (133) | REPLACE | element `to_vector()` | private `_framing_lift`, `_engine_element` | 36e286e |
| `_span_embedding` (152) | MIGRATE | research#401: `RestrictedScalarsModules.ParentMethods.subobject_on(elements)` (Res(V) lacks spans) | HNF private | 36e286e |
| `ArithmeticSubgroup` (193), `.supergroup`; `GeneratedSubgroup` (352), `.selected_generators` | REPLACE / DROP | `GeneratedSubgroups`, `Subgroups.supergroup`, `selected_subgroup_generators()` | generators never checked to lie in the supergroup | 36e286e |
| `FinitePermutationRepresentation` (208): `__init__`, `group`, `points`, `image_order`, `_permutation_of`, `_image_subgroup` | REPLACE (carrier) | `FiniteGSets(G)(points, action)` (`group/g_sets.py:61, 672`), `permutation_representation()`, `GroupMorphism.image` | GAP handles private | 36e286e |
| `FinitePermutationRepresentation.is_faithful` (250) | MIGRATE | research#396 `is_faithful_action()` | swallows TypeError/ValueError/OverflowError → False; may not terminate on infinite groups | 36e286e |
| `._generator_word_matrix` (267), `.generator_word_lift` (281), `._gap_word_lift`, `._gap_word_matrix` (314-322) | MIGRATE | research#396 word lifting on `GeneratedSubgroups` (next to `GroupsWithChosenFreeBasis.reduced_word`, `groups.py:3640`) | duplicated word parsing; private `_row_action_matrix`, `_isometry_from_column_matrix` | 36e286e |
| `._preimage_generators` (324), `._subgroup_index` (344) | MIGRATE (core capability) | research#396: generator resolution of `PreimageSubgroups` via `GroupMorphism.preimage_subgroup` (`groups.py:1945`); `index()` | no reduction of Schreier generators | 36e286e |
| `FinitePreimageSubgroup` (367), `KernelSubgroup` (393), `StabilizerSubgroup` (414), `CentralizerSubgroup` (450), `IntersectionSubgroup` (484) and accessors | MIGRATE (generators) into the preamble shells | research#396: `PreimageSubgroups`, `KernelSubgroups`, `StabilizerSubgroups`, `CentralizerSubgroups`, `IntersectionSubgroups` (`predicate_subgroups.py:326-562`) gain generators | centralizer and intersection need a faithful finite image, so they cover finite G only; not used in src | 36e286e |
| `RightCosetDecomposition` (525), `DoubleCosetIntersection` (554), `DoubleCosetDecomposition` (577) | MIGRATE | research#396: cosets and double cosets with a finite-index side (generalising `groups.py:3183-3210` and `orthogonal_quotients.py:257`) | representatives are left cosets g·G_L despite the name | 36e286e |
| `FiniteIntegralRepresentation` (631) and `orbit`, `orbit_witnesses`, `_compute_submodule_orbit`, `_lattice_stabilizer`, `transporter_to(_key)`, `right_cosets`, `double_cosets`, key and permutation helpers | MIGRATE | research#396, merged into `FiniteCommensurabilityQuotient` (`rational_integral_stabilizers.py:238`); the orbit with transversal goes to research#396 | hand-written GF(p) row reduction (should be Sage `echelon_form`); GAP `fail` detected by string comparison; left generators not checked to lie in G; `PreImage` run twice | 36e286e |
| `.orbit_images` (686), `._evaluate_free_word` (914) | DROP | trivial wrappers | — | 36e286e |
| `IntegralStructureAction` (1056) and `rational_group`, `lattice_inclusion`, `preserves_selected_lattice` | MIGRATE / REPLACE | research#396: merged into `rational_integral_stabilizers.IntegralStructureAction` (`:129`); the port callback at `:118-126` is deleted | takes a preamble generated subgroup | 36e286e |
| `._ambient_action_matrix` (1094), `._scaling_morphism` (1257), `.finite_module`, `.finite_projection`, `.intermediate_image`, `._carries_lattice` | REPLACE | `_row_action_matrix().transpose()`; `FiniteCommensurabilityQuotient` `scaling_morphism`, `quotient_module`, `quotient_projection`, `intermediate_image` (`:286-344`); `is_in_image` | `_carries_lattice` returns False instead of raising | 36e286e |
| `._invariant_overlattice`, `.invariant_overlattice` (1140, 1205) | MIGRATE | research#401 invariant overlattice | unbounded loop if the precondition fails | 36e286e |
| `._quotient_exponent`, `.quotient_exponent`, `._restricted_matrix*`, `._invariant_basis_matrix`, `._selected_invariant_coordinates`, `._intermediate_image_key`, `._finite_representation`, `.finite_representation`, `.lattice_stabilizer`, `.right_cosets`, `.double_cosets`, `._repr_` | MIGRATE | research#396 | — | 36e286e |
| `.transporter` (1389) | MIGRATE | research#396 `transporter(target)` (`:188`) | three `except ValueError: return None` blocks | 36e286e |
| `._selected_submodule`, `.selected_submodule`, `.double_cosets_from_generators`, `integral_structure_action`, `integral_structure_action_for_group` | DROP | unused, or the callback target that the migration removes | — | 36e286e |

## backends/, groups/groupoid_cache.py, invariants.py, small modules

| Symbol | Disposition | Destination or replacement | Adaptations and defects | Commit |
| --- | --- | --- | --- | --- |
| `canonization.CellConfiguration` (37) and `__post_init__` | MIGRATE | research#399: generalise `VectorConfigurations` to spanning, coloured families (colours as a mark family, as in `polyhedral_cones.py:552`) | optional `roles=()` | 36e286e |
| `CellConfiguration.from_coordinate_rows` (55) | DROP | callers build preamble elements | raw rows | 36e286e |
| `_BlissCanonization` (70), `_bliss_canonization` (298), `_configuration_graph` (327), `_coordinate_matrix` (351), `_configuration_row_action` (361), `_verified_integral_isometry` (391), `_configuration_rational_automorphism_group` (435) | MIGRATE (private) | research#399, merged with `_canonical_pairing_data` and `_pairing_graph` (`vector_configurations.py:71, 114`) | private `_isometry_from_column_matrix`; unused parameter `canon`; dead guards | 36e286e |
| `cell_stabilizer` (176) | MIGRATE | research#399, built on research#401 (invariant overlattice) and research#396 (stabilizer of an integral structure) | inline overlattice loop and raw libGAP orbit and preimage; `libgap.Group([])` | 36e286e |
| `cell_transporter` (112) | MIGRATE | research#399 configuration transporter | private `_row_action_matrix` | 36e286e |
| `_configuration_integral_structure_action` (92) | MIGRATE | research#396 | — | 36e286e |
| `_rational_lattice_with_integral_structure` (457) | MIGRATE | research#396: `(L⊗QQ).integral_structure(basis_vectors)` | imported across modules though private | 36e286e |
| `_gram_matrix` (425) | REPLACE | `L.gram_matrix()` (`lattices.py:1874`) | — | 36e286e |
| `presentation_bucket_key` (79), `_identity_rows` (488), `_rational_rows` (492) | DROP | duplicates `lattice_prefilter` / trivial conversions | — | 36e286e |
| `polyhedral.configuration_cone` (20) | DROP | `RationalPolyhedralCones(L).from_rays(vectors)` (`polyhedral_cones.py:162`) | — | 36e286e |
| `polyhedral.configuration_facets` (28) | REPLACE | `facet_covectors()` (`:570`) with `evaluate_covector(c, v) == 0` (`:334`) | private `_engine_polyhedron` | 36e286e |
| `polyhedral.facet_orbits` (43) | REPLACE | orbits of a group on a finite set (`group/groups.py:3689`, `g_sets.py:218`), once the stabilizer is a preamble group | raw libGAP | 36e286e |
| `polyhedral._configuration_permutation` (68) | MIGRATE | research#399 action of isometries on a configuration | private `_row_action_matrix` | 36e286e |
| `definite.definite_orthogonal_group` (12), `definite.definite_isometry` (21) | REPLACE | `L.O()`, `source.Isom(target).an_element()` | wrappers | 36e286e |
| `groupoid_cache.IsometryGroupoidCache` (16-82) | DROP (unused); the operation migrates | research#397: transporting group generators along an isometry | conflates "unknown" with "non-isometric" | deb8572 |
| `invariants.LatticePrefilter`, `lattice_prefilter` (23, 33) | KEEP-I | bucketing for the port's deduplication, composed of research invariants: genus [P] where it applies, otherwise signature, parity and discriminant form [P] | private `_engine_element`, `_engine_component_matrix` → public invariants | 36e286e |
| `invariants.AttackProfile` (52), `from_lattice` | KEEP-I | sign normalization | `positive_index` holds min(n+, n−) | 36e286e |
| `invariants.vector_content` (71) | MIGRATE | research#401 `Element.content()` | duplicated at 97-100 | 36e286e |
| `invariants.VectorPrefilter` (79), `from_vector` | KEEP-I | bucketing for vector deduplication, composed of `q`, `div`, `divided_discriminant_class`, `isotropic_reduction` [P] and `content` (research#401) | the orbit-key defect above; boolean mode flag; `repr` keys | 36e286e |
| `presentations.ReducedLatticePresentation` (13), `from_lattice` | DROP | preamble `LatticeReduction`; indefinite reduction is research#400 | raises on indefinite input | deb8572 |
| `isotropic_reductions` (FlagType re-export), `recursive_provider` (all), `groups/__init__`, `indefinite/__init__` re-exports, `__init__.__version__` | DROP | unused facades, raw-Gram adapters and re-exports | `__version__` duplicates pyproject; package docstring claims libGAP delegation | deb8572 |
| `readiness` (all) | KEEP | process infrastructure | — | — |

## indefinite/lorentzian_cells.py (1069 lines)

| Symbol | Disposition | Destination or replacement | Adaptations and defects | Commit |
| --- | --- | --- | --- | --- |
| `PerfectMode` (56) | KEEP-I | — | — | — |
| `IndefiniteOrthogonalAlgorithm.__init__`, `.orthogonal_group` case 1 (76-88), `.isometry` Lorentzian branch (138-149) | KEEP-I | become the Lorentzian leaf functions `lorentzian_orthogonal_group(L)` / `lorentzian_isometry(L, M)` that `recursive.py` calls; no base class | the complex cache is ignored by `orthogonal_group` | 4b6b63e |
| `.attack_profile` (66), `.orthogonal_group` case 0, `.isometry` definite branch (132-135) | REPLACE | `AttackProfile`; `L.O()`; `Isom(...).an_element()` | — | 4b6b63e |
| `._two_u_cover_model` (92), `._higher_witt_orthogonal_group` (98), `.isometry` 2U branch (150-157), `.vector_stabilizer` (159), `.vector_transporter` (178), `.orthogonal_group` default case | DROP | unreachable behind `recursive.py`'s overrides | `vector_stabilizer` is also wrong (requires every generator to lift) | 4b6b63e |
| `TraversalObjectRecord`, `TraversalAdjacencyData`, `TraversalAdjacencyRecord`, `TraversalRecord` (218-235), `LorentzianPerfectComplex.records` (487) | DROP | the boundary returns preamble objects (research#403) | private `_row_action_matrix` | 4b6b63e |
| `LorentzianPerfectCell` (239), `LorentzianCellAdjacency` (263), `LorentzianPerfectLocalBackend` (270), `.initial_cell` (271), `.flip_across` (313) | KEEP-I | `.configuration()` becomes a research#399 configuration; facet normals come from `facet_covectors` | raw matrices → preamble objects (#15); dead guards | 4b6b63e |
| `.cell_stabilizer` (306), `.cell_transporter` (351), `.facet_orbits` (310) | REPLACE | research#399 stabilizer and transporter; cone facets [P] plus orbits on a finite set [P]/research#396 | `@cache` keeps cells alive | 4b6b63e |
| `LorentzianPerfectComplex` (357) and traversal methods, `.isometry_to` | KEEP-I | — | fail-open mode mismatch; optional `local_backend` fallback | 43056f8 |
| `.component_preserving_group` (475), `.full_orthogonal_group` (482) | KEEP-I collection | research#400 generated subgroup | — | 4b6b63e |
| `_perfect_form_hash_key` (525) | KEEP-I | upstream `f_hash` bucketing, computed with `L.b` [P] | raw Gram rows | dfc1f04 |
| `perfect_domain_traversal` (543) | KEEP-I entry point | frozen: `(L, mode) -> LorentzianPerfectDomainTraversal` | contract mismatch (defects table) | e4a2279 |
| `NoLocalMarkTheoremError` (549), `MarkedCellOrbitAlgorithm` (553) with `.local_marks`, `.transport_marks`, `.global_orbits` | KEEP-I | — | no caller or test yet | — |
| `MarkedCellOrbitAlgorithm.local_stabilizer_orbits` (568) | REPLACE | orbits on a finite set [P] | hand-written search | 4b6b63e |
| `_normalizing_sign` (676), `_normalized_lattice` (687) | REPLACE | `signature_pair()`, `twist(-1)`; merged with `AttackProfile` | module-level `@cache` | 43056f8 |
| `_coordinate_row` (667), `_same_coordinates` (693), `_ambient_element` (703), `_vector_rows` (1013) | REPLACE / DROP | `to_vector()`, `L.linear_combination` (`lattices.py:1182`) | private `_element_from_coordinates`, `_engine_element` | 43056f8 |
| `_positive_direction` (707) | MIGRATE (positive vector) + KEEP-I (greedy reduction) | research#401 `L.positive_vector()` | private `_rational_positive_vector` | 43056f8 |
| `_primitive_integral_direction` (750) | MIGRATE | research#401 `primitive_part` / `content` | duplicate of three other copies | 43056f8 |
| `_bezout_partner` (763) | DROP (unused); the operation migrates | research#401 `bezout_partner()` | one of four copies | 43056f8 |
| `_negative_direction_in_span` (1003) | MIGRATE | research#401 vector of a given sign in a subspace | — | 43056f8 |
| `_PositiveVectorEnumerator` (786) and methods, `_find_positive_vectors` (891), `_search_initial_vectors` (921), `_source_mid_value` (932), `_upper_bound` (958), `_kernel_flipping` (1017) | KEEP-I | uses `orthogonal_complement`, `twist`, `LLL` [P] and the new public first-shell operation (research#401) | private `_ExactCVPEngine`; unbounded loops; dead `only_shortest=False` branch | 43056f8 |

## indefinite/recursive.py (934 lines) and indefinite/isotropic_flags.py

| Symbol | Disposition | Destination or replacement | Adaptations and defects | Commit |
| --- | --- | --- | --- | --- |
| `IndefiniteOrthogonalAlgorithm.__init__` (49) | KEEP-I | `cached_method` on objects instead of `id()`-keyed caches | caches live for one call | 36e286e |
| `._section` (107) | DROP cache; section MIGRATES | research#401 | — | 36e286e |
| `._lattice_prefilter` (118) | REPLACE | genus comparison [P]; the migrated `isometry_prefilter` for odd lattices | — | dfc1f04 |
| `._orbit_cover_model` (129) | KEEP-I | the port's `OrbitCoverModel` | — | 36e286e |
| `._simple_indefinite_reduction_data` (141), `._simple_indefinite_reduction` (224) | MIGRATE / REPLACE | research#400 `indefinite_reduction() -> LatticeReduction` | private `_engine_component_matrix`; `id()` caches | 36e286e |
| `._retarget_isometry` (245) | REPLACE + MIGRATE | `twist_functor(-1)` (`lattices.py:635`), `base_change` (`lattice_morphisms.py:598`); QQ→ZZ descent → research#401 | integrality signalled by exceptions | 36e286e |
| `._selected_generators` (288) | REPLACE | `group_generators()` | type-sniffing | 36e286e |
| `.orthogonal_group` (296) | KEEP-I recursion; definite branch REPLACE; assembly KEEP-I as the single owner | research#400 generated subgroup | dead assignment; impossible guard | 36e286e |
| `.isometry` (339) | KEEP-I | genus [P], reduction research#400, covers [P] | bounded splitting-vector search (three copies); private API | 36e286e |
| `._integral_stabilizer_from_rational_lifts` (448) | REPLACE | research#396 `stabilizer()` | swallows TypeError/ValueError | 36e286e |
| `._rational_stabilizer_lifts` (481) | KEEP-I composition | lifts research#401; unipotent generators [P]; odd kernel research#398 | the isotropic-case concern above | 36e286e |
| `._split_unimodular_isotropic_transporter` (501) | KEEP-I | research#401 split plus research#397 direct sum | private API | 36e286e |
| `._split_unimodular_isotropic_data` (555) | MIGRATE; Bezout part REPLACE | research#401 `integral_hyperbolic_splitting()`; `rational_witt_decomposition().bezout_partner()` (`lattices.py:4484`) | private `_engine_component_matrix` | 36e286e |
| `.vector_stabilizer` (627), `.vector_transporter` (633), `.vector_orbit_representatives` (763) | KEEP-I | research#401 extension and content; research#396 transporter | discards transporters; NotImplementedError for index < 2 | 36e286e |
| entry points (891-924), six `isotropic_flags` functions | KEEP-I | frozen in #33 | — | e4a2279 |

## indefinite/eichler.py, vector_sections.py, isotropic_lifts.py

| Symbol | Disposition | Destination or replacement | Adaptations and defects | Commit |
| --- | --- | --- | --- | --- |
| `eichler.OrbitCover` (36) | REPLACE | `finite_indexed_family` from `covering_vector_representatives` | — | 36e286e |
| `OrbitCoverModel` (48): `.from_lattice`, `.lattice`, `._subgroup_column_matrices`, `.subgroup`, `.covering_representatives` | KEEP-I | the approximate model of the recursion (map §9), composed of research primitives: SL₂(Z)×SL₂(Z) on U⊕U, Eichler transvections, O(K), discriminant classes, Eichler's criterion; A_L and its cosets through research#396, restriction through research#397. Research's `TwoUEichlerModel` approximate-model and covering parts are themselves algorithm glue and move here (audit pending) | the isotropic coverage gap; private API | 36e286e |
| `.one_representative`, `.refined_by`, `EichlerOrbitCover.one_representative` | DROP | unused | — | 36e286e |
| `.choose_splitting_vector` (181, 376) | KEEP-I (one copy) | a real recursive-complexity score per the design | bounded search | 36e286e |
| `TwoHyperbolicPlaneDecomposition` (214) | REPLACE | the `TwoUEichlerModel` accessors (`:325, 334, 337`) | — | 36e286e |
| `EichlerEnvelope` (236), `build_eichler_envelope` (399), `.similarity_scale` | KEEP-I | Dutour Sikirić's hyperplane-basis step (`GetEichlerHyperplaneBasis`), composed of research primitives: isotropic vector with witness and hyperbolic partner (research#401), orthogonal complement and similarity morphisms [P] | a similarity embedding, not an overlattice; raw matrices → preamble objects | 36e286e |
| `EichlerEnvelope._embedding_matrix`, `._embedding_inverse`, `.pullback_element`, `._pullback_coordinate_row`, `.pullback_isometry`, `.stabilizer_of_original_lattice`, `.right_cosets` | DROP / REPLACE | `lift` / `is_in_image` (`module_morphisms.py:1236, 1255`); research#397, #396 | unused | 36e286e |
| `EichlerOrbitCover` (341): `.lattice`, `.subgroup`, `.covering_representatives` (primitive) | REPLACE | `TwoUEichlerModel` `lattice()`, `approximate_generating_family()`, `covering_vector_representatives` | — | 36e286e |
| `EichlerOrbitCover.covering_representatives` (non-primitive), `square_divisors` (512), `InfiniteLocusError` (31) | KEEP-I | the cover of non-primitive vectors as multiples c·w of primitive covers at norm n/c², part of the recursion | boolean mode flag → two functions | 36e286e |
| `eichler_transvection` (495) | DROP | `lattices.py:3269` | stricter (requires even) | 36e286e |
| `find_hyperbolic_pair` (520) | MIGRATE (witness) + REPLACE (partner) | research#401 `isotropic_vector_witness()`, `hyperbolic_partner()`; `rational_witt_decomposition` | bounded search; non-minimal pairing; dead parameters | 36e286e |
| `vector_sections` `NonIsotropicLiftContext`, `IsotropicLiftContext`, `.reduced_object`, `._rational_lift_context`, `VectorOrthogonalSection`, `orthogonal_section` | DROP / REPLACE | `v.orthogonal_complement()` (`lattices.py:4198`), `v.isotropic_reduction()` (`:4206`) | — | 36e286e |
| `NonIsotropicVectorSection.rational_lift` (99), `.integral_lift` (206) | MIGRATE | research#401 `rational_extension_across_perpendicular`, `integral_extension_across_perpendicular`; absorbs the private `_definite_complement_extensions` (`vector_orbits.py:257`) | private API; duplicated paths | 36e286e |
| `IsotropicVectorSection.rational_lift` (333) | MIGRATE | research#398/#401: `IsotropicReductions.rational_lift_to(target, ψ)` returning the full n−2 dimensional family | the missing-directions defect | 36e286e |
| `isotropic_lifts.MatrixEquationSolution` (36), `solve_isotropic_extension_equation` (461) | MIGRATE | research#401 linear map X ↦ XA + AᵀXᵀ with `kernel` / `solve_right` | — | 36e286e |
| `CodimensionOneIsotropicExtension` (54) | MIGRATE | research#401 `witt_extension` on isometries | private API; mathematics checked correct | 36e286e |
| `CodimensionOneIsotropicExtensionResult` (45), `NoIntegralExtensionError` (162) | DROP | return the isometry | — | 36e286e |
| `IntegralParameterCoset` (166), `IsometryExtensionTorsor` (174) and methods | MIGRATE | research#401 `affine_integral_locus(particular, directions)` | temporal coupling; endpoints only rank-checked | 36e286e |
| `PointwisePerpendicularKernel` (316) and methods, `pointwise_perpendicular_kernel` (454) | MIGRATE | research#398, as the body of `IsotropicReductions.pointwise_perpendicular_kernel` (`lattices.py:5133`) | trivial for rank 1 (defects table) | 36e286e |

## indefinite/edgewalk.py, edgewalk_rank2.py, edgewalk_extensions.py

| Symbol | Disposition | Destination or replacement | Adaptations and defects | Commit |
| --- | --- | --- | --- | --- |
| `edgewalk_rank2` `_gram`, `_zvector`, type aliases | DROP | rank-2 preamble lattices and elements | — | 6444c29 |
| `quadratic_eval`, `scalar_eval` | REPLACE | `L.q`, `L.b` | — | 6444c29 |
| `oriented_determinant` | KEEP-I | — | — | 6444c29 |
| `_primitive_integral` | MIGRATE | research#401 `primitive_part` | three copies | 6444c29 |
| `isotropic_factorization`, `primitive_isotropic_vectors`, `fixed_norm_vectors_isotropic`, `_floor_sqrt`, `canonical_companion`, `promised_step`, `shorter_pair`, `reduced_start_pair`, `oriented_complement`, `anisotropic_cycle` | MIGRATE | research#401 binary-form representations, beside the preamble's binary reduction (`lattice_morphisms.py:93-221`); the automorph as an isometry | inconsistent output at n = 0; unguarded loop; `nonintegral cycle` raise | 6444c29 |
| `positive_primitive_vector`, `first_next_vector` | DROP | unused | — | 6444c29 |
| `edgewalk_extensions.COXETER_INFINITY`, `compute_coxeter_matrix` | REPLACE | Sage `Infinity`; `CoxeterDiagrams().from_roots` (`coxeter_diagrams.py:1190`, `_coxeter_entry` `:49`) | — | 6444c29 |
| `NormedDynkinExtension`, `_exact_negative_pairing`, `_norm_ratio_allows`, `compute_possible_extensions` | KEEP-I | — | boolean mode flag; dead deduplication; unwritten field | 6444c29 |
| `_as_qq_matrix`, `_as_qq_roots` | DROP | preamble objects | — | 6444c29 |
| `edgewalk.EdgewalkRecord` | REPLACE | return `AllcockEdgewalkReport` (`hyperbolic_lattices.py:373`) | — | 6444c29 |
| `_signature_pair`, `_Vertex`, `_integer_gram`, `_q`, `_pair`, `_row_of_element`, `_matrix_rows` | REPLACE / DROP | `signature_pair()`, cones `from_rays(..., wall_roots=)`, `HyperbolicLattices(ZZ)`, `L.q`, `L.b`, isometries | `int()` truncation | 6444c29 |
| `_primitive_row` | MIGRATE | research#401 `primitive_part` | — | 6444c29 |
| `_corner_from_roots` | REPLACE | `fundamental_chamber` (`hyperbolic_lattices.py:667`) rays | vertex not certified | 6444c29 |
| `_initial_vertex` | KEEP-I (choice); internals REPLACE | public `vinberg_simple_roots` (`:500`) | private `_vinberg_search`; bounded search | 6444c29 |
| `_constraint_lattice` | REPLACE | dual lattice, scaling, intersection [P] | — | 6444c29 |
| `_class_action_order` | MIGRATE | research#396 order of an induced action | arbitrary bound | 6444c29 |
| `_full_root_orbit` | MIGRATE | research#396 orbit closure from seeds | — | 4b6b63e |
| `_rational_gcd_pair`, `_resolve_cusp_lattice_equation` (lattice part) | MIGRATE | research#401 points of an affine line in a lattice | fail-open None returns | 6444c29 |
| `_positive_vector_in_plane`, `_projection_data`, `_orientation_compare`, `_Candidate`, `_extension_roots`, `_signed_sqrt_compare`, `_candidate_compare`, `_vertex_from_root`, `_CuspCandidate`, `_cusp_roots`, `_isotropic_next_vertex`, `_edge_step`, `_edge_directions`, `_vertex_configuration`, `_update_finiteness`, `edgewalk_fundamental_domain` | KEEP-I | cones, roots, invariant lattices, `multiplicative_order` [P]; configurations research#399 | the unhandled cases in the defects table | 6444c29 |

</details>
