# source-to-sage-translation-map

Source-to-Sage translation map for the indefinite lattice port (received 2026-09-02, external report; citation markers of the generating tool removed). The plan card PLAN-INDEFINITE-EXTRACTION-AND-DECOMPOSITION is derived from this text.

## Source authority and versioning

Use two pinned source lines:

- `MathieuDutSik/Indefinite.jl@374a5ebb…` for the readable GAP formulation of the high-level algorithm.
- `MathieuDutSik/polyhedral_common@a55fcb7b…` for the later C++ reimplementation, including the newer isotropic-subspace double-coset machinery, memoization, and optimized Lorentzian implementation.

The Julia code itself is not an algorithmic source: it exposes five thin Julia-to-GAP entry points and provides subprocess wrappers for compiled executables. The substantive legacy implementations are in `Indefinite.jl/indef/lib/IndefiniteForms.g`, `Lorentzian.g`, and `GroupAction.g`; the current C++ is a later reimplementation and extension of that corpus.

The implementation rule should be:

- read the GAP code first for mathematical control flow;
- use the C++ to identify later refinements, exact corner cases, and the true performance-critical leaves;
- do not translate the Julia wrappers, command-line programs, temporary-file protocols, or custom matrix framework.

The module names below are proposed logical boundaries, not assumptions about the existing repository layout.

---

# 1. Semantic carriers and public API

Adopt the following carriers.

```python
L: IntegralQuadraticLattice
v, w: LatticeElement
I, J: LatticeSubobject
f: LatticeIsometry
G, H: OrthogonalArithmeticGroup | ArithmeticSubgroup
```

Use the convention

O^+(L) = ker( O(L) --> O(A_L) ),

and write O^Ω(L) for the component-preserving subgroup.

## 1.1 Lattices and elements

```python
B = L.gram_tensor()               # symmetric (0,2)-tensor
M = L.gram_matrix(basis=None)     # coordinate representation only

L.b(v, w)
L.q(v)                            # b(v,v), no factor 1/2
beta = L.metric_map()             # L -> L.linear_dual()

v.to_vector(basis=None)
v.to_covector()                   # exactly beta(v), not a transpose operation

L.linear_dual()                   # Hom_Z(L,Z)
L.dual_lattice()                  # L^vee inside L_Q
L.discriminant_module()           # A_L = L^vee/L

v.content()                       # largest d with v in dL
v.divisor()                       # positive generator of b(v,L)
v.discriminant_class()            # [v/div(v)] in A_L
v.is_isotropic()
v.is_primitive()
```

For a non-unimodular lattice, `v.to_covector()` lands in the linear dual, and `L.metric_map()` is an injective finite-cokernel map rather than an identification. The metric dual lattice and the linear dual must remain distinct objects.

## 1.2 Sublattices

```python
I = L.sublattice_from([v1, ..., vr], saturate=False)
Iprim = I.saturation()

v.sublattice()                    # Z v, not silently saturated
v.primitive_line()                # saturation of Z v

I.inclusion()
I.ambient_lattice()
I.perp()
I.is_primitive()
I.is_totally_isotropic()
```

For a primitive isotropic vector or primitive totally isotropic sublattice:

```python
R = v.isotropic_reduction()
R = I.isotropic_reduction()

R.isotropic_sublattice()          # I
R.orthogonal_complement()         # I^perp
R.quotient_lattice()              # I^perp/I
R.inclusion()                     # I -> I^perp
R.projection()                    # I^perp -> I^perp/I
```

`v.isotropic_reduction()` should require `v.is_primitive()`. For a nonprimitive isotropic vector, v^⊥/Zv has torsion and is not an integral lattice in the same category; the general operation should instead be named `v.isotropic_quotient()` and return a formed module.

## 1.3 Isometries and groups

```python
G = L.O()
Gplus = L.O_plus()                # kernel on A_L
Gomega = L.O_component()          # O^Omega(L), when defined

g = G.element(matrix)
g(v)
g.restrict(I)
g.induced_on(R.quotient_lattice())

G.gens()
G.stabilizer(v)
G.stabilizer(I, action="setwise")
G.stabilizer(I, action="pointwise")
G.centralizer(f)
G.intersection(H)

rho_A = G.discriminant_representation()
G.kernel(rho_A)
G.preimage(rho_A, H_finite)
```

A subgroup object must retain how it was constructed:

```text
GeneratedSubgroup
KernelSubgroup
FinitePreimageSubgroup
StabilizerSubgroup
CentralizerSubgroup
IntersectionSubgroup
```

This is required because generator filtering does not compute a predicate-defined subgroup.

## 1.4 Orbit domains

Infinite domains should be symbolic loci:

```python
X = L.vector_locus(norm=m, primitive=True)
Y = L.isotropic_sublattice_locus(rank=k)
F = L.isotropic_flag_locus(ranks=(1, 2))

D = G.orbit_decomposition(X)
```

An orbit decomposition should expose:

```python
D.representatives()
D.stabilizer(x)
D.transporter(x, representative)
D.acting_group()
D.domain()
```

The representatives should be lattice elements or subobjects, not raw row matrices.

---

# 2. Proposed implementation modules

```text
quadratic_lattices/
    core.py
    subobjects.py
    morphisms.py
    finite_forms.py

    groups/
        arithmetic.py
        finite_representations.py
        integral_structures.py
        groupoid_cache.py

    indefinite/
        invariants.py
        vector_sections.py
        isotropic_reductions.py
        isotropic_lifts.py
        eichler.py
        recursive.py
        lorentzian_cells.py
        isotropic_flags.py
        arithmetic_subgroups.py
        equivariant.py

    backends/
        definite.py
        polyhedral.py
        lorentzian_native.py
```

The dependency direction should be:

```text
core/subobjects/morphisms
        ↓
finite_forms + finite_representations
        ↓
integral_structures + isotropic_lifts
        ↓
definite backend + Lorentzian cell backend
        ↓
Eichler orbit covers
        ↓
recursive O(L), isometry, vector orbits
        ↓
isotropic planes/flags
        ↓
finite-index Γ and equivariant/centralizer variants
```

---

# 3. `indefinite/invariants.py`

## 3.1 Port `AttackScheme` as a signed view

Reference:

```text
src_indefinite/IndefiniteFormFundamental.h
    AttackScheme
    INDEF_FORM_GetAttackScheme
```

The source computes h = min(n_+, n_-) and replaces Q by -Q when necessary so that h becomes the positive index.

Translate this as:

```python
@dataclass(frozen=True)
class AttackProfile:
    lattice: IntegralQuadraticLattice
    signed_view: SignedFormView
    sign: int                    # +1 or -1
    positive_index: int          # min(n_+, n_-)
    negative_index: int
```

```python
profile = L.attack_profile()
```

`SignedFormView` should share the underlying module with L and replace b by -b. Do not construct a supposedly different mathematical lattice and silently transfer objects between them.

## 3.2 Replace hash invariants by structured prefilters

Reference:

```text
INDEF_FORM_Invariant
INDEF_FORM_InvariantVector
INDEF_FORM_Invariant_IsotropicKplane_Raw
```

The source combines rank, signature, determinant, parity, vector data, and reduced-form data into `size_t` hashes. These are valid candidate-bucketing devices but not mathematical invariants whose hash equality may be treated as equivalence.

Use:

```python
@dataclass(frozen=True)
class LatticePrefilter:
    rank: int
    signature: tuple[int, int, int]
    parity: str
    determinant: Integer
    discriminant_elementary_divisors: tuple[Integer, ...]


@dataclass(frozen=True)
class VectorPrefilter:
    norm: Integer
    content: Integer
    divisor: Integer
    discriminant_class_orbit_key: object
    orthogonal_reduction_prefilter: LatticePrefilter


@dataclass(frozen=True)
class IsotropicSubspacePrefilter:
    rank: int
    flag_type: FlagType
    quotient_prefilter: LatticePrefilter
    embedding_elementary_divisors: tuple[Integer, ...]
```

There is an important source-name correction. `INDEF_FORM_InvariantVector` calls one scalar extracted from v `divisor` and one extracted from Qv `index`. Semantically these should be separated as:

content(v) = max{ d : v ∈ dL },

and

div(v) = gcd b(v, L) = content(v^♭).

The source code should not determine the public terminology.

## 3.3 Do not port `ExpandMatrix` as a public primitive

`ExpandMatrix` merely inserts a reduced matrix into a chosen block basis. Its replacement is a direct-sum morphism:

```python
f = phi.direct_sum(identity_map(rank_one_factor))
```

The basis matrix used to express this morphism belongs in the coordinate adapter.

---

# 4. `indefinite/vector_sections.py`

## 4.1 Translate `INDEF_FORM_GetVectorStructure`

Reference:

```text
CombinedAlgorithms.h
    INDEF_FORM_GetVectorStructure
    MapOrthogonalSublatticeEndomorphism
    MapOrthogonalSublatticeGroup
```

The source computes the integral kernel of b(v,-), calls it `NSP`, restricts the Gram form, and then uses different lift constructions according to whether q(v) is zero.

The semantic replacement is:

```python
@dataclass(frozen=True)
class VectorOrthogonalSection:
    vector: LatticeElement
    perpendicular: LatticeSubobject
    inclusion: LatticeMorphism
    reduction: IntegralQuadraticLattice | IsotropicReduction

    def rational_lift(
        self,
        reduced_isometry,
        *,
        target: "VectorOrthogonalSection | None" = None,
    ) -> "IsometryExtensionTorsor":
        ...

    def integral_stabilizer(
        self,
        reduced_group,
    ) -> ArithmeticSubgroup:
        ...
```

Construction:

```python
def orthogonal_section(v):
    L = v.parent()
    P = kernel(L.metric_map().then(evaluation_at(v)))
    if not v.is_isotropic():
        return NonIsotropicVectorSection(v, P)
    return IsotropicVectorSection(v, P, v.isotropic_reduction())
```

More directly, `P` is the kernel of x ↦ b(x,v).

## 4.2 Nonisotropic branch

For q(v) ≠ 0, L_Q = Qv ⊥ v_Q^⊥.

Given an isometry φ: v_1^⊥ → v_2^⊥ with q(v_1) = q(v_2), the rational extension is uniquely specified by v_1 ↦ v_2, x ↦ φ(x) for x ∈ v_1^⊥.

This is what the source implements by constructing `Pmat`, inserting the reduced isometry with `ExpandMatrix`, and conjugating back.

The new method should be:

```python
extension = section1.rational_lift(phi, target=section2)
g = extension.integral_point(source=L1, target=L2)
```

The integrality correction must be a separate operation. A rational block extension is not automatically an integral lattice isometry.

## 4.3 Isotropic branch

The source does not quotient by Zv. Instead, it recursively handles the degenerate formed module v^⊥, then calls

```text
LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1
```

and changes sign when the resulting lift sends v to -v.

The semantic implementation should instead expose K_v = v^⊥/Zv as the primary reduced lattice. The source's degenerate recursion may remain as an interim backend, but it should not define the public object.

Use:

```python
R = v.isotropic_reduction()
P = R.parabolic_data()

Gv = P.point_stabilizer_group()
rho = P.action_on_quotient()      # Gv -> O(K_v)
U = P.unipotent_kernel()
```

This makes explicit that the stabilizer is not simply O(K_v).

---

# 5. `indefinite/isotropic_reductions.py`

## 5.1 Translate `INDEF_FORM_Rec_IsotropicKplane`

This C++ structure is the main source object to reinterpret.

| C++ field | Mathematical meaning |
|---|---|
| `Qmat` | ambient lattice form L |
| `Plane` | basis matrix for I |
| `NSP_T` | a basis for I^⊥ |
| `GramMatRed` | degenerate form on I^⊥ |
| `PlaneExpr` | coordinates of I ↪ I^⊥ |
| `TheCompl` | a noncanonical complement representing I^⊥/I |
| `QmatRed` | Gram matrix of K_I = I^⊥/I in that complement |
| `FullBasis` | the chosen splitting basis of I^⊥ |
| `LiftToFullAutomorphism` | one rational lift of a quotient isometry |
| `ComputeRelevantKernel` | kernel acting trivially on I^⊥ |
| `ComputeInvariantSublattice` | attempted common integral structure for lifts |

The constructor computes I^⊥, expresses I inside it, chooses an integral complement, and obtains a representative Gram matrix for I^⊥/I.

The replacement is:

```python
@dataclass(frozen=True)
class IsotropicReduction:
    ambient: IntegralQuadraticLattice
    isotropic: PrimitiveLatticeSubobject
    perpendicular: LatticeSubobject
    quotient: IntegralQuadraticLattice
    inclusion: LatticeMorphism
    projection: LatticeMorphism
    coordinate_frame: "IsotropicReductionFrame"

    def extension_problem(self, levi_element) -> IsometryExtensionTorsor:
        ...

    def parabolic_data(
        self,
        *,
        flag_type: FlagType | None = None,
    ) -> IntegralParabolicDatum:
        ...
```

`coordinate_frame` stores the equivalent of `PlaneExpr`, `TheCompl`, and `FullBasis`; it is explicitly a choice and is excluded from equality of `IsotropicReduction` objects.

## 5.2 Replace `SeqDims` by `FlagType`

Reference:

```text
SeqDims
seq_dims_plane
seq_dims_flag
seq_dims_reduced
seq_dims_append_one
f_get_list_spaces
GetAutomorphismOfFlag
```

The source encodes a plane as `dims={k}` and a complete flag as `dims={1,...,1}`; it then constructs the successive prefix spaces and a block-triangular group preserving them.

Use:

```python
@dataclass(frozen=True)
class FlagType:
    dimensions: tuple[int, ...]

    @classmethod
    def plane(cls, k):
        return cls((k,))

    @classmethod
    def complete(cls, k):
        return cls(tuple(range(1, k + 1)))

    @property
    def block_sizes(self):
        ...
```

An actual flag is:

```python
@dataclass(frozen=True)
class IsotropicFlag:
    sublattices: tuple[PrimitiveLatticeSubobject, ...]

    def flag_type(self) -> FlagType:
        ...
```

`GetAutomorphismOfFlag` and `ExtendIsometryGroup_Triangular` should not survive as raw block-matrix generators. Replace them with:

```python
P = GL(I).stabilizer(flag_on_I)
```

represented as the integral parabolic subgroup determined by `FlagType`.

## 5.3 Port the radical extension as a semidirect product

Reference:

```text
ExtendIsometryGroup
ExtendIsometryGroup_Triangular
ExtendIsometryGroup_IsotropicOrth
INDEF_FORM_AutomorphismGroup_Reduced
INDEF_FORM_TestEquivalence_Reduced
```

For a degenerate free formed module M with radical R, choose the quotient M̄ = M/R. After a module splitting,

Aut(M, b) ≅ Hom_Z(M̄, R) ⋊ ( O(M̄) × GL(R) ).

For a flag in R, replace GL(R) by its flag stabilizer.

This should become:

```python
D = L.radical_reduction()
G = D.automorphism_group(flag=None)
```

The source's elementary shear generators are then generated from the ordinary module Hom_Z(M̄, R), rather than inserted by coordinate loops.

---

# 6. `indefinite/isotropic_lifts.py`

This is one of the best pieces to reconstruct directly in Sage.

## 6.1 Port the codimension-one extension

Reference:

```text
Lorentzian.g
    LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1

lorentzian_linalg.h
    LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1_Basis
    LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1_Kernel
```

The algorithm chooses a complementary vector, solves the prescribed pairings against the target subspace, and adjusts it along the one-dimensional isotropic radical to obtain the required norm. The GAP implementation states this calculation directly.

Translate it as:

```python
class CodimensionOneIsotropicExtension:
    def __init__(self, source_subspace, target_subspace, partial_isometry):
        ...

    def rational_extension(self) -> LatticeIsometryOverQQ:
        ...
```

The return value should include the source and target subspace morphisms, not only a matrix.

## 6.2 Rewrite `SpecialEquationSolving` as a linear-map problem

Reference:

```text
SpecialEquationSolving
LORENTZ_ExtendOrthogonalIsotropicIsomorphism
```

The source flattens r×r matrices and solves B = XA + A^T X^T.

Define the linear map

Φ_A : Mat_r(Q) → Sym_r(Q), X ↦ XA + A^T X^T.

Then use Sage kernels and affine solution spaces:

```python
@dataclass(frozen=True)
class MatrixEquationSolution:
    particular: Matrix
    homogeneous_space: VectorSpace
    homogeneous_lattice: FreeModule | None


def solve_isotropic_extension_equation(A, B):
    phi = matrix_linear_map(
        domain=MatrixSpace(QQ, r, r),
        codomain=SymmetricMatrices(QQ, r),
        function=lambda X: X*A + A.T*X.T,
    )
    return phi.solve_affine(B)
```

Do not port the manual index map `f(i,j)` or the flattened dense matrix construction unless profiling shows that Sage's generic linear-map implementation is too slow.

## 6.3 Replace denominator elimination by affine-lattice intersection

The source computes a rational affine family T_0 + Σ_i λ_i T_i and calls

```text
EliminateSuperfluousPrimeDenominators_Matrix
```

to find an integral member.

The correct semantic operation is:

```python
E = IsometryExtensionTorsor(T0, [T1, ..., Tr])
E.integral_parameters(source_lattice=L1, target_lattice=L2)
E.one_integral_extension()
```

Integrality is a system of affine congruences in the λ_i. Clear denominators once and solve the resulting affine lattice problem by Smith normal form and Chinese remaindering. This gives either:

- an empty integral locus;
- or a coset λ_0 + Λ ⊂ Q^r.

No heuristic prime-denominator cancellation is required.

## 6.4 Port the pointwise kernel, not the entire unipotent radical

Reference:

```text
IntegralKernelSpecialEquation
GetOrthogonalTotallyIsotropicKernelSubspace
```

The latter computes isometries acting identically on I^⊥. It reduces the condition to HU^T + UH^T = 0 and constructs generators from the integral kernel. The source explicitly notes that this particular kernel is commutative.

Translate this as:

```python
K = R.pointwise_perpendicular_kernel()

K.parameter_lattice()
K.embedding_into(R.ambient.O())
K.gens()
```

Do not call this object the full unipotent radical. The full parabolic unipotent group may include additional shear directions and need not be commutative.

---

# 7. Replace `ComputeInvariantSublattice` by intrinsic parabolic data

This is the most important deliberate deviation.

The source methods are:

```text
INDEF_FORM_Rec_IsotropicKplane.LiftToFullAutomorphism
ComputeInvariantSublattice_method1
ComputeInvariantSublattice_method2
ComputeRelevantKernel
MapOrthogonalSublatticeGroupUsingSublattice
```

The source itself states that its attempted invariant-sublattice construction has no theoretical guarantee, although it had worked on tested inputs. Method 2 lifts all quotient generators and then feeds those rational matrices to `MatrixIntegral_GetInvariantSpace`.

Do not port either method.

## 7.1 Construct an `IntegralParabolicDatum`

For a primitive totally isotropic I ⊂ L, form P = I^⊥, K = P/I.

The pairing induces an injection of free rank-k modules

L/P ↪ I^∨, [x] ↦ b(x,-)|_I,

with finite cokernel. Record:

1. I;
2. P = I^⊥;
3. K = P/I;
4. the lattice J_I = im(L/P → I^∨);
5. the extension 0 → I → P → K → 0;
6. the finite gluing data recovering L inside a rational Witt decomposition;
7. any requested flag on I.

```python
@dataclass(frozen=True)
class IntegralParabolicDatum:
    ambient: IntegralQuadraticLattice
    isotropic: PrimitiveLatticeSubobject
    perpendicular: LatticeSubobject
    quotient: IntegralQuadraticLattice
    dual_pairing_lattice: LatticeSubobject
    extension_class: ModuleExtension
    gluing: FiniteGluingDatum
    flag_type: FlagType | None
```

## 7.2 Compute the Levi image intrinsically

There is a natural homomorphism

Stab_{O(L)}(I) → GL(I) × O(K).

Its image consists of pairs preserving the finite extension and gluing data.

```python
P = R.parabolic_data(flag_type=flag_type)

M = P.levi_candidate_group()      # GL(I) × O(K), or flag parabolic × O(K)
F = P.finite_compatibility_action()
M_integral = M.stabilizer(F.gluing_object())
```

This is a finite-module stabilizer problem and belongs to the integral-structure/libGAP layer.

## 7.3 Lift generators by exact affine solving

For each generator (a, g) of `M_integral`:

```python
E = P.extension_problem(a, g)
lift = E.one_integral_extension()
```

Then adjoin the relevant unipotent generators.

```python
P_group = generated_group(
    lifted_levi_generators + P.unipotent_generators()
)
```

This produces the actual parabolic stabilizer without guessing a common "helping lattice."

---

# 8. `groups/integral_structures.py`

This module should reconstruct the mathematics of `GroupAction.g` and `MatrixGroup.h`, but not their custom permutation-group implementation.

## 8.1 Target API

```python
QG = RationalMatrixGroup(V_Q, generators)

A = IntegralStructureAction(
    rational_group=QG,
    lattice=L,
)

M = A.invariant_overlattice()
rho = A.finite_representation()

H = A.lattice_stabilizer()
t = A.transporter(L1, L2)
C = A.right_cosets()

D = A.double_cosets(left_subgroup=V)
```

The double-coset result should state its sides:

```python
@dataclass(frozen=True)
class DoubleCosetDecomposition:
    left_subgroup: RationalMatrixGroup
    ambient_group: RationalMatrixGroup
    right_subgroup: ArithmeticSubgroup
    representatives: tuple
    intersections: tuple
```

For the source call used in isotropic-plane induction, this is V \ G / G_L, with G_L = G ∩ GL(L).

## 8.2 Port `MatrixIntegral_GetInvariantSpace` as a lattice closure

The source starts from the standard lattice and repeatedly replaces it by the Z-span of its images under every generator and inverse until the determinant stabilizes.

Translate this as:

```python
def smallest_invariant_overlattice(G, L):
    M = L
    while True:
        M_next = saturation(sum_of_sublattices(
            [M] + [g(M) for g in G.gens_and_inverses()]
        ))
        if M_next == M:
            return M
        M = M_next
```

The method should explicitly have the precondition that G stabilizes some lattice commensurable with L. It should not be advertised as a terminating procedure for every finitely generated subgroup of GL_n(Q).

## 8.3 Replace `LinearSpace_GetDivisor`

The GAP implementation increments d = 1, 2, … until it finds dZ^n ⊆ L.

Replace this by the exponent of the finite quotient. If L ⊆ M and the Smith invariants of M/L are d_1 | ⋯ | d_r, then the least d satisfying dM ⊆ L is d_r.

```python
d = (M / L).exponent()
```

## 8.4 Replace residue-vector enumeration by a finite-module action

Given the invariant over-lattice M and d with dM ⊆ L, define F = M/dM, S = L/dM ⊆ F.

Then G_L = ρ^{-1}( Stab_{ρ(G)}(S) ).

The source implements this through:

```text
LinearSpace_ModStabilizer
MapToPermutationOrderedOrbit
LinearSpace_Stabilizer_Kernel
LinearSpace_Stabilizer_RightCoset
LinearSpace_ModEquivalence
```

It repeatedly discovers a violating residue vector, enumerates its finite orbit, constructs a permutation action, and stabilizes the corresponding subset.

The Sage implementation should instead construct an explicit finite-module object:

```python
F = M.quotient(d * M)
S = F.submodule(L)

rho = QG.action_on(F)
P = rho.image_as_permutation_group(orbit_of=S)
stab = libgap.Stabilizer(P, S)
H = rho.preimage(stab)
```

For large F, refine prime by prime using F = ⊕_{p | d} F_(p) and the elementary-divisor filtration. This recovers the source's incremental modular refinement without treating residue vectors as the public abstraction.

## 8.5 Map the `MatrixIntegral_*` family directly

| Source symbol | Target method |
|---|---|
| `MatrixIntegral_Stabilizer_General` | `QG.stabilizer_of_lattice(L)` |
| `MatrixIntegral_Equivalence_General` | `QG.transporter_of_lattices(L1,L2)` |
| `MatrixIntegral_Equivalence_Bis_General` | `QG.integral_element_in_coset(t*QG)` |
| `MatrixIntegral_RightCosets_General` | `QG.right_cosets(QG.stabilizer_of_lattice(L))` |
| `MatrixIntegral_DoubleCosets_General` | `QG.double_cosets(V, QG.stabilizer_of_lattice(L))` |
| `LinearSpace_Stabilizer` | `G.stabilizer(LatticeSubobject)` |
| `LinearSpace_Equivalence` | `G.transporter(L1,L2)` |

The current C++ routines first find an invariant lattice, conjugate the rational group into integral coordinates, perform a finite modular action, and conjugate the results back.

## 8.6 Delegate all finite group work to libGAP

Do not port:

- `PersoGroup`;
- custom Schreier routines;
- custom permutation types;
- `SmallGeneratingSet`;
- direct double-coset enumeration.

Maintain a homomorphism ρ: G → P to the finite permutation group together with the correspondence between matrix generators and permutation generators. Lift GAP words back by evaluating them in the original matrices.

The parity corpus for this module is `CI_tests/01_RatIntAutomorphy/ProcessExamples.g`, which explicitly verifies that the returned double cosets are disjoint and exhaust the finite ambient group.

---

# 9. `indefinite/eichler.py`

## 9.1 Port `INDEF_FORM_Eichler_Transvection`

Reference:

```text
IndefiniteForms.g
ApproximateModels.h
    INDEF_FORM_Eichler_Transvection
```

For an even lattice, isotropic f, and x ∈ f^⊥,

E_{f,x}(y) = y + b(y,x) f − (q(x)/2) b(y,f) f − b(y,f) x.

Use:

```python
g = L.eichler_transvection(f, x)
```

The constructor checks:

```python
L.is_even()
f.is_isotropic()
L.b(f, x) == 0
```

and returns a `LatticeIsometry`, not a matrix.

## 9.2 Replace `ApproximateModel` callbacks by an immutable object

The source structure contains four callbacks:

```text
GetApproximateGroup
SetListClassesOrbitwise
GetCoveringOrbitRepresentatives
GetOneOrbitRepresentative
```

Use:

```python
class OrbitCoverModel:
    def lattice(self): ...
    def subgroup(self) -> ArithmeticSubgroup: ...
    def covering_representatives(
        self,
        norm,
        *,
        primitive,
    ) -> OrbitCover: ...

    def one_representative(self, norm, *, primitive): ...
    def refined_by(self, subgroup) -> "OrbitCoverModel": ...
```

`OrbitCover` must not be named `OrbitDecomposition`: its representatives may contain several members of the same full O(L)-orbit.

## 9.3 Translate the two-U model

Reference:

```text
INDEF_FORM_EichlerCriterion_TwoHyperplanesEven
GeneratorsSL2Z
InternalEichler
EnumerateVectorOverDiscriminant
SetListClassesOrbitwise
GetApproximateGroup
```

The source:

1. identifies L = 2U ⊕ K;
2. constructs left and right SL_2(Z)-actions on the 2U-block;
3. adds Eichler transvections;
4. enumerates discriminant classes of K;
5. refines those classes under known easy isometries;
6. constructs vector representatives satisfying the norm congruences.

The target object is:

```python
@dataclass(frozen=True)
class TwoHyperbolicPlaneDecomposition:
    lattice: IntegralQuadraticLattice
    U1: LatticeSubobject
    U2: LatticeSubobject
    complement: LatticeSubobject
    sum_isometry: LatticeIsometry
```

```python
class EichlerOrbitCover(OrbitCoverModel):
    decomposition: TwoHyperbolicPlaneDecomposition
```

`SetListClassesOrbitwise` should become:

```python
A = L.discriminant_module()
image = extra_group.image_on(A)
class_orbits = image.orbits(A)
model2 = model.refined_by(extra_group)
```

Do not mutate an internal list of coordinate representatives.

## 9.4 Preserve the primitive/nonprimitive decomposition

For m ≠ 0, the source runs through positive integers c such that c^2 | m and multiplies primitive representatives of norm m/c^2 by c. For m = 0, it returns only primitive isotropic vectors.

Expose this distinction:

```python
model.covering_representatives(norm=m, primitive=True)
model.covering_representatives(norm=m, primitive=False)
```

For `norm=0`, `primitive=False` is not a finite orbit problem and should not silently call the primitive routine.

## 9.5 Generalize `EichlerReduction` to an envelope

Reference:

```text
ResultHyperbolicPlane
GetHyperbolicPlane
EichlerReduction
GetEichlerHyperplaneBasis
INDEF_FORM_GetApproximateModel
```

`GetEichlerHyperplaneBasis` finds two scaled hyperbolic pairs and constructs an embedding of the input lattice into an over-lattice having a literal 2U-summand. `INDEF_FORM_GetApproximateModel` then computes the subgroup of the envelope model preserving the original embedded lattice by `LinearSpace_Stabilizer_RightCoset`.

Use:

```python
@dataclass(frozen=True)
class EichlerEnvelope:
    lattice: IntegralQuadraticLattice
    envelope: IntegralQuadraticLattice
    inclusion: LatticeMorphism
    two_U_decomposition: TwoHyperbolicPlaneDecomposition

    def stabilizer_of_original_lattice(self, group):
        ...
```

Then:

```python
envelope_model = EichlerOrbitCover(E.envelope)
A = envelope_model.subgroup()
A_L = A.stabilizer(E.inclusion.image())
cosets = A.right_cosets(A_L)
```

This is the semantic content of the source's `Embed`, `EmbedInv_T`, `QmatEichler`, and right-coset manipulations.

## 9.6 Rewrite the randomized hyperbolic-plane search

`GetHyperbolicPlane` repeatedly applies random integral perturbations, asks for isotropic vectors, and uses heuristic stopping counters to find a pair with small pairing.

Replace it with:

```python
L.find_hyperbolic_pair(
    *,
    primitive=True,
    method="auto",
)
```

Backend order:

1. use an already known orthogonal decomposition;
2. use an exact rational isotropy solver and integral saturation;
3. solve for w with q(w) = 0 and prescribed b(v,w);
4. use randomized reduction only as a candidate accelerator.

Every result returns actual elements v, w ∈ L and verifies q(v) = q(w) = 0, b(v,w) > 0.

## 9.7 Replace `GetFirstNorm`

The source seeks the first positive integer represented by the model and uses a vector of that norm as the splitting vector.

The recursion only requires a vector v with positive norm in the sign-normalized form, because then the positive index of v^⊥ drops by one. The smallest represented positive norm is a complexity heuristic, not a mathematical requirement.

Use:

```python
v = model.choose_splitting_vector(
    objective="minimize_recursive_complexity",
)
```

Possible scores include:

- |q(v)|;
- determinant of v^⊥;
- discriminant-module size of v^⊥;
- estimated quotient-cover size.

---

# 10. `indefinite/recursive.py`

## 10.1 Translate the core recursion closely

Reference:

```text
IndefiniteForms.g
    INDEF_FORM_Machinery_AllFct

CombinedAlgorithms.h
    IndefiniteCombinedAlgo
    INDEF_FORM_AutomorphismGroup_Kernel
    INDEF_FORM_TestEquivalence_Kernel
    INDEF_FORM_GetOrbitRepresentative_Reduced
    INDEF_FORM_StabilizerVector_Reduced
    INDEF_FORM_EquivalenceVector_Reduced
```

The current automorphism routine has exactly the following structure:

1. dispatch definite/Lorentzian/higher-Witt-index by `AttackScheme`;
2. construct an approximate model;
3. choose a represented positive norm and vector v;
4. insert generators of the approximate subgroup;
5. recursively compute Stab_{O(L)}(v);
6. obtain a finite covering list on the norm shell;
7. add every exact transporter from v to a candidate in its full orbit.

Translate it as:

```python
class IndefiniteOrthogonalAlgorithm:
    def orthogonal_group(self, L):
        profile = L.attack_profile()

        match profile.positive_index:
            case 0:
                return self.definite_backend.orthogonal_group(L)

            case 1:
                return self.lorentzian_backend.full_orthogonal_group(L)

            case _:
                model = self.orbit_cover_model(L)
                v = model.choose_splitting_vector()

                A = model.subgroup()
                Pv = self.vector_stabilizer(v)

                transporters = []
                for w in model.covering_representatives(
                    norm=L.q(v),
                    primitive=v.is_primitive(),
                ):
                    t = self.vector_transporter(v, w)
                    if t is not None:
                        transporters.append(t)

                return L.O_from_generators(
                    A.gens() + Pv.gens() + transporters
                )
```

This is one of the places where the source control flow should be preserved almost verbatim.

## 10.2 Translate lattice equivalence

The source:

1. compares attack indices;
2. chooses a splitting vector v_1 in L_1;
3. enumerates a covering list of vectors of the same norm in L_2;
4. calls the recursive vector-transporter routine;
5. optionally simplifies the resulting matrix by multiplying on both sides by approximate groups.

Use:

```python
def isometry(self, L1, L2):
    if not basic_invariants_match(L1, L2):
        return None

    model1 = self.orbit_cover_model(L1)
    v1 = model1.choose_splitting_vector()

    model2 = self.orbit_cover_model(L2)
    for v2 in model2.covering_representatives(
        norm=L1.q(v1),
        primitive=v1.is_primitive(),
    ):
        f = self.vector_transporter(v1, v2, target_lattice=L2)
        if f is not None:
            return f
    return None
```

Do not initially port `ExhaustiveMatrixDoubleCosetSimplifications`. That routine only shortens the matrix witness after equivalence has already been proved.

## 10.3 Translate vector stabilizers

The current source does:

```text
INDEF_FORM_GetVectorStructure
    ↓
O(restricted or degenerate Gram form)
    ↓
MapOrthogonalSublatticeGroup
    ↓
MatrixIntegral_Stabilizer_General
```

The target is:

```python
def vector_stabilizer(self, v):
    section = v.orthogonal_section()

    Gred = self.orthogonal_group(section.reduced_object())
    QG = section.rational_lift_group(Gred)

    return QG.stabilizer_of_lattice(v.parent())
```

For isotropic v, replace `reduced_object()` by the explicit parabolic datum as that implementation becomes available.

## 10.4 Translate vector transporters

The target sequence is:

```python
def vector_transporter(self, v1, v2, target_lattice=None):
    S1 = v1.orthogonal_section()
    S2 = v2.orthogonal_section()

    phi = self.isometry(S1.reduced_object(), S2.reduced_object())
    if phi is None:
        return None

    torsor = S1.rational_lift(phi, target=S2)
    return torsor.one_integral_extension()
```

This corresponds to the source's reduced-form equivalence, rational lift, and `MatrixIntegral_Equivalence_Bis_General` correction.

## 10.5 Translate vector orbit decomposition

The source obtains a finite approximate cover, simplifies candidates under the approximate subgroup, and then deduplicates by repeated exact `INDEF_FORM_EquivalenceVector` calls.

Use a structured bucket map:

```python
def vector_orbits(self, G, locus):
    model = self.orbit_cover_model(G.ambient_lattice())
    candidates = model.covering_representatives(
        norm=locus.norm,
        primitive=locus.primitive,
    )

    buckets = group_by(candidates, key=vector_prefilter)

    representatives = []
    for bucket in buckets.values():
        representatives.extend(
            exact_orbit_deduplication(
                bucket,
                transporter=G.transporter,
            )
        )

    return OrbitDecomposition(G, locus, representatives)
```

Store the transporters found during deduplication instead of discarding them.

---

# 11. Reduction and memoization

## 11.1 Keep `IndefiniteReduction` as an exact presentation change

The public C++ entry points first call `IndefiniteReduction`, work with a reduced Gram matrix, and conjugate the result back.

Use:

```python
@dataclass(frozen=True)
class ReducedLatticePresentation:
    source: IntegralQuadraticLattice
    reduced: IntegralQuadraticLattice
    isometry: LatticeIsometry
```

```python
R = L.reduced_presentation()
```

No reduced Gram matrix should be detached from the isometry relating it to L.

## 11.2 Do not expose `ApproxCanonicalIndefiniteForm`

`IndefApproxCanonical.h`:

- reduces connected blocks;
- canonically orders absolute-value patterns by graph canonicalization;
- changes coordinate signs;
- orders blocks by coefficient-size heuristics.

This is not a canonical form for integral-lattice isometry. Retain it only as a private cache-bucketing and matrix-size heuristic:

```python
L.presentation_bucket_key()
```

An exact equivalence test is still mandatory after a key collision.

## 11.3 Generalize `DatabaseResultEquiStab` to an isometry-groupoid cache

The source stores:

- known isometries;
- known nonisometries;
- known stabilizer generators;

and transports stabilizers by conjugation along known isometries.

This is mathematically useful. Translate it as:

```python
class IsometryGroupoidCache:
    def remember_isometry(self, f: LatticeIsometry): ...
    def remember_nonisometric(self, L1, L2): ...
    def lookup_isometry(self, L1, L2): ...

    def remember_orthogonal_group(self, L, G): ...
    def lookup_orthogonal_group(self, L): ...
```

If f: L → M and O(M) is known, return f^{-1} O(M) f as O(L), with actual conjugated `LatticeIsometry` elements.

---

# 12. `indefinite/lorentzian_cells.py`

The Lorentzian perfect-domain code should initially remain partly native, but the native boundary should be much smaller than `LORENTZ_GetGeneratorsAutom`.

## 12.1 Use `DataPerfectLorentzianFunc` as the extraction seam

The existing C++ object already has approximately the desired local protocol:

| Existing method | Target method |
|---|---|
| `f_init` | `initial_cell()` |
| `f_hash` | `cell_bucket_key(cell)` |
| `f_repr` | `cell_transporter(C1,C2)` |
| `f_adj` | `adjacent_cells(C)` |
| `f_adji_obj` | `adjacency_target(edge)` |
| `LORENTZ_ComputeStabilizer` | `cell_stabilizer(C)` |
| `LORENTZ_DoFlipping` | `flip_across(C,facet)` |

`f_adj` computes the cell stabilizer, enumerates facet orbits by dual description, and flips across each facet. `f_repr` performs an exact integral configuration isomorphism test.

Define:

```python
@dataclass(frozen=True)
class LorentzianPerfectCell:
    lattice: IntegralQuadraticLattice
    vector_configuration: tuple[LatticeElement, ...]
    mode: str                     # "total" or "isotropic"


@dataclass(frozen=True)
class LorentzianCellAdjacency:
    source: LorentzianPerfectCell
    facet: PolyhedralFace
    target: LorentzianPerfectCell
    transporter: LatticeIsometry
```

```python
class LorentzianPerfectLocalBackend:
    def initial_cell(self, L, mode): ...
    def cell_stabilizer(self, C): ...
    def facet_orbits(self, C): ...
    def flip_across(self, C, facet): ...
    def cell_transporter(self, C1, C2): ...
```

## 12.2 Keep these functions native initially

The minimal native residue is:

```text
LORENTZ_FindPositiveVectorsKernel
GetUpperBound
LORENTZ_Kernel_Flipping
LORENTZ_GetOnePerfect
LORENTZ_DoFlipping
LORENTZ_ComputeStabilizer
LORENTZ_TestEquivalence
```

These perform exact CVP-like enumeration, wall movement, and cell construction.

Do not keep native:

- `EnumerateAndStore_Serial`;
- final group extraction;
- global memoization/database management;
- command-line configuration and serialization.

Python should own the quotient-cell traversal.

## 12.3 Rewrite the quotient traversal as a graph-of-groups computation

```python
class LorentzianPerfectComplex:
    def quotient_cells(self):
        ...

    def component_preserving_group(self):
        cells = self.quotient_cells()
        return generated_group(
            all_cell_stabilizers(cells)
            + all_side_pairings(cells)
        )
```

This is the semantic version of `LORENTZ_ExtractGeneratorsFromObjList`, which collects adjacency transporters and cell stabilizers.

The result of the cell traversal is O^Ω(L). The full group method separately adjoins a verified component-flipping isometry. In Lorentzian signature, −1_L is canonical. The legacy GAP implementation explicitly begins with −I because the perfect-domain generators do not flip the cone.

```python
def full_orthogonal_group(self):
    Gomega = self.component_preserving_group()
    return generated_group(Gomega.gens() + [-self.L.identity_isometry()])
```

`L.O_plus()` is subsequently computed as the kernel of the discriminant action; it is not the output of the Lorentzian cell traversal.

## 12.4 Generalize the norm-zero orbit logic to marked cells

`LORENTZ_GetOrbitRepresentative_Kernel` currently:

1. computes local orbits of isotropic vertices under each cell stabilizer;
2. transports them across adjacency edges;
3. builds a graph on the local orbit labels;
4. takes connected components.

For nonzero norm it stops with an explicit "some code needs to be written" error.

Extract the implemented logic as:

```python
class MarkedCellOrbitAlgorithm:
    def local_marks(self, cell): ...
    def local_stabilizer_orbits(self, cell): ...
    def transport_marks(self, edge): ...
    def global_orbits(self): ...
```

Then isotropic-vector orbits are one instance:

```python
marks = lambda C: C.isotropic_vertices()
```

Nonzero fixed-norm vectors require a new finite local-mark theorem or enumeration, not removal of the source's exception.

The Lorentzian parity corpus is `CI_tests/28B_LorentzianPerfStabEqui`, whose test driver compares exact perfect-domain counts in ranks 3, 4, 5, and > 5 with a frozen result file.

---

# 13. `indefinite/isotropic_flags.py`

This is the most complex high-level translation.

## 13.1 Source symbols

```text
CombinedAlgorithms.h
    f_stab
    f_stab_plane_v
    f_double_cosets

    INDEF_FORM_Stabilizer_IsotropicKstuff_Reduced
    INDEF_FORM_Equivalence_IsotropicKstuff_Reduced
    INDEF_FORM_RightCosets_IsotropicKstuff_Reduced
    INDEF_FORM_GetOrbit_IsotropicKstuff_Method
    INDEF_FORM_GetOrbit_IsotropicKstuff_Kernel
```

The current C++ starts from primitive isotropic-vector orbits and inductively extends rank-(k−1) sublattices by isotropic vectors in the quotient. It uses the double-coset path by default.

## 13.2 Express the induction intrinsically

Let I represent a G-orbit of primitive isotropic rank-(k−1) sublattices. Put

K_I = I^⊥/I, P_I = Stab_G(I), H_I = im(P_I → O(K_I)).

Let u represent an O(K_I)-orbit of primitive isotropic vectors and set Q_u = Stab_{O(K_I)}(u).

Then the H_I-orbits inside the O(K_I)-orbit of u are represented by Q_u \ O(K_I) / H_I in the source's right-action convention.

Translate this directly:

```python
class IsotropicFlagOrbitAlgorithm:
    def extend_one_rank(self, G, previous_orbits, flag_type):
        result = OrbitAccumulator(G)

        for I in previous_orbits.representatives():
            R = I.isotropic_reduction()
            K = R.quotient_lattice()

            P = G.stabilizer(I, action="setwise")
            H = P.induced_image(K)

            OK = K.O()

            for u in OK.primitive_isotropic_vector_orbits().representatives():
                Qu = OK.stabilizer(u)

                for d in OK.double_cosets(left=Qu, right=H):
                    u_d = d.representative()(u)
                    J = R.lift_isotropic_line(u_d.sublattice()).saturation()
                    result.insert_with_exact_transporter(J)

        return result.finish()
```

The group-action convention should be normalized in this method; callers should never reverse a double coset by hand.

## 13.3 Interpret the source helper methods

`f_stab(eRec, sd)` computes a rationally extended group preserving the existing isotropic plane or flag.

Translate as:

```python
P = R.parabolic_data(flag_type=sd)
P.rational_candidate_group()
```

`f_stab_plane_v(eRec, v, sd)` constructs the subgroup additionally preserving the next quotient isotropic vector, including transformations adding elements of the old isotropic space.

Translate as:

```python
P.stabilizer_of_quotient_line(v.primitive_line())
```

`f_double_cosets` computes the double cosets between:

- the rational candidate group;
- the rational subgroup stabilizing the extension vector;
- the integral subgroup preserving the ambient lattice.

The source comments explicitly identify this as an orbit-splitting problem and invoke `MatrixIntegral_DoubleCosets_General`.

Translate it as:

```python
D = IntegralStructureAction(G_rational, L).double_cosets(
    left_subgroup=V_rational
)
```

Once `IntegralParabolicDatum` is available, this can be simplified further by computing H_I directly and performing the quotient-lattice double cosets there.

## 13.4 Preserve final exact deduplication

After lifting candidate quotient lines, the source uses:

```text
INDEF_FORM_Equivalence_IsotropicKstuff_Reduced
```

to merge candidates only after an actual ambient integral isometry has been found.

Retain that structure:

```python
accumulator.insert_with_exact_transporter(J)
```

The prefilter may use:

- rank;
- flag type;
- J^⊥/J;
- discriminant gluing data;

but it may not decide equivalence.

## 13.5 Do not expose arbitrary k until the rewrite is complete

The source says that the initial-set computation "works by kind of chance" because only k = 1, 2 had been used.

Therefore:

1. implement rank-two planes and flags first;
2. validate the exact parabolic/gluing construction;
3. only then enable arbitrary `rank=k`.

For signature (2,n), rank two is already the maximal isotropic rank and recovers the Baily–Borel/Tits-building applications.

---

# 14. `indefinite/arithmetic_subgroups.py`

This is a generalization above the upstream algorithm.

## 14.1 Finite-preimage groups

```python
rho = G.discriminant_representation()
Gamma = G.preimage(rho, H)
Oplus = G.kernel(rho)
```

For any finite representation ρ: G → F and H ≤ F, Γ = ρ^{-1}(H) is computed by a finite subgroup-preimage calculation in libGAP, with words lifted to isometries.

This covers:

- O^+(L);
- prescribed subgroups of O(A_L);
- component kernels;
- congruence groups;
- intersections of these conditions.

## 14.2 Orbit splitting

For a full G-orbit representative x, let P_x = G_x. Then

P_x \ G / Γ ≅ ρ(P_x) \ ρ(G) / H.

Use:

```python
D_gamma = G.split_orbit(
    full_orbit=D_G[x],
    subgroup=Gamma,
)
```

The method returns:

- the split representatives;
- lifts of the finite double-coset representatives;
- stabilizers in Γ.

This is a separate use of double cosets from the internal integralization double cosets in `f_double_cosets`; both can share the same finite-group abstraction.

---

# 15. `indefinite/equivariant.py`

Centralizers are new functionality, not a small modification of the source recursion.

## 15.1 Public API

```python
Lf = L.with_isometry(f)
C = L.O().centralizer(f)

C = Lf.centralizer_group()
C.stabilizer(v)
C.stabilizer(I)
C.intersection(L.O_plus())
```

## 15.2 Involution implementation

For f^2 = 1, construct L_± = L ∩ ker(f ∓ 1). Let M = L_+ ⊕ L_- ⊆ L.

The finite-index overlattice L/M is the gluing object. Then

O(L, f) ≅ Stab_{O(L_+) × O(L_-)}(L/M).

Implementation:

```python
class EquivariantLattice:
    def eigensublattices(self): ...
    def split_sublattice(self): ...
    def gluing_datum(self): ...

    def centralizer_group(self):
        Lp, Lm = self.eigensublattices()
        G0 = Lp.O().direct_product(Lm.O())

        H = self.gluing_datum()
        Hstab = G0.stabilizer(H, action="setwise")

        return Hstab.extend_to_overlattice(self.lattice)
```

This reuses:

- orthogonal-group computation on the two eigensublattices;
- finite discriminant/gluing actions;
- subgroup preimages;
- exact extension to an overlattice.

No generator filtering occurs.

## 15.3 Finite-order generalization

For finite-order f, replace the ±1-decomposition by the rational cyclotomic decomposition

L_Q = ⊕_{d | ord(f)} ker Φ_d(f).

Intersect each rational isotypic component with L, compute the product of componentwise centralizers, and impose the finite equivariant gluing condition.

A general infinite-order centralizer should not be claimed by this initial implementation.

---

# 16. Exact symbol-level disposition

## Port semantically

```text
IndefiniteFormFundamental.h
    AttackScheme
    INDEF_FORM_GetAttackScheme
    INDEF_FORM_IsEven

IndefiniteForms.g / ApproximateModels.h
    INDEF_FORM_Eichler_Transvection
    INDEF_FORM_EichlerCriterion_TwoHyperplanesEven
    GetSquareDivisors
    INDEF_FORM_GetApproximateModel

CombinedAlgorithms.h
    INDEF_FORM_AutomorphismGroup_Kernel
    INDEF_FORM_TestEquivalence_Kernel
    INDEF_FORM_GetOrbitRepresentative_Reduced
    INDEF_FORM_StabilizerVector_Reduced
    INDEF_FORM_EquivalenceVector_Reduced

Lorentzian.g / lorentzian_linalg.h
    LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1
    SpecialEquationSolving
    LORENTZ_ExtendOrthogonalIsotropicIsomorphism
    IntegralKernelSpecialEquation
    GetOrthogonalTotallyIsotropicKernelSubspace

GroupAction.g / MatrixGroupBasic.h
    MatrixIntegral_GetInvariantSpace
```

## Generalize into typed abstractions

```text
ApproximateModel
    -> OrbitCoverModel

INDEF_FORM_GetVectorStructure
    -> VectorOrthogonalSection

INDEF_FORM_Rec_IsotropicKplane
    -> IsotropicReduction + IntegralParabolicDatum

SeqDims
    -> FlagType

DataPerfectLorentzianFunc
    -> LorentzianPerfectLocalBackend

DatabaseResultEquiStab
    -> IsometryGroupoidCache

MatrixIntegral_* and LinearSpace_*
    -> IntegralStructureAction
```

## Rewrite rather than port

```text
ComputeInvariantSublattice_method1
ComputeInvariantSublattice_method2
    -> intrinsic extension/gluing compatibility

LinearSpace_GetDivisor
    -> exponent from Smith normal form

GetHyperbolicPlane
    -> exact hyperbolic-pair solver with optional randomized acceleration

GetFirstNorm
    -> splitting-vector selection strategy

INDEF_FORM_Invariant*
    -> structured prefilters, not size_t hashes

INDEF_FORM_GetOrbit_IsotropicKstuff_Method
    -> quotient/parabolic/double-coset induction using actual subobjects

LORENTZ_GetOrbitRepresentative_Kernel
    -> generic marked-cell orbit traversal
```

## Delegate

```text
PositiveNegative.h
LatticeStabEquiCan.h
LATT_Automorphism
LATT_Isomorphism
short-vector/CVP leaves
    -> Sage/Oscar/PARI/fplll-backed definite lattice backend

finite permutation groups, stabilizers, transporters, double cosets
    -> libGAP

generic facets, rays, redundancy, LP
    -> Normaliz/cddlib/PPL/Sage polyhedral backends

HNF, SNF, saturation, kernels, quotient modules
    -> Sage/FLINT
```

## Retain native temporarily

```text
LORENTZ_FindPositiveVectorsKernel
GetUpperBound
LORENTZ_Kernel_Flipping
LORENTZ_GetOnePerfect
LORENTZ_DoFlipping
LORENTZ_ComputeStabilizer
LORENTZ_TestEquivalence
```

## Drop

```text
Indefinite.jl/src/Functions.jl
Indefinite.jl/src/ExternalCalls.jl
gap_polyhedral/lib/stubs.g
INDEF_FORM_*.cpp command-line programs
SystemNamelist and output-format code
temporary-file serialization
custom global debug/timing macros
custom matrix/vector carrier classes
ExhaustiveMatrixDoubleCosetSimplifications
ExhaustiveReductionComplexityGroupMatrix
```

The first three are wrappers rather than mathematical implementations.

---

# 17. Implementation order and concrete recovery points

## Tranche 1: subobjects, reductions, and exact lifts

Implement:

```text
VectorOrthogonalSection
IsotropicReduction
FlagType
IsometryExtensionTorsor
solve_isotropic_extension_equation
pointwise_perpendicular_kernel
```

Recover v^⊥/Zv for primitive isotropic vectors in U ⊕ E_8(−1), together with exact quotient isometries and rational/integral lifts.

## Tranche 2: integral structures and double cosets

Implement:

```text
IntegralStructureAction
smallest_invariant_overlattice
finite module action M/dM
lattice stabilizer
lattice transporter
right cosets
double cosets
```

Run the complete `CI_tests/01_RatIntAutomorphy` and `CI_tests/DoubleCosets/DBL` corpora. The former checks rational versus integral automorphisms and exact disjoint exhaustion by double cosets.

## Tranche 3: Lorentzian cell backend

Expose only the local cell protocol from C++, with Python owning quotient traversal and group extraction.

Recover:

- O^Ω(U ⊕ E_8(−1));
- full O(U ⊕ E_8(−1));
- primitive isotropic-vector orbits;
- their stabilizers;
- the frozen Lorentzian perfect-domain counts.

## Tranche 4: Eichler model and recursive full group

Implement:

```text
EichlerEnvelope
EichlerOrbitCover
IndefiniteOrthogonalAlgorithm
vector stabilizer
vector transporter
lattice isometry
full O(L)
```

Recover, for N = U ⊕ U(2) ⊕ E_8(−2), full O(N), explicit isometries, primitive isotropic-vector orbits, and the action on A_N.

The upstream broad regression file exercises automorphisms, equivalence, norm-0 and norm-2 vector orbits, and isotropic planes/flags on six lattice families including the Enriques case.

## Tranche 5: exact parabolics and rank-two isotropic planes

Implement:

```text
IntegralParabolicDatum
parabolic Levi image
unipotent generators
isotropic rank extension
quotient double-coset splitting
```

Recover the published line, plane, and flag orbits for the Enriques and polarized K3 lattices from the preceding oracle addendum.

## Tranche 6: finite-index arithmetic groups

Implement:

```text
G.kernel(rho)
G.preimage(rho,H)
G.split_orbit(...)
```

Recover O^+(L), O^{+,Ω}(L), and the published arithmetic-subgroup cusp buildings.

## Tranche 7: centralizers and intersections

Implement:

```text
EquivariantLattice
G.centralizer(f)
construction-aware intersections
```

Recover the centralizer of an Enriques involution from eigensublattices and gluing, then intersect it with stable and polarization-stabilizer subgroups.

The resulting Sage code owns the mathematical recursion and object semantics. The only initially retained C++ is the local Lorentzian perfect-cell engine; even there, traversal, orbit assembly, subgroup semantics, and result construction move into the Sage layer.
