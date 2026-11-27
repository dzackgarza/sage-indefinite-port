# Allcock edgewalk source correspondence

These notes preserve the source correspondence and mathematical rationale from parent revision `1cba545b64637f873bf0b801bd6157ff1a042130`, with original source line numbers. Their relocation does not certify the whole algorithm.

## `src/sage_indefinite_port/indefinite/edgewalk.py`

### `edgewalk_fundamental_domain`

Lines 88–89:

A prescribed corner belongs to the original presentation and must
be used as the start, including in the special rank-two cycle.

### `edgewalk_isometry`

Lines 147–151:

The rank-two edgewalk has a binary cycle rather than the
higher-dimensional finite vertex-orbit graph. A chosen initial
wall need not match one from a different integral presentation.
The exact binary Lorentzian isometry owner handles both the
isotropic and anisotropic cases, including the cycle period.

Lines 158–160:

Absence of an edgewalk starting corner is not an isometry
obstruction. The complete Lorentzian perfect-cell equivalence
algorithm remains available for this non-vertex branch.

Lines 182–185:

A nonreflective edge has no next chamber vertex. The
vertex-orbit graph is then not a complete isometry
obstruction: use the complete perfect-cell equivalence
algorithm rather than returning a false negative.

### `_Edgewalk.__init__`

Lines 406–409:

edgewalk.h:1985-1999 raise NonReflectivityException from
get_simple_cone; StandardEdgewalkAnalysis (edgewalk.h:2175-2183)
reports the reason. A starting corner that certifies nonreflectivity
therefore has to carry that reason to the report.

Lines 412–414:

The congruence matrix and its homogeneous integral kernel in
Get_DataMapping depend on the root norm, not the affine coset.
They are reused for every height of one Vinberg search.

### `_Edgewalk.initial_isotropic_vertex`

Lines 502–509:

k^perp/k on engine data (edgewalk.h:1955-1968, LatticeProjectionFramework):
the rows K of an integral basis of k^perp (the integral left kernel
of G k), and the Smith transform U of the column c^T of the
coordinates of k in K, U c^T = (+-1, 0, ..., 0)^T, as Sage frames
the cokernel of k -> k^perp (the owned isotropic reduction's Smith
generators).  A vector x K of k^perp projects to the last
coordinates of U x^T; the columns of U^-1 after the first lift
the quotient generators.

Lines 523–524:

edgewalk.h:1951-1954: Latt_inter_NSP = IntersectionLattice(Latt, NSP),
the integral left kernel of B G k on the rows B of a basis of C_N.

Lines 531–532:

fr.BasisProj (edgewalk.h:1955-1960): the projected generators,
read through the perpendicular basis and the quotient framing.

Lines 540–543:

The norm shell of the image lattice (edgewalk.h:1969-2000) on its
Gram matrix in the basis image_basis: PARI qfminim, the vectors
of square exactly N, each followed by its negative, as the
definite shell enumeration orders them.

Lines 554–554:

fr.GetOnePreimage (edgewalk.h:2010): one preimage in C_N cap k^perp.

Lines 561–564:

get_simple_cone raises NonReflectivityException here
(edgewalk.h:1985-1988) and StandardEdgewalkAnalysis turns it
into the report's reason (edgewalk.h:2175-2183): "The list of
vectors is empty. Cannot be reflective."

Lines 569–570:

edgewalk.h:1990-1999: rnk < G.rows() - 2 means the roots at the
cusp cannot be a reflective cusp's simple system.

Lines 574–580:

GetFacetOneDomain_ListIdx (lorentzian_linalg.h:267–339):
regular splitting vector, positive half and its extreme rays.
A regular splitting is outside the finitely many hyperplanes
annihilating the nonzero projected roots. On the moment curve
(1,t,...,t^(rank-1)), each forbidden hyperplane has at most
rank-1 integer parameters t. Exhausting t=0,1,... therefore
constructs a regular splitting without probabilistic retries.

Lines 590–591:

The extreme rays of the cone of the positive projected roots in
k^perp/k (pplpy dual description, as the owned cone's engine).

Lines 597–598:

get_one_root (edgewalk.h:2001-2022): the admissible root on this
ray, i.e. a positive multiple of it (the quotient is definite).

Lines 606–607:

Source get_simple_cone accepts non-simplicial faces. Test rank in
the isotropic quotient, not in the ambient lattice of their lifts.

Lines 611–613:

edgewalk.h:2020–2049: choose a negative kP orthogonal to the
simple walls, orient by b(kP,k)<0 and complete the cusp by the
source's affine-line equation fibres.

### `_Edgewalk._vinberg_polytope_roots`

Lines 632–635:

LLL_gram returns U with U^t G U reduced. Source Pmat = U^t.
Both this change of coordinates and FullGramInfo's exact
Bareiss data depend on the root-norm kernel, not on the
affine shift, height, or current system of facet constraints.

### `_Edgewalk._vinberg_polytope_roots.last_true`

Lines 665–665:

Shvec_LastTrue: exact doubling and bisection.

### `_Edgewalk._vinberg_polytope_roots.spherical_bounds`

Lines 683–683:

Exact Bound_Ceil / Bound_Floor, with no inexact seed.

### `_Edgewalk._vinberg_polytope_roots`

Lines 695–700:

The section at level i substitutes precisely the coordinates
j > i already chosen by the spherical recursion.  Carry these
constant terms down the search tree, just as the fraction-free
Shvec recursion carries its affine offset.  Constructing them
afresh from every suffix for each LP node redundantly recomputes
all earlier coordinate assignments.

### `_Edgewalk._vinberg_polytope_roots.enumerate_level`

Lines 723–728:

A nonempty section need not be pointed.  For example, a
bounded strip times a line has no extreme points, although
its intersection with the positive-definite sphere is finite.
Such a section cannot supply a vertex-derived coordinate
bound.  Retain the exhaustive spherical bounds instead of
discarding its valid integral solutions.

### `_Edgewalk._vinberg_decomposed_roots`

Lines 795–798:

A = [2GM | -kI] has full row rank, and its Smith
transformations are independent of the affine offset a.
The last n-1 columns of the right transformation give an
integral kernel basis for every congruence fibre.

Lines 815–818:

Sage smith_form returns (D,P,Q) with P*A*Q=D.  Since
-k*I is a block of A, its row rank is n.  The final n-1
columns of Q project to the *integral* kernel in w, not
merely the rational nullspace of the congruence.

### `_Edgewalk.prescribed_vertex`

Lines 874–876:

Both FileVertexRoots (:2107) and the FileVertex/common return
(:2143) apply RemoveFractionVector. The roots of FileVertexRoots are
data and are not primitive-normalised with the generator.

### `_Edgewalk.initial_vertex`

Lines 913–916:

ordinary_vertex_roots has already enumerated the complete
perpendicular root family and selected its cone rays.
The require_full path would enumerate all fixed-norm vectors
again, although its sole remaining condition is this count.

Lines 920–931:

vinberg.h FindRoots_Kernel / FindOneInitialRay: examine the accepted
root system after its initial finite cone and after each new root,
stopping on the first feasible nonpositive ray.

GetVinbergAux (vinberg.h:213-307) takes the saturated integral
perpendicular and ComputeTranslationClasses for
L / (v0^perp + Z*v0).  Its pinned basic_common_cpp@995fd625
MAT_MatrixInt.h:2881-2958 traverses the classes breadth first,
starting at zero and adding the ambient basis vectors in order.
Preserve that order through the actual module cokernel: the root
iterator starts at index 1 and consumes W sequentially.  Neither
ComputeTranslationClasses nor GetVinbergAux sorts by height.

Lines 954–955:

GetVinbergAux uses integer QuoInt (vinberg.h:258-260).
Division in QQ computes the field quotient, including for //.

Lines 958–959:

GetVinbergAux uses QuoInt/ResInt to place each class in
the fundamental timelike strip, 0 <= -<w,v0> < -v0^2.

Lines 964–968:

vinberg.h FundCone -> GetFacetOneDomain: select one member of
each perpendicular root pair by a regular linear functional,
and retain the original normed roots on the extreme cone rays.
Lexicographic sign is a deterministic regular functional for this
finite integral root family.

Lines 970–972:

GetFacetOneDomain first reduces the root coordinate columns to
their actual span, and GetFacetOneDomain_ListIdx requires at least
twice that reduced dimension in the nonempty root family.

Lines 978–980:

GetFacetOneDomain_ListIdx allocates exactly half the roots on
the assumption that the family is centrally symmetric. A count
alone cannot establish the required opposite-root correspondence.

Lines 991–992:

GetFacetOneDomain_ListIdx selects nonredundant *indices*;
parallel root vectors cannot both survive its reduction.

### `_Edgewalk.initial_vertex.accepted_vertex`

Lines 1022–1026:

GetOneInteriorVertex: exact dual description of -R*G,
tested on each extreme ray in the selected time component.
GetOneInteriorVertex uses the covectors -G*r. The
lattice's algebraic correlation produces these directly,
without assuming that primal and dual generator labels agree.

Lines 1029–1033:

GetOneInteriorVertex (vinberg.h:833-909) accepts the first
nonspacelike extreme ray of the feasible cone, regardless
of its time orientation.  The perpendicular-root chamber
was selected by an independent regular functional, so it
need not lie in the component containing the initial v0.

### `_Edgewalk.initial_vertex`

Lines 1044–1049:

vinberg.h FindRoots_Kernel tests its exit only after every root
belonging to one decomposition request has been inserted, so an
intermediate root cannot end the search.
vinberg.h IterateRootDecompositions::get_k (pinned 1592b246)
prioritizes each length's next coset by val/sqrt(val), with
val = -<v0, a>.

Lines 1063–1064:

vinberg.h FindRoot_filter/fct_CVP retains every norm-k
candidate satisfying all accepted wall inequalities.

Lines 1073–1075:

vinberg.h ReduceListRoot removes generators redundant in the
accepted root cone. Keep the original normed root vectors,
using primitive rays only to identify their directions.

Lines 1078–1080:

ReduceListRoot/get_non_redundant_indices preserves the order
of the original family. Keep its first normed representative
of each extremal direction rather than reordering by rays.

Lines 1084–1085:

Compare primitive directions once per accepted root rather
than testing every root against every extremal ray.

### `_Edgewalk.ordinary_vertex_roots`

Lines 1106–1113:

edgewalk.h:1845-1866: one FindFixedDistVectors per norm on the
Gram matrix of the definite complement k^perp.  The complement's
basis rows and Gram matrix are its engine data; PARI qfminim (the
engine of the complement's vectors_of_square) returns one vector
of each sign pair, in the order that shell lists them, and only
the accepted roots are wrapped as lattice vectors.
The rows of k^perp: the echelon basis of the integral left kernel
of G k, the basis of the owned orthogonal complement.

Lines 1130–1131:

The upstream ordinary-cone constructor rejects a root
family whose facet count is not dimSpace (including zero).

Lines 1134–1137:

edgewalk.h get_simple_cone_from_lattice selects actual normed
roots from FindFixedDistVectors, not primitive multiples of them.
Use geometric rays only to identify which original roots span
the selected cone, retaining each root's prescribed norm.

Lines 1140–1141:

edgewalk.h get_simple_cone_from_lattice requires exactly
dimSpace selected walls for an ordinary starting vertex.

### `_Edgewalk.constraint_rows`

Lines 1167–1170:

C = {v in Z^n : G v in d Z^n}.  The integral kernel of [G | -d I]
is {(v, w) : G v = d w}; G is nondegenerate, so w is determined
by v and the projection to v is an isomorphism onto C.  Its
Hermite form is the echelon basis of the ZZ-module intersection.

### `_Edgewalk.edge_orientation`

Lines 1215–1221:

A vector of positive square in the Lorentzian plane, on its
basis rows b_0, b_1 with Gram g: b_0 if g_00 > 0; else, for
g_00 < 0, b_1 - (g_01/g_00) b_0 of square det(g)/g_00 > 0; and for
g_00 = 0, b_1 + t b_0 of square g_11 + 2 t g_01 = 1.  Every such
vector with b(., k) < 0 lies in one component of the positive
cone, which k bounds, so the orientation of (positive, k) and
the isotropy of the plane do not depend on the choice.

### `_Edgewalk.edge_plane_for_divisor`

Lines 1287–1289:

The ZZ span of the projected rows C K / e is (1/e) times the span of
the integer rows C K; the Hermite normal form is unique and scales
with a positive factor, so this is the echelon basis of Sage's ZZ span.

Lines 1294–1296:

The basis (h_0, h_1) / e is positively oriented against (o_0, o_1)
exactly when the 2 x 2 minors of both on one pair of columns i, j
(with nonzero minor of o) have one sign (Cramer's rule).

### `_Edgewalk.edge_plane_for_divisor.coordinates`

Lines 1315–1315:

e * row = x h_0 + y h_1 on columns i, j.

### `_Edgewalk.edge_plane_for_divisor`

Lines 1321–1325:

edgewalk.h:763-764, Basis_P_inter_LN = IntersectionLattice_VectorSpace(Latt, Pplane):
C cap P on the same coordinate rows as the projection above. A
vector c * B of C (B the rows of a basis of C) lies in P exactly
when it pairs to zero with every edge wall, so the coefficient
lattice of C cap P is the integral left kernel of B G R^T.

Lines 1327–1329:

The Hermite normal form of a basis of the integral left kernel is the
echelon basis Sage's left_kernel() carries (the form is unique for
the kernel lattice).

### `_Edgewalk._anisotropic_edge_candidate_family`

Lines 1385–1387:

edge.form is scale * Gwork and the bound is scale * N, while
Anisotropic receives Gwork and res_norm itself (edgewalk.h:765-777);
pass the factor so Canonical's floor is taken at the source's scale.

Lines 1404–1414:

lorentzian_linalg.h:70-143,
GetMatrixExponentSublattice_TrivClass: first find the least
positive power preserving C cap Pi, then the order of its
action on Pi/(C cap Pi). The source tests integral membership
in the sublattice and then the induced finite quotient action,
not the generally larger order in GL(2, Z/exp(Pi/S)).
On coordinate rows of Pi the period acts by x -> x * A, A the rows
of the images of the basis (edgewalk.h:833, TheMat.transpose() * v).
With M the rows of a basis of S = C cap Pi, h preserves S exactly
when M h M^{-1} is integral, and acts trivially on Pi/S exactly
when (h - 1) M^{-1} is integral.

### `_Edgewalk.binary_candidates`

Lines 1450–1452:

The pinned fixed-norm factorization admits negative residual
squares as well as positive ones. Only the zero-square case lacks
a finite signed factor-pair enumeration.

Lines 1463–1468:

Pinned edgewalk.h:921-961 uses the *entire* signed factor-pair
enumeration, including imprimitive vectors, before ordering
roots and checking their ambient integral reconstructions.
A vector-orbit representative family is not that candidate set.
The source factors the form on integers (Factor_GP_LN); so does
binary_split_form_vectors_of_square, with no plane lattice.

### `_Edgewalk.extension_roots`

Lines 1498–1500:

edgewalk.h:905,988 get_next_anisotropic/get_next_isotropic:
SolutionIntMat(Latt, v_T), every integral candidate must
belong to the norm-specific lattice.

Lines 1506–1511:

edgewalk.h get_next_anisotropic/get_next_isotropic returns the
FIRST integral member of this ordered family. Later members
are not inspected, even if vertex_from_root subsequently
rejects the resulting endpoint. Preserve that stopping rule
rather than reconstructing the entire family for each Dynkin
extension (the family itself is cached by binary_candidates).

### `_Edgewalk.vertex_from_root`

Lines 1546–1549:

A Lorentzian chamber corner is a nonspacelike ray.  The root
congruence and side-pairing tests alone do not imply this: a
positive-norm intersection can satisfy both and is not a vertex
of the hyperbolic polyhedron.

Lines 1558–1561:

edgewalk.h:1032-1058 stores alpha itself in MatRoot. A scalar
multiple defines the same hyperplane, but not the same normed root
configuration; replacing alpha by its primitive part here changes
the reference algorithm's subsequent extensions and stabilizers.

### `_Edgewalk.cusp_root`

Lines 1581–1582:

IntersectionLattice_VectorSpace(Latt, {u, k}): Z^n cap span in the
coordinates of the basis of C.

Lines 1586–1586:

ComputePairGcd(k1, k2): a unimodular P with (k1, k2) P = (g, 0).

Lines 1590–1590:

P = [[s, -right], [t, left]] has determinant s*left + t*right = 1.

### `_Edgewalk.cusp_roots`

Lines 1630–1631:

edgewalk.h:393-408, DetermineRootsCuspidalCase, appends the
actual chosen vector, without primitive normalization.

### `_Edgewalk.isotropic_next_vertex`

Lines 1636–1638:

The two null lines of the split edge plane L cap P on engine rows:
P is the integral left kernel of G R^T (R the rows of the walls),
and each line is primitive in L cap P, which is saturated in L.

Lines 1640–1642:

A binary edge has no retained wall roots (l_ui = []). Preserve the
ambient column count explicitly: matrix([]) is 0 x 0, not the
0 x rank wall matrix whose transpose multiplies G.

Lines 1650–1654:

edgewalk.h:1094-1120, get_can_gen orients each isotropic
line by b(k, v) < 0 *before* testing the discarded wall.
Primitive-isotropic-line enumerators need not return the
representative in this time component; rejecting the wrong
sign loses the unique forward cusp altogether.

Lines 1657–1657:

In the split edge plane this is the old isotropic line.

### `_Edgewalk.edge_step`

Lines 1677–1679:

edgewalk.h:966-1061 calls get_next(extension) once: it takes
the first integral candidate in the binary search ordering.
Only that candidate is tested as a possible next vertex.

Lines 1689–1702:

With no finite root extension, an anisotropic rational edge plane
has no rational ideal endpoint either.  The chamber edge therefore
has no vertex in the projective closure, so the reflection chamber
cannot be a finite-volume fundamental polyhedron.  This is a
nonreflective terminal edge, not a failed cusp computation.

Substitution, recorded deliberately: the source has no output for
this configuration.  Its isotropic fallback returns an empty
generator list on an anisotropic edge plane and it aborts with
TerminalException{"Rethink needed"} (edgewalk.h:1117-1122), after
which no verdict exists.  The negative verdict comes from Allcock's
dichotomy instead (main.tex:1971-1985, 2082-2089): a finite-volume
chamber has every edge shared with a neighbouring chamber, so a
terminal anisotropic edge rules reflectivity out.

### `_Edgewalk.edge_directions`

Lines 1723–1732:

edgewalk.h:525-539 GetAdjacencyDirection: v_disc is the last
root outside the facet. An ideal vertex with several affine
components has facets that omit several walls (U+<2>+<2>
under K3 norms). Every omitted wall gives the same edge: the
roots vanish at k, so a point x = a k + c d of the edge plane
P = (included)^perp pairs with an omitted root r as
b(r, x) = c b(r, d), where d spans the edge ray and
b(r, d) < 0 because r is not a wall of that edge. The tests
in edge_plane, vertex_from_root and isotropic_next_vertex
read only this sign, which is -sign(c) for every omitted r.

### `_Edgewalk.coloured_family`

Lines 1746–1747:

The stabilizer and transporter read each value many times; each
point of L x {0, 1} is built on its first reading only.

### `_Edgewalk.root_chamber_certificate`

Lines 1832–1834:

DualDescriptionStandard of ListIneq = -LVect * G (edgewalk.h:1721-1726)
on the Gram matrix and the root coordinates: the cone
{x : -b(root, x) >= 0} and its extreme rays.

Lines 1843–1849:

The dual-description norm check alone is vacuous on the zero
cone.  A hyperbolic chamber must have full span and contain a
timelike vector.  Moreover, different nonspacelike rays must
belong to the same closed time cone.  Otherwise their positive
span contains spacelike vectors despite each ray having q <= 0.
Pairwise nonpositive inner products make the entire positive span
nonspacelike by expanding the square of any positive combination.

### `_Edgewalk.binary_walk`

Lines 1932–1936:

The DEFAULT isotropic_vinberg start is a null ray. In rank two
its perpendicular quotient is zero-dimensional, so the source's
get_simple_cone raises NonReflectivityException("list of roots
is empty"). A FileVertex/FileVertexRoots start takes the other
branch of get_initial_vertex and must not be intercepted here.

Lines 1942–1944:

As in the higher-rank walk, a bare vertex has its simple
cone computed while a (vertex, roots) pair is taken as
given and checked below (edgewalk.h:2100-2115).

Lines 1967–1979:

Choose primitive INTEGRAL generators e,f of the two saturated
rational null lines. An integral isometry maps each line to
itself or the other, so its matrix is diagonal or antidiagonal
with multipliers a,d in Z (ge lies in the target rank-one
saturated submodule), and preserving b(e,f) != 0 gives a*d = 1.
Hence a,d = +-1 and O(L) is FINITE. A
finite group acting on H^1 = R has quotient R or a half-line,
both of infinite length, so no reflection subgroup can be
reflective. This exact
arithmetic argument covers prescribed starts too; unlike the
default null start it does not pretend that their wall list was
empty. In particular an ideal neighbour has no wall root and
must not be unpacked as another reflective endpoint.

### `_Edgewalk.walk`

Lines 1996–2001:

edgewalk.h:2100-2143 FileVertex / FileVertexRoots. A bare
vertex has its simple cone computed by get_simple_cone,
including the cusp route for a null ray and the check that
the perpendicular roots are a full-rank corner (:1848-1873);
a (vertex, roots) pair is taken as given, exactly as the
source trusts FileVertexRoots (:2100-2110).

Lines 2016–2018:

The corner search can end in the source's own nonreflectivity
decision (edgewalk.h:1985-1999); report its reason rather than a
negative verdict with no mathematical justification.

Lines 2027–2031:

A rational row x is fixed by a chamber symmetry with row matrix M
precisely when x(M-I)=0. A finite cone-preserving group fixes a
timelike vector by averaging; conversely the integral stabilizer of
a timelike vector is discrete in a compact orthogonal group and is
therefore finite. Store the complete common fixed space exactly.

### `_Edgewalk.walk.certified_roots`

Lines 2058–2063:

LORENTZ_RunEdgewalkAlgorithm::f_try_terminate: close the
accumulated walls under every discovered chamber symmetry.
compute_full_root_orbit_iter (edgewalk.h:1216-1258): the
orbits of the simple roots under the generators, on the
coordinates of the roots and the integral matrices of the
isometries; each root is wrapped once per walk.

### `_Edgewalk.walk.try_reflectivity_certificate`

Lines 2088–2089:

LORENTZ_RunEdgewalkAlgorithm::f_try_terminate tests the dual
description only when the closed root orbit has grown.

### `_Edgewalk.walk`

Lines 2101–2102:

The source checks after insertion of the initial vertex and after
every newly discovered vertex, not only once the quotient closes.

Lines 2162–2165:

Allcock, Theorem (ComputesChamber): the complete corner walk
terminates exactly in finite volume. The closed quotient graph
and finite chamber-symmetry group unfold to finitely many
corners with every edge treated. No wall orbit is needed here.

## `witnesses/test_edgewalk_source_options.sage`

### `test_binary_file_vertex_and_file_vertex_roots`

Lines 32–33:

The displayed root is reflective on the integral lattice: reflection is
-1 on r=e0+e1 and +1 on k=e0-e1, so it maps (e0,e1) to (-e1,-e0).

Lines 42–45:

FileVertex and FileVertexRoots must not be confused with the default
isotropic_vinberg start. O(U) is finite (its null-line scalings are
integral units), so the exact nonreflective verdict has infinite
covolume, not a fictitious empty wall list at the prescribed finite k.

Lines 56–57:

A bare null FileVertex follows get_simple_cone, including ray
normalisation, and returns the source's cusp nonreflectivity report.

Lines 64–65:

The empty-wall isotropic fallback used above must not accept bad roots
merely because the ambient binary lattice is isotropic.

### `test_exhaustive_mode_skips_early_finiteness_and_preserves_chamber`

Lines 87–89:

The transporter is explicit, so this checks the ENTIRE integral root
multiset and its lattice gluing rather than merely a Coxeter graph.
It is a signed swap of the last two coordinates, with determinant -1.

### `test_anisotropic_cycle_period_and_input_guards`

Lines 107–108:

For q=2*x^2-y^2, P sends (e0,e1) to (3,4),(2,3). Direct evaluation
gives q(3,4)=2, q(2,3)=-1 and b((3,4),(2,3))=0; det(P)=1.

Lines 120–120:

split: discriminant one

Lines 122–122:

not Lorentzian
