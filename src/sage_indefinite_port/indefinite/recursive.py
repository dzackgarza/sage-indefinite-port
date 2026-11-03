r"""Dutour Sikirić's Eichler recursion over preamble operations (map §10, references/port-design.md).

For L of positive index at least 2 after the sign normalization of AttackProfile:

- generation step: if the A_L-translates of finitely many vectors r_i = t_i(v) cover the
  O(L)-orbit of a vector v, then O(L) = <A_L, Aut((L, v)), t_i>;
- Aut((L, v)): for q(v) != 0, the extensions across v^perp of generators of O(v^perp), then the
  subgroup of the group they generate preserving L; for q(v) = 0, the pointwise stabilizer of
  Zv in its parabolic subgroup;
- an element of Isom((L, v), (M, w)): for q(v) != 0, an isometry of the complements by
  recursion, its extension, and its correction by an element of the group above carrying the
  image of L onto M; for q(v) = 0, an integral lift of an isometry of the reductions, searched
  over the cosets of the Levi image;
- isometry(L, M): a splitting vector v of L, the covering vectors w of M of the same square
  and content, and an element of Isom((L, v), (M, w));
- primitive vector classes: covering vectors deduplicated by deciding Isom exactly.

Definite lattices are the preamble's; signature (1, n) is the Lorentzian leaf (#15).
"""

import itertools
from dataclasses import dataclass
from math import gcd

from sage_extensions.all import GObjects, Lattices, Modules, Sets, ZZ
from sage_extensions.categories.group.groups import OwnedGroups
from sage_extensions.categories.lattice_morphisms import LatticeIsometryMethods
from sage_extensions.categories.modules.orbit_spans import (
    OrbitSpanConstructions,
)
from sage_extensions.categories.modules.pure.modules import (
    Modules,
    ModulesOverIntegralDomains,
)
from sage_extensions.categories.rational_integral_stabilizers import (
    IntegralStructureAction,
)
from sage.groups.matrix_gps.finitely_generated import MatrixGroup
from sage.libs.gap.libgap import libgap
from sage.groups.matrix_gps.linear import GL
from sage.matrix.constructor import identity_matrix, matrix
from sage.matrix.special import diagonal_matrix
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.rings.integer_ring import ZZ as SageZZ
from sage.structure.parent import Parent

from sage_indefinite_port.indefinite.edgewalk_scalars import _integer
from sage_indefinite_port.readiness import UnfinishedCapability
from sage_indefinite_port.indefinite.eichler import OrbitCoverModel
from sage_indefinite_port.indefinite.groupoid_cache import IsometryGroupoidCache, _gram_in_label_order
from sage_indefinite_port.indefinite.indefinite_reduction import (
    indefinite_reduced_presentation,
)
from sage_indefinite_port.indefinite.lorentzian_cells import (
    lorentzian_isometry,
    lorentzian_orthogonal_group_generators,
)
from sage_indefinite_port.invariants import (
    AttackProfile,
    LatticePrefilter,
    VectorPrefilter,
)

_UNTWIST = Lattices(ZZ).twist_functor(-1)


def _untwisted(isometry: LatticeIsometryMethods, profile: AttackProfile) -> LatticeIsometryMethods:
    r"""An isometry of the signed view, as an isometry of the lattice itself."""
    return isometry if profile.sign == 1 else _UNTWIST(isometry)


def _extend_perpendicular_isometry(
    source_reduction,
    target_reduction,
    source_basis: tuple[Lattices.ElementMethods, ...],
    images: tuple[Lattices.ElementMethods, ...],
) -> tuple[LatticeIsometryMethods, int]:
    r"""Extend a supplied isometry of I^perp using the preamble Witt operation.

    source: Lorentzian.g:4-50 LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1
    (rank-one I, called from IndefiniteForms.g:1029 and :1638) and
    Lorentzian.g:115-249 LORENTZ_ExtendOrthogonalIsotropicIsomorphism
    get_one_transformation (rank-two I, IndefiniteForms.g:1097-1106, :1685-1686).
    The source's Subspace1/Subspace2 are a basis of N=I^perp and its images;
    here they are ``source_basis`` and ``images``. Its extension to L_Q -> M_Q, the Euclidean-nullspace
    complement (Fundamental.g:509-516), the pairing equations
    (SpecialEquationSolving, Lorentzian.g:53-99) and the prime-denominator
    elimination (Fundamental.g:751-794) belong to the preamble's
    ``witt_extension``. For rank-one I the extension is unique (an isometry
    fixing the hyperplane v^perp pointwise of a nondegenerate space with
    isotropic normal is trivial), so the choice of complement does not
    change the result. ``GetDenominatorMatrix`` (IndefiniteForms.g:1105) is
    the denominator ideal of the extension from L to M.
    """
    source = source_reduction.isotropic_embedding().codomain()
    target = target_reduction.isotropic_embedding().codomain()
    fraction_map = source.base_ring().fraction_field_map()
    field = fraction_map.codomain()
    rank = int(source_reduction.orthogonal_complement().inclusion().domain().module_rank())
    if len(source_basis) != rank or len(images) != rank:
        raise ValueError("the perpendicular partial isometry needs images of a full chosen basis")
    rational_source = source.base_change(fraction_map)
    rational_target = target.base_change(fraction_map)
    source_rationalize = ModulesOverIntegralDomains(source.base_ring())(source).generic_fibre_map()
    target_rationalize = ModulesOverIntegralDomains(target.base_ring())(target).generic_fibre_map()
    # The source's Subspace1 and Subspace2 are the rows of the two families.
    # Here they are the maps D -> L_Q and D -> M_Q out of the rational
    # lattice D on the basis, with the Gram pulled back along the first; the
    # partial isometry is the identity of D. No coordinate solve is needed.
    frame = field.free_module(rank)
    labels = tuple(frame.module_generating_set())
    gram = rational_source.gram_tensor().pullback(
        Modules(field).Mor(frame, rational_source)(dict(zip(labels, (source_rationalize(point).underlying_element() for point in source_basis), strict=True)))
    )
    span = Lattices(field)([[gram[row, column] for column in range(rank)] for row in range(rank)])
    span_labels = tuple(span.module_generating_set())
    source_inclusion = Modules(field).Mor(span, rational_source)(
        dict(zip(span_labels, (source_rationalize(point).underlying_element() for point in source_basis), strict=True))
    )
    target_inclusion = Modules(field).Mor(span, rational_target)(dict(zip(span_labels, (target_rationalize(point).underlying_element() for point in images), strict=True)))
    partial = span.Isom(span).identity()
    extension = partial.witt_extension(source_inclusion, target_inclusion)
    (denominator,) = tuple(extension.denominator_ideal(source, target).ideal_generators())
    return extension, int(denominator)


def _column_matrix(isometry: LatticeIsometryMethods) -> Matrix_integer_dense:
    r"""The integer matrix whose column ``j`` holds the coordinates of the image of generator ``j``."""
    domain, codomain = isometry.domain(), isometry.codomain()
    codomain_labels = tuple(codomain.module_generating_set())
    domain_labels = tuple(domain.module_generating_set())
    domain_basis = tuple(domain.module_generator(label) for label in domain_labels)
    images = tuple(isometry(generator).to_vector() for generator in domain_basis)
    return matrix(SageZZ, [[_integer(image(label)) for image in images] for label in codomain_labels])


def _from_column(lattice: Lattices.ParentMethods, column) -> Lattices.ElementMethods:
    r"""The element of ``lattice`` with the integer coordinates ``column``."""
    return lattice.linear_combination({label: int(entry) for label, entry in zip(lattice.module_generating_set(), column, strict=True) if entry})


def _isometry_from_row_matrix(source: Lattices.ParentMethods, target: Lattices.ParentMethods, rows: Matrix_integer_dense) -> LatticeIsometryMethods:
    r"""Build a morphism from image rows in coordinate-label order.

    The pinned source's degenerate extensions and block actions use row
    matrices, while the public constructor consumes images in module-generator
    order. Convert explicitly at this boundary.
    """
    source_labels = tuple(source.module_generating_set())
    target_labels = tuple(target.module_generating_set())
    if source.base_ring() != target.base_ring():
        raise ValueError("row isometry matrix requires a common base ring")
    if rows.nrows() != len(source_labels) or rows.ncols() != len(target_labels):
        raise ValueError("row isometry matrix must have shape (source rank, target rank)")
    source_gram = _gram_in_label_order(source)
    target_gram = _gram_in_label_order(target)
    if rows * target_gram * rows.transpose() != source_gram:
        raise ValueError("row matrix does not define an isometry between the supplied Gram forms")
    images = []
    for generator in source.module_generators():
        coordinates = generator.to_vector()
        source_coordinates = tuple(_integer(coordinates(label)) for label in source_labels)
        target_coordinates = tuple(
            sum(source_coordinates[row] * int(rows[row, column]) for row in range(len(source_labels)))
            for column in range(len(target_labels))
        )
        images.append(_from_column(target, target_coordinates))
    morphism_parent = source.O() if source is target else source.Isom(target)
    return morphism_parent(tuple(images))


def _isometry_from_column_matrix(source: Lattices.ParentMethods, target: Lattices.ParentMethods, columns: Matrix_integer_dense) -> LatticeIsometryMethods:
    r"""Build a morphism from a label-coordinate matrix at the owned API boundary.

    Matrix columns are indexed by ``source.module_generating_set()`` labels;
    the public Isom constructor consumes images in ``source.module_generators()``
    order. Translate explicitly between those conventions.
    """
    source_labels = tuple(source.module_generating_set())
    target_labels = tuple(target.module_generating_set())
    if source.base_ring() != target.base_ring():
        raise ValueError("column isometry matrix requires a common base ring")
    if columns.nrows() != len(target_labels) or columns.ncols() != len(source_labels):
        raise ValueError("column isometry matrix must have shape (target rank, source rank)")
    source_gram = _gram_in_label_order(source)
    target_gram = _gram_in_label_order(target)
    if columns.transpose() * target_gram * columns != source_gram:
        raise ValueError("column matrix does not define an isometry between the supplied Gram forms")
    images = []
    for generator in source.module_generators():
        coordinates = generator.to_vector()
        source_coordinates = tuple(_integer(coordinates(label)) for label in source_labels)
        target_coordinates = tuple(
            sum(int(columns[row, column]) * source_coordinates[column] for column in range(len(source_labels)))
            for row in range(len(target_labels))
        )
        images.append(_from_column(target, target_coordinates))
    morphism_parent = source.O() if source is target else source.Isom(target)
    return morphism_parent(tuple(images))


def _coordinate_box(lattice: Lattices.ParentMethods, height: int):
    r"""The integer coordinate tuples of a lattice with `max |entry| = height`, in a fixed order.

    Used only to walk candidate vectors for :meth:`EichlerRecursion.marked_orbit_representatives`. The order
    is by height and then lexicographically, so the walk is exhaustive on the box and never materialises more
    than one box at a time.
    """
    rank = int(lattice.rank())
    if height == 0:
        yield tuple(0 for _ in range(rank))
        return
    ranges = [range(-height, height + 1)] * rank
    for coordinates in itertools.product(*ranges):
        if max(abs(entry) for entry in coordinates) == height:
            yield coordinates


def _conjugated(transport: LatticeIsometryMethods, generators: tuple[LatticeIsometryMethods, ...]) -> tuple[LatticeIsometryMethods, ...]:
    r"""The generators ``h`` of O(R) carried to L by ``transport: R -> L``, each built once from ``P H P^-1``.

    source: IndefiniteForms.g:1487 ``Inverse(RecRed.B) * eGen1 * RecRed.B``.
    Each conjugate is constructed from the images of the generators of L, so
    later evaluations do not compose the three maps again.
    """
    lattice = transport.codomain()
    change = _column_matrix(transport)
    inverse = change.inverse().change_ring(SageZZ)
    return tuple(_isometry_from_column_matrix(lattice, lattice, change * _column_matrix(generator) * inverse) for generator in generators)


def _transported(target_transport: LatticeIsometryMethods, witness: LatticeIsometryMethods, source_transport: LatticeIsometryMethods) -> LatticeIsometryMethods:
    r"""``target_transport * witness * source_transport^-1``, built once from its column matrix.

    source: IndefiniteForms.g:1389 ``Inverse(RecRed2.B) * test * RecRed1.B`` and
    :1419 (the degenerate extension ``Inverse(FullBasis2) * TheEquivTest * FullBasis1``).
    """
    source, target = source_transport.codomain(), target_transport.codomain()
    columns = _column_matrix(target_transport) * _column_matrix(witness) * _column_matrix(source_transport).inverse().change_ring(SageZZ)
    return _isometry_from_column_matrix(source, target, columns)


@dataclass(frozen=True)
class VectorClass:
    r"""An O(L)-orbit of primitive vectors: a representative, and the isometries of L found carrying it to the other covering vectors of its orbit."""

    representative: Lattices.ElementMethods
    witnesses: tuple[LatticeIsometryMethods, ...]


@dataclass(frozen=True)
class MarkedOrbitCanonical:
    r"""The canonical representative of a pointed orbit, with a transporter back to the input.

    ``representative`` is the least element, in the declared coordinate order, of the orbit of the
    canonicalized vector among the lattice's elements; ``transporter`` is the actual isometry with
    ``transporter(representative) == vector``; ``height`` is ``max |coordinate|`` of the representative.

    The representative is canonical **relative to the declared basis and coordinate order**, not a
    functorial object: a different basis presents the same orbit relation with different minima. The orbit
    relation itself is untouched.
    """

    representative: Lattices.ElementMethods
    transporter: LatticeIsometryMethods
    height: int


@dataclass(frozen=True)
class MarkedOrbitDecision:
    r"""One candidate's verdict in the duplicate-free orbit enumeration (issue #62, part B).

    ``candidate`` is the vector under test, ``emitted`` says whether it is the first member of its orbit in
    the enumeration order -- and therefore that orbit's canonical representative -- and ``representative``
    with ``transporter`` carries the verdict's witness: the canonical representative of the candidate's orbit
    and the isometry with ``transporter(representative) == candidate``. For an emitted candidate the
    representative is the candidate itself and the transporter is the identity.

    A comparison that cannot be decided is never read as inequivalence: ``UnfinishedCapability`` raised by
    the equivalence backend propagates out of the enumeration, because turning "undecided" into "distinct"
    would emit one class twice.
    """

    candidate: Lattices.ElementMethods
    emitted: bool
    representative: Lattices.ElementMethods
    transporter: LatticeIsometryMethods


@dataclass(frozen=True)
class RadicalCongruenceImage:
    r"""The **actual** image of `GL_k(Z)` in `GL_k(Z/modulus)`, as a maintained matrix group plus its lifts.

    `group` is Sage's GAP-backed `MatrixGroup` over `Zmod(modulus)`, generated by the reductions of the integral
    matrices in `lifts`; each lift is integral, of determinant `+-1`, and reduces to its generator. Order,
    membership, orbits, transporters, stabilizers and the word section below are all the maintained GAP operations
    reached through `libgap` -- nothing here reimplements matrix arithmetic over `Zmod`, group closure or an order
    formula.

    **The group is not `GL_k(Z/modulus)` in general.** A lift is integral and unimodular, so its determinant is
    `+-1` modulo the modulus, and the image is exactly `det^{-1}({+-1})`, of index `phi(modulus)/2`. That is proper
    as soon as `phi(modulus) > 2` -- `modulus = 5` already has index two, and `diag(2, 1) mod 5` has determinant
    `2` and no integral unimodular lift at all. The image is the whole group exactly for `modulus` in
    `{1, 2, 3, 4, 6}`, where `phi(modulus) <= 2`, which is why the distinction is invisible on small moduli.

    **Vector orbits are a different question from surjectivity of the group map.** The image is proper at rank two
    and modulus five, and the specimens below still transport `(5,1,0) -> (5,2,0)`; the carrier is an explicit
    integral `D`, and no transitivity is read off the group map.

    Rows multiply on the right in this port, so `D` acts by `r |-> r D`; the corresponding GAP action is `OnRight`
    on row vectors over the residue ring, which is what `orbit_of`/`transporter`/`stabilizer_order` use.
    """

    modulus: int
    radical_rank: int
    group: object
    generators: tuple[Matrix_integer_dense, ...]
    lifts: tuple[Matrix_integer_dense, ...]
    general_linear_order: int

    @property
    def order(self) -> int:
        r"""`|im|`, from the maintained matrix group."""
        return int(self.group.order())

    @property
    def is_proper(self) -> bool:
        r"""Whether the image is a proper subgroup of `GL_k(Z/modulus)`."""
        return self.order < self.general_linear_order

    @property
    def index(self) -> int:
        return self.general_linear_order // self.order

    @property
    def determinant_image(self) -> set[int]:
        r"""The determinants of the **generators**, modulo the modulus; multiplicativity carries it to the group."""
        return {int(generator.det()) % self.modulus for generator in self.generators}

    def contains(self, residue) -> bool:
        r"""Membership of a `k x k` residue matrix in the image, decided by the maintained group.

        This is `matrix in group`; the "shadow" arithmetic an earlier revision used is gone.
        """
        entries = [[int(residue[row, column]) % self.modulus for column in range(self.radical_rank)]
                   for row in range(self.radical_rank)]
        return matrix(Zmod(self.modulus), entries) in self.group

    def _point(self, row):
        r"""The residue row as a GAP vector over the ring the matrices live in, so the action is well defined.

        The ring is taken from the group's own coefficients: `ZmodnZ(m)`, which is `GF(p)` when the modulus is
        prime. Building points in that ring -- rather than as integer lists -- is what makes GAP's `OnRight`
        action well defined for these matrices.
        """
        ring = libgap.ZmodnZ(self.modulus)
        one = libgap.One(ring)
        entries = tuple(int(entry) for entry in row)
        if len(entries) != self.radical_rank:
            raise ValueError(f"a radical row must have {self.radical_rank} entries, got {len(entries)}")
        return libgap([libgap(entry % self.modulus) * one for entry in entries])

    def _orbit_gap(self, row):
        point = self._point(row)
        return point, libgap.Orbit(self.group.gap(), point, libgap.OnRight)

    def orbit_of(self, row) -> tuple[tuple[int, ...], ...]:
        r"""The orbit of a residue row under the image, as integer residues -- GAP's `Orbit`."""
        _, orbit = self._orbit_gap(row)
        return tuple(tuple(int(libgap.Int(entry)) for entry in vector_) for vector_ in orbit)

    def stabilizer_order(self, row) -> int:
        r"""`|{D in im : r D = r}|`, from GAP's `Stabilizer` -- the congruence part of a pointwise stabilizer."""
        point, _ = self._orbit_gap(row)
        stabilizer = libgap.Stabilizer(self.group.gap(), point, libgap.OnRight)
        return int(libgap.Size(stabilizer))

    def stabilizer_generators(self, row):
        r"""Integral lifts of generators of the **residue** stabilizer `{D in im : r D = r (mod modulus)}`.

        Each generator of GAP's `Stabilizer` subgroup of the image is pulled back to a word in the integral
        generating set, composed from `lifts`, and verified to be unimodular and to fix the residue -- the same
        lift-word obligation `transporter` discharges, on the same row convention.

        **Scope, exactly.** These matrices generate the residue stabilizer, i.e. the image of the stabilizer in
        `GL_k(Z/modulus)`. They do **not** generate the integral inverse image
        `{D in GL_k(Z) : r D = r (mod modulus)}`, which is infinite and contains the congruence kernel
        `{D : D = 1 (mod modulus), r D = r}`. Agreement with `stabilizer_order` modulo the modulus is precisely a
        statement about that finite image and says nothing about the kernel, so no downstream claim of full
        integral stabilizer generation may rest on this method. Generating the integral inverse image needs the
        kernel as well, and that is a separate obligation which is **not** supplied here.
        """
        point, _ = self._orbit_gap(row)
        group = self.group.gap()
        stabilizer = libgap.Stabilizer(group, point, libgap.OnRight)
        epimorphism = libgap.EpimorphismFromFreeGroup(self.group)
        residue = tuple(int(entry) % self.modulus for entry in row)
        lifts = []
        for generator in libgap.GeneratorsOfGroup(stabilizer):
            word = libgap.PreImagesRepresentative(epimorphism, generator)
            syllables = [int(entry) for entry in libgap.ExtRepOfObj(word)]
            lift = identity_matrix(SageZZ, self.radical_rank)
            for generator_index, exponent in zip(syllables[0::2], syllables[1::2], strict=True):
                lift = lift * self.lifts[generator_index - 1] ** exponent
            if abs(int(lift.det())) != 1:
                raise ArithmeticError("a stabilizer lift is not unimodular")
            image = tuple(
                int(entry) % self.modulus
                for entry in (matrix(SageZZ, [list(residue)]) * lift).list()
            )
            if image != residue:
                raise ArithmeticError(
                    f"a stabilizer lift carries {residue} to {image}: the word order and the lift order disagree"
                )
            lifts.append(lift)
        return tuple(lifts)

    def transporter(self, source_row, target_row):
        r"""An integral unimodular `D` with `source_row D = target_row` modulo the modulus, or `None`.

        The chain is GAP's documented one, with every argument in the ring of the matrices: orbit the source point,
        form the `ActionHomomorphism` on that orbit, ask `RepresentativeAction` for a permutation carrying the
        source position to the target's, pull it back to an element of the matrix group, and pull *that* back to a
        word in the generators via `EpimorphismFromFreeGroup`. The integral lift is then the same word composed
        from `lifts` -- the lift-word obligation, which is the port's own data -- and it is **verified** to be
        unimodular and to carry the source residue onto the target.

        `None` is a decision, not a failure: the orbit is exhaustive, so a target outside the source's orbit has no
        carrier in the image, and none in `GL_k(Z)` either.
        """
        point, orbit = self._orbit_gap(source_row)
        target = self._point(target_row)
        target_position = _gap_position(orbit, target)
        if target_position is None:
            # Exhaustive nonmembership: the orbit list is complete, so no element of the image carries the source
            # residue onto the target -- and none of `GL_k(Z)` does either.
            return None
        source_position = _gap_position(orbit, point)
        if source_position is None:
            raise ArithmeticError("the source residue is not a row of the module the image acts on")
        group = self.group.gap()
        action = libgap.ActionHomomorphism(group, orbit, libgap.OnRight)
        permutation = libgap.RepresentativeAction(
            libgap.Image(action), source_position, target_position, libgap.OnPoints
        )
        if permutation == libgap.fail:
            # Both positions exist and the image acts transitively on its orbit, so this is an unfinished
            # execution, not a decision: it must never be reported as "no carrier".
            raise ArithmeticError(
                "RepresentativeAction returned fail although both residues lie in the same orbit; the transport "
                "is undecided rather than absent"
            )
        element = libgap.PreImagesRepresentative(action, permutation)
        word = libgap.PreImagesRepresentative(libgap.EpimorphismFromFreeGroup(self.group), element)
        syllables = [int(entry) for entry in libgap.ExtRepOfObj(word)]
        lift = identity_matrix(SageZZ, self.radical_rank)
        for generator_index, exponent in zip(syllables[0::2], syllables[1::2], strict=True):
            lift = lift * self.lifts[generator_index - 1] ** exponent
        if abs(int(lift.det())) != 1:
            raise ArithmeticError("the assembled congruence lift is not unimodular")
        moved = tuple(
            int(entry) % self.modulus
            for entry in (matrix(SageZZ, [[int(e) for e in source_row]]) * lift).list()
        )
        expected = tuple(int(entry) % self.modulus for entry in target_row)
        if moved != expected:
            raise ArithmeticError(
                f"the assembled lift moved {tuple(int(e) for e in source_row)} to {moved}, not {expected}: the "
                "free-group letter order and the lift order disagree"
            )
        return lift


def _gap_position(sequence, element):
    r"""GAP's `Position` with its real sentinel: `fail` when absent, an integer when present.

    GAP returns `fail`, not `false`, for an absent entry, and libgap hands that back as a `GapElement_Boolean`
    whose string is `'fail'`. Comparing `str(value) == "false"` is therefore never true and silently leaves the
    sentinel in place, which is what an earlier revision did: the missing-target case then fell through to
    `RepresentativeAction` with a `fail` argument, whose own `fail` return was read as "no carrier". Those are
    different facts -- an exhaustive nonmembership versus an unfinished execution -- so they are separated here.
    """
    position = libgap.Position(sequence, element)
    if position == libgap.fail:
        return None
    return int(position)


def _coordinate_height(lattice: Lattices.ParentMethods, vector: Lattices.ElementMethods) -> int:
    """The height of a vector: ``max |coordinate|`` in the lattice's declared label order."""
    coordinates = vector.to_vector()
    return max((abs(_integer(coordinates(label))) for label in lattice.module_generating_set()), default=0)


def _coordinate_element(lattice: Lattices.ParentMethods, coordinates) -> Lattices.ElementMethods:
    """The lattice element with the given coordinates in declared label order."""
    entries = {
        label: int(entry)
        for label, entry in zip(lattice.module_generating_set(), coordinates, strict=True)
        if entry
    }
    return lattice.linear_combination(entries) if entries else lattice.zero()


class EichlerRecursion:
    r"""The recursion; every general step is a preamble operation."""

    def __init__(self) -> None:
        # IndefiniteForms.g:1173-1310: the groupoid database is shared by
        # every recursive orthogonal-group and lattice-equivalence request.
        self._groupoid = IsometryGroupoidCache()
        self._models: dict[Lattices.ParentMethods, OrbitCoverModel] = {}

    def _degenerate_group_generators(self, lattice: Lattices.ParentMethods) -> tuple[LatticeIsometryMethods, ...]:
        r"""IndefiniteForms.g:853–901, 1496–1511: extend O(C) over the radical.

        The presentation of ``lattice`` is ``C ⊕ R`` with Gram ``C ⊕ 0``.
        Every integral isometry has block rows ``[A B; 0 D]``, where
        ``A`` preserves C, ``D`` belongs to GL(R), and B is arbitrary.
        The source extends the generators of O(C) and GL(R) and adjoining
        the first-radical-coordinate shears e_i -> e_i + z_1: the GL(R)
        action supplies the other coordinates of B.
        """
        gram = _gram_in_label_order(lattice)
        size, nonsingular_rank = gram.nrows(), gram.rank()
        assert nonsingular_rank < size
        if any(gram[i, j] for i in range(nonsingular_rank, size) for j in range(size)):
            raise ArithmeticError("the radical must occupy the final Gram coordinates")
        matrices = []
        if nonsingular_rank:
            complement_gram = gram.matrix_from_rows_and_columns(range(nonsingular_rank), range(nonsingular_rank))
            complement = lattice.lattice_category()([[int(entry) for entry in row] for row in complement_gram.rows()])
            complement_labels = tuple(complement.module_generating_set())
            complement_basis = tuple(complement.module_generator(label) for label in complement_labels)
            for generator in self.orthogonal_group_generators(complement):
                # `gram` is expressed in module_generating_set() order. Build
                # the extended block in that same order rather than assuming
                # module_generators() enumerates the identical basis.
                top = tuple(
                    tuple(int(generator(point).to_vector()(label)) for label in complement_labels)
                    for point in complement_basis
                )
                block = identity_matrix(SageZZ, size)
                for row in range(nonsingular_rank):
                    for column in range(nonsingular_rank):
                        block[row, column] = top[row][column]
                matrices.append(block)
        radical_rank = size - nonsingular_rank
        general_linear = OwnedGroups().GL(radical_rank, lattice.base_ring())
        for generator in general_linear.group_generators():
            action = general_linear._to_engine(generator).matrix()
            block = identity_matrix(SageZZ, size)
            for i in range(radical_rank):
                for j in range(radical_rank):
                    block[nonsingular_rank + i, nonsingular_rank + j] = SageZZ(action[i, j])
            matrices.append(block)
        for row in range(nonsingular_rank):
            block = identity_matrix(SageZZ, size)
            block[row, nonsingular_rank] = 1
            matrices.append(block)
        return tuple(_isometry_from_row_matrix(lattice, lattice, block) for block in matrices)

    def _degenerate_isometry(self, source: Lattices.ParentMethods, target: Lattices.ParentMethods) -> LatticeIsometryMethods | None:
        r"""IndefiniteForms.g:1387–1423: split radicals, identify them, and extend the nondegenerate witness."""
        source_reduction = indefinite_reduced_presentation(source)
        target_reduction = indefinite_reduced_presentation(target)
        left, right = source_reduction.reduced, target_reduction.reduced
        left_gram = _gram_in_label_order(left)
        right_gram = _gram_in_label_order(right)
        size, rank = left_gram.nrows(), left_gram.rank()
        if rank != right_gram.rank():
            return None
        cacheable = left.base_ring() == SageZZ and right.base_ring() == SageZZ
        if cacheable:
            known, cached = self._groupoid.cached_isometry(left, right)
            if known:
                return None if cached is None else _transported(target_reduction.isometry, cached, source_reduction.isometry)
        extension = identity_matrix(SageZZ, size)
        if rank:
            left_complement_gram = left_gram.matrix_from_rows_and_columns(range(rank), range(rank))
            right_complement_gram = right_gram.matrix_from_rows_and_columns(range(rank), range(rank))
            left_complement = left.lattice_category()([[int(entry) for entry in row] for row in left_complement_gram.rows()])
            right_complement = right.lattice_category()([[int(entry) for entry in row] for row in right_complement_gram.rows()])
            complement_witness = self.isometry(left_complement, right_complement)
            if complement_witness is None:
                if cacheable:
                    self._groupoid.record_isometry(left, right, None)
                return None
            left_labels = tuple(left_complement.module_generating_set())
            right_labels = tuple(right_complement.module_generating_set())
            left_basis = tuple(left_complement.module_generator(label) for label in left_labels)
            for row, point in enumerate(left_basis):
                image_coordinates = complement_witness(point).to_vector()
                for column, label in enumerate(right_labels):
                    extension[row, column] = int(image_coordinates(label))
        reduced_witness = _isometry_from_row_matrix(left, right, extension)
        if cacheable:
            self._groupoid.record_isometry(left, right, reduced_witness)
        return _transported(target_reduction.isometry, reduced_witness, source_reduction.isometry)

    def marked_orbit_action(self, lattice: Lattices.ParentMethods):
        r"""The action of `O(L)` on the elements of `L` by evaluation, as the retained G-set.

        Built through the provider's own action constructor rather than from the group's attribute
        surface: `GObjects(G, Sets()).on_set(points, action, orbit_relation)` is the supported
        construction for an action on a set, and it yields the native `_GSetOnPoints` with `orbits()`,
        `stabilizer()` and `in_same_orbit()`.

        The orbit relation is supplied by :meth:`pointed_isometry`, which decides whether an isometry
        carries one vector onto the other (and is exact for degenerate `L` as well). This matters for
        issue #41/#53: `points()` of an orbit is a predicate-cut set, so it is never materialized and an
        infinite orbit stays a mathematical set instead of becoming a truncated list.
        """
        return GObjects(lattice.O(), Sets()).on_set(
            lattice,
            lambda isometry, point: isometry(point),
            orbit_relation=lambda left, right: self.pointed_isometry(left, right) is not None,
        )

    def marked_orbit(self, vector: Lattices.ElementMethods):
        r"""The `O(L)`-orbit of a marked vector, as the native lazy set.

        The result is the provider's `condition_set` cut out by the orbit relation. It supports membership
        decisions and representative transport, and it makes no finiteness claim: asking it `is_finite()`
        returns a proposition, not `True`. Consumers must therefore never coerce it to a finite ordered
        set -- `L[n]/O(L)` may be infinite, and a list would be a truncation.
        """
        return self.marked_orbit_action(vector.parent()).orbits().orbit_of(vector).points()

    def marked_square_slice(self, lattice: Lattices.ParentMethods, square):
        r"""The set `L[n] = {v in L : q(v) = n}` as the provider's predicate-cut set.

        This is a mathematical set, not a list: it is cut out by a predicate on the lattice, so membership
        is decided exactly and the set is never materialized.
        """
        target = int(square)
        return Sets().condition_set(lattice, lambda candidate: int(candidate.q()) == target)

    def marked_orbit_quotient(self, lattice: Lattices.ParentMethods, square):
        r"""`L[n] / O(L)`: the quotient set of orbit classes of square-`n` vectors (issues #41 and #53).

        This is the object the orbit consumer needs, and it is *not* the points of a single orbit:
        `action.orbits()` returns the provider's `_OrbitSet`, whose points are orbit classes, obtained
        through `orbit_of`. It is lazy: nothing is enumerated until asked, `orbit_of` returns the class of
        one vector, and `is_finite()` stays a proposition -- `L[n]/O(L)` may be infinite (see #53's
        example), so neither a finite list nor a finiteness claim is made here.
        """
        action = GObjects(lattice.O(), Sets()).on_set(
            self.marked_square_slice(lattice, square),
            lambda isometry, point: isometry(point),
            orbit_relation=lambda left, right: self.pointed_isometry(left, right) is not None,
        )
        return action.orbits()

    def marked_orbit_canonical(self, vector: Lattices.ElementMethods) -> MarkedOrbitCanonical:
        r"""The canonical representative of `vector`'s pointed orbit, with its transporter (issue #62, part C).

        ``can_Gamma(x) = min_{y ~ x} y`` in the declared coordinate order: height ``max |coordinate|`` first,
        then lexicographic within one height. The search walks only the **finite initial segment** of that
        order up to `vector`'s own height and returns the first `y` with an isometry carrying `y` onto
        `vector`; it terminates because `vector` lies in its own box, and the first hit is the global minimum
        of the orbit, since any smaller equivalent member would be inside the segment.

        Three properties follow from minimality and are what orbit consumers use:

        * stability, ``x ~ x'`` implies ``can(x) = can(x')``;
        * idempotence, ``can(can(x)) = can(x)``;
        * ``can(x)`` together with its transporter is exactly the orbit-class equality check: two vectors are
          in one orbit iff their canonical representatives are equal as lattice elements.

        The equivalence decided here is :meth:`pointed_isometry`. On the ``#41`` domain this is exact and
        terminating for degenerate `L` too, which is what makes the search a decision procedure rather than a
        search that might give up. It is **not** total everywhere: a backend that cannot decide raises
        ``UnfinishedCapability``, and that propagates -- an undecided comparison is never read as
        inequivalence, because that would split one orbit into two.
        """
        lattice = vector.parent()
        target_height = _coordinate_height(lattice, vector)
        for height in range(target_height + 1):
            for coordinates in _coordinate_box(lattice, height):
                candidate = _coordinate_element(lattice, coordinates)
                transporter = self.pointed_isometry(candidate, vector)
                if transporter is None:
                    continue
                return MarkedOrbitCanonical(candidate, transporter, height)
        raise ArithmeticError(
            "the coordinate box of a vector's own height contains the vector, so the canonical search cannot "
            "fail to terminate; reaching this point means the equivalence backend denied a vector to itself"
        )

    def marked_orbit_is_canonical(self, vector: Lattices.ElementMethods) -> bool:
        """Whether `vector` is the canonical representative of its own pointed orbit."""
        canonical = self.marked_orbit_canonical(vector)
        return canonical.representative == vector

    def marked_orbit_enumeration(self, lattice: Lattices.ParentMethods, square, limit: int | None = None, height_bound: int | None = None):
        r"""Exhaustive, duplicate-free enumeration of the orbits of square-`square` vectors, with witnesses.

        Each amplification step walks the next coordinate box of :func:`_coordinate_box` and decides every
        square-`square` candidate against the representatives already emitted. A candidate is **emitted** iff
        it is equivalent to none of them, and otherwise yields its witness -- the isometry carrying an earlier
        representative onto it. Since candidates are reached in the declared order, an emitted vector is the
        first member of its orbit to be reached, i.e. that orbit's canonical representative; two emitted
        vectors therefore cannot be equivalent (both would be minimal, hence equal). That is the
        duplicate-free and exhaustive argument of issue #62 part B, and it needs no finiteness of the
        quotient.

        What this does **not** do is signal that a finite quotient has ended. The candidate walk is unbounded,
        so a quotient whose classes are exhausted stops producing representatives while the search continues:
        absence of further output is not a completion certificate. `limit` bounds the number of
        representatives yielded, not the search.

        ``height_bound`` restricts the **search** to candidates of height at most that value and is issue #62 part
        D's bounded-height classification: one representative is emitted per orbit *meeting* the bound. The
        equivalence test stays global -- a transporter may carry a vector out of the box, which is why classes
        are never split by the bound. With no bound the search is unbounded, as above.

        ``UnfinishedCapability`` from the equivalence backend propagates: an undecided pair never becomes two
        classes.
        """
        target = int(square)
        emitted: list[Lattices.ElementMethods] = []
        height = 0
        while (limit is None or len(emitted) < limit) and (height_bound is None or height <= int(height_bound)):
            for coordinates in _coordinate_box(lattice, height):
                candidate = _coordinate_element(lattice, coordinates)
                if int(candidate.q()) != target:
                    continue
                witness = None
                for representative in emitted:
                    witness = self.pointed_isometry(representative, candidate)
                    if witness is not None:
                        break
                if witness is not None:
                    yield MarkedOrbitDecision(
                        candidate=candidate,
                        emitted=False,
                        representative=representative,
                        transporter=witness,
                    )
                    continue
                emitted.append(candidate)
                yield MarkedOrbitDecision(
                    candidate=candidate,
                    emitted=True,
                    representative=candidate,
                    # An emitted candidate is carried onto itself by the identity; `pointed_isometry(v, v)`
                    # returns exactly that, through the provider's own identity constructor.
                    transporter=self.pointed_isometry(candidate, candidate),
                )
                if limit is not None and len(emitted) >= limit:
                    return
            height += 1

    def marked_orbit_representatives(self, lattice: Lattices.ParentMethods, square, limit: int | None = None, height_bound: int | None = None):
        r"""One canonical representative per `O(L)`-orbit among the square-`square` vectors.

        The emitted vectors are the ones :meth:`marked_orbit_enumeration` marks ``emitted``: the first member
        of each orbit in the declared coordinate order, hence canonical relative to that order and
        duplicate-free. The provider's own ``_OrbitSet`` refuses to iterate unless the underlying G-set is
        placed finite, which is the right discipline but not an enumeration; this walks candidates instead.

        Nothing is materialized beyond the representatives already yielded, and `limit` is a **maximum on the
        number of representatives yielded, not a termination guarantee for the search**. A square whose
        classes are exhausted keeps searching for a next representative and may never find one, so this is not
        an infinite stream merely because the candidate walk is infinite: new classes appear as candidates are
        reached, and no finiteness or completeness claim is made.
        """
        for decision in self.marked_orbit_enumeration(lattice, square, limit, height_bound):
            if decision.emitted:
                yield decision.candidate

    def shear_kernel(self, complement: Lattices.ParentMethods, radical: Lattices.ParentMethods, marked):
        r"""The kernel of `Hom_Z(C,R) -> R`, `B |-> B x`, through the provider's own contracts.

        This is the shear kernel of issue #41: for `v = (x, r)` the pointwise stabilizer contains every
        `B` with `B x = 0`, and because `R` is torsion-free that is the same as `B x_0 = 0` for the
        primitive direction `x_0` of `x`. It is built from the maintained operations only -- the Hom module
        `Modules(Z).Mor(C, R)`, its `module_generating_set`/`module_generator`, `evaluation`, the morphism
        constructor, `kernel()` and the kernel's `inclusion` -- so no coordinate matrix is manufactured
        here; the provider's `_kernel_spanning_family` is what decides the linear algebra.

        Returns the kernel's elements as morphisms `C -> R` (the images of the kernel's generators under its
        inclusion), which callers evaluate on the generators of `C`.
        """
        hom = Modules(ZZ).Mor(complement, radical)
        labels = tuple(hom.module_generating_set())
        if not labels:
            raise ValueError("Hom(C,R) has no generators to evaluate")
        images = {
            label: hom.evaluation(hom.module_generator(label), marked)
            for label in labels
        }
        evaluation = Modules(ZZ).Mor(hom, radical)(images)
        kernel = evaluation.kernel()
        inclusion = kernel.inclusion()
        return tuple(inclusion(kernel.module_generator(label)) for label in kernel.module_generating_set())

    def orthogonal_group_generators(self, lattice: Lattices.ParentMethods) -> tuple[LatticeIsometryMethods, ...]:
        if lattice.is_definite():
            # IndefiniteForms.g:1275-1310 INDEF_FORM_AutomorphismGroup_Memoize
            # covers the definite kernel (:1435) too: a group computed
            # for an isometric lattice is transported along the recorded
            # isometry instead of being recomputed.
            cacheable = lattice.base_ring() == ZZ
            if cacheable:
                cached = self._groupoid.cached_group(lattice)
                if cached is not None:
                    return cached
            definite_generators = tuple(lattice.O().group_generators())
            if cacheable:
                self._groupoid.record_group(lattice, definite_generators)
            return definite_generators
        # IndefiniteForms.g:1474-1494: reduce before the memoized group
        # computation, then conjugate its generators back to the input form.
        # The reduction's isometry is directed from the reduced lattice to L.
        reduction = indefinite_reduced_presentation(lattice)
        cacheable = reduction.reduced.base_ring() == ZZ
        if cacheable:
            cached = self._groupoid.cached_group(reduction.reduced)
            if cached is not None:
                return _conjugated(reduction.isometry, cached)
        if not reduction.reduced.is_nondegenerate():
            # IndefiniteForms.g:1496–1511 applies ExtendIsometryGroup to
            # the nondegenerate quotient before transporting the generators.
            reduced_generators = self._degenerate_group_generators(reduction.reduced)
        else:
            profile = AttackProfile.from_lattice(reduction.reduced)
            match profile.positive_index:
                case 1:
                    generators = lorentzian_orthogonal_group_generators(profile.signed_view)
                case _:
                    generators = self._generation_step(profile.signed_view)
            reduced_generators = tuple(_untwisted(generator, profile) for generator in generators)
        if cacheable:
            self._groupoid.record_group(reduction.reduced, reduced_generators)
        return _conjugated(reduction.isometry, reduced_generators)

    def _model(self, lattice: Lattices.ParentMethods) -> OrbitCoverModel:
        r"""The approximate model of ``lattice``, built once per lattice (IndefiniteForms.g:1338, :1346, :1450, :1530)."""
        if lattice not in self._models:
            self._models[lattice] = OrbitCoverModel(lattice)
        return self._models[lattice]

    def _generation_step(self, lattice: Lattices.ParentMethods) -> tuple[LatticeIsometryMethods, ...]:
        model = self._model(lattice)
        vector = model.splitting_vector()
        transporters = tuple(
            isometry for candidate in model.covering_vectors(int(vector.q()), int(vector.content())) if (isometry := self.pointed_isometry(vector, candidate)) is not None
        )
        return model.approximate_generators() + self.pointed_automorphism_generators(vector) + transporters

    def isometry(self, source: Lattices.ParentMethods, target: Lattices.ParentMethods) -> LatticeIsometryMethods | None:
        # IndefiniteForms.g:1387–1423 splits the radical *before* the
        # nondegenerate genus/discriminant prefilters and memoized recursion.
        if source.is_nondegenerate() != target.is_nondegenerate():
            return None
        if not source.is_nondegenerate():
            if source.module_rank() != target.module_rank():
                return None
            return self._degenerate_isometry(source, target)
        if LatticePrefilter.from_lattice(source) != LatticePrefilter.from_lattice(target):
            return None
        if source.is_even() and source.genus() != target.genus():
            return None
        if not source.is_definite():
            source_reduction = indefinite_reduced_presentation(source)
            target_reduction = indefinite_reduced_presentation(target)
            reduced_witness = self._indefinite_isometry(source_reduction.reduced, target_reduction.reduced)
            if reduced_witness is None:
                return None
            return _transported(target_reduction.isometry, reduced_witness, source_reduction.isometry)
        # IndefiniteForms.g:1239-1273 INDEF_FORM_TestEquivalence_Memoize
        # records every exact verdict, the definite kernel's (:1326)
        # included, so later O(K) requests are transported along it.
        cacheable = source.base_ring() == ZZ and target.base_ring() == ZZ
        if cacheable:
            known, cached = self._groupoid.cached_isometry(source, target)
            if known:
                return cached
        homset = source.Isom(target)
        witness = None if homset.is_empty() else homset.an_element()
        if cacheable:
            self._groupoid.record_isometry(source, target, witness)
        return witness

    def _indefinite_isometry(self, source: Lattices.ParentMethods, target: Lattices.ParentMethods) -> LatticeIsometryMethods | None:
        cacheable = source.base_ring() == ZZ and target.base_ring() == ZZ
        if cacheable:
            known, cached = self._groupoid.cached_isometry(source, target)
            if known:
                return cached
        source_profile = AttackProfile.from_lattice(source)
        target_profile = AttackProfile.from_lattice(target)
        signed_source, signed_target = source_profile.signed_view, target_profile.signed_view
        match source_profile.positive_index:
            case 1:
                witness = lorentzian_isometry(signed_source, signed_target)
            case _:
                source_model = self._model(signed_source)
                target_model = self._model(signed_target)
                vector = source_model.splitting_vector()
                square, content = int(vector.q()), int(vector.content())
                # IndefiniteForms.g:1320-1359 and
                # ApproximateModels.h:GetCoveringOrbitRepresentatives: enumerate
                # the full target covering family at the source's first norm,
                # testing its pointed isometries in order. A separate
                # GetOneOrbitRepresentative trial is not part of this kernel:
                # it can demand an additional exact pointed equivalence and
                # does not replace any covering case.
                covering = target_model.covering_vectors(square, content)
                witness = next(
                    (isometry for candidate in covering if (isometry := self.pointed_isometry(vector, candidate)) is not None),
                    None,
                )
        result = None if witness is None else _untwisted(witness, source_profile)
        if cacheable:
            self._groupoid.record_isometry(source, target, result)
        return result

    def _lattice_action(self, vector: Lattices.ElementMethods) -> IntegralStructureAction:
        r"""The rational extensions acting at L, with a finite stable containing module.

        source: IndefiniteForms.g:988-1052 INDEF_FORM_GetVectorStructure for q(v) != 0:
        each generator of O(v^perp) is extended by the identity on v
        (``PmatInv * ExpandMatrix(eEndoRed) * Pmat``, :1019-1020).
        N=Zv+v-perp has finite index in L and is stable under every extension.
        If d annihilates L/N, then L is contained in the stable module N/d.
        This is the boundedness hypothesis for the finite residue computation.
        """
        lattice = vector.parent()
        perpendicular = vector.orthogonal_complement().inclusion()
        complement = perpendicular.domain()
        extensions = tuple(generator.extension_across(vector, vector) for generator in self.orthogonal_group_generators(complement))
        rational = lattice.base_change(lattice.base_ring().fraction_field_map())
        inclusion = ModulesOverIntegralDomains(lattice.base_ring())(lattice).generic_fibre_map()
        space = inclusion.codomain()
        split = lattice.subobject_on((vector, *(perpendicular(point) for point in complement.module_generators())))
        denominator = split.inclusion().cokernel().exponent()
        scale = rational.base_ring().one() / rational.base_ring()(denominator)
        containing = Modules(lattice.base_ring()).Mono(split, space)(
            lambda label: space.wrap(rational.scalar_multiple(scale, inclusion(split.inclusion()(split.module_generator(label))).underlying_element()))
        )
        return IntegralStructureAction(rational.O().subgroup(extensions), inclusion).with_stable_containing_bound(containing)

    def isotropic_stabilizer_action(
        self, isotropic: tuple[Lattices.ElementMethods, ...], radical_images: tuple[tuple[Lattices.ElementMethods, ...], ...]
    ) -> IntegralStructureAction:
        r"""The stabilizer problem of the primitive totally isotropic J = span(isotropic), at L.

        source: IndefiniteForms.g:1555-1568 INDEF_FORM_StabilizerVector with :988-1052
        INDEF_FORM_GetVectorStructure (isotropic v), and :1060-1131
        INDEF_FORM_GetRec_IsotropicKplane MapOrthogonalSublatticeGroup with :1710-1723
        INDEF_FORM_Stabilizer_IsotropicKstuff_Kernel (k-planes and flags); the stabilizer
        of L is MatrixIntegral_Stabilizer.
        Aut(N) for N = J-perp is Hom(K, J) x| (O(K) x Gamma) with K = J-perp/J and Gamma the
        allowed action on J (trivial for a pointed vector, +-1 for a line, GL(J) for a plane, the
        flag stabilizer for a flag), given as images of the vectors of J. Its generators are
        extended by the preamble Witt operation, the kernel generators at the common
        denominator of those extensions are adjoined, and the stabilizer of L in the rational group
        they generate is the stabilizer of J (with its marking). Every element preserves N, so the
        orbit span of L is computed inside L/D for D = 2 d^2 e^2 a, where d is that denominator,
        e = |disc K| and a = |coker(L -> J^*)| (issue #18); leaving that module raises.
        """
        lattice = isotropic[0].parent()
        reduction = lattice.primitive_isotropic_subobject(*isotropic).isotropic_reduction()
        quotient = reduction.quotient_lattice()
        perpendicular = reduction.orthogonal_complement().inclusion()
        frame = reduction.coordinate_frame()
        labels = tuple(quotient.module_generating_set())
        lifts = {label: perpendicular(frame(label)) for label in labels}
        basis = (*isotropic, *(lifts[label] for label in labels))

        def lift(element):
            coordinates = element.to_vector()
            return sum((lattice.scalar_multiple(coordinates(label), lifts[label]) for label in labels if coordinates(label)), lattice.zero())

        # Images of the basis (J, then the frame lifts), in that order.
        frame = tuple(lifts[label] for label in labels)
        image_families = []
        for isometry in self.orthogonal_group_generators(quotient):
            image_families.append((*isotropic, *(lift(isometry(quotient.module_generator(label))) for label in labels)))
        for radical in radical_images:
            image_families.append((*radical, *frame))
        for sheared in range(len(labels)):
            for point in isotropic:
                image_families.append((*isotropic, *(lifted + point if position == sheared else lifted for position, lifted in enumerate(frame))))
        extended = [_extend_perpendicular_isometry(reduction, reduction, basis, images) for images in image_families]
        denominator = SageZZ.one()
        for _isometry, entry_denominator in extended:
            denominator = denominator.lcm(entry_denominator)
        # Lorentzian.g:115-249 supplies 1 + T/d for every generator T of the
        # saturated integral homogeneous module. The preamble supplies 1 + T
        # through the complete pointwise-perpendicular kernel. Since T kills
        # N=J^perp and has image in J contained in N, T_i T_j = 0: the scaled
        # generators give exactly the source's (1/d)-parameter module.
        field_map = lattice.base_ring().fraction_field_map()
        rational = lattice.base_change(field_map)
        inverse_denominator = rational.base_ring().one() / rational.base_ring()(int(denominator))
        inclusion = ModulesOverIntegralDomains(lattice.base_ring())(lattice).generic_fibre_map()
        kernel_generators = []
        for integral_generator in reduction.pointwise_perpendicular_kernel().group_generators():
            # T = generator - 1 is read on the integral basis and rationalized;
            # 1 + T/d is the rational map with those scaled images.
            def scaled_image(label, generator=integral_generator):
                point = lattice.module_generator(label)
                direction = inclusion(generator(point) - point).underlying_element()
                return rational.module_generator(label) + rational.scalar_multiple(inverse_denominator, direction)

            kernel_generators.append(rational.Isom(rational)(scaled_image))
        generators = (*(isometry for isometry, _entry_denominator in extended), *kernel_generators)
        space = inclusion.codomain()
        # Unlike the rational coordinate complement, the pairing L -> J^*
        # determines the invariant index a in D=2d^2e^2a intrinsically.
        isotropic_inclusion = reduction.isotropic_embedding()
        # The pairing L -> I^* is i^* composed with the metric map,
        # where i: I -> L. Dualization is represented on the opposite
        # category, so its input is the reversed arrow op(L) -> op(I).
        dualization = Modules(lattice.base_ring()).dualization()
        opposite = dualization.domain()
        opposite_inclusion = opposite.Mor(opposite.object(lattice), opposite.object(isotropic_inclusion.domain()))(isotropic_inclusion)
        pairing = dualization(opposite_inclusion) * lattice.metric_map()
        # Here pairing is a map between finite free ZZ-modules L -> J*. Its
        # cokernel is finite, and its order is the product of the nonzero Smith
        # entries of its coordinate matrix. Compute that exact index from these
        # presentations: the generic provider's index() currently delegates to
        # cardinality() on a represented torsion module whose construction does
        # not yet supply a decided finite cardinality.
        pairing_domain_labels = tuple(pairing.domain().module_generating_set())
        pairing_codomain_labels = tuple(pairing.codomain().module_generating_set())
        pairing_rows = tuple(
            tuple(pairing(pairing.domain().module_generator(label)).to_vector()(target_label) for target_label in pairing_codomain_labels)
            for label in pairing_domain_labels
        )
        pairing_matrix = matrix(SageZZ, pairing_rows)
        pairing_rank = int(pairing_matrix.rank())
        if pairing_rank != len(pairing_codomain_labels):
            raise ArithmeticError("the isotropic pairing L -> J* must have finite cokernel")
        smith = pairing_matrix.smith_form()[0]
        pairing_index = SageZZ.one()
        for position in range(pairing_rank):
            pairing_index *= abs(smith[position, position])
        bound = 2 * int(denominator) ** 2 * abs(int(quotient.discriminant())) ** 2 * int(pairing_index)
        scale = rational.base_ring().one() / rational.base_ring()(bound)
        containing = Modules(lattice.base_ring()).Mono(lattice, space)(
            lambda label: space.wrap(rational.scalar_multiple(scale, inclusion(lattice.module_generator(label)).underlying_element()))
        )
        action = IntegralStructureAction(rational.O().subgroup(tuple(generators)), inclusion)
        span = OrbitSpanConstructions(lattice.base_ring())(action.ambient_action_functor(), inclusion, containing=containing)
        return action.with_stable_containing_bound(span.inclusion())

    def line_stabilizer_generators(self, vector: Lattices.ElementMethods) -> tuple[LatticeIsometryMethods, ...]:
        r"""Generators of Aut((L, Zv)) for primitive isotropic v: Gamma = {+-1} on the line."""
        lattice = vector.parent()
        action = self.isotropic_stabilizer_action((vector,), ((-vector,),))
        return tuple(lattice.integral_isometry(generator) for generator in action.stabilizer().group_generators())

    def pointwise_stabilizer_generators(self, vector: Lattices.ElementMethods) -> tuple[LatticeIsometryMethods, ...]:
        r"""Generators of `Aut((L, v))`, the **pointwise** stabilizer `{g : g(v) = v}` (issues #53 and #62).

        This is the group a quotient such as `L[n]/Gamma` runs over when `Gamma` is a pointwise stabilizer, and it
        is deliberately distinct from :meth:`line_stabilizer_generators`, which is `Aut((L, Zv))` and admits
        `-1` on the line. For a primitive isotropic `v` the two differ, and the difference is visible: the line
        group contains an element with `g(v) = -v`, the pointwise group contains none.

        The generators are the ones :meth:`pointed_automorphism_generators` computes, re-verified here: each is
        applied to `v` and must return `v` exactly, because the consumer's quotient is only the pointwise one if
        that holds for every generator.

        What this does *not* supply, and what decidability of equivalence does not give: a claim that these
        generators are complete for `Aut((L, v))` beyond the source's own statement, or any way to decide
        membership in the subgroup they generate. Deciding `Gamma`-equivalence of two vectors is a different
        problem from having generators, and this port decides it only as far as
        :meth:`pointwise_stabilizer_invariant` reaches.
        """
        generators = self.pointed_automorphism_generators(vector)
        for generator in generators:
            if generator(vector) != vector:
                raise ArithmeticError(
                    "a claimed pointwise stabilizer generator does not fix the marked vector; the quotient would "
                    "be over a larger group than the pointwise stabilizer"
                )
        return generators

    def pointwise_stabilizer_invariant(
        self, marked: tuple[Lattices.ElementMethods, ...], vector: Lattices.ElementMethods
    ) -> tuple[int, ...]:
        r"""An exact `Gamma`-invariant of `vector`, for `Gamma = {g : g(u) = u for u in marked}` (issue #53).

        The entry is `(q(v), content(v), b(u_1, v), ..., b(u_k, v))`. Every component is invariant: `q` and the
        content are preserved by any isometry, and `b(u_i, .)` is preserved because `g` fixes `u_i` and preserves
        the pairing, so `b(u_i, g v) = b(g u_i, g v) = b(u_i, v)`.

        Two vectors with different entries therefore lie in **different** `Gamma`-orbits -- a proved negative, not
        a search that failed. Equal entries prove nothing: the invariant is not claimed to be complete, and this
        port does not decide `Gamma`-equivalence beyond what the invariant separates.
        """
        lattice = vector.parent()
        for one in marked:
            if one.parent() is not lattice:
                raise ValueError(
                    "the pointwise stabilizer is of vectors of one lattice, and this vector belongs to another"
                )
        return (int(vector.q()), int(vector.content()), *(int(lattice.b(one, vector)) for one in marked))

    def pointwise_stabilizer_separates(
        self, marked: tuple[Lattices.ElementMethods, ...], source: Lattices.ElementMethods, target: Lattices.ElementMethods
    ) -> bool:
        r"""Whether :meth:`pointwise_stabilizer_invariant` **proves** `source` and `target` inequivalent.

        `True` is a decision: the two lie in different orbits of the pointwise stabilizer of `marked`. `False`
        means only that the invariant did not separate them, and must never be read as equivalence.
        """
        if source.parent() is not target.parent():
            raise ValueError("the invariant compares vectors of one lattice presentation")
        return self.pointwise_stabilizer_invariant(marked, source) != self.pointwise_stabilizer_invariant(
            marked, target
        )

    def verify_pointwise_stabilizer_transport(
        self,
        marked: tuple[Lattices.ElementMethods, ...],
        source: Lattices.ElementMethods,
        target: Lattices.ElementMethods,
        candidate: LatticeIsometryMethods,
    ) -> None:
        r"""Certify that `candidate` lies in the pointwise stabilizer of `marked` and carries `source` to `target`.

        A certificate *checker*, not a search: it runs in the time of a few evaluations and makes no completeness
        claim. Schematically, given the generator set, "does some word give a witness" is membership in a subgroup
        of infinite index, which this port does not decide -- so a caller supplies a candidate element (for
        example one built from the structure of the lattice) and this decides whether it really is a Gamma-element
        with the required action.
        """
        lattice = source.parent()
        if target.parent() is not lattice:
            raise ValueError("a pointwise stabilizer transport runs inside one lattice presentation")
        if candidate.domain() is not lattice or candidate.codomain() is not lattice:
            raise ValueError("the candidate must be an isometry of the lattice of the two markings")
        for one in marked:
            if candidate(one) != one:
                raise ArithmeticError(f"the candidate does not fix the marked vector {one} pointwise")
        if candidate(source) != target:
            raise ArithmeticError("the candidate does not carry the source marking onto the target marking")

    def levi_image(self, vector: Lattices.ElementMethods) -> Parent:
        r"""H, the image of Aut((L, Zv)) in O(K), K = v-perp/v, generated by the images of its computed generators."""
        reduction = vector.parent().primitive_isotropic_subobject(vector).isotropic_reduction()
        levi = reduction.levi_action()
        parabolic = reduction.parabolic_subgroup()
        return reduction.quotient_lattice().O().subgroup(tuple(levi(parabolic(generator)) for generator in self.line_stabilizer_generators(vector)))

    def _radical_pointed_isometry(self, source: Lattices.ElementMethods, target: Lattices.ElementMethods) -> LatticeIsometryMethods | None:
        r"""An isometry of `L` carrying the radical vector `source` onto `target`, or `None` (issues #41 and #62).

        **Conventions.** The port works in the row convention -- an isometry is a matrix `M` whose rows are the
        images of the basis vectors, and a vector `v` maps to `v M`. In the reduced presentation `L = C + R` with
        `R = rad(L)`, the retained block form (`_degenerate_group_generators`, `IndefiniteForms.g:853-901`) is

            M = [[A, B],
                 [0, D]]

        with the zero block **bottom-left**, `A` in `O(C)`, `D` in `GL(R)` and `B` in `Hom(C, R)` arbitrary. So

            (x, y) M = (x A, x B + y D),

        and a radical vector `(0, r)` maps to `(0, r D)`: it **stays in the radical**, because the radical is
        characteristic -- an isometry preserves the pairing, hence preserves `{v : b(v, L) = 0}`. `B` never moves a
        radical vector; it moves the nondegenerate part into the radical direction, which is exactly the shear
        `e_i -> e_i + z_1` of the source. (An earlier version of this docstring had the block layout transposed,
        wrote `(0, r) M = (r B, r D)`, and concluded that arbitrary `B` could push a radical vector out of the
        radical. Neither the layout nor the consequence is right.)

        **The decision.** Since a radical marking is `(0, r)` and its image is `(0, r D)`, transport is decided by
        `D` alone: `r D = s` with `D` unimodular, and `GL_k(Z)` acts on a column by its **content** and is
        transitive on primitive vectors, so `r ~ s` iff `gcd(r) = gcd(s)`. The zero vector is the content-0 case
        and is fixed by everything.

        **The witness built here** takes `A = 1` and `B = 0` -- a choice, not a requirement, and the simplest one
        that leaves `D` free: Smith normal form gives `V_r r_0 = e_1` and `V_s s_0 = e_1` for the primitive parts,
        `D = V_r V_s^{-1}` is unimodular, and `D r = s` is checked as it is built. The block is conjugated back to
        `L` with the retained `_conjugated` and the result is verified to carry the marking onto its image.

        Everything happens in the **reduced** presentation, so a permuted or mixed basis is handled by the
        retained reduction rather than by assuming the radical occupies the last coordinates of the given one.

        Cross-lattice radical transport is refused rather than guessed: it would need an isometry between the two
        reduced presentations carrying one radical coordinate onto the other, and this port does not decide that
        yet.
        """
        source_lattice, target_lattice = source.parent(), target.parent()
        if source_lattice is not target_lattice:
            return self._transported_radical_pointed_isometry(source, target)
        lattice = source_lattice
        reduction = indefinite_reduced_presentation(lattice)
        reduced = reduction.reduced
        transport = reduction.isometry
        gram = _gram_in_label_order(reduced)
        size, nonsingular_rank = int(gram.nrows()), int(gram.rank())
        radical_rank = size - nonsingular_rank
        if radical_rank == 0:
            raise ArithmeticError("a nondegenerate lattice has no radical vectors to transport")

        inverse = transport.inverse()
        labels = tuple(reduced.module_generating_set())
        def radical_coordinates(vector):
            coordinates = tuple(_integer(inverse(vector).to_vector()(label)) for label in labels)
            if any(coordinates[:nonsingular_rank]):
                raise ArithmeticError(
                    "a radical vector must vanish on the nondegenerate complement in the reduced presentation"
                )
            return coordinates[nonsingular_rank:]

        left = radical_coordinates(source)
        right = radical_coordinates(target)
        left_content = 0
        for entry in left:
            left_content = gcd(left_content, entry)
        right_content = 0
        for entry in right:
            right_content = gcd(right_content, entry)
        if left_content != right_content:
            return None

        if left_content == 0:
            block = identity_matrix(SageZZ, size)
        else:
            def normalize(row):
                """A unimodular `V` with `row V = e_1`, from Smith form of the primitive row."""
                primitive = [entry // left_content for entry in row]
                _, _, change = matrix(SageZZ, [primitive]).smith_form()
                if abs(int(change.det())) != 1:
                    raise ArithmeticError("Smith form must return a unimodular change of basis")
                if tuple(int(entry) for entry in (matrix(SageZZ, [primitive]) * change).list()) != (1,) + (0,) * (radical_rank - 1):
                    raise ArithmeticError("Smith reduction must carry a primitive row to e_1")
                return change
            from_left = normalize(left)
            from_right = normalize(right)
            radical_block = from_left * from_right.inverse().change_ring(SageZZ)
            if abs(int(radical_block.det())) != 1:
                raise ArithmeticError("the constructed radical block must be unimodular")
            if tuple(int(entry) for entry in (matrix(SageZZ, [left]) * radical_block).list()) != tuple(right):
                raise ArithmeticError("the constructed radical block must carry one content onto the other")
            block = identity_matrix(SageZZ, size)
            for row in range(radical_rank):
                for column in range(radical_rank):
                    block[nonsingular_rank + row, nonsingular_rank + column] = radical_block[row, column]
        witness = _conjugated(transport, (_isometry_from_row_matrix(reduced, reduced, block),))[0]
        if witness(source) != target:
            raise ArithmeticError("the constructed radical isometry does not carry the marking onto its image")
        return witness

    def radical_congruence_image(
        self, lattice: Lattices.ParentMethods, modulus: int
    ) -> RadicalCongruenceImage:
        r"""The image of `GL(R)` in `GL_k(Z/modulus)`, as a maintained matrix group (issue #41).

        The mixed marking condition of #41 needs `s - D r` in `c R`, i.e. `r D = s` modulo `c` for some `D` in
        `GL_k(Z)` -- the port multiplies rows on the right. The group that acts there is the **image** of `GL_k(Z)`,
        and this method builds it from integral generator *lifts*, so the object it returns cannot be mistaken for
        all of `GL_k(Z/c)`: a lift is integral and unimodular, hence has determinant `+-1`, and the image is
        `det^{-1}({+-1})` of index `phi(c)/2`, proper whenever `phi(c) > 2`.

        Everything here is the maintained machinery, not a parallel one: the generators are Sage integer matrices
        (`identity_matrix`, `diagonal_matrix`), their reductions are matrices over `Zmod`, `MatrixGroup` is Sage's
        GAP-backed matrix group, and `order`/`contains`/the general linear order all come from `MatrixGroup` and
        `GL`. There is no element enumeration, no hand-rolled group closure and no order formula, so no size bound
        is needed for what this object answers.

        What it does **not** answer, because no maintained contract supplies it yet, is which residue pairs
        transport: that needs a transporter in the image together with an integral word section so the lift can be
        assembled. That precise missing operation is routed in `.tmp/handoff-owner-congruence-transporter.md` rather
        than reimplemented here -- the port's mixed-marking decision itself stays on
        :meth:`pointed_isometry`, which decides every specimen.
        """
        modulus = int(modulus)
        if modulus < 2:
            raise ValueError(
                f"a congruence modulus of {modulus} makes the residue condition vacuous; no image is needed"
            )
        reduction = indefinite_reduced_presentation(lattice)
        gram = _gram_in_label_order(reduction.reduced)
        radical_rank = int(gram.nrows()) - int(gram.rank())
        if radical_rank == 0:
            raise ArithmeticError("a nondegenerate lattice has no radical, so there is no congruence image")

        if radical_rank == 1:
            lifts = (matrix(SageZZ, [[1]]), matrix(SageZZ, [[-1]]))
        else:
            elementary = []
            for row in range(radical_rank):
                for column in range(radical_rank):
                    if row == column:
                        continue
                    unit = identity_matrix(SageZZ, radical_rank)
                    unit[row, column] = 1
                    elementary.append(unit)
            sign = diagonal_matrix(SageZZ, [-1] + [1] * (radical_rank - 1))
            if tuple(sign.rows()[0]) != (-1,) + (0,) * (radical_rank - 1):
                raise ArithmeticError(
                    f"the sign generator must be diag(-1, 1, ...), got first row {tuple(sign.rows()[0])}"
                )
            lifts = tuple(elementary) + (sign,)

        residue_ring = Zmod(modulus)
        generators = tuple(
            matrix(
                residue_ring,
                [[int(lift[row, column]) % modulus for column in range(radical_rank)]
                 for row in range(radical_rank)],
            )
            for lift in lifts
        )
        for lift, generator in zip(lifts, generators, strict=True):
            if abs(int(lift.det())) != 1:
                raise ArithmeticError("a generator lift must be integral unimodular")
            if lift.change_ring(residue_ring) != generator:
                raise ArithmeticError("a lift does not reduce to its generator")
        group = MatrixGroup(list(generators))
        return RadicalCongruenceImage(
            modulus=modulus,
            radical_rank=radical_rank,
            group=group,
            generators=generators,
            lifts=lifts,
            general_linear_order=int(GL(radical_rank, residue_ring).order()),
        )

    def _mixed_marking_residue_obstruction(self, source: Lattices.ElementMethods, target: Lattices.ElementMethods) -> bool:
        r"""Whether the radical-residue condition admits a carrier; `False` is a **proved** obstruction (issue #41).

        In the reduced presentation `L = C + R` write `v = (x, r)` and `w = (y, s)`, with `c = content(x)`. Every
        isometry `g(v) = w` has a block form whose radical part `D` satisfies

            r D = s (mod c R),      i.e.  r D = s (mod c),

        so if no integral unimodular `D` does that, no isometry carries `v` onto `w` and the caller may answer
        `None` without any search. That is exactly the fast negative this returns `False` for.

        It is a **necessary condition only**: `True` reports "no obstruction found", never "an isometry exists" --
        the complementary condition `A x = y` with `A` in `O(C)` is not tested here, and `C`'s own orbit problem is
        the nondegenerate backend's. `c = 0` is the radical case, decided by `_radical_pointed_isometry`, and
        `c = 1` makes the congruence vacuous; both return `True`.

        The coordinates come from the lattice's cached reduced presentation, the same one
        `radical_congruence_image` reads, so the rows handed to the image are in one convention throughout: rows
        of length `radical_rank`, multiplied on the right.
        """
        lattice = source.parent()
        reduction = indefinite_reduced_presentation(lattice)
        reduced = reduction.reduced
        gram = _gram_in_label_order(reduced)
        nonsingular_rank = int(gram.rank())
        radical_rank = int(gram.nrows()) - nonsingular_rank
        if radical_rank == 0:
            return True
        labels = tuple(reduced.module_generating_set())
        inverse = reduction.isometry.inverse()

        def split(vector):
            coordinates = tuple(_integer(inverse(vector).to_vector()(label)) for label in labels)
            return coordinates[:nonsingular_rank], coordinates[nonsingular_rank:]

        source_complement, source_radical = split(source)
        _, target_radical = split(target)
        content = 0
        for entry in source_complement:
            content = gcd(content, entry)
        if content <= 1:
            return True
        image = self.radical_congruence_image(lattice, content)
        return image.transporter(source_radical, target_radical) is not None

    def _transported_radical_pointed_isometry(
        self, source: Lattices.ElementMethods, target: Lattices.ElementMethods
    ) -> LatticeIsometryMethods | None:
        r"""Radical pointed transport between two presentations (issues #41 and #62), or `None`.

        Two distinct parent objects are never identified: the result is an isometry from `source`'s parent to
        `target`'s parent, built by composing maintained operations rather than by assuming the two presentations
        coincide.

        The argument is short and it is also the reason the answer does not depend on which isometry is chosen.
        Suppose `phi : L -> M` is an isometry. Then every isometry `g : L -> M` with `g(v) = w` factors as
        `g = sigma . phi` with `sigma = g . phi^{-1}` in `O(M)`, so

            such a `g` exists  <=>  `phi(v)` and `w` lie in one `O(M)`-orbit,

        and if `phi_1`, `phi_2` are two isometries `L -> M` then `phi_2 phi_1^{-1}` is in `O(M)`, so `phi_1(v)` and
        `phi_2(v)` are in the same orbit: the class tested is independent of the choice. Both `phi(v)` and `w` lie
        in `rad(M)`, because an isometry preserves the radical, so the second question is the **same-lattice**
        radical decision already implemented above.

        Provable inequivalence is separated from the undecided case:

        * if no isometry `L -> M` exists, `self.isometry` returns `None` and this returns `None` -- the markings
          cannot be transported because the lattices themselves are not isometric;
        * if `isometry` cannot decide, it raises `UnfinishedCapability` and that propagates untouched;
        * if `phi(v)` and `w` differ in radical content, the same-lattice decision returns `None`.
        """
        source_lattice, target_lattice = source.parent(), target.parent()
        transport = self.isometry(source_lattice, target_lattice)
        if transport is None:
            return None
        image = transport(source)
        if not self.vector_lies_in_radical(image):
            raise ArithmeticError(
                "an isometry must carry the radical onto the radical, so a radical marking cannot have a "
                "non-radical image"
            )
        correction = self._radical_pointed_isometry(image, target)
        if correction is None:
            return None
        columns = _column_matrix(correction) * _column_matrix(transport)
        composed = _isometry_from_column_matrix(source_lattice, target_lattice, columns)
        if composed(source) != target:
            raise ArithmeticError("the composed radical isometry does not carry the marking onto its image")
        if composed.domain() is not source_lattice or composed.codomain() is not target_lattice:
            raise ArithmeticError("the composed witness must run between the two given parents, not one of them")
        return composed

    def _radical_stabilizer_of_row(self, row) -> tuple[Matrix_integer_dense, ...]:
        r"""Generators of `{D in GL_k(Z) : r D = r}` for a primitive row `r` (issue #41).

        The pointwise stabilizer of a vector inside the radical is `r D = r` on the radical block, so this
        is the `GL(R)` factor of it. For a primitive row, Smith reduction gives a unimodular `V` with
        `r V = e_1`, and in that basis the stabilizer of `e_1` is exactly `[[1, 0], [*, M]]` with
        `M` in `GL_{k-1}(Z)` and the column `*` arbitrary: generated by the shears `e_i -> e_i + e_1` for
        `i >= 2` together with `GL_{k-1}(Z)` on the remaining coordinates. Conjugating by `V` carries them
        back. `k = 1` has only `D = 1`, since `r D = r` forces `D = 1` in `GL_1(Z)`.
        """
        entries = tuple(int(entry) for entry in row)
        rank = len(entries)
        identity = identity_matrix(SageZZ, rank)
        content = 0
        for entry in entries:
            content = gcd(content, entry)
        if content == 0:
            raise ValueError("the zero vector is fixed by every isometry; it has no radical stabiliser")
        # R is torsion-free, so r D = r holds exactly when (r/content) D = r/content: the stabiliser of a
        # non-primitive row is the stabiliser of its primitive direction.
        entries = tuple(entry // content for entry in entries)
        if rank == 1:
            return (identity,)

        _, _, change = matrix(SageZZ, [list(entries)]).smith_form()
        assert change.det() in (1, -1)
        larger = OwnedGroups().GL(rank - 1, ZZ)

        def embedded(action) -> Matrix_integer_dense:
            block = identity_matrix(SageZZ, rank)
            for row_index in range(rank - 1):
                for column_index in range(rank - 1):
                    block[1 + row_index, 1 + column_index] = SageZZ(action[row_index, column_index])
            return block

        stabilizer_of_first = []
        for index in range(1, rank):
            shear = identity_matrix(SageZZ, rank)
            shear[index, 0] = 1
            stabilizer_of_first.append(shear)
        for generator in larger.group_generators():
            stabilizer_of_first.append(embedded(larger._to_engine(generator).matrix()))

        inverse = change.inverse().change_ring(SageZZ)
        return tuple(change * block * inverse for block in stabilizer_of_first)

    def _radical_blocks(self, lattice: Lattices.ParentMethods, vector: Lattices.ElementMethods) -> tuple[Matrix_integer_dense, ...]:
        r"""Row-matrix blocks of `Aut((L, v))` for `v = (0, r)` in a *reduced*, radical-last presentation.

        The retained reduction is what supplies the split: `indefinite_reduced_presentation` puts the
        radical in the final coordinates and provides the coordinate transport, so the block form
        `[[A,B],[0,D]]` is available here and nowhere else. `A` runs over `O(C)`, `D` over the `GL(R)`
        stabilizer of `r`, and `B` over all of `Hom(C,R)`, independently.
        """
        gram = _gram_in_label_order(lattice)
        size, complement_rank = gram.nrows(), gram.rank()
        radical_rank = size - complement_rank
        labels = tuple(lattice.module_generating_set())
        coordinates = tuple(_integer(vector.to_vector()(label)) for label in labels)
        if any(coordinates[:complement_rank]):
            raise ArithmeticError(
                "the retained reduction must place the radical in the final coordinates, but this vector "
                "has a nonzero nondegenerate-complement coordinate"
            )

        blocks = []
        if complement_rank:
            complement_gram = gram.matrix_from_rows_and_columns(range(complement_rank), range(complement_rank))
            complement = lattice.lattice_category()([[int(entry) for entry in row] for row in complement_gram.rows()])
            complement_labels = tuple(complement.module_generating_set())
            complement_basis = tuple(complement.module_generator(label) for label in complement_labels)
            for generator in self.orthogonal_group_generators(complement):
                top = tuple(
                    tuple(_integer(generator(point).to_vector()(label)) for label in complement_labels)
                    for point in complement_basis
                )
                block = identity_matrix(SageZZ, size)
                for row_index in range(complement_rank):
                    for column_index in range(complement_rank):
                        block[row_index, column_index] = top[row_index][column_index]
                blocks.append(block)

        if radical_rank:
            for action in self._radical_stabilizer_of_row(coordinates[complement_rank:]):
                block = identity_matrix(SageZZ, size)
                for row_index in range(radical_rank):
                    for column_index in range(radical_rank):
                        block[complement_rank + row_index, complement_rank + column_index] = action[row_index, column_index]
                blocks.append(block)

        for row_index in range(complement_rank):
            for radical_index in range(radical_rank):
                block = identity_matrix(SageZZ, size)
                block[row_index, complement_rank + radical_index] = 1
                blocks.append(block)
        return tuple(blocks)

    def _degenerate_radical_pointed_generators(self, vector: Lattices.ElementMethods) -> tuple[LatticeIsometryMethods, ...]:
        r"""Generators of `Aut((L, v))` for `v` inside the radical of a degenerate `L` (issue #41).

        The radical is *not* assumed to occupy the final coordinates of the caller's basis. The retained
        reduction `indefinite_reduced_presentation` supplies the split and its exact coordinate transport,
        so a permuted or mixed basis is carried into the radical-last presentation, the block generators are
        built there, and the results are conjugated back to the caller's lattice. `v = 0` fixes every
        isometry, so the answer is simply the whole group.
        """
        lattice = vector.parent()
        labels = tuple(lattice.module_generating_set())
        if not any(_integer(vector.to_vector()(label)) for label in labels):
            return self.orthogonal_group_generators(lattice)

        reduction = indefinite_reduced_presentation(lattice)
        reduced = reduction.reduced
        transport = reduction.isometry
        reduced_vector = transport.inverse()(vector)
        blocks = self._radical_blocks(reduced, reduced_vector)
        reduced_generators = tuple(_isometry_from_row_matrix(reduced, reduced, block) for block in blocks)
        return _conjugated(transport, reduced_generators)

    def vector_lies_in_radical(self, vector: Lattices.ElementMethods) -> bool:
        r"""Whether the vector pairs trivially with the whole lattice, i.e. lies in `rad(L)`.

        `q(v) = 0` is necessary but not sufficient -- an isotropic vector off the radical is isotropic too
        -- so the pairing against the basis is what decides it.
        """
        lattice = vector.parent()
        if lattice.is_nondegenerate():
            return False
        gram = _gram_in_label_order(lattice)
        labels = tuple(lattice.module_generating_set())
        coordinates = [_integer(vector.to_vector()(label)) for label in labels]
        return not any(
            sum(gram[row_index, column_index] * coordinates[column_index] for column_index in range(len(coordinates)))
            for row_index in range(len(coordinates))
        )

    def pointed_automorphism_generators(self, vector: Lattices.ElementMethods) -> tuple[LatticeIsometryMethods, ...]:
        r"""Generators of Aut((L, v)): the stabilizer of L in the rational extensions of Aut(v-perp) (IndefiniteForms.g:1533-1545).

        A vector inside the radical of a degenerate lattice takes the split-radical construction of issue
        #41: the isotropic reduction refuses a rank-deficient ambient form, and the radical case has its own
        exact description.
        """
        lattice = vector.parent()
        if self.vector_lies_in_radical(vector):
            return self._degenerate_radical_pointed_generators(vector)
        action = self.isotropic_stabilizer_action((vector,), ()) if vector.q() == lattice.base_ring().zero() else self._lattice_action(vector)
        return tuple(lattice.integral_isometry(generator) for generator in action.stabilizer().group_generators())

    def pointed_isometry(self, source: Lattices.ElementMethods, target: Lattices.ElementMethods) -> LatticeIsometryMethods | None:
        r"""An element of Isom((L, v), (M, w)), or None."""
        if VectorPrefilter.from_vector(source) != VectorPrefilter.from_vector(target):
            return None
        source_lattice, target_lattice = source.parent(), target.parent()
        if source_lattice is target_lattice and source == target:
            return source_lattice.O().identity()
        if source_lattice is not target_lattice and LatticePrefilter.from_lattice(source_lattice) != LatticePrefilter.from_lattice(target_lattice):
            return None
        if not source_lattice.is_nondegenerate() and self.vector_lies_in_radical(source):
            # A radical marking is not an isotropic *line*: the source's isotropic machinery needs the line
            # primitive and non-degenerate, and it asserts for a radical vector of content > 1. The radical has
            # its own exact description (see `_radical_pointed_isometry`), so route it before that machinery.
            return self._radical_pointed_isometry(source, target)
        if (
            source_lattice is target_lattice
            and not source_lattice.is_nondegenerate()
            and not self._mixed_marking_residue_obstruction(source, target)
        ):
            # A proved obstruction: the radical residues cannot be matched modulo `content(x)`, so no isometry
            # carries one marking onto the other. Answering here is exact and skips the search below.
            return None
        if source.q() == source_lattice.base_ring().zero():
            if source_lattice.is_even() and source.div() == source_lattice.base_ring().one():
                return self._split_isotropic_isometry(source, target)
            return self._isotropic_pointed_isometry(source, target)
        source_complement = source.orthogonal_complement().inclusion().domain()
        target_complement = target.orthogonal_complement().inclusion().domain()
        reduced = self.isometry(source_complement, target_complement)
        if reduced is None:
            return None
        rational_transport = reduced.extension_across(source, target)
        if rational_transport.is_integral_on(source_lattice, target_lattice):
            return target_lattice.integral_isometry(rational_transport, source=source_lattice)
        fraction_map = source_lattice.base_ring().fraction_field_map()
        restriction = Modules(fraction_map.codomain()).restriction_of_scalars(fraction_map)
        image = restriction(rational_transport) * ModulesOverIntegralDomains(source_lattice.base_ring())(source_lattice).generic_fibre_map()
        correction = self._lattice_action(target).transporter(image)
        return None if correction is None else target_lattice.integral_isometry(correction.inverse() * rational_transport, source=source_lattice)

    def _reduced_isotropic_extension(self, source_reduction, target_reduction, source, target, reduced) -> LatticeIsometryMethods:
        r"""The rational extension L_Q -> M_Q of the integral isometry N -> N' that ``reduced`` determines.

        source: IndefiniteForms.g:1623-1641 (vectors) and :1673-1686 (k-planes and
        flags). The source's ``test`` is an integral isometry of the perpendicular
        lattices N=J^perp -> N'=J'^perp; here it sends the basis ``source`` of J to
        ``target`` and the frame lift of each k in K=N/J to the frame lift of
        ``reduced(k)``. Both families are Z-bases of N and N', so this is an
        integral isometry onto N' and carries the marking as given (the sign
        normalization of :1639-1641 is then the identity). Its extension is
        ``LORENTZ_ExtendOrthogonalIsotropicIsomorphism``.
        """
        target_lattice = target[0].parent()
        source_quotient = source_reduction.quotient_lattice()
        target_quotient = target_reduction.quotient_lattice()
        source_labels = tuple(source_quotient.module_generating_set())
        target_labels = tuple(target_quotient.module_generating_set())
        source_frame = source_reduction.coordinate_frame()
        target_frame = target_reduction.coordinate_frame()
        source_perpendicular = source_reduction.orthogonal_complement().inclusion()
        target_perpendicular = target_reduction.orthogonal_complement().inclusion()
        source_lifts = tuple(source_perpendicular(source_frame(label)) for label in source_labels)
        target_lifts = {label: target_perpendicular(target_frame(label)) for label in target_labels}

        def target_lift(element):
            coordinates = element.to_vector()
            return sum((target_lattice.scalar_multiple(coordinates(label), target_lifts[label]) for label in target_labels if coordinates(label)), target_lattice.zero())

        images = (*target, *(target_lift(reduced(source_quotient.module_generator(label))) for label in source_labels))
        rational_transport, _denominator = _extend_perpendicular_isometry(source_reduction, target_reduction, (*source, *source_lifts), images)
        return rational_transport

    def _isotropic_correction(
        self,
        rational_transport: LatticeIsometryMethods,
        source_lattice,
        target: tuple[Lattices.ElementMethods, ...],
        radical_images: tuple[tuple[Lattices.ElementMethods, ...], ...],
    ) -> LatticeIsometryMethods | None:
        r"""source: IndefiniteForms.g:1569-1604 Kernel_Equivalence_Qmat.

        An integral ``EquivRat`` is returned as it is (:1578-1580). Otherwise
        the source searches its extended stabilizer group of J (:1581-1590,
        MatrixIntegral_Equivalence_Bis) for g with ``EquivRat * g`` integral.
        Here the same search is the integral-structure transporter of the
        conjugate group at J' (:meth:`isotropic_stabilizer_action`), with the
        allowed action ``radical_images`` on J'. The transport must carry N onto
        N' integrally (:meth:`_reduced_isotropic_extension`); every group element
        preserves N', so a correction exists exactly when an isometry with the
        marking exists.
        """
        target_lattice = target[0].parent()
        if rational_transport.is_integral_on(source_lattice, target_lattice):
            return target_lattice.integral_isometry(rational_transport, source=source_lattice)
        fraction_map = source_lattice.base_ring().fraction_field_map()
        restriction = Modules(fraction_map.codomain()).restriction_of_scalars(fraction_map)
        image = restriction(rational_transport) * ModulesOverIntegralDomains(source_lattice.base_ring())(source_lattice).generic_fibre_map()
        correction = self.isotropic_stabilizer_action(tuple(target), radical_images).transporter(image)
        return None if correction is None else target_lattice.integral_isometry(correction.inverse() * rational_transport, source=source_lattice)

    def isotropic_subspace_isometry(
        self, source: tuple[Lattices.ElementMethods, ...], target: tuple[Lattices.ElementMethods, ...], radical_images: tuple[tuple[Lattices.ElementMethods, ...], ...]
    ) -> LatticeIsometryMethods | None:
        r"""An isometry L -> M carrying J = span(source) onto J' = span(target), compatible with the marking.

        source: IndefiniteForms.g:1660-1694 INDEF_FORM_Equivalence_IsotropicKstuff_Kernel.
        An isometry psi of the reductions K -> K' (``f_equiv``, :1673) gives an
        isometry of the perpendicular lattices, its rational extension
        (:meth:`_reduced_isotropic_extension`), and Kernel_Equivalence_Qmat
        (:meth:`_isotropic_correction`), whose group's allowed action
        ``radical_images`` on J' (images of ``target``) absorbs the choice of
        basis of J'.
        """
        source_lattice = source[0].parent()
        target_lattice = target[0].parent()
        source_reduction = source_lattice.primitive_isotropic_subobject(*source).isotropic_reduction()
        target_reduction = target_lattice.primitive_isotropic_subobject(*target).isotropic_reduction()
        reduced = self.isometry(source_reduction.quotient_lattice(), target_reduction.quotient_lattice())
        if reduced is None:
            return None
        rational_transport = self._reduced_isotropic_extension(source_reduction, target_reduction, source, target, reduced)
        return self._isotropic_correction(rational_transport, source_lattice, target, radical_images)

    def _isotropic_pointed_isometry(self, source: Lattices.ElementMethods, target: Lattices.ElementMethods) -> LatticeIsometryMethods | None:
        r"""Isom of primitive isotropic marked vectors.

        source: IndefiniteForms.g:1605-1651 INDEF_FORM_EquivalenceVector, isotropic
        branch. ``test`` (:1623) is a recursive isometry of the reductions; the
        rational transport extends the integral isometry of perpendicular
        lattices it determines and is integral already in the Ansatz case
        (:1578). Otherwise the preamble's marked-line lift torsor supplies the
        complete integral fibre over this reduction isometry (port-design
        "Isometry-extension fibres"). An empty fibre for this one reduction map
        does not make the pointed Isom empty: Kernel_Equivalence_Qmat then
        ranges over the full stabilizer of the marked vector, which acts on K
        through all of O(K).
        """
        source_lattice, target_lattice = source.parent(), target.parent()
        source_reduction = source_lattice.primitive_isotropic_subobject(source).isotropic_reduction()
        target_reduction = target_lattice.primitive_isotropic_subobject(target).isotropic_reduction()
        reduced = self.isometry(source_reduction.quotient_lattice(), target_reduction.quotient_lattice())
        if reduced is None:
            return None
        rational_transport = self._reduced_isotropic_extension(source_reduction, target_reduction, (source,), (target,), reduced)
        if rational_transport.is_integral_on(source_lattice, target_lattice):
            return target_lattice.integral_isometry(rational_transport, source=source_lattice)

        torsor = source_reduction.rational_lifts(target_reduction, reduced)
        if torsor.integral_members().is_empty() is False:
            integral = torsor.integral_base_point()
            # The torsor fixes the reductions' chosen line generators; when one
            # is the negative of the supplied marking, the global sign of
            # IndefiniteForms.g:1639-1641 is applied against the actual vectors.
            fraction_map = source_lattice.base_ring().fraction_field_map()
            source_vector = ModulesOverIntegralDomains(source_lattice.base_ring())(source_lattice).generic_fibre_map()(source).underlying_element()
            target_vector = ModulesOverIntegralDomains(target_lattice.base_ring())(target_lattice).generic_fibre_map()(target).underlying_element()
            mapped = integral(source_vector)
            if mapped == -target_vector:
                target_rational = target_lattice.base_change(fraction_map)
                integral = target_rational.Isom(target_rational)(lambda label: -target_rational.module_generator(label)) * integral
            elif mapped != target_vector:
                raise ArithmeticError("the integral reduction lift does not carry the marked isotropic line onto the requested target")
            return target_lattice.integral_isometry(integral, source=source_lattice)
        return self._isotropic_correction(rational_transport, source_lattice, (target,), ())

    def _split_isotropic_isometry(self, source: Lattices.ElementMethods, target: Lattices.ElementMethods) -> LatticeIsometryMethods | None:
        r"""For e of divisibility one in an even lattice, L = U + e^perp/e; an isometry is id_U plus one of the reductions."""
        source_splitting = source.parent().primitive_isotropic_subobject(source).isotropic_reduction().integral_hyperbolic_splitting()
        target_splitting = target.parent().primitive_isotropic_subobject(target).isotropic_reduction().integral_hyperbolic_splitting()
        hyperbolic_plane, source_complement = source_splitting.codomain().biproduct_factors()
        _hyperbolic_plane, target_complement = target_splitting.codomain().biproduct_factors()
        reduced = self.isometry(source_complement, target_complement)
        if reduced is None:
            return None
        split_map = hyperbolic_plane.identity_morphism().biproduct_map(reduced, source=source_splitting.codomain(), target=target_splitting.codomain())
        # The biproduct operation returns a module map.  Its summands are
        # isometries, so retain that structure before composing the splittings.
        split = source_splitting.codomain().Isom(target_splitting.codomain())(
            lambda label: split_map(source_splitting.codomain().module_generator(label))
        )
        return target_splitting.inverse() * split * source_splitting

    def primitive_vector_classes(self, lattice: Lattices.ParentMethods, square: int) -> tuple[VectorClass, ...]:
        r"""One class per O(L)-orbit of primitive vectors of the square, for L of positive index at least 2."""
        profile = AttackProfile.from_lattice(lattice)
        if profile.positive_index < 2:
            raise ValueError(
                f"{lattice} has positive index {profile.positive_index}: vector classes here need positive index at least 2; the Lorentzian case is marked cells (#15)"
            )
        # Real positive and negative index at least two does not imply
        # rational isotropy: diag(2, 2, -6, -6), for example, is anisotropic
        # by descent modulo 3.  Its zero-square fibre has no primitive
        # vectors, so decide that empty case before requesting the two
        # rational hyperbolic pairs of the Eichler envelope.
        if square == 0 and not lattice.is_isotropic():
            return ()
        buckets: dict[VectorPrefilter, list[tuple[Lattices.ElementMethods, list[LatticeIsometryMethods]]]] = {}
        for candidate in self._model(profile.signed_view).primitive_covering_vectors(profile.sign * int(square)):
            bucket = buckets.setdefault(VectorPrefilter.from_vector(candidate), [])
            match next(((witnesses, isometry) for kept, witnesses in bucket if (isometry := self.pointed_isometry(kept, candidate)) is not None), None):
                case None:
                    bucket.append((candidate, []))
                case (witnesses, isometry):
                    witnesses.append(isometry)
        return tuple(
            VectorClass(
                lattice.linear_combination({label: kept.to_vector()(label) for label in lattice.module_generating_set()}),
                tuple(_untwisted(isometry, profile) for isometry in witnesses),
            )
            for bucket in buckets.values()
            for kept, witnesses in bucket
        )


__all__ = ["EichlerRecursion", "MarkedOrbitCanonical", "MarkedOrbitDecision", "VectorClass"]
