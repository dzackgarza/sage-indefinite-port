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
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from math import gcd, isqrt

import numpy as np

from sage_extensions.all import GObjects, Lattices, Modules, Sets, ZZ
from sage_extensions.categories.group.groups import OwnedGroups
from sage_extensions.categories.lattice_morphisms import LatticeIsometryMethods
from sage_extensions.categories.modules.orbit_spans import (
    OrbitSpanConstructions,
)
from sage_extensions.categories.modules.pure.modules import (
    Modules,
    ModulesOverIntegralDomains,
    _engine_matrix,
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
from sage.modules.free_module_element import vector
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

type _IsotropicActionKey = tuple[
    int, tuple[tuple[int, ...], ...], tuple[tuple[tuple[int, ...], ...], ...],
    tuple[tuple[tuple[int, ...], ...], ...] | None, int,
]


def _untwisted(isometry: LatticeIsometryMethods, profile: AttackProfile) -> LatticeIsometryMethods:
    r"""An isometry of the signed view, as an isometry of the lattice itself."""
    return isometry if profile.sign == 1 else _UNTWIST(isometry)


def _integral_stabilizer_family(action: IntegralStructureAction) -> Iterable[LatticeIsometryMethods]:
    r"""Realize the provider's full Schreier family from its supplied factors.

    GroupAction.g:869-886 supplies the integral stabilizer before the source
    generator simplification. The provider retains its orbit, letters and
    transversals and specifies ``t_(s.p)^-1 s t_p`` in orbit-major order.
    The provider's ambient-group inclusion retains each factor's actual
    rational isometry. Composites use owned isometry multiplication and
    preserve every source family entry.
    """
    letters, transversals = action.stabilizer_schreier_factors()
    g_set = action.finite_representation().g_set()
    if g_set.cardinality() == 1:
        yield from letters
        return
    points = tuple(g_set)
    permutations = g_set.orbit_generator_permutations(points[0])
    positions = {point: position for position, point in enumerate(points)}
    inclusion = action.ambient_group_inclusion()
    realized_letters = tuple(inclusion(letter) for letter in letters)
    realized_transversals = tuple(inclusion(transversal) for transversal in transversals)
    inverse_transversals = tuple(~transversal for transversal in realized_transversals)
    for position, point in enumerate(points):
        for letter, permutation in zip(realized_letters, permutations, strict=True):
            image_position = positions[permutation(point)]
            yield inverse_transversals[image_position] * letter * realized_transversals[position]


def _rational_action_columns(lattice, generators, *, inclusion=None):
    r"""Read the supplied maps in the lattice's retained rational framing."""
    labels = tuple(lattice.module_generating_set())
    if inclusion is None:
        inclusion = ModulesOverIntegralDomains(lattice.base_ring())(lattice).generic_fibre_map()
    if inclusion.domain() is not lattice:
        raise ValueError("the rational framing must retain the supplied integral lattice")
    basis = tuple(inclusion(lattice.module_generator(label)).underlying_element() for label in labels)
    field = inclusion.codomain().extension_ring()
    frame = field.free_module(lattice.module_generating_set())
    comparisons = {}
    inverse_comparisons = {}
    for generator in generators:
        domain, codomain = generator.domain(), generator.codomain()
        for endpoint in (domain, codomain):
            if endpoint not in comparisons:
                comparison = Modules(field).Mor(frame, endpoint)({
                    label: endpoint(point) for label, point in zip(labels, basis, strict=True)
                })
                comparisons[endpoint] = _engine_matrix(comparison)
        if codomain not in inverse_comparisons:
            inverse_comparisons[codomain] = comparisons[codomain].inverse()
        yield (
            inverse_comparisons[codomain]
            * _engine_matrix(generator.module_morphism())
            * comparisons[domain]
        )


def _integral_stabilizer_action_columns(lattice, action):
    r"""Every full integral Schreier action, with factors in one rational frame.

    GroupAction.g:869-886 constructs the integral stabilizer before
    MatrixGroupSimplification.h reduces its full matrix generators.
    The provider supplies t_(s.p)^-1 s t_p and the exact orbit action.
    Rational factors are composed before their product is admitted over ZZ.
    """
    g_set = action.finite_representation().g_set()
    if g_set.cardinality() == 1:
        for columns in _rational_action_columns(
            lattice, _integral_stabilizer_family(action), inclusion=action.lattice_inclusion()
        ):
            yield columns.change_ring(SageZZ)
        return
    letters, transversals = action.stabilizer_schreier_factors()
    inclusion = action.ambient_group_inclusion()
    factors = itertools.chain(
        (inclusion(letter) for letter in letters),
        (inclusion(transversal) for transversal in transversals),
    )
    columns = tuple(_rational_action_columns(lattice, factors, inclusion=action.lattice_inclusion()))
    letter_columns = columns[:len(letters)]
    transversal_columns = columns[len(letters):]
    inverse_transversals = tuple(transversal.inverse() for transversal in transversal_columns)
    points = tuple(g_set)
    permutations = g_set.orbit_generator_permutations(points[0])
    positions = {point: position for position, point in enumerate(points)}
    for position, point in enumerate(points):
        for letter, permutation in zip(letter_columns, permutations, strict=True):
            image_position = positions[permutation(point)]
            yield (
                inverse_transversals[image_position] * letter * transversal_columns[position]
            ).change_ring(SageZZ)


def _integral_generator_set(
    lattice: Lattices.ParentMethods,
    generators: Iterable[LatticeIsometryMethods] | IntegralStructureAction,
) -> tuple[LatticeIsometryMethods, ...]:
    r"""Reduce full integral actions before constructing integral maps.

    Pinned MatrixGroupSimplification.h:650-751 reduces pairs by products
    and inverses; :1386-1487 applies this to the entrywise L1 norm.
    These are invertible Nielsen transformations, so the integral kernel
    survives even when several generators have the same finite image.
    Native matrices are temporary exact basis actions; the result consists
    of owned lattice isometries. A supplied integral-structure action
    realizes its complete Schreier family in the common frame first.
    """
    labels = tuple(lattice.module_generating_set())
    columns = (
        _integral_stabilizer_action_columns(lattice, generators)
        if isinstance(generators, IntegralStructureAction)
        else _rational_action_columns(lattice, generators)
    )
    identity = tuple(int(row == column) for row in range(len(labels)) for column in range(len(labels)))
    seen = {identity}
    actions = {}
    norm_counts = Counter()

    def retire(key):
        norm = actions.pop(key)[0]
        norm_counts[norm] -= 1
        if norm_counts[norm] == 0:
            del norm_counts[norm]

    # MatrixGroupSimplification.h:1495-1605 selects bounded integer
    # arithmetic. Entrywise L1 is submultiplicative. If both directions
    # of each pair have norm <= B, every product and partial dot sum has
    # absolute sum <= B^2 <= 2^63-1. Larger pairs use Python integers.
    machine_bound = isqrt(int(np.iinfo(np.int64).max))

    def action_data(action):
        # Fixed-size row-major keys preserve equality and lexicographic order.
        key = tuple(action.ravel().tolist())
        return key, int(np.add.reduce(np.abs(action), axis=None))

    def normalized(action, inverse, direct_norm, inverse_norm):
        dtype = np.dtype(np.int64) if max(direct_norm, inverse_norm) <= machine_bound else np.dtype("object")
        action = np.asarray(action, dtype=dtype)
        inverse = np.asarray(inverse, dtype=dtype)
        # Choose a fixed representative of {g,g^-1} at equal complexity.
        if inverse_norm < direct_norm:
            return tuple(inverse.ravel().tolist()), (inverse_norm, inverse, action)
        direct_key = tuple(action.ravel().tolist())
        if inverse_norm == direct_norm:
            inverse_key = tuple(inverse.ravel().tolist())
            if inverse_key < direct_key:
                return inverse_key, (inverse_norm, inverse, action)
        return direct_key, (direct_norm, action, inverse)

    for action_columns in columns:
        engine_action = action_columns.transpose().change_ring(SageZZ)
        action = np.asarray(tuple(map(int, engine_action.list())), dtype="object").reshape(len(labels), len(labels))
        key, direct_norm = action_data(action)
        if key in seen:
            continue
        seen.add(key)
        engine_inverse = engine_action.inverse().change_ring(SageZZ)
        inverse = np.asarray(tuple(map(int, engine_inverse.list())), dtype="object").reshape(len(labels), len(labels))
        # Nielsen updates preserve the subgroup: both g and g^-1 remain
        # words in the retained family after either action is retired.
        seen.add(tuple(inverse.ravel().tolist()))
        key, entry = normalized(action, inverse, direct_norm, int(np.add.reduce(np.abs(inverse), axis=None)))
        if key in actions:
            continue
        actions[key] = entry
        norm_counts[entry[0]] += 1
        pending = [key]
        # A newly inserted action must meet every retained action. An
        # unchanged pair needs no second visit. Replacements strictly lower
        # total L1 complexity, including when a repeated action disappears.
        # Thus this worklist reaches the source's exhaustive pair reduction
        # without constructing all pending pairs of a Schreier family.
        while pending:
            current = pending.pop()
            if current not in actions:
                continue
            for other in tuple(actions):
                if current not in actions:
                    break
                if other == current:
                    continue
                a_norm, a, a_inverse = actions[current]
                b_norm, b, b_inverse = actions[other]
                # The same four source products and their inverses, evaluated
                # together. The L1 bound above certifies every integer product
                # and reduction independently of this batching.
                products = np.asarray((a, a_inverse, a, a_inverse, b_inverse, b_inverse, b, b)) @ np.asarray(
                    (b, b, b_inverse, b_inverse, a_inverse, a, a_inverse, a)
                )
                norms = np.add.reduce(np.abs(products), axis=(1, 2))
                # Equal integral actions have equal L1 norms. A product
                # can change this source step only by shortening the pair
                # or matching a retained action, so other full matrix keys
                # need no materialization.
                candidates = tuple(
                    normalized(products[index], products[index + 4], int(norms[index]), int(norms[index + 4]))
                    for index in range(4)
                    if (
                        min(int(norms[index]), int(norms[index + 4])) < max(a_norm, b_norm)
                        or min(int(norms[index]), int(norms[index + 4])) in norm_counts
                    )
                )
                if not candidates:
                    continue
                # The same source pair products can certify redundancy
                # even when no L1 norm decreases. If c = a^eps*b^eta
                # (up to inversion) is already retained independently of
                # a and b, either a or b is a word in the other and c.
                # Delete one using equality of the full integral actions;
                # equality only in a finite quotient would lose kernels.
                if any(key in actions and key not in (current, other) for key, _entry in candidates):
                    retire(current if a_norm >= b_norm else other)
                    continue
                candidate_key, candidate = min(candidates, key=lambda item: (item[1][0], item[0]))
                if candidate[0] >= max(a_norm, b_norm):
                    continue
                replaced = current if a_norm >= b_norm else other
                retire(replaced)
                if candidate_key != identity and candidate_key not in actions:
                    actions[candidate_key] = candidate
                    norm_counts[candidate[0]] += 1
                    pending.append(candidate_key)

    result = []
    for key in actions:
        # The stabilizer supplies an automorphism preserving this integral
        # structure. Its full basis action has just been read exactly; use
        # those images in the owned isometry homset instead of evaluating
        # the rational expression again through integral restriction.
        integral_images = {}
        for position, source_label in enumerate(labels):
            row = key[position * len(labels):(position + 1) * len(labels)]
            integral_images[source_label] = lattice.linear_combination(
                {
                    label: lattice.base_ring()(coefficient)
                    for label, coefficient in zip(labels, row, strict=True)
                }
            )
        result.append(lattice.Isom(lattice)(integral_images))
    return tuple(result)


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
    extension = _extend_rational_perpendicular(
        rational_source, rational_target,
        tuple(source_rationalize(point).underlying_element() for point in source_basis),
        tuple(target_rationalize(point).underlying_element() for point in images),
    )
    (denominator,) = tuple(extension.denominator_ideal(source, target).ideal_generators())
    return extension, int(denominator)


def _extend_rational_perpendicular(
    rational_source: Lattices.ParentMethods, rational_target: Lattices.ParentMethods,
    source_basis: tuple[Lattices.ElementMethods, ...], images: tuple[Lattices.ElementMethods, ...],
) -> LatticeIsometryMethods:
    r"""Invoke the pinned Witt extension in the supplied rational framings.

    CombinedAlgorithms.h:314-332 LiftToFullAutomorphism changes the
    ambient framing before invoking the same perpendicular extension.
    """
    field = rational_source.base_ring()
    rank = len(source_basis)
    # The source's Subspace1 and Subspace2 are the rows of the two families.
    # Here they are the maps D -> L_Q and D -> M_Q out of the rational
    # lattice D on the basis, with the Gram pulled back along the first; the
    # partial isometry is the identity of D. No coordinate solve is needed.
    frame = field.free_module(rank)
    labels = tuple(frame.module_generating_set())
    gram = rational_source.gram_tensor().pullback(
        Modules(field).Mor(frame, rational_source)(dict(zip(labels, source_basis, strict=True)))
    )
    span = Lattices(field)([[gram[row, column] for column in range(rank)] for row in range(rank)])
    span_labels = tuple(span.module_generating_set())
    source_inclusion = Modules(field).Mor(span, rational_source)(
        dict(zip(span_labels, source_basis, strict=True))
    )
    target_inclusion = Modules(field).Mor(span, rational_target)(dict(zip(span_labels, images, strict=True)))
    partial = span.Isom(span).identity()
    return partial.witt_extension(source_inclusion, target_inclusion)


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



def _fair_unimodular_matrices(n: int):
    r"""Enumerate GL_n(Z) exactly by increasing elementary-generator word length.

    Elementary row additions, basis swaps and a sign change generate
    GL_n(Z).  Breadth-first enumeration is therefore exhaustive.  This
    procedure must be called only after an independent *exact positive*
    isometry decision; otherwise it need not terminate.
    """
    from collections import deque

    identity = identity_matrix(SageZZ, n)
    if n == 0:
        # GL_0(ZZ) is the trivial group. In particular there is no first
        # basis coordinate on which to apply the sign generator below.
        yield identity
        return
    generators = []
    for i in range(n):
        for j in range(n):
            if i != j:
                for sign in (-1, 1):
                    elementary = identity_matrix(SageZZ, n)
                    elementary[i, j] = sign
                    generators.append(elementary)
    for i in range(n - 1):
        swap = identity_matrix(SageZZ, n)
        swap.swap_rows(i, i + 1)
        generators.append(swap)
    sign = identity_matrix(SageZZ, n)
    sign[0, 0] = -1
    generators.append(sign)
    frontier = deque((identity,))
    seen = {tuple(identity.list())}
    while frontier:
        candidate = frontier.popleft()
        yield candidate
        for generator in generators:
            next_candidate = candidate * generator
            key = tuple(next_candidate.list())
            if key not in seen:
                seen.add(key)
                frontier.append(next_candidate)



def _fair_unimodular_isometry_witness(source: Lattices.ParentMethods, target: Lattices.ParentMethods):
    r"""First exact integral witness, after an independently certified positive decision."""
    source_gram = _gram_in_label_order(source)
    target_gram = _gram_in_label_order(target)
    for candidate in _fair_unimodular_matrices(int(source_gram.nrows())):
        if candidate * target_gram * candidate.transpose() == source_gram:
            return _isometry_from_row_matrix(source, target, candidate)


def iter_exact_integral_orthogonal_isometries(lattice: Lattices.ParentMethods):
    r"""Fair, nontruncated enumeration of all integral isometries of a finite-rank lattice.

    The output is a stream of owned O(L)-elements, not a finite generating
    family.  Explicit arithmetic parametrizations and breadth-first
    enumeration of complete group generators cover the specialized cases;
    otherwise the generic GL_n(ZZ) matrix stream is exhaustive.  Infinite
    orthogonal groups produce infinite streams, while the finite definite
    and split binary branches terminate after their last element.
    """
    gram = _gram_in_label_order(lattice)
    if int(gram.nrows()) == 0:
        yield lattice.O().identity()
        return
    if int(gram.nrows()) == 1:
        # GL_1(ZZ)={+1,-1}; both preserve every integral rank-one
        # quadratic form, including the degenerate zero form.
        yield _isometry_from_row_matrix(lattice, lattice, identity_matrix(SageZZ, 1))
        yield _isometry_from_row_matrix(lattice, lattice, -identity_matrix(SageZZ, 1))
        return
    if not gram:
        # For the zero form, O(L)=GL(L).  Use the maintained integral
        # general-linear generating family, with inverses, instead of
        # constructing the exponentially branching GL_n word search.
        # This is a genuine group enumeration: all words and their
        # inverses occur at finite breadth-first distance.
        from collections import deque

        linear = OwnedGroups().GL(int(lattice.rank()), lattice.base_ring())
        forward = tuple(linear._to_engine(g).matrix().change_ring(SageZZ)
                        for g in linear.group_generators())
        generators = forward + tuple(g.inverse().change_ring(SageZZ) for g in forward)
        identity = identity_matrix(SageZZ, int(lattice.rank()))
        frontier = deque((identity,))
        seen = {tuple(identity.list())}
        while frontier:
            current = frontier.popleft()
            yield _isometry_from_row_matrix(lattice, lattice, current)
            for generator in generators:
                successor = current * generator
                key = tuple(successor.list())
                if key not in seen:
                    seen.add(key)
                    frontier.append(successor)
        return
    if not lattice.is_nondegenerate() and int(lattice.rank()) >= 3:
        # Reduce the radical to the final coordinates and extend the
        # nondegenerate complement group by all integral radical shears
        # and GL(R) actions.  This is the exact degenerate extension
        # algorithm, not an ambient GL_n enumeration.  Its generators
        # are complete whenever the complement group generators are.
        from collections import deque

        recursion = EichlerRecursion()
        forward = tuple(_column_matrix(g) for g in recursion.orthogonal_group_generators(lattice))
        generators = forward + tuple(g.inverse().change_ring(SageZZ) for g in forward)
        identity = identity_matrix(SageZZ, int(lattice.rank()))
        frontier = deque((identity,))
        seen = {tuple(identity.list())}
        while frontier:
            current = frontier.popleft()
            yield _isometry_from_column_matrix(lattice, lattice, current)
            for generator in generators:
                successor = current * generator
                key = tuple(successor.list())
                if key not in seen:
                    seen.add(key)
                    frontier.append(successor)
        return
    if int(gram.nrows()) == 2 and lattice.is_nondegenerate() and not lattice.is_definite() and lattice.is_isotropic():
        # Every rational isometry of a split binary plane permutes its two
        # rational null lines, scaling them by mutually inverse factors.
        # For an integral unimodular isometry each primitive integral ray
        # maps to a signed primitive ray; hence there are at most eight
        # possibilities.  The pinned binary algorithm tests all eight.
        from sage_indefinite_port.indefinite.lorentzian_cells import _binary_split_automorphisms

        yield from _binary_split_automorphisms(lattice)
        return
    signature = tuple(map(int, lattice.signature_pair()))
    if lattice.is_nondegenerate() and not lattice.is_definite() and (
        int(lattice.rank()) == 2
        or (int(lattice.witt_index()) > 0 and (min(signature) == 1 or int(lattice.witt_index()) >= 2))
    ):
        # The pinned indefinite generation step supplies the full integral
        # orthogonal group in every supported signature: Lorentzian leaves
        # use the perfect-domain algorithm; higher Witt index uses the
        # approximate-model covering plus exact pointed stabilizer.  In
        # either case enumerate the resulting Cayley graph, not GL_n(ZZ).
        # The signed-view convention is undone on every generator below.
        from collections import deque

        # A group generating family is not generally a *monoid*
        # generating family.  In particular an infinite-order Pell period
        # and its negative powers must both occur in the enumeration.
        # Include inverse matrices explicitly before traversing words.
        if int(lattice.rank()) == 2:
            # The rank-two arithmetic edgewalk is the exact provider for
            # both split and anisotropic indefinite binary forms; the
            # higher-rank Eichler generation step requires a rational
            # isotropic vector and is not applicable to anisotropic planes.
            forward = tuple(_column_matrix(g) for g in lorentzian_orthogonal_group_generators(lattice))
        else:
            forward = tuple(_column_matrix(g)
                            for g in EichlerRecursion().orthogonal_group_generators(lattice))
        generators = forward + tuple(g.inverse().change_ring(SageZZ) for g in forward)
        identity = identity_matrix(SageZZ, int(lattice.rank()))
        frontier = deque((identity,))
        seen = {tuple(identity.list())}
        while frontier:
            current = frontier.popleft()
            yield _isometry_from_column_matrix(lattice, lattice, current)
            for generator in generators:
                successor = current * generator
                key = tuple(successor.list())
                if key not in seen:
                    seen.add(key)
                    frontier.append(successor)
        return
    if int(gram.nrows()) == 2 and gram.rank() == 1:
        # A rank-one radical in a binary lattice is primitive.  Choose a
        # unimodular basis with the radical last; every isometry in that
        # basis is [[epsilon, b], [0, delta]] (row convention), where
        # epsilon,delta=+/-1 and b ranges over ZZ.  Unlike a generic
        # GL_2 word search this gives every isometry once, by an explicit
        # arithmetic parametrization, and reaches every shear fairly.
        reduction = indefinite_reduced_presentation(lattice)
        reduced = reduction.reduced
        reduced_gram = _gram_in_label_order(reduced)
        if reduced_gram[0, 1] or reduced_gram[1, 0] or reduced_gram[1, 1]:
            raise ArithmeticError("binary radical reduction did not put the radical last")
        from itertools import count

        transport = _column_matrix(reduction.isometry)
        inverse_transport = transport.inverse().change_ring(SageZZ)
        for height in count(0):
            shears = (0,) if height == 0 else (-height, height)
            for shear in shears:
                for first_sign in (-1, 1):
                    for radical_sign in (-1, 1):
                        rows = matrix(SageZZ, [[first_sign, shear], [0, radical_sign]])
                        reduced_map = _isometry_from_row_matrix(reduced, reduced, rows)
                        yield _isometry_from_column_matrix(
                            lattice, lattice,
                            transport * _column_matrix(reduced_map) * inverse_transport,
                        )
        return
    if lattice.is_definite():
        # Closure of the exact definite O(L) generators is finite.  Matrix
        # multiplication takes place in the same labelled integral basis
        # as the owned isometry constructor; no heuristic cutoff is used.
        from collections import deque

        generators = tuple(_column_matrix(g) for g in lattice.O().group_generators())
        identity = identity_matrix(SageZZ, int(gram.nrows()))
        frontier = deque((identity,))
        seen = {tuple(identity.list())}
        while frontier:
            current = frontier.popleft()
            yield _isometry_from_column_matrix(lattice, lattice, current)
            for generator in generators:
                successor = current * generator
                key = tuple(successor.list())
                if key not in seen:
                    seen.add(key)
                    frontier.append(successor)
        return
    for candidate in _fair_unimodular_matrices(int(gram.nrows())):
        if candidate * gram * candidate.transpose() == gram:
            yield _isometry_from_row_matrix(lattice, lattice, candidate)


def _norm_fibre_has_local_obstruction(gram, square: int, *, primitive: bool) -> bool:
    r"""Certify norm-fibre emptiness by integral-value and local congruences.

    For primitive vectors the residue must remain nonzero modulo the
    underlying prime.  Exhaustive tests are bounded to 10000 tuples.
    """
    rank = int(gram.nrows())
    value_ideal = 0
    for i in range(rank):
        value_ideal = gcd(value_ideal, _integer(gram[i, i]))
        for j in range(i + 1, rank):
            value_ideal = gcd(value_ideal, 2 * _integer(gram[i, j]))
    if value_ideal and int(square) % value_ideal:
        return True
    for modulus, prime in ((2, 2), (3, 3), (4, 2), (8, 2), (9, 3), (5, 5), (7, 7), (25, 5)):
        if modulus ** rank > 10000:
            continue
        if not any(
            sum(int(gram[i, j]) * x[i] * x[j] for i in range(rank) for j in range(rank)) % modulus == int(square) % modulus
            for x in itertools.product(range(modulus), repeat=rank)
            if not primitive or any(value % prime for value in x)
        ):
            return True
    return False


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
    Owned composition is the product of the framed matrices, so each
    conjugate retains the matrix P H P^-1.
    """
    inverse = ~transport
    return tuple(transport * generator * inverse for generator in generators)


def _transported(target_transport: LatticeIsometryMethods, witness: LatticeIsometryMethods, source_transport: LatticeIsometryMethods) -> LatticeIsometryMethods:
    r"""``target_transport * witness * source_transport^-1``, built once from its column matrix.

    source: IndefiniteForms.g:1389 ``Inverse(RecRed2.B) * test * RecRed1.B`` and
    :1419 (the degenerate extension ``Inverse(FullBasis2) * TheEquivTest * FullBasis1``).
    """
    return target_transport * witness * ~source_transport


def _witt_zero_generators(lattice):
    r"""Return ``O(L)`` generators for ``L`` isometric to the norm lattice of ``ZZ<1,i,j,k>`` of ``(-1, 3)``, or ``None``."""
    from sage.matrix.constructor import matrix as sage_matrix

    from sage_indefinite_port.indefinite.li_margulis import is_isometric_in_box
    from sage_indefinite_port.indefinite.quaternion_lattices import QuaternionNormLattice
    from sage_indefinite_port.indefinite.quaternion_voronoi import QuaternionVoronoi, orthogonal_group_generators

    model = QuaternionNormLattice(-1, 3, 2)
    source = model.lattice()
    transport = is_isometric_in_box(source, lattice, 3)
    if transport is None:
        return None
    frame = tuple(source.module_generators())
    columns = sage_matrix([[int(c) for c in transport(v)] for v in frame]).transpose()
    inverse = columns.inverse()
    labels = tuple(lattice.module_generating_set())
    target = tuple(lattice.module_generators())
    result = []
    for generator in orthogonal_group_generators(model, QuaternionVoronoi(-1, 3), model.normalizer_classes_in_box(3)):
        on_model = sage_matrix([[int(c) for c in generator(v)] for v in frame]).transpose()
        moved = columns * on_model * inverse
        result.append(
            lattice.O()(
                {
                    label: sum(
                        (lattice.scalar_multiple(int(c), g) for c, g in zip(moved.column(j), target, strict=True)),
                        lattice.zero(),
                    )
                    for j, label in enumerate(labels)
                }
            )
        )
    return tuple(result)


def _witt_one_generators(lattice):
    r"""Return isometries generating ``O(L)`` from the Hilbert modular model, or ``None``.

    For every real quadratic maximal order: Vaserstein's generation of ``SL_2(O)`` by elementary matrices plus
    the exact cusp links (:func:`recognized_generators`).
    """
    from sage_indefinite_port.indefinite.witt_one_hilbert import recognized_generators

    matrices = recognized_generators(lattice)
    if matrices is None:
        return None
    frame = tuple(lattice.module_generators())
    labels = tuple(lattice.module_generating_set())
    return tuple(
        lattice.O()(
            {
                label: sum(
                    (lattice.scalar_multiple(int(c), g) for c, g in zip(m.column(j), frame, strict=True)),
                    lattice.zero(),
                )
                for j, label in enumerate(labels)
            }
        )
        for m in matrices
    )


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
        self._isotropic_actions: dict[_IsotropicActionKey, IntegralStructureAction] = {}

    def congruence_image(self, lattice: Lattices.ParentMethods, modulus: int):
        r"""Actual finite image of ``O(lattice)`` on ``L / modulus L`` (#42).

        Uses the complete constructive integral orthogonal generators and the
        lattice's coordinate-label order. Composite moduli are supported.
        """
        from sage_indefinite_port.indefinite.congruence_images import IntegralCongruenceImage

        generators = self.orthogonal_group_generators(lattice)
        return IntegralCongruenceImage(
            tuple(_column_matrix(g) for g in generators),
            modulus,
            rank=len(tuple(lattice.module_generating_set())),
        )

    def congruence_stabilizer_generators(self, lattice: Lattices.ParentMethods, modulus: int, residue):
        r"""Integral isometries fixing a coordinate residue modulo ``modulus``."""
        image = self.congruence_image(lattice, modulus)
        index, matrices = image.residue_vector_stabilizer(residue)
        return index, tuple(_isometry_from_column_matrix(lattice, lattice, g) for g in matrices)

    def congruence_submodule_stabilizer_generators(self, lattice: Lattices.ParentMethods, modulus: int, generators):
        r"""Integral setwise stabilizer of a submodule of ``L / modulus L``."""
        image = self.congruence_image(lattice, modulus)
        index, matrices = image.residue_submodule_stabilizer(generators)
        return index, tuple(_isometry_from_column_matrix(lattice, lattice, g) for g in matrices)

    def congruence_flag_stabilizer_generators(self, lattice: Lattices.ParentMethods, modulus: int, stages):
        r"""Integral setwise stabilizer of a nested flag in ``L / modulus L``."""
        image = self.congruence_image(lattice, modulus)
        index, matrices = image.residue_flag_stabilizer(stages)
        return index, tuple(_isometry_from_column_matrix(lattice, lattice, g) for g in matrices)

    def joint_determinant_congruence_image(self, lattice: Lattices.ParentMethods, modulus: int):
        r"""Actual joint determinant and level image of the integral orthogonal group."""
        from sage_indefinite_port.indefinite.congruence_images import JointDeterminantCongruenceImage

        generators = self.orthogonal_group_generators(lattice)
        return JointDeterminantCongruenceImage(
            tuple(_column_matrix(g) for g in generators),
            modulus,
            rank=len(tuple(lattice.module_generating_set())),
        )

    def joint_determinant_congruence_preimage_generators(self, lattice: Lattices.ParentMethods, modulus: int, allowed):
        r"""Integral generators for a subgroup of the actual joint level/sign image."""
        image = self.joint_determinant_congruence_image(lattice, modulus)
        index, matrices = image.preimage_generators(allowed)
        return index, tuple(_isometry_from_column_matrix(lattice, lattice, g) for g in matrices)

    def finite_index_sublattice_stabilizer_generators(self, lattice: Lattices.ParentMethods, modulus: int, basis):
        r"""Integral isometries preserving the specified finite-index sublattice.

        ``basis`` consists of coordinate rows, with ``modulus L`` contained
        in their span. The exact congruence image supplies the stabilizer.
        """
        image = self.congruence_image(lattice, modulus)
        index, matrices = image.finite_index_sublattice_stabilizer(basis)
        return index, tuple(_isometry_from_column_matrix(lattice, lattice, g) for g in matrices)

    def joint_real_spinor_congruence_image(self, lattice: Lattices.ParentMethods, modulus: int):
        r"""Actual joint image of level reduction and normalized real spinor sign.

        Uses the provider's ``O_plus`` character with the Gritsenko--Hulek--
        Sankaran ``-b/2`` reflection normalization, not determinant or
        positive-cone component preservation.
        """
        from sage_indefinite_port.indefinite.congruence_images import JointSignCongruenceImage

        character = lattice.O_plus().preimage_morphism()

        def sign(columns):
            isometry = _isometry_from_column_matrix(lattice, lattice, columns)
            return 1 if character(isometry).is_one() else -1

        generators = self.orthogonal_group_generators(lattice)
        return JointSignCongruenceImage(
            tuple(_column_matrix(g) for g in generators),
            modulus, sign, rank=len(tuple(lattice.module_generating_set())),
        )

    def joint_real_spinor_congruence_preimage_generators(self, lattice: Lattices.ParentMethods, modulus: int, allowed):
        r"""Integral generators of a joint level/real-spinor inverse image."""
        image = self.joint_real_spinor_congruence_image(lattice, modulus)
        index, matrices = image.preimage_generators(allowed)
        return index, tuple(_isometry_from_column_matrix(lattice, lattice, g) for g in matrices)

    def joint_determinant_real_spinor_congruence_image(self, lattice: Lattices.ParentMethods, modulus: int):
        r"""Actual image of level, determinant and normalized real-spinor signs."""
        from sage_indefinite_port.indefinite.congruence_images import JointTwoSignsCongruenceImage

        character = lattice.O_plus().preimage_morphism()

        def spinor(columns):
            isometry = _isometry_from_column_matrix(lattice, lattice, columns)
            return 1 if character(isometry).is_one() else -1

        generators = self.orthogonal_group_generators(lattice)
        return JointTwoSignsCongruenceImage(
            tuple(_column_matrix(g) for g in generators), modulus,
            lambda g: int(g.det()), spinor,
            rank=len(tuple(lattice.module_generating_set())),
        )

    def joint_determinant_real_spinor_congruence_preimage_generators(self, lattice: Lattices.ParentMethods, modulus: int, allowed):
        r"""Integral preimage generators for a subgroup of the actual triple image."""
        image = self.joint_determinant_real_spinor_congruence_image(lattice, modulus)
        index, matrices = image.preimage_generators(allowed)
        return index, tuple(_isometry_from_column_matrix(lattice, lattice, g) for g in matrices)

    def congruence_colouring_stabilizer_generators(self, lattice: Lattices.ParentMethods, modulus: int, colour_classes):
        r"""Integral stabilizer preserving each named residue-vector colour."""
        image = self.congruence_image(lattice, modulus)
        index, matrices = image.residue_colouring_stabilizer(colour_classes)
        return index, tuple(_isometry_from_column_matrix(lattice, lattice, g) for g in matrices)

    def joint_field_spinor_congruence_image(self, lattice: Lattices.ParentMethods, modulus: int, field_map=None, form_multiplier=None):
        r"""Actual finite generated image of level and field-spinor square classes.

        The field-spinor normalization is exactly that of the provider's
        ``spinor_norm``; unlike the real-spinor sign this retains full
        square classes, including the finite image in an infinite target.
        """
        from sage_indefinite_port.indefinite.congruence_images import JointFieldSpinorCongruenceImage

        character = lattice.spinor_norm(field_map, form_multiplier)

        def spinor(columns):
            isometry = _isometry_from_column_matrix(lattice, lattice, columns)
            return character(isometry)

        generators = self.orthogonal_group_generators(lattice)
        return JointFieldSpinorCongruenceImage(
            tuple(_column_matrix(g) for g in generators), modulus,
            spinor, character.codomain().one(),
            rank=len(tuple(lattice.module_generating_set())),
        )

    def joint_field_spinor_congruence_preimage_generators(self, lattice: Lattices.ParentMethods, modulus: int, allowed, field_map=None, form_multiplier=None):
        r"""Integral generators of a subgroup of the actual joint field-spinor image."""
        image = self.joint_field_spinor_congruence_image(lattice, modulus, field_map, form_multiplier)
        index, matrices = image.preimage_generators(allowed)
        return index, tuple(_isometry_from_column_matrix(lattice, lattice, g) for g in matrices)

    def field_spinor_kernel_generators(self, lattice: Lattices.ParentMethods, field_map=None, form_multiplier=None):
        r"""Exact integral generators of the field-spinor kernel.

        The finite generated spinor image is computed using the joint level
        model at modulus two and then *all* reduction classes with trivial
        spinor are retained. Thus no level-two restriction is imposed.
        """
        image = self.joint_field_spinor_congruence_image(lattice, 2, field_map, form_multiplier)
        one = image.identity[1]
        allowed = (key for key in image.image_keys() if key[1] == one)
        index, matrices = image.preimage_generators(allowed)
        return index, tuple(_isometry_from_column_matrix(lattice, lattice, g) for g in matrices)

    def principal_congruence_subgroup_generators(self, lattice: Lattices.ParentMethods, modulus: int):
        r"""Generators of ``ker(O(L) -> Aut(L/modulus*L))`` (principal congruence subgroup)."""
        image = self.congruence_image(lattice, modulus)
        index, matrices = image.preimage_generators((image.image_one,))
        return index, tuple(_isometry_from_column_matrix(lattice, lattice, g) for g in matrices)

    def joint_positive_cone_congruence_image(self, lattice: Lattices.ParentMethods, modulus: int):
        r"""Actual joint image of congruence reduction and positive-cone component."""
        from sage_indefinite_port.indefinite.congruence_images import JointSignCongruenceImage

        # Calling the provider predicate enforces its exact signature (1,n)
        # hypothesis and selected positive-cone convention.
        lattice.O_component()  # Enforce the provider's exact signature hypothesis.

        def sign(columns):
            isometry = _isometry_from_column_matrix(lattice, lattice, columns)
            return 1 if isometry.preserves_positive_cone() else -1

        generators = self.orthogonal_group_generators(lattice)
        return JointSignCongruenceImage(
            tuple(_column_matrix(g) for g in generators), modulus, sign,
            rank=len(tuple(lattice.module_generating_set())),
        )

    def joint_positive_cone_congruence_preimage_generators(self, lattice: Lattices.ParentMethods, modulus: int, allowed):
        r"""Integral generators of a subgroup of the actual level/component image."""
        image = self.joint_positive_cone_congruence_image(lattice, modulus)
        index, matrices = image.preimage_generators(allowed)
        return index, tuple(_isometry_from_column_matrix(lattice, lattice, g) for g in matrices)

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
        if int(lattice.rank()) == 0:
            # O(0) is trivial. In particular, PARI's definite qfauto
            # normalization divides by the gcd of the empty Gram matrix,
            # so it is not an algorithm for this boundary object.
            return ()
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
            # PARI qfauto is sensitive to the Gram presentation: a badly
            # conditioned perpendicular complement can overflow its stack
            # even though the lattice's finite orthogonal group is small.
            # Reduce the *integral* definite Gram before invoking the
            # provider's exact qfauto kernel, then conjugate its generators
            # through the certified unimodular LLL isometry.
            reduction = lattice.lll_reduction()
            reduced = reduction.reduced
            reduced_generators = self._groupoid.cached_group(reduced) if cacheable else None
            if reduced_generators is None:
                reduced_generators = tuple(reduced.O().group_generators())
                if cacheable:
                    self._groupoid.record_group(reduced, reduced_generators)
            definite_generators = _conjugated(reduction.isometry, reduced_generators)
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
                case _ if int(profile.signed_view.module_rank()) == 4 and int(profile.signed_view.witt_index()) == 0 and (
                    witt_zero := _witt_zero_generators(profile.signed_view)
                ) is not None:
                    generators = witt_zero
                case _ if int(profile.signed_view.module_rank()) == 4 and int(profile.signed_view.witt_index()) == 1 and (
                    witt_one := _witt_one_generators(profile.signed_view)
                ) is not None:
                    generators = witt_one
                case _ if int(profile.signed_view.witt_index()) < 2:
                    # The approximate model needs two hyperbolic planes; the source
                    # stops below rational Witt index two (IndefiniteForms.g:830).
                    raise UnfinishedCapability(
                        f"O(L) generators for {profile.signed_view}: positive index {profile.positive_index}, "
                        f"rational Witt index {profile.signed_view.witt_index()} and no recognized model; "
                        "work unit #69 owns this case"
                    )
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
        r"""An isometry between reduced indefinite lattices, or None.

        source: IndefiniteForms.g:1320-1359 and ApproximateModels.h:1305-1390.

        Documented departure (#68): over ZZ at rank at least 3 the verdict is
        decided before the covering search by genus and spinor-genus equality.
        Eichler-Kneser (Kneser 1956; Cassels, Rational Quadratic Forms, Ch. 11,
        Thm. 1.3 and Cor.; Conway-Sloane, SPLAG Ch. 15 Sec. 9): an indefinite
        genus of rank at least 3 has one class per spinor genus. The source
        reaches a negative verdict only by exhausting the covering family; for
        a same-genus pair of large 2-adic index (diag(1,-2,-64,128) and the
        other class of its genus) its coset expansion exceeds 10^6 cosets and
        does not finish. The positive witness keeps the source covering route,
        which must then succeed.
        """
        cacheable = source.base_ring() == ZZ and target.base_ring() == ZZ
        spinor_decided = cacheable and source.module_rank() >= 3
        if cacheable:
            known, cached = self._groupoid.cached_isometry(source, target)
            if known:
                return cached
        source_profile = AttackProfile.from_lattice(source)
        target_profile = AttackProfile.from_lattice(target)
        signed_source, signed_target = source_profile.signed_view, target_profile.signed_view
        if source_profile.positive_index >= 2 and (int(source.witt_index()) < 2 or int(target.witt_index()) < 2):
            decision = source.Isom(target).is_empty()
            if decision is True:
                result = None
            elif decision is False:
                result = _fair_unimodular_isometry_witness(source, target)
            else:
                raise ArithmeticError("low-Witt integral isometry decision is not exact")
            if cacheable:
                self._groupoid.record_isometry(source, target, result)
            return result
        if spinor_decided:
            # Departure from IndefiniteForms.g:1320-1359 and
            # ApproximateModels.h:1305-1390 (#68): the source decides a
            # same-genus pair only by exhausting the covering family, whose
            # coset expansion exceeds 10^6 cosets at large 2-adic index.
            # Eichler-Kneser (Kneser 1956; Cassels, Rational Quadratic Forms,
            # Ch. 11, Thm. 1.3 and Cor.; SPLAG Ch. 15 Sec. 9): an indefinite
            # genus of rank >= 3 has one class per spinor genus. So distinct
            # genera or spinor genera give None, and otherwise the covering
            # search below must find a witness.
            from sage_extensions.categories.lattice_engines import isometric_within_genus

            if source.genus() != target.genus() or not isometric_within_genus(source, target):
                self._groupoid.record_isometry(source, target, None)
                return None
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
        if spinor_decided:
            assert witness is not None, "Eichler-Kneser: one spinor genus is one class, so the covering search must find a witness"
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
        self, isotropic: tuple[Lattices.ElementMethods, ...], radical_images: tuple[tuple[Lattices.ElementMethods, ...], ...],
        *, quotient_generators: tuple[LatticeIsometryMethods, ...] | None = None, extension_denominator: int = 1,
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

        Supplied quotient generators replace O(K); the empty family fixes K
        pointwise. Lorentzian.g:225-249 permits the common denominator d in
        its saturated kernel module. Include a prescribed rational lift's
        denominator when this action corrects that lift in its quotient fibre.
        """
        lattice = isotropic[0].parent()
        labels = tuple(lattice.module_generating_set())
        if extension_denominator < 1:
            raise ValueError("the isotropic extension denominator must be positive")
        reduction = lattice.primitive_isotropic_subobject(*isotropic).isotropic_reduction()
        quotient = reduction.quotient_lattice()
        quotient_labels = tuple(quotient.module_generating_set())
        if quotient_generators is not None and any(
            generator.domain() is not quotient or generator.codomain() is not quotient
            for generator in quotient_generators
        ):
            raise ValueError("the supplied quotient generators must act on this isotropic reduction")
        quotient_key = None if quotient_generators is None else tuple(
            tuple(tuple(_integer(generator(point).to_vector()(label)) for label in quotient_labels) for point in quotient.module_generators())
            for generator in quotient_generators
        )

        def coordinates(point: Lattices.ElementMethods) -> tuple[int, ...]:
            if point.parent() is not lattice:
                raise ValueError("an isotropic stabilizer action needs all basis images in its ambient lattice")
            row = point.to_vector()
            return tuple(_integer(row(label)) for label in labels)

        key = (
            id(lattice), tuple(coordinates(point) for point in isotropic),
            tuple(tuple(coordinates(point) for point in images) for images in radical_images),
            quotient_key, extension_denominator,
        )
        cached = self._isotropic_actions.pop(key, None)
        if cached is not None:
            self._isotropic_actions[key] = cached
            return cached
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
        reduced_generators = self.orthogonal_group_generators(quotient) if quotient_generators is None else quotient_generators
        for isometry in reduced_generators:
            image_families.append((*isotropic, *(lift(isometry(quotient.module_generator(label))) for label in labels)))
        for radical in radical_images:
            image_families.append((*radical, *frame))
        for sheared in range(len(labels)):
            for point in isotropic:
                image_families.append((*isotropic, *(lifted + point if position == sheared else lifted for position, lifted in enumerate(frame))))
        extended = [_extend_perpendicular_isometry(reduction, reduction, basis, images) for images in image_families]
        denominator = SageZZ(extension_denominator)
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
        # The pairing L -> I^* has matrix (b(e_i, i(f_j))) in integral
        # generator bases. Construct this matrix directly: the opposite
        # module category need not accept a lattice embedding as one of its
        # own arrows, whereas these pairings are intrinsic to the lattice.
        isotropic_basis = tuple(isotropic_inclusion(point)
                               for point in isotropic_inclusion.domain().module_generators())
        ambient_basis = tuple(lattice.module_generators())
        pairing_matrix = matrix(SageZZ, [
            [_integer(lattice.b(point, isotropic)) for isotropic in isotropic_basis]
            for point in ambient_basis
        ])
        # Here pairing is a map between finite free ZZ-modules L -> J*. Its
        # cokernel is finite, and its order is the product of the nonzero Smith
        # entries of its coordinate matrix. Compute that exact index from these
        # presentations: the generic provider's index() currently delegates to
        # cardinality() on a represented torsion module whose construction does
        # not yet supply a decided finite cardinality.
        pairing_rank = int(pairing_matrix.rank())
        if pairing_rank != len(isotropic_basis):
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
        # The source's scalar dilation bounds the orbit span but need not
        # itself be invariant under nonintegral rational generators. Close
        # the orbit inside that finite bound before supplying a stable module.
        span = OrbitSpanConstructions(lattice.base_ring())(action.ambient_action_functor(), inclusion, containing=containing)
        result = action.with_stable_containing_bound(span.inclusion())
        # Invariant comparison, extension and integral correction consume
        # the same source f_stab datum. Retain the two current actions, with
        # their exact bases and allowed radical images, rather than rebuild
        # their Witt lifts and stable modules. Each action retains L through
        # its inclusion, so a live key cannot be reused for another parent.
        self._isotropic_actions[key] = result
        while len(self._isotropic_actions) > 2:
            del self._isotropic_actions[next(iter(self._isotropic_actions))]
        return result

    def line_stabilizer_generators(self, vector: Lattices.ElementMethods) -> tuple[LatticeIsometryMethods, ...]:
        r"""Generators of Aut((L, Zv)) for nonzero isotropic v, via its primitive saturation."""
        lattice = vector.parent()
        if vector == lattice.zero() or vector.q() != lattice.base_ring().zero():
            raise ValueError("a line stabilizer requires a nonzero isotropic vector")
        if int(vector.content()) > 1:
            # Stab(Z(c*v)) = Stab(Zv) for every nonzero integer c.
            return self.line_stabilizer_generators(vector.primitive_part())
        if self.vector_lies_in_radical(vector):
            # The radical is characteristic. Its exact pointed stabilizer
            # is already the split C+R block group; central -id_L adds
            # the missing sign and gives the entire moving-line group.
            minus = _isometry_from_row_matrix(
                lattice, lattice, -identity_matrix(SageZZ, int(lattice.rank()))
            )
            return (*self.pointed_automorphism_generators(vector), minus)
        if not lattice.is_nondegenerate():
            reduction = indefinite_reduced_presentation(lattice)
            reduced = reduction.reduced
            gram = _gram_in_label_order(reduced)
            n, size = int(gram.rank()), int(gram.nrows())
            labels = tuple(reduced.module_generating_set())
            coordinates = (~reduction.isometry)(vector).to_vector()
            x = tuple(_integer(coordinates(label)) for label in labels[:n])
            r = tuple(_integer(coordinates(label)) for label in labels[n:])
            if n and gcd(*x) > 1:
                # Write v=(c*x0,r), x0 primitive. Pointwise fixation is
                # x0*A=x0, r*D=r (mod c), x0*B=(r-r*D)/c.
                # The GL(R) congruence stabilizer has index at most c^k.
                # Use the owned *integral* finite-index Schreier algorithm,
                # not the finite residue-image stabilizer: the latter
                # omits the infinite congruence kernel. The free shear
                # kernel is Hom(C/Zx0,R); -id_L adds the line reversal.
                c = gcd(*x)
                x0 = tuple(entry // c for entry in x)
                k = size - n
                complement = reduced.lattice_category()([
                    [int(gram[i, j]) for j in range(n)] for i in range(n)
                ])
                c_labels = tuple(complement.module_generating_set())
                projected = complement.linear_combination({
                    label: complement.base_ring()(entry)
                    for label, entry in zip(c_labels, x0, strict=True) if entry
                })
                _, _, change = matrix(SageZZ, [x0]).smith_form()
                blocks = []
                for generator in self.pointed_automorphism_generators(projected):
                    block = identity_matrix(SageZZ, size)
                    for i, point in enumerate(complement.module_generators()):
                        image = generator(point).to_vector()
                        for j, label in enumerate(c_labels):
                            block[i, j] = _integer(image(label))
                    blocks.append(block)
                linear = OwnedGroups().GL(k, reduced.base_ring())

                def radical_matrix(generator):
                    return linear._to_engine(generator).matrix()

                def preserves_residue(generator):
                    action = radical_matrix(generator)
                    return all((sum(r[i] * _integer(action[i, j]) for i in range(k)) - r[j]) % c == 0
                               for j in range(k))

                stabilizer = linear.predicate_subgroup(
                    preserves_residue, lambda: f"stabilizer of radical residue {r} modulo {c}"
                )
                for generator in stabilizer.finite_index_schreier_generators(
                    index_bound=c ** k, generators=linear.group_generators()
                ):
                    action = radical_matrix(generator)
                    block = identity_matrix(SageZZ, size)
                    for i in range(k):
                        for j in range(k):
                            block[n + i, n + j] = _integer(action[i, j])
                    for j in range(k):
                        delta = (r[j] - sum(r[i] * _integer(action[i, j]) for i in range(k))) // c
                        for row in range(n):
                            block[row, n + j] = int(change[row, 0]) * delta
                    blocks.append(block)
                for i in range(1, n):
                    for j in range(k):
                        block = identity_matrix(SageZZ, size)
                        for row in range(n):
                            block[row, n + j] = int(change[row, i])
                        blocks.append(block)
                blocks.append(-identity_matrix(SageZZ, size))
                return _conjugated(reduction.isometry, tuple(
                    _isometry_from_row_matrix(reduced, reduced, block) for block in blocks
                ))
            if n and gcd(*x) == 1:
                # For L=C+R, every isometry is [A B; 0 D]. If x is
                # primitive, choose u with x*u=1. The shear T=[I u*r;
                # 0 I] carries (x,0) to (x,r). Its line stabilizer is
                # T^-1 (Stab_O(C)(Zx) x GL(R) semidirect Hom(C/Zx,R)) T.
                # This is exact: xB=0 is the sole constraint on B for
                # (x,0) to span its own line.
                complement = reduced.lattice_category()([
                    [int(gram[i, j]) for j in range(n)] for i in range(n)
                ])
                c_labels = tuple(complement.module_generating_set())
                projected = complement.linear_combination({
                    label: complement.base_ring()(entry)
                    for label, entry in zip(c_labels, x, strict=True) if entry
                })
                _, _, change = matrix(SageZZ, [x]).smith_form()
                shear = identity_matrix(SageZZ, size)
                for i in range(n):
                    for j, entry in enumerate(r):
                        shear[i, n + j] = int(change[i, 0]) * entry
                blocks = []
                for generator in self.line_stabilizer_generators(projected):
                    block = identity_matrix(SageZZ, size)
                    for i, point in enumerate(complement.module_generators()):
                        image = generator(point).to_vector()
                        for j, label in enumerate(c_labels):
                            block[i, j] = _integer(image(label))
                    blocks.append(block)
                radical_rank = size - n
                general_linear = OwnedGroups().GL(radical_rank, reduced.base_ring())
                for generator in general_linear.group_generators():
                    action = general_linear._to_engine(generator).matrix()
                    block = identity_matrix(SageZZ, size)
                    for i in range(radical_rank):
                        for j in range(radical_rank):
                            block[n + i, n + j] = _integer(action[i, j])
                    blocks.append(block)
                for i in range(1, n):
                    for j in range(radical_rank):
                        block = identity_matrix(SageZZ, size)
                        for row in range(n):
                            block[row, n + j] = int(change[row, i])
                        blocks.append(block)
                reduced_generators = tuple(_isometry_from_row_matrix(
                    reduced, reduced, shear.inverse() * block * shear
                ) for block in blocks)
                return _conjugated(reduction.isometry, reduced_generators)
        action = self.isotropic_stabilizer_action((vector,), ((-vector,),))
        # The source supplies a stable containing module, so this is the
        # finite residue action (GroupAction.g:813-900), independently of
        # whether the rational arithmetic group itself is finite.
        return _integral_generator_set(lattice, action)

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
        r"""H, the actual integral parabolic image on K = v-perp/v.

        The isotropic-reduction owner computes the image, its membership
        predicate, and its lift section. Do not replace this with a subgroup
        generated by a potentially incomplete finite stabilizer resolution.
        """
        reduction = vector.parent().primitive_isotropic_subobject(vector).isotropic_reduction()
        return reduction.levi_image()

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

    def radical_stabilizer_congruence_image(self, lattice, modulus: int, fixed_row):
        r"""Actual mod-c image of the integral GL(R) stabilizer of a nonzero radical row.

        Its finite generators are reductions of complete integral stabilizer
        generators, with exact lifts retained for transporter words.  This
        is generally smaller than the stabilizer of the row modulo c in
        GL_k(Z/c), so using the latter would invent integral witnesses.
        """
        rank = len(fixed_row)
        lifts = self._radical_stabilizer_of_row(fixed_row)
        residue_ring = Zmod(modulus)
        generators = tuple(matrix(residue_ring, [[int(g[i, j]) % modulus for j in range(rank)]
            for i in range(rank)]) for g in lifts)
        full = self.radical_congruence_image(lattice, modulus)
        return RadicalCongruenceImage(
            modulus=modulus, radical_rank=rank,
            group=MatrixGroup(list(generators)), generators=generators,
            lifts=lifts, general_linear_order=full.general_linear_order,
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
        composed = correction * transport
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
        if vector == lattice.zero():
            # The zero vector imposes no stabilizer condition, including
            # in the zero lattice and in degenerate ambient forms.
            return self.orthogonal_group_generators(lattice)
        # Stab(c*v) = Stab(v) for every nonzero integer c: the
        # underlying free module is torsion-free.  In particular, do not
        # pass an imprimitive isotropic vector to a primitive line action.
        if int(vector.content()) > 1:
            return self.pointed_automorphism_generators(vector.primitive_part())
        if self.vector_lies_in_radical(vector):
            return self._degenerate_radical_pointed_generators(vector)
        if not lattice.is_nondegenerate() and vector.q() == lattice.base_ring().zero():
            # Central -id_L reverses every nonzero line generator. Thus
            # Stab(Zv) = Stab(v) x < -id_L >, and correcting each complete
            # line generator by its sign yields a complete generating
            # family for the pointwise stabilizer, without attempting a
            # Witt extension across the ambient radical.
            minus = _isometry_from_row_matrix(
                lattice, lattice, -identity_matrix(SageZZ, int(lattice.rank()))
            )
            result = []
            for generator in self.line_stabilizer_generators(vector):
                image = generator(vector)
                if image not in (vector, -vector):
                    raise ArithmeticError("a claimed line stabilizer generator does not preserve the marked line")
                result.append(generator if image == vector else minus * generator)
            return tuple(result)
        action = self.isotropic_stabilizer_action((vector,), ()) if vector.q() == lattice.base_ring().zero() else self._lattice_action(vector)
        return _integral_generator_set(lattice, action)

    def _mixed_degenerate_pointed_isometry(self, source: Lattices.ElementMethods, target: Lattices.ElementMethods, *, radical_action=None) -> LatticeIsometryMethods | None:
        r"""Construct a mixed radical/shear transporter in the reduced split presentation.

        For row coordinates ``(x,r)``, ``(y,s)`` on ``C+R``, every isometry
        has block form ``[[A,B],[0,D]]``.  First solve ``xA=y`` in O(C),
        then the exact finite-image congruence ``rD=s mod c``, where
        ``c=content(x)``.  A Bezout column ``u`` with ``x*u=c`` gives
        ``B=u*(s-rD)/c``.  The resulting matrix is integral unimodular,
        and its full form and marked-vector action are verified.
        """
        source_lattice, target_lattice = source.parent(), target.parent()
        if source_lattice is not target_lattice:
            transport = self.isometry(source_lattice, target_lattice)
            if transport is None:
                return None
            correction = self._mixed_degenerate_pointed_isometry(transport(source), target)
            if correction is None:
                return None
            return correction * transport
        lattice = source_lattice
        reduction = indefinite_reduced_presentation(lattice)
        reduced = reduction.reduced
        gram = _gram_in_label_order(reduced)
        n = int(gram.rank())
        k = int(gram.nrows()) - n
        if not k:
            raise ArithmeticError("mixed radical transport requires a nonzero radical")
        labels = tuple(reduced.module_generating_set())
        inverse = reduction.isometry.inverse()
        x = tuple(int(inverse(source).to_vector()(label)) for label in labels)
        y = tuple(int(inverse(target).to_vector()(label)) for label in labels)
        x_c, r = x[:n], x[n:]
        y_c, t = y[:n], y[n:]
        c = 0
        for value in x_c:
            c = gcd(c, value)
        if not c:
            return self._radical_pointed_isometry(source, target)
        complement = reduced.lattice_category()([[int(gram[i, j]) for j in range(n)] for i in range(n)])
        c_labels = tuple(complement.module_generating_set())
        x_vector = complement.linear_combination({label: complement.base_ring()(value)
            for label, value in zip(c_labels, x_c, strict=True) if value})
        y_vector = complement.linear_combination({label: complement.base_ring()(value)
            for label, value in zip(c_labels, y_c, strict=True) if value})
        a_isometry = self.pointed_isometry(x_vector, y_vector)
        if a_isometry is None:
            return None
        a = matrix(SageZZ, [[int(a_isometry(complement.module_generator(label)).to_vector()(key))
            for key in c_labels] for label in c_labels])
        d = (radical_action if radical_action is not None else
             (self.radical_congruence_image(lattice, c).transporter(r, t)
              if c > 1 else identity_matrix(SageZZ, k)))
        if d is None:
            return None
        residual = matrix(SageZZ, [list(t)]) - matrix(SageZZ, [list(r)]) * d
        if any(int(value) % c for value in residual.list()):
            if radical_action is not None:
                return None
            raise ArithmeticError("the radical congruence witness does not satisfy its residue equation")
        _, _, change = matrix(SageZZ, [list(x_c)]).smith_form()
        bezout = change.column(0)
        if sum(x_c[i] * int(bezout[i]) for i in range(n)) == -c:
            bezout = -bezout
        if sum(x_c[i] * int(bezout[i]) for i in range(n)) != c:
            raise ArithmeticError("Smith form did not produce the required Bezout column")
        block = identity_matrix(SageZZ, n + k)
        for i in range(n):
            for j in range(n):
                block[i, j] = a[i, j]
            for j in range(k):
                block[i, n + j] = int(bezout[i]) * (int(residual[0, j]) // c)
        for i in range(k):
            for j in range(k):
                block[n + i, n + j] = d[i, j]
        witness = _conjugated(reduction.isometry, (_isometry_from_row_matrix(reduced, reduced, block),))[0]
        if witness(source) != target:
            raise ArithmeticError("the mixed radical/shear witness does not transport the marking")
        return witness

    def pointed_isometry(self, source: Lattices.ElementMethods, target: Lattices.ElementMethods) -> LatticeIsometryMethods | None:
        r"""An element of Isom((L, v), (M, w)), or None."""
        if VectorPrefilter.from_vector(source) != VectorPrefilter.from_vector(target):
            return None
        # An integral isometry preserves content, and cancelling the
        # common nonzero content is equivalent over torsion-free lattices.
        source_content, target_content = int(source.content()), int(target.content())
        if source_content != target_content:
            return None
        if source_content > 1:
            return self.pointed_isometry(source.primitive_part(), target.primitive_part())
        source_lattice, target_lattice = source.parent(), target.parent()
        if source == source_lattice.zero():
            # Both vectors have content zero. Their pointed isometry
            # problem is precisely the unmarked lattice isometry problem.
            return self.isometry(source_lattice, target_lattice)
        if source_lattice is target_lattice and source == target:
            return source_lattice.O().identity()
        if source_lattice is target_lattice:
            # If q(x)=q(y), then for d=x-y one has
            # 2 b(x,d)=q(d), hence s_d(x)=y whenever q(d)!=0.
            # An integral reflection in the primitive line of d is
            # therefore a complete transporter certificate.  Test it
            # before constructing O(x^perp): in high-rank definite
            # complements the latter can exhaust PARI's qfauto stack even
            # when a single root reflection is already the answer.
            difference = source - target
            if difference.q() != source_lattice.base_ring().zero():
                root = difference.primitive_part()
                if root.is_integral_reflection_vector():
                    reflection = source_lattice.reflection(root)
                    if reflection(source) != target:
                        raise ArithmeticError("the integral difference reflection failed to transport the marking")
                    return reflection
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
        if not source_lattice.is_nondegenerate():
            return self._mixed_degenerate_pointed_isometry(source, target)
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
        target_lattice = target_reduction.isotropic_embedding().codomain()
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

    def isotropic_quotient_lift(
        self, reduction, reduced: LatticeIsometryMethods,
        *, radical_images: tuple[tuple[Lattices.ElementMethods, ...], ...] | None = None,
    ) -> LatticeIsometryMethods | None:
        r"""An integral parabolic lift of a prescribed quotient isometry, or None.

        Lorentzian.g:115-249 extends the integral map of the perpendicular.
        CombinedAlgorithms.h:779-877 supplies its triangular correction group
        with identity quotient action. The full GL(I), integral Hom(K,I)
        and saturated pointwise-perpendicular kernel preserve the prescribed
        reduction map throughout the finite integral-structure transporter.
        Supplied radical images restrict the correction action on I; the
        empty family gives pointwise lifts.
        """
        from sage_indefinite_port.indefinite.isotropic_induction import IsotropicInduction

        embedding = reduction.isotropic_embedding()
        lattice = embedding.codomain()
        isotropic = tuple(embedding(point) for point in embedding.domain().module_generators())
        rational_transport = self._reduced_isotropic_extension(reduction, reduction, isotropic, isotropic, reduced)
        if rational_transport.is_integral_on(lattice, lattice):
            return lattice.integral_isometry(rational_transport)
        if radical_images is None:
            subspace = lattice.primitive_isotropic_subobject(*isotropic)
            radical_images = IsotropicInduction.sublattice_action_images(subspace)
        (denominator,) = tuple(rational_transport.denominator_ideal(lattice, lattice).ideal_generators())
        correction_action = self.isotropic_stabilizer_action(
            isotropic, radical_images,
            quotient_generators=(), extension_denominator=_integer(denominator),
        )
        fraction_map = lattice.base_ring().fraction_field_map()
        restriction = Modules(fraction_map.codomain()).restriction_of_scalars(fraction_map)
        image = restriction(rational_transport) * ModulesOverIntegralDomains(lattice.base_ring())(lattice).generic_fibre_map()
        correction = correction_action.transporter(image)
        return None if correction is None else lattice.integral_isometry(correction.inverse() * rational_transport)

    def isotropic_quotient_kernel_generators(
        self, reduction, *, radical_images: tuple[tuple[Lattices.ElementMethods, ...], ...] | None = None,
    ) -> tuple[LatticeIsometryMethods, ...]:
        r"""Integral generators fixing the perpendicular quotient pointwise.

        CombinedAlgorithms.h:779-877 with the trivial quotient group retains
        GL(I), Hom(K,I) and the full pointwise-perpendicular kernel. The
        source integral-structure stabilizer selects their integral subgroup.
        A supplied radical action replaces GL(I); the empty family fixes I
        pointwise and gives the joint unipotent kernel.
        """
        from sage_indefinite_port.indefinite.isotropic_induction import IsotropicInduction

        embedding = reduction.isotropic_embedding()
        lattice = embedding.codomain()
        isotropic = tuple(embedding(point) for point in embedding.domain().module_generators())
        if not isotropic:
            return ()
        if radical_images is None:
            subspace = lattice.primitive_isotropic_subobject(*isotropic)
            radical_images = IsotropicInduction.sublattice_action_images(subspace)
        action = self.isotropic_stabilizer_action(
            isotropic, radical_images, quotient_generators=(),
        )
        return _integral_generator_set(lattice, action)

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
        integral = torsor.integral_base_point_or_none()
        if integral is not None:
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

    def exact_vector_orbit_iterator(self, lattice: Lattices.ParentMethods, square: int, *, primitive_only: bool = False):
        r"""Enumerate fixed-norm orbit representatives fairly without assuming finiteness.

        Every integral vector occurs in one finite coordinate-height shell.
        Exact pointed isometry decisions discard precisely those in an
        earlier orbit.  The iterator may run forever after its last yield
        when the orbit set is finite; it is not a finite-class certificate.
        """
        from itertools import count

        # These are decidably empty or singleton fibres.  In particular,
        # do not scan infinitely many nonzero vectors of an anisotropic
        # lattice when the prescribed norm is zero.
        if int(lattice.rank()) == 0:
            if int(square) == 0 and not primitive_only:
                yield lattice.zero()
            return
        if int(lattice.rank()) >= 3 and _gram_in_label_order(lattice).rank() == 1:
            # For a rank-one nondegenerate quotient and radical rank >= 2,
            # the reduced form is diag(a,0,...,0).  The full group is
            # {+/-1} x GL(R) extended by every shear Hom(Z,R).
            # For fixed |x|>0, the radical coordinates modulo |x| form
            # GL(R)-orbits indexed exactly by gcd(|x|,y_1,...,y_k),
            # i.e. the positive divisors of |x|.  When x=0 the radical
            # content is the complete orbit invariant.
            reduction = indefinite_reduced_presentation(lattice)
            reduced = reduction.reduced
            reduced_gram = _gram_in_label_order(reduced)
            coefficient = int(reduced_gram[0, 0])
            if coefficient == 0 or any(reduced_gram[i, j] for i in range(int(lattice.rank()))
                                         for j in range(int(lattice.rank())) if (i, j) != (0, 0)):
                raise ArithmeticError("rank-one quotient reduction did not isolate the radical")
            radical_zeros = (0,) * (int(lattice.rank()) - 2)
            if int(square) == 0:
                if primitive_only:
                    yield reduction.isometry(_from_column(reduced, (0, 1, *radical_zeros)))
                    return
                yield lattice.zero()
                for content in count(1):
                    yield reduction.isometry(_from_column(reduced, (0, content, *radical_zeros)))
                return
            quotient = int(square) // coefficient if int(square) % coefficient == 0 else 0
            if quotient <= 0 or isqrt(quotient) ** 2 != quotient:
                return
            x = isqrt(quotient)
            # Divisor enumeration is output-sensitive: scanning 1..x is
            # prohibitively expensive for a large represented square with
            # only a few divisor classes (e.g. x a large prime).
            for divisor in SageZZ(x).divisors():
                if primitive_only and divisor != 1:
                    continue
                yield reduction.isometry(_from_column(reduced, (x, divisor, *radical_zeros)))
            return
        if int(lattice.rank()) == 2 and _gram_in_label_order(lattice).rank() == 1:
            # In a reduced basis the form is diag(a,0), a != 0.  Its
            # integral orthogonal group consists exactly of
            # (x,y) -> (eps*x, b*x + delta*y), eps,delta=+/-1, b in ZZ.
            # For x != 0, orbits are |x| together with y mod |x| up to
            # sign; for x=0 they are |y|.  This gives a complete finite
            # nonzero-norm classification and an ordered null-content stream.
            reduction = indefinite_reduced_presentation(lattice)
            reduced = reduction.reduced
            reduced_gram = _gram_in_label_order(reduced)
            coefficient = int(reduced_gram[0, 0])
            if coefficient == 0 or any(reduced_gram[i, j] for i in range(2) for j in range(2) if (i, j) != (0, 0)):
                raise ArithmeticError("binary radical reduction did not produce diag(a,0)")

            def original_point(x: int, y: int):
                return reduction.isometry(_from_column(reduced, (x, y)))

            if int(square) == 0:
                if primitive_only:
                    yield original_point(0, 1)
                    return
                yield lattice.zero()
                for content in count(1):
                    yield original_point(0, content)
                return
            quotient = int(square) // coefficient if int(square) % coefficient == 0 else 0
            if quotient <= 0 or isqrt(quotient) ** 2 != quotient:
                return
            x = isqrt(quotient)
            for residue in range(x // 2 + 1):
                y = min(residue, (-residue) % x)
                if primitive_only and gcd(x, y) != 1:
                    continue
                yield original_point(x, y)
            return
        if int(lattice.rank()) == 2 and lattice.is_nondegenerate() and not lattice.is_definite() and lattice.is_isotropic():
            # The pinned binary factorization algorithm gives the complete
            # finite primitive orbit set, including the two null rays.
            # For nonprimitive vectors of nonzero square, content c obeys
            # c^2 | square; each primitive class at square/c^2 scales to
            # precisely one class of content c.  Null-vector contents are
            # unbounded and are enumerated in increasing order.
            from sage_indefinite_port.indefinite.lorentzian_cells import binary_split_primitive_vector_classes

            if primitive_only:
                yield from binary_split_primitive_vector_classes(lattice, int(square))
                return
            if int(square) == 0:
                yield lattice.zero()
                for content in count(1):
                    for representative in binary_split_primitive_vector_classes(lattice, 0):
                        yield lattice.scalar_multiple(lattice.base_ring()(content), representative)
                return
            for content in range(1, isqrt(abs(int(square))) + 1):
                if int(square) % (content * content):
                    continue
                for representative in binary_split_primitive_vector_classes(lattice, int(square) // (content * content)):
                    yield lattice.scalar_multiple(lattice.base_ring()(content), representative)
            return
        if int(lattice.rank()) == 2 and lattice.is_nondegenerate() and not lattice.is_definite():
            # The maintained binary quadratic-form orbit set computes a
            # *finite complete* transversal even for nonsplit indefinite
            # binary forms (Pell arithmetic), including imprimitive
            # vectors.  Unlike the generic coordinate walk, this branch
            # terminates after its last orbit representative.  The null
            # case was already decided by rational anisotropy above.
            orbit_set = lattice.binary_fixed_norm_representatives(int(square))
            section = orbit_set.section()
            for orbit in orbit_set:
                representative = section(orbit)
                if not primitive_only or int(representative.content()) == 1:
                    yield representative
            return
        profile = AttackProfile.from_lattice(lattice)
        if (lattice.is_nondegenerate() and profile.positive_index == 1
                and int(lattice.rank()) >= 3 and int(lattice.witt_index()) > 0
                and profile.sign * int(square) >= 0):
            # In the Lorentzian positive cone the perfect-cell algorithm
            # returns the complete finite primitive norm-orbit set.  The
            # null fibre is the induction algorithm's finite primitive
            # line-class set.  Recover all nonprimitive classes by content:
            # c^2 divides the norm, and scaling preserves orbit bijectively.
            if int(square) == 0:
                from sage_indefinite_port.indefinite.isotropic_induction import IsotropicInduction

                primitive = IsotropicInduction(self).isotropic_vector_classes(lattice)
                if primitive_only:
                    yield from primitive
                    return
                yield lattice.zero()
                if not primitive:
                    return
                for content in count(1):
                    for representative in primitive:
                        yield lattice.scalar_multiple(lattice.base_ring()(content), representative)
                return
            from sage_indefinite_port.indefinite.lorentzian_cells import (
                NONNEGATIVE_SQUARES,
                MarkedCellOrbitAlgorithm,
                perfect_complex,
            )

            reduction = indefinite_reduced_presentation(profile.signed_view)
            signed_norm = profile.sign * int(square)
            for content in range(1, isqrt(abs(int(square))) + 1):
                if int(square) % (content * content):
                    continue
                if primitive_only and content != 1:
                    continue
                primitive_norm = signed_norm // (content * content)
                orbits = MarkedCellOrbitAlgorithm(
                    perfect_complex(reduction.reduced, NONNEGATIVE_SQUARES), primitive_norm
                ).global_orbits()
                for orbit in orbits:
                    reduced_vector = orbit[0]
                    transported = reduction.isometry(reduced_vector)
                    coordinates = tuple(int(transported.to_vector()(label))
                                         for label in lattice.module_generating_set())
                    representative = _from_column(lattice, coordinates)
                    if not primitive_only or int(representative.content()) == 1:
                        yield lattice.scalar_multiple(lattice.base_ring()(content), representative)
            return
        if (lattice.is_nondegenerate() and profile.positive_index >= 2
                and int(lattice.witt_index()) >= 2):
            # The higher-Witt Eichler covering model computes a complete
            # finite transversal of primitive vector orbits at every
            # represented norm.  Lift it to all vectors by the exact
            # content decomposition: c^2 divides a nonzero norm, while
            # isotropic contents are unbounded and are streamed in order.
            if int(square) == 0:
                primitive = tuple(
                    vector_class.representative
                    for vector_class in self.primitive_vector_classes(lattice, 0)
                )
                if primitive_only:
                    yield from primitive
                    return
                yield lattice.zero()
                if not primitive:
                    return
                for content in count(1):
                    for representative in primitive:
                        yield lattice.scalar_multiple(lattice.base_ring()(content), representative)
                return
            for content in range(1, isqrt(abs(int(square))) + 1):
                if int(square) % (content * content):
                    continue
                if primitive_only and content != 1:
                    continue
                for vector_class in self.primitive_vector_classes(
                    lattice, int(square) // (content * content)
                ):
                    yield lattice.scalar_multiple(
                        lattice.base_ring()(content), vector_class.representative
                    )
            return
        if int(square) != 0 and not lattice.is_nondegenerate():
            # For L=C (+) R with nondegenerate C and zero radical
            # R, O(L) is the semidirect product of O(C) x GL(R) with
            # Hom(C,R).  A shear changes the radical component of (v,r)
            # by any multiple of content(v).  For rank(R)=1, residual
            # classes are r mod c up to sign; for rank(R)>=2, GL(R)
            # orbits of residual vectors are indexed by gcd(c,r), hence
            # by divisors of c.  This produces a *terminating* complete
            # norm-fibre classification whenever the quotient C admits
            # the existing finite exact fixed-norm classification.
            reduction = indefinite_reduced_presentation(lattice)
            reduced = reduction.reduced
            reduced_gram = _gram_in_label_order(reduced)
            quotient_rank = int(reduced_gram.rank())
            radical_rank = int(lattice.rank()) - quotient_rank
            if quotient_rank == 0:
                # A zero bilinear form represents no nonzero square;
                # constructing its rank-zero complement from [] would
                # also violate the owned lattice-constructor contract.
                return
            if any(reduced_gram[i, j] for i in range(quotient_rank, int(lattice.rank()))
                   for j in range(int(lattice.rank()))):
                raise ArithmeticError("degenerate reduction did not isolate the radical")
            quotient_gram = reduced_gram.matrix_from_rows_and_columns(
                range(quotient_rank), range(quotient_rank)
            )
            complement = lattice.lattice_category()(
                [[int(value) for value in row] for row in quotient_gram.rows()]
            )
            complement_profile = AttackProfile.from_lattice(complement)
            for content in range(1, isqrt(abs(int(square))) + 1):
                if int(square) % (content * content):
                    continue
                quotient_square = int(square) // (content * content)
                if complement.is_definite():
                    quotient_classes = tuple(
                        point for point in self.exact_vector_orbit_iterator(complement, quotient_square)
                        if int(point.content()) == 1
                    )
                elif int(complement.witt_index()) == 0:
                    # A real-indefinite but rationally anisotropic quotient
                    # has no Eichler envelope.  Its primitive vector orbits
                    # are primitive rank-one embeddings; use the exact
                    # genus/gluing embedding classifier rather than the
                    # approximate-model search or an unbounded coordinate
                    # walk.
                    source = lattice.lattice_category()([[quotient_square]])
                    embeddings = source.Emb(complement).primitive_embedding_class_representatives(
                        classification="emb"
                    )
                    source_generator = next(iter(source.module_generators()))
                    quotient_classes = tuple(embedding(source_generator) for embedding in embeddings)
                else:
                    # The public lattice operation owns the full primitive
                    # vector-class algorithm in every indefinite signature,
                    # including rationally anisotropic higher-rank quotients.
                    quotient_classes = tuple(
                        complement.O().vector_orbit_representatives(quotient_square, content=1)
                    )
                for quotient_point in quotient_classes:
                    coordinates = tuple(content * int(quotient_point.to_vector()(label))
                                         for label in complement.module_generating_set())
                    if radical_rank == 1:
                        residues = range(content // 2 + 1)
                    else:
                        residues = SageZZ(content).divisors()
                    for residue in residues:
                        if primitive_only and gcd(content, int(residue)) != 1:
                            continue
                        yield reduction.isometry(_from_column(
                            reduced, coordinates + (int(residue),) + (0,) * (radical_rank - 1)
                        ))
            return
        # On a zero-form lattice every vector has norm zero, and GL_n(ZZ)
        # acts transitively on primitive vectors: Bezout completion extends
        # each primitive column to an integral unimodular basis.  Thus the
        # primitive null fibre has exactly one orbit.  For arbitrary
        # vectors, content c>=0 is a complete GL_n(ZZ)-orbit invariant:
        # Bezout completion sends each vector to c times the first basis
        # vector.  Enumerate those classes directly rather than repeatedly
        # invoking the pointed transporter on an infinite coordinate walk.
        if all(
            lattice.b(x, y) == 0
            for x in lattice.module_generators()
            for y in lattice.module_generators()
        ):
            if int(square) != 0:
                return
            if primitive_only:
                yield lattice.module_generator(next(iter(lattice.module_generating_set())))
                return
            basis_vector = lattice.module_generator(next(iter(lattice.module_generating_set())))
            yield lattice.zero()
            for content in count(1):
                yield lattice.scalar_multiple(lattice.base_ring()(content), basis_vector)
            return
        if int(square) == 0 and not lattice.is_nondegenerate():
            # If the nondegenerate quotient C is rationally anisotropic,
            # every null vector lies in the radical R.  The exact
            # degenerate group contains GL(R), which is transitive on
            # vectors of each positive content, and fixes content.
            # Thus the entire null orbit space is {0} union {c*r:c>=1}
            # for one primitive radical generator r.  No infinite
            # coordinate search or pointed equivalence test is needed.
            reduction = indefinite_reduced_presentation(lattice)
            reduced = reduction.reduced
            gram_reduced = _gram_in_label_order(reduced)
            quotient_rank = int(gram_reduced.rank())
            quotient_gram = gram_reduced.matrix_from_rows_and_columns(
                range(quotient_rank), range(quotient_rank)
            )
            if any(gram_reduced[i, j] for i in range(quotient_rank, int(lattice.rank()))
                   for j in range(int(lattice.rank()))):
                raise ArithmeticError("degenerate reduction did not isolate the radical")
            complement = lattice.lattice_category()(
                [[int(value) for value in row] for row in quotient_gram.rows()]
            )
            if complement.is_definite() or not complement.is_isotropic():
                radical_point = reduction.isometry(_from_column(
                    reduced, (0,) * quotient_rank + (1,) + (0,) * (int(lattice.rank()) - quotient_rank - 1)
                ))
                if primitive_only:
                    yield radical_point
                    return
                yield lattice.zero()
                for content in count(1):
                    yield lattice.scalar_multiple(lattice.base_ring()(content), radical_point)
                return
            # For L=C (+) R and q(v,r)=0, the quotient vector v has
            # an invariant content c.  For c>0 the shears Hom(C,R)
            # reduce r modulo cR, and GL(R) gives residue classes
            # modulo sign when rank R=1, or gcd(c,r) when rank R>=2.
            # For c=0, only the radical content survives.  Both c
            # and the radical-only content are unbounded; interleave
            # their finite strata to obtain a fair complete stream.
            quotient_classes = tuple(
                complement.O().vector_orbit_representatives(0, content=1)
            )
            radical_rank = int(lattice.rank()) - quotient_rank
            radical_zeros = (0,) * (radical_rank - 1)
            radical_point = reduction.isometry(_from_column(
                reduced, (0,) * quotient_rank + (1,) + radical_zeros
            ))
            if not primitive_only:
                yield lattice.zero()
            for content in count(1):
                if not primitive_only or content == 1:
                    yield lattice.scalar_multiple(lattice.base_ring()(content), radical_point)
                for quotient_point in quotient_classes:
                    quotient_coordinates = tuple(int(quotient_point.to_vector()(label))
                                                 for label in complement.module_generating_set())
                    if radical_rank == 1:
                        residues = range(content // 2 + 1)
                    else:
                        residues = SageZZ(content).divisors()
                    for residue in residues:
                        if primitive_only and gcd(content, int(residue)) != 1:
                            continue
                        yield reduction.isometry(_from_column(
                            reduced,
                            tuple(content * value for value in quotient_coordinates)
                            + (int(residue),) + radical_zeros,
                        ))
            return
        if int(square) == 0 and lattice.is_nondegenerate() and not lattice.is_isotropic():
            if not primitive_only:
                yield lattice.zero()
            return
        # For an integral symmetric Gram matrix G, the integer values of
        # q(x)=xGx generate the ideal (G_ii, 2G_ij : i<j).  A requested
        # square outside that ideal has no integral representatives, even
        # when the real quadratic space is indefinite and isotropic.
        gram = _gram_in_label_order(lattice)
        if _norm_fibre_has_local_obstruction(gram, square, primitive=primitive_only):
            return
        # A definite quadratic form has a finite norm fibre.  Cauchy--
        # Schwarz for the positive definite form H=+/-G gives
        # x_i^2 <= |square|*(H^-1)_ii.  This provides an exact finite
        # search bound, rather than an iterator which hangs after its last
        # representative.  The bound is computed in QQ, without floats.
        if lattice.is_definite():
            sign = 1 if int(gram[0, 0]) > 0 else -1
            if sign * int(square) < 0:
                return
            inverse = (sign * gram).inverse()
            bounds = tuple(isqrt(int((abs(int(square)) * inverse[i, i]).floor()))
                           for i in range(int(lattice.rank())))
            # The finite group O(L) acts on the finite norm fibre.  Mark
            # whole orbits using its exact generators instead of invoking
            # a separate pointed-isometry decision for every pair of
            # candidates (quadratically many expensive group problems).
            from collections import deque

            generators = tuple(_column_matrix(g) for g in lattice.O().group_generators())
            seen: set[tuple[int, ...]] = set()
            for coordinates in itertools.product(*(range(-b, b + 1) for b in bounds)):
                if coordinates in seen:
                    continue
                point = _from_column(lattice, coordinates)
                if int(point.q()) != int(square) or (primitive_only and int(point.content()) != 1):
                    continue
                seen.add(coordinates)
                frontier = deque((coordinates,))
                while frontier:
                    current = frontier.popleft()
                    for generator in generators:
                        image = tuple(int(value) for value in generator * vector(SageZZ, current))
                        if image not in seen:
                            seen.add(image)
                            frontier.append(image)
                yield point
            return
        seen = []
        for height in count(0):
            for coordinates in _coordinate_box(lattice, height):
                point = _from_column(lattice, coordinates)
                if int(point.q()) != int(square) or (primitive_only and int(point.content()) != 1):
                    continue
                if all(self.pointed_isometry(previous, point) is None for previous in seen):
                    seen.append(point)
                    yield point

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
        # Reject impossible integral norms before requesting an Eichler
        # envelope.  In particular a rationally anisotropic rank-four
        # lattice may have no two hyperbolic planes at all, yet its odd
        # norm fibres can be certified empty by the integral value ideal.
        gram = _gram_in_label_order(lattice)
        if _norm_fibre_has_local_obstruction(gram, square, primitive=True):
            return ()

        # Rational Witt index one has no second hyperbolic pair for the
        # Eichler envelope (#40). A recognized rank-four Hilbert model decides
        # its classes, isotropic ones included, by the cusps of SL_2(O); every
        # other Witt-index-one lattice uses the proved Li-Margulis box, whose
        # isotropic orbits are its cusps.
        if int(lattice.witt_index()) == 1:
            from sage_indefinite_port.indefinite.li_margulis import primitive_vector_orbits
            from sage_indefinite_port.indefinite.witt_one_fixed_norm_classes import witt_one_fixed_norm_vector_classes
            from sage_indefinite_port.indefinite.witt_one_hilbert import recognized_model

            if int(lattice.module_rank()) == 4 and recognized_model(lattice) is not None:
                vectors = witt_one_fixed_norm_vector_classes(lattice, int(square), True)
            else:
                vectors = primitive_vector_orbits(lattice, int(square))
            return tuple(VectorClass(v, ()) for v in vectors)
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
