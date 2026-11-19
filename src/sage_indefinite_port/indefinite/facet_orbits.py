"""Orbit-representative facet enumeration for a perfect cell's cone.

`cell_cone` asks an engine for the **complete** H-description of the cell's ray cone and then keeps one facet per
stabilizer orbit (``lorentzian_cells.py:2834-2840``); on the A7 specimen that single ``SupportHyperplanes`` call
does not finish in 1 800 s. This module computes what the traversal consumes: the facets of ``cone(cell.rows)``,
grouped into orbits of the cell stabilizer, with exact normals, incidences and transport words. The reference's
``lorentzian_perfect.h:1433-1472`` calls ``DualDescriptionRecordFullDim(EXT_T, GRPperm, rddo)``. Its splitting
rules select direct exact dual description or adjacency decomposition: seed one facet, then walk the facet
graph across ridges, using the same dispatch on each face.

Every decision uses a documented exact operator:

* **seed** -- a vertex of the dual slice ``{u : <u,r> >= 0, <h,u> = 1}`` where ``h = sum(rays)`` is an interior
  point of the cone. This is a **linear-programming engine call** (cddlib's GMP backend), not a hand computation:
  the route removes the *complete H-description* call, not every engine call, and the seed is one LP per
  coordinate. The slice is bounded exactly because `h` is interior to the *primal* cone (a functional
  positive on the rays is not enough -- it leaves the slice unbounded), and its vertices are exactly the facet
  normals. Computed by lexicographic minimisation with cddlib's GMP backend: exact rational, no perturbation.
* **rank / kernel** -- ``matrix(...).rank()``, ``.right_kernel()``, ``.left_kernel_matrix()``,
  ``.row_space().basis_matrix()``.
* **ridges** -- the ridges of a facet are the facets of its cone inside its own span. The source's splitting
  rules select direct Normaliz computation or recursion on an injective coordinate projection of its rows.
* **ridge pencil** -- the covectors vanishing on a ridge form a 2-dimensional space, and the normal cone in it is
  a pointed 2-cone whose extreme rays are exactly the two facet normals through the ridge. Found by exact integer
  comparison of the constraint directions; no linear program, no tolerance.
* **certification** -- ``facet_certificates``: every incidence is checked with `is_facet_certificate`.

Completeness uses ridge adjacency and the pinned Balinski termination certificate. After one orbit has been
expanded, the remaining known facets cannot separate a processed region from an undiscovered region when their
full count is below the dual graph's connectivity, or when they all contain a common ray. Every discovered orbit
is returned, including those whose further expansion the certificate makes unnecessary. Otherwise the walk
expands an orbit of minimum facet incidence.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from fractions import Fraction
from heapq import heappop, heappush
from math import gcd

from sage.matrix.constructor import matrix
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.modules.free_module_element import vector
from sage.all import Permutation, PermutationGroup
from sage.libs.gap.libgap import libgap

from sage_indefinite_port.indefinite.facet_certificates import (
    _is_facet_certificate_from_pairings,
    is_facet_certificate,
    orient_facet_normal,
    pointedness_functional,
)

__all__ = [
    "FacetOrbitRecord",
    "cell_facet_orbits",
    "cone_facets",
    "orbit_representatives",
    "ray_permutations_from_stabilizer",
]


@dataclass(frozen=True)
class FacetOrbitRecord:
    """A facet representative and exact orbit size, with optional full transport words."""

    normal: tuple[int, ...]
    incidence: frozenset[int]
    orbit_size: int
    words: tuple[tuple[frozenset[int], tuple[int, ...]], ...]


def _primitive(entries) -> tuple[int, ...]:
    """Clear denominators, divide by the content, sign-normalise the first nonzero entry."""
    rational = [QQ(entry) for entry in entries]
    denominator = 1
    for entry in rational:
        denominator = denominator * int(entry.denominator()) // gcd(denominator, int(entry.denominator()))
    integral = [int(entry * denominator) for entry in rational]
    content = 0
    for entry in integral:
        content = gcd(content, entry)
    if content == 0:
        raise ArithmeticError("a zero covector has no primitive representative")
    primitive = [entry // content for entry in integral]
    first = next((entry for entry in primitive if entry), 0)
    if first < 0:
        primitive = [-entry for entry in primitive]
    return tuple(primitive)


def _pair(covector, ray) -> int:
    return sum(int(a) * int(b) for a, b in zip(covector, ray, strict=True))


def _same_ray(left, right) -> bool:
    r"""Whether two integer vectors are **positive** rational multiples of one another.

    The contract is proportionality against one nonzero pivot with a positive ratio, not a scan of adjacent
    minors. Scanning cyclically adjacent minors is not enough in dimension 4: ``e1`` and ``e3`` pass all four of
    those tests and are not proportional, and a vector and its negative pass every minor test while spanning
    *opposite* rays. For a pointed cone those are different rays, and confusing them corrupts the incidence
    bookkeeping.

    Comparison is over `Q` and exact: ``left = lambda * right`` with ``lambda > 0``.
    """
    left = [int(entry) for entry in left]
    right = [int(entry) for entry in right]
    if len(left) != len(right):
        return False
    pivot = next((index for index, entry in enumerate(right) if entry), None)
    if pivot is None:
        return False
    if all(entry == 0 for entry in left):
        return False
    # lambda = left[pivot] / right[pivot]; it must be positive and must reproduce every entry.
    numerator, denominator = left[pivot], right[pivot]
    if numerator == 0:
        return False
    if (numerator > 0) != (denominator > 0):
        return False
    return all(left[index] * denominator == numerator * right[index] for index in range(len(left)))


def _orient(rays, normal):
    """The non-negative normal, its incidence and all its exact ray pairings."""
    oriented, values = orient_facet_normal(normal, rays)
    return oriented, frozenset(index for index, value in enumerate(values) if value == 0), values


def _span_basis(rays):
    r"""An integral basis of the rational span, without integer-lattice reduction.

    Adjacency decomposition uses rational subspaces, not the module
    generated by the rays. A rational echelon basis followed by clearing
    each row's denominators gives the same span. Restriction to this basis
    is an invertible linear chart, so facet incidences and their supplied
    permutation actions are unchanged (lorentzian_perfect.h:1433-1472).
    """
    rational_basis = matrix(QQ, [list(ray) for ray in rays]).row_space().basis_matrix()
    return matrix(ZZ, [_primitive(row) for row in rational_basis.rows()])


def _seed_facet(rays, dimension: int) -> tuple[tuple[int, ...], frozenset[int], list[int]]:
    r"""A facet normal of ``cone(rays)``, from a vertex of the dual slice.

    Take ``h = sum(rays)``, which lies in the **interior** of the cone (the sum of all generators of a pointed
    cone spans its interior). Then ``S = {u : <u,r> >= 0 for all r, <h,u> = 1}`` is a bounded polytope -- this is
  the aperture condition, and it is why the normalising vector must be an interior point of the *primal* cone
    and not a functional positive on the rays: the latter leaves ``S`` unbounded, which is exactly how the first
    version of this routine failed (PPL: "Solution is unbounded"). The vertices of `S` are the extreme rays of
    the dual cone, i.e. exactly the facet normals. Lexicographic minimisation finds one vertex, in exact rational
    arithmetic through cddlib's GMP backend.
    """
    import cdd.gmp as cdd

    interior = [sum(int(ray[index]) for ray in rays) for index in range(dimension)]
    if all(entry == 0 for entry in interior):
        raise ArithmeticError("the sum of the rays vanished, so the cone cannot be pointed here")
    # cdd's H rows encode b + A*x >= 0. Mark the normalization and
    # subsequent lexicographic constraints as equalities. The GMP module
    # retains exact rational inputs and outputs:
    # https://pycddlib.readthedocs.io/en/stable/cdd.html#linear-programming
    rows = [[0, *(int(entry) for entry in ray)] for ray in rays]
    rows.append([-1, *interior])
    equalities = {len(rows) - 1}
    span_basis = _span_basis(rays)
    if span_basis.nrows() < dimension:
        annihilator = matrix(QQ, span_basis).right_kernel().basis_matrix()
        for row in range(annihilator.nrows()):
            rows.append([Fraction(0), *(Fraction(int(entry.numerator()), int(entry.denominator())) for entry in annihilator.row(row))])
            equalities.add(len(rows) - 1)
    coordinates = []
    for index in range(dimension):
        objective = [0] * (dimension + 1)
        objective[index + 1] = 1
        data = cdd.matrix_from_array(rows, lin_set=equalities, rep_type=cdd.RepType.INEQUALITY, obj_type=cdd.LPObjType.MIN, obj_func=objective)
        program = cdd.linprog_from_matrix(data)
        cdd.linprog_solve(program)
        if program.status != cdd.LPStatusType.OPTIMAL:
            raise ArithmeticError(f"the exact dual-slice program has status {program.status}")
        value = program.primal_solution[index]
        coordinates.append(QQ(value.numerator) / QQ(value.denominator))
        equality = [0] * (dimension + 1)
        equality[0], equality[index + 1] = -value, 1
        rows.append(equality)
        equalities.add(len(rows) - 1)
    return _orient(rays, _primitive(coordinates))


def _normal_cone_through_ridge(rays, ridge, span_basis):
    r"""The two facet normals whose incidence contains ``ridge``, by the 2-dimensional normal cone.

    The covectors annihilating the ridge span a 2-dimensional space `W`; the normal cone
    ``{w in W : <w,r_i> >= 0}`` is a pointed 2-cone whose extreme rays are exactly the two facet normals through
    the ridge. Every extreme ray of a pointed 2-cone is orthogonal to one of its boundary constraints, so the
    candidate set is complete and each candidate is tested exactly.
    """
    dimension = span_basis.ncols()
    if ridge:
        gram = span_basis * matrix(ZZ, [list(rays[index]) for index in sorted(ridge)]).transpose()
        # The pencil is a rational subspace; an integral kernel would
        # additionally saturate a lattice that this operation never uses.
        kernel = gram.change_ring(QQ).left_kernel_matrix()
        parameters = [[QQ(entry) for entry in kernel.row(index)] for index in range(kernel.nrows())]
    else:
        parameters = [[QQ(1 if index == position else 0) for position in range(span_basis.nrows())]
                      for index in range(span_basis.nrows())]
    if len(parameters) != 2:
        raise ArithmeticError(f"a ridge must have a two-dimensional pencil, got {len(parameters)}")
    basis = [
        _primitive([sum(theta[i] * QQ(span_basis[i, j]) for i in range(span_basis.nrows()))
                    for j in range(dimension)])
        for theta in parameters
    ]
    left, right = basis
    projections = []
    for ray in rays:
        first, second = _pair(left, ray), _pair(right, ray)
        if first or second:
            projections.append((first, second))
    if not projections:
        raise ArithmeticError("a ridge pencil must see a non-orthogonal ray")
    angular_max = angular_min = projections[0]
    for candidate in projections[1:]:
        if angular_max[0] * candidate[1] - angular_max[1] * candidate[0] > 0:
            angular_max = candidate
        if angular_min[0] * candidate[1] - angular_min[1] * candidate[0] < 0:
            angular_min = candidate
    found: list[tuple[int, ...]] = []
    for a, b in ((angular_max[1], -angular_max[0]), (-angular_min[1], angular_min[0])):
        if all(a * c + b * d >= 0 for c, d in projections):
            candidate = _primitive([a * QQ(p) + b * QQ(q) for p, q in zip(left, right, strict=True)])
            if not any(_same_ray(candidate, known) for known in found):
                found.append(candidate)
    if len(found) != 2:
        raise ArithmeticError(f"a ridge must lie in exactly two facets, found {len(found)}")
    return found


def _local_coordinates(rays, incidence):
    """The facet's rows in an injective coordinate projection of its span.

    Returns ``(local_rays, local_dimension)``.
    """
    face_rows = matrix(ZZ, [list(rays[index]) for index in sorted(incidence)])
    # The row-echelon basis has an invertible submatrix in its pivot
    # columns. Projection to those columns is therefore an isomorphism on
    # the rational span: it preserves all faces, incidences and supplied
    # ray permutations in the recursive adjacency decomposition
    # (lorentzian_perfect.h:1453-1454). Its coordinates are entries of the
    # original rays, so recursive charts introduce no coefficient growth.
    pivots = face_rows.pivots()
    coordinates = face_rows.matrix_from_columns(pivots)
    return [tuple(int(entry) for entry in row) for row in coordinates.rows()], len(pivots)


def cone_facets(rays) -> dict[frozenset[int], tuple[int, ...]]:
    r"""All facets of ``cone(rays)``, as ``incidence -> primitive normal``, by ridge closure.

    The cone must be full-dimensional and pointed (a perfect cell's ray cone is: see
    `facet_certificates.pointedness_functional`). The result is complete by the closure argument: it contains a
    certified seed and is closed under ridge adjacency.
    """
    rays = [tuple(int(entry) for entry in ray) for ray in rays]
    dimension = len(rays[0])
    if dimension == 1:
        raise ArithmeticError("a one-dimensional cone has no facet in this sense")

    # DualDescriptionRecordFullDim retains one cone representation while
    # walking its facets (lorentzian_perfect.h:1453-1454). The rational
    # span depends only on these rays, so all ridge pencils share it.
    span_basis = _span_basis(rays)
    if dimension == 2:
        normals = _normal_cone_through_ridge(rays, frozenset(), span_basis)
        facets = {}
        for candidate in normals:
            normal, incidence, values = _orient(rays, candidate)
            if not _is_facet_certificate_from_pairings(normal, incidence, rays, values):
                raise ArithmeticError("a two-dimensional facet failed certification")
            facets[incidence] = normal
        return facets

    seed_normal, seed_incidence, values = _seed_facet(rays, dimension)
    if not _is_facet_certificate_from_pairings(seed_normal, seed_incidence, rays, values):
        raise ArithmeticError("the seed facet failed certification")
    facets: dict[frozenset[int], tuple[int, ...]] = {seed_incidence: seed_normal}
    pending = [seed_incidence]
    while pending:
        current = pending.pop()
        normal = facets[current]
        for ridge in _ridges(rays, current):
            for candidate in _normal_cone_through_ridge(rays, ridge, span_basis):
                if _same_ray(candidate, normal):
                    continue
                oriented, incidence, values = _orient(rays, candidate)
                if incidence in facets:
                    continue
                if not _is_facet_certificate_from_pairings(oriented, incidence, rays, values):
                    raise ArithmeticError(
                        "a facet adjacent across a ridge failed certification (incidence size "
                        f"{len(incidence)}); excluding it is not proved, so the walk refuses to continue"
                    )
                facets[incidence] = oriented
                pending.append(incidence)
    return facets


def _ridges(rays, incidence):
    """The ridges of one facet: the facets of its cone inside its own span, mapped back to ambient incidences."""
    local_rays, local_dimension = _local_coordinates(rays, incidence)
    if local_dimension < 2:
        return ()
    local_facets = cone_facets(local_rays)
    positions = sorted(incidence)
    return tuple(frozenset(positions[index] for index in ridge) for ridge in local_facets)


def _integral_coordinate(value) -> int:
    r"""The exact integer value of an owned coordinate, refusing a rational one.

    The generator acts by an integral matrix, so every coordinate of an image ray must be an integer. Converting
    with ``int(QQ(numerator)/QQ(denominator))`` truncates instead: a rational coordinate such as ``3/2`` becomes
    ``1`` and could coincide with a row it has no right to, fabricating a permutation of the rays.
    """
    numerator = int(value.numerator())
    denominator = int(value.denominator())
    if denominator == 0:
        raise ArithmeticError("an owned coordinate has a zero denominator")
    if numerator % denominator:
        raise ArithmeticError(
            f"a stabilizer generator produced the non-integral coordinate {numerator}/{denominator}; the "
            "generators must act integrally on the rays"
        )
    return numerator // denominator


def ray_permutations_from_stabilizer(cell, stabilizer) -> tuple[tuple[int, ...], ...]:
    """Position permutations of the cell's rows induced by the stabilizer generators; raises unless they permute."""
    rays = [tuple(int(entry) for entry in row) for row in cell.rows]
    positions = {ray: index for index, ray in enumerate(rays)}
    if len(positions) != len(rays):
        raise ArithmeticError("the cell configuration contains repeated rays")
    labels = tuple(cell.lattice.module_generating_set())
    permutations = []
    for generator in stabilizer:
        image = []
        for point in cell.vector_configuration:
            moved = generator(point).to_vector()
            row = tuple(_integral_coordinate(moved(label)) for label in labels)
            if row not in positions:
                raise ArithmeticError("a stabilizer generator did not permute the cell's rays")
            image.append(positions[row])
        if sorted(image) != list(range(len(rays))):
            raise ArithmeticError("a stabilizer generator is not a permutation of the rays")
        permutations.append(tuple(image))
    return tuple(permutations)


def _orbit(incidence, permutations):
    r"""Exact orbit and Schreier words, with bitset actions on large ray sets.

    Every generator is an actual permutation of a finite ray universe.
    Replace an incidence F by its characteristic integer bitset and each
    generator by a support mask plus the image bit of each moved position.
    Applying it to F changes ONLY the occupied moved positions, retaining
    fixed positions in one integer AND. No abstract group equivalence or
    approximate hash replaces the orbit membership test.

    Breadth-first insertion and generator order are identical to the
    ordinary subset BFS; the returned map of incidences to transporter
    words is therefore exactly the same. For a sparse-support generator
    (such as a coordinate swap of a high-dimensional symmetric cell),
    this costs O(number of moved rays in F), not O(|F|).
    """
    if permutations and len(permutations) >= 8:
        degree = len(permutations[0])
        if all(len(generator) == degree for generator in permutations):
            images = []
            for generator in permutations:
                support_mask = 0
                bit_images = [0] * degree
                for position, image in enumerate(generator):
                    if not 0 <= int(image) < degree:
                        raise ArithmeticError("a ray stabilizer action contains an invalid permutation image")
                    if image != position:
                        support_mask |= 1 << position
                        bit_images[position] = 1 << int(image)
                images.append((support_mask, tuple(bit_images)))

            initial_mask = sum(1 << int(position) for position in incidence)
            words = {initial_mask: ()}
            pending = deque([initial_mask])
            while pending:
                current = pending.popleft()
                word = words[current]
                for index, (support_mask, generator_images) in enumerate(images):
                    moved = current & ~support_mask
                    affected = current & support_mask
                    while affected:
                        singleton = affected & -affected
                        moved |= generator_images[singleton.bit_length()-1]
                        affected ^= singleton
                    if moved not in words:
                        words[moved] = word + (index,)
                        pending.append(moved)

            def subset(mask):
                positions = []
                while mask:
                    singleton = mask & -mask
                    positions.append(singleton.bit_length() - 1)
                    mask ^= singleton
                return frozenset(positions)
            return {subset(mask): word for mask, word in words.items()}

    members = {incidence: ()}
    queue = deque([incidence])
    while queue:
        current = queue.popleft()
        word = members[current]
        for index, permutation in enumerate(permutations):
            moved = frozenset(permutation[position] for position in current)
            if moved not in members:
                members[moved] = word + (index,)
                queue.append(moved)
    return members


def _orbit_members(incidence, permutations):
    r"""The orbit of a facet incidence under the ray permutations, as a set."""
    members = {incidence}
    queue = deque([incidence])
    while queue:
        current = queue.popleft()
        for permutation in permutations:
            moved = frozenset(permutation[position] for position in current)
            if moved not in members:
                members.add(moved)
                queue.append(moved)
    return members


def _stabilizer_on_facet(group, incidence, orbit_size, *, orbit_words=None, permutations=()):
    r"""The induced action of the exact set stabilizer on its facet.

    Sage's PermutationGroup.stabilizer(..., action="OnSets") delegates
    the finite set stabilizer to GAP. Restricting its generators gives
    the full image on the facet, including when the restriction has a
    kernel. Global transport words still use the supplied ray generators.
    """
    # The supplied orbit words are already an exact Schreier transversal.
    # They determine the full subgroup without asking GAP to calculate
    # the order of the possibly enormous ambient permutation group.
    if orbit_words is not None and not permutations:
        return ()
    positions = tuple(sorted(incidence))
    local_position = {position + 1: index for index, position in enumerate(positions)}
    if orbit_words is not None and permutations:
        # Schreier's lemma for the action on subsets: if t_x takes the
        # reference facet F to x, then t_x s t_(xs)^-1 fixes F. The BFS
        # transversal is already available from the exact facet-orbit
        # computation. Restrict these Schreier generators to F instead of
        # invoking GAP's second, potentially costly OnSets stabilizer search.
        # Evaluate each transversal word on the facet only once while it
        # remains in a bounded cache. A Schreier edge then takes one
        # generator application and one inverse lookup per facet position.
        transport_cache = {}
        cache_limit = max(2, 250_000 // max(1, len(positions)))
        local_generators = set()
        identity = tuple(range(len(positions)))
        for member, word in orbit_words.items():
            source_data = transport_cache.get(member)
            if source_data is None:
                source_images = list(positions)
                for step in word:
                    generator = permutations[step]
                    source_images = [generator[position] for position in source_images]
                if set(source_images) != member:
                    raise ArithmeticError("a facet Schreier word failed to transport its incidence")
                source_data = (tuple(source_images), {image: index for index, image in enumerate(source_images)})
                transport_cache[member] = source_data
                if len(transport_cache) > cache_limit:
                    transport_cache.pop(next(iter(transport_cache)))
            source_images = source_data[0]
            for generator_index, generator in enumerate(permutations):
                next_member = frozenset(generator[position] for position in member)
                target_data = transport_cache.get(next_member)
                if target_data is None:
                    target_images = list(positions)
                    for step in orbit_words[next_member]:
                        target_generator = permutations[step]
                        target_images = [target_generator[position] for position in target_images]
                    if set(target_images) != next_member:
                        raise ArithmeticError("a facet Schreier target word failed to transport its incidence")
                    target_data = (tuple(target_images), {image: index for index, image in enumerate(target_images)})
                    transport_cache[next_member] = target_data
                    if len(transport_cache) > cache_limit:
                        transport_cache.pop(next(iter(transport_cache)))
                action = tuple(target_data[1][generator[position]] for position in source_images)
                if action != identity:
                    local_generators.add(action)
        return tuple(sorted(local_generators))
    # The GAP orbit-stabilizer shortcut is needed only on the
    # non-Schreier path. No group order was computed on the exact
    # word-transversal path above.
    if orbit_size == int(group.order()):
        return ()
    stabilizer = group.stabilizer([position + 1 for position in positions], action="OnSets")
    return tuple(
        tuple(local_position[int(generator(position + 1))] for position in positions)
        for generator in stabilizer.gens()
    )


def orbit_representatives(rays, permutations, *, level: int = 1):
    r"""The certified incidences and normals of the complete facet-orbit walk."""
    return tuple(
        (record.incidence, record.normal)
        for record in _facet_orbit_records(rays, permutations, level=level)
    )


def _crosspolytope_facet_orbit(rays, permutations, *, retain_words):
    r"""Certify the *complete* facet orbit of an affine crosspolytope.

    Suppose the configuration has 2m+1 distinct rows in rank m+1,
    consisting of a centre c and m pairs c+v_i, c-v_i. If the v_i are
    linearly independent and c is outside their span, the affine section
    is exactly an m-dimensional crosspolytope. Its facets are precisely
    the 2^m choices of one endpoint from each pair.

    Require the supplied integral cell-stabilizer permutations to contain
    each independent antipodal-pair exchange, fixing every other ray.
    These verified generators act transitively on all those facets.
    Consequently ONE representative and its 2^m-member orbit are the
    *entire* facet set modulo the actual supplied group, with each
    Schreier word built from the verified exchanges. No facet can be
    omitted and no unrelated rational symmetry is assumed.

    This avoids a 65,536-facet Normaliz dual description and the
    subsequent recursive/GAP backtracking for the initial rank-17
    IsotropicCases_8819 perfect cell.
    """
    dimension = len(rays[0])
    count = len(rays)
    if dimension < 8 or count != 2 * dimension - 1:
        return None
    if len(set(rays)) != count:
        return None
    if any(len(row) != dimension for row in rays):
        return None
    if not permutations:
        return None

    total = tuple(sum(row[j] for row in rays) for j in range(dimension))
    if any(value % count for value in total):
        return None
    centre = tuple(value // count for value in total)
    positions = {row: i for i, row in enumerate(rays)}
    centre_index = positions.get(centre)
    if centre_index is None:
        return None
    remaining = set(positions) - {centre}
    pairs = []
    while remaining:
        first = min(remaining)
        second = tuple(2 * c - x for c, x in zip(centre, first, strict=True))
        if second == first or second not in remaining:
            return None
        remaining.remove(first)
        remaining.remove(second)
        pairs.append((positions[first], positions[second]))
    if len(pairs) != dimension - 1:
        return None
    differences = matrix(QQ, [
        [rays[first][j] - centre[j] for j in range(dimension)]
        for first, _second in pairs
    ])
    if differences.rank() != dimension - 1:
        return None
    if matrix(QQ, [centre, *differences.rows()]).rank() != dimension:
        return None

    pair_swaps = []
    for first, second in pairs:
        index = next(
            (
                index for index, permutation in enumerate(permutations)
                if len(permutation) == count
                and permutation[first] == second
                and permutation[second] == first
                and all(permutation[position] == position
                        for position in range(count) if position not in (first, second))
            ),
            None,
        )
        if index is None:
            return None
        pair_swaps.append(index)

    incidence = frozenset(first for first, _second in pairs)
    kernel = matrix(QQ, [rays[position] for position in sorted(incidence)]).right_kernel_matrix()
    if kernel.nrows() != 1:
        return None
    normal, certified_incidence, values = _orient(rays, _primitive(kernel.row(0)))
    if certified_incidence != incidence or not _is_facet_certificate_from_pairings(
        normal, incidence, rays, values,
    ):
        raise ArithmeticError("the crosspolytope's purported facet failed its exact support certificate")

    if not retain_words:
        return FacetOrbitRecord(normal, incidence, 1 << len(pairs), ())
    # Independent, disjoint pair swaps commute. Their indexed words
    # send the chosen facet to every other sign choice exactly once.
    words = []
    for mask in range(1 << len(pairs)):
        member = frozenset(
            second if mask & (1 << i) else first
            for i, (first, second) in enumerate(pairs)
        )
        word = tuple(pair_swaps[i] for i in range(len(pairs)) if mask & (1 << i))
        words.append((member, word))
    return FacetOrbitRecord(normal, incidence, len(words), tuple(words))


def _symmetric_chamber_facet_orbits(rays, permutations, *, retain_words):
    r"""All facet orbits from a certified symmetric fundamental chamber.

    Write K=cone(rays) for a full-dimensional pointed rational cone in Q^d.
    An integral permutation of coordinate positions that permutes the
    supplied rays preserves K and acts dually on its cone K^vee of
    nonnegative supporting covectors. Only coordinate transpositions
    whose *entire ray permutation occurs among the supplied generators*
    are used; this does not enlarge the acting subgroup.

    The graph of these transpositions on coordinate positions generates
    the full symmetric group W on each connected component. For one
    component I choose the closed fundamental chamber
    C={n_i<=n_j for adjacent sorted i,j in I}. Every W-orbit of K^vee
    meets C. An extreme ray rho of K^vee that lies in C is also an
    extreme ray of K^vee intersect C (any decomposition in the
    intersection is a decomposition in K^vee). Hence the extreme rays
    of K^vee intersect C include a representative of *every* facet
    orbit of K under the supplied group.

    Normaliz enumerates those candidate rays exactly from the inequalities
    given by the original rays and chamber walls. An artificial extreme
    ray cut out by a chamber wall need not be extreme in K^vee, so test
    that its zero incidence has exact rank d-1. The accepted candidates
    are genuine facets, and their complete orbits and Schreier words are
    computed under the original supplied generator action, never under
    a guessed full orthogonal group. Thus no facet is omitted or wrongly
    identified, while the chamber may shrink an enormous dual
    description to a small exact list.

    The rank-17 205-ray source configuration has W=S_14 and its chamber
    has just 57 candidate normals, all certified facets; the unrestricted
    dual has far more. The procedure applies to arbitrary source
    configurations whenever the stated group certificate exists.
    """
    dimension = len(rays[0])
    if dimension < 8 or len(rays) < 3 * dimension or len(permutations) < 8:
        return None
    rows = tuple(tuple(int(entry) for entry in row) for row in rays)
    positions = {row: i for i, row in enumerate(rows)}
    if len(positions) != len(rows):
        return None
    supplied = {tuple(int(position) for position in permutation) for permutation in permutations}
    if any(len(permutation) != len(rows) for permutation in supplied):
        return None

    # For an edge i-j in the swap graph the corresponding position
    # permutation is literally a supplied generator. Connected edge
    # transpositions generate the entire symmetric group on each
    # component; choosing a sorted chamber is therefore legitimate
    # even when the edge set is not the complete graph.
    neighbors = [set() for _ in range(dimension)]
    for i in range(dimension):
        for j in range(i + 1, dimension):
            images = []
            for row in rows:
                swapped = row[:i] + (row[j],) + row[i + 1:j] + (row[i],) + row[j + 1:]
                image = positions.get(swapped)
                if image is None:
                    break
                images.append(image)
            if len(images) == len(rows) and tuple(images) in supplied:
                neighbors[i].add(j)
                neighbors[j].add(i)
    remaining = set(range(dimension))
    components = []
    while remaining:
        start = min(remaining)
        remaining.remove(start)
        component = {start}
        queue = deque([start])
        while queue:
            i = queue.popleft()
            for j in sorted(neighbors[i] & remaining):
                remaining.remove(j)
                component.add(j)
                queue.append(j)
        components.append(tuple(sorted(component)))
    largest = max(components, key=len)
    if len(largest) < 8:
        return None

    # The dual inequalities are dot(normal,ray)>=0; the chamber
    # inequalities sort its entries. This cone is pointed because K
    # spans Q^d, so its exact extreme rays are a complete finite set.
    from sage.geometry.polyhedron.backend_normaliz import NmzCone, NmzResult

    chamber_inequalities = []
    for left, right in zip(largest, largest[1:]):
        inequality = [0] * dimension
        inequality[left] = -1
        inequality[right] = 1
        chamber_inequalities.append(inequality)

    # If one coordinate sign change on this connected symmetric block
    # is a *supplied* generator, conjugation by the already certified
    # exchanges supplies every independent sign change. We have the
    # hyperoctahedral reflection group (Z/2)^m ⋊ S_m, and can further
    # restrict the dual normal to 0 <= n_1 <= ... <= n_m. This signed
    # Weyl chamber is a fundamental domain, so the earlier extreme-ray
    # completeness proof remains valid without enumerating all 2^m
    # sign copies. The frozen 2049-ray IsotropicCases_8812 configuration
    # has exactly this reflection subgroup.
    for index in largest:
        images = []
        for row in rows:
            reversed_row = row[:index] + (-row[index],) + row[index+1:]
            image = positions.get(reversed_row)
            if image is None:
                break
            images.append(image)
        if len(images) == len(rows) and tuple(images) in supplied:
            inequality = [0] * dimension
            inequality[largest[0]] = 1
            chamber_inequalities.append(inequality)
            break

    cone = NmzCone(inequalities=[*rows, *chamber_inequalities])
    candidates = NmzResult(cone, "ExtremeRays")
    del cone
    if not candidates:
        raise ArithmeticError("a pointed full-dimensional cone has no dual extreme rays")
    records = []
    covered = set()
    for candidate in candidates:
        normal, incidence, values = _orient(rows, _primitive(candidate))
        if not _is_facet_certificate_from_pairings(normal, incidence, rows, values):
            continue
        if incidence in covered:
            continue
        members = _orbit(incidence, permutations)
        covered.update(members)
        records.append(FacetOrbitRecord(
            normal, incidence, len(members),
            tuple(members.items()) if retain_words else (),
        ))
    if not records:
        raise ArithmeticError("the chamber-ray enumeration certified no facet of the pointed cone")
    return tuple(records)


def _facet_orbit_records(rays, permutations, *, level: int = 1, retain_words: bool = True):
    r"""One certified facet record per orbit, with transport words when requested.

    The source's splitting rules select direct dual description or adjacency decomposition. In the latter,
    the worklist holds one representative per orbit, and a
    representative contributes only the orbits of its own ridges, computed by the same routine one dimension down
    on the facet's cone modulo the facet's own stabilizer (Schreier). A new facet enters only if its orbit is
    disjoint from every orbit already found, so no orbit is duplicated. The source Balinski certificate or an
    empty worklist establishes completeness. Every discovered representative is yielded before further expansion.
    """
    rays = [tuple(int(entry) for entry in ray) for ray in rays]
    dimension = len(rays[0])
    if dimension < 2:
        raise ArithmeticError("a cone of dimension below two has no facets in this sense")
    certified_crosspolytope = _crosspolytope_facet_orbit(
        rays, permutations, retain_words=retain_words,
    )
    if certified_crosspolytope is not None:
        yield certified_crosspolytope
        return
    chamber_records = _symmetric_chamber_facet_orbits(
        rays, permutations, retain_words=retain_words,
    )
    if chamber_records is not None:
        yield from chamber_records
        return
    span_basis = _span_basis(rays)
    if retain_words:
        # Schreier words on the *actual* facet orbits give exact local
        # stabilizers. In this branch, constructing an ambient GAP group
        # and computing its gigantic order are mathematically unnecessary.
        # Point-orbit connected components need only the supplied action
        # on the finite ray index set, computed by ordinary generator
        # closure, not an ambient group presentation.
        group = None
        gap_group = None
        order = None
        remaining_points = set(range(len(rays)))
        point_orbits_list = []
        while remaining_points:
            start = min(remaining_points)
            orbit = {start}
            frontier = deque([start])
            while frontier:
                point = frontier.popleft()
                for permutation in permutations:
                    image = int(permutation[point])
                    if image not in orbit:
                        orbit.add(image)
                        frontier.append(image)
            remaining_points.difference_update(orbit)
            point_orbits_list.append(frozenset(orbit))
        point_orbits = tuple(point_orbits_list)
    else:
        group = PermutationGroup(
            [Permutation([int(position) + 1 for position in permutation]) for permutation in permutations],
            domain=tuple(range(1, len(rays) + 1)),
        )
        order = int(group.order())
        gap_group = libgap(group)
        point_orbits = tuple(frozenset(int(point) - 1 for point in orbit) for orbit in group.orbits())
    records: list[FacetOrbitRecord] = []
    covered: set[frozenset[int]] = set()
    # RepresentativeAction is asked only between facets of equal incidence
    # size; an incidence met again in the walk keeps its earlier decision.
    by_size: dict[int, list[FacetOrbitRecord]] = {}
    decided: set[frozenset[int]] = set()

    def register(normal, incidence):
        if retain_words:
            if incidence in covered:
                return None
            members = _orbit(incidence, permutations)
            covered.update(members)
            record = FacetOrbitRecord(normal, incidence, len(members), tuple(members.items()))
        else:
            # SetFunctionality.g:174-201,234-260: compare orbit
            # representatives by RepresentativeAction(..., OnSets), and
            # retain |G|/|Stab_G(F)| rather than enumerating the orbit.
            if incidence in decided or incidence in covered:
                return None
            decided.add(incidence)
            positions = sorted(position + 1 for position in incidence)
            key = len(incidence)
            for previous in by_size.get(key, ()):
                target = sorted(position + 1 for position in previous.incidence)
                if libgap.RepresentativeAction(gap_group, positions, target, libgap.OnSets) != libgap.eval("fail"):
                    return None
            # Only |Stab_G(F)| is read, so GAP's Stabilizer is not raised to
            # a Sage permutation group.
            stabilizer_order = int(libgap.Stabilizer(gap_group, positions, libgap.OnSets).Size())
            record = FacetOrbitRecord(normal, incidence, order // stabilizer_order, ())
            by_size.setdefault(key, []).append(record)
            # A listed orbit recognizes its later members by lookup instead
            # of by RepresentativeAction; either test decides the same orbit.
            # An orbit is listed when it has at most n d members, as many as
            # the n x d ray matrix has entries, so each list holds at most
            # n^2 d integers, polynomial in the input; a larger orbit keeps
            # the source's |G|/|Stab_G(F)| record.
            if record.orbit_size <= len(rays) * dimension:
                covered.update(_orbit_members(incidence, permutations))
        records.append(record)
        return record

    if dimension == 2:
        for candidate in _normal_cone_through_ridge(rays, frozenset(), span_basis):
            normal, incidence, values = _orient(rays, candidate)
            if not _is_facet_certificate_from_pairings(normal, incidence, rays, values):
                raise ArithmeticError("a two-dimensional facet failed certification")
            record = register(normal, incidence)
            if record is not None:
                yield record
        return
    # The pinned StandardHeuristicSplitting switches back to direct dual
    # description at deeper levels or for small groups. On the frozen
    # Lorentzian configurations that repeatedly submits large recursive
    # faces to Normaliz, even though the exact ridge walk is already
    # available. Use direct enumeration only for near-simplicial faces;
    # otherwise continue the certified recursive walk regardless of the
    # group's order or the depth. Each recursive face lowers dimension,
    # and the two-dimensional base case terminates exactly. This changes
    # the engine dispatch, not the facet set or its orbit equivalence.
    # The absolute excess-ray threshold misses the rank-17 crosspolytope
    # (33 rays, dimension 17), whose exact dual has 2^16 facets although
    # its excess ray count is only 16.  Direct Normaliz enumerates every
    # facet before quotienting by the full signed-permutation cell group.
    # An orbit-pruned recursive ridge walk instead proves completeness by
    # the same adjacency/Balinski certificate without expanding the whole
    # dual facet set. The choice changes only the algorithm, not its answer.
    split = (
        len(rays) - dimension >= 20
        or (dimension >= 8 and len(rays) >= 2 * dimension - 1)
    )
    if not split:
        # POLY_RecursiveDualDesc.h:1465-1474, DirectFacetOrbitComputation.
        # Use the exact Normaliz cone backend already owned by the local
        # perfect-cell engine boundary. It replaces the source's choice
        # between LRS and Normaliz while retaining every supporting row.
        from sage.geometry.polyhedron.backend_normaliz import NmzCone, NmzResult

        cone = NmzCone(cone=rays, subspace=())
        if NmzResult(cone, "Equations"):
            raise ArithmeticError("direct facet enumeration requires a full-dimensional cone")
        hyperplanes = NmzResult(cone, "SupportHyperplanes")
        # POLY_DirectDualDesc.h:334-337 rejects an empty direct result.
        if not hyperplanes:
            raise ArithmeticError("the direct facet computation returned no facets")
        direct_facets = {}
        for row in hyperplanes:
            normal, incidence, values = _orient(rays, _primitive(row))
            if not _is_facet_certificate_from_pairings(normal, incidence, rays, values):
                raise ArithmeticError("a direct facet failed certification")
            if not retain_words:
                direct_facets.setdefault(incidence, normal)
                continue
            record = register(normal, incidence)
            if record is not None:
                yield record
        if not retain_words:
            # Direct dual description has already supplied *every* facet.
            # The ray permutations therefore act on this finite set of
            # incidences. Enumerate its generated orbits by closure, not by
            # repeated RepresentativeAction/Stabilizer searches in the
            # potentially enormous group of ray permutations. This performs
            # one incidence image per (facet, supplied generator), regardless
            # of the order of the ambient group.
            #
            # A generator is a permutation of the rays, so its inverse is a
            # positive power; forward-generator closure is the full group
            # orbit. The complete normalised facet domain must be invariant.
            # Keep Normaliz's first facet for each orbit to preserve the
            # historical representative and traversal order.
            covered_facets: set[frozenset[int]] = set()
            for incidence, normal in direct_facets.items():
                if incidence in covered_facets:
                    continue
                orbit = {incidence}
                pending = deque((incidence,))
                while pending:
                    current = pending.popleft()
                    for permutation in permutations:
                        moved = frozenset(permutation[position] for position in current)
                        if moved not in direct_facets:
                            raise ArithmeticError("the stabilizer action did not preserve the complete facet domain")
                        if moved not in orbit:
                            orbit.add(moved)
                            pending.append(moved)
                covered_facets.update(orbit)
                yield FacetOrbitRecord(normal, incidence, len(orbit), ())
        return
    seed_normal, seed, values = _seed_facet(rays, dimension)
    if not _is_facet_certificate_from_pairings(seed_normal, seed, rays, values):
        raise ArithmeticError("the seed facet failed certification")
    seed_record = register(seed_normal, seed)
    assert seed_record is not None
    # SetFunctionality.g:304-323 and POLY_RecursiveDualDesc.h:837-855
    # expand a pending facet of minimum incidence, rather than the last one
    # discovered. Equal incidences retain their discovery order.
    pending = [(len(seed), 0)]
    processed_any = False
    yield seed_record
    while pending:
        if processed_any:
            # DualDescriptionAndADM.g:133-145; Balinski_basic.h:55-60.
            # Count full disjoint facet orbits. A common ray must belong
            # to every translate of every pending incidence, hence its
            # entire point orbit must lie in their intersection.
            if sum(records[index].orbit_size for _size, index in pending) <= dimension - 2:
                return
            common = set(records[pending[0][1]].incidence)
            for _size, index in pending[1:]:
                common.intersection_update(records[index].incidence)
                if not common:
                    break
            if common and any(orbit <= common for orbit in point_orbits):
                return
        _size, index = heappop(pending)
        record = records[index]
        incidence, normal = record.incidence, record.normal
        positions = sorted(incidence)
        local_permutations = _stabilizer_on_facet(
            group, incidence, record.orbit_size,
            orbit_words=dict(record.words) if retain_words else None,
            permutations=permutations,
        )
        local_rays, local_dimension = _local_coordinates(rays, incidence)
        if local_dimension < 2:
            processed_any = True
            continue
        for local_record in _facet_orbit_records(
            local_rays, local_permutations, level=level + 1, retain_words=retain_words
        ):
            ridge = frozenset(positions[offset] for offset in sorted(local_record.incidence))
            for candidate in _normal_cone_through_ridge(rays, ridge, span_basis):
                if _same_ray(candidate, normal):
                    continue
                oriented, adjacent, values = _orient(rays, candidate)
                if adjacent in covered:
                    continue
                if not _is_facet_certificate_from_pairings(oriented, adjacent, rays, values):
                    # Silence here would be a completeness claim built on a dropped facet: the worklist could
                    # empty with a genuine facet left unfound and the closure certificate would be false.
                    raise ArithmeticError(
                        "a facet adjacent across a ridge failed certification (incidence size "
                        f"{len(adjacent)}); excluding it is not proved, so the walk refuses to continue"
                    )
                adjacent_record = register(oriented, adjacent)
                if adjacent_record is not None:
                    heappush(pending, (len(adjacent), len(records) - 1))
                    yield adjacent_record
        processed_any = True


def cell_facet_orbits(cell, stabilizer) -> tuple[FacetOrbitRecord, ...]:
    r"""One certified record per orbit of ``stabilizer`` on the facets of ``cone(cell.rows)``.

    ``stabilizer`` is the cell's stabilizer as isometries (``cell_stabilizer(cell)``). The records come from
    `orbit_representatives`, with the source's direct or recursive dispatch. Each record carries a
    representative normal and incidence, its orbit size, and a word of generators transporting the representative
    onto every member of its orbit -- the transport witness a consumer needs to reconstruct any facet of the
    orbit exactly.
    """
    permutations = ray_permutations_from_stabilizer(cell, stabilizer)
    return _cell_facet_orbits_from_permutations(cell, permutations)


def _cell_facet_orbits_from_permutations(cell, permutations, *, retain_words: bool = True) -> tuple[FacetOrbitRecord, ...]:
    r"""Walk facet orbits using the certified ray action supplied by the caller.

    Lorentzian.g:1318-1322 passes GRPperm to adjacency decomposition.
    The local backend computes this action on the ambient basis and applies
    it to the ray rows; the public caller computes it on the ray elements.
    """
    rays = [tuple(int(entry) for entry in row) for row in cell.rows]
    _height, pairings = pointedness_functional(cell)
    if not all(value > 0 for value in pairings):
        raise ArithmeticError("the cell must be pointed: the sum of its rays has to be an interior point")

    records = []
    for record in _facet_orbit_records(rays, permutations, retain_words=retain_words):
        if not is_facet_certificate(record.normal, record.incidence, rays):
            raise ArithmeticError(f"the walk produced an uncertified facet with incidence of size {len(record.incidence)}")
        records.append(record)
    return tuple(records)
