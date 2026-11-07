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

Completeness is the closure argument, not a sample: the set produced is closed under ridge adjacency, and the
facet graph of a full-dimensional pointed cone is connected, so a closed non-empty set containing a seed is every
facet. Orbit grouping is then exact set-stabiliser bookkeeping.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from fractions import Fraction
from math import gcd

from sage.all import QQ, ZZ, Permutation, PermutationGroup, matrix, vector

from sage_indefinite_port.indefinite.facet_certificates import (
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
    """One stabilizer orbit of facets: a representative and its transport witnesses."""

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


def _incidence(rays, normal) -> frozenset[int]:
    return frozenset(index for index, ray in enumerate(rays) if _pair(normal, ray) == 0)


def _orient(rays, normal):
    """The sign of ``normal`` that is non-negative on every ray, plus its incidence."""
    oriented, values = orient_facet_normal(normal, rays)
    return oriented, frozenset(index for index, value in enumerate(values) if value == 0)


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


def _seed_facet(rays, dimension: int) -> tuple[int, ...]:
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
    normal, _incidence = _orient(rays, _primitive(coordinates))
    return normal


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
            normal, incidence = _orient(rays, candidate)
            if not is_facet_certificate(normal, incidence, rays):
                raise ArithmeticError("a two-dimensional facet failed certification")
            facets[incidence] = normal
        return facets

    seed_normal = _seed_facet(rays, dimension)
    seed_incidence = _incidence(rays, seed_normal)
    if not is_facet_certificate(seed_normal, seed_incidence, rays):
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
                oriented, incidence = _orient(rays, candidate)
                if incidence in facets:
                    continue
                if not is_facet_certificate(oriented, incidence, rays):
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


def _stabilizer_on_facet(group, incidence):
    r"""The induced action of the exact set stabilizer on its facet.

    Sage's PermutationGroup.stabilizer(..., action="OnSets") delegates
    the finite set stabilizer to GAP. Restricting its generators gives
    the full image on the facet, including when the restriction has a
    kernel. Global transport words still use the supplied ray generators.
    """
    positions = tuple(sorted(incidence))
    local_position = {position + 1: index for index, position in enumerate(positions)}
    stabilizer = group.stabilizer([position + 1 for position in positions], action="OnSets")
    return tuple(
        tuple(local_position[int(generator(position + 1))] for position in positions)
        for generator in stabilizer.gens()
    )


def orbit_representatives(rays, permutations, *, level: int = 1):
    r"""One certified ``(incidence, normal)`` per orbit, by the orbit-pruned adjacency decomposition.

    The source's splitting rules select direct dual description or adjacency decomposition. In the latter,
    the worklist holds one representative per orbit, and a
    representative contributes only the orbits of its own ridges, computed by the same routine one dimension down
    on the facet's cone modulo the facet's own stabilizer (Schreier). A new facet enters only if its orbit is
    disjoint from every orbit already found, so no orbit is duplicated, and termination with an empty worklist is
    the closure certificate, so none is missed.
    """
    rays = [tuple(int(entry) for entry in ray) for ray in rays]
    dimension = len(rays[0])
    if dimension < 2:
        raise ArithmeticError("a cone of dimension below two has no facets in this sense")
    span_basis = _span_basis(rays)
    if dimension == 2:
        representatives = []
        covered: set[frozenset[int]] = set()
        for candidate in _normal_cone_through_ridge(rays, frozenset(), span_basis):
            normal, incidence = _orient(rays, candidate)
            if not is_facet_certificate(normal, incidence, rays):
                raise ArithmeticError("a two-dimensional facet failed certification")
            if incidence in covered:
                continue
            covered.update(_orbit(incidence, permutations))
            representatives.append((incidence, normal))
        return tuple(representatives)

    group = PermutationGroup(
        [Permutation([int(position) + 1 for position in permutation]) for permutation in permutations],
        domain=tuple(range(1, len(rays) + 1)),
    )
    # POLY_Heuristics.h:320-328, StandardHeuristicSplitting, in rule order.
    # LORENTZ_GetGeneratorsAutom_Reduced sets dimEXT = n + 1, so the
    # initial n-dimensional cone has level one. Each facet lowers rank by
    # one. Both source dispatches compute the complete exact facet set;
    # these rules select computation, not a search bound.
    order = int(group.order())
    if len(rays) - dimension < 20:
        split = False
    elif order > 5000 and level < 3:
        split = True
    else:
        split = len(rays) >= 40 and level <= 1 and order >= 100
    if not split:
        # POLY_RecursiveDualDesc.h:1465-1474, DirectFacetOrbitComputation.
        # Use the exact Normaliz cone backend already owned by the local
        # perfect-cell engine boundary. It replaces the source's choice
        # between LRS and Normaliz while retaining every supporting row.
        from sage.geometry.polyhedron.backend_normaliz import NmzCone, NmzResult

        cone = NmzCone(cone=rays, subspace=())
        if NmzResult(cone, "Equations"):
            raise ArithmeticError("direct facet enumeration requires a full-dimensional cone")
        representatives = []
        covered = set()
        for row in NmzResult(cone, "SupportHyperplanes"):
            normal, incidence = _orient(rays, _primitive(row))
            if not is_facet_certificate(normal, incidence, rays):
                raise ArithmeticError("a direct facet failed certification")
            if incidence in covered:
                continue
            covered.update(_orbit(incidence, permutations))
            representatives.append((incidence, normal))
        return tuple(representatives)
    seed_normal = _seed_facet(rays, dimension)
    seed = _incidence(rays, seed_normal)
    if not is_facet_certificate(seed_normal, seed, rays):
        raise ArithmeticError("the seed facet failed certification")
    representatives: list[tuple[frozenset[int], tuple[int, ...]]] = [(seed, seed_normal)]
    covered = set(_orbit(seed, permutations))
    pending = [0]
    while pending:
        index = pending.pop()
        incidence, normal = representatives[index]
        positions = sorted(incidence)
        local_permutations = _stabilizer_on_facet(group, incidence)
        local_rays, local_dimension = _local_coordinates(rays, incidence)
        if local_dimension < 2:
            continue
        for ridge_local, _ridge_normal in orbit_representatives(local_rays, local_permutations, level=level + 1):
            ridge = frozenset(positions[offset] for offset in sorted(ridge_local))
            for candidate in _normal_cone_through_ridge(rays, ridge, span_basis):
                if _same_ray(candidate, normal):
                    continue
                oriented, adjacent = _orient(rays, candidate)
                if adjacent in covered:
                    continue
                if not is_facet_certificate(oriented, adjacent, rays):
                    # Silence here would be a completeness claim built on a dropped facet: the worklist could
                    # empty with a genuine facet left unfound and the closure certificate would be false.
                    raise ArithmeticError(
                        "a facet adjacent across a ridge failed certification (incidence size "
                        f"{len(adjacent)}); excluding it is not proved, so the walk refuses to continue"
                    )
                members = _orbit(adjacent, permutations)
                if members.keys() & covered:
                    continue
                covered.update(members)
                representatives.append((adjacent, oriented))
                pending.append(len(representatives) - 1)
    return tuple(representatives)


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


def _cell_facet_orbits_from_permutations(cell, permutations) -> tuple[FacetOrbitRecord, ...]:
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
    for incidence, normal in orbit_representatives(rays, permutations):
        if not is_facet_certificate(normal, incidence, rays):
            raise ArithmeticError(f"the walk produced an uncertified facet with incidence of size {len(incidence)}")
        members = _orbit(incidence, permutations)
        records.append(
            FacetOrbitRecord(
                normal=normal,
                incidence=incidence,
                orbit_size=len(members),
                words=tuple((member, word) for member, word in members.items()),
            )
        )
    return tuple(records)
