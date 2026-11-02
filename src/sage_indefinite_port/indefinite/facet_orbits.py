"""Orbit-representative facet enumeration for a perfect cell's cone.

`cell_cone` asks an engine for the **complete** H-description of the cell's ray cone and then keeps one facet per
stabilizer orbit (``lorentzian_cells.py:2834-2840``); on the A7 specimen that single ``SupportHyperplanes`` call
does not finish in 1 800 s. This module computes what the traversal consumes: the facets of ``cone(cell.rows)``,
grouped into orbits of the cell stabilizer, with exact normals, incidences and transport words. The facet set is
produced by the reference's adjacency decomposition (``polyhedral_common`` ``lorentzian_perfect.h:1433-1472``:
``f_adj`` -> ``DualDescriptionRecordFullDim(EXT_T, GRPperm, rddo)``): seed one facet, then walk the facet graph
across ridges.

Every decision uses a documented exact operator:

* **seed** -- a vertex of the dual slice ``{u : <u,r> >= 0, <h,u> = 1}`` where ``h = sum(rays)`` is an interior
  point of the cone. The slice is bounded exactly because `h` is interior to the *primal* cone (a functional
  positive on the rays is not enough -- it leaves the slice unbounded), and its vertices are exactly the facet
  normals. Computed by lexicographic minimisation with Sage's PPL backend: exact rational, no perturbation.
* **rank / kernel** -- ``matrix(...).rank()``, ``.right_kernel()``, ``.left_kernel_matrix()``,
  ``.row_space().basis_matrix()``.
* **ridges** -- the ridges of a facet are the facets of that facet's cone inside its own span, so the routine
  recurses one dimension down, on the facet's rows written in a basis of its span.
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
from math import gcd

from sage.all import QQ, ZZ, matrix, vector

from sage_indefinite_port.indefinite.facet_certificates import (
    is_facet_certificate,
    orient_facet_normal,
    pointedness_functional,
)

__all__ = ["FacetOrbitRecord", "cell_facet_orbits", "cone_facets", "ray_permutations_from_stabilizer"]


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
    return all(int(left[i]) * int(right[(i + 1) % len(left)]) == int(left[(i + 1) % len(left)]) * int(right[i])
               for i in range(len(left)))


def _incidence(rays, normal) -> frozenset[int]:
    return frozenset(index for index, ray in enumerate(rays) if _pair(normal, ray) == 0)


def _orient(rays, normal):
    """The sign of ``normal`` that is non-negative on every ray, plus its incidence."""
    oriented, values = orient_facet_normal(normal, rays)
    return oriented, frozenset(index for index, value in enumerate(values) if value == 0)


def _seed_facet(rays, dimension: int) -> tuple[int, ...]:
    r"""A facet normal of ``cone(rays)``, from a vertex of the dual slice.

    Take ``h = sum(rays)``, which lies in the **interior** of the cone (the sum of all generators of a pointed
    cone spans its interior). Then ``S = {u : <u,r> >= 0 for all r, <h,u> = 1}`` is a bounded polytope -- this is
    the aperture condition, and it is why the normalising vector must be an interior point of the *primal* cone
    and not a functional positive on the rays: the latter leaves ``S`` unbounded, which is exactly how the first
    version of this routine failed (PPL: "Solution is unbounded"). The vertices of `S` are the extreme rays of
    the dual cone, i.e. exactly the facet normals. Lexicographic minimisation finds one vertex, in exact rational
    arithmetic through Sage's PPL backend.
    """
    from sage.numerical.mip import MixedIntegerLinearProgram

    interior = [sum(int(ray[index]) for ray in rays) for index in range(dimension)]
    if all(entry == 0 for entry in interior):
        raise ArithmeticError("the sum of the rays vanished, so the cone cannot be pointed here")
    program = MixedIntegerLinearProgram(maximization=False, solver="PPL")
    variables = program.new_variable(real=True, indices=range(dimension))
    program.add_constraint(sum(ZZ(interior[index]) * variables[index] for index in range(dimension)) == 1)
    for ray in rays:
        program.add_constraint(sum(ZZ(ray[index]) * variables[index] for index in range(dimension)) >= 0)
    span_basis = matrix(ZZ, [list(ray) for ray in rays]).row_space().basis_matrix()
    if span_basis.nrows() < dimension:
        annihilator = matrix(QQ, span_basis).right_kernel().basis_matrix()
        for row in range(annihilator.nrows()):
            program.add_constraint(
                sum(QQ(annihilator[row, index]) * variables[index] for index in range(dimension)) == 0
            )
    coordinates = []
    for index in range(dimension):
        program.set_objective(variables[index])
        program.solve()
        value = QQ(program.get_values(variables)[index])
        coordinates.append(value)
        program.add_constraint(variables[index] == value)
    normal, _incidence = _orient(rays, _primitive(coordinates))
    return normal


def _normal_cone_through_ridge(rays, ridge, dimension: int):
    r"""The two facet normals whose incidence contains ``ridge``, by the 2-dimensional normal cone.

    The covectors annihilating the ridge span a 2-dimensional space `W`; the normal cone
    ``{w in W : <w,r_i> >= 0}`` is a pointed 2-cone whose extreme rays are exactly the two facet normals through
    the ridge. Every extreme ray of a pointed 2-cone is orthogonal to one of its boundary constraints, so the
    candidate set is complete and each candidate is tested exactly.
    """
    span_basis = matrix(ZZ, [list(ray) for ray in rays]).row_space().basis_matrix()
    if ridge:
        gram = span_basis * matrix(ZZ, [list(rays[index]) for index in sorted(ridge)]).transpose()
        kernel = gram.left_kernel_matrix()
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
    """The facet's rows in a basis of its own span, with the parent's positive functional carried along.

    Returns ``(local_rays, local_dimension)``.
    """
    face_rows = matrix(ZZ, [list(rays[index]) for index in sorted(incidence)])
    span_basis = face_rows.row_space().basis_matrix()
    local_dimension = span_basis.nrows()
    coordinates = span_basis * face_rows.transpose()
    local_rays = [
        tuple(int(coordinates[row, column]) for row in range(local_dimension))
        for column in range(coordinates.ncols())
    ]
    return local_rays, local_dimension


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

    if dimension == 2:
        normals = _normal_cone_through_ridge(rays, frozenset(), 2)
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
            for candidate in _normal_cone_through_ridge(rays, ridge, dimension):
                if _same_ray(candidate, normal):
                    continue
                oriented, incidence = _orient(rays, candidate)
                if incidence in facets:
                    continue
                if not is_facet_certificate(oriented, incidence, rays):
                    continue
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
            row = tuple(
                int(QQ(int(moved(label).numerator())) / QQ(int(moved(label).denominator()))) for label in labels
            )
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


def cell_facet_orbits(cell, stabilizer) -> tuple[FacetOrbitRecord, ...]:
    r"""One certified record per orbit of ``stabilizer`` on the facets of ``cone(cell.rows)``.

    ``stabilizer`` is the cell's stabilizer as isometries (``cell_stabilizer(cell)``); the record set carries a
    representative normal and incidence plus, for every member of the orbit, a word of generators transporting the
    representative onto it. The facet set itself is complete by ridge closure, so the orbit set is the whole orbit
    set -- not a sample.
    """
    rays = [tuple(int(entry) for entry in row) for row in cell.rows]
    height, pairings = pointedness_functional(cell)
    if not all(value > 0 for value in pairings):
        raise ArithmeticError("the dual slice requires a functional strictly positive on every ray")
    permutations = ray_permutations_from_stabilizer(cell, stabilizer)
    facets = cone_facets(rays)

    seen: set[frozenset[int]] = set()
    records = []
    for incidence in sorted(facets, key=lambda face: (len(face), sorted(face))):
        if incidence in seen:
            continue
        members = _orbit(incidence, permutations)
        if not set(members) <= set(facets):
            raise ArithmeticError("the stabilizer did not preserve the facet set, so the orbit grouping is invalid")
        seen.update(members)
        records.append(
            FacetOrbitRecord(
                normal=facets[incidence],
                incidence=incidence,
                orbit_size=len(members),
                words=tuple((member, word) for member, word in members.items()),
            )
        )
    if len(seen) != len(facets):
        raise ArithmeticError("the orbits did not cover the facet set")
    return tuple(records)
