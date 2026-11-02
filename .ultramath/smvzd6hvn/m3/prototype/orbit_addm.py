"""Orbit-wise adjacency-decomposition dual description (exact prototype).

Given a pointed full-dimensional rational cone ``C = cone(r_1, ..., r_m)`` in
``Q^n`` and a group ``G <= Sym(m)`` that permutes the rays (each generator is
induced by an integral isometry preserving the configuration, but only the
permutation is used here), enumerate the facets of ``C`` **modulo** ``G``: one
structural representative per orbit, with the exact Schreier word of every
orbit member.  Nothing in here calls a polyhedral engine.

The routine is the adjacency decomposition of the polyhedral_common /
Indefinite.jl reference (``DualDescriptionAndADM.g``::
``__ListFacetByAdjacencyDecompositionMethod``): a breadth-first search on
facet orbits, where the neighbours of a facet ``F`` are obtained from the
facets of the facet cone ``cone(F)`` (its ridges) through the two-facet pencil
of a ridge (``Fundamental.g``::``__ProjectionLiftingFramework_Rational``::
``FuncLift``).

Everything is exact integer / rational arithmetic.
"""

from __future__ import annotations

from collections import deque
from math import gcd

from sage.all import QQ, ZZ, Polyhedron, matrix, vector


# --------------------------------------------------------------------------
# small exact primitives
# --------------------------------------------------------------------------

def primitive(vec):
    """Divide an integer vector by its content; fix no sign."""
    content = 0
    for entry in vec:
        content = gcd(content, abs(int(entry)))
    if content == 0:
        raise ArithmeticError("the zero vector is not a ray of a pointed cone")
    return tuple(int(entry) // content for entry in vec)


def pair(covector, ray):
    return sum(int(a) * int(b) for a, b in zip(covector, ray, strict=True))


def same_ray(left, right):
    """True iff ``left = lambda * right`` for some rational ``lambda > 0``."""
    scale = None
    for a, b in zip(left, right, strict=True):
        if b == 0:
            if a != 0:
                return False
            continue
        if scale is None:
            candidate = QQ(a) / QQ(b)
            if candidate <= 0:
                return False
            scale = candidate
        elif QQ(a) / QQ(b) != scale:
            return False
    return scale is not None


def incidence(rows, covector):
    """The exact zero set of a covector on a ray list."""
    return frozenset(index for index, row in enumerate(rows) if pair(covector, row) == 0)


def rank_of(rows):
    if not rows:
        return 0
    return int(matrix(QQ, rows).rank())


def row_basis(rows):
    """An integer basis of the row space of ``rows`` (rows of maximal rank kept)."""
    if not rows:
        return []
    transposed = matrix(QQ, rows).transpose()
    pivots = tuple(int(p) for p in transposed.pivots())
    return [tuple(int(x) for x in rows[position]) for position in pivots]


def kernel_basis(rows, n):
    """Integer row basis of ``{v in Q^n : v . r = 0 for r in rows}``."""
    if not rows:
        return [tuple(1 if i == j else 0 for i in range(n)) for j in range(n)]
    basis = matrix(ZZ, rows).right_kernel_matrix()
    return [tuple(int(x) for x in row) for row in basis.rows()]


def facet_normal(rows, face, n):
    """The primitive integer supporting covector of a facet, from its incidence.

    Raises unless ``face`` has rank ``n - 1`` in the cone's ambient span and one
    of the two kernel directions is nonnegative on every ray (which is what
    "facet" means, and which is unique on a pointed cone).
    """
    basis = kernel_basis([rows[index] for index in sorted(face)], n)
    if len(basis) != 1:
        raise ArithmeticError(f"a facet incidence must have a one-dimensional kernel, got {len(basis)}")
    for sign in (1, -1):
        candidate = primitive(tuple(sign * entry for entry in basis[0]))
        values = [pair(candidate, row) for row in rows]
        if all(value >= 0 for value in values) and incidence(rows, candidate) == face:
            return candidate
    raise ArithmeticError("a facet incidence has no supporting normal")


# --------------------------------------------------------------------------
# the two-facet pencil of a ridge (adjacent facet operator)
# --------------------------------------------------------------------------

def facet_normals_through_ridge(rows, span_basis, ridge, n, *, ray_transpose=None):
    """The two primitive facet normals of ``cone(rows)`` whose incidence contains ``ridge``.

    ``ridge`` is a codimension-two face of ``cone(rows)``, so the covectors that
    vanish on it *and lie in the span of the configuration* form a two
    dimensional space ``W``; the normal cone ``N_ridge = {w in W : w . r_i >= 0}``
    is a pointed two dimensional cone whose extreme rays are precisely the two
    facet normals containing the ridge.  They are found by intersecting the
    half-planes ``{y : y . a_i >= 0}`` in a two dimensional coordinate system
    ``W = <p, q>``: every extreme ray of a pointed two dimensional cone is
    orthogonal to one of the boundary half-planes, so the candidate set
    ``{+/- (-a_i2, a_i1)}`` is complete, and each candidate is tested exactly.
    """
    ridge_rows = [rows[index] for index in sorted(ridge)]
    span = matrix(ZZ, span_basis)
    if ridge_rows:
        gram = span * matrix(ZZ, ridge_rows).transpose()
        kernel = gram.left_kernel_matrix()
        parameters = [tuple(int(x) for x in row) for row in kernel.rows()]
    else:
        parameters = [tuple(1 if i == j else 0 for j in range(span.nrows())) for i in range(span.nrows())]
    if len(parameters) != 2:
        raise ArithmeticError(
            f"a ridge of this cone must have a two-dimensional pencil, got {len(parameters)}"
        )
    basis = [tuple(int(x) for x in (vector(ZZ, theta) * span)) for theta in parameters]
    left, right = basis
    if ray_transpose is None:
        ray_transpose = matrix(ZZ, rows).transpose()
    projection_matrix = matrix(ZZ, 2, n, [*left, *right]) * ray_transpose
    projections = [
        (int(projection_matrix[0, index]), int(projection_matrix[1, index]))
        for index in range(len(rows))
        if projection_matrix[0, index] != 0 or projection_matrix[1, index] != 0
    ]
    if not projections:
        raise ArithmeticError("a ridge pencil must see at least one non-orthogonal ray")
    # The constraints span a pointed cone, hence an angular sector of width
    # below 180 degrees; on such a set "cross(b, a) < 0" is a total order, so
    # the angular minimum and maximum are one pass each.  With
    # theta_min/theta_max the extreme constraint angles, the normal cone is
    # the angular interval [theta_max - 90, theta_min + 90].
    angular_max = projections[0]
    for candidate in projections[1:]:
        if angular_max[0] * candidate[1] - angular_max[1] * candidate[0] > 0:
            angular_max = candidate
    angular_min = projections[0]
    for candidate in projections[1:]:
        if angular_min[0] * candidate[1] - angular_min[1] * candidate[0] < 0:
            angular_min = candidate
    candidates: list[tuple[int, ...]] = []
    for a, b in (
        (angular_max[1], -angular_max[0]),
        (-angular_min[1], angular_min[0]),
    ):
        if all(a * c + b * d >= 0 for c, d in projections):
            candidate = primitive(
                tuple(a * p + b * q for p, q in zip(left, right, strict=True))
            )
            if not any(same_ray(candidate, known) for known in candidates):
                candidates.append(candidate)
    if len(candidates) != 2:
        raise ArithmeticError(f"a ridge must lie in exactly two facets, found {len(candidates)}")
    return candidates


# --------------------------------------------------------------------------
# seed facet: one vertex of the dual slice
# --------------------------------------------------------------------------

def seed_facet(rows, height, n):
    """One facet incidence of a pointed cone, from a vertex of the dual slice.

    ``height`` is an integer covector that is strictly positive on every ray.
    The slice ``S = {u : u . r_i >= 0, height . u = 1}`` is a bounded polytope
    (``height`` is interior to the dual cone), and its lexicographic minimum is
    a vertex of ``S``.  Vertices of ``S`` are exactly the extreme rays of the
    dual cone, i.e. exactly the facet normals.  No genericity assumption is
    needed: the lexicographic order is total, so the optimum is a vertex, and it
    is obtained by ``n`` exact linear programs over the same feasible set (each
    one pinning the previous coordinate with an exact equality).
    """
    from sage.numerical.mip import MixedIntegerLinearProgram

    program = MixedIntegerLinearProgram(maximization=False, solver="PPL")
    variables = program.new_variable(real=True, indices=range(n))
    program.add_constraint(
        sum(ZZ(height[index]) * variables[index] for index in range(n)) == 1
    )
    for row in rows:
        program.add_constraint(
            sum(ZZ(row[index]) * variables[index] for index in range(n)) >= 0
        )
    # A cone that is not full dimensional has an unbounded dual: it contains
    # the annihilator of the cone's span.  Restrict to covectors inside the
    # span, which is where the intrinsic facet normals of the sub-cone live.
    for annihilator in kernel_basis(rows, n):
        program.add_constraint(
            sum(ZZ(annihilator[index]) * variables[index] for index in range(n)) == 0
        )
    coordinates = []
    for index in range(n):
        program.set_objective(variables[index])
        program.solve()
        value = QQ(program.get_values(variables)[index])
        coordinates.append(value)
        program.add_constraint(variables[index] == value)
    vertex = coordinates
    denominators = [entry.denominator() for entry in vertex]
    scale = 1
    for denominator in denominators:
        scale = scale * denominator // gcd(scale, denominator)
    covector = primitive(tuple(int(entry * scale) for entry in vertex))
    values = [pair(covector, row) for row in rows]
    if not all(value >= 0 for value in values):
        raise ArithmeticError("the dual-slice vertex is not a supporting covector")
    face = incidence(rows, covector)
    if rank_of([rows[index] for index in sorted(face)]) != rank_of(rows) - 1:
        raise ArithmeticError("the dual-slice vertex did not produce a facet")
    return face


# --------------------------------------------------------------------------
# the group action: orbit closure with exact Schreier words
# --------------------------------------------------------------------------

def permute_face(face, generator):
    return frozenset(generator[index] for index in face)


def orbit_closure(seed, generators):
    """Breadth-first orbit of ``seed`` with a word carrying the seed to each member."""
    members = {seed: ()}
    queue = deque([seed])
    while queue:
        face = queue.popleft()
        word = members[face]
        for index, generator in enumerate(generators):
            image = permute_face(face, generator)
            if image not in members:
                members[image] = word + (index,)
                queue.append(image)
    return members


def apply_word(face, generators, word):
    for index in word:
        face = permute_face(face, generators[index])
    return face


def setwise_stabilizer_generators(generators, face, point_count):
    """Generators of the setwise stabilizer of ``face``, restricted to ``face``.

    Schreier's lemma applied to the orbit of the set ``face``: with ``t_h`` a
    transversal element carrying ``face`` to the orbit member ``h``, the
    elements ``t_{s h}^{-1} . s . t_h`` (``s`` a generator) fix ``face`` and
    generate its setwise stabilizer.  This is exactly how the production code
    derives them from the transversal words that ``_iter_facet_mask_orbits``
    already retains; no permutation-group machinery is needed.
    """
    positions = tuple(sorted(face))
    if not generators:
        return [], positions
    members = orbit_closure(face, generators)
    elements = {}
    inverses = {}
    for member, word in members.items():
        permutation = list(range(point_count))
        for generator_index in word:
            generator = generators[generator_index]
            permutation = [generator[entry] for entry in permutation]
        elements[member] = permutation
        inverse = [0] * point_count
        for source, target in enumerate(permutation):
            inverse[target] = source
        inverses[member] = inverse
    local = {point: index for index, point in enumerate(positions)}
    restricted = []
    for member in members:
        transversal = elements[member]
        for generator in generators:
            image = permute_face(member, generator)
            inverse_image = inverses[image]
            permutation = []
            for point in positions:
                moved = generator[transversal[point]]
                permutation.append(local[inverse_image[moved]])
            candidate = tuple(permutation)
            if candidate != tuple(range(len(positions))) and candidate not in restricted:
                restricted.append(candidate)
    return restricted, positions


def restrict_generators(generators, face, point_count):
    """Backwards-compatible name for the setwise stabilizer routine above."""
    return setwise_stabilizer_generators(generators, face, point_count)


# --------------------------------------------------------------------------
# the adjacency decomposition
# --------------------------------------------------------------------------

class Stats:
    def __init__(self):
        self.nodes = 0
        self.seeds = 0
        self.inherited_seeds = 0
        self.ridges = 0
        self.pencil_calls = 0
        self.inserts = 0
        self.members = 0

    def __repr__(self):
        return (
            f"nodes={self.nodes} seeds={self.seeds} inherited={self.inherited_seeds} "
            f"ridges={self.ridges} pencils={self.pencil_calls} inserts={self.inserts} "
            f"members={self.members}"
        )


def facet_orbits(
    rows, generators, n, *, height, seed=None, seed_ridge=None, memo=None, stats=None,
    want_normals=True,
):
    """Facets of ``cone(rows)`` modulo ``<generators>``, one representative per orbit.

    Returns a list of records ``(face, normal, members)`` where ``face`` is the
    representative incidence (a frozenset of positions in ``rows``), ``normal``
    is its exact primitive supporting covector (``None`` when ``want_normals``
    is false, or in the rank-one base case) and ``members`` maps every orbit
    member to a word in ``generators`` carrying ``face`` to it.

    ``seed`` is one facet incidence of ``cone(rows)`` when the caller already
    knows one (an inherited seed); ``seed_ridge`` is a codimension-two face
    inside that seed facet, which seeds the first inner run.  Both are exact
    certificates, never guesses: they are used only as starting points, so a
    wrong one could only make the run fail its own rank/incidence checks, never
    silently change the result.

    Correctness needs: ``rows`` spans a ``d >= 1`` dimensional pointed cone,
    ``height . r > 0`` for every ray, and every generator permutes ``rows``.
    """
    if stats is None:
        stats = Stats()
    if memo is None:
        memo = {}
    rows = tuple(tuple(int(entry) for entry in row) for row in rows)
    generators = tuple(tuple(int(entry) for entry in generator) for generator in generators)
    key = (rows, generators, want_normals)
    if key in memo:
        return memo[key]
    dimension = rank_of(rows)
    if dimension == 0:
        raise ArithmeticError("a facet enumeration needs a nonzero cone")
    stats.nodes += 1
    if dimension == 1:
        result = [{"face": frozenset(), "normal": None, "members": {frozenset(): ()}}]
        memo[key] = result
        return result

    span_basis = row_basis(rows)
    ray_transpose = matrix(ZZ, rows).transpose()
    if seed is None:
        seed = seed_facet(rows, height, n)
        stats.seeds += 1
        inherited_ridge = None
    else:
        seed = frozenset(seed)
        if rank_of([rows[index] for index in sorted(seed)]) != dimension - 1:
            raise ArithmeticError("an inherited seed facet does not have codimension one")
        stats.inherited_seeds += 1
        inherited_ridge = None if seed_ridge is None else frozenset(seed_ridge)

    orbits: list[dict] = []
    where: dict[frozenset, int] = {}
    pending: deque[int] = deque()

    def insert(face, discovery_ridge=None):
        if face in where:
            return
        members = orbit_closure(face, generators)
        index = len(orbits)
        for member in members:
            where[member] = index
        stats.inserts += 1
        stats.members += len(members)
        orbits.append(
            {"face": face, "normal": None, "members": members, "ridge": discovery_ridge}
        )
        pending.append(index)

    insert(seed, inherited_ridge)
    while pending:
        current = pending.popleft()
        record = orbits[current]
        face = record["face"]
        local_generators, positions = restrict_generators(generators, face, len(rows))
        local_index = {point: index for index, point in enumerate(positions)}
        discovery_ridge = record["ridge"]
        inner_seed = None
        if discovery_ridge is not None:
            # The ridge through which this facet was found lies inside it, so it
            # is a facet of this facet cone: a free seed for the inner run.
            inner_seed = frozenset(local_index[point] for point in discovery_ridge)
        inner = facet_orbits(
            [rows[index] for index in positions],
            local_generators,
            n,
            height=height,
            seed=inner_seed,
            memo=memo,
            stats=stats,
            want_normals=False,
        )
        for ridge_record in inner:
            ridge = frozenset(positions[index] for index in ridge_record["face"])
            stats.ridges += 1
            stats.pencil_calls += 1
            for candidate in facet_normals_through_ridge(
                rows, span_basis, ridge, n, ray_transpose=ray_transpose
            ):
                neighbour = incidence(rows, candidate)
                if neighbour == face:
                    continue
                if rank_of([rows[index] for index in sorted(neighbour)]) != dimension - 1:
                    raise ArithmeticError("the ridge pencil produced a non-facet")
                insert(neighbour, ridge)
    if want_normals:
        for record in orbits:
            if record["normal"] is None:
                record["normal"] = facet_normal(rows, record["face"], n)
    memo[key] = orbits
    return orbits


# --------------------------------------------------------------------------
# validation helpers
# --------------------------------------------------------------------------

def full_facet_incidences(rows):
    """All facet incidences of ``cone(rows)`` from an exact engine, for cross-checking."""
    polyhedron = Polyhedron(rays=[list(row) for row in rows])
    facets = []
    for inequality in polyhedron.Hrepresentation():
        if not inequality.is_inequality():
            continue
        covector = tuple(int(x) for x in inequality.A())
        value = QQ(inequality.b())
        # b + a . x >= 0 with the cone (0 in the polyhedron): b >= 0, and b = 0
        # for a genuine facet of the cone.
        if value != 0:
            continue
        facets.append(incidence(rows, covector))
    return facets
