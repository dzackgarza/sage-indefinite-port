"""Standalone prototype of the reference's orbit-reduced facet enumeration for CONES.

Reconstructed from the vendored sources (see report.md):
  * Indefinite.jl@374a5ebb/indef/lib/DualDescriptionAndADM.g:77-188
      __ListFacetByAdjacencyDecompositionMethod: orbit work list over facet
      incidences, recursion on the incidence's own row subconfiguration with the
      set stabilizer, lift of each returned ridge to the adjacent facet.
  * Indefinite.jl@374a5ebb/indef/lib/Fundamental.g:605-690
      __ProjectionLiftingFramework_Rational: the lift is a 2-dimensional
      annihilator computation (NullspaceMat(TransposedMat(ridge rows))), an
      angular extreme-ray sweep in that 2-plane, and an incidence recovery by
      filtering rows that pair to zero.  No linear programme, no bound.
  * Indefinite.jl@374a5ebb/indef/lib/SetFunctionality.g:204-330
      OrbitGroupFormalism: insert / mark-done / minimal-undone work list.

This file is a *validation* artefact: it runs on synthetic symmetric cones where
the complete facet list can be obtained independently from Sage's polyhedron
backend, and asserts that
  (1) every enumerated incidence is a genuine facet (exact normal, nonnegative
      everywhere, rank n-1 incidence),
  (2) the orbit set equals the orbit set of the complete facet list,
  (3) the orbit set does not depend on the seed facet,
  (4) the cost is one subcone call per facet orbit, one lift per ridge orbit.

Deliberately NOT a port patch and NOT imported by the port.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from typing import Sequence

Q = Fraction


# ---------------------------------------------------------------- linear algebra
def rank(rows: Sequence[Sequence[int]]) -> int:
    return len(_echelon(rows)[1])


def _echelon(rows: Sequence[Sequence[int]]) -> tuple[list[list[Q]], list[int]]:
    """Row echelon form over Q; returns (nonzero rows, pivot columns)."""
    mat = [[Q(entry) for entry in row] for row in rows if any(row)]
    if not mat:
        return [], []
    width = len(mat[0])
    pivots: list[int] = []
    r = 0
    for c in range(width):
        pivot = next((i for i in range(r, len(mat)) if mat[i][c] != 0), None)
        if pivot is None:
            continue
        mat[r], mat[pivot] = mat[pivot], mat[r]
        lead = mat[r][c]
        mat[r] = [value / lead for value in mat[r]]
        for i in range(len(mat)):
            if i != r and mat[i][c] != 0:
                factor = mat[i][c]
                mat[i] = [a - factor * b for a, b in zip(mat[i], mat[r])]
        pivots.append(c)
        r += 1
        if r == len(mat):
            break
    return mat[:r], pivots


def pivot_columns(rows: Sequence[Sequence[int]]) -> list[int]:
    """Columns selecting an injective projection of span(rows) (reference __ProjectionFrame)."""
    return _echelon(rows)[1]


def project(rows: Sequence[Sequence[int]]) -> tuple[list[tuple[int, ...]], tuple[int, ...]]:
    cols = tuple(pivot_columns(rows))
    return [tuple(row[c] for c in cols) for row in rows], cols


def nullspace(rows: Sequence[Sequence[int]]) -> list[list[Q]]:
    """Basis of {a : a . r = 0 for every r in rows}, as rows of length = len(r)."""
    if not rows:
        raise ValueError("nullspace of an empty row family is not used here")
    width = len(rows[0])
    reduced, pivots = _echelon(rows)
    free = [c for c in range(width) if c not in pivots]
    basis: list[list[Q]] = []
    for column in free:
        vector = [Q(0)] * width
        vector[column] = Q(1)
        for index, pivot in enumerate(pivots):
            vector[pivot] = -reduced[index][column]
        basis.append(vector)
    return basis


def dot(left: Sequence[Q | int], right: Sequence[Q | int]) -> Q:
    return sum((Q(a) * Q(b) for a, b in zip(left, right)), Q(0))


# ---------------------------------------------------------------- 2-plane geometry
def _perp(point: Sequence[Q]) -> tuple[Q, Q]:
    return (Q(-point[1]), Q(point[0]))


def _primitive_direction(vector: Sequence[Q]) -> tuple[Q, Q]:
    """Canonical representative of the POSITIVE scaling class of a nonzero 2-vector.

    Dividing by the first nonzero entry keeps the sign: on a pointed cone v and -v
    are different rays, and in the dual 2-plane the halfplane y . p >= 0 is not
    invariant under y -> -y.
    """
    x, y = Q(vector[0]), Q(vector[1])
    if x != 0:
        scale = abs(x)          # positive scaling: the sign of the ray is data
        return (x / scale, y / scale)
    if y != 0:
        scale = abs(y)
        return (x / scale, y / scale)
    raise ValueError("zero vector has no direction")


def cone_extreme_rays_2d(points: Sequence[Sequence[Q]], want: int) -> list[tuple[Q, Q]]:
    """Extreme rays (primitive, deduplicated) of the cone spanned by nonzero 2-vectors."""
    rays: dict[tuple[Q, Q], None] = {}
    for p in points:
        key = _primitive_direction(p)
        if key in rays:
            continue
        # p is extreme iff the closed halfplane through p supporting the cone
        # contains every point: the supporting normal is perpendicular to p.
        normal = _perp(key)
        side = [Q(dot(normal, q)) for q in points]
        if all(s >= 0 for s in side) or all(s <= 0 for s in side):
            rays[key] = None
    if len(rays) != want:
        raise ArithmeticError(f"expected {want} extreme rays in the 2-plane, found {len(rays)}")
    return list(rays)


def dual_cone_extreme_rays_2d(points: Sequence[Sequence[Q]]) -> list[tuple[Q, Q]]:
    """Extreme rays of {y in Q^2 : y . p >= 0 for every p}.

    Every extreme ray is perpendicular to some input point (a 2-dim cone's
    extreme rays are intersections of two of its facet lines), so the candidates
    are +-perp(p); testing each exactly is O(m^2) and needs no bound.
    """
    found: dict[tuple[Q, Q], None] = {}
    for p in points:
        if p[0] == 0 and p[1] == 0:
            continue
        for candidate in (_perp(p), _perp((-p[0], -p[1]))):
            if all(dot(candidate, q) >= 0 for q in points):
                found[_primitive_direction(candidate)] = None
    if len(found) != 2:
        raise ArithmeticError(f"a genuine ridge supports exactly 2 facets, found {len(found)}")
    return list(found)


# ---------------------------------------------------------------- configuration
class Configuration:
    """An integer ray family with a permutation action on the rows.

    rows are projected to pivot columns so that rank(rows) == ambient dimension,
    exactly as __ProjectionFrame/ColumnReduction does in the reference.
    """

    def __init__(self, rows: Sequence[Sequence[int]], gens: Sequence[Sequence[int]]):
        projected, cols = project(rows)
        self.rows = projected
        self.columns = cols
        self.count = len(projected)
        self.dimension = len(projected[0]) if projected else 0
        if rank(self.rows) != self.dimension:
            raise ValueError("configuration must be projected to full rank")
        self.gens = tuple(tuple(g) for g in gens)

    # -- group closure: |G| is small in every specimen used here
    def group(self) -> tuple[tuple[int, ...], ...]:
        cached = self.__dict__.get("_group")
        if cached is not None:
            return cached
        identity = tuple(range(self.count))
        seen = {identity}
        frontier = [identity]
        while frontier:
            element = frontier.pop()
            for generator in self.gens:
                image = tuple(generator[element[i]] for i in range(self.count))
                if image not in seen:
                    seen.add(image)
                    frontier.append(image)
        self.__dict__["_group"] = tuple(seen)
        return self.__dict__["_group"]

    def pointed(self) -> bool:
        """True iff the cone has no lineality: some strictly positive form exists."""
        from sage.all import QQ, Polyhedron, vector

        cone = Polyhedron(rays=[vector(QQ, row) for row in self.rows])
        return cone.codim() == 0 and len(cone.lines()) == 0

    def stabilizer(self, incidence: Sequence[int]) -> tuple[tuple[int, ...], ...]:
        target = tuple(sorted(incidence))
        return tuple(
            element
            for element in self.group()
            if tuple(sorted(element[i] for i in target)) == target
        )


# ---------------------------------------------------------------- the algorithm
class Stats:
    def __init__(self) -> None:
        self.subcone_calls = 0
        self.lift_calls = 0
        self.facet_orbits = 0
        self.max_incidence = 0
        self.expansions = 0
        self.heights: list[int] = []
        self.cache_hits = 0
        self.level_lifts: dict[int, int] = {}
        self.level_expansions: dict[int, int] = {}
        self.balinski_hits = 0
        self.shared_hits = 0


def normal_and_incidence(rows: Sequence[tuple[int, ...]], incidence: Sequence[int]) -> tuple[Q, ...]:
    """The primitive-ish supporting form of an incidence, oriented positive on rows."""
    space = nullspace([rows[i] for i in incidence])
    if len(space) != 1:
        raise ArithmeticError("a facet incidence must have a one-dimensional annihilator")
    normal = space[0]
    values = [dot(normal, row) for row in rows]
    if any(v < 0 for v in values):
        normal = [-v for v in normal]
    return tuple(normal)


def is_facet(rows: Sequence[tuple[int, ...]], incidence: Sequence[int]) -> bool:
    if not incidence or len(incidence) == len(rows):
        return False
    if rank([rows[i] for i in incidence]) != len(rows[0]) - 1:
        return False
    normal = normal_and_incidence(rows, incidence)
    values = [dot(normal, row) for row in rows]
    zeros = tuple(i for i, v in enumerate(values) if v == 0)
    return all(v >= 0 for v in values) and zeros == tuple(sorted(incidence))


def seed_facet(rows: Sequence[tuple[int, ...]], rng, stats: Stats) -> tuple[int, ...]:
    """One facet, found by exact candidate-and-verify.

    Any facet is an acceptable seed: completeness comes from orbit saturation,
    never from the seed.  Selection here is random, verification is exact.
    """
    width = len(rows[0])
    indices = list(range(len(rows)))
    for _attempt in range(4000):
        rng.shuffle(indices)
        basis = indices[: width - 1]
        if rank([rows[i] for i in basis]) != width - 1:
            continue
        space = nullspace([rows[i] for i in basis])
        normal = space[0]
        values = [dot(normal, row) for row in rows]
        if any(v < 0 for v in values):
            normal = [-v for v in normal]
            values = [-v for v in values]
        if any(v < 0 for v in values):
            continue
        incidence = tuple(i for i, v in enumerate(values) if v == 0)
        if rank([rows[i] for i in incidence]) == width - 1:
            assert is_facet(rows, incidence)
            return incidence
    raise ArithmeticError("no facet found by candidate-and-verify; is the cone pointed?")


def lift(rows: Sequence[tuple[int, ...]], facet: Sequence[int], ridge: Sequence[int]) -> tuple[int, ...]:
    """The unique facet != facet containing the ridge (Fundamental.g:614-688)."""
    ridge_rows = [rows[i] for i in ridge]
    space = nullspace(ridge_rows)
    if len(space) != 2:
        raise ArithmeticError(f"a ridge must be codimension 2, got annihilator {len(space)}")
    points = [(dot(space[0], row), dot(space[1], row)) for row in rows]
    candidates = dual_cone_extreme_rays_2d([p for p in points if p != (Q(0), Q(0))])
    normals = []
    for y in candidates:
        normals.append([y[0] * space[0][i] + y[1] * space[1][i] for i in range(len(space[0]))])
    facet_set = set(facet)
    keep = []
    for normal in normals:
        zeros = {i for i, row in enumerate(rows) if dot(normal, row) == 0}
        if zeros >= facet_set:
            continue  # this is the facet we came from
        keep.append(normal)
    if len(keep) != 1:
        raise ArithmeticError(f"expected exactly one adjacent facet, got {len(keep)}")
    normal = keep[0]
    values = [dot(normal, row) for row in rows]
    if any(v < 0 for v in values):
        raise ArithmeticError("the lifted form is not a supporting form")
    adjacent = tuple(i for i, v in enumerate(values) if v == 0)
    if not is_facet(rows, adjacent):
        raise ArithmeticError("the lift did not produce a facet")
    return adjacent


_RIDGE_CACHE: dict = {}


def facet_orbits(
    configuration: Configuration,
    *,
    stats: Stats | None = None,
    rng=None,
    depth: int = 0,
) -> list[tuple[int, ...]]:
    """One facet incidence per orbit of cone(rows) under configuration.gens.

    Memoized on (row multiset, generator tuple): the answer is a function of the
    subconfiguration and its group alone.  This is the reference's
    BankRecording (DualDescriptionAndADM.g:12-74, 79-82); it is what keeps the
    recursion a DAG instead of a tree.
    """
    key = (tuple(sorted(configuration.rows)), tuple(sorted(configuration.gens)))
    cached = _RIDGE_CACHE.get(key)
    if cached is not None:
        stats_ = stats if stats is not None else Stats()
        stats_.cache_hits += 1
        return list(cached)
    stats = stats if stats is not None else Stats()
    stats.subcone_calls += 1
    stats.heights.append(depth)
    rows = configuration.rows
    dimension = configuration.dimension
    group = configuration.group()
    if dimension <= 2:
        if dimension == 1:
            return []
        points = [(Q(row[0]), Q(row[1])) for row in rows]
        rays = cone_extreme_rays_2d([p for p in points if p != (Q(0), Q(0))], 2)
        facets = []
        for ray in rays:
            incidence = tuple(i for i, p in enumerate(points) if p == (Q(0), Q(0)) or _primitive_direction(p) == ray)
            facets.append(incidence)
        return _reduce_orbits(facets, group)
    seed = seed_facet(rows, rng, stats)
    result = _orbit_representatives([seed], group, rows, stats, configuration=configuration, rng=rng, depth=depth)
    _RIDGE_CACHE[(tuple(sorted(rows)), tuple(sorted(configuration.gens)))] = tuple(result)
    return result


def _reduce_orbits(seeds, group):
    """One representative per orbit of the given incidences (base case)."""
    seen = set()
    representatives = []
    for seed in seeds:
        target = tuple(sorted(seed))
        if target in seen:
            continue
        orbit = {tuple(sorted(element[i] for i in target)) for element in group}
        seen |= orbit
        representatives.append(min(orbit))
    return representatives


def _orbit_representatives(
    seeds: Sequence[tuple[int, ...]],
    group: Sequence[tuple[int, ...]],
    rows: Sequence[tuple[int, ...]],
    stats: Stats,
    *,
    configuration: Configuration | None = None,
    rng=None,
    depth: int = 0,
) -> list[tuple[int, ...]]:
    if configuration is None:
        # base case: the seeds already exhaust the facets
        return list(seeds)

    def canonical(incidence: Sequence[int]) -> tuple[int, ...]:
        target = tuple(sorted(incidence))
        return min(tuple(sorted(element[i] for i in target)) for element in group)

    registry: dict[tuple[int, ...], tuple[int, ...]] = {}
    done: set[tuple[int, ...]] = set()
    order: list[tuple[int, ...]] = []
    for seed in seeds:
        key = canonical(seed)
        if key not in registry:
            registry[key] = key
            order.append(key)
    while True:
        pending = [key for key in order if key not in done]
        if not pending:
            break
        if CUTOFF:
            # DualDescriptionAndADM.g:131-142, verbatim: quit expanding when the
            # number of *undone facet elements* drops to TheDim-1, or when all
            # undone facet orbits share a group orbit of rows.
            size = len(group)
            undone = sum(size // _stabilizer_order(key, group) for key in pending)
            if done and undone <= (len(rows[0]) - 1) - 1:
                stats.balinski_hits += 1
                break
            shared = _shared_orbit_rows(pending, group) if RULE2 else set()
            if done and shared:
                stats.shared_hits += 1
                break
        representative = min(pending, key=lambda key: (len(key), key))
        done.add(representative)
        stats.expansions += 1
        stats.max_incidence = max(stats.max_incidence, len(representative))
        # ridge orbits of the representative: facets of the facet cone,
        # modulo the set stabilizer of the representative.
        stabilizer = tuple(
            element
            for element in group
            if tuple(sorted(element[i] for i in representative)) == representative
        )
        sub_rows = [rows[i] for i in representative]
        sub_configuration = Configuration(
            sub_rows,
            [_restrict(element, representative) for element in stabilizer],
        )
        ridges = facet_orbits(sub_configuration, stats=stats, rng=rng, depth=depth + 1)
        stats.level_expansions[depth] = stats.level_expansions.get(depth, 0) + 1
        for ridge in ridges:
            stats.lift_calls += 1
            stats.level_lifts[depth] = stats.level_lifts.get(depth, 0) + 1
            adjacent_local = lift(rows, representative, [representative[i] for i in ridge])
            adjacent = tuple(sorted(adjacent_local))
            assert is_facet(rows, adjacent), "enumerated incidence is not a facet"
            key = canonical(adjacent)
            if key not in registry:
                registry[key] = key
                order.append(key)
    stats.facet_orbits += len(registry)
    return list(registry.values())


CUTOFF = False
RULE2 = True


def _stabilizer_order(incidence, group):
    target = tuple(sorted(incidence))
    return sum(1 for element in group if tuple(sorted(element[i] for i in target)) == target)


def _shared_orbit_rows(pending, group):
    """Rows lying in a full group orbit inside every pending incidence (OrbitIntersection)."""
    if not pending:
        return set()
    shared = None
    for key in pending:
        fixed = set(key)
        for element in group:
            fixed &= {element[i] for i in key}
        shared = fixed if shared is None else shared & fixed
        if not shared:
            return set()
    return shared


def _restrict(element: Sequence[int], incidence: Sequence[int]) -> tuple[int, ...]:
    """Permutation induced by `element` on the positions inside `incidence`."""
    positions = {value: index for index, value in enumerate(incidence)}
    return tuple(positions[element[i]] for i in incidence)


# ---------------------------------------------------------------- oracle and checks
def oracle_facets(rows: Sequence[Sequence[int]]) -> list[tuple[int, ...]]:
    """The complete facet list, from the maintained polyhedron backend."""
    from sage.all import QQ, Polyhedron, vector

    cone = Polyhedron(rays=[vector(QQ, row) for row in rows])
    if cone.dimension() != len(rows[0]):
        raise AssertionError("the oracle specimen must be full-dimensional")
    facets = []
    for facet in cone.Hrepresentation():
        if not facet.is_inequality():
            continue
        normal = [QQ(entry) for entry in facet.A()]
        values = [sum(a * b for a, b in zip(normal, row)) for row in rows]
        if any(v < 0 for v in values):
            normal = [-v for v in normal]
            values = [-v for v in values]
        if any(v < 0 for v in values):
            raise AssertionError("backend facet is not a supporting form")
        facets.append(tuple(i for i, v in enumerate(values) if v == 0))
    return facets


def orbit_partition(facets, group):
    seen = set()
    orbits = []
    for facet in facets:
        target = tuple(sorted(facet))
        if target in seen:
            continue
        orbit = {tuple(sorted(element[i] for i in target)) for element in group}
        seen |= orbit
        orbits.append(min(orbit))
    return orbits


def check_specimen(name: str, rows, gens, *, verbose: bool = True, cutoff: bool = False):
    global CUTOFF
    CUTOFF = cutoff
    import random as _random

    configuration = Configuration(rows, gens)
    group = configuration.group()
    reference = oracle_facets(configuration.rows)
    expected = orbit_partition(reference, group)
    stats = Stats()
    produced = facet_orbits(configuration, stats=stats, rng=_random.Random(20240607))
    produced_keys = {tuple(sorted(facet)) for facet in produced}
    expected_keys = {tuple(sorted(facet)) for facet in expected}
    assert produced_keys == expected_keys, (
        f"{name}: orbit sets differ; missing={sorted(expected_keys - produced_keys)[:3]} "
        f"extra={sorted(produced_keys - expected_keys)[:3]}"
    )
    # seed independence: two other seeds must give the same orbit set again
    for trial in range(2):
        other = Stats()
        again = facet_orbits(configuration, stats=other, rng=_random.Random(1000 + trial))
        assert {tuple(sorted(f)) for f in again} == expected_keys, f"{name}: seed-dependent result"
    if verbose:
        print(
            f"{name}: rays={configuration.count} dim={configuration.dimension} "
            f"|G|={len(group)} facets={len(reference)} facet_orbits={len(expected)} "
            f"ray_orbits={len({min(tuple(sorted(e[i] for i in range(configuration.count))) for e in group) for _ in [0]}) } "
            f"expansions={stats.expansions} subcone_calls={stats.subcone_calls} "
            f"lifts={stats.lift_calls} max_incidence={stats.max_incidence} "
            f"depth={max(stats.heights)} cache_hits={stats.cache_hits} "
            f"L0:{stats.level_expansions.get(0, 0)}exp/{stats.level_lifts.get(0, 0)}lifts "
            f"L1:{stats.level_expansions.get(1, 0)}exp/{stats.level_lifts.get(1, 0)}lifts "
            f"balinski_stops={stats.balinski_hits} shared_orbit_stops={stats.shared_hits}"
        )
    return stats


# ---------------------------------------------------------------- specimens
def _apply(row, matrix):
    return tuple(sum(row[k] * matrix[k][j] for k in range(len(row))) for j in range(len(row)))


def orbit_cone(seed_rays, matrices, name=""):
    """Cone generated by the union of the matrix-group orbits of the given seeds.

    Pointed and full-dimensional by construction of the caller, with the induced
    permutation action on the rows -- the exact data an orbit-reduced facet
    enumeration receives.
    """
    rows: list[tuple[int, ...]] = []
    seen: set[tuple[int, ...]] = set()
    frontier = [tuple(int(entry) for entry in ray) for ray in seed_rays]
    while frontier:
        ray = frontier.pop()
        if ray in seen or not any(ray):
            continue
        seen.add(ray)
        rows.append(ray)
        for matrix in matrices:
            image = _apply(ray, matrix)
            if image not in seen:
                frontier.append(image)
    rows.sort()
    position = {row: index for index, row in enumerate(rows)}
    gens = [tuple(int(position[_apply(row, matrix)]) for row in rows) for matrix in matrices]
    return rows, gens


def identity(dimension: int):
    return [[1 if i == j else 0 for j in range(dimension)] for i in range(dimension)]


def transposition(dimension: int, i: int, j: int):
    matrix = identity(dimension)
    matrix[i][i] = matrix[j][j] = 0
    matrix[i][j] = matrix[j][i] = 1
    return matrix


def cycle_matrix(dimension: int):
    return [[1 if j == (i + 1) % dimension else 0 for j in range(dimension)] for i in range(dimension)]


def sign_flip(dimension: int, index: int):
    matrix = identity(dimension)
    matrix[index][index] = -1
    return matrix


def simplex_cone(rank_):
    """Cone over the standard `rank_`-simplex: `rank_+1` rays, S_{rank_+1} action."""
    rows = []
    for index in range(rank_ + 1):
        row = [0] * (rank_ + 1)
        row[index] = 1
        rows.append(tuple(row))
    dimension = rank_ + 1

    def swap(i, j, size=dimension):
        return tuple(j if k == i else i if k == j else k for k in range(size))

    return rows, [swap(0, i) for i in range(1, dimension)]


def specimens():
    spec = [("simplex4", simplex_cone(4)), ("simplex5", simplex_cone(5))]
    # cones in the positive orthant: pointed, non-simplicial facets, nontrivial stabilizers
    s4 = [transposition(4, 0, 1), cycle_matrix(4)]
    for seeds in (
        [(1, 1, 0, 0), (1, 1, 1, 0)],
        [(1, 2, 0, 0), (1, 1, 1, 1)],
        [(2, 1, 1, 0), (0, 1, 1, 1)],
        [(1, 0, 1, 0), (1, 2, 2, 0)],
    ):
        spec.append((f"S4{seeds}", orbit_cone(seeds, s4)))
    b5 = [transposition(5, 0, 1), transposition(5, 1, 2), transposition(5, 3, 4)]
    spec.append(("S3xS2", orbit_cone([(1, 1, 0, 1, 1), (1, 0, 1, 1, 0)], b5)))
    c6 = [cycle_matrix(6)]
    spec.append(("cyclic6", orbit_cone([(1, 2, 0, 0, 0, 0), (2, 1, 1, 0, 0, 0)], c6)))
    spec.append(("cyclic6b", orbit_cone([(1, 1, 1, 0, 0, 0), (1, 1, 0, 1, 0, 0)], c6)))
    spec.append(("cyclic6c", orbit_cone([(1, 1, 1, 1, 0, 0), (2, 1, 1, 1, 0, 0)], c6)))
    # a Lorentzian-flavoured cone: rays on the light cone's positive shell in signature (3,1),
    # plus the spacelike coordinate directions, with the S3 action on the spacelike part
    lor = [
        (1, 1, 0, 1),
        (1, 0, 1, 1),
        (0, 1, 1, 1),
        (1, 1, 1, 2),
        (2, 1, 1, 2),
        (1, 0, 0, 1),
        (0, 1, 0, 1),
        (0, 0, 1, 1),
    ]
    spec.append(("lorentzian_S3", orbit_cone(lor, [transposition(4, 0, 1), transposition(4, 1, 2)])))
    return spec


def main() -> None:
    import random as _random

    import sys as _sys

    global CUTOFF
    with_cutoff = "--cutoff" in _sys.argv or "--balinski" in _sys.argv
    global RULE2
    RULE2 = "--balinski" not in _sys.argv
    print("== upstream cutoff emulation (nbUndone <= TheDim-1 or shared undone row orbit)" if with_cutoff else "== full saturation")
    incomplete = []
    for name, data in specimens():
        rows, gens = data
        configuration = Configuration(rows, gens)
        if not configuration.pointed():
            print(f"{name}: skipped (not pointed)")
            continue
        try:
            check_specimen(name, rows, gens, cutoff=with_cutoff)
        except AssertionError as error:
            incomplete.append((name, str(error)[:120]))
            print(f"{name}: INCOMPLETE UNDER CUTOFF: {error}")
    if incomplete:
        print("cutoff made the enumeration incomplete on:", [name for name, _ in incomplete])


if __name__ == "__main__":
    main()
