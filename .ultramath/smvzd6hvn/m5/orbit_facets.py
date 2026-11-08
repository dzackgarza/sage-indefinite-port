"""Exact orbit-reduced facet enumeration by adjacency decomposition (standalone prototype).

Self-contained reference implementation of the construction specified in SOLUTION.md.  No
dependency on the port, SageMath, Normaliz or any LP solver: integer and `fractions.Fraction`
arithmetic throughout, with sympy only for rank and nullspace of the small integer matrices that
occur (at most (n-2) x n).

Objects
-------
A *level* is a ray configuration `rows`: m nonzero integer tuples of length d, of rank d, whose
cone is pointed, together with a *positive covector* `g` -- a rational covector with g.r > 0 for
every ray r.  For a perfect cell the top level always has one: the grading covector, which pairs
to a positive constant with every configuration vector.

A *group action* is `perms`: position permutations generating a subgroup G of Sym(m) that carries
the ray multiset to itself.

A *facet incidence* is the sorted tuple of positions on which a facet normal vanishes; a *facet
normal* is a primitive integer covector a with <a,r> >= 0 for all rays whose zero set spans a
hyperplane.  Every normal is certified before it leaves a routine.

    one_facet_by_descent(rows, g)        -- one facet normal, from the positive covector g
    pivot_across_ridge(rows, a, ridge)   -- the unique other facet normal through a ridge
    facet_orbits(rows, perms, g)         -- one facet normal per stabilizer orbit, recursively
"""
from __future__ import annotations

from collections import deque
from fractions import Fraction
from math import gcd

from sympy import Matrix

from two_dim import extremes_2d


# ---------------------------------------------------------------------------------------------
# exact helpers
# ---------------------------------------------------------------------------------------------
def dot(u, v):
    return sum(x * y for x, y in zip(u, v))


def primitive(v) -> tuple[int, ...]:
    """Scale a nonzero rational/integer vector to primitive integer form (sign preserved)."""
    entries = [Fraction(x) for x in v]
    scale = 1
    for entry in entries:
        scale = scale * entry.denominator // gcd(scale, entry.denominator)
    integers = [int(entry * scale) for entry in entries]
    content = 0
    for entry in integers:
        content = gcd(content, abs(entry))
    if content == 0:
        raise ArithmeticError("the zero covector is not a ray normal")
    return tuple(entry // content for entry in integers)


def rank_int(rows) -> int:
    if not rows:
        return 0
    return int(Matrix([list(row) for row in rows]).rank())


def nullspace_int(rows, ncols: int):
    """An integer basis of {x in Q^ncols : <row, x> = 0 for every row}."""
    matrix = Matrix([list(row) for row in rows]) if rows else Matrix.zeros(0, ncols)
    return [primitive(vector) for vector in matrix.nullspace()]


def pivot_columns(rows):
    """Column indices whose projection is injective on the row span (an invertible block)."""
    matrix = Matrix([list(row) for row in rows])
    return tuple(int(column) for column in matrix.rref()[1])


def independent_subset(rows, size: int):
    """A greedily selected independent subset of the given row vectors, of the requested size."""
    chosen: list[int] = []
    span_rank = 0
    for index in range(len(rows)):
        new_rank = rank_int([rows[i] for i in chosen + [index]])
        if new_rank > span_rank:
            chosen = chosen + [index]
            span_rank = new_rank
        if span_rank == size:
            return tuple(chosen)
    raise ArithmeticError("the configuration does not contain an independent subset of that size")


def solve_rational(square, target):
    """Solve square * x = target exactly over the rationals."""
    solution = Matrix([list(row) for row in square]).solve(Matrix([list(target)]).T)
    return tuple(Fraction(value) for value in solution)


# ---------------------------------------------------------------------------------------------
# facet bookkeeping
# ---------------------------------------------------------------------------------------------
def incidence_of(rows, a) -> tuple[int, ...]:
    """The zero set of an exact covector, as sorted ray positions."""
    return tuple(index for index, row in enumerate(rows) if dot(a, row) == 0)


def normal_of_incidence(rows, incidence):
    """The unique primitive covector vanishing on `incidence`, oriented nonnegative on all rays."""
    candidates = nullspace_int([rows[i] for i in incidence], len(rows[0]))
    if len(candidates) != 1:
        raise ArithmeticError(
            f"a facet incidence must annihilate a hyperplane, got nullspace dimension {len(candidates)}"
        )
    a = candidates[0]
    if any(dot(a, row) < 0 for row in rows):
        a = tuple(-x for x in a)
    return a


def certify_facet(rows, a) -> tuple[int, ...]:
    """Check the facet contract exactly and return the incidence."""
    d = len(rows[0])
    if not any(a):
        raise ArithmeticError("the zero covector defines no facet")
    values = [dot(a, row) for row in rows]
    if any(value < 0 for value in values):
        raise ArithmeticError(f"candidate {a} is not nonnegative on the ray configuration")
    incidence = tuple(index for index, value in enumerate(values) if value == 0)
    zero_rank = rank_int([rows[i] for i in incidence])
    if zero_rank != d - 1:
        raise ArithmeticError(f"candidate {a} has a zero set of rank {zero_rank}, not {d - 1}")
    return incidence


# ---------------------------------------------------------------------------------------------
# the pivot: the unique other facet normal through a ridge
# ---------------------------------------------------------------------------------------------
def pivot_across_ridge(rows, a, ridge):
    """The other facet normal of cone(rows) whose incidence contains `ridge`.

    The set W = ann(span ridge) is two dimensional and contains the current facet normal `a`.
    Inside W the supporting covectors form the normal cone of the ridge, a pointed two-dimensional
    cone whose two extreme rays are the two facets through the ridge.  The computation reduces W
    to coordinates (a, e) and takes the polar of the two extreme rays of the primal image cone.
    """
    d = len(rows[0])
    ridge_rows = [rows[i] for i in ridge]
    if rank_int(ridge_rows) != d - 2:
        raise ArithmeticError(f"a ridge must span rank {d - 2}, got {rank_int(ridge_rows)}")
    annihilator = nullspace_int(ridge_rows, d)
    if len(annihilator) != 2:
        raise ArithmeticError("the ridge annihilator must be two dimensional")
    a = primitive(a)
    if any(dot(a, row) != 0 for row in ridge_rows):
        raise ArithmeticError("the current facet normal does not vanish on the ridge")
    if any(dot(a, row) < 0 for row in rows):
        raise ArithmeticError("the current normal is not a supporting covector")
    e = next((vector for vector in annihilator if rank_int([a, vector]) == 2), None)
    if e is None:
        raise ArithmeticError("the ridge annihilator is spanned by the facet normal itself")

    # Images of the rays in the two-dimensional coordinates (a, e).
    images = [(dot(a, row), dot(e, row)) for row in rows]
    slopes = [(Fraction(second, first), index) for index, (first, second) in enumerate(images) if first > 0]
    if not slopes:
        raise ArithmeticError("the facet incidence cannot cover every ray")
    candidates = [(1, min(slopes)[0]), (1, max(slopes)[0])]
    if any(first == 0 and second > 0 for first, second in images):
        candidates.append((0, 1))
    if any(first == 0 and second < 0 for first, second in images):
        candidates.append((0, -1))
    primal_extremes = extremes_2d(candidates)

    def polar_normal(extreme):
        normal = (extreme[1], -extreme[0])
        if any(normal[0] * first + normal[1] * second < 0 for first, second in images):
            normal = (-normal[0], -normal[1])
        return normal

    normals = []
    for extreme in primal_extremes:
        first, second = polar_normal(extreme)
        candidate = primitive([Fraction(first) * Fraction(a[k]) + Fraction(second) * Fraction(e[k])
                               for k in range(d)])
        if all(dot(candidate, row) >= 0 for row in rows):
            normals.append(candidate)
    distinct = [candidate for candidate in normals if rank_int([a, candidate]) == 2]
    if len(distinct) != 1:
        raise ArithmeticError(
            "the prescribed set is not a ridge of the current facet: its normal cone does not "
            f"have the current facet as one of two extremes (found {len(distinct)} other extremes)"
        )
    b = distinct[0]
    if any(dot(b, rows[index]) != 0 for index in ridge):
        raise ArithmeticError("the pivoted normal does not vanish on the ridge")
    certify_facet(rows, b)
    return b


# ---------------------------------------------------------------------------------------------
# facet finding: iterated two-dimensional reduction from a positive covector
# ---------------------------------------------------------------------------------------------
def one_facet_by_descent(rows, positive_covector, *, stats=None):
    """A facet normal of cone(rows), by strictly increasing the rank of a zero set.

    Each step picks `e` in the annihilator of the current zero set outside span(a) and moves to
    the min-slope extreme of the two-dimensional image cone, which zeroes at least one further
    ray outside the current span.  The zero-set rank therefore grows by at least one per step, so
    at most d-1 steps reach a hyperplane, i.e. a facet.
    """
    d = len(rows[0])
    a = primitive(positive_covector)
    if not any(a) or any(dot(a, row) < 0 for row in rows):
        raise ArithmeticError("the descent requires a nonzero covector nonnegative on every ray")
    steps = 0
    while True:
        incidence = incidence_of(rows, a)
        if rank_int([rows[i] for i in incidence]) == d - 1:
            if stats is not None:
                stats["descent_steps"] = stats.get("descent_steps", 0) + steps
            return a
        annihilator = nullspace_int([rows[i] for i in incidence], d)
        e = next((vector for vector in annihilator if rank_int([a, vector]) == 2), None)
        if e is None:
            raise ArithmeticError("the descent cannot leave a zero set that is already maximal")
        slopes = [
            (Fraction(dot(e, row), dot(a, row)), index)
            for index, row in enumerate(rows)
            if dot(a, row) > 0
        ]
        if not slopes:
            raise ArithmeticError("the descent found no ray with a strictly positive pairing")
        slope, witness = min(slopes)
        a = primitive([Fraction(e[k]) - slope * Fraction(a[k]) for k in range(d)])
        steps += 1
        if steps > d + 2:
            raise ArithmeticError("the facet descent failed to terminate")
        if dot(a, rows[witness]) != 0:
            raise ArithmeticError("the facet descent did not enlarge the zero set")


# ---------------------------------------------------------------------------------------------
# permutation-group helpers (position permutations, exact and explicit)
# ---------------------------------------------------------------------------------------------
def compose(left, right):
    """The product g*h as position maps, applying `right` first."""
    return tuple(left[right[i]] for i in range(len(right)))


def invert(permutation):
    inverse = [0] * len(permutation)
    for position, image in enumerate(permutation):
        inverse[image] = position
    return tuple(inverse)


def word_as_permutation(permutations, word, degree=None):
    if degree is None:
        if not permutations:
            raise ArithmeticError("a word needs an explicit degree when the generator list is empty")
        degree = len(permutations[0])
    permutation = tuple(range(degree))
    for generator_index in word:
        permutation = compose(permutations[generator_index], permutation)
    return permutation


def orbit_closure(incidence, permutations):
    """Every group image of `incidence`, each with a word carrying `incidence` onto it."""
    start = tuple(sorted(incidence))
    closure = {start: ()}
    queue = deque([start])
    steps = 0
    while queue:
        current = queue.popleft()
        word = closure[current]
        for generator_index, permutation in enumerate(permutations):
            image = tuple(sorted(permutation[position] for position in current))
            steps += 1
            if image not in closure:
                closure[image] = word + (generator_index,)
                queue.append(image)
    return closure, steps


def stabilizer_generators(incidence, permutations, degree):
    """Generators of the set stabilizer of `incidence`, by Schreier's lemma, as permutations.

    For each orbit member x with a transporter t_x (t_x(incidence) = x) and each generator s
    fixing x, the conjugate t_x^-1 s t_x stabilizes the incidence; those conjugates generate the
    whole stabilizer (orbit-stabilizer with Schreier's lemma).
    """
    if not permutations:
        return ()
    closure, _ = orbit_closure(incidence, permutations)
    identity = tuple(range(degree))
    generators = []
    seen = set()
    for image, word in closure.items():
        transporter = word_as_permutation(permutations, word, degree)
        transporter_inverse = invert(transporter)
        for permutation in permutations:
            if tuple(sorted(permutation[position] for position in image)) != image:
                continue
            element = compose(compose(transporter_inverse, permutation), transporter)
            if element in seen or element == identity:
                continue
            seen.add(element)
            generators.append(element)
    return tuple(generators)


def restrict_generators(permutations, subset):
    """The action of a set stabilizer on positions of `subset`, re-indexed 0..len(subset)-1."""
    positions = {position: index for index, position in enumerate(subset)}
    restricted = []
    for permutation in permutations:
        image = tuple(permutation[position] for position in subset)
        if any(position not in positions for position in image):
            raise ArithmeticError("a purported set stabilizer leaves the set")
        restricted.append(tuple(positions[position] for position in image))
    return tuple(restricted)


# ---------------------------------------------------------------------------------------------
# the recursion: a subcone is described in the coordinates of a projection frame
# ---------------------------------------------------------------------------------------------
def subcone_level(rows, incidence, normal, positive_covector):
    """Rows and positive covector for the facet cone of the given incidence.

    The facet cone lives in normal^perp, of dimension d-1.  A set of coordinate directions on
    which that hyperplane projects injectively keeps the data integral; the positive covector is
    transported by solving the (d-1) x (d-1) system on a basis of the hyperplane.
    """
    d = len(rows[0])
    hyperplane_basis = nullspace_int([normal], d)             # basis of normal^perp, dim d-1
    if len(hyperplane_basis) != d - 1:
        raise ArithmeticError("the facet normal must be nonzero")
    frame = pivot_columns(hyperplane_basis)
    if len(frame) != d - 1:
        raise ArithmeticError("the hyperplane projection frame must have d-1 coordinates")
    coordinates = tuple(tuple(rows[index][column] for column in frame) for index in incidence)
    if rank_int(coordinates) != d - 1:
        raise ArithmeticError("the facet cone must be full dimensional in its hyperplane")
    projected_basis = [tuple(vector[column] for column in frame) for vector in hyperplane_basis]
    target = [Fraction(dot(positive_covector, vector)) for vector in hyperplane_basis]
    transported = solve_rational(projected_basis, target)
    if any(dot(transported, row) <= 0 for row in coordinates):
        raise ArithmeticError("the transported covector must stay strictly positive on the subcone")
    return coordinates, transported


# ---------------------------------------------------------------------------------------------
# the enumeration
# ---------------------------------------------------------------------------------------------
def facet_orbits(rows, permutations, positive_covector, *, stats=None, depth=0):
    """One facet normal per orbit of cone(rows) under the group generated by `permutations`.

    Returns one record per orbit with keys `incidence`, `normal`, `orbit_size`, `depth`; when
    `stats` is a dict, accumulates `pivots`, `dual_description_calls`, `orbit_steps`,
    `descent_steps`, `max_subcone_rays`, `pivots_by_depth`.
    """
    if stats is None:
        stats = {}
    for key in ("pivots", "dual_description_calls", "orbit_steps", "descent_steps", "max_subcone_rays"):
        stats.setdefault(key, 0)
    stats.setdefault("pivots_by_depth", {})
    stats["dual_description_calls"] += 1
    stats["max_subcone_rays"] = max(stats["max_subcone_rays"], len(rows))

    d = len(rows[0])
    degree = len(rows)
    orbit_records: list[dict] = []
    visited: set[tuple[int, ...]] = set()
    queue: deque[int] = deque()

    def count_pivot():
        stats["pivots"] += 1
        stats["pivots_by_depth"][depth] = stats["pivots_by_depth"].get(depth, 0) + 1

    def insert(incidence):
        key = tuple(sorted(incidence))
        if key in visited:
            return None
        closure, steps = orbit_closure(key, permutations)
        stats["orbit_steps"] += steps
        record = {"incidence": key, "normal": None, "orbit_size": len(closure), "depth": depth}
        visited.update(closure)
        orbit_records.append(record)
        queue.append(len(orbit_records) - 1)
        return len(orbit_records) - 1

    # Seed: one facet, found from the positive covector by the two-dimensional descent.
    seed_normal = one_facet_by_descent(rows, positive_covector, stats=stats)
    insert(certify_facet(rows, seed_normal))

    while queue:
        index = queue.popleft()
        record = orbit_records[index]
        incidence = record["incidence"]
        normal = record["normal"] if record["normal"] is not None else normal_of_incidence(rows, incidence)
        if certify_facet(rows, normal) != incidence:
            raise ArithmeticError("the stored incidence is not the zero set of the stored normal")
        record["normal"] = normal

        if d == 2:
            ridges = [()]                      # a two-dimensional pointed cone has one ridge: the apex
        else:
            stabilizer = stabilizer_generators(incidence, permutations, degree)
            coordinates, transported = subcone_level(rows, incidence, normal, positive_covector)
            restricted = restrict_generators(stabilizer, incidence)
            ridge_records = facet_orbits(
                coordinates, restricted, transported, stats=stats, depth=depth + 1
            )
            ridges = [ridge["incidence"] for ridge in ridge_records]

        for ridge in ridges:
            global_ridge = tuple(incidence[position] for position in ridge)
            other = pivot_across_ridge(rows, normal, global_ridge)
            count_pivot()
            insert(incidence_of(rows, other))
    return orbit_records
