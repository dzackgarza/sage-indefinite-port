r"""Orbits of a pointwise stabilizer on vectors, by prescribed pairing slices (#53).

For a fixed tuple `s` of vectors, `\Gamma = Stab_{O(L)}(s)` preserves every `b(s_i, x)` and
`q(x)`. These invariants give exact negative decisions. On one slice, `x` and `y` are in one
`\Gamma`-orbit exactly when the tuples `(s, x)` and `(s, y)` are `O(L)`-equivalent as
markings of the same Gram lattice, which is a decision on coslice objects.
"""

from itertools import count, product

from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage_extensions.categories.lattices import Lattices
from sage_indefinite_port.indefinite.marked_lattices import isometry


def _marking(lattice, vectors):
    gram = [[lattice.b(left, right) for right in vectors] for left in vectors]
    lattices = Lattices(lattice.base_ring())
    source = lattices(gram)
    labels = tuple(source.module_generating_set())
    images = dict(zip(labels, vectors, strict=True))
    category = lattices.CosliceUnder(source)
    return category, category.object(lattices.Mor(source, lattice)(lambda label: images[label]))


def _saturated_bases(lattice, source, target):
    r"""Return bases of the saturation of the source span and its image, or ``None``.

    A map fixing `source_i -> target_i` is rational-linear on the span, extends uniquely to
    the saturation, and exists only when the relations match and the image is a saturated
    integral basis, since an isometry carries a saturated sublattice to a saturated one.
    """
    rows = matrix(QQ, [[int(c) for c in v] for v in source])
    images = matrix(QQ, [[int(c) for c in v] for v in target])
    pivots = rows.pivot_rows()
    independent = rows[list(pivots)]
    transfer = independent.solve_left(rows)
    if transfer * images[list(pivots)] != images:
        return None
    saturated = matrix(ZZ, independent.denominator() * independent).saturation().change_ring(QQ)
    moved = independent.solve_left(saturated) * images[list(pivots)]
    if moved.denominator() != 1 or any(d != 1 for d in matrix(ZZ, moved).elementary_divisors()):
        return None
    generators = tuple(lattice.module_generators())

    def assemble(row):
        return sum(
            (lattice.scalar_multiple(int(c), g) for c, g in zip(row, generators, strict=True)),
            lattice.zero(),
        )

    return tuple(assemble(r) for r in saturated.rows()), tuple(assemble(r) for r in moved.rows())


def slice_invariant(lattice, fixed, vector):
    r"""Return the pairings with the fixed tuple and the square of ``vector``."""
    return tuple(lattice.b(s, vector) for s in fixed), lattice.q(vector)


def marked_transporter(lattice, source, target):
    r"""Return `g \in O(L)` with `g(source_i) = target_i` for tuples, or ``None``."""
    gram = lambda vectors: [[lattice.b(x, y) for y in vectors] for x in vectors]
    if gram(source) != gram(target):
        return None
    bases = _saturated_bases(lattice, source, target)
    if bases is None:
        return None
    mixed = _mixed_marking_transporter(lattice, *bases)
    if mixed is not NotImplemented:
        if mixed is not None:
            assert all(mixed(x) == y for x, y in zip(source, target, strict=True))
        return mixed
    category, marked_source = _marking(lattice, bases[0])
    _category, marked_target = _marking(lattice, bases[1])
    marked_target = category.object(_rebase(category, marked_target))
    witness = isometry(category, marked_source, marked_target)
    if witness is None:
        return None
    assert all(witness(x) == y for x, y in zip(source, target, strict=True))
    return witness


def _mixed_marking_transporter(lattice, source, target):
    r"""Decide a marking whose span has a nonzero radical and a nonzero nondegenerate part.

    Choose a nonsingular principal block `S` of the Gram matrix. An element of `O(L)` moves
    the block of ``source`` to that of ``target``, and any solution differs from it by the
    pointwise stabilizer of `S'`. When `L = S' + K`, that stabilizer is `O(K)`, and the
    remaining vectors project to a totally isotropic tuple of `K`. Returns
    ``NotImplemented`` when the marking is not mixed.
    """
    gram = matrix(ZZ, [[int(lattice.b(x, y)) for y in source] for x in source])
    rank = gram.rank()
    if rank == 0 or rank == len(source):
        return NotImplemented
    pivots = tuple(gram.pivot_rows())
    block = tuple(source[i] for i in pivots)
    block_target = tuple(target[i] for i in pivots)
    move = marked_transporter(lattice, block, block_target)
    if move is None:
        return None
    frame = tuple(lattice.module_generators())
    coordinates = lambda v: vector(ZZ, [int(c) for c in v])
    ambient = matrix(ZZ, [[int(lattice.b(x, y)) for y in frame] for x in frame])
    block_matrix = matrix(ZZ, [coordinates(v) for v in block_target])
    complement = (block_matrix * ambient).right_kernel().basis_matrix()
    together = block_matrix.stack(complement)
    assert abs(together.det()) == 1, (
        f"{block_target} does not split off {lattice}: its stabilizer needs the glue of the complement"
    )
    assemble = lambda row: sum(
        (lattice.scalar_multiple(int(c), g) for c, g in zip(row, frame, strict=True)), lattice.zero()
    )
    complement_vectors = tuple(assemble(r) for r in complement.rows())
    inverse = together.change_ring(QQ).inverse()
    tail = len(block_target)

    def in_complement(v):
        return vector(ZZ, list(coordinates(v) * inverse)[tail:])

    rest_source = tuple(move(source[i]) for i in range(len(source)) if i not in pivots)
    rest_target = tuple(target[i] for i in range(len(source)) if i not in pivots)
    complement_lattice = Lattices(lattice.base_ring())(
        [[int(lattice.b(x, y)) for y in complement_vectors] for x in complement_vectors]
    )
    complement_frame = tuple(complement_lattice.module_generators())
    build = lambda row: sum(
        (complement_lattice.scalar_multiple(int(c), g) for c, g in zip(row, complement_frame, strict=True)),
        complement_lattice.zero(),
    )
    image_source = tuple(build(in_complement(v)) for v in rest_source)
    image_target = tuple(build(in_complement(v)) for v in rest_target)
    inner = marked_transporter(complement_lattice, image_source, image_target)
    if inner is None:
        return None
    labels = tuple(lattice.module_generating_set())
    images = {}
    for label, generator in zip(labels, frame, strict=True):
        parts = coordinates(generator) * inverse
        head = sum((int(c) * together.row(i) for i, c in enumerate(parts[:tail])), vector(ZZ, len(frame)))
        moved = inner(build(parts[tail:]))
        body = sum(
            (int(c) * complement.row(j) for j, c in enumerate(moved)),
            vector(ZZ, len(frame)),
        )
        images[label] = assemble(head + body)
    return lattice.O()(images) * move


def pointwise_transporter(lattice, fixed, left, right):
    r"""Return `g \in Stab_{O(L)}(fixed)` with `g(left) = right`, or ``None``.

    ``left`` and ``right`` are vectors or tuples of vectors of equal length.
    """
    left = left if isinstance(left, tuple) else (left,)
    right = right if isinstance(right, tuple) else (right,)
    return marked_transporter(lattice, (*fixed, *left), (*fixed, *right))


def line_transporter(lattice, generator, left, right):
    r"""Return `g` in the setwise stabilizer of `Z generator` with `g(left) = right`, or ``None``.

    `g` sends `generator` to its positive or its negative, so two markings decide.
    """
    for image in (generator, -generator):
        witness = marked_transporter(lattice, (generator, left), (image, right))
        if witness is not None:
            return witness
    return None


def setwise_transporter(lattice, fixed, left, right):
    r"""Return `g` preserving the span of ``fixed`` with `g(left) = right`, or ``None``.

    The tuple ``fixed`` is a basis of a primitive definite sublattice `S`. Restriction sends the
    setwise stabilizer into the finite group `O(S)`, so running over `O(S)` decides.
    """
    gram = [[int(lattice.b(x, y)) for y in fixed] for x in fixed]
    source = Lattices(lattice.base_ring())(gram)
    assert source.is_definite(), f"{fixed} does not span a definite sublattice"
    frame = tuple(source.module_generators())
    for automorphism in source.O():
        images = tuple(
            sum(
                (lattice.scalar_multiple(int(c), f) for c, f in zip(automorphism(v), fixed, strict=True)),
                lattice.zero(),
            )
            for v in frame
        )
        witness = marked_transporter(lattice, (*fixed, left), (*images, right))
        if witness is not None:
            return witness
    return None


def _rebase(category, marked):
    r"""Return the arrow of ``marked`` from the base object of ``category``."""
    base = category.base_object()
    arrow = marked.arrow()
    lattice = arrow.codomain()
    labels = tuple(arrow.domain().module_generating_set())
    images = {label: arrow(arrow.domain().module_generator(label)) for label in labels}
    return category.base_category().Mor(base, lattice)(lambda label: images[label])


class PointwiseStabilizerOrbits:
    r"""The orbit set of `Stab_{O(L)}(fixed)` on primitive vectors of a given square.

    The set can be infinite, so it has no list of representatives. Iteration visits all
    vectors by sup-norm height and yields the first vector of each orbit. Every primitive
    vector of the square has a height, and the transporter decides equivalence, so the
    stream is exhaustive and free of duplicates.
    """

    def __init__(self, lattice, fixed, square):
        self.lattice = lattice
        self.fixed = tuple(fixed)
        self.square = square

    def __contains__(self, vector):
        return self.lattice.q(vector) == self.square and vector.content() == 1

    def equivalent(self, left, right):
        return pointwise_transporter(self.lattice, self.fixed, left, right) is not None

    def transporter(self, left, right):
        return pointwise_transporter(self.lattice, self.fixed, left, right)

    def _shell(self, height):
        generators = tuple(self.lattice.module_generators())
        for coordinates in product(range(-height, height + 1), repeat=len(generators)):
            if max((abs(c) for c in coordinates), default=0) == height:
                yield sum(
                    (self.lattice.scalar_multiple(c, g) for c, g in zip(coordinates, generators, strict=True)),
                    self.lattice.zero(),
                )

    def __iter__(self):
        known = []
        for height in count(0):
            for vector in self._shell(height):
                if vector == self.lattice.zero() or vector not in self:
                    continue
                if not any(self.equivalent(representative, vector) for representative in known):
                    known.append(vector)
                    yield vector


class SubgroupVectorOrbits:
    r"""The orbit set of a caller-fixed subgroup `\Gamma` on primitive vectors of a given square.

    `\Gamma` needs a ``vector_equivalence_witness``, which is its exact decision. The stream visits
    vectors by height and keeps the first of each orbit, so it is exhaustive and free of
    duplicates, with no finiteness assumed.
    """

    def __init__(self, lattice, group, square):
        self.lattice = lattice
        self.group = group
        self.square = square

    def __contains__(self, vector):
        return int(self.lattice.q(vector)) == int(self.square) and int(vector.content()) == 1

    def transporter(self, left, right):
        return self.group.vector_equivalence_witness(left, right)

    def equivalent(self, left, right):
        return self.transporter(left, right) is not None

    def __iter__(self):
        frame = tuple(self.lattice.module_generators())
        known = []
        for height in count(0):
            for coordinates in product(range(-height, height + 1), repeat=len(frame)):
                if max((abs(c) for c in coordinates), default=0) != height or not any(coordinates):
                    continue
                vector = sum(
                    (self.lattice.scalar_multiple(c, g) for c, g in zip(coordinates, frame, strict=True)),
                    self.lattice.zero(),
                )
                if vector in self and not any(self.equivalent(representative, vector) for representative in known):
                    known.append(vector)
                    yield vector

