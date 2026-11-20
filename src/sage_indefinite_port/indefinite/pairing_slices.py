r"""Orbits of a pointwise stabilizer on vectors, by prescribed pairing slices (#53).

For a fixed tuple `s` of vectors, `\Gamma = Stab_{O(L)}(s)` preserves every `b(s_i, x)` and
`q(x)`. These invariants give exact negative decisions. On one slice, `x` and `y` are in one
`\Gamma`-orbit exactly when the tuples `(s, x)` and `(s, y)` are `O(L)`-equivalent as
markings of the same Gram lattice, which is a decision on coslice objects.

For the centralizer `\Gamma = C_{O(L)}(F_1, \dots, F_s)` of arbitrary isometries, which has infinite
index in general, :class:`CentralizerSubgroup` decides vector and tuple equivalence through the
equivariant marked transporter of `(L, F, v)` (:class:`CentralizerTupleOrbit`, over #54).
"""

from itertools import count, product

from sage.matrix.constructor import identity_matrix, matrix
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage_extensions.categories.lattices import Lattices
from sage_indefinite_port.indefinite.equivariant_markings import CentralizerTupleOrbit
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
    source, target = tuple(source), tuple(target)
    if len(source) != len(target):
        return None
    if not source or all(point == lattice.zero() for point in source):
        # The empty marking has the entire orthogonal group as stabilizer;
        # the all-zero marking has the same stabilizer, but cannot be sent
        # to a tuple containing a nonzero vector.
        return lattice.O().identity() if all(point == lattice.zero() for point in target) else None
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
    if not together.is_square() or abs(together.det()) != 1:
        # A primitive nondegenerate summand need not split off integrally.
        # Its orthogonal complement can meet it in an overlattice with
        # nontrivial discriminant glue (e.g. <e+f> in U has index two).
        # The direct-product O(S) x O(S^perp) construction below is exact
        # only at index one. Pass a genuinely glued marking to the full
        # based Coslice Isom algorithm, which checks the integral gluing
        # conditions rather than asserting or discarding the orbit.
        return NotImplemented
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


def _full_radical_setwise_transporter(lattice, fixed, left, right):
    r"""Exact setwise transporter when S=P+R contains R=rad(L), P nondegenerate.

    The reduced lattice L=C+R is an integral direct sum of a nondegenerate
    quotient and a free zero-form radical. Since S contains the whole
    radical, its image P in C is a primitive sublattice and S=P+R;
    it has no additional graph/gluing datum over R.

    Every integral isometry has the unique row block form
        [[A,B],[0,D]], A in O(C), B in Hom_Z(C,R), D in GL(R).
    It preserves S precisely when A(P)=P. For x=(c,r), y=(c',s),
    first decide c A=c' by the complete nondegenerate setwise
    transporter of P. Then solve rD+cB=s. If c has content d>0,
    the set of possible cB is exactly dR by Bezout, so a solution
    exists iff the *actual* integral GL(R) image carries r mod d
    to s mod d. The maintained RadicalCongruenceImage decides this
    and returns an integral lift D; one Bezout column constructs B.
    At d=0, the shear disappears and GL(R) must carry r exactly
    to s, decided by the complete pure-radical pointed transporter.

    Every decision is necessary and sufficient, and the constructed
    map is checked on the original ambient lattice and marked span.
    Return NotImplemented unless the full-radical, nondegenerate-P
    hypothesis holds, so the caller may choose another valid route.
    """
    from math import gcd

    from sage_indefinite_port.indefinite.indefinite_reduction import indefinite_reduced_presentation
    from sage_indefinite_port.indefinite.recursive import (
        EichlerRecursion, _gram_in_label_order, _isometry_from_row_matrix,
    )

    reduction = indefinite_reduced_presentation(lattice)
    reduced = reduction.reduced
    gram = _gram_in_label_order(reduced)
    size, n = int(gram.nrows()), int(gram.rank())
    k = size - n
    if not n or not k:
        return NotImplemented
    labels = tuple(reduced.module_generating_set())
    inclusion = lattice.subobject_on(tuple(fixed)).inclusion()
    radical = tuple(
        reduction.isometry(reduced.module_generator(label))
        for label in labels[n:]
    )
    if not inclusion.is_primitive() or not all(inclusion.is_in_image(point) for point in radical):
        return NotImplemented

    complement = Lattices(lattice.base_ring())([
        [int(gram[i,j]) for j in range(n)] for i in range(n)
    ])
    c_labels = tuple(complement.module_generating_set())
    inverse = ~reduction.isometry

    def chart(point):
        row = inverse(point).to_vector()
        return tuple(int(row(label)) for label in labels)

    def c_point(coords):
        return complement.linear_combination({
            label: complement.base_ring()(coefficient)
            for label, coefficient in zip(c_labels, coords, strict=True) if coefficient
        })

    projected = tuple(c_point(chart(point)[:n]) for point in fixed)
    plane = complement.subobject_on(projected)
    if not plane.inclusion().is_primitive() or not plane.is_nondegenerate():
        return NotImplemented
    if int(plane.rank()) + k != int(inclusion.domain().rank()):
        raise ArithmeticError("the marked full-radical lattice is not the pullback of its projected sublattice")

    x, y = chart(left), chart(right)
    c, c_target = x[:n], y[:n]
    r, t = x[n:], y[n:]
    if not any(c) and not any(c_target):
        a = identity_matrix(ZZ, n)
    else:
        projected_basis = tuple(plane.inclusion()(point) for point in plane.module_generators())
        a_witness = setwise_transporter(complement, projected_basis, c_point(c), c_point(c_target))
        if a_witness is None:
            return None
        a = matrix(ZZ, [
            [int(a_witness(point).to_vector()(label)) for label in c_labels]
            for point in complement.module_generators()
        ])

    d = 0
    for entry in c:
        d = gcd(d, entry)
    recursion = EichlerRecursion()
    if d > 1:
        D = recursion.radical_congruence_image(lattice, d).transporter(r, t)
        if D is None:
            return None
    elif d == 1 or r == t:
        D = identity_matrix(ZZ, k)
    else:
        pure_source = reduced.linear_combination({
            labels[n+i]: reduced.base_ring()(entry)
            for i, entry in enumerate(r) if entry
        })
        pure_target = reduced.linear_combination({
            labels[n+i]: reduced.base_ring()(entry)
            for i, entry in enumerate(t) if entry
        })
        radical_witness = recursion.pointed_isometry(pure_source, pure_target)
        if radical_witness is None:
            return None
        D = matrix(ZZ, [
            [int(radical_witness(reduced.module_generator(labels[n+i]))
                 .to_vector()(labels[n+j])) for j in range(k)]
            for i in range(k)
        ])

    residual = matrix(ZZ, [list(t)]) - matrix(ZZ, [list(r)]) * D
    if d:
        if any(int(entry) % d for entry in residual.list()):
            raise ArithmeticError("the radical GL congruence lift violates the prescribed residue")
        _, _, smith_right = matrix(ZZ, [list(c)]).smith_form()
        bezout = smith_right.column(0)
        if sum(c[i]*int(bezout[i]) for i in range(n)) == -d:
            bezout = -bezout
        if sum(c[i]*int(bezout[i]) for i in range(n)) != d:
            raise ArithmeticError("the projected content lacks its integral Bezout column")
        B = matrix(ZZ, n, k, lambda i,j: int(bezout[i])*int(residual[0,j])//d)
    else:
        if any(residual.list()):
            return None
        B = matrix(ZZ, n, k, 0)

    block = identity_matrix(ZZ, size)
    block[:n,:n] = a
    block[:n,n:] = B
    block[n:,n:] = D
    reduced_witness = _isometry_from_row_matrix(reduced, reduced, block)
    witness = reduction.isometry * reduced_witness * inverse
    if witness(left) != right:
        raise ArithmeticError("the mixed setwise transporter failed the prescribed vector equation")
    if any(not inclusion.is_in_image(witness(point)) for point in fixed):
        raise ArithmeticError("the mixed setwise transporter failed its marked lattice equation")
    return witness


def setwise_transporter(lattice, fixed, left, right):
    r"""Return `g` preserving the span of ``fixed`` with `g(left) = right`, or ``None``.

    The tuple ``fixed`` is a basis of a primitive nondegenerate sublattice `S`. For definite
    `S`, restriction sends the setwise stabilizer into the finite group `O(S)`, so running over
    `O(S)` decides. For indefinite `S`, the decision splits along `S + S^\perp`.
    """
    gram = [[int(lattice.b(x, y)) for y in fixed] for x in fixed]
    source = Lattices(lattice.base_ring())(gram)
    if not source.is_nondegenerate():
        found = _full_radical_setwise_transporter(lattice, fixed, left, right)
        if found is NotImplemented:
            raise NotImplementedError(
                "setwise transport of a degenerate marking not containing the entire "
                "ambient radical with nondegenerate quotient projection requires "
                "the general equivariant gluing algorithm"
            )
        return found
    if not source.is_definite():
        return _glued_setwise_transporter(lattice, fixed, left, right)
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


def _glued_setwise_transporter(lattice, fixed, left, right):
    r"""Decide the setwise transporter for a nondegenerate primitive sublattice `S`.

    Let `K = S^\perp`. Then `S + K \subseteq L \subseteq S^\vee + K^\vee`, and `L` is the
    graph of an anti-isometry of subgroups of `A_S` and `A_K` (Nikulin, *Integral symmetric
    bilinear forms*, Prop. 1.5.1). An isometry `g` with `g(S) = S` also preserves `K`, so
    `g = h + k` with `h \in O(S)` and `k \in O(K)`; conversely `h + k` lies in `O(L)` exactly
    when it maps `L` into `L`. Write `x = a + c` with `a \in S^\vee_Q` and `c \in K^\vee_Q`.
    Then `g(x) = y` exactly when `h(a) = a'` and `k(c) = c'`. The solutions `h` form a coset
    `h_0 Stab_{O(S)}(a)`, and only its image on the finite group `S^\vee/S` affects whether
    `h + k` preserves `L`. The same holds for `k`. So finitely many pairs of classes decide.

    Matrices act on coefficient rows from the right: `p \mapsto p M`.
    """
    frame = tuple(lattice.module_generators())
    coordinates = lambda v: vector(QQ, [int(c) for c in v])
    ambient = matrix(QQ, [[int(lattice.b(x, y)) for y in frame] for x in frame])
    fixed_rows = matrix(QQ, [coordinates(v) for v in fixed])
    complement_rows = matrix(ZZ, fixed_rows * ambient).right_kernel().basis_matrix().change_ring(QQ)
    together = fixed_rows.stack(complement_rows)
    inverse = together.inverse()
    tail = fixed_rows.nrows()
    left_parts = coordinates(left) * inverse
    right_parts = coordinates(right) * inverse
    cosets = []
    for rows, source, target in (
        (fixed_rows, left_parts[:tail], right_parts[:tail]),
        (complement_rows, left_parts[tail:], right_parts[tail:]),
    ):
        if rows.nrows() == 0:
            cosets.append((matrix(QQ, 0, 0),))
            continue
        coset = _dual_vector_coset(rows * ambient * rows.transpose(), source, target)
        if coset is None:
            return None
        cosets.append(coset)
    for first in cosets[0]:
        for second in cosets[1]:
            moved = inverse * matrix.block_diagonal([first, second], subdivide=False) * together
            if moved.denominator() != 1:
                continue
            images = {
                label: sum(
                    (lattice.scalar_multiple(int(c), g) for c, g in zip(moved.row(i), frame, strict=True)),
                    lattice.zero(),
                )
                for i, label in enumerate(lattice.module_generating_set())
            }
            witness = lattice.O()(images)
            assert witness(left) == right
            return witness
    return None


def _dual_vector_coset(gram, source, target):
    r"""Return matrices of `\{h \in O(S) : h(a) = a'\}`, one for each class on `S^\vee/S`.

    ``source`` and ``target`` are coefficient rows of `a, a' \in S^\vee_Q`. Returns ``None``
    when no `h` exists. Two matrices `M_1, M_2` act equally on `S^\vee/S` exactly when
    `G^{-1}(M_1 - M_2)` is integral, since the rows of `G^{-1}` span `S^\vee`.
    """
    if source.denominator() != target.denominator():
        return None
    abstract = Lattices(ZZ)([[int(c) for c in row] for row in gram])
    generators = tuple(abstract.module_generators())
    rank = len(generators)
    scale = source.denominator()

    def element(row):
        return sum(
            (abstract.scalar_multiple(int(c), g) for c, g in zip(scale * row, generators, strict=True)),
            abstract.zero(),
        )

    def as_matrix(isometry):
        return matrix(QQ, [[int(c) for c in isometry(g)] for g in generators])

    group = abstract.O()
    if source.is_zero():
        if not target.is_zero():
            return None
        start = matrix.identity(QQ, rank)
        stabilizer = tuple(group.group_generators())
    else:
        move = group.vector_equivalence_witness(element(source), element(target))
        if move is None:
            return None
        start = as_matrix(move)
        stabilizer = tuple(group.vector_stabilizer_generators(element(source)))
    dual = gram.inverse()
    key = lambda m: tuple(c - c.floor() for c in (dual * m).list())
    steps = tuple(as_matrix(s) for s in stabilizer)
    identity = matrix.identity(QQ, rank)
    seen = {key(identity): identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for step in steps:
            following = current * step
            if key(following) not in seen:
                seen[key(following)] = following
                frontier.append(following)
    return tuple(m * start for m in seen.values())


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



class CentralizerSubgroup:
    r"""`C_{O(L)}(F_1, \dots, F_s)` for isometries `F_i` of `L`, with its exact vector and sublattice decisions.

    Elements are elements of ``lattice.O()``; internally an isometry is the integer matrix whose
    columns are the coordinates of the images of the module generators.
    """

    def __init__(self, lattice, isometries):
        self.lattice = lattice
        self.frame = tuple(lattice.module_generators())
        self.labels = tuple(lattice.module_generating_set())
        self.gram = matrix(ZZ, [[int(lattice.b(x, y)) for y in self.frame] for x in self.frame])
        self.isometries = tuple(isometries)
        self.matrices = tuple(self._matrix(F) for F in self.isometries)

    def _matrix(self, element):
        return matrix(ZZ, [[int(c) for c in element(g)] for g in self.frame]).transpose()

    def _columns(self, vectors):
        return matrix(ZZ, [[int(c) for c in v] for v in vectors]).transpose()

    def _element(self, columns):
        images = {
            label: sum(
                (self.lattice.scalar_multiple(int(c), g) for c, g in zip(columns.column(j), self.frame, strict=True)),
                self.lattice.zero(),
            )
            for j, label in enumerate(self.labels)
        }
        return self.lattice.O()(images)

    def __contains__(self, element):
        g = self._matrix(element)
        return g.transpose() * self.gram * g == self.gram and all(g * F == F * g for F in self.matrices)

    def _orbit(self, vectors):
        return CentralizerTupleOrbit(self.gram, self.matrices, self._columns(vectors))

    def tuple_equivalence_witness(self, left, right):
        r"""Return `g \in C` with `g(left_i) = right_i`, or ``None`` as a proof that none exists."""
        found = self._orbit(left).transporter(self._columns(right))
        if found is None:
            return None
        witness = self._element(found)
        assert all(witness(x) == y for x, y in zip(left, right, strict=True))
        return witness

    def vector_equivalence_witness(self, left, right):
        r"""Return `g \in C` with `g(left) = right`, or ``None`` as a proof that none exists."""
        return self.tuple_equivalence_witness((left,), (right,))

    def vector_stabilizer_generators(self, vector):
        return tuple(self._element(g) for g in self._orbit((vector,)).stabilizer_generators())

    def pointwise_stabilizer_generators(self, vectors):
        return tuple(self._element(g) for g in self._orbit(tuple(vectors)).stabilizer_generators())

    def setwise_stabilizer_generators(self, basis):
        r"""Return generators of `\{g \in C : g \Lambda = \Lambda\}` for the sublattice `\Lambda` with ``basis``."""
        return tuple(self._element(g) for g in self._orbit(tuple(basis)).setwise_stabilizer_generators())
