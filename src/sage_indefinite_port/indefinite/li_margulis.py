r"""Bounded searches of Li-Margulis, with exact verification (#55, #56).

Theorem 4 and Lemma 20 of Li-Margulis bound one member of every orbit of primitive embeddings
of a Gram matrix `A` into `L`. Given a box containing such members, this module lists the
primitive tuples in the box with Gram matrix `A` and keeps one per `O(L)`-orbit, by the exact
marked transporter. The result is a statement about the supplied box: the bound is the input.
"""

from itertools import product

from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ

from sage_indefinite_port.indefinite.pairing_slices import marked_transporter


def _vectors_in_box(lattice, bound):
    generators = tuple(lattice.module_generators())
    for coordinates in product(range(-bound, bound + 1), repeat=len(generators)):
        yield sum(
            (lattice.scalar_multiple(c, g) for c, g in zip(coordinates, generators, strict=True)), lattice.zero()
        )


def primitive_tuples_in_box(lattice, gram, bound):
    r"""Yield the primitive tuples with Gram matrix ``gram`` and coordinates of absolute value at most ``bound``."""
    gram = matrix(ZZ, gram)
    size = gram.nrows()
    shell = {}
    for index in range(size):
        target = gram[index, index]
        shell[index] = [v for v in _vectors_in_box(lattice, bound) if int(lattice.b(v, v)) == int(target)]

    def extend(prefix):
        index = len(prefix)
        if index == size:
            rows = matrix(ZZ, [[int(c) for c in v] for v in prefix])
            if rows.rank() == size and all(d == 1 for d in rows.elementary_divisors()):
                yield tuple(prefix)
            return
        for candidate in shell[index]:
            if all(int(lattice.b(prefix[i], candidate)) == int(gram[i, index]) for i in range(index)):
                yield from extend([*prefix, candidate])

    yield from extend([])


def embedding_orbits_in_box(lattice, gram, bound):
    r"""Return one tuple per ``O(L)``-orbit among the primitive tuples in the box."""
    representatives = []
    for candidate in primitive_tuples_in_box(lattice, gram, bound):
        if all(marked_transporter(lattice, known, candidate) is None for known in representatives):
            representatives.append(candidate)
    return tuple(representatives)


def is_isometric_in_box(source, target, bound):
    r"""Return an isometry `source -> target` with images in the box, or ``None`` if the box has none."""
    frame = tuple(source.module_generators())
    gram = [[int(source.b(x, y)) for y in frame] for x in frame]
    labels = tuple(source.module_generating_set())
    for images in primitive_tuples_in_box(target, gram, bound):
        rows = matrix(ZZ, [[int(c) for c in v] for v in images])
        if abs(rows.det()) == 1:
            table = dict(zip(labels, images, strict=True))
            return source.Hom(target)(lambda label: table[label])
    return None
