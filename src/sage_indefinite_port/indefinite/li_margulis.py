r"""Bounded searches of Li-Margulis, with exact verification (#55, #56).

Theorem 4 and Lemma 20 of Li-Margulis bound one member of every orbit of primitive embeddings
of a Gram matrix `A` into `L`. Given a box containing such members, this module lists the
primitive tuples in the box with Gram matrix `A` and keeps one per `O(L)`-orbit, by the exact
marked transporter. The result is a statement about the supplied box: the bound is the input.
"""

from itertools import product

from sage.matrix.constructor import matrix
from sage.functions.other import ceil
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from sage_indefinite_port.indefinite.pairing_slices import marked_transporter


def generator_bound(gram, constant):
    r"""Return the integer box size of Theorem 2, `C_r |\det B|^{r^6} \|B\|^{3 r^4}`, rounded up and made strict.

    ``constant`` is a certified majorant of the dimension constant `C_r`. Section 1.7 of the
    source computes it from the proof and appendix, and gives no numerical value.
    """
    gram = matrix(ZZ, gram)
    rank = gram.nrows()
    if rank < 3 or gram.det() == 0:
        raise ValueError("Theorem 2 needs a nonsingular form of rank at least three")
    norm = max(abs(entry) for entry in gram.list())
    return ceil(constant * abs(gram.det()) ** (rank**6) * norm ** (3 * rank**4)) + 1


def isometry_box_bound(left, right, constant, anisotropic=False):
    r"""Return an upper bound on the entries of `g` with `g^t B_2 g = B_1`, from Theorem 1.

    Equation (3) applies in general and has the negative determinant exponent of the source.
    Equation (4) is the anisotropic case, for a caller who has established anisotropy.
    The right-hand norm factors are replaced by the monotone majorant `\|B\|^{r^4}`, so
    the box is larger than the source's and remains valid.
    """
    left, right = matrix(ZZ, left), matrix(ZZ, right)
    rank = right.nrows()
    if rank < 3 or left.det() == 0 or right.det() == 0:
        raise ValueError("Theorem 1 needs nonsingular forms of rank at least three")
    norm = max(abs(entry) for entry in (*left.list(), *right.list()))
    determinant = abs(right.det())
    numerator = 13 * rank**2 + 8 * rank - 26 if anisotropic else -(13 * rank**3 - 39 * rank**2 + 36 * rank + 52)
    denominator = 10 * rank if anisotropic else 20 * rank
    exponent = max(QQ(numerator) / denominator, 0)
    return ceil(constant * determinant ** ceil(exponent) * norm ** (rank**4)) + 1


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
