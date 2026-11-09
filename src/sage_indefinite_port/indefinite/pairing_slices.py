r"""Orbits of a pointwise stabilizer on vectors, by prescribed pairing slices (#53).

For a fixed tuple `s` of vectors, `\Gamma = Stab_{O(L)}(s)` preserves every `b(s_i, x)` and
`q(x)`. These invariants give exact negative decisions. On one slice, `x` and `y` are in one
`\Gamma`-orbit exactly when the tuples `(s, x)` and `(s, y)` are `O(L)`-equivalent as
markings of the same Gram lattice, which is a decision on coslice objects.
"""

from sage_extensions.categories.lattices import Lattices
from sage_indefinite_port.indefinite.marked_lattices import isometry


def _marking(lattice, vectors):
    gram = [[lattice.b(left, right) for right in vectors] for left in vectors]
    lattices = Lattices(lattice.base_ring())
    source = lattices(gram)
    labels = tuple(source.module_generating_set())
    images = dict(zip(labels, vectors, strict=True))
    category = lattices.CosliceUnder(source)
    return category, category.object(source.Mor(lattice)(lambda label: images[label]))


def slice_invariant(lattice, fixed, vector):
    r"""Return the pairings with the fixed tuple and the square of ``vector``."""
    return tuple(lattice.b(s, vector) for s in fixed), lattice.q(vector)


def pointwise_transporter(lattice, fixed, left, right):
    r"""Return `g \in Stab_{O(L)}(fixed)` with `g(left) = right`, or ``None``."""
    if slice_invariant(lattice, fixed, left) != slice_invariant(lattice, fixed, right):
        return None
    category, source = _marking(lattice, (*fixed, left))
    _category, target = _marking(lattice, (*fixed, right))
    target = category.object(_rebase(category, target))
    witness = isometry(category, source, target)
    if witness is None:
        return None
    assert all(witness(s) == s for s in fixed) and witness(left) == right
    return witness


def _rebase(category, marked):
    r"""Return the arrow of ``marked`` from the base object of ``category``."""
    base = category.base_object()
    arrow = marked.arrow()
    lattice = arrow.codomain()
    labels = tuple(arrow.domain().module_generating_set())
    images = {label: arrow(arrow.domain().module_generator(label)) for label in labels}
    return base.Mor(lattice)(lambda label: images[label])
