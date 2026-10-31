"""Isometry prefilters and the sign normalization of the Eichler recursion.

Replaces the ``size_t`` hash of ``INDEF_FORM_Invariant`` (polyhedral_common
``src_indefinite/IndefiniteFormFundamental.h``) by structured records of preamble invariants.
A prefilter may reject or bucket candidates; it never decides equivalence.
"""

from __future__ import annotations

from dataclasses import dataclass

from sage_extensions.categories.lattices import Lattices
from sage_extensions.tensors.tensor import _engine_component_matrix
from sage.matrix.constructor import matrix as SageMatrix
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.modules.free_module_element import vector as SageVector
from sage.rings.integer_ring import ZZ as SageZZ


@dataclass(frozen=True)
class LatticePrefilter:
    """Rank, signature, parity and discriminant: equal for isometric lattices."""

    rank: int
    signature: tuple[int, int]
    even: bool
    discriminant: int

    @classmethod
    def from_lattice(cls, lattice: Lattices.ParentMethods) -> LatticePrefilter:
        positive, negative = lattice.signature_pair()
        return cls(int(lattice.rank()), (int(positive), int(negative)), bool(lattice.is_even()), int(lattice.discriminant()))


@dataclass(frozen=True)
class AttackProfile:
    """L or L(-1), whichever has positive index min(n_+, n_-)."""

    lattice: Lattices.ParentMethods
    signed_view: Lattices.ParentMethods
    sign: int
    positive_index: int
    negative_index: int

    @classmethod
    def from_lattice(cls, lattice: Lattices.ParentMethods) -> AttackProfile:
        positive, negative = (int(index) for index in lattice.signature_pair())
        if positive <= negative:
            return cls(lattice, lattice, 1, positive, negative)
        return cls(lattice, lattice.twist(-lattice.base_ring().one()), -1, negative, positive)


_GRAMS: dict[int, tuple[Lattices.ParentMethods, Matrix_integer_dense]] = {}


def _gram(lattice: Lattices.ParentMethods) -> Matrix_integer_dense:
    r"""The Gram matrix of a lattice on its framing, lowered once per lattice."""
    if id(lattice) not in _GRAMS:
        _GRAMS[id(lattice)] = (lattice, _engine_component_matrix(lattice.gram_tensor()).change_ring(SageZZ))
    return _GRAMS[id(lattice)][1]


@dataclass(frozen=True)
class VectorPrefilter:
    r"""O(L)-invariants of a vector v: those of INDEF_FORM_InvariantVector.

    # source: IndefiniteFormFundamental.h:124-138 INDEF_FORM_InvariantVector (IndefiniteFormsFundamental.g:91),
    # compared first by CombinedAlgorithms.h:1816-1817 / IndefiniteForms.g:1607 INDEF_FORM_EquivalenceVector.
    The source records the rank of L, q(v), the index of b(v, L) (``div``), the content
    (``RemoveFractionVectorPlusCoeff(v)``) and INDEF_FORM_Invariant of the Gram matrix of
    v^perp (IndefiniteFormFundamental.h:75-88: rank, the determinant of the nondegenerate part,
    the signature with its zero count, and parity). For a vector of L these signature and
    determinant entries are functions of L, q(v), the content and div(v); the parity of v^perp
    is the one independent entry, read here from an integral basis of v^perp.
    """

    norm: int
    content: int
    divisibility: int
    perpendicular_even: bool

    @classmethod
    def from_vector(cls, vector: Lattices.ElementMethods) -> VectorPrefilter:
        lattice = vector.parent()
        gram = _gram(lattice)
        coordinates = vector.to_vector()
        row = SageVector(SageZZ, [int(coordinates(label)) for label in lattice.module_generating_set()])
        perpendicular = SageMatrix(SageZZ, [gram * row]).right_kernel_matrix()
        perpendicular_even = all((basis * gram * basis) % 2 == 0 for basis in perpendicular.rows())
        return cls(int(vector.q()), int(vector.content()), int(vector.div()), perpendicular_even)
