"""Integral-isometry oracles from indefinite lattice classification.

Original inputs (NOT verdicts): Indefinite.jl pinned at 374a5ebb,
TestLor/{U_I3,U_E8,U_2U_2I3}_mat{1,2}. Its test_gap.jl
CaseEquiv1/2/3 are all disabled; no positive verdict is inferred.

Independent theorem: indefinite odd unimodular integral lattices are
classified by signature and parity; indefinite even unimodular lattices
are uniquely determined by signature (Milnor--Husemoller, Symmetric
Bilinear Forms, chapter IV). For even nonunimodular indefinite lattices,
Nikulin, Integral symmetric bilinear forms and some of their applications,
Math. USSR Izv. 14 (1980), Theorem 1.14.2: an even indefinite lattice
with rank >= length(discriminant group)+2 is unique in its genus.

For the first pair the Gram forms are odd, unimodular, signature (1,4).
For the second they are even, unimodular, signature (1,9). For the third
they are even, signature (2,5), with elementary divisors
(1,1,2,2,2,2,2), so length=5 and rank=7=length+2. The third pair's
local genus agreement is independently checked using exact genus data;
the theorem then proves their integral isometry. Every returned witness
is separately checked by its defining integral Gram congruence.
"""

from pathlib import Path

from sage.all import ZZ, matrix

from tests.acceptance.consumer import MARKED_ISOMETRY, lattice, require, rows_of


_SOURCE = Path(__file__).resolve().parents[1] / "references/vendor/Indefinite.jl@374a5ebb/TestLor"


def _source_gram(stem, suffix):
    tokens = (_SOURCE / f"{stem}_mat{suffix}").read_text().split()
    n, m = map(int, tokens[:2])
    assert n == m
    assert len(tokens) == 2 + n * n
    return matrix(ZZ, n, n, [int(x) for x in tokens[2:]])


def _check_pair(stem, *, rank, determinant, even, length):
    require(MARKED_ISOMETRY)
    a, b = (_source_gram(stem, j) for j in (1, 2))
    assert a.nrows() == b.nrows() == rank
    for gram in (a, b):
        assert abs(gram.det()) == determinant
        assert all(gram[i, i] % 2 == 0 for i in range(rank)) == even
        assert sum(int(d) > 1 for d in gram.elementary_divisors()) == length
    left, right = lattice(a.rows()), lattice(b.rows())
    assert left.genus() == right.genus()
    assert left.is_isometric(right)
    witness = left.Isom(right).an_element()
    images = matrix(ZZ, rows_of(witness))
    assert abs(images.det()) == 1
    assert images * b * images.transpose() == a


def test_odd_unimodular_signature_one_four_unique_class():
    _check_pair("U_I3", rank=5, determinant=1, even=False, length=0)


def test_even_unimodular_signature_one_nine_unique_class():
    _check_pair("U_E8", rank=10, determinant=1, even=True, length=0)


def test_even_rank_seven_genus_unique_by_nikulin():
    _check_pair("U_2U_2I3", rank=7, determinant=32, even=True, length=5)
