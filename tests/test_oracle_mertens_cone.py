"""Mertens's cone stabilizer is not the normalized real-spinor kernel.

Primary source: Mertens, arXiv:1303.3478, ``hyperbolic.tex``
lines 125--130 (row action), 259--269 (Omega is the stabilizer of a
negative-cone sheet), and 575--579 (the connecting element c_{1,1}).
Only this single source-printed element is tested; no completeness of the
listed generating system is inferred.
"""

from sage.all import ZZ, identity_matrix, matrix

from tests.acceptance.consumer import isometry_from_rows, lattice


def test_mertens_c11_preserves_time_cone_but_not_real_spinor_kernel():
    """The positive-norm reflection has negative normalized spinor sign.

    For row vectors and A below, v=(30,-5,10) has vAv^t=50.
    The reflection is R=I-2 A v^t v/50, precisely the printed c_{1,1}.
    Consequently det(R)=-1 and R^2=I. Its normalized real spinor
    class is -q(v)/2=-25<0. On the other hand e1 and e1 R are
    timelike with B(e1,e1 R)=-26<0, so they belong to the same
    component of the negative cone (signature (2,1)).
    """
    gram = ((-1, -3, -1), (-3, 14, 8), (-1, 8, 11))
    rows = ((31, -5, 10), (96, -15, 32), (-48, 8, -15))
    A = matrix(ZZ, gram)
    R = matrix(ZZ, rows)
    v = matrix(ZZ, [[30, -5, 10]])
    assert A.det() == -155
    assert (v * A * v.transpose())[0, 0] == 50
    assert R == identity_matrix(ZZ, 3) - A * v.transpose() * v / 25
    assert R * A * R.transpose() == A
    assert R.det() == -1
    assert R * R == identity_matrix(ZZ, 3)
    assert A[0, 0] == -1
    assert (R * A)[0, 0] == -26

    L = lattice(gram)
    g = isometry_from_rows(L, rows)
    assert int(g.determinant()) == -1
    assert not L.O_plus().preimage_morphism()(g).is_one()
