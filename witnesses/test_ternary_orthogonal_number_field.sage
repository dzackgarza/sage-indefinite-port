r"""`O(M)` for the anisotropic ternary `M = \langle 1, -11, -19 \rangle` over `K = Q(\sqrt 5)` (#54).

`C_0(M) \cong (11, 19)_K`. Over `Q`, `(11, 19)` ramifies at `11` (`(19/11) = (8/11) = -1`) and `2` (`(-1)^{5 \cdot 9} = -1`), not at `19`
(`(11/19) = -(19/11) = 1`) or `\infty` (both positive). In `K`, `2` is inert (residue degree 2 kills the local invariant `1/2`) and `11 = 4^2 - 5`
splits, so `(11, 19)_K` ramifies exactly at the two primes over `11`: a division algebra unramified at both real places. So `M` is
anisotropic and indefinite at both real places, `SO(M \otimes R_\sigma) \cong SO(1, 2)` is noncompact, and the arithmetic group `O(M)` is
infinite. An element of finite order of `GL_3(K)` has eigenvalues `\zeta_n` with `[K(\zeta_n) : K] \le 2`, so `\varphi(n) \le 4` and `g^{120} = 1`;
some generator therefore has `g^{120} \ne 1`."""

from sage_extensions.all import *
from sage.matrix.constructor import matrix as sage_matrix
from sage.matrix.special import identity_matrix as sage_identity
from sage.rings.number_field.number_field import QuadraticField
from sage_indefinite_port.indefinite.genus_number_field import standard_z_basis
from sage_indefinite_port.indefinite.ternary_orthogonal_number_field import TernaryOrthogonalGroup


def test_ternary_orthogonal_group_of_one_minus_eleven_minus_nineteen_over_sqrt5():
    field = QuadraticField(int(5), "s")
    gram = sage_matrix(field, [[int(1), int(0), int(0)], [int(0), int(-11), int(0)], [int(0), int(0), int(-19)]])
    group = TernaryOrthogonalGroup(field, gram, standard_z_basis(field, int(3)))
    finite, infinite = group.algebra.ramified_places(inf=True)
    assert len(finite) == int(2) and all(p.norm() == int(11) for p in finite) and len(infinite) == int(0)
    generators = group.generators()
    one = sage_identity(field, int(3))
    for g in generators:
        assert g.transpose() * gram * g == gram
        assert all(c.is_integral() for c in g.list())
    assert any(g ** int(120) != one for g in generators)
