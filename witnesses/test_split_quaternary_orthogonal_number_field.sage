r"""`O(M)` for `M = \langle 1, -11, -19, 209 \rangle` over `K = Q(\sqrt 5)`, the norm form of `B = (11, 19)_K` (#54).

`\mathrm{nrd}(t + x i + y j + z k) = t^2 - 11 x^2 - 19 y^2 + 209 z^2`, so the discriminant `11 \cdot 19 \cdot 209 = 209^2` is a square and `C_0(M) = B \times B`.
`B` ramifies exactly at the two primes over `11` (see the ternary test: `(11, 19)_Q` ramifies at `\{2, 11\}`, `2` is inert and `11` splits in
`K`), so `B` is a division algebra unramified at both real places: `M` is anisotropic of signature `(2, 2)` at both places, `SO(2, 2)` is
noncompact, and `O(M)` is infinite. An element of `GL_4(K)` of finite order `n` has `\varphi(n) \le 8`, so `n \mid 5040` and `g^{5040} = 1`; some
generator has `g^{5040} \ne 1`. The standard lattice is stable under `x \mapsto \bar x` (`t, x, y, z \mapsto t, -x, -y, -z`), an integral improper
isometry, so some generator has determinant `-1`."""

from sage_extensions.all import *
from sage.matrix.constructor import matrix as sage_matrix
from sage.matrix.special import identity_matrix as sage_identity
from sage.rings.number_field.number_field import QuadraticField
from sage_indefinite_port.indefinite.genus_number_field import standard_z_basis
from sage_indefinite_port.indefinite.quaternary_orthogonal_number_field import SplitQuaternaryOrthogonalGroup


def test_split_quaternary_orthogonal_group_of_the_norm_form_of_eleven_nineteen_over_sqrt5():
    field = QuadraticField(int(5), "s")
    gram = sage_matrix(field, [[int(1), int(0), int(0), int(0)], [int(0), int(-11), int(0), int(0)], [int(0), int(0), int(-19), int(0)], [int(0), int(0), int(0), int(209)]])
    group = SplitQuaternaryOrthogonalGroup(field, gram, standard_z_basis(field, int(4)))
    finite, infinite = group.algebra.ramified_places(inf=True)
    assert len(finite) == int(2) and all(p.norm() == int(11) for p in finite) and len(infinite) == int(0)
    generators = group.generators()
    one = sage_identity(field, int(4))
    for g in generators:
        assert g.transpose() * gram * g == gram
        assert all(c.is_integral() for c in g.list())
    assert any(g.det() == int(-1) for g in generators)
    assert any(g ** int(5040) != one for g in generators)
