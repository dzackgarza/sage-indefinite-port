"""Independently derived low-Witt mathematical test cases.

Primary source: H. Li and G. A. Margulis, *Effective estimates on integral
quadratic forms*, author's PDF, https://li-han.weebly.com/uploads/5/1/3/3/51338781/equivalence.pdf,
printed pp. 2--3, Theorems 1 (equations (3), (4)) and 2.
The determinant exponent in (3) was checked against the RENDERED p. 2,
not against extracted text or the port's documentation. The author PDF
prints a NEGATIVE determinant exponent. For nonsingular integral B,
|det B| >= 1, hence its negative-power factor is at most one.

The Hilbert-form examples below are direct proofs from the definition of
field norm and the displayed Hermitian determinant; they do not assert
unverified completeness of an orbit or generating set.
"""

from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ

from sage_indefinite_port.indefinite.li_margulis import (
    generator_bound,
    isometry_box_bound,
    isometries_in_box,
)
from sage_indefinite_port.indefinite.witt_one_hilbert import HilbertWittOneLattice
from sage_extensions.all import Lattices


def test_theorem_one_equation_three_negative_determinant_majorant():
    r"""Li--Margulis, Theorem 1(3), printed p. 2: negative determinant power.

    The actual determinant exponent is
    -(13n^3-39n^2+36n+52)/(20n), equal to -160/60=-8/3
    at n=3. Since B is nonsingular integral, |det B|>=1 and
    |det B|^(-8/3)<=1. Replacing that factor by one is a
    valid (weaker) upper bound. The positive norm exponents
    410/40=41/4 and 430/40=43/4 round upward to 11 each.
    For A=B=diag(2,-2,6), both sup norms are 6.
    The multiplier 1 is only a formal arithmetic input, NOT
    the unknown theorem constant C_3.
    """
    B = matrix(ZZ, [[2, 0, 0], [0, -2, 0], [0, 0, 6]])
    assert isometry_box_bound(B, B, ZZ(1)) == 6**22 + 1


def test_theorem_one_equation_four_anisotropic_exponents():
    r"""Li--Margulis, printed p. 2, Theorem 1, equation (4).

    B=diag(1,1,-3) is anisotropic over Q: if x²+y²=3z²,
    reduction modulo 3 forces 3|x,y, hence 3|z, and infinite descent
    excludes nonzero integral solutions. For n=3, the determinant exponent
    (13n²+8n-26)/(10n)=115/30=23/6 rounds to 4; the A-norm
    exponent 1/2 rounds to 1 and B-norm exponent (n-1)/2 equals 1.
    Therefore the routine's integer majorant for A=B, formal multiplier
    one, is 3^4 * 3 * 3 + 1. This does not certify C_3=1.
    """
    B = matrix(ZZ, [[1, 0, 0], [0, 1, 0], [0, 0, -3]])
    assert isometry_box_bound(B, B, ZZ(1), anisotropic=True) == 3**4 * 3**2 + 1


def test_theorem_two_generator_bound_exponent_arithmetic():
    r"""Li--Margulis, printed p. 3, Theorem 2.

    The printed factor is |det B|^(n^6) ||B||^(3n^4).
    For B=diag(1,1,-2), n=3, determinant magnitude 2, norm 2,
    hence with formal multiplier one and the API's strict integer
    enlargement the result is 2^(729+243)+1. No value of C_3 is claimed.
    """
    B = matrix(ZZ, [[1, 0, 0], [0, 1, 0], [0, 0, -2]])
    assert generator_bound(B, ZZ(1)) == 2**972 + 1


def test_hilbert_gram_from_quadratic_field_norm():
    r"""Direct derivation from q(a,d,beta)=-2ad+2N(beta).

    For D=2, beta=x+y sqrt(2), N(beta)=x²-2y²; polarizing q
    gives Gram U(-1) + diag(2,-4). For D=5, omega=(1+sqrt(5))/2,
    omega+bar(omega)=1 and N(omega)=-1; hence the trace pairing
    on (1,omega) is [[2,1],[1,-2]]. All entries follow by expansion,
    independently of the implementation.
    """
    assert HilbertWittOneLattice(2).gram() == matrix(
        ZZ, [[0, -1, 0, 0], [-1, 0, 0, 0],
             [0, 0, 2, 0], [0, 0, 0, -4]])
    assert HilbertWittOneLattice(5).gram() == matrix(
        ZZ, [[0, -1, 0, 0], [-1, 0, 0, 0],
             [0, 0, 2, 1], [0, 0, 1, -2]])


def test_hilbert_elementary_shear_and_swap_explicit_action():
    r"""Exact Hermitian multiplication, independent of generator enumeration.

    For g=[[1,1],[0,1]], X=[[a,beta],[bar(beta),d]],
    gXg*= [[a+beta+bar(beta)+d,beta+d],[bar(beta)+d,d]].
    At D=2, Tr(beta)=2x, so (a,d,x,y) maps to
    (a+d+2x,d,x+d,y). For w=[[0,1],[-1,0]], wXw* maps
    (a,d,x,y) to (d,a,-x,y). Both preserve -2det(X) since
    det(g)=det(w)=1; this also proves their integral isometry status.
    """
    H = HilbertWittOneLattice(2)
    E = H.field
    shear = matrix(E, [[1, 1], [0, 1]])
    swap = matrix(E, [[0, 1], [-1, 0]])
    S = matrix(ZZ, [[1, 1, 2, 0], [0, 1, 0, 0],
                    [0, 1, 1, 0], [0, 0, 0, 1]])
    W = matrix(ZZ, [[0, 1, 0, 0], [1, 0, 0, 0],
                    [0, 0, -1, 0], [0, 0, 0, 1]])
    assert H.action(shear) == S
    assert H.action(swap) == W
    B = H.gram()
    assert S.transpose() * B * S == B
    assert W.transpose() * B * W == B
    assert S.det() == W.det() == 1


def test_hyperbolic_plane_box_exact_four_isometries():
    r"""Direct Diophantine classification, not an orbit-count oracle.

    In U with Gram [[0,1],[1,0]], a matrix [[a,b],[c,d]]
    preserves the Gram iff ac=bd=0 and ad+bc=1. Thus the
    only integral solutions are diag(1,1), diag(-1,-1),
    [[0,1],[1,0]], [[0,-1],[-1,0]]. All fit radius one.
    The public bounded search must return exactly these owned arrows.
    """
    L = Lattices(ZZ)("U")
    expected = {((1, 0), (0, 1)), ((-1, 0), (0, -1)),
                ((0, 1), (1, 0)), ((0, -1), (-1, 0))}
    found = tuple(isometries_in_box(L, 1))
    frame = tuple(L.module_generators())
    actual = {
        tuple(tuple(ZZ(c) for c in g(v)) for v in frame)
        for g in found
    }
    assert actual == expected
    assert len(found) == 4
    assert all(g.parent() is L.O() for g in found)
