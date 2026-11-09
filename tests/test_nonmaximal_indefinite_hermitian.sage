r"""U(Lambda, J) for Lambda = Z[2i]^2 inside Lambda' = Z[i]^2, J the hyperbolic plane.

Lambda has basis (1, 0), (2i, 0), (0, 1), (0, 2i). The saturation group U(1,1; Z[i]) has exact
generators, and U(Lambda) is the stabilizer of Lambda / 2 Lambda'. The index must equal the
orbit size of that submodule under the finite image, which the test enumerates directly.
"""

from sage_extensions.all import *
from sage.matrix.constructor import block_matrix, matrix as sage_matrix
from sage.matrix.special import identity_matrix as sage_identity
from sage.rings.finite_rings.finite_field_constructor import GF as SageGF
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational_field import QQ as SageQQ
from sage_extensions.algebras.order_hermitian_modules import OrderHermitianModule
from sage_indefinite_port.indefinite.binary_hyperbolic_unitary import generators, is_element


def _trace_matrix(g):
    r"""Return ``g`` on the `Z`-basis (1, 0), (i, 0), (0, 1), (0, i) of `Z[i]^2`."""
    rows = []
    for i in range(int(2)):
        real = []
        imaginary = []
        for j in range(int(2)):
            a, b = g[i, j].vector()
            real += [a, -b]
            imaginary += [b, a]
        rows += [real, imaginary]
    return sage_matrix(SageZZ, rows)


def _module():
    half = SageQQ(int(1)) / SageQQ(int(2))
    i_on_second = sage_matrix(SageQQ, int(2), int(2), [SageQQ(int(0)), SageQQ(int(-2)), half, SageQQ(int(0))])
    zero = sage_matrix(SageQQ, int(2), int(2))
    i_action = block_matrix(SageQQ, [[i_on_second, zero], [zero, i_on_second]], subdivide=False)
    gram = sage_matrix(SageZZ, int(4), int(4), [int(v) for v in (0, 0, 2, 0, 0, 0, 0, 8, 2, 0, 0, 0, 0, 8, 0, 0)])
    return OrderHermitianModule(gram, [SageQQ(int(2)) * i_action], [sage_identity(SageQQ, int(4)), i_action])


def test_stabilizer_index_matches_the_finite_orbit():
    d = int(1)
    module = _module()
    assert int(module.conductor_exponent) == int(2)
    saturated = sage_matrix(SageQQ, module.saturation_basis)
    half = SageQQ(int(1)) / SageQQ(int(2))
    frame = sage_matrix.diagonal(SageQQ, [SageQQ(int(1)), half, SageQQ(int(1)), half])
    change = saturated.inverse() * frame
    moves = [change * _trace_matrix(g).change_ring(SageQQ) * change.inverse() for g in generators(d).values()]
    index, schreier = module.stabilizer_of_module(moves)
    assert all(module.is_unitary(g) for g in schreier)

    field = SageGF(int(2))
    one = sage_identity(field, int(4))
    one.set_immutable()
    letters = [sage_matrix(SageZZ, m).change_ring(field) for m in moves]
    known = set([one])
    frontier = [one]
    while frontier:
        current = frontier.pop()
        for letter in letters:
            image = letter * current
            image.set_immutable()
            if image not in known:
                known.add(image)
                frontier.append(image)
    inclusion = module._inclusion.change_ring(field)
    span = lambda m: (m * inclusion).transpose().row_module()
    base = span(one)
    assert int(index) == len(set(span(g) for g in known))
