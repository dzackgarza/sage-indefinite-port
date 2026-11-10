r"""Level structure on the integral parabolic of the Gaussian line in Z[i]^3, J = antidiag(1, 1, 1).

The Schreier generators of Gamma(2) in the setwise stabilizer reduce to the identity and preserve
the lattice and the form, and the index is the order of the finite image. The pointwise
stabilizer fixes the marked vector, and the setwise stabilizer contains a unit action on it.
"""

from sage_extensions.all import *
from sage.matrix.constructor import matrix as sage_matrix
from sage.matrix.special import identity_matrix as sage_identity
from sage.rings.number_field.number_field import CyclotomicField as SageCyclotomicField
from sage_extensions.categories.hermitian_lattices import HermitianLattice, HermitianWittFrame
from sage_indefinite_port.indefinite.picard_congruence import congruence_generators, finite_image, reduction


def _frame():
    field = SageCyclotomicField(int(4))
    zero, one = field.zero(), field.one()
    gram = sage_matrix(field, int(3), int(3), [zero, zero, one, zero, one, zero, one, zero, zero])
    lattice = HermitianLattice((one, field.gen()), gram)
    return field, HermitianWittFrame(lattice, sage_identity(field, int(3)), sage_matrix(field, int(1), int(1), [one]), int(1))


def test_level_two_subgroup_of_the_setwise_stabilizer():
    field, frame = _frame()
    generators = tuple(frame.stabilizer_generators(False))
    one = sage_identity(field, int(3))
    reduce = reduction(int(2))
    image = finite_image(generators, reduce, one)
    index, schreier = congruence_generators(generators, reduce, one, [one])
    assert int(index) == len(image) and len(image) > int(1)
    assert all(reduce(g) == reduce(one) and frame.preserves_lattice(g) for g in schreier)
    assert all(g.conjugate_transpose() * frame.j() * g == frame.j() for g in schreier)


def test_pointwise_stabilizer_fixes_the_marked_vector():
    field, frame = _frame()
    e = sage_matrix(field, int(3), int(1), [field.one(), field.zero(), field.zero()])
    assert all(g * e == e for g in frame.stabilizer_generators(True))
    units = [g[int(0), int(0)] for g in frame.stabilizer_generators(False)]
    assert field.gen() in units or -field.gen() in units


def test_gaussian_line_stabilizer_equals_the_source_cusp_group_up_to_scalars():
    r"""Z[i]^3 with the source form C of Falbel-Francsics-Lax-Parker and the null line of e_3.

    The cusp group P of arXiv:0911.1104 is the stabilizer of that line in SU(2,1; Z[i]), generated
    by N0, N1 and M (Proposition 3.2). Scalar i has determinant -i, so U = SU x <i> on this rank:
    an element of the stabilizer in U is a power of i times an element of P.

    Forward: each generator, normalized to determinant 1, has a word in N0, N1, M.
    Reverse: each of N0, N1, M is a lift of a Levi pair times a product of kernel generators.
    """
    from sage.modules.free_module_element import vector as sage_vector
    from sage.rings.integer_ring import ZZ as SZZ
    from sage.rings.rational_field import QQ as SQQ
    from sage_indefinite_port.indefinite.gaussian_picard import (
        CUSP_STABILIZER_GENERATORS, FIELD, I, ONE, ZERO, evaluate, is_element, lattice, word,
    )
    hermitian = lattice()
    frame = hermitian.witt_frame([[ZERO, ZERO, ONE]])
    e3 = sage_matrix(FIELD, int(3), int(1), [ZERO, ZERO, ONE])
    units = [I**int(k) for k in range(int(4))]
    generators = frame.stabilizer_generators(False)
    assert generators
    for g in generators:
        scalar = next(I**int(a) for a in range(int(4)) if (I**int(a) * g).det() == ONE)
        element = scalar * g
        assert is_element(element)
        assert element * e3 in [u * e3 for u in units]
        expression = word(element)
        assert all(name in ("N0", "N1", "M") for name, _ in expression)
        assert evaluate(expression) == element

    kernel, centre = frame.integral_unipotent_kernel()
    index, pairs = frame.levi_image_schreier()
    assert int(index) == int(1)
    lifted = [(frame.induced_action(frame.in_lattice(frame.levi(a, u) * frame.unipotent(*frame.levi_lift(a, u)))), frame.in_lattice(frame.levi(a, u) * frame.unipotent(*frame.levi_lift(a, u)))) for a, u in pairs]
    identity = sage_identity(FIELD, int(3))
    closure = {(ONE, ONE): identity}
    frontier = [(ONE, ONE)]
    while frontier:
        key = frontier.pop()
        for (a, u), g in lifted:
            image = (key[int(0)] * a[int(0), int(0)], key[int(1)] * u[int(0), int(0)])
            if image not in closure:
                closure[image] = closure[key] * g
                frontier.append(image)

    flatten = lambda value: list(value.vector())
    rows = sage_matrix(SQQ, [flatten(pair[int(0)][int(0), int(0)]) for pair in kernel])
    scale = rows.denominator()
    echelon, transform = (scale * rows).change_ring(SZZ).echelon_form(transformation=True)
    rank = echelon.rank()
    hermite = echelon[:rank].change_ring(SQQ) / scale
    kernel_matrices = [frame.in_lattice(frame.unipotent(*pair)) for pair in kernel]
    central = frame.in_lattice(frame.unipotent(*centre[int(0)]))
    for cusp in CUSP_STABILIZER_GENERATORS:
        a, u = frame.induced_action(cusp)
        levi_element = closure[(a[int(0), int(0)], u[int(0), int(0)])]
        remainder = ~levi_element * cusp
        witt = frame._frame.inverse() * remainder * frame._frame
        Z = witt[int(1), int(2)]
        coordinates = hermite.solve_left(sage_vector(SQQ, flatten(Z)))
        assert all(c in SZZ for c in coordinates)
        exponents = sage_vector(SZZ, [SZZ(c) for c in coordinates]) * transform[:rank]
        product = identity
        for matrix_, exponent in zip(kernel_matrices, exponents, strict=True):
            product = product * matrix_ ** int(exponent)
        rest = ~product * remainder
        multiples = [central**int(s) for s in range(int(-6), int(7))]
        assert rest in multiples
