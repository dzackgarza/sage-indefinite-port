r"""U + (Z[sqrt 2], Tr) has Gram blocks U(-1) and diag(2, -4), which is U + <2, -4> of issue 40.

Every primitive null vector reaches e, by the Euclidean reduction. So a product g of generators
factors as g = h^-1 s with s = h g fixing e.
"""

from random import Random

from sage_extensions.all import *
from sage.matrix.constructor import matrix as sage_matrix
from sage.modules.free_module_element import vector as sage_vector
from sage.rings.integer_ring import ZZ as SageZZ
from sage_indefinite_port.indefinite.witt_one_hilbert import HilbertWittOneLattice


def test_gram_is_the_issue_lattice_for_d_two():
    B = HilbertWittOneLattice(int(2)).gram()
    expected = sage_matrix(SageZZ, 4, 4, [0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, -4])
    assert B == expected


def test_generators_are_integral_isometries():
    for D in (int(2), int(3), int(5)):
        lattice = HilbertWittOneLattice(D)
        B = lattice.gram()
        for g in lattice.generators():
            assert g.transpose() * B * g == B
            assert int(abs(g.det())) == int(1)


def test_null_vectors_reach_e_and_products_factor():
    rng = Random(int(3))
    for D in (int(2), int(5)):
        lattice = HilbertWittOneLattice(D)
        generators = lattice.generators()
        e = sage_vector(SageZZ, [SageZZ(int(1)), SageZZ(int(0)), SageZZ(int(0)), SageZZ(int(0))])
        for _trial in range(int(10)):
            product = generators[int(0)] ** int(2)
            for _letter in range(int(10)):
                product = product * generators[rng.randrange(len(generators))]
            v = product * e
            h = lattice.transporter_to_e(v)
            assert h * v == e
            assert (h * product) * e == e
