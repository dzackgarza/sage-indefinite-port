r"""Picard groups act on primitive null vectors with a single orbit up to units (FFLP, Zhao).

Images of the standard cusp under random words reduce to it, by construction of the source reduction.
"""

from random import Random

from sage_extensions.all import *
from sage.modules.free_module_element import vector as sage_vector
from sage_indefinite_port.indefinite import euclidean_picard as ep
from sage_indefinite_port.indefinite import gaussian_picard as gp


def test_gaussian_images_of_the_cusp_reduce():
    rng = Random(int(21))
    names = tuple(gp.GENERATORS)
    base = sage_vector(gp.FIELD, [gp.ZERO, gp.ZERO, gp.ONE])
    for _trial in range(int(8)):
        product = gp.evaluate([(rng.choice(names), int(1)) for _ in range(int(12))])
        v = product * base
        g = gp.null_transporter(v)
        assert gp.is_element(g)
        image = g * v
        assert image[int(0)] == gp.ZERO and image[int(1)] == gp.ZERO and image[int(2)].norm() == gp.ONE


def test_zhao_images_of_the_cusp_reduce():
    rng = Random(int(22))
    for d in (int(2), int(7), int(11)):
        table = ep.generators(d)
        names = tuple(table)
        E = ep.field(d)
        base = sage_vector(E, [E.one(), E.zero(), E.zero()])
        for _trial in range(int(6)):
            product = table["R1"] ** int(2)
            for _letter in range(int(10)):
                product = product * table[rng.choice(names)]
            v = product * base
            g = ep.null_transporter(d, v)
            assert ep.is_element(d, g)
            image = g * v
            assert image[int(1)] == E.zero() and image[int(2)] == E.zero() and image[int(0)] in (E.one(), -E.one())
