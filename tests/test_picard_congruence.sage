r"""Principal and proper congruence subgroups of SU(2,1; Z[i]) and of Zhao's d = 7 group.

Gamma(m) has index the order of the finite image, every Schreier generator reduces to the
identity, and every one is a group element. The preimage of the image of a generator's cyclic
group has that image order fewer cosets: index = |image| / |<g> mod m|.
"""

from sage_extensions.all import *
from sage.matrix.special import identity_matrix as sage_identity
from sage_indefinite_port.indefinite import euclidean_picard, gaussian_picard
from sage_indefinite_port.indefinite.picard_congruence import congruence_generators, finite_image, reduction


def test_gaussian_principal_subgroup_of_level_two():
    generators = tuple(gaussian_picard.GENERATORS.values())
    one = sage_identity(gaussian_picard.FIELD, int(3))
    reduce = reduction(int(2))
    image = finite_image(generators, reduce, one)
    index, schreier = congruence_generators(generators, reduce, one, [one])
    assert int(index) == len(image) and len(image) > int(1)
    assert all(reduce(g) == reduce(one) and gaussian_picard.is_element(g) for g in schreier)


def test_gaussian_preimage_of_a_cyclic_image():
    generators = tuple(gaussian_picard.GENERATORS.values())
    one = sage_identity(gaussian_picard.FIELD, int(3))
    reduce = reduction(int(2))
    image = finite_image(generators, reduce, one)
    rotation = gaussian_picard.GENERATORS["M"]
    cyclic = finite_image((rotation,), reduce, one)
    index, schreier = congruence_generators(generators, reduce, one, list(cyclic.values()))
    assert int(index) * len(cyclic) == len(image)
    assert all(reduce(g) in cyclic for g in schreier)


def test_zhao_principal_subgroup_of_level_two():
    d = int(7)
    generators = tuple(euclidean_picard.generators(d).values())
    one = sage_identity(euclidean_picard.field(d), int(3))
    basis = euclidean_picard.order_basis(d)
    # The coordinates in the maximal order are integral on 1, omega.
    from sage.matrix.constructor import matrix as sage_matrix
    from sage.rings.rational_field import QQ as SageQQ

    change = sage_matrix(SageQQ, [list(b.vector()) for b in basis]).inverse()

    def reduce(element):
        return tuple(int(c) % int(2) for entry in element.list() for c in entry.vector() * change)

    image = finite_image(generators, reduce, one)
    index, schreier = congruence_generators(generators, reduce, one, [one])
    assert int(index) == len(image)
    assert all(reduce(g) == reduce(one) and euclidean_picard.is_element(d, g) for g in schreier)
