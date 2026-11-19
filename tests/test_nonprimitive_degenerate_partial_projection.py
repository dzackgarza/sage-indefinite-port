"""Finite-index markings of a primitive mixed sublattice with proper quotient projection."""

from sage_extensions.all import Lattices, ZZ
from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_sublattice_classes


def test_homothetic_index_four_submodule_of_partial_radical_marking_is_unique():
    r"""In L=U+<0>^2, S=<e+f,r1> is primitive and has
    proper nondegenerate projection <e+f> in U, with
    S intersect rad(L)=<r1> a proper radical summand.
    Any rank-two submodule J of S of Smith type (2,2)
    equals 2S: indeed J is contained in 2S and both
    have index four in S. The orbit classification is
    therefore one submodule irrespective of its group.
    """
    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0,0],[0,0]])
    e, f, r1, r2 = tuple(ambient.module_generators())
    source = ambient.subobject_on((e+f,r1))
    assert source.inclusion().is_primitive()
    classes = tuple(nonprimitive_sublattice_classes(source, (2,2)))
    assert len(classes) == 1
    inclusion = classes[0]
    expected = ambient.subobject_on((2*(e+f),2*r1)).inclusion()
    assert all(inclusion.is_in_image(expected(point))
               for point in expected.domain().module_generators())
    assert all(expected.is_in_image(inclusion(point))
               for point in inclusion.domain().module_generators())


def test_homothetic_index_eight_submodule_of_full_radical_proper_projection():
    r"""Let L=U+R, R zero of rank two, and S=<e+f>+R.
    The quotient projection P=<e+f> is proper inside U, although
    rank(S)-rank(U)=1. Every index-eight submodule of Smith
    type (2,2,2) is necessarily 2S: it is contained in 2S
    and has the same index eight. Thus the finite-index
    classification must not assume the full-rank P=U branch.
    """
    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0,0],[0,0]])
    e, f, r1, r2 = tuple(ambient.module_generators())
    source = ambient.subobject_on((e+f,r1,r2))
    assert source.inclusion().is_primitive()
    classes = tuple(nonprimitive_sublattice_classes(source,(2,2,2)))
    assert len(classes)==1
    inclusion = classes[0]
    expected = ambient.subobject_on((2*(e+f),2*r1,2*r2)).inclusion()
    assert all(inclusion.is_in_image(expected(point))
               for point in expected.domain().module_generators())
    assert all(expected.is_in_image(inclusion(point))
               for point in inclusion.domain().module_generators())
