"""Mixed degenerate markings across a genuine index-two orthogonal glue."""

from sage_extensions.all import Lattices, ZZ
from sage_indefinite_port.indefinite.pairing_slices import marked_transporter


def test_mixed_glued_plane_transporter_preserves_integral_overlattice():
    r"""In L=U+<0>, S=<e+f,r> has Gram diag(2,0), and
    S+S^perp=<e+f,e-f,r> has index 2 in L. The map acting as
    identity on U and negation on the radical is an integral
    isometry carrying (e+f,r) to (e+f,-r). A direct-sum-only
    mixed transporter cannot reject this admissible glued input.
    """
    lattice = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, r = tuple(lattice.module_generators())
    source, target = (e+f,r), (e+f,-r)
    witness = marked_transporter(lattice, source, target)
    assert witness is not None
    assert tuple(witness(point) for point in source) == target
    assert (~witness)(target[1]) == source[1]
    assert all(lattice.b(witness(x), witness(y)) == lattice.b(x, y)
               for x in (e,f,r) for y in (e,f,r))
