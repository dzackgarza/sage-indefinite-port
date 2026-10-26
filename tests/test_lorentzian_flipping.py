"""Supporting faces of the integral cone x squared minus y squared >= 0."""

from sage.all import QQ, vector

from sage_indefinite_port.indefinite.lorentzian_cells import NONNEGATIVE_SQUARES, _kernel_flipping
from tests.acceptance.consumer import lattice


def test_flip_completes_the_supporting_face_at_x_equal_one() -> None:
    """The nonnegative integral vectors at x=1 are (1,-1), (1,0), (1,1)."""
    form = lattice([[1, 0], [0, -1]])
    for transverse in ([-1, 0, 1], [-QQ(3) / 2, QQ(1) / 2, 1]):
        rows, normal, direction, bound = _kernel_flipping(
            form,
            ((1, 1),),
            vector(QQ, [-1, 2, -1]),
            vector(QQ, transverse),
            NONNEGATIVE_SQUARES,
        )
        assert set(rows) == {(1, -1), (1, 0), (1, 1)}
        assert normal[1] > 0 and normal[0] == -normal[1] and normal[2] == 0
        assert direction[0] == bound and direction[1] == 0
