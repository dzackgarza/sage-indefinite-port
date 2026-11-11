r"""Isometry test and Kneser neighbours over `O_K`, `K = Q(\sqrt 5)` (#54)."""

from sage_extensions.all import *
from sage.matrix.constructor import matrix as sage_matrix
from sage.matrix.special import block_diagonal_matrix, identity_matrix as sage_identity
from sage.modules.free_module_element import vector as sage_vector
from sage.rings.number_field.number_field import QuadraticField
from sage.combinat.root_system.cartan_matrix import CartanMatrix
from sage_indefinite_port.indefinite.genus_number_field import Cusps, lattice_isometry, standard_z_basis
from sage_indefinite_port.indefinite.witt_one_number_field import WittOneLattice


def _e8(field):
    return sage_matrix(field, CartanMatrix(["E", int(8)]))


def test_odd_unimodular_rank_nine_genus_has_two_inequivalent_lattices_over_sqrt5():
    r"""`I_9` and `E_8 \perp \langle 1 \rangle` over `O_K` are odd unimodular of rank 9, in one genus (base change of one `Z`-genus), and
    not isometric: `a \in O_K` with `a^2 \preceq 1` at both places is `0, \pm 1` (`\varphi = 1.618\ldots`), so `I_9 \otimes O_K` has the
    18 vectors `\pm \epsilon_i` of norm 1, while `h(y) + c^2 = 1` with `h(y) \in 2 O_K` on `E_8 \otimes O_K` forces `y = 0`, `c = \pm 1`:
    2 vectors. A unimodular basis change of `I_9` is isometric, with the returned images reproducing the Gram matrix."""
    field = QuadraticField(int(5), "s")
    nine = sage_identity(field, int(9))
    other = block_diagonal_matrix(_e8(field), sage_matrix(field, [[int(1)]]))
    basis = standard_z_basis(field, int(9))
    assert lattice_isometry(field, int(1), nine, basis, other, basis) is None
    phi = (int(1) + field.gen()) / int(2)
    T = sage_identity(field, int(9))
    T[int(0), int(1)] = phi
    T[int(1), int(2)] = int(1)
    changed = T.transpose() * nine * T
    images = lattice_isometry(field, int(1), nine, basis, changed, basis)
    assert images is not None
    assert len(images) == int(18)


def test_sqrt5_neighbour_of_three_squares_is_three_squares():
    r"""`\Lambda = U \perp I_3` over `O_K`, `\mathfrak p = (\sqrt 5)`, `x = m_1 + 2 m_2` (`h(x, x) = 5 \in \mathfrak p^2`). The cusp
    `u = \sqrt 5 e + \lambda f + x'` (`x' \equiv x \bmod \mathfrak p`) is isotropic of divisor `O_K` and `u^\perp/u` is the neighbour
    `N = M_x + O_K x/\sqrt 5`. By hand `N \cong I_3`: `m_3 \in N`, and `N \cap \langle m_1, m_2 \rangle_K` has the index of `O_K^2` (index
    `5` down to `M_x`, `5` up to `N`), so it is unimodular of rank 2 containing `x/\sqrt 5` of norm 1, whose complement is `\langle u \rangle` with
    `u` a totally positive unit, a square `\varphi^{2j}`; so `N \cong \langle 1, 1, 1 \rangle`. The isometry test confirms it."""
    field = QuadraticField(int(5), "s")
    root = field.gen()
    gram = block_diagonal_matrix(sage_matrix(field, [[int(0), int(1)], [int(1), int(0)]]), sage_identity(field, int(3)))
    lattice = WittOneLattice(field, gram, [int(1)] * int(5))
    cusps = Cusps(lattice)
    prime = field.ideal(root)
    x0 = sage_vector(field, [int(0), int(0), int(1), int(2), int(0)])
    e = sage_vector(field, [int(1), int(0), int(0), int(0), int(0)])
    f = sage_vector(field, [int(0), int(1), int(0), int(0), int(0)])
    u = cusps.neighbour_vector(e, f, list(lattice.z_basis), x0, prime, root)
    assert u is not None and lattice.h(u, u) == int(0) and u[int(0)] == root
    partner, complement, _even = cusps.frame(u)
    assert lattice.h(u, partner) == int(1)
    assert lattice_isometry(field, int(1), gram, list(lattice.z_basis), gram, complement) is not None
