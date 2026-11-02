"""Issue #41: the retained `<2> + 0` pointed/vector checks, as executable evidence.

`L = <2> + 0` with coordinates `(x, y)` and `q(x, y) = 2 x^2` is degenerate: its radical is the
`y`-coordinate line. The pinned GAP and C++ public vector routines reject a rank-deficient ambient form
before their nonzero-norm branch (IndefiniteForms.g:1555-1564, :1605-1622; CombinedAlgorithms.h:1750-1756,
:1813-1826), so these cases are an extension, not a translation. The pinned upstream is therefore not the
authority for the expected values here; the criterion is.

For `L = C + R` with `C = <2>`, `R = Z` (column convention `g = [[A,0],[B,D]]`, `A` in `O(C)`, `D` in
`GL(R)`, `B` in `Hom(C,R)`), a transporter `g(v) = w` for `v = (x, r)`, `w = (y, s)` exists if and only if

    A x = y  and  s - D r is in cR,   c = content_C(x),

followed by the linear equation `B x = s - D r`. Those two conditions are what the checks below exercise:
the positive ones exhibit the integral unimodular witness, and the negative one fails the congruence rather
than any norm, content or divisibility invariant (all three agree there).

The unmarked radical `Aut`/`Isom` route is implemented separately and must keep working; it is checked here
so a change to the marked path cannot quietly take it out.
"""
import pytest

from sage.all import MatrixGroup, ZZ, matrix, vector

from sage_indefinite_port.indefinite.recursive import EichlerRecursion, _column_matrix, shear_kernel_basis
from tests.acceptance.consumer import (
    automorphisms,
    element,
    is_isometry,
    isomorphisms,
    lattice,
    pointed_category,
    rows_of,
    underlying_isometry,
)

_RECURSION = EichlerRecursion()
_GRAM = [[2, 0], [0, 0]]


def _degenerate_lattice():
    return lattice(_GRAM)


def _vector(coordinates):
    lattice_ = _degenerate_lattice()
    return lattice_, element(lattice_, coordinates)


def _marking(category, lattice_, vector):
    return category.object(category.base_object().hom([vector], codomain=lattice_))


def test_transports_one_zero_to_one_one_by_an_integral_unimodular_shear() -> None:
    """`(1,0) -> (1,1)` is realised by exactly the shear `g_{1,1,1}`, of determinant 1."""
    lattice_, source = _vector([1, 0])
    _, target = _vector([1, 1])
    witness = _RECURSION.pointed_isometry(source, target)

    assert witness is not None
    assert witness(source) == target
    assert is_isometry(rows_of(witness), _GRAM)
    columns = _column_matrix(witness)
    assert columns == matrix(ZZ, [[1, 0], [1, 1]])
    assert int(columns.det()) == 1


def test_transports_five_one_to_five_six_by_the_same_shear() -> None:
    """`(5,1) -> (5,6)`: same shear. Here `c = 5` and `6 - 1 = 5` is in `5R`, so the congruence holds."""
    lattice_, source = _vector([5, 1])
    _, target = _vector([5, 6])
    witness = _RECURSION.pointed_isometry(source, target)

    assert witness is not None
    assert witness(source) == target
    assert _column_matrix(witness) == matrix(ZZ, [[1, 0], [1, 1]])


def test_equal_norm_content_and_divisibility_do_not_transport_five_one_to_five_two() -> None:
    """`(5,1)` is not equivalent to `(5,2)`: the condition `2 - delta r` in `5Z` fails for `delta = +-1`.

    The three invariants a naive test would use all agree, so the refusal must come from the congruence.
    """
    lattice_, source = _vector([5, 1])
    _, target = _vector([5, 2])

    assert int(source.q()) == int(target.q())
    assert int(source.content()) == int(target.content())
    assert int(source.div()) == int(target.div())

    assert _RECURSION.pointed_isometry(source, target) is None


def test_pointwise_stabilizer_of_one_zero_is_the_order_two_reflection_of_the_radical() -> None:
    """`O(L)_{(1,0)} = {g_{1,0,1}, g_{1,0,-1}}`, of order 2 and fixing the vector pointwise."""
    lattice_, vector = _vector([1, 0])
    generators = tuple(_RECURSION.pointed_automorphism_generators(vector))

    assert generators
    assert all(generator(vector) == vector for generator in generators)
    forms = [matrix(ZZ, rows_of(generator)) for generator in generators]
    assert MatrixGroup(forms).order() == 2
    # The non-identity element negates the radical coordinate and fixes nothing else.
    assert matrix(ZZ, [[1, 0], [0, -1]]) in forms


def test_public_pointed_consumers_reach_the_same_verdicts() -> None:
    """The owned pointed category, not the private recursion, is what a consumer sees."""
    lattice_, source = _vector([1, 0])
    _, target = _vector([1, 1])
    _, unrelated = _vector([5, 2])
    category = pointed_category(2)
    source_object = _marking(category, lattice_, source)
    target_object = _marking(category, lattice_, target)

    assert not isomorphisms(category, source_object, target_object).is_empty()
    morphism = isomorphisms(category, source_object, target_object).an_element()
    assert underlying_isometry(morphism)(source) == target

    automorphism_group = automorphisms(category, source_object)
    generators = [underlying_isometry(generator) for generator in automorphism_group.group_generators()]
    assert generators and all(generator(source) == source for generator in generators)
    assert MatrixGroup([matrix(ZZ, rows_of(g)) for g in generators]).order() == 2

    content_five = pointed_category(50)
    source_object = _marking(content_five, lattice_, element(lattice_, [5, 1]))
    assert isomorphisms(content_five, source_object,
                        _marking(content_five, lattice_, unrelated)).is_empty()


def test_marking_transports_in_a_presentation_where_the_radical_is_not_last() -> None:
    """The construction is about the lattice, not about a Gram block: permute the coordinates and it holds."""
    permuted = lattice([[0, 0], [0, 2]])
    source = element(permuted, [0, 1])
    target = element(permuted, [1, 1])
    assert int(source.q()) == int(target.q()) == 2

    witness = _RECURSION.pointed_isometry(source, target)
    assert witness is not None
    assert witness(source) == target
    assert is_isometry(rows_of(witness), [[0, 0], [0, 2]])


def test_unmarked_radical_automorphism_route_is_preserved() -> None:
    """The separately implemented unmarked radical `Aut`/`Isom` route must not regress."""
    lattice_, _ = _vector([1, 0])
    generators = tuple(lattice_.O().group_generators())

    assert generators
    forms = [matrix(ZZ, rows_of(generator)) for generator in generators]
    assert all(is_isometry(rows, _GRAM) for rows in [rows_of(g) for g in generators])
    # The unmarked group is generated by O(C) = {+-1}, GL(R) = {+-1} and the shear, so it contains
    # the radical negation and the shear as actual generators.
    assert matrix(ZZ, [[1, 0], [0, -1]]) in forms
    assert matrix(ZZ, [[1, 1], [0, 1]]) in forms or matrix(ZZ, [[1, 0], [1, 1]]) in forms
    assert MatrixGroup(forms).order() > 2


# ---------------------------------------------------------------------------------------------------
# Rank-two radical: the actual integral GL(R) image, and the nontrivial shear kernel.
#
# For rank(R) >= 2 the integral group GL(R) genuinely moves residues: on the specimen below the witness
# needs D = [[2,-1],[5,-2]] (determinant 1), not the identity, because s - r is not in 5R. The decision is
# the action of the *integral* generators on R/5R; over Z the image of GL(R) is what must be used, and at
# rank one it is a proper subgroup of Aut(R/cR), which is the (5,1) -/-> (5,2) refusal above.
# ---------------------------------------------------------------------------------------------------

_RANK_TWO_GRAM = [[2, 0, 0], [0, 0, 0], [0, 0, 0]]


def _split_blocks(columns, c_rank, r_rank):
    """A, B, D of the column matrix [[A, 0], [B, D]] of g = [[A,0],[B,D]]."""
    A = columns.submatrix(0, 0, c_rank, c_rank)
    B = columns.submatrix(c_rank, 0, r_rank, c_rank)
    D = columns.submatrix(c_rank, c_rank, r_rank, r_rank)
    return A, B, D


def test_rank_two_radical_transport_needs_a_nonidentity_integral_D() -> None:
    """`(5,1,0) -> (5,2,0)`: `s - r` is not in `5R`, so some `D` in GL_2(Z) must move `r` modulo 5."""
    lattice_ = lattice(_RANK_TWO_GRAM)
    source, target = element(lattice_, [5, 1, 0]), element(lattice_, [5, 2, 0])
    witness = _RECURSION.pointed_isometry(source, target)

    assert witness is not None
    assert witness(source) == target
    assert is_isometry(rows_of(witness), _RANK_TWO_GRAM)

    columns = _column_matrix(witness)
    A, shear, D = _split_blocks(columns, 1, 2)
    residue = vector(ZZ, [1, 0])          # the radical part of the source
    assert D * residue != residue, "this specimen is only interesting when D moves integer r"

    # The radical block must be an integral automorphism of the radical lattice: det(D) = +-1.
    assert D.change_ring(ZZ) == D
    assert int(D.det()) in (1, -1)
    # ... and the whole isometry must be unimodular, not just the D block.
    assert int(columns.det()) in (1, -1)
    assert A == matrix(ZZ, [[1]])

    # The criterion: s - D r lies in content_C(x) * R = 5R, and B x = s - D r.
    difference = vector(ZZ, [2, 0]) - D * residue
    assert all(int(entry) % 5 == 0 for entry in difference)
    assert shear * vector(ZZ, [5]) == difference


def test_rank_two_radical_refuses_when_r_vanishes_modulo_content_but_s_does_not() -> None:
    """`(5,5,0) -/-> (5,0,1)`: here `r = (5,0) = 0` in `R/5R`, so `D r = 0` for every `D`, and `s != 0`."""
    lattice_ = lattice(_RANK_TWO_GRAM)
    source, target = element(lattice_, [5, 5, 0]), element(lattice_, [5, 0, 1])

    assert int(source.q()) == int(target.q()) == 50
    assert _RECURSION.pointed_isometry(source, target) is None


def test_pure_shear_transport_when_d_minus_identity_already_lies_in_the_content_module() -> None:
    """`(5,1,0) -> (5,6,5)` needs no `D` at all: `s - r = (5,5)` is in `5R`, so `D = 1` and `B` carries it."""
    lattice_ = lattice(_RANK_TWO_GRAM)
    source, target = element(lattice_, [5, 1, 0]), element(lattice_, [5, 6, 5])
    witness = _RECURSION.pointed_isometry(source, target)

    assert witness is not None and witness(source) == target
    columns = _column_matrix(witness)
    A, shear, D = _split_blocks(columns, 1, 2)
    # A "pure shear" claim is only that if BOTH diagonal blocks are the identity: a nonzero B with a
    # non-identity A or D would be transporting through the radical, not shearing.
    assert A == matrix(ZZ, [[1]])
    assert D == matrix(ZZ, [[1, 0], [0, 1]])
    assert shear != 0
    assert shear * vector(ZZ, [5]) == vector(ZZ, [5, 5])


def test_pointwise_stabilizer_contains_the_nontrivial_shear_kernel() -> None:
    """With `C` of rank two the kernel `{B : B x = 0}` is nontrivial, so the stabilizer is infinite.

    For `C` of rank one the kernel is trivial: `B` is determined by `B(x_0)` and `c * B(x_0) = 0` forces
    `B = 0` in a torsion-free `R`. Rank two is therefore the first place the kernel shows up. The claimed
    kernel elements must have identity diagonal blocks -- that is what makes them shears rather than
    transports through the radical -- and they must span the *whole* kernel, not one nonzero example.
    """
    gram = [[2, 0, 0], [0, 2, 0], [0, 0, 0]]
    lattice_ = lattice(gram)
    marked = element(lattice_, [1, 0, 0])
    generators = tuple(_RECURSION.pointed_automorphism_generators(marked))

    assert generators
    assert all(generator(marked) == marked for generator in generators)

    identity_c = matrix(ZZ, [[1, 0], [0, 1]])
    shears = []
    for generator in generators:
        A, shear, D = _split_blocks(_column_matrix(generator), 2, 1)
        if shear == 0:
            continue
        assert A == identity_c, "a shear with a non-identity A is not in the kernel"
        assert D == matrix(ZZ, [[1]]), "a shear with a non-identity D is not in the kernel"
        assert shear * vector(ZZ, [1, 0]) == 0, "a kernel shear must annihilate the marked C-coordinate"
        shears.append(shear)
    assert shears, "the kernel {B : Bx = 0} should be nontrivial for rank(C) = 2"

    # The criterion's own basis: rank(C) - 1 = 1 generator per radical coordinate, saturated.
    expected = shear_kernel_basis(2, 1, (1, 0))
    assert len(expected) == 1
    expected_rows = matrix(ZZ, [[int(entry) for row in block.rows() for entry in row] for block in expected])
    assert expected_rows.rank() == 1

    # The port's shears lie in that kernel and generate all of it: same rank, and primitive (the row
    # lattice is saturated, so it is the kernel itself rather than a proper sublattice).
    found = matrix(ZZ, [[int(entry) for row in shear.rows() for entry in row] for shear in shears])
    assert found.rank() == 1
    invariants = [int(entry) for entry in found.smith_form()[0].diagonal() if int(entry) != 0]
    assert invariants == [1], f"the shears must span a saturated lattice, got invariants {invariants}"
    # Each found shear is an integral multiple of the criterion's generator.
    generator_row = expected_rows[0]
    for row in found.rows():
        assert (row[0] * generator_row[1] - row[1] * generator_row[0]) == 0


# ---------------------------------------------------------------------------------------------------
# The full kernel, as a general construction rather than one example.
# ---------------------------------------------------------------------------------------------------

@pytest.mark.parametrize(
    ("c_rank", "radical_rank", "x"),
    [
        (2, 1, (1, 0)),
        (2, 1, (5, 0)),          # content 5: the same primitive direction, so the same kernel
        (2, 1, (5, 15)),         # content 5, primitive direction (1, 3)
        (2, 2, (3, 6)),          # two radical coordinates
        (3, 2, (3, 4, 5)),       # rank(C) = 3: two kernel dimensions per radical coordinate
    ],
)
def test_shear_kernel_basis_generates_the_whole_kernel(c_rank, radical_rank, x) -> None:
    """The generators span the kernel itself: the right rank, and saturated (all Smith invariants 1)."""
    basis = shear_kernel_basis(c_rank, radical_rank, x)
    marked = vector(ZZ, x)

    assert len(basis) == radical_rank * (c_rank - 1)
    assert all(block * marked == 0 for block in basis)

    rows = [[int(entry) for row in block.rows() for entry in row] for block in basis]
    stacked = matrix(ZZ, rows)
    assert stacked.rank() == radical_rank * (c_rank - 1)
    # Saturated: a proper sublattice of the kernel would leave a Smith invariant above 1.
    invariants = [int(entry) for entry in stacked.smith_form()[0].diagonal() if int(entry) != 0]
    assert invariants == [1] * (radical_rank * (c_rank - 1)), invariants


def test_shear_kernel_of_a_vanishing_quotient_coordinate_is_all_of_hom() -> None:
    """For `x = 0` the condition is vacuous: the kernel is the whole `Hom(C,R)` of rank `c_rank*radical_rank`."""
    c_rank, radical_rank = 2, 2
    basis = shear_kernel_basis(c_rank, radical_rank, (0, 0))

    assert len(basis) == c_rank * radical_rank
    rows = [[int(entry) for row in block.rows() for entry in row] for block in basis]
    stacked = matrix(ZZ, rows)
    assert stacked.rank() == c_rank * radical_rank
    assert abs(int(stacked.det())) == 1, "the standard basis is unimodular, so this is all of Hom(C,R)"
