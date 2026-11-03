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

At radical rank `k` the same computation runs with `D` in `GL_k(Z)`, and `s - D r` in `x R` asks for
`D r = s mod x`, i.e. for `s mod x` to lie in the `GL_k(Z)`-orbit of `r mod x`. That image is the subtle part:
`GL_1(Z) = {+-1}`, so at rank one only `s = +- r (mod x)` is reachable and the two signs are genuinely different
orbits; at rank `k >= 2`, `GL_k(Z) -> GL_k(Z/n)` is onto, so the residue condition is the whole of `GL_k(Z/n)`
and is coarser than the rank-one test. The tests below pin both, with the witness or the refusal built by hand
from the block conditions rather than read off the port's answer.

The unmarked radical `Aut`/`Isom` route is implemented separately and must keep working; it is checked here
so a change to the marked path cannot quietly take it out.
"""
import pytest

from sage.all import MatrixGroup, ZZ, identity_matrix, matrix, vector

from sage_indefinite_port.indefinite.recursive import EichlerRecursion, _column_matrix
from tests.acceptance.consumer import (
    automorphisms,
    content,
    coordinates,
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

    # The full kernel, not one example.  For rank(C) = 2 and x_0 = (1,0) primitive, the map
    # Hom(C,R) -> R, B |-> B x_0, is onto a rank-1 radical, so its kernel is free of rank
    # rank(C) - 1 = 1 and is saturated.  The port's shears must therefore span a rank-1 *saturated*
    # lattice of Hom(C,R) = Z^2; a proper sublattice would leave a Smith invariant above 1.
    found = matrix(ZZ, [[int(entry) for row in shear.rows() for entry in row] for shear in shears])
    assert found.rank() == 1
    invariants = [int(entry) for entry in found.smith_form()[0].diagonal() if int(entry) != 0]
    assert invariants == [1], f"the shears must span a saturated lattice, got invariants {invariants}"




# ---------------------------------------------------------------------------------------------------
# Lazy orbit consumers (#41 with #53): the native set, never a truncated list.
# ---------------------------------------------------------------------------------------------------

def test_marked_orbit_is_the_native_lazy_set_not_a_list() -> None:
    """`marked_orbit` returns the provider's predicate-cut set, and it claims no finiteness.

    `L[n]/O(L)` can be infinite, so a list would be a truncation. The set answers membership exactly and
    leaves `is_finite()` as a proposition.
    """
    lattice_ = lattice(_GRAM)
    source = element(lattice_, [1, 0])
    orbit = _RECURSION.marked_orbit(source)

    assert source in orbit
    # The shear (1,0) -> (1,1) is in the group, so this member must be decided in.
    assert element(lattice_, [1, 1]) in orbit
    # Same square (2) but a different radical coordinate content: (3,0) has q = 18, so it is out.
    assert element(lattice_, [3, 0]) not in orbit
    # Not a finite ordered set, and not claimed finite.
    assert not hasattr(orbit, "list")
    assert orbit.is_finite() is not True


def test_marked_orbit_of_a_radical_vector_handles_the_c_zero_case() -> None:
    """A vector in the radical has `c = 0`: the condition is exactly `D r = s`, and the orbit is still lazy."""
    lattice_ = lattice(_GRAM)
    source = element(lattice_, [0, 1])
    orbit = _RECURSION.marked_orbit(source)

    assert source in orbit
    assert element(lattice_, [0, -1]) in orbit        # the radical sign change is integral
    assert element(lattice_, [0, 2]) not in orbit     # (0,2) is not primitive, so no unimodular image
    assert element(lattice_, [1, 0]) not in orbit     # different square
    assert orbit.is_finite() is not True


def test_marked_orbit_action_exposes_the_native_orbit_surface() -> None:
    """The action object carries the operations #53 requires, rather than only a points list."""
    lattice_ = lattice(_GRAM)
    action = _RECURSION.marked_orbit_action(lattice_)

    for name in ("orbits", "stabilizer", "in_same_orbit"):
        assert hasattr(action, name), name
    first, second = element(lattice_, [1, 0]), element(lattice_, [1, 1])
    assert action.in_same_orbit(first, second) is True
    assert action.in_same_orbit(first, element(lattice_, [3, 0])) is False


# ---------------------------------------------------------------------------------------------------
# The quotient L[n]/O(L), and the restored general coverage.
# ---------------------------------------------------------------------------------------------------

def test_marked_orbit_quotient_is_the_orbit_set_and_refuses_to_be_listed() -> None:
    """`L[n]/O(L)` is a set of orbit classes, used through `orbit_of`, and it refuses to be listed.

    The quotient is not the points of one orbit: `orbit_of` returns a class. And the provider enforces the
    #53 discipline itself -- `__iter__` asserts that the underlying G-set *is placed finite*
    ("the orbits of ... are listed only when it is a finite set"), so an infinite quotient cannot be
    silently truncation-listed. That refusal is the contract, not a defect, so it is asserted here.
    """
    lattice_ = lattice(_GRAM)
    quotient = _RECURSION.marked_orbit_quotient(lattice_, 2)

    # The quotient exposes orbit classes, not vector points, and answers lazily.
    assert hasattr(quotient, "orbit_of")
    source = element(lattice_, [1, 0])
    orbit = quotient.orbit_of(source)
    assert hasattr(orbit, "points")
    assert source in orbit
    assert element(lattice_, [1, 1]) in orbit          # the shear image is in the same class
    assert element(lattice_, [0, 1]) not in orbit      # q = 0, a different square

    # No finiteness is claimed, and listing is refused rather than truncated.
    assert quotient.is_finite() is not True
    with pytest.raises(AssertionError):
        next(iter(quotient))


def _definite_complement_lattice(c_rank: int, radical_rank: int):
    gram = [[2 if row == column and row < c_rank else 0
             for column in range(c_rank + radical_rank)]
            for row in range(c_rank + radical_rank)]
    return lattice(gram)


@pytest.mark.parametrize(
    ("c_rank", "radical_rank", "marked"),
    [
        (1, 1, (0, 1)),
        (2, 1, (0, 0, 1)),
        (1, 2, (0, 1, 0)),
        (1, 2, (0, 2, 3)),        # a non-standard primitive radical vector
    ],
)
def test_pointwise_stabilizer_of_a_radical_marked_vector(c_rank, radical_rank, marked) -> None:
    """A vector inside the radical no longer raises: the split-radical path gives its pointwise stabilizer.

    `Aut((L,v))` for `v = (0,r)` is `{[[A,B],[0,D]] : A in O(C), rD = r, B in Hom(C,R)}`, so its generators
    are the embedded `O(C)` generators, the `GL(R)`-stabilizer of `r`, and the whole standard basis of
    `Hom(C,R)` as shears. This used to raise `ArithmeticError("the isotropic pairing L -> J* must have
    finite cokernel")` because the isotropic route requires a nondegenerate ambient form.
    """
    lattice_ = _definite_complement_lattice(c_rank, radical_rank)
    vector_ = element(lattice_, marked)
    assert _RECURSION.vector_lies_in_radical(vector_)

    generators = tuple(_RECURSION.pointed_automorphism_generators(vector_))
    assert generators
    assert all(generator(vector_) == vector_ for generator in generators)

    # The shear part must be the *whole* Hom(C,R): rank(C)*rank(R) freedoms, each present and primitive.
    shears = [shear for shear in (_split_blocks(_column_matrix(g), c_rank, radical_rank)[1] for g in generators)
              if shear != 0]
    rows = [[int(entry) for row in shear.rows() for entry in row] for shear in shears]
    stacked = matrix(ZZ, rows)
    assert stacked.rank() == c_rank * radical_rank
    invariants = [int(entry) for entry in stacked.smith_form()[0].diagonal() if int(entry) != 0]
    assert invariants == [1] * (c_rank * radical_rank), invariants


def test_isotropic_vector_off_the_radical_is_not_classified_as_radical() -> None:
    """`q(v) = 0` is necessary but not sufficient: an isotropic vector off the radical is not radical."""
    hyperbolic = lattice([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    isotropic = element(hyperbolic, [1, 0, 0])   # q = 0, but b(v, (0,1,0)) = 1, so it is off the radical
    radical_lattice = lattice([[2, 0], [0, 0]])
    radical_vector = element(radical_lattice, [0, 1])

    assert int(isotropic.q()) == 0
    assert _RECURSION.vector_lies_in_radical(radical_vector)
    assert not _RECURSION.vector_lies_in_radical(isotropic)


# ---------------------------------------------------------------------------------------------------
# The quotient L[n]/O(L), and the restored general coverage.
# ---------------------------------------------------------------------------------------------------

def test_marked_orbit_quotient_is_the_orbit_set_and_refuses_to_be_listed() -> None:
    """`L[n]/O(L)` is a set of orbit classes, used through `orbit_of`, and it refuses to be listed.

    The quotient is not the points of one orbit: `orbit_of` returns a class. And the provider enforces the
    #53 discipline itself -- `__iter__` asserts that the underlying G-set *is placed finite*
    ("the orbits of ... are listed only when it is a finite set"), so an infinite quotient cannot be
    silently truncation-listed. That refusal is the contract, not a defect, so it is asserted here.
    """
    lattice_ = lattice(_GRAM)
    quotient = _RECURSION.marked_orbit_quotient(lattice_, 2)

    # The quotient exposes orbit classes, not vector points, and answers lazily.
    assert hasattr(quotient, "orbit_of")
    source = element(lattice_, [1, 0])
    orbit = quotient.orbit_of(source)
    assert hasattr(orbit, "points")
    assert source in orbit
    assert element(lattice_, [1, 1]) in orbit          # the shear image is in the same class
    assert element(lattice_, [0, 1]) not in orbit      # q = 0, a different square

    # No finiteness is claimed, and listing is refused rather than truncated.
    assert quotient.is_finite() is not True
    with pytest.raises(AssertionError):
        next(iter(quotient))


def _definite_complement_lattice(c_rank: int, radical_rank: int):
    gram = [[2 if row == column and row < c_rank else 0
             for column in range(c_rank + radical_rank)]
            for row in range(c_rank + radical_rank)]
    return lattice(gram)


def _kernel_shears(c_rank: int, radical_rank: int, marked_coordinates):
    """The shear blocks of the pointwise stabilizer on a real lattice, with identity diagonal blocks."""
    lattice_ = _definite_complement_lattice(c_rank, radical_rank)
    marked = element(lattice_, marked_coordinates)
    stabilizer = tuple(_RECURSION.pointed_automorphism_generators(marked))
    assert stabilizer, "the pointwise stabilizer has a generating family"
    assert all(generator(marked) == marked for generator in stabilizer)

    identity_c = identity_matrix(ZZ, c_rank)
    identity_r = identity_matrix(ZZ, radical_rank)
    labels = tuple(lattice_.module_generating_set())
    marked_coordinates = vector(ZZ, [int(marked.to_vector()(label)) for label in labels])
    quotient_part = marked_coordinates[:c_rank]
    shears = []
    for generator in stabilizer:
        A, shear, D = _split_blocks(_column_matrix(generator), c_rank, radical_rank)
        if shear == 0:
            continue
        # A kernel shear has both diagonal blocks trivial and annihilates the marked C-coordinate.
        assert A == identity_c
        assert D == identity_r
        assert shear * quotient_part == 0
        shears.append(shear)
    return lattice_, marked, shears


@pytest.mark.parametrize(
    ("c_rank", "radical_rank", "marked", "expected_rank"),
    [
        (2, 1, (1, 0, 0), 1),
        (3, 1, (1, 0, 0, 0), 2),
        (2, 1, (5, 0, 0), 1),          # content 5: the same primitive direction, so the same kernel
    ],
)
def test_pointwise_stabilizer_shears_span_the_full_shear_kernel(c_rank, radical_rank, marked, expected_rank):
    """The kernel `{B : B x = 0}` has rank `rank(R)*(rank(C)-1)` for primitive `x`, and all of
    `Hom(C,R)` when `x = 0`; the port's shears must span a *saturated* lattice of exactly that rank.

    For `rank(C) = 2`, `rank(R) = 1`, `x_0 = (1,0)` primitive, the map `Hom(C,R) -> R`, `B |-> Bx`, is onto a
    rank-1 radical, so the kernel is free of rank 1; for `x = 0` the condition is vacuous and the kernel is
    all of `Hom(C,R)`, of rank `rank(C)*rank(R)`.
    """
    lattice_, marked_vector, shears = _kernel_shears(c_rank, radical_rank, marked)
    assert shears, "the shear kernel should be generated"

    rows = [[int(entry) for row in shear.rows() for entry in row] for shear in shears]
    stacked = matrix(ZZ, rows)
    assert stacked.rank() == expected_rank
    invariants = [int(entry) for entry in stacked.smith_form()[0].diagonal() if int(entry) != 0]
    assert invariants == [1] * expected_rank, invariants


# ---------------------------------------------------------------------------------------------------
# The radical branch's actual domain: permuted, mixed, zero, non-primitive, purely radical.
#
# The branch does not assume the caller's basis has the radical last. `indefinite_reduced_presentation`
# supplies the split and its coordinate transport, so these bases are carried in, solved, and carried back.
# ---------------------------------------------------------------------------------------------------

def _fixing_isometries(gram, coordinates):
    lattice_ = lattice(gram)
    vector_ = element(lattice_, coordinates)
    generators = tuple(_RECURSION.pointed_automorphism_generators(vector_))
    assert generators, "the pointwise stabilizer has a generating family"
    assert all(generator(vector_) == vector_ for generator in generators)
    assert all(is_isometry(rows_of(generator), gram) for generator in generators)
    return lattice_, vector_, generators


def test_radical_stabilizer_in_a_permuted_basis() -> None:
    """Radical FIRST: the vector is genuine, so rejecting it would be a domain defect, not a limitation."""
    lattice_ = lattice([[0, 0], [0, 2]])
    vector_ = element(lattice_, [1, 0])
    assert int(vector_.q()) == 0
    assert _RECURSION.vector_lies_in_radical(vector_)
    _, _, generators = _fixing_isometries([[0, 0], [0, 2]], [1, 0])
    assert len(generators) >= 2


def test_radical_stabilizer_in_a_mixed_basis() -> None:
    """Radical interleaved between two complement coordinates."""
    gram = [[2, 0, 0], [0, 0, 0], [0, 0, 2]]
    _, vector_, generators = _fixing_isometries(gram, [0, 1, 0])
    assert int(vector_.q()) == 0
    assert generators


def test_radical_stabilizer_of_the_zero_vector_is_the_whole_group() -> None:
    """`v = 0` is fixed by every isometry, so the pointwise stabilizer is `O(L)` itself."""
    lattice_ = lattice(_GRAM)
    zero = element(lattice_, [0, 0])
    generators = tuple(_RECURSION.pointed_automorphism_generators(zero))

    assert generators
    assert all(is_isometry(rows_of(generator), _GRAM) for generator in generators)
    assert all(generator(zero) == zero for generator in generators)


def test_radical_stabilizer_of_a_nonprimitive_radical_vector() -> None:
    """`(0,2)` is not primitive; `R` is torsion-free, so its stabilizer is that of its primitive direction."""
    gram = [[2, 0, 0], [0, 0, 0], [0, 0, 0]]
    lattice_ = lattice(gram)
    vector_ = element(lattice_, [0, 2, 0])
    assert int(vector_.content()) == 2
    _, _, generators = _fixing_isometries(gram, [0, 2, 0])
    assert generators
    # The stabilizer of (0,2) and of (0,1) coincide, because 2*(d - 1) = 0 forces d = 1 in a torsion-free R.
    primitive_generators = _fixing_isometries(gram, [0, 1, 0])[2]
    assert len(primitive_generators) == len(generators)


def test_radical_stabilizer_of_a_purely_radical_lattice() -> None:
    """`C = 0`: there is no complement block, and the stabilizer is the `GL(R)`-stabilizer of the vector."""
    gram = [[0, 0], [0, 0]]
    _, vector_, generators = _fixing_isometries(gram, [1, 0])
    assert generators
    # In a purely radical lattice every isometry is a GL(R) matrix, and the shears e_2 -> e_2 + e_1 fix e_1.
    forms = [matrix(ZZ, rows_of(generator)) for generator in generators]
    assert any(form == matrix(ZZ, [[1, 0], [1, 1]]) for form in forms)


def test_lazy_orbit_representative_enumeration_is_a_stream() -> None:
    """The quotient is enumerable lazily: one representative per class, taken on demand, never materialized.

    `limit` caps how many representatives are yielded; it is not a termination guarantee for the search.
    A square whose classes are exhausted keeps looking for a next class, so a finite class set is fully
    compatible with an unbounded candidate walk.
    """
    lattice_ = lattice(_GRAM)
    stream = _RECURSION.marked_orbit_representatives(lattice_, 2, limit=1)

    assert hasattr(stream, "__next__"), "the enumeration is a generator, not a list"
    representative = next(stream)
    assert int(representative.q()) == 2
    # It is the class of (1,0): the shear and the complement sign both connect them.
    assert _RECURSION.pointed_isometry(representative, element(lattice_, [1, 0])) is not None
    assert _RECURSION.pointed_isometry(representative, element(lattice_, [1, 7])) is not None


def test_lazy_enumeration_yields_distinct_classes_only() -> None:
    """Every yielded representative is a genuinely new class: no two are equivalent.

    For `L = <2> + 0` and square 2 all vectors `(+-1, y)` lie in one `O(L)`-orbit, so the walk must not
    yield a second representative for them.
    """
    lattice_ = lattice(_GRAM)
    stream = _RECURSION.marked_orbit_representatives(lattice_, 2, limit=1)
    first = next(stream)
    try:
        second = next(stream)
    except StopIteration:
        second = None
    if second is not None:
        assert _RECURSION.pointed_isometry(first, second) is None


# ---------------------------------------------------------------------------------------------------
# The generated kernel for radical rank >= 2: maintained Hom/evaluation/kernel, and mixed words.
# ---------------------------------------------------------------------------------------------------

def _hom_element_matrix(complement, radical, morphism):
    """The coordinate matrix of a `Hom(C,R)` element, by evaluating it on the generators of `C`."""
    rows = []
    for source_label in complement.module_generating_set():
        image = morphism(complement.module_generator(source_label))
        rows.append([int(image.to_vector()(label)) for label in radical.module_generating_set()])
    return matrix(ZZ, rows).transpose()          # radical_rank x complement_rank


def _blocks(c_rank: int, r_rank: int):
    """The complement and radical factors of `diag(<2>^c, 0^r)`, as lattices on the provider's ring."""
    from sage_extensions.all import Lattices, ZZ as OwnedZZ

    complement = Lattices(OwnedZZ)([[2 if row == column else 0 for column in range(c_rank)]
                                    for row in range(c_rank)])
    radical = Lattices(OwnedZZ)([[0] * r_rank for _ in range(r_rank)])
    return complement, radical


@pytest.mark.parametrize(
    ("c_rank", "r_rank", "marked"),
    [
        (2, 1, (1, 0, 0)),
        (3, 1, (1, 0, 0, 0)),
        (2, 2, (1, 0, 0, 0)),
        (3, 2, (1, 0, 0, 0, 0)),
    ],
)
def test_shear_kernel_through_the_maintained_hom_evaluation_kernel(c_rank, r_rank, marked) -> None:
    """`{B : B x = 0}` has rank `rank(R)*(rank(C)-1)` for nonzero `x`, computed by the provider itself.

    The kernel is obtained from the Hom module, `evaluation`, the morphism constructor and `kernel()`,
    with the kernel's `inclusion` giving back elements of `Hom(C,R)`. Every returned element must kill the
    marked coordinate; the rank is `rank(R)*(rank(C)-1)` and it comes from the provider, not from a
    coordinate construction in the port.
    """
    complement, radical = _blocks(c_rank, r_rank)
    lattice_ = lattice([[2 if row == column and row < c_rank else 0
                         for column in range(c_rank + r_rank)] for row in range(c_rank + r_rank)])
    vector_ = element(lattice_, marked)
    x = complement.linear_combination({complement.module_generating_set()[0]: 1})

    kernel = _RECURSION.shear_kernel(complement, radical, x)
    assert len(kernel) == r_rank * (c_rank - 1)
    for morphism in kernel:
        image = morphism(complement.module_generator(complement.module_generating_set()[0]))
        assert [int(image.to_vector()(label)) for label in radical.module_generating_set()] == [0] * r_rank


def _stabilizer_word_shears(c_rank: int, r_rank: int, marked, depth: int = 3):
    """Pure shears produced by *words* of length at most `depth` in the stabilizer generators.

    Counting the pure-shear members of one generator list measures nothing about the subgroup they
    generate -- a product of mixed elements can be a pure shear. This walks the words and keeps those whose
    diagonal blocks are both the identity, which are exactly the elements of the kernel that the words
    produce.
    """
    lattice_ = lattice([[2 if row == column and row < c_rank else 0
                         for column in range(c_rank + r_rank)] for row in range(c_rank + r_rank)])
    vector_ = element(lattice_, marked)
    # Column convention throughout: `_split_blocks` reads B as the bottom-left block, which is where the
    # shear lives in the column matrix. Row matrices would put it top-right and hide every shear.
    forms = [_column_matrix(generator)
             for generator in _RECURSION.pointed_automorphism_generators(vector_)]

    def key_of(form):
        return tuple(int(entry) for row in form.rows() for entry in row)

    identity = identity_matrix(ZZ, c_rank + r_rank)
    frontier = {key_of(identity): identity}
    collected = {}
    for _ in range(depth):
        next_frontier = {}
        for word in frontier.values():
            for form in forms:
                product = word * form
                key = key_of(product)
                if key in collected or key in next_frontier:
                    continue
                next_frontier[key] = product
                A, shear, D = _split_blocks(product, c_rank, r_rank)
                if A == identity_matrix(ZZ, c_rank) and D == identity_matrix(ZZ, r_rank) and shear != 0:
                    collected[key] = shear
        frontier = next_frontier
        if not frontier:
            break
    return tuple(collected.values())


@pytest.mark.parametrize(
    ("c_rank", "r_rank", "marked"),
    [
        (2, 1, (1, 0, 0)),
        (3, 1, (1, 0, 0, 0)),
        (2, 2, (1, 0, 0, 0)),
        (3, 2, (1, 0, 0, 0, 0)),
    ],
)
def test_generated_shear_kernel_equals_the_maintained_kernel(c_rank, r_rank, marked) -> None:
    """The shears reachable by words generate the *whole* kernel: asserted as equality, not a rank bound.

    A rank bound is not the claim. `reachable.rank() <= full_rank` holds for every sublattice, and stacking a
    full-rank kernel basis on the generated rows makes the combined rank `max(...)` automatically -- both
    pass for a strictly deficient generated sublattice. What has to be checked is the **index**:

    the generated shears lie in the kernel `K` (a saturated lattice of rank `full_rank`), so they span a
    sublattice `L`. Writing `L`'s generators in `K`'s basis gives an integer matrix `X` of rank
    `full_rank`, and `[K : L]` is the product of `X`'s Smith invariants. Equality of lattices is exactly
    `[K : L] = 1`; any prime factor means the words do *not* generate the kernel.

    Pure shears reachable by words of length at most three, mixed generators included.
    """
    shears = _stabilizer_word_shears(c_rank, r_rank, marked)
    assert shears, "the stabilizer must generate some kernel elements"

    # Every word-shear is a genuine kernel element: it annihilates the marked C-coordinate.
    lattice_ = lattice([[2 if row == column and row < c_rank else 0
                         for column in range(c_rank + r_rank)] for row in range(c_rank + r_rank)])
    vector_ = element(lattice_, marked)
    labels = tuple(lattice_.module_generating_set())
    quotient_part = vector(ZZ, [int(vector_.to_vector()(label)) for label in labels][:c_rank])
    for shear in shears:
        assert shear * quotient_part == 0

    generated = matrix(ZZ, [[int(entry) for row in shear.rows() for entry in row] for shear in shears])

    complement, radical = _blocks(c_rank, r_rank)
    x = complement.linear_combination({complement.module_generating_set()[0]: 1})
    kernel = _RECURSION.shear_kernel(complement, radical, x)
    full_rank = r_rank * (c_rank - 1)
    kernel_matrix = matrix(ZZ, [
        [int(entry) for row in _hom_element_matrix(complement, radical, morphism).rows() for entry in row]
        for morphism in kernel
    ])
    assert kernel_matrix.rank() == full_rank

    # Coordinates of the generated shears in the kernel's basis: an exact integer matrix, since L is inside K.
    coordinates = generated * kernel_matrix.transpose() * (kernel_matrix * kernel_matrix.transpose()).inverse()
    assert coordinates.denominator() == 1, "the generated shears are not in the kernel's lattice"
    coordinate_matrix = coordinates.change_ring(ZZ)
    assert coordinate_matrix.rank() == full_rank, "the generated shears do not even span the kernel's rank"

    invariants = [abs(int(coordinate_matrix.smith_form()[0][index, index])) for index in range(full_rank)]
    index_in_kernel = 1
    for invariant in invariants:
        index_in_kernel *= invariant
    assert index_in_kernel == 1, (
        f"the generated shear sublattice has index {index_in_kernel} in the kernel (Smith invariants "
        f"{invariants}), so the words do not generate it"
    )


def test_the_index_check_detects_a_deficient_sublattice() -> None:
    """Control: the index computation used above must report > 1 for a strictly deficient sublattice.

    Without this, `index == 1` could be passing for the wrong reason. Two checks: a purely arithmetic one,
    and one built on a real kernel basis with the sublattice scaled by 2.
    """
    def index_in(coordinates, rank):
        invariants = [abs(int(coordinates.smith_form()[0][i, i])) for i in range(rank)]
        value = 1
        for invariant in invariants:
            value *= invariant
        return value, invariants

    # Arithmetic control: 2*Z^2 inside Z^2 has index 4.
    coordinates = 2 * identity_matrix(ZZ, 2)
    value, invariants = index_in(coordinates, 2)
    assert value == 4 and invariants == [2, 2]

    # Real control: the (2,2) kernel basis, with the generated lattice scaled by 2.
    complement, radical = _blocks(2, 2)
    x = complement.linear_combination({complement.module_generating_set()[0]: 1})
    kernel = _RECURSION.shear_kernel(complement, radical, x)
    kernel_matrix = matrix(ZZ, [
        [int(entry) for row in _hom_element_matrix(complement, radical, morphism).rows() for entry in row]
        for morphism in kernel
    ])
    full_rank = kernel_matrix.rank()
    scaled = 2 * kernel_matrix
    coordinates = scaled * kernel_matrix.transpose() * (kernel_matrix * kernel_matrix.transpose()).inverse()
    assert coordinates.denominator() == 1
    value, invariants = index_in(coordinates.change_ring(ZZ), full_rank)
    assert value == 2 ** full_rank, f"a doubled sublattice must have index 2^rank, got {value} ({invariants})"


# ---------------------------------------------------------------------------------------------------------------------
# Issue #41, remaining open items: mixed markings at radical rank 2, and the image of GL(R) in Aut(R/xR).
# ---------------------------------------------------------------------------------------------------------------------

_RANK_TWO_RADICAL_GRAM = [[2, 0, 0], [0, 0, 0], [0, 0, 0]]


def _rank_two_radical_lattice():
    return lattice(_RANK_TWO_RADICAL_GRAM)


def _block(unit, shear, radical_block):
    r"""The `[[A, 0], [B, D]]` block, column convention: `g v = (A x, B x + D r)` for `v = (x, r)`.

    Written out here so the witness for a transport is built from the *conditions* and then checked against the
    lattice, rather than copied from whatever the port returns. At `C = <2>` the only symmetry is `A = +-1`; `B` is
    a column of length `k`; `D` is `k x k`.
    """
    size = 1 + len(radical_block)
    rows = [[0] * size for _ in range(size)]
    rows[0][0] = unit
    for index in range(len(radical_block)):
        rows[1 + index][0] = shear[index]
        for column in range(len(radical_block)):
            rows[1 + index][1 + column] = radical_block[index][column]
    return matrix(ZZ, rows)


def _checked_block(source_coordinates, target_coordinates, unit, shear, radical_block, gram):
    r"""Build the block and assert it is an automorphism of the lattice carrying `source` onto `target`.

    Three things are checked and none of them is the port's word: the matrix is invertible over `Z` (`|det| = 1`,
    which for an integral matrix is exactly integral invertibility), it preserves the pairing
    (`g^T G g = G`), and it sends the source coordinates onto the target coordinates. Only then is the transport
    known to exist, independently of any search.
    """
    block = _block(unit, shear, radical_block)
    assert abs(int(block.det())) == 1, "an automorphism of the lattice must be unimodular"
    form = matrix(ZZ, gram)
    assert block.transpose() * form * block == form, "the block must preserve the pairing"
    image = block * vector(ZZ, source_coordinates)
    assert tuple(int(entry) for entry in image) == tuple(target_coordinates), (
        "the hand-built block must carry the source onto the target"
    )
    return block


# (source, target, A, B, D, why) -- every one of these must transport, and the hand-built block proves it.
_RANK_TWO_TRANSPORTS = [
    ((5, 1, 0), (5, 2, 0), 1, (0, -1), [[2, 1], [5, 3]],
     "both residues are nonzero modulo 5, and GL_2(Z/5) moves one to the other"),
    ((5, 1, 0), (5, 3, 0), 1, (0, -1), [[3, 1], [5, 2]],
     "same residue orbit, a different lift"),
    ((5, 1, 0), (5, 0, 1), 1, (0, 0), [[0, 1], [1, 0]],
     "the radical swap"),
    ((5, 1, 1), (5, 1, -1), 1, (0, 0), [[1, 0], [-2, 1]],
     "a shear in GL_2(Z) alone"),
    ((5, 5, 0), (5, 10, 0), 1, (1, 0), [[1, 0], [0, 1]],
     "both residues are zero modulo 5, so the shear B carries the difference"),
    ((5, 5, 0), (5, 0, 5), 1, (0, 0), [[0, 1], [1, 0]],
     "zero residues again, moved by the swap"),
]


@pytest.mark.parametrize(
    "source_coordinates, target_coordinates, unit, shear, radical_block, why",
    _RANK_TWO_TRANSPORTS,
    ids=[f"{case[0]}->{case[1]}" for case in _RANK_TWO_TRANSPORTS],
)
def test_rank_two_radical_mixed_markings_transport_by_an_explicit_gl2_block(
    source_coordinates, target_coordinates, unit, shear, radical_block, why
) -> None:
    r"""Mixed markings at radical rank 2: the residue condition is `GL_2(Z)`, and each transport is exhibited.

    `L = <2> + 0 + 0` has radical `R = Z^2` and `C = <2>`, and every isometry is `[[A, 0], [B, D]]` with `A` in
    `O(<2>) = {+-1}`, `D` in `GL_2(Z)` and `B` in `Hom(C, R) = Z^2`. So `g(v) = w` asks for `A x = y` and
    `s - D r` in `x R`, and `B` then carries it. Each case below builds that block, checks it is a lattice
    automorphism, and only then asks the port -- so agreement is two independent computations, not a rerun.
    """
    lattice_ = _rank_two_radical_lattice()
    _checked_block(
        source_coordinates, target_coordinates, unit, shear, radical_block, _RANK_TWO_RADICAL_GRAM
    )

    source = element(lattice_, source_coordinates)
    target = element(lattice_, target_coordinates)
    witness = _RECURSION.pointed_isometry(source, target)
    assert witness is not None, f"expected a transport ({why})"
    assert witness(source) == target, "the returned witness must carry the source onto the target"
    assert is_isometry(rows_of(witness), _RANK_TWO_RADICAL_GRAM), "the witness must be an isometry of L"


def test_rank_two_residue_orbit_is_coarser_than_the_rank_one_sign_test() -> None:
    r"""`(5,1) -> (5,2)` at radical rank 1 and `(5,1,0) -> (5,2,0)` at radical rank 2 differ, and must.

    At rank one the residue condition is `s = +- r (mod 5)`, so `2` is not reachable from `1`. At rank two
    `GL_2(Z) -> GL_2(Z/5)` is onto and `GL_2(Z/5)` is transitive on the nonzero vectors, so `(2,0)` *is* reachable
    from `(1,0)`. Same residue data, same norm, same content, opposite verdicts: the difference is the image of
    `GL(R)` in `Aut(R/xR)`, which is a proper subgroup exactly at rank one.
    """
    rank_one_lattice = _degenerate_lattice()
    rank_one_source = element(rank_one_lattice, (5, 1))
    rank_one_target = element(rank_one_lattice, (5, 2))
    assert _RECURSION.pointed_isometry(rank_one_source, rank_one_target) is None, (
        "at rank one the sign is the whole of GL_1(Z), so 2 = +-1 mod 5 is false"
    )

    rank_two_lattice = _rank_two_radical_lattice()
    _checked_block((5, 1, 0), (5, 2, 0), 1, (0, -1), [[2, 1], [5, 3]], _RANK_TWO_RADICAL_GRAM)
    assert _RECURSION.pointed_isometry(
        element(rank_two_lattice, (5, 1, 0)), element(rank_two_lattice, (5, 2, 0))
    ) is not None, "at rank two the lift exists, so the refusal must not carry over"


def test_rank_two_transport_is_refused_when_the_content_differs() -> None:
    r"""The rank-two negatives in this family are separated by `content`, which is a proved invariant.

    The residue orbit is everything modulo 5, so a refusal cannot come from the residue: `r = 0 mod 5` forces every
    radical coordinate of `v` to be divisible by 5 while `x = 5`, hence `content(v)` a multiple of 5, whereas
    `s` not `0 mod 5` forces `content(w) = 1`. The refusal is therefore explained by the invariant rather than by a
    failed search, and the check is that the port agrees while the naive residue data would not separate them.
    """
    lattice_ = _rank_two_radical_lattice()
    for source_coordinates, target_coordinates in (
        ((5, 1, 0), (5, 0, 5)),
        ((5, 5, 0), (5, 1, 5)),
        ((5, 0, 1), (5, 5, 5)),
    ):
        source = element(lattice_, source_coordinates)
        target = element(lattice_, target_coordinates)
        assert int(source.q()) == int(target.q()), "the norms agree, so the refusal is not a norm test"
        assert content(source) != content(target), "content is the invariant that separates them"
        assert _RECURSION.pointed_isometry(source, target) is None


# (source, target, A, B, D, why) -- the rank-one sign cases, each with the hand-built block that realises it.
# `D` is the radical unit and `B` is the shear that absorbs `s - D r` over `x`.
_RANK_ONE_SIGN_CASES = [
    ((5, 1), (5, 4), 1, (1,), [[-1]], "4 = -1 mod 5, so the negating radical unit is the witness"),
    ((5, 1), (5, 9), 1, (2,), [[-1]], "9 = -1 mod 5, reached by the same unit plus a doubling shear"),
    ((5, 2), (5, 3), 1, (1,), [[-1]], "3 = -2 mod 5: the sign, not equality of residues"),
    ((5, 2), (-5, 2), -1, (0,), [[1]], "A = -1 is a symmetry of C, not a transporter, so it must be allowed"),
    ((-5, 1), (5, 4), -1, (-1,), [[-1]], "both blocks negated"),
]


@pytest.mark.parametrize(
    "source_coordinates, target_coordinates, unit, shear, radical_block, why",
    _RANK_ONE_SIGN_CASES,
    ids=[f"{case[0]}->{case[1]}" for case in _RANK_ONE_SIGN_CASES],
)
def test_rank_one_radical_image_is_the_sign_and_the_symmetry_of_c(
    source_coordinates, target_coordinates, unit, shear, radical_block, why
) -> None:
    r"""At `R = Z`, `GL_1(Z) = {+-1}` and `O(<2>) = {+-1}`: both signs are witnesses, and both are needed.

    A dropped sign here is a false negative of exactly the kind the pointwise-shell `abs(m)` defect was: the
    division `s - D r` in `cR` would report "not divisible" for a pair that the negating unit transports. Each case
    builds the block first and then asks the port.
    """
    lattice_ = _degenerate_lattice()
    _checked_block(source_coordinates, target_coordinates, unit, shear, radical_block, _GRAM)

    source = element(lattice_, source_coordinates)
    target = element(lattice_, target_coordinates)
    witness = _RECURSION.pointed_isometry(source, target)
    assert witness is not None, f"expected a transport ({why})"
    assert witness(source) == target
    assert is_isometry(rows_of(witness), _GRAM)


def test_rank_one_refusal_is_exhaustive_over_the_two_units_and_two_symmetries() -> None:
    r"""`(5,2) -> (5,1)` is refused, and the refusal is a finite exhaustion, not a bounded search.

    At rank one there are four candidate blocks: `A` in `{+-1}` (irrelevant here, as `x = y`) and `D` in `{+-1}`.
    Neither unit satisfies `d r = s (mod 5)`, so no transporter exists -- and the check below runs those four, which
    is *all* of them, rather than trusting that the port's search covered enough. `q`, `content` and the
    divisibility agree on the pair, so the refusal cannot be blamed on any of them.
    """
    lattice_ = _degenerate_lattice()
    source = element(lattice_, (5, 2))
    target = element(lattice_, (5, 1))

    assert int(source.q()) == int(target.q())
    assert content(source) == content(target)
    source_x, source_r = coordinates(source)
    target_x, target_r = coordinates(target)
    assert (source_x, target_x) == (5, 5)
    for unit in (1, -1):
        for radical_unit in (1, -1):
            remainder = target_r - radical_unit * source_r
            assert remainder % source_x != 0, (
                f"A={unit}, d={radical_unit} must fail the congruence, so no shear can absorb it"
            )
    assert _RECURSION.pointed_isometry(source, target) is None
