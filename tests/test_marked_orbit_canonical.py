"""Canonical representatives and duplicate-free orbit enumeration (issue #62, parts B and C).

The domain is the one the port already supports: `O(L)` acting on the vectors of an integral lattice, with
`pointed_isometry` as the equivalence backend -- exact and terminating on the `#41` domain, degenerate `L`
included. The order is the declared one: height `max |coordinate|` first, then lexicographic within one height.

Two small specimens carry the interesting behaviour, and both come from the port's own `#41` record:

* `L = <2> + 0`, the degenerate rank-two lattice with `q(x, y) = 2x^2`, radical last. In the reduced
  presentation `O(L)` is `[[A, 0], [B, D]]` with `A, D = ±1` and `B` arbitrary, so `(x, y) |-> (±x, Bx ± y)`.
  On the square-50 slice that forces `x = ±5` and makes `y mod 5` the class invariant, giving **exactly three
  orbits**, with canonical representatives derivable by hand: `(-5, -5)` (y = 0 mod 5), `(-5, -4)` (y = ±1) and
  `(-5, -3)` (y = ±2). Note `(5, 1)` and `(5, 6)` are in the same orbit, and `(5, 2)` is in another one although
  norm, content and divisibility all agree with `(5, 1)` -- the `#41` finding a canonicalizer must respect.
* the square-0 slice of the same lattice: `q(0, c) = 0` for every `c` and `(0, c) |-> (0, ±c)`, so `|c|` is an
  invariant and the quotient is **infinite**. The enumeration is exhaustive and duplicate-free with no
  termination signal -- issue #62's central distinction.

`-v` always lies in the orbit of `v`, so the minimum of an orbit is frequently a negative vector; the expected
values below are derived from the orbit description above, not guessed from the input.
"""
import pytest

from sage_indefinite_port.indefinite.recursive import (
    EichlerRecursion,
    MarkedOrbitCanonical,
    MarkedOrbitDecision,
    _coordinate_box,
    _coordinate_element,
)
from sage_indefinite_port.indefinite.groupoid_cache import _gram_in_label_order
from sage_indefinite_port.readiness import UnfinishedCapability
from tests.acceptance.consumer import element, lattice

_RECURSION = EichlerRecursion()
_DEGENERATE = [[2, 0], [0, 0]]
_HYPERBOLIC = [[0, 1], [1, 0]]


_LATTICES: dict[tuple, object] = {}


def _presentation(gram):
    """A lattice object per presentation, shared across tests: the reduction is cached per object."""
    key = tuple(tuple(int(entry) for entry in row) for row in gram)
    if key not in _LATTICES:
        _LATTICES[key] = lattice(gram)
    return _LATTICES[key]


def _degenerate():
    return _presentation(_DEGENERATE)


_CANONICAL_CACHE: dict[tuple, MarkedOrbitCanonical] = {}


def _canonical(lattice_, coordinates) -> MarkedOrbitCanonical:
    """`marked_orbit_canonical`, memoized: each call scans a coordinate box and each `pointed_isometry` is costly."""
    gram = _gram_in_label_order(lattice_)
    key = (tuple(tuple(int(gram[row, column]) for column in range(gram.ncols())) for row in range(gram.nrows())),
           tuple(coordinates))
    if key not in _CANONICAL_CACHE:
        _CANONICAL_CACHE[key] = _RECURSION.marked_orbit_canonical(element(lattice_, coordinates))
    return _CANONICAL_CACHE[key]


def _column_det(isometry) -> int:
    """The determinant of an isometry, in the port's column convention."""
    from sage_indefinite_port.indefinite.recursive import _column_matrix

    return int(_column_matrix(isometry).det())


def _coordinates(lattice_, vector) -> tuple[int, ...]:
    values = vector.to_vector()
    return tuple(int(values(label)) for label in lattice_.module_generating_set())


def test_canonical_representative_is_the_least_member_of_its_orbit() -> None:
    r"""`(5, 1) ~ (5, 6)`, and the minimum of that orbit is `(-5, -4)`: derived, then checked.

    With `x = ±5` forced and `y |-> Bx ± y`, the members are `(±5, y)` for `y = ±1 mod 5`; the least in the
    declared order is `x = -5` and then the smallest admissible `y`, which is `-4`.
    """
    L = _degenerate()
    v = element(L, (5, 6))

    canonical = _canonical(L, (5, 6))

    assert isinstance(canonical, MarkedOrbitCanonical)
    assert _coordinates(L, canonical.representative) == (-5, -4)
    assert canonical.height == 5
    # The transporter is a real isometry carrying the representative onto the input, not a relabelling.
    assert canonical.transporter(canonical.representative) == v


def test_canonical_representatives_distinguish_equal_invariant_vectors() -> None:
    r"""`(5, 1)` and `(5, 2)`: equal norm, content and divisibility, different orbits, different canonicals."""
    L = _degenerate()
    first = _canonical(L, (5, 1))
    second = _canonical(L, (5, 2))

    assert _coordinates(L, first.representative) == (-5, -4)
    assert _coordinates(L, second.representative) == (-5, -3)
    assert first.representative != second.representative
    assert _RECURSION.pointed_isometry(first.representative, second.representative) is None


def test_canonicalization_is_stable_and_idempotent() -> None:
    """`x ~ x'` implies equal canonicals; canonicalizing a canonical changes nothing."""
    L = _degenerate()
    representative = _canonical(L, (5, 1)).representative

    assert _canonical(L, (5, 6)).representative == representative
    assert _RECURSION.marked_orbit_canonical(representative).representative == representative
    assert _RECURSION.marked_orbit_is_canonical(representative)
    assert not _RECURSION.marked_orbit_is_canonical(element(L, (5, 6)))


def test_enumeration_of_a_finite_slice_is_exact_and_duplicate_free() -> None:
    r"""Square 50 of `<2> + 0` has exactly three orbits, with the derived canonical representatives.

    Three is provable: `q = 50` forces `x = ±5`, and `y mod 5` is the class invariant, so the classes are
    `y = 0`, `y = ±1` and `y = ±2`. The enumeration must emit exactly one representative per class, each of them
    canonical, and each absorbed candidate must carry a witness.
    """
    L = _degenerate()
    # A bounded scan (issue #62 part D): the search stops at height 5, so the absorbed candidates that follow
    # the third emission are reachable instead of the generator waiting for a fourth class that does not exist.
    decisions = list(_RECURSION.marked_orbit_enumeration(L, 50, height_bound=5))
    emitted = [decision for decision in decisions if decision.emitted]

    assert [_coordinates(L, decision.candidate) for decision in emitted] == [(-5, -5), (-5, -4), (-5, -3)]
    for decision in emitted:
        assert isinstance(decision, MarkedOrbitDecision)
        assert decision.representative == decision.candidate
        assert decision.transporter(decision.representative) == decision.candidate
        assert _canonical(L, _coordinates(L, decision.candidate)).representative == decision.candidate

    absorbed = [decision for decision in decisions if not decision.emitted]
    assert absorbed, "the slice contains absorbed candidates, so the witness path is exercised"
    emitted_set = {decision.candidate for decision in emitted}
    for decision in absorbed:
        # the witness carries the retained representative exactly onto the absorbed candidate
        assert decision.transporter(decision.representative) == decision.candidate
        assert decision.representative in emitted_set
    # ... and for a sample of them, the retained representative is also the candidate's own canonical one, so the
    # enumerator's bookkeeping agrees with the independent canonicalizer
    for decision in absorbed[:3]:
        assert _canonical(L, _coordinates(L, decision.candidate)).representative == decision.representative

    # distinct emissions are inequivalent, which is the duplicate-free property stated abstractly
    representatives = [decision.candidate for decision in emitted]
    for index, earlier in enumerate(representatives):
        for later in representatives[index + 1:]:
            assert _RECURSION.pointed_isometry(earlier, later) is None


def test_the_stream_is_prefix_stable_and_reports_no_completion() -> None:
    r"""Asked for fewer representatives, the same stream is a prefix of the longer one.

    That is what "an enumerator, not a count" means in issue #62's sense: the consumer decides when to stop, and
    nothing in the protocol says the classes are exhausted. Square 50 of `<2> + 0` has exactly three, which is
    provable from `x = ±5` and the invariant `y mod 5`; the stream cannot know that and does not claim it.
    """
    L = _degenerate()
    short = list(_RECURSION.marked_orbit_representatives(L, 50, limit=2))
    long = list(_RECURSION.marked_orbit_representatives(L, 50, limit=3))

    assert [tuple(_coordinates(L, vector)) for vector in short] == [(-5, -5), (-5, -4)]
    assert long[:2] == short
    assert len(long) == 3


def test_radical_transport_is_decided_for_every_content_with_an_integral_witness() -> None:
    r"""A radical marking is decided by its content, and the witness is an integral invertible isometry.

    In the reduced presentation `O(L)` is `[[A, 0], [B, D]]` with `D` unimodular, so a radical vector `(0, r)`
    goes to `(0, r D)` and the class is the content `gcd(r)`. `(0, 2)` and `(0, -2)` are therefore equivalent --
    the case that previously reached `primitive_isotropic_subobject` and asserted -- while `(0, 2)` and `(0, 3)`
    are not, and neither is `(0, 1)` against `(0, 2)`.
    """
    L = _degenerate()

    for pair in (((0, 2), (0, -2)), ((0, 1), (0, -1)), ((0, 6), (0, -6))):
        source, target = element(L, pair[0]), element(L, pair[1])
        witness = _RECURSION.pointed_isometry(source, target)
        assert witness is not None, f"{pair[0]} and {pair[1]} have the same content and must be transported"
        assert witness(source) == target, "the witness must carry the marking onto its image"
        assert abs(int(_column_det(witness))) == 1, "the witness must be an invertible integral isometry"

    for pair in (((0, 1), (0, 2)), ((0, 2), (0, 3)), ((5, 1), (0, 1)), ((0, 2), (0, 0))):
        source, target = element(L, pair[0]), element(L, pair[1])
        assert _RECURSION.pointed_isometry(source, target) is None, f"{pair[0]} and {pair[1]} differ in content"


def test_the_zero_vector_is_its_own_radical_class() -> None:
    """Content 0 is the zero vector, which every isometry fixes, and it is not conflated with a nonzero class."""
    L = _degenerate()
    zero = element(L, (0, 0))

    witness = _RECURSION.pointed_isometry(zero, zero)
    assert witness is not None and witness(zero) == zero
    assert _RECURSION.pointed_isometry(zero, element(L, (0, 2))) is None
    assert _coordinates(L, _canonical(L, (0, 2)).representative) == (0, -2)


def test_the_unbounded_content_family_is_a_constructive_infinite_quotient() -> None:
    r"""Square 0 of `<2> + 0` has one class per `|c|`, so the quotient is infinite -- and now decidable.

    `q(0, c) = 0` and `(0, c) |-> (0, +-c)`, so `|c|` separates the classes; the least member of the class of `c`
    is `(0, -|c|)` because `-c` precedes `c` in the declared order. Nothing reports that this quotient is
    finished, and nothing may: it is not finite. What is evidenced here is that the stream is prefix-stable and
    keeps producing classes past any fixed bound, each with an exact transporter to its canonical representative.
    """
    L = _degenerate()
    representatives = list(_RECURSION.marked_orbit_representatives(L, 0, limit=4))
    assert [_coordinates(L, vector) for vector in representatives] == [(0, 0), (0, -1), (0, -2), (0, -3)]

    # prefix stability: the same stream, asked for fewer
    assert list(_RECURSION.marked_orbit_representatives(L, 0, limit=2)) == representatives[:2]

    # each emitted class is distinct, and each transports to a non-canonical member of its class
    for earlier, later in zip(representatives, representatives[1:]):
        assert _RECURSION.pointed_isometry(earlier, later) is None
    for vector in representatives[1:3]:
        positive = element(L, (0, -_coordinates(L, vector)[1]))
        canonical = _canonical(L, _coordinates(L, positive))
        assert canonical.transporter(canonical.representative) == positive

    # a class beyond the ones yielded is already decidable, which is what "unbounded content" means here: the
    # number is not bounded by anything the enumerator knows, only by how far a test chooses to walk
    far = element(L, (0, 7))
    assert _RECURSION.pointed_isometry(far, element(L, (0, -7))) is not None
    assert _coordinates(L, _canonical(L, (0, 7)).representative) == (0, -7)


def test_radical_transport_in_a_permuted_and_mixed_basis() -> None:
    r"""The decision goes through the retained reduction, so the radical need not be the last coordinates.

    Three presentations of the same lattice: the radical last, the radical first, and a mixed basis in which no
    coordinate is purely radical. In each, content decides the class and the witness transports.
    """
    for gram in ([[2, 0], [0, 0]], [[0, 0], [0, 2]], [[2, 2], [2, 2]]):
        L = _presentation(gram)
        # in the mixed presentation `(1, -1)` spans the radical and `(1, 1)` does not
        radical = (0, 1) if gram == [[2, 0], [0, 0]] else ((1, 0) if gram == [[0, 0], [0, 2]] else (1, -1))
        doubled = tuple(2 * entry for entry in radical)
        negated = tuple(-entry for entry in doubled)

        assert _RECURSION.vector_lies_in_radical(element(L, doubled)), f"expected a radical vector in {gram}"
        witness = _RECURSION.pointed_isometry(element(L, doubled), element(L, negated))
        assert witness is not None, f"content does not depend on the presentation ({gram})"
        assert witness(element(L, doubled)) == element(L, negated)
        assert abs(int(_column_det(witness))) == 1


def test_the_quotient_is_unchanged_under_a_permuted_presentation() -> None:
    r"""Canonicalization is basis-relative; the orbit relation it induces is not.

    The representative's *coordinates* depend on the declared order, which is why no functorial claim is made.
    What must not depend on it is the equivalence itself: two vectors are in one orbit in the standard
    presentation exactly when their coordinate-swapped images are in one orbit in the permuted one.
    """
    standard = _degenerate()
    permuted = _presentation([[0, 0], [0, 2]])  # the same lattice, radical first

    canonical = _canonical(permuted, (6, 5))
    assert _RECURSION.marked_orbit_is_canonical(canonical.representative)
    assert canonical.transporter(canonical.representative) == element(permuted, (6, 5))

    assert _RECURSION.pointed_isometry(element(standard, (5, 1)), element(standard, (5, 6))) is not None
    assert _RECURSION.pointed_isometry(element(permuted, (1, 5)), element(permuted, (6, 5))) is not None
    assert _RECURSION.pointed_isometry(element(standard, (5, 1)), element(standard, (5, 2))) is None
    assert _RECURSION.pointed_isometry(element(permuted, (1, 5)), element(permuted, (2, 5))) is None


def test_an_undecided_equivalence_is_never_read_as_inequivalence(monkeypatch) -> None:
    r"""A backend that cannot decide must propagate, never be read as "distinct".

    Issue #62 is explicit: turning "unknown" into inequivalence would emit one orbit twice. Both the
    canonicalizer and the enumerator let `UnfinishedCapability` escape instead.
    """
    L = _degenerate()

    def refuse(*_args, **_kwargs):
        raise UnfinishedCapability("the equivalence backend is unfinished: work unit #54 owns it")

    monkeypatch.setattr(EichlerRecursion, "pointed_isometry", refuse)

    with pytest.raises(UnfinishedCapability):
        list(_RECURSION.marked_orbit_representatives(L, 0, limit=2))
    with pytest.raises(UnfinishedCapability):
        _RECURSION.marked_orbit_canonical(element(L, (0, 1)))


def test_canonicalizer_recognizes_an_already_canonical_vector() -> None:
    r"""On a non-degenerate lattice the minimum is recognised through the identity path.

    `L = U` with `q(x, y) = 2xy`: square 4 is one orbit `{(+-1,+-2),(+-2,+-1)}`, whose least member is `(-2, -1)`.
    Canonicalizing that member must return it unchanged, which the search does on its first comparison. Transport
    *between* distinct members of this orbit on `U` is the costly non-integral correction branch (measured at 57 s
    for `(2, 1) -> (-2, -1)`), a property of that branch and not of canonicalization, so this test does not
    exercise it.
    """
    L = _presentation(_HYPERBOLIC)
    canonical = _canonical(L, (-2, -1))

    assert _coordinates(L, canonical.representative) == (-2, -1)
    assert canonical.transporter(canonical.representative) == element(L, (-2, -1))
    assert _RECURSION.marked_orbit_is_canonical(element(L, (-2, -1)))
