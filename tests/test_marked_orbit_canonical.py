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


def _degenerate():
    return lattice(_DEGENERATE)


_CANONICAL_CACHE: dict[tuple, MarkedOrbitCanonical] = {}


def _canonical(lattice_, coordinates) -> MarkedOrbitCanonical:
    """`marked_orbit_canonical`, memoized: each call scans a coordinate box and each `pointed_isometry` is costly."""
    gram = _gram_in_label_order(lattice_)
    key = (tuple(tuple(int(gram[row, column]) for column in range(gram.ncols())) for row in range(gram.nrows())),
           tuple(coordinates))
    if key not in _CANONICAL_CACHE:
        _CANONICAL_CACHE[key] = _RECURSION.marked_orbit_canonical(element(lattice_, coordinates))
    return _CANONICAL_CACHE[key]


def _coordinates(lattice_, vector) -> tuple[int, ...]:
    values = vector.to_vector()
    return tuple(int(values(label)) for label in lattice_.module_generating_set())


def _earlier_equivalent_exists(lattice_, vector) -> bool:
    """Whether some lattice element strictly before `vector` in the declared order is equivalent to it."""
    target_height = max((abs(entry) for entry in _coordinates(lattice_, vector)), default=0)
    for height in range(target_height + 1):
        for coordinates in _coordinate_box(lattice_, height):
            candidate = _coordinate_element(lattice_, coordinates)
            if candidate == vector:
                return False
            if _RECURSION.pointed_isometry(candidate, vector) is not None:
                return True
    return False


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
    assert not _earlier_equivalent_exists(L, canonical.representative)


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


def test_primitive_radical_vectors_are_decided_and_non_primitive_ones_are_refused() -> None:
    r"""The undecided subdomain is exactly the non-primitive part of the radical, and it is named.

    Measured on `<2> + 0`, whose radical is the `y`-axis: `(0, 1)` and `(0, -1)` are one class and the transport
    is decided, while `(0, 2)` reaches `primitive_isotropic_subobject`, which asserts that the span of `(2 e_1)`
    has torsion. The enumeration must not survive that by reading it as "different orbit" -- it raises, naming the
    gap -- and it must not refuse the primitive part either, which *is* decided.

    This is also why the content family of issue #62's test 4 (`c e_1`, `c = 1, 2, ...`) is not available through
    this route: those classes become radical vectors exactly when `q` vanishes on them.
    """
    L = _degenerate()

    assert _RECURSION.vector_lies_in_radical(element(L, (0, 1)))
    assert _RECURSION.vector_lies_in_radical(element(L, (0, 2)))
    assert not _RECURSION.vector_lies_in_radical(element(L, (5, 1)))

    # the primitive radical class is decided on both sides
    canonical = _canonical(L, (0, 1))
    assert _coordinates(L, canonical.representative) == (0, -1)
    assert canonical.transporter(canonical.representative) == element(L, (0, 1))

    with pytest.raises(UnfinishedCapability):
        _RECURSION.marked_orbit_canonical(element(L, (0, 2)))
    # the enumeration keeps going to height 1 and refuses when it reaches the non-primitive class
    assert [tuple(_coordinates(L, vector)) for vector in
            _RECURSION.marked_orbit_representatives(L, 0, height_bound=1)] == [(0, 0), (0, -1)]
    with pytest.raises(UnfinishedCapability):
        list(_RECURSION.marked_orbit_representatives(L, 0, height_bound=2))


def test_the_quotient_is_unchanged_under_a_permuted_presentation() -> None:
    r"""Canonicalization is basis-relative; the orbit relation it induces is not.

    The representative's *coordinates* depend on the declared order, which is why no functorial claim is made.
    What must not depend on it is the equivalence itself: two vectors are in one orbit in the standard
    presentation exactly when their coordinate-swapped images are in one orbit in the permuted one.
    """
    standard = _degenerate()
    permuted = lattice([[0, 0], [0, 2]])  # the same lattice, radical first

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


def test_canonicalizer_on_the_hyperbolic_plane() -> None:
    r"""`L = U` with `q(x, y) = 2xy`: square 4 is one orbit `{(+-1,+-2),(+-2,+-1)}` with minimum `(-2, -1)`."""
    L = lattice(_HYPERBOLIC)
    canonical = _canonical(L, (2, 1))

    assert _coordinates(L, canonical.representative) == (-2, -1)
    assert canonical.transporter(canonical.representative) == element(L, (2, 1))
    assert list(_RECURSION.marked_orbit_representatives(L, 4, limit=1)) == [canonical.representative]
    # square 2 holds only `(+-1, +-1)`, whose minimum is `(-1, -1)`
    assert _coordinates(L, list(_RECURSION.marked_orbit_representatives(L, 2, limit=1))[0]) == (-1, -1)
