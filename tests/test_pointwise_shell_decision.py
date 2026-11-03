r"""Exact `Gamma`-equivalence on the whole isotropic cone of `U + U` (issues #53 and #62).

`Gamma = Stab_{O(L)}(e1)` **pointwise**. The decision is by classification -- for `m = b(e1, v) != 0` the class of
`(c, d)` in `(Z/m)^2` up to `m Z^2` and four symmetries; for `m = 0` the class of `(c, d)` plus `a mod gcd(c, d)`;
and in general the pair `(content, primitive class)`, since an isometry is a unimodular change of basis that is
`Z`-linear -- and every witness is an explicit element of `Gamma`, validated by the port's checked row-matrix
boundary.

The classification is cross-checked against the **retained generator set** rather than against itself: each
retained `Gamma` generator is applied to a sample of cone vectors, and the class key must agree. That would fail if
the key were too fine (which is the plausible error, since a too-coarse key would not be caught by it). The other
direction is checked as a law against the proved separating invariant.
"""
from itertools import product

import pytest

from sage_indefinite_port.indefinite.pointwise_shell import (
    U_PLUS_U_GRAM,
    gamma_action,
    gamma_class,
    gamma_decompose,
    gamma_element,
    gamma_equivalent,
    gamma_witness,
    isotropic_content,
    pointwise_shell_action,
    pointwise_shell_group,
    pointwise_shell_quotient,
    shell_coordinates,
)
from sage_indefinite_port.indefinite.recursive import EichlerRecursion
from sage_indefinite_port.readiness import UnfinishedCapability
from tests.acceptance.consumer import element, lattice

_RECURSION = EichlerRecursion()
_LATTICE = None


def _L():
    global _LATTICE
    if _LATTICE is None:
        _LATTICE = lattice([list(row) for row in U_PLUS_U_GRAM])
    return _LATTICE


def _e1():
    return element(_L(), (1, 0, 0, 0))


def _v(a, m, c, d):
    return element(_L(), (a, m, c, d))


def _cone_box(radius):
    r"""Every isotropic vector of the box `[-radius, radius]^4`, the zero vector included.

    Satisfying `a m + c d = 0` is the isotropic condition in this basis, so this enumerates the cone exactly rather
    than by testing `q` -- an independent statement of the same hypothesis, and it includes the imprimitive vectors.
    """
    return [
        _v(a, m, c, d)
        for a, m, c, d in product(range(-radius, radius + 1), repeat=4)
        if a * m + c * d == 0
    ]


def _shell_sample():
    """Primitive isotropic vectors: `a m + c d = 0` and content one."""
    sample = [_v(0, 1, 1, 0), _v(-6, 1, 2, 3), _v(0, 1, 3, 0), _v(3, 1, -3, 1),
              _v(0, 2, 1, 0), _v(0, 0, 1, 0), _v(1, 0, 2, 0), _v(2, 0, 3, 0),
              _v(1, 0, 0, 0), _v(-1, 0, 0, 0), _v(0, 0, 0, 1),
              # negative `m` and opposite vectors: `b(e1, .)` is signed, so these are not the negatives' classes
              _v(0, -1, 1, 0), _v(0, -2, 1, 0), _v(0, 1, 0, 0), _v(0, -1, 0, 0),
              _v(0, -1, -1, 0), _v(6, -2, 3, 4)]
    for vector in sample:
        assert int(vector.q()) == 0 and int(vector.content()) == 1, f"bad sample member {vector}"
    return sample


def test_the_cone_hypotheses_are_checked_not_assumed() -> None:
    """A non-isotropic vector is refused, an imprimitive one is not, and `a m + c d = 0` holds for every cone vector.

    Content is deliberately *not* a hypothesis: imprimitive isotropic vectors are classified by splitting content
    off, so refusing them here would hide a decision the module makes.
    """
    L = _L()
    assert shell_coordinates(L, _v(-6, 1, 2, 3)) == (-6, 1, 2, 3)
    with pytest.raises(UnfinishedCapability):
        shell_coordinates(L, element(L, (1, 0, 0, 0)) + element(L, (0, 1, 0, 0)))  # q = 2
    assert shell_coordinates(L, _v(0, 2, 2, 0)) == (0, 2, 2, 0), "content must not be a refusal"
    assert shell_coordinates(L, _v(0, 0, 0, 0)) == (0, 0, 0, 0)
    for vector in _cone_box(2):
        a, m, c, d = shell_coordinates(L, vector)
        assert a * m + c * d == 0
        assert int(vector.q()) == 0


def test_the_class_key_is_invariant_under_the_retained_gamma_generators() -> None:
    r"""Every retained generator of `Aut((L, e1))` preserves the class key on the whole sample.

    This is the cross-check that matters: the key is computed from the classification, not from the generators, so
    agreement on 6 generators against 11 vectors is independent evidence that the key is not too fine.
    """
    L = _L()
    generators = _RECURSION.pointwise_stabilizer_generators(_e1())
    assert generators
    for vector in _shell_sample():
        for generator in generators:
            image = generator(vector)
            assert int(image.q()) == 0 and int(image.content()) == 1, "a Gamma image must stay in the shell"
            assert gamma_class(L, image) == gamma_class(L, vector), (
                f"the class key split a Gamma-orbit: {vector} and its image {image}"
            )


def test_the_decision_is_an_equivalence_relation_and_agrees_with_witnesses() -> None:
    """Reflexive, symmetric on the sample, and a `None` witness exactly when the classes differ."""
    L = _L()
    sample = _shell_sample()
    for source in sample:
        assert gamma_equivalent(L, source, source)
        for target in sample:
            same = gamma_equivalent(L, source, target)
            assert same == gamma_equivalent(L, target, source)
            witness = gamma_witness(L, source, target)
            if same:
                assert witness is not None, f"class equality must produce a witness: {source} -> {target}"
                assert witness(source) == target, "the witness must carry the source onto the target exactly"
                assert witness(_e1()) == _e1(), "the witness must lie in Gamma, i.e. fix e1"
            else:
                assert witness is None


def test_the_family_v_m_is_separated_by_m_and_reaches_across_by_symmetry() -> None:
    r"""`v_m = e2 + m f1` for distinct `m` are inequivalent; within one `m` the symmetries act transitively on the
    residue classes they are allowed to."""
    L = _L()
    family = [element(L, (0, m, 1, 0)) for m in range(1, 5)]  # v_m, primitive isotropic with b(e1, v_m) = m
    for m, source in enumerate(family, start=1):
        assert shell_coordinates(L, source) == (0, m, 1, 0)
    for index, source in enumerate(family):
        for target in family[index + 1:]:
            assert not gamma_equivalent(L, source, target), "distinct m must be distinct Gamma-orbits"
            assert gamma_witness(L, source, target) is None

    # m = 2: (c,d) = (1,0) and (3,0) agree modulo 2, so they are equivalent and a witness exists
    equivalent_pair = (element(L, (0, 2, 1, 0)), element(L, (-6, 2, 3, 4)))
    assert shell_coordinates(L, equivalent_pair[1]) == (-6, 2, 3, 4)
    assert gamma_equivalent(L, *equivalent_pair)
    witness = gamma_witness(L, *equivalent_pair)
    assert witness is not None and witness(equivalent_pair[0]) == equivalent_pair[1]


def test_e1_and_minus_e1_are_different_orbits() -> None:
    r"""The pointwise stabilizer fixes `e1`, so `-e1` is a different class -- based vector, not unoriented line.

    In the classification this is the `m = 0`, `c = d = 0` case, where `a` itself is the invariant because
    `gcd(0, 0) = 0` admits no shift.
    """
    L = _L()
    assert gamma_class(L, _v(1, 0, 0, 0)) != gamma_class(L, _v(-1, 0, 0, 0))
    assert gamma_witness(L, _v(1, 0, 0, 0), _v(-1, 0, 0, 0)) is None
    # and the line stabilizer does carry one to the other, which is the distinction being preserved
    line = _RECURSION.line_stabilizer_generators(_e1())
    assert any(generator(_e1()) == -_e1() for generator in line)


def test_the_unipotent_part_moves_a_inside_the_shell() -> None:
    r"""For `m = 0` the unipotent `u(rho, gam)` shifts `a` by the ideal `gcd(c, d) Z`; that is what makes
    `e1 + 2 e2` and `e1 + 4 e2` inequivalent while `e2` and `e1 + e2` are equivalent."""
    L = _L()
    # gcd(1, 0) = 1, so a is unconstrained: e2 ~ e1 + e2
    assert gamma_equivalent(L, _v(0, 0, 1, 0), _v(1, 0, 1, 0))
    witness = gamma_witness(L, _v(0, 0, 1, 0), _v(1, 0, 1, 0))
    assert witness is not None and witness(_v(0, 0, 1, 0)) == _v(1, 0, 1, 0)

    # gcd(3, 0) = 3, so a moves by multiples of 3: a = 1 and a = 4 agree, a = 1 and a = 2 do not
    assert gamma_equivalent(L, _v(1, 0, 3, 0), _v(4, 0, 3, 0))
    assert not gamma_equivalent(L, _v(1, 0, 3, 0), _v(2, 0, 3, 0))
    witness = gamma_witness(L, _v(1, 0, 3, 0), _v(4, 0, 3, 0))
    assert witness is not None and witness(_v(1, 0, 3, 0)) == _v(4, 0, 3, 0)
    # and (3,0) is not (6,0) under the symmetries, so those differ whatever a is
    assert not gamma_equivalent(L, _v(1, 0, 3, 0), _v(1, 0, 6, 0))


def test_every_witness_is_a_genuine_gamma_element_of_the_retained_generators_group() -> None:
    r"""Each constructed witness fixes `e1` and preserves the pairing -- an isometry of `L`, checked directly."""
    L = _L()
    sample = _shell_sample()
    for source in sample:
        for target in sample:
            witness = gamma_witness(L, source, target)
            if witness is None:
                continue
            assert witness(_e1()) == _e1()
            for probe in (_v(1, 0, 0, 0), _v(0, 1, 0, 0), _v(0, 0, 1, 0), _v(2, -3, 1, 4)):
                assert int(witness(probe).q()) == int(probe.q()), "the witness must preserve norms"
                assert int(L.b(witness(probe), witness(source))) == int(L.b(probe, source)), (
                    "the witness must preserve the pairing"
                )


def test_the_classification_covers_the_whole_cone_and_only_the_cone() -> None:
    r"""Every isotropic vector gets a class -- imprimitive and zero included -- and a non-isotropic one is refused.

    The refusal is the honest boundary of the decision: `q` alone does not classify the non-isotropic vectors, so
    none of them is given a class it has no business having.
    """
    L = _L()
    for vector in _cone_box(2):
        assert gamma_class(L, vector) is not None
    assert gamma_class(L, _v(0, 0, 0, 0)) == ("linear", (0, 0), 0), "the zero vector is its own class"
    assert gamma_class(L, _v(0, 2, 2, 0)) == ("content", 2, gamma_class(L, _v(0, 1, 1, 0)))
    with pytest.raises(UnfinishedCapability):
        gamma_class(L, element(L, (1, 0, 0, 0)) + element(L, (0, 1, 0, 0)))  # q = 2, not isotropic


def test_content_splits_the_cone_into_scaled_primitive_orbits() -> None:
    r"""`t w` and `t' w'` are equivalent exactly when `t = t'` and `w ~ w'` -- the content corollary, executed.

    `content` is preserved because an isometry of `L` is a unimodular change of basis, and `g(t w) = t g(w)` because
    it is `Z`-linear. So the cone's orbits are the primitive orbits scaled by `t`, and the check is two-sided:
    same-content witnesses are certified at the scaled endpoints, and scaling is shown **not** to be an orbit map.
    """
    L = _L()
    sample = _shell_sample()
    scaled_by_content = {}
    for t in (1, 2, 3):
        scaled = [_v(*(t * entry for entry in shell_coordinates(L, w))) for w in sample]
        scaled_by_content[t] = scaled
        for w, tw in zip(sample, scaled, strict=True):
            assert int(tw.content()) == t, "the scaled sample must have content t"
            assert isotropic_content(shell_coordinates(L, tw)) == t
            expected_class = gamma_class(L, w) if t == 1 else ("content", t, gamma_class(L, w))
            assert gamma_class(L, tw) == expected_class
        for index, source in enumerate(scaled):
            for other in range(index + 1, len(scaled)):
                target = scaled[other]
                expected = gamma_equivalent(L, sample[index], sample[other])
                assert gamma_equivalent(L, source, target) == expected, (
                    f"scaling by {t} changed the orbit relation: {source} vs {target}"
                )
                witness = gamma_witness(L, source, target)
                if expected:
                    assert witness is not None and witness(source) == target
                    assert witness(_e1()) == _e1()
                else:
                    assert witness is None

    # different contents are never equivalent: `content` is a proved Gamma-invariant
    for t, scaled in scaled_by_content.items():
        for other_t, other_scaled in scaled_by_content.items():
            if other_t <= t:
                continue
            for source in scaled[:4]:
                for target in other_scaled[:4]:
                    assert not gamma_equivalent(L, source, target), (
                        f"content {t} and {other_t} must not be merged"
                    )
                    assert gamma_witness(L, source, target) is None


def test_the_cone_is_decided_exhaustively_on_a_box() -> None:
    r"""On every box of the cone: same class means a certified witness, and the proved invariant never splits a class.

    Same class is checked **constructively**, and that direction is not a formality: every element of `Gamma` is
    `u(rho, gam) . S` by the module's theorem, so if the key merged two vectors no such composite carries one onto
    the other and `_certified` raises rather than passing. So on this box, class equality implies an exhibited
    element of `Gamma` -- a key that is too coarse cannot survive the loop.

    Different class is the direction no search can settle, and the class key is its authority. What can be checked
    independently is the **proved** separating invariant as a law: whenever it separates a pair, the classification
    must separate that pair too. That catches a key that is too fine; the constructive loop above catches one that is
    too coarse.
    """
    L = _L()
    box = _cone_box(2)
    by_class = {}
    for vector in box:
        by_class.setdefault(gamma_class(L, vector), []).append(vector)

    witnessed = 0
    for members in by_class.values():
        for index, source in enumerate(members):
            for target in members[index + 1:]:
                witness = gamma_witness(L, source, target)
                assert witness is not None, f"same class but no witness: {source} -> {target}"
                assert witness(source) == target, "the witness must carry the source onto the target"
                assert witness(_e1()) == _e1(), "the witness must fix e1"
                witnessed += 1
    assert witnessed > 100, f"the box must exercise many same-class pairs, saw {witnessed}"

    marked = (_e1(),)
    for index, source in enumerate(box):
        for target in box[index + 1:]:
            if _RECURSION.pointwise_stabilizer_separates(marked, source, target):
                assert gamma_class(L, source) != gamma_class(L, target), (
                    f"the proved invariant separates {source} and {target} but the classification merged them"
                )


def test_retained_generators_preserve_content_and_class_on_the_imprimitive_cone() -> None:
    r"""Content and class survive the retained generators on imprimitive vectors, not only on primitive ones.

    This is the too-fine direction of the same cross-check, carried past primitivity: an isometry cannot change the
    gcd of a vector's coordinates, and it cannot merge content classes.
    """
    L = _L()
    generators = _RECURSION.pointwise_stabilizer_generators(_e1())
    assert generators
    imprimitive = [vector for vector in _cone_box(2) if int(vector.content()) > 1]
    assert len(imprimitive) == 32, f"the radius-two box has 32 imprimitive cone vectors, got {len(imprimitive)}"
    for vector in imprimitive:
        for generator in generators:
            image = generator(vector)
            assert int(image.content()) == int(vector.content()), "content is an isometry invariant"
            assert int(image.q()) == 0, "an isometry maps the cone to itself"
            assert gamma_class(L, image) == gamma_class(L, vector), (
                f"the class key split a Gamma-orbit: {vector} and its image {image}"
            )


def test_the_gamma_set_reaches_the_imprimitive_cone() -> None:
    r"""The widened G-set: `in_same_orbit` is the exact decision on imprimitive vectors and on the zero vector.

    The point set is now the whole cone `q = 0` rather than the primitive slice, and the relation is the same exact
    one, driven through the provider's retained action constructor.
    """
    L = _L()
    action = pointwise_shell_action(L)
    pairs = [
        (_v(0, 2, 2, 0), _v(0, 2, 6, 0)),   # 2*(0,1,1,0) and 2*(0,1,3,0): m = 1, so equivalent
        (_v(0, 2, 2, 0), _v(0, 3, 3, 0)),   # content 2 against content 3
        (_v(0, 0, 0, 0), _v(0, 0, 0, 0)),   # the zero vector is a fixed point
        (_v(0, 0, 0, 0), _v(0, 2, 2, 0)),
        (_v(-12, 2, 4, 6), _v(-6, 1, 2, 3)),  # 2*(-6,1,2,3) against (-6,1,2,3): content 2 against 1
    ]
    for left, right in pairs:
        assert action.in_same_orbit(left, right) == gamma_equivalent(L, left, right)
    assert action.orbits().orbit_of(_v(0, 2, 2, 0)) is not None
    assert pointwise_shell_quotient(L).is_finite() is not True, (
        "the quotient is infinite and must not be reported finite"
    )


def test_the_signed_pairing_defect_is_fixed() -> None:
    r"""`f1` and `-f1` are primitive isotropic, and they are **not** `Gamma`-equivalent.

    The published revision keyed the affine class on `abs(m)` and lost `sign b(e1, v)`, which every pointwise
    stabilizer preserves. `f1 = (0,1,0,0)` and `-f1 = (0,-1,0,0)` therefore shared a key, `gamma_equivalent`
    returned `True`, and `gamma_witness` returned the identity -- which does not carry one onto the other. Both
    halves are pinned here.
    """
    L = _L()
    source, target = _v(0, 1, 0, 0), _v(0, -1, 0, 0)

    assert int(L.b(_e1(), source)) == 1 and int(L.b(_e1(), target)) == -1
    assert gamma_class(L, source) != gamma_class(L, target), "signed m must be part of the class key"
    assert not gamma_equivalent(L, source, target)
    assert gamma_witness(L, source, target) is None

    # and every returned witness must actually reach its target, which the identity did not
    for candidate in (_v(0, 1, 1, 0), _v(0, -1, 1, 0), _v(-6, 1, 2, 3)):
        for other in _shell_sample():
            witness = gamma_witness(L, candidate, other)
            if witness is not None:
                assert witness(candidate) == other, "a successful witness must carry source onto target"
                assert witness(_e1()) == _e1()


def test_the_invariant_is_never_finer_than_the_classification() -> None:
    r"""Whenever the retained invariant separates two vectors, the classification must separate them too.

    The invariant `(q, content, b(e1, .))` is proved `Gamma`-invariant, so a class key that merges a pair the
    invariant separates is wrong by construction -- which is exactly the `abs(m)` defect. This law is checked over
    the whole sample in both directions of the pair, so it would have caught the defect and will catch its
    relatives.
    """
    L = _L()
    marked = (_e1(),)
    sample = _shell_sample()
    for index, source in enumerate(sample):
        for target in sample[index + 1:]:
            if _RECURSION.pointwise_stabilizer_separates(marked, source, target):
                assert not gamma_equivalent(L, source, target), (
                    f"the proved invariant separates {source} and {target} but the classification merged them"
                )


def test_negative_m_classes_behave_like_positive_ones() -> None:
    r"""`m < 0` follows the same rule with `m` kept signed: `(Z/|m|)^2` up to symmetries, `a = -cd/m`."""
    L = _L()
    # m = -2, (c,d) = (1,0) and (3,4) agree modulo 2 and are equivalent, with a witnessed transport
    source, target = _v(0, -2, 1, 0), _v(6, -2, 3, 4)
    assert shell_coordinates(L, source) == (0, -2, 1, 0)
    assert shell_coordinates(L, target) == (6, -2, 3, 4)
    assert gamma_equivalent(L, source, target)
    witness = gamma_witness(L, source, target)
    assert witness is not None and witness(source) == target and witness(_e1()) == _e1()

    # the same (c,d) data with the opposite sign of m is a different class
    mirrored = _v(0, 2, 1, 0)
    assert not gamma_equivalent(L, source, mirrored)
    assert gamma_witness(L, source, mirrored) is None

    # and within one sign the symmetries act: (1,0) ~ (-1,0)
    assert gamma_equivalent(L, _v(0, -1, 1, 0), _v(0, -1, -1, 0))


def test_opposite_vectors_are_never_equivalent_unless_the_class_says_so() -> None:
    r"""`-v` need not lie in the `Gamma`-orbit of `v`; each case is decided, not assumed."""
    L = _L()
    decided = []
    for vector in _shell_sample():
        negated = element(L, tuple(-entry for entry in shell_coordinates(L, vector)))
        decided.append((vector, gamma_equivalent(L, vector, negated)))
    # `-e1` is fixed by no element of Gamma, and neither is `-f1` reachable from `f1`; but `-` does preserve some
    # classes, so the point is that each is decided rather than that all are the same
    assert any(equivalent for _, equivalent in decided)
    assert any(not equivalent for _, equivalent in decided)

    for vector in _shell_sample():
        negated = element(L, tuple(-entry for entry in shell_coordinates(L, vector)))
        witness = gamma_witness(L, vector, negated)
        if witness is not None:
            assert witness(vector) == negated
            assert witness(_e1()) == _e1()


def test_the_proved_generating_set_decomposes_exactly_and_generates_gamma() -> None:
    r"""`Gamma = <sigma, tau, u(1,0), u(0,1)>`, and each generator's decomposition is decided, not sampled.

    The four elements decompose to `(1,0,0)`, `(2,0,0)`, `(0,1,0)` and `(0,0,1)`, the structural statement that the
    generating set is the `N`-part times the `Theta`-part. Every generator returned by the retained
    `pointed_automorphism_generators` also decomposes, which decides its membership in `Gamma` exactly instead of
    comparing a few images.
    """
    L = _L()
    assert gamma_decompose(L, gamma_action(L, "sigma")) == (1, 0, 0)
    assert gamma_decompose(L, gamma_action(L, "tau")) == (2, 0, 0)
    assert gamma_decompose(L, gamma_action(L, "u", 1, 0)) == (0, 1, 0)
    assert gamma_decompose(L, gamma_action(L, "u", 0, 1)) == (0, 0, 1)
    # composition lands where the semidirect product says
    assert gamma_decompose(L, gamma_action(L, "u", 1, 1)) == (0, 1, 1)

    retained = _RECURSION.pointwise_stabilizer_generators(_e1())
    assert retained
    for generator in retained:
        assert gamma_decompose(L, generator) is not None, (
            "a retained generator failed the exact membership test in Gamma"
        )

    assert len(tuple(pointwise_shell_group(L).group_generators())) == 4


def test_an_isometry_outside_gamma_does_not_decompose() -> None:
    r"""The decomposition decides membership: an isometry moving `e1` returns `None`."""
    from sage.all import ZZ, matrix

    from sage_indefinite_port.indefinite.recursive import _isometry_from_row_matrix

    L = _L()
    not_gamma = _isometry_from_row_matrix(
        L, L, matrix(ZZ, [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    )
    assert not_gamma(_e1()) != _e1()
    assert gamma_decompose(L, not_gamma) is None


def test_the_retained_action_constructor_carries_the_pointwise_shell() -> None:
    r"""The G-set is the provider's, driven by `Gamma`, with the exact relation as its orbit relation.

    Reusing `GObjects(...).on_set(...)` rather than a second framework: `in_same_orbit` agrees with
    `gamma_equivalent` on the shell, `orbit_of` returns the provider's orbit object, and the quotient's
    `is_finite()` stays a proposition, because the quotient here is infinite.
    """
    L = _L()
    action = pointwise_shell_action(L)

    pairs = [(_v(0, 2, 1, 0), _v(-6, 2, 3, 4)), (_v(0, 2, 1, 0), _v(0, 3, 1, 0)),
             (_v(0, 1, 0, 0), _v(0, -1, 0, 0)), (_v(1, 0, 3, 0), _v(4, 0, 3, 0))]
    for left, right in pairs:
        assert action.in_same_orbit(left, right) == gamma_equivalent(L, left, right)

    assert action.orbits().orbit_of(_v(0, 2, 1, 0)) is not None
    assert pointwise_shell_quotient(L).is_finite() is not True, (
        "the quotient is infinite and must not be reported finite"
    )


def test_mixed_symmetry_unipotent_round_trips_across_all_four_symmetries() -> None:
    r"""Every `u(rho, gam) . S` decomposes back to `(S, rho, gam)` -- for all four `S` and mixed `(rho, gam)`.

    The published `gamma_decompose` read `rho` off the `e1` coefficient of `g(e2)`, which is only correct for
    `S = identity`; for `S = sigma` that coefficient is `-rho`, so `u(1,0) sigma` -- an element that fixes `e1` --
    came back as `None`. The extraction now solves `A . (rho, gam)^T = (e1-coefficient of g(e2), ... of g(f2))`,
    where the rows of `A` are `S(1,0)` and `S(0,1)`.
    """
    L = _L()
    cases = [(0, 0), (1, 0), (0, 1), (1, 1), (-1, 0), (0, -1), (2, -3), (-2, 5)]
    for index in range(4):
        for rho, gam in cases:
            element_ = gamma_element(L, index, rho, gam)
            assert element_(_e1()) == _e1(), "a Gamma element must fix e1"
            assert gamma_decompose(L, element_) == (index, rho, gam), (
                f"round trip failed for S index {index}, rho {rho}, gam {gam}"
            )


def test_the_reported_mixed_composition_case_is_decided() -> None:
    r"""`u(1,0) . sigma`: `g(e2) = -e2 - e1`, `g(f2) = -f2`, and the decomposition is `(1, 1, 0)`.

    This is the exact false negative reported against the previous revision: the element fixes `e1` and was
    reported outside `Gamma`.
    """
    L = _L()
    labels = tuple(L.module_generating_set())
    element_ = gamma_element(L, 1, 1, 0)

    assert element_(_e1()) == _e1()

    def image_coordinates(position):
        image = element_(L.module_generator(labels[position])).to_vector()
        return [int(image(label)) for label in labels]

    assert image_coordinates(2) == [-1, 0, -1, 0], "g(e2) must be -e2 - e1"
    assert image_coordinates(3) == [0, 0, 0, -1], "g(f2) must be -f2"
    assert gamma_decompose(L, element_) == (1, 1, 0)


def test_decompose_is_total_exactly_on_the_isometries_fixing_e1() -> None:
    r"""`Gamma` *is* `{g : g(e1) = e1}`, so `None` means precisely "moves `e1`" -- and nothing else.

    An earlier attempt to build an isometry that fixes `e1` yet has a non-symmetry `Q`-part was abandoned: `O(Q)` has
    exactly four elements, so no such isometry exists, and that is exactly why the decomposition is total on `Gamma`
    rather than a partial extraction. The refusals are therefore the elements that move `e1`.
    """
    from sage.all import ZZ, matrix

    from sage_indefinite_port.indefinite.recursive import _isometry_from_row_matrix

    L = _L()
    labels = tuple(L.module_generating_set())
    moving = _isometry_from_row_matrix(
        L, L, matrix(ZZ, [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    )
    assert moving(_e1()) != _e1()
    assert gamma_decompose(L, moving) is None

    # every mixed element fixes e1 and decomposes; that is the other half of the same statement
    for index in range(4):
        candidate = gamma_element(L, index, 1, -2)
        assert candidate(_e1()) == _e1()
        assert gamma_decompose(L, candidate) is not None
