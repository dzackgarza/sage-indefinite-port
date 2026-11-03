r"""Exact `Gamma`-equivalence on the primitive isotropic shell of `U + U` (issues #53 and #62).

`Gamma = Stab_{O(L)}(e1)` **pointwise**. The decision is by classification -- for `m = b(e1, v) != 0` the class of
`(c, d)` in `(Z/m)^2` up to `m Z^2` and four symmetries; for `m = 0` the class of `(c, d)` plus `a mod gcd(c, d)`
-- and every witness is an explicit element of `Gamma`, validated by the port's checked row-matrix boundary.

The classification is cross-checked against the **retained generator set** rather than against itself: each
retained `Gamma` generator is applied to a sample of shell vectors, and the class key must agree. That would fail if
the key were too fine (which is the plausible error, since a too-coarse key would not be caught by it).
"""
import pytest

from sage_indefinite_port.indefinite.pointwise_shell import (
    U_PLUS_U_GRAM,
    gamma_action,
    gamma_class,
    gamma_equivalent,
    gamma_witness,
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


def _shell_sample():
    """Primitive isotropic vectors: `a m + c d = 0` and content one."""
    sample = [_v(0, 1, 1, 0), _v(-6, 1, 2, 3), _v(0, 1, 3, 0), _v(3, 1, -3, 1),
              _v(0, 2, 1, 0), _v(0, 0, 1, 0), _v(1, 0, 2, 0), _v(2, 0, 3, 0),
              _v(1, 0, 0, 0), _v(-1, 0, 0, 0), _v(0, 0, 0, 1)]
    for vector in sample:
        assert int(vector.q()) == 0 and int(vector.content()) == 1, f"bad sample member {vector}"
    return sample


def test_the_shell_hypotheses_are_checked_not_assumed() -> None:
    """Non-isotropic and imprimitive vectors are refused, and the coordinates satisfy `a m + c d = 0`."""
    L = _L()
    assert shell_coordinates(L, _v(-6, 1, 2, 3)) == (-6, 1, 2, 3)
    with pytest.raises(UnfinishedCapability):
        shell_coordinates(L, element(L, (1, 0, 0, 0)) + element(L, (0, 1, 0, 0)))  # q = 2
    with pytest.raises(UnfinishedCapability):
        shell_coordinates(L, _v(0, 2, 2, 0))


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


def test_the_classification_covers_the_whole_shell() -> None:
    r"""A vector outside the shell is refused rather than classed; inside it, every vector gets a class.

    The refusal is the honest boundary: the classification uses `a m + c d = 0` and primitivity, and neither is
    meaningful outside. Nothing is silently given a class it has no business having.
    """
    L = _L()
    for vector in _shell_sample():
        assert gamma_class(L, vector) is not None
    outside = element(L, (2, 0, 0, 0))  # content 2
    with pytest.raises(UnfinishedCapability):
        gamma_class(L, outside)
