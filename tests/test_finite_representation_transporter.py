"""The finite-representation transporter is reachable and works on a rational-nonintegral specimen.

`EichlerRecursion.pointed_isometry` extends an isometry of `v^perp` across `v`; when the extended rational map
is **not** integral on the lattices, it corrects it through
`self._lattice_action(target).transporter(image)` — the provider's `IntegralStructureAction.transporter`, whose
finite representation is built by `_FiniteIntegralRepresentation` from the residue module `M/dM`
(`rational_integral_stabilizers.py:112-125`, `:505-520`).

A backlog note recorded that this transporter "still fails when it requests cardinality of a represented
finite torsion module", written before the provider's cardinality-admission fix. The first specimen I tried
(divisibility-four isotropic `U(4)+U(4)`) never enters the branch at all — it finds its witness on the
integral path — so it could not settle the note either way. The specimen below does enter it, and at provider
`42d7361c` the transporter returns a witness.
"""
import pytest

from sage_extensions.categories.rational_integral_stabilizers import IntegralStructureAction

from sage_indefinite_port.indefinite.recursive import EichlerRecursion
from tests.acceptance.consumer import element, isometry_from_rows, lattice, rows_of


@pytest.fixture()
def transporter_calls(monkeypatch):
    """Count calls to the provider transporter without changing what it does."""
    calls = []
    original = IntegralStructureAction.transporter

    def counted(self, *args, **kwargs):
        try:
            result = original(self, *args, **kwargs)
        except Exception as error:          # recorded, then re-raised: the caller must still see it
            calls.append(("raised", type(error).__name__, str(error)[:160]))
            raise
        calls.append(("ok", result is not None))
        return result

    monkeypatch.setattr(IntegralStructureAction, "transporter", counted)
    return calls


def test_pointed_isometry_reaches_and_survives_the_finite_representation_transporter(transporter_calls) -> None:
    """`(-2,-2,-1) -> (2,2,-1)` in `U + <-2>`: the rational transport is non-integral, the correction works."""
    gram = [[0, 1, 0], [1, 0, 0], [0, 0, -2]]
    lattice_ = lattice(gram)
    source = element(lattice_, [-2, -2, -1])
    target = element(lattice_, [2, 2, -1])

    # Same square and both primitive: the prefilter cannot separate them.
    assert int(source.q()) == int(target.q()) == 6
    assert int(source.content()) == int(target.content()) == 1

    witness = EichlerRecursion().pointed_isometry(source, target)

    assert witness is not None, "the transporter is the correction for a non-integral rational transport"
    assert witness(source) == target
    assert isometry_from_rows(lattice_, rows_of(witness)) is not None

    # The claim that matters: the correction branch was actually entered, and it did not fail.
    assert transporter_calls, "the rational transport was integral, so this specimen cannot test the transporter"
    assert all(event[0] == "ok" for event in transporter_calls), transporter_calls
