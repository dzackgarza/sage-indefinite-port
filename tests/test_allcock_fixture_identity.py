"""Identity invariants of the rank-3 Allcock fixture.

The acceptance node id for `test_root_systems_match_allcock` is composed from the case's
retained source identity `(weyl_group_id, lattice_id)`, because `case["id"]` elides the
Weyl-group id and names up to 374 different lattices. These checks are the guard on that
composition: they assert the fixture really does carry a unique retained identity, so the
node id stays unambiguous, and that no case is lost or duplicated along the way.

They are fixture-identity checks only. Nothing here asserts or substitutes for the
mathematical assertions of the frozen acceptance case.
"""
import hashlib
import json
from collections import Counter
from pathlib import Path

import pytest

FIXTURE = Path(__file__).parent / "fixtures" / "allcock_rank3_reflective.json"
ENTRIES = 8595
DISTINCT_IDS = 176
WEYL_GROUPS = 374


def load():
    with open(FIXTURE, encoding="utf-8") as handle:
        return json.load(handle)


def content_digest(case):
    """The case's mathematical content, independent of its reporting label."""
    body = {key: value for key, value in case.items() if key != "id"}
    return hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def composed_id(case):
    """The same composition the acceptance parametrization uses."""
    return f"{case['id']}_w{case['weyl_group_id']}_l{case['lattice_id']}"


def test_fixture_keeps_every_case():
    assert len(load()) == ENTRIES


def test_id_is_the_ambiguous_lattice_only_label():
    """Document the defect the composed id repairs: `id` alone is not an identity."""
    cases = load()
    assert all(case["id"] == f"allcock_rk3_{case['lattice_id']}" for case in cases)
    ids = Counter(case["id"] for case in cases)
    assert len(ids) == DISTINCT_IDS
    assert ids.most_common(1)[0][1] > 1


def test_retained_source_identity_is_unique():
    cases = load()
    pairs = [(case["weyl_group_id"], case["lattice_id"]) for case in cases]
    assert len(set(pairs)) == ENTRIES
    assert len({case["weyl_group_id"] for case in cases}) == WEYL_GROUPS


def test_composed_id_is_injective_and_covers_every_case():
    cases = load()
    composed = [composed_id(case) for case in cases]
    assert len(set(composed)) == len(composed) == ENTRIES


def test_no_two_cases_share_content_so_no_occurrence_counter_is_needed():
    cases = load()
    digests = Counter(content_digest(case) for case in cases)
    assert len(digests) == ENTRIES
    assert all(count == 1 for count in digests.values())


def test_every_case_cites_its_source():
    for case in load():
        source = case["source"]
        assert source["kind"] == "published_classification"
        assert "GeometryDatabase_Rank3_Lorentzian_lattices" in source["file"]


def test_mathematical_data_is_preserved_per_case():
    """The identity fix must not touch the data it identifies."""
    for case in load():
        assert len(case["gram"]) == 3
        assert len(case["simple_roots"]) == case["num_simple_roots"]
        assert len(case["elementary_divisors"]) == 3


@pytest.mark.parametrize(
    "label,weyl,lattice",
    [
        ("allcock_rk3_4", 0, 4),
        ("allcock_rk3_2", 2, 2),
    ],
)
def test_the_recorded_sample_cases_resolve_to_one_lattice(label, weyl, lattice):
    """The #33 branch sample names these two by id *and* Gram; the id alone does not pin them."""
    cases = load()
    matching = [
        case
        for case in cases
        if case["id"] == label and (case["weyl_group_id"], case["lattice_id"]) == (weyl, lattice)
    ]
    assert len(matching) == 1, f"{label} w{weyl} l{lattice} does not resolve to one case"
    assert [case for case in cases if case["id"] == label].__len__() > 1, (
        f"{label} is expected to be ambiguous; if it is unique the composed id is unnecessary"
    )
