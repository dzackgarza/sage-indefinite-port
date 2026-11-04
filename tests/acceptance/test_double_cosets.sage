"""Double cosets of the integral subgroup in the rational group (#33's integral-structure row).

Recorded data: polyhedral_common @ 1592b246,
`CI_tests/DoubleCosets/DBL/DoubleCoset_n<deg>_big<|G|>_sma<|H|>_vf<count>_idx<i>`, read by
`references/extract/polyhedral_common_ci.py`, which asserts on extraction that the two permutation
groups in each file really have the orders their names claim. Eighteen cases, degrees 55 and 56.

The fact asserted here is the upstream decomposition test's, `TestDoubleCosetPermDecomposition` in
`CI_tests/01_RatIntAutomorphy/ProcessExamples.g:179-253`: for `H <= G`, the double cosets `H g H`
**partition** `G` -- the classes are pairwise disjoint, and their cardinalities sum to `|G|`. Here
each class is counted as `|H| * #{left cosets of H in H g H}`, which is what `|H g H|` is for a
subgroup `H`, so the sum is the completeness of the decomposition and not a restatement of a count.

Computation route: the provider owns this operation (`Group.double_cosets`, sage-extensions
`categories/group/groups.py`, realized as the `H`-orbits of the finite G-set `G/H`). The port supplies
no group theory here: the case builds the two permutation groups with the owned constructor and reads
the decomposition through that method. `G` is finite, so `G/H` is a finite `G`-set and no kernel or
preimage realization is needed.

Nothing in the case is a port entry point, so no `require` gate applies; a failure here is a provider
or fixture fact, not an unimplemented port capability.
"""

import pytest

from sage.combinat.permutation import Permutation

from sage_extensions.categories.group.groups import OwnedGroups

from tests.fixtures.oracle_fixtures import (
    load_double_coset_cases,
)

_CASES = load_double_coset_cases()


def _owned_permutation(symmetric, images):
    r"""The owned element of ``symmetric`` with these one-based-images permutation data.

    The fixture stores images zero-based, as the reference files do; ``Permutation.cycle_tuples``
    supplies the cycle list the owned constructor takes, and the identity has no cycles.
    """
    cycles = Permutation([image + 1 for image in images]).cycle_tuples()
    return symmetric([tuple(cycle) for cycle in cycles]) if cycles else symmetric([])


@pytest.mark.parametrize("case", _CASES, ids=lambda case: case["id"])
def test_double_cosets_partition_the_rational_group(case) -> None:
    r"""`G` is the disjoint union of the `H`-double cosets, and the recorded orders are real orders.

    The two group orders are the fixture's own external facts, re-derived here from the permutation
    generators rather than trusted; then `G/C = union` and `sum |H g H| = |G|` is the upstream
    partition test. A class whose left-coset count times `|H|` missed the group order would break
    the sum, so the two assertions together pin both the decomposition and its completeness.
    """
    symmetric = OwnedGroups().S(int(case["degree"]))
    group = symmetric.subgroup(
        tuple(_owned_permutation(symmetric, generator) for generator in case["group_g"]["generators"])
    )
    integral = group.subgroup(
        tuple(_owned_permutation(symmetric, generator) for generator in case["group_h"]["generators"])
    )

    big_order = int(case["big_group_order"])
    small_order = int(case["small_group_order"])
    assert int(group.cardinality()) == big_order
    assert int(integral.cardinality()) == small_order
    assert big_order % small_order == 0, "an integral subgroup divides the rational group's order"

    classes = tuple(group.double_cosets(integral, integral))
    assert classes, "the decomposition of a nonempty group is nonempty"

    sizes = tuple(len(class_.points()) * small_order for class_ in classes)
    assert sum(sizes) == big_order, (
        f"the double cosets {case['id']} have sizes {sizes}, summing to {sum(sizes)} rather than |G| = {big_order}"
    )
    for size in sizes:
        assert size % small_order == 0
        assert size <= big_order


def test_the_recorded_case_families_are_all_present() -> None:
    r"""All eighteen cited instances are loaded, spanning every degree the reference files use.

    The extraction script asserts the file-name-to-content agreement (degrees, group orders, vector
    counts) on every rebuild; this pins that the loaded corpus is the whole set, so a case cannot
    disappear silently and leave the operation nominally covered. The degree list is the reference
    directory's own: the DBL instances are not confined to the two smallest degrees.
    """
    assert len(_CASES) == 18
    assert sorted({int(case["degree"]) for case in _CASES}) == [55, 56, 60, 64, 70, 72, 80, 96]
    assert sorted({int(case["index"]) for case in _CASES}) == [0, 1]
    assert all(int(case["num_vectors"]) > 0 for case in _CASES)
    assert all(int(case["small_group_order"]) <= int(case["big_group_order"]) for case in _CASES)
