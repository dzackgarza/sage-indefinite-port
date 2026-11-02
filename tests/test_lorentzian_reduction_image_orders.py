"""The LORENTZ reduction images: the acceptance property, and identification where the form is nondegenerate.

`test_lorentzian_generators_match_polyhedral_common` certifies a generating set by comparing image orders
modulo every prime in `comparison_primes(gram)` = {2, 3, 5, 7}. At p=5 and p=7 that comparison does not
terminate in the resource envelope (see `.tmp/handoff-owner-lorentz-image-order.md`), so the acceptance node
stays open. This regression pins what *is* affordable, at p=2 and p=3:

* always: the computed and the recorded generating set have the **same** image order -- the property the
  acceptance assertion checks, and the only one that is meaningful when `p` divides `det L`;
* when `p` does not divide `det L`: that order is exactly the order of the full even-type orthogonal group

      |O^+_{2n}(q)| = 2 q^{n(n-1)} (q^n - 1) prod_{i=1}^{n-1} (q^{2i} - 1),   n = 5,

  computed from the classical formula rather than from GAP, so the pin does not just restate the
  computation it checks.

`LORENTZ_OrthogonalGroup1` has determinant -256, so at p = 2 the reduction is degenerate and the image is not
a full orthogonal group -- its measured order there is 2^30 * 3^5 * 5^2 * 7. That is why the identification
step is conditional and why it is not asserted at p=2 for that case.
"""
import pytest

from sage.all import ZZ, matrix

from tests.acceptance.consumer import lattice, reduction_image_order, rows_of
from tests.fixtures.oracle_fixtures import load_lorentzian_stabilizer_cases


def orthogonal_plus_order(n: int, q: int) -> int:
    """|O^+_{2n}(q)| from the classical order formula."""
    value = 2 * q ** (n * (n - 1))
    for index in range(1, n):
        value *= q ** (2 * index) - 1
    return value * (q ** n - 1)


_CASES = {case["id"]: case for case in load_lorentzian_stabilizer_cases()}
_IDS = ["LORENTZ_OrthogonalGroup0", "LORENTZ_OrthogonalGroup1", "LORENTZ_OrthogonalGroup2"]


@pytest.mark.parametrize("case_id", _IDS)
def test_reduction_images_agree_and_are_orthogonal_groups_where_nondegenerate(case_id) -> None:
    case = _CASES[case_id]
    determinant = int(matrix(ZZ, case["gram"]).det())
    lattice_ = lattice(case["gram"])
    computed_rows = [rows_of(isometry) for isometry in lattice_.O().group_generators()]

    for prime in (2, 3):
        computed_order = reduction_image_order(computed_rows, prime)
        recorded_order = reduction_image_order(case["generators"], prime)
        assert computed_order == recorded_order, (case_id, prime)
        if determinant % prime:
            assert computed_order == orthogonal_plus_order(5, prime), (case_id, prime)
