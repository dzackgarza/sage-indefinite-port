"""Exact finite-field coverage of the quaternion rank-one residue search."""

from sage.rings.finite_rings.finite_field_constructor import GF
from sage_indefinite_port.indefinite.quaternion_orders_number_field import _projective_plane_representatives


def test_projective_residue_representatives_cover_all_nonzero_lines():
    for field in (GF(2), GF(3), GF(4, "a"), GF(5)):
        representatives = list(_projective_plane_representatives(field))
        q = field.cardinality()
        assert len(representatives) == q * q + q + 1
        assert len(set(representatives)) == len(representatives)
        for x in field:
            for y in field:
                for z in field:
                    if (x, y, z) == (0, 0, 0):
                        continue
                    first = next(c for c in (x, y, z) if c)
                    assert tuple(c / first for c in (x, y, z)) in representatives
