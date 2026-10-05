import pytest
from sage.matrix.constructor import matrix
from sage.rings.rational_field import QQ

from sage_indefinite_port.indefinite.edgewalk_extensions import (
    COXETER_INFINITY,
    compute_coxeter_matrix,
    compute_possible_extensions,
)


def test_compute_coxeter_matrix_recovers_crystallographic_labels():
    gram = matrix(QQ, [[2, -1, 0], [-1, 2, -2], [0, -2, 4]])
    coxeter, scalar = compute_coxeter_matrix(gram, ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    assert scalar == gram
    assert coxeter == matrix(QQ, [[2, 3, 2], [3, 2, 4], [2, 4, 4]])


def test_compute_coxeter_matrix_recovers_infinite_edge():
    coxeter, _ = compute_coxeter_matrix(
        matrix(QQ, [[2, -2], [-2, 2]]), ((1, 0), (0, 1))
    )
    assert coxeter[0, 1] == COXETER_INFINITY


def test_A1_extensions_match_source_spherical_and_affine_candidates():
    gram = matrix(QQ, [[2]])
    spherical = compute_possible_extensions(gram, ((1,),), (2,), only_spherical=True)
    assert {x.coxeter_entries for x in spherical} == {(QQ(2),), (QQ(3),)}
    affine = compute_possible_extensions(gram, ((1,),), (2,))
    assert {x.coxeter_entries for x in affine} == {
        (QQ(2),),
        (QQ(3),),
        (COXETER_INFINITY,),
    }
    assert {x.coxeter_entries: x.residual_norm for x in affine}[
        (COXETER_INFINITY,)
    ] == 0


def test_norm_ratios_filter_double_and_triple_edges_source_faithfully():
    gram = matrix(QQ, [[2]])
    labels4 = {
        x.coxeter_entries for x in compute_possible_extensions(gram, ((1,),), (4,))
    }
    assert (QQ(4),) in labels4 and (QQ(3),) not in labels4 and (QQ(6),) not in labels4
    labels6 = {
        x.coxeter_entries for x in compute_possible_extensions(gram, ((1,),), (6,))
    }
    assert (QQ(6),) in labels6 and (QQ(3),) not in labels6 and (QQ(4),) not in labels6


def test_non_dynkin_root_pair_is_rejected():
    with pytest.raises(ValueError, match="unsupported Coxeter quotient"):
        compute_coxeter_matrix(matrix(QQ, [[2, -1], [-1, 4]]), ((1, 0), (0, 1)))
