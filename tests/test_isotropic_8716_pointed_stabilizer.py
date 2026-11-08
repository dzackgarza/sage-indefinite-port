"""The complete pointed null stabilizer in the rank-six IsotropicCases_8716 lattice.

Its positive-definite quotient has full group O(<1,1> + 5 A2) and its
translation kernel is exactly the even-norm sublattice. The proof is in
references/isotropic-8716-null-stabilizer.md.
"""

from sage.all import QQ, ZZ, identity_matrix, matrix, vector
from sage.libs.gap.libgap import libgap
from sage_extensions.categories.orthogonal_quotients import _GapWordSection

from tests.acceptance.consumer import (
    automorphisms, element, lattice, pointed, pointed_category, rows_of, underlying_isometry,
)
from tests.fixtures.oracle_fixtures import load_isotropic_cases


def test_pointed_null_stabilizer_has_the_full_levi_group_and_even_translation_kernel():
    case = next(case for case in load_isotropic_cases() if case["id"] == "IsotropicCases_8716")
    space = lattice(case["gram"])
    v = element(space, (1, 0, 0, 1, 0, 0))
    marked = pointed(space, v)
    isometries = tuple(underlying_isometry(g) for g in automorphisms(pointed_category(v.q()), marked).group_generators())
    assert all(g(v) == v and (~g)(v) == v for g in isometries)
    data = {"gram": case["gram"], "vector": (1, 0, 0, 1, 0, 0), "actions": tuple(rows_of(g) for g in isometries)}
    gram=matrix(ZZ,data['gram'])
    actions=tuple(matrix(ZZ,rows) for rows in data['actions'])
    indices=(1,2,4,5)
    quotient_gram=gram.matrix_from_rows_and_columns(indices,indices)
    null=vector(ZZ,data['vector'])
    quotients=[]
    for action in actions:
        assert action*gram*action.transpose()==gram
        assert action.det() in (-1,1) and null*action==null and null*action.inverse()==null
        quotient=action.matrix_from_rows_and_columns(indices,indices)
        assert quotient*quotient_gram*quotient.transpose()==quotient_gram
        quotients.append(quotient)
    section=_GapWordSection(tuple(libgap(q.change_ring(QQ)) for q in quotients))
    assert section.order()==96
    permutation_group=section.permutation_group
    degree=int(libgap.LargestMovedPoint(permutation_group))
    def key(point):
        return tuple(int(entry) for entry in libgap.ListPerm(point,degree))
    lifts={}
    inverses={}
    for point in libgap.AsList(permutation_group):
        action=identity_matrix(ZZ,6)
        for letter in section.letters_of_permutation(point):
            action=action*(actions[letter-1] if letter>0 else actions[-letter-1].inverse())
        assert libgap.Image(section._to_permutation, libgap(action.matrix_from_rows_and_columns(indices, indices).change_ring(QQ))) == point
        index=key(point)
        lifts[index]=action
        inverses[index]=action.inverse()
    parameters=[]
    letters=tuple(libgap.Image(section._to_permutation,libgap(q.change_ring(QQ))) for q in quotients)
    for point in libgap.AsList(permutation_group):
        for action,letter in zip(actions,letters,strict=True):
            kernel=lifts[key(point)]*action*inverses[key(point*letter)]
            parameter=vector(ZZ,[kernel[0,index] for index in indices])
            square=parameter*quotient_gram*parameter
            assert square%2==0
            expected=identity_matrix(ZZ,6)
            for position,index in enumerate(indices):
                expected[0,index]=parameter[position]
                expected[3,index]=-parameter[position]
                pairing=(quotient_gram*parameter)[position]
                expected[index,0]=-pairing
                expected[index,3]=-pairing
            expected[0,0]-=square//2
            expected[0,3]-=square//2
            expected[3,0]+=square//2
            expected[3,3]+=square//2
            assert kernel==expected
            parameters.append(parameter)
    actual=matrix(ZZ,parameters).row_module()
    expected=matrix(ZZ,[[1,1,0,0],[0,2,0,0],[0,0,1,0],[0,0,0,1]]).row_module()
    assert actual==expected
