from dzack_research.preamble.all import Lattices, ZZ

from sage_indefinite_port.indefinite.lorentzian_cells import IndefiniteOrthogonalAlgorithm
from sage_indefinite_port.indefinite.vector_sections import orthogonal_section
from sage_indefinite_port.invariants import AttackProfile, VectorPrefilter, vector_content


def test_attack_profile_normalizes_sign() -> None:
    negative = Lattices(ZZ)("A2")
    direct = AttackProfile.from_lattice(negative)
    assert direct.positive_index == 0
    assert direct.sign == 1

    positive = negative.twist(-ZZ.one())
    flipped = AttackProfile.from_lattice(positive)
    assert flipped.positive_index == 0
    assert flipped.sign == -1


def test_vector_content_and_prefilter_are_distinct_from_pairing_divisor() -> None:
    plane = Lattices(ZZ)("U")
    e, f = plane.module_generators()
    vector = plane.scalar_multiple(ZZ.one() + ZZ.one(), e) + f
    assert vector_content(vector) == 1

    isotropic = VectorPrefilter.from_vector(e)
    assert isotropic.norm == 0
    assert isotropic.content == 1
    assert isotropic.divisor == 1
    assert isotropic.orthogonal_reduction_prefilter is not None


def test_recursive_dispatch_uses_definite_leaf() -> None:
    lattice = Lattices(ZZ)("A2")
    group = IndefiniteOrthogonalAlgorithm().orthogonal_group(lattice)
    assert group.domain() is lattice


def test_vector_stabilizer_lifts_reduced_group_and_fixes_vector() -> None:
    plane = Lattices(ZZ)("U")
    e, f = plane.module_generators()
    vector = e + f
    stabilizer = IndefiniteOrthogonalAlgorithm().vector_stabilizer(vector)
    assert all(generator(vector) == vector for generator in stabilizer.generators())


def test_vector_transporter_lifts_exact_reduced_isometry() -> None:
    plane = Lattices(ZZ)("U")
    e, f = plane.module_generators()
    source = plane.scalar_multiple(ZZ.one() + ZZ.one(), e) + f
    target = e + plane.scalar_multiple(ZZ.one() + ZZ.one(), f)
    transporter = IndefiniteOrthogonalAlgorithm().vector_transporter(source, target)
    assert transporter is not None
    assert transporter(source) == target



def test_higher_witt_branch_selects_two_u_model_and_rank_drop() -> None:
    zero = Lattices(ZZ)(ZZ.free_module(0))
    lattice = zero.two_u_eichler_model().lattice()
    algorithm = IndefiniteOrthogonalAlgorithm()
    model = algorithm._two_u_cover_model(lattice)
    vector = model.choose_splitting_vector()
    section = orthogonal_section(vector)
    assert algorithm.attack_profile(lattice).positive_index == 2
    assert algorithm.attack_profile(section.reduced_object()).positive_index == 1


def test_ci_indefinite_2u_public_orthogonal_group() -> None:
    lattice = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    group = IndefiniteOrthogonalAlgorithm().orthogonal_group(lattice)
    assert all(generator.parent() is lattice.Aut() for generator in group.generators())
