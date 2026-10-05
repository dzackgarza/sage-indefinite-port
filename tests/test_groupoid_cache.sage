from dzack_research.preamble.all import Lattices, ZZ

from sage_indefinite_port.indefinite.presentations import ReducedLatticePresentation, presentation_bucket_key
from sage_indefinite_port.groups.groupoid_cache import IsometryGroupoidCache
from tests.fixtures.oracle_fixtures import load_conway_sloane_cases


def test_reduced_presentation_retains_verified_isometry() -> None:
    lattice = Lattices(ZZ)("A2")
    presentation = ReducedLatticePresentation.from_lattice(lattice)

    assert presentation.source is presentation.isometry.domain()
    assert presentation.reduced is presentation.isometry.codomain()


def test_splag_51a_51b_share_bucket_but_exactly_fail_isometry() -> None:
    pair = load_conway_sloane_cases()["spinor_genus_pair_determinant_minus_128"]
    first = Lattices(ZZ)(pair["form_a"]["gram"])
    second = Lattices(ZZ)(pair["form_b"]["gram"])

    assert presentation_bucket_key(first) == presentation_bucket_key(second)
    assert pair["integrally_equivalent"] is False


def test_cached_orthogonal_group_transports_by_verified_conjugation() -> None:
    lattice = Lattices(ZZ)("A2")
    presentation = ReducedLatticePresentation.from_lattice(lattice)
    cache = IsometryGroupoidCache()
    target = presentation.reduced
    generators = tuple(target.Aut().framing().group_generators())

    cache.remember_isometry(presentation.isometry)
    cache.remember_orthogonal_group(target, generators)
    transported = cache.lookup_orthogonal_group(presentation.source)

    assert transported is not None
    assert all(generator.parent() is presentation.source.Aut() for generator in transported)
