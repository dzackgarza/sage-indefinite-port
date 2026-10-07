"""The consumer side of the acceptance suite: corpus data in, preamble objects out.

Every case reaches this port the way the research preamble's consumers do: a lattice is
``Lattices(ZZ)(gram)`` (or a named preamble lattice), and every group, orbit,
stabilizer, witness and cusp is a method of that lattice, of its orthogonal group, or of
a subgroup of it. These helpers only translate corpus coordinates into preamble elements
and read results back as coordinates, so that results can be compared with recorded data.

Matrices follow the corpus's row convention: row ``i`` holds the coordinates of the image
of the ``i``-th basis vector, so an isometry ``R`` of a lattice with Gram matrix ``G``
satisfies ``R G R^T = G``, and a witness ``R`` from ``L_1`` to ``L_2`` satisfies
``R G_2 R^T = G_1``.
"""

import importlib
from functools import reduce
from math import gcd

from dzack_research.preamble.all import ZZ as PreambleZZ
from dzack_research.preamble.all import Lattices
from sage.all import GF, ZZ, MatrixGroup, matrix, prime_divisors

from sage_indefinite_port.readiness import UNFINISHED, UnfinishedCapability

O_L = "sage_indefinite_port.indefinite.recursive:orthogonal_group_generators"
ISOMETRY = "sage_indefinite_port.indefinite.recursive:isometry"
VECTOR_WITNESS = "sage_indefinite_port.indefinite.recursive:vector_equivalence_witness"
VECTOR_STABILIZER = "sage_indefinite_port.indefinite.recursive:vector_stabilizer_generators"
VECTOR_ORBITS = "sage_indefinite_port.indefinite.recursive:vector_orbit_representatives"
ISOTROPIC_ORBITS = "sage_indefinite_port.indefinite.isotropic_flags:isotropic_sublattice_orbit_representatives"
ISOTROPIC_WITNESS = "sage_indefinite_port.indefinite.isotropic_flags:isotropic_sublattice_equivalence_witness"
ISOTROPIC_STABILIZER = "sage_indefinite_port.indefinite.isotropic_flags:isotropic_sublattice_stabilizer_generators"
FLAG_ORBITS = "sage_indefinite_port.indefinite.isotropic_flags:isotropic_flag_orbit_representatives"
SPLIT_ORBIT = "sage_indefinite_port.groups.finite_index:split_orbit"
SUBGROUP_GENERATORS = "sage_indefinite_port.groups.finite_index:subgroup_generators"
EQUIVARIANT_LATTICE = "sage_indefinite_port.groups.equivariant:EquivariantLattice"
EDGEWALK = "sage_indefinite_port.indefinite.edgewalk:edgewalk_fundamental_domain"
PERFECT_DOMAINS = "sage_indefinite_port.indefinite.lorentzian_cells:perfect_domain_traversal"


def require(*entry_points: str) -> None:
    """The first line of a case: each port entry point it needs exists and is finished.

    A missing module, a missing name, or an entry point still marked ``@unfinished`` raises
    ``UnfinishedCapability`` at once, before any of the preamble work that precedes the missing
    capability. Every other exception, such as a module that exists but does not import,
    propagates as itself, so the case's xfail (``raises=UnfinishedCapability``) does not
    absorb it.
    """
    for entry_point in entry_points:
        module_name, _, name = entry_point.partition(":")
        try:
            module = importlib.import_module(module_name)
        except ModuleNotFoundError as error:
            if error.name != module_name:
                raise
            raise UnfinishedCapability(f"{entry_point} does not exist yet: its module is not written") from error
        if not hasattr(module, name):
            raise UnfinishedCapability(f"{entry_point} does not exist yet")
        if entry_point in UNFINISHED:
            raise UnfinishedCapability(f"{entry_point} is unfinished: work unit #{UNFINISHED[entry_point]} owns it")


def lattice(gram):
    return Lattices(PreambleZZ)([[int(entry) for entry in row] for row in gram])


def basis(lattice_):
    return tuple(lattice_.module_generators().values())


def element(lattice_, coordinates):
    return lattice_([int(entry) for entry in coordinates])


def coordinates(vector) -> tuple[int, ...]:
    return tuple(int(entry) for entry in vector.to_vector())


def content(vector) -> int:
    return reduce(gcd, coordinates(vector), 0)


def primitive(vectors):
    return tuple(vector for vector in vectors if content(vector) == 1)


def isometry_from_rows(lattice_, rows):
    return lattice_.Aut()([element(lattice_, row) for row in rows])


def rows_of(morphism) -> list[list[int]]:
    return [list(coordinates(morphism(vector))) for vector in basis(morphism.domain())]


def is_isometry(rows, gram) -> bool:
    images = matrix(ZZ, rows)
    form = matrix(ZZ, gram)
    return images * form * images.transpose() == form


def is_witness(rows, source_gram, target_gram) -> bool:
    images = matrix(ZZ, rows)
    return images * matrix(ZZ, target_gram) * images.transpose() == matrix(ZZ, source_gram)


def discriminant_image_order(lattice_, isometries) -> int:
    target = lattice_.discriminant_group().orthogonal_group()
    return int(target.subgroup_on(tuple(isometry.discriminant_morphism() for isometry in isometries)).cardinality())


def reduction_image_order(generator_rows, prime: int) -> int:
    """The order of the image in GL_n(F_p) of the group generated by the given row matrices."""
    return int(MatrixGroup([matrix(GF(prime), rows) for rows in generator_rows]).order())


def comparison_primes(gram) -> tuple[int, ...]:
    """The primes dividing 2 det G, and 3, 5 and 7: the reductions on which two generating sets are compared."""
    determinant = int(matrix(ZZ, gram).determinant())
    return tuple(sorted(set(prime_divisors(2 * determinant)) | {3, 5, 7}))


def assert_same_finite_images(gram, computed, recorded_rows) -> None:
    """Assert that computed isometries and recorded generators have the same finite images.

    The images compared are the image in O(q_L) and the reductions modulo the primes of
    ``comparison_primes``. Recorded generators come from the reference implementation or
    the literature as row matrices.
    """
    lattice_ = computed[0].domain()
    computed_rows = [rows_of(isometry) for isometry in computed]
    assert all(is_isometry(rows, gram) for rows in computed_rows)
    assert all(is_isometry(rows, gram) for rows in recorded_rows)
    recorded = [isometry_from_rows(lattice_, rows) for rows in recorded_rows]
    assert discriminant_image_order(lattice_, computed) == discriminant_image_order(lattice_, recorded)
    for prime in comparison_primes(gram):
        assert reduction_image_order(computed_rows, prime) == reduction_image_order(recorded_rows, prime), prime


def transpose(rows):
    return [list(column) for column in zip(*rows, strict=True)]
