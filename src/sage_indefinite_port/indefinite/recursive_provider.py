"""Engine-boundary adapters for recursive indefinite lattice operations."""

from dzack_research.preamble.all import ZZ, Lattices

from sage_indefinite_port.indefinite.lorentzian_cells import IndefiniteOrthogonalAlgorithm


def _lattice(gram):
    return Lattices(ZZ)(gram)


def _vector(lattice, coordinates):
    labels = tuple(lattice.module_generating_set())
    return lattice.linear_combination({label: lattice.base_ring()(int(coefficient)) for label, coefficient in zip(labels, coordinates, strict=True) if coefficient})


def _row_matrix(isometry):
    return [[int(entry) for entry in row] for row in isometry.domain().Aut()._row_action_matrix(isometry).rows()]


def indefinite_automorphism_group(gram):
    lattice = _lattice(gram)
    group = IndefiniteOrthogonalAlgorithm().orthogonal_group(lattice)
    return tuple(_row_matrix(generator) for generator in group.generators())


def indefinite_isometry_witness(source_gram, target_gram):
    source = _lattice(source_gram)
    target = _lattice(target_gram)
    witness = IndefiniteOrthogonalAlgorithm().isometry(source, target)
    return None if witness is None else _row_matrix(witness)


def indefinite_vector_isometry_witness(gram, source_coordinates, target_coordinates):
    lattice = _lattice(gram)
    source = _vector(lattice, source_coordinates)
    target = _vector(lattice, target_coordinates)
    witness = IndefiniteOrthogonalAlgorithm().vector_transporter(source, target)
    return None if witness is None else _row_matrix(witness)


def indefinite_vector_stabilizer(gram, coordinates):
    lattice = _lattice(gram)
    vector = _vector(lattice, coordinates)
    subgroup = IndefiniteOrthogonalAlgorithm().vector_stabilizer(vector)
    return tuple(_row_matrix(generator) for generator in subgroup.generators())


__all__ = [
    "indefinite_automorphism_group",
    "indefinite_isometry_witness",
    "indefinite_vector_isometry_witness",
    "indefinite_vector_stabilizer",
]
