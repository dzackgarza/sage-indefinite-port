"""Show the pre-fix expression raises on the same object the fixed helper accepts."""
from sage_indefinite_port.indefinite.isotropic_induction import IsotropicInduction
from sage_indefinite_port.indefinite.edgewalk_scalars import _integer
from tests.acceptance.consumer import lattice

lattice_ = lattice([[0, 1], [1, 0]])
vector = next(iter(lattice_.module_generators()))
rational = lattice_.base_change(lattice_.base_ring().fraction_field_map())
point = rational.linear_combination({label: rational.base_ring()(1 if index == 0 else 0)
                                     for index, label in enumerate(lattice_.module_generating_set())})
moved = point.to_vector()
label = lattice_.module_generating_set()[0]
entry = moved(label)
print("coordinate container :", type(entry).__module__ + "." + type(entry).__name__)
try:
    print("pre-fix  int(entry)  ->", int(entry))
except TypeError as error:
    print("pre-fix  int(entry)  -> TypeError:", error)
print("fixed    _integer(entry) ->", _integer(entry))
