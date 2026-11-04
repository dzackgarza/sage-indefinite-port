"""Does IntegralStructureAction.stabilizer's composite disjunct actually hold?

If it does, `rational_group.is_finite()` sits in a dead `or` operand and the same reorder that
00bcc86 applied elsewhere would skip the GAP query. Checked on the hyperbolic plane for both
pointed branches: isotropic (q = 0, the branch square[0] takes) and non-isotropic.
"""
from sage_extensions.all import Lattices, ZZ

from sage_indefinite_port.indefinite.recursive import EichlerRecursion

lattice = Lattices(ZZ)([[0, 1], [1, 0]])
labels = tuple(lattice.module_generating_set())
recursion = EichlerRecursion()


def report(name, action):
    group = action.rational_group()
    bound = action.stable_containing_bound()
    print(f"--- {name}")
    print(f"   stable_containing_bound is not None : {bound is not None}")
    print(f"   has_selected_group_resolution       : {group.has_selected_group_resolution()}")
    print(f"   has_selected_finite_group_generating_set: {group.has_selected_finite_group_generating_set()}")
    disjunct = (
        group.has_selected_group_resolution()
        and group.has_selected_finite_group_generating_set()
        and bound is not None
    )
    print(f"   => composite disjunct               : {disjunct}")
    print(f"   => is_finite() is True              : {group.is_finite() is True}")
    print(f"   => query is DEAD when disjunct holds: {disjunct}")


timelike = lattice.linear_combination({labels[0]: 1, labels[1]: 1})
print("timelike v q =", timelike.q())
report("non-isotropic branch (_lattice_action)", recursion._lattice_action(timelike))

isotropic = lattice.module_generator(labels[0])
print("\nisotropic v q =", isotropic.q())
report("isotropic branch (isotropic_stabilizer_action)", recursion.isotropic_stabilizer_action((isotropic,), ()))
