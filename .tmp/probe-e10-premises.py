"""Measure the composite-disjunct premises on the ACTUAL E10 square-0 specimen.

Builds the same action the failing case builds -- via vector_classes, then the port's pointed
path -- and stops before `.stabilizer()`, which is the call that enters the expensive GAP
IsFinite branch. Prints the three premises and the or-condition they sit in, without ever
evaluating is_finite().
"""
import time

from tests.acceptance.consumer import lattice, marked_vector, vector_classes
from tests.fixtures.oracle_fixtures import load_e10_fundamental_domain

from sage_indefinite_port.indefinite.recursive import EichlerRecursion

_E10 = load_e10_fundamental_domain()
started = time.monotonic()
lattice_ = lattice(_E10["gram"])
print("E10 built in", round(time.monotonic() - started, 1), "s", flush=True)

started = time.monotonic()
marked = next(iter(vector_classes(lattice_, 0)))
vector = marked_vector(marked)
print(f"square-0 marked class obtained in {time.monotonic() - started:.1f}s; v.q() = {vector.q()}", flush=True)

recursion = EichlerRecursion()
started = time.monotonic()
if int(vector.q()) == 0:
    action = recursion.isotropic_stabilizer_action((vector,), ())
    branch = "isotropic_stabilizer_action"
else:
    action = recursion._lattice_action(vector)
    branch = "_lattice_action"
print(f"{branch} built in {time.monotonic() - started:.1f}s  (is_finite NOT evaluated)", flush=True)

group = action.rational_group()
bound = action.stable_containing_bound()
print("--- composite disjunct premises on the E10 square-0 specimen ---", flush=True)
print("   stable_containing_bound is not None      :", bound is not None, flush=True)
print("   has_selected_group_resolution            :", group.has_selected_group_resolution(), flush=True)
print("   has_selected_finite_group_generating_set :", group.has_selected_finite_group_generating_set(), flush=True)
print("=> composite disjunct                       :",
      (group.has_selected_group_resolution()
       and group.has_selected_finite_group_generating_set()
       and bound is not None), flush=True)
print("is_finite() deliberately NOT called: it is the expensive branch under investigation", flush=True)
