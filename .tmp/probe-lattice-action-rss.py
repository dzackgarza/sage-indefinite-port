"""Which construction step inside _lattice_action allocates?

Reimplements the six steps of EichlerRecursion._lattice_action on the rank-2 hyperbolic plane
for v = e0 + e1, printing process RSS after each, and stopping before it can threaten the host.
No port code is modified.
"""
import sys
import time

import sage_indefinite_port.indefinite.recursive as R
from sage_extensions.all import Lattices, ZZ

LIMIT_MIB = 2560


def rss_mib():
    with open("/proc/self/status") as handle:
        for line in handle:
            if line.startswith("VmRSS:"):
                return int(line.split()[1]) // 1024
    return -1


STARTED = time.monotonic()


def step(name, call):
    before = rss_mib()
    started = time.monotonic()
    value = call()
    after = rss_mib()
    print(f"{name:38s} {time.monotonic() - started:6.2f}s  RSS {before:5d} -> {after:5d} MiB "
          f"(+{after - before})", flush=True)
    if after > LIMIT_MIB:
        print(f"STOP: RSS {after} MiB exceeds {LIMIT_MIB} MiB", flush=True)
        sys.exit(2)
    return value


lattice = Lattices(ZZ)([[0, 1], [1, 0]])
labels = tuple(lattice.module_generating_set())
vector = lattice.linear_combination({labels[0]: 1, labels[1]: 1})
recursion = R.EichlerRecursion()
print(f"rank {int(lattice.module_rank())}, v.q() = {vector.q()}, start RSS {rss_mib()} MiB", flush=True)

perpendicular = step("1 perpendicular inclusion", lambda: vector.orthogonal_complement().inclusion())
complement = step("2 complement domain", lambda: perpendicular.domain())
extensions = step("3 orthogonal_group_generators(complement)  [.extension_across]",
                  lambda: tuple(generator.extension_across(vector, vector)
                                for generator in recursion.orthogonal_group_generators(complement)))
rational = step("4 lattice.base_change(QQ)", lambda: lattice.base_change(lattice.base_ring().fraction_field_map()))
inclusion = step("5 generic_fibre_map", lambda: R.ModulesOverIntegralDomains(lattice.base_ring())(lattice).generic_fibre_map())
space = inclusion.codomain()
split = step("6 lattice.subobject_on(v, perp generators)",
             lambda: lattice.subobject_on((vector, *(perpendicular(point) for point in complement.module_generators()))))
denominator = step("7 split cokernel exponent", lambda: split.inclusion().cokernel().exponent())
scale = rational.base_ring().one() / rational.base_ring()(denominator)


def build_bound():
    return R.Modules(lattice.base_ring()).Mono(split, space)(
        lambda label: space.wrap(rational.scalar_multiple(
            scale, inclusion(split.inclusion()(split.module_inject(label) if False else split.module_generator(label))).underlying_element()))
    )


containing = step("8 stable containing bound (Mono)", build_bound)
action = step("9 IntegralStructureAction(...)", 
              lambda: R.IntegralStructureAction(rational.O().subgroup(extensions), inclusion).with_stable_containing_bound(containing))
print("FINAL RSS", rss_mib(), "MiB after", round(time.monotonic() - STARTED, 1), "s", flush=True)
print("action:", action, flush=True)
