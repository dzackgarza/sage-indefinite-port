"""Which accessor on a degenerate lattice exposes the group whose orbit machinery is native?"""
from tests.acceptance.consumer import element, lattice

L = lattice([[2, 0], [0, 0]])
v = element(L, [1, 0])
for name in ("O", "Aut", "orthogonal_group", "bilinear_orthogonal_group", "O_plus"):
    if not hasattr(L, name):
        print(f"  {name}: absent", flush=True)
        continue
    try:
        group = getattr(L, name)()
    except Exception as error:
        print(f"  {name}(): raised {type(error).__name__}: {str(error)[:80]}", flush=True)
        continue
    surface = [attr for attr in ("orbit", "orbits", "natural_g_set", "action_on", "subgroup",
                                 "stabilizer", "group_generators") if hasattr(group, attr)]
    print(f"  {name}() -> {type(group).__name__} | {surface}", flush=True)
    if hasattr(group, "natural_g_set"):
        try:
            gset = group.natural_g_set()
            print(f"      natural_g_set -> {type(gset).__name__} | orbits? {hasattr(gset, 'orbits')}", flush=True)
            orbit = gset.orbits().orbit_of(v)
            print(f"      orbit_of -> {type(orbit).__name__} | points? {hasattr(orbit, 'points')}", flush=True)
        except Exception as error:
            print(f"      natural_g_set raised {type(error).__name__}: {str(error)[:110]}", flush=True)
