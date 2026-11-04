"""The maintained route to {B in Hom(C,R) : B x = 0}: Hom, evaluation morphism, kernel."""
from sage_extensions.all import Lattices, Modules, ZZ as OwnedZZ
from tests.acceptance.consumer import element

# C = <2> + <2> (rank 2), R = zero form rank 1 (the radical)
C = Lattices(OwnedZZ)([[2, 0], [0, 2]])
R = Lattices(OwnedZZ)([[0]])
modules = Modules(OwnedZZ)
print("C:", C, "| R:", R, flush=True)

hom = modules.Mor(C, R)
print("Hom(C,R):", type(hom).__name__, "| evaluation:", hasattr(hom, "evaluation"), flush=True)
print("  rank/生成 info:", [m for m in ("module_rank", "generators", "basis", "an_element") if hasattr(hom, m)],
      flush=True)

x = C.module_generator(C.module_generating_set()[0])      # (1,0), primitive
print("x:", x, flush=True)
try:
    evaluation_morphism = modules.Mor(hom, R)(lambda map_element: hom.evaluation(map_element, x))
    print("evaluation morphism:", type(evaluation_morphism).__name__, flush=True)
    kernel = evaluation_morphism.kernel()
    print("kernel:", type(kernel).__name__, flush=True)
    print("  surface:", [m for m in ("generators", "basis", "inclusion", "module_rank", "domain")
                          if hasattr(kernel, m)], flush=True)
    for name in ("module_rank",):
        if hasattr(kernel, name):
            print(f"  {name}: {getattr(kernel, name)()}", flush=True)
    if hasattr(kernel, "generators"):
        generators = tuple(kernel.generators())
        print("  generators:", len(generators), [str(g)[:40] for g in generators], flush=True)
except Exception as error:
    import traceback
    print(f"raised {type(error).__name__}: {str(error)[:200]}", flush=True)
    traceback.print_exc()
