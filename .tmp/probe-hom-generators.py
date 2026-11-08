"""The Hom module's public generator contract: module_generating_set / module_generator."""
from sage_extensions.all import Lattices, Modules, ZZ as OwnedZZ

C = Lattices(OwnedZZ)([[2, 0], [0, 2]])
R = Lattices(OwnedZZ)([[0]])
hom = Modules(OwnedZZ).Mor(C, R)
print("Hom:", hom, "| type", type(hom).__name__, flush=True)
for name in ("module_generating_set", "module_generator", "module_generators", "module_rank", "evaluation"):
    print(f"  has {name}: {hasattr(hom, name)}", flush=True)
try:
    labels = tuple(hom.module_generating_set())
    print("  module_generating_set():", labels, flush=True)
    elements = tuple(hom.module_generator(label) for label in labels)
    print("  module_generator(label) for each:", [str(e)[:30] for e in elements], flush=True)
    x = C.module_generator(C.module_generating_set()[0])
    print("  evaluation of each generator at x:", [hom.evaluation(e, x) for e in elements], flush=True)
except Exception as error:
    import traceback
    print(f"  raised {type(error).__name__}: {str(error)[:200]}", flush=True)
    traceback.print_exc()
