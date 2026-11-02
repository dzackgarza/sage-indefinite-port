"""Build ev_x: Hom(C,R) -> R from the generator images, then take its kernel()."""
from sage_extensions.all import Lattices, Modules, ZZ as OwnedZZ

C = Lattices(OwnedZZ)([[2, 0], [0, 2]])
R = Lattices(OwnedZZ)([[0, 0], [0, 0]])          # radical rank 2 for a non-trivial kernel
hom = Modules(OwnedZZ).Mor(C, R)
modules = Modules(OwnedZZ)
labels = tuple(hom.module_generating_set())
generators = tuple(hom.module_generator(label) for label in labels)
x = C.module_generator(C.module_generating_set()[0])       # (1,0), primitive
print("Hom rank:", hom.module_rank(), "| labels:", labels, flush=True)
images = tuple(hom.evaluation(generator, x) for generator in generators)
print("ev_x images:", [tuple(int(image.to_vector()(l)) for l in R.module_generating_set()) for image in images],
      flush=True)

for label, build in (
    ("dict of images", lambda: modules.Mor(hom, R)(dict(zip(labels, images)))),
    ("list of images", lambda: modules.Mor(hom, R)(list(images))),
):
    try:
        morphism = build()
        print(f"[{label}] -> {type(morphism).__name__}", flush=True)
        kernel = morphism.kernel()
        print(f"    kernel: {type(kernel).__name__}", flush=True)
        print(f"    surface: {[m for m in ('module_generating_set','module_generator','module_rank','inclusion','domain') if hasattr(kernel, m)]}", flush=True)
        if hasattr(kernel, "module_rank"):
            print(f"    kernel module_rank: {kernel.module_rank()}", flush=True)
        if hasattr(kernel, "module_generating_set"):
            klabels = tuple(kernel.module_generating_set())
            print(f"    kernel generators: {klabels}", flush=True)
            inclusion = kernel.inclusion()
            for kl in klabels:
                element_ = inclusion(kernel.module_generator(kl))
                image = hom.evaluation(element_, x)
                print(f"      {kl} -> ev_x image {tuple(int(image.to_vector()(l)) for l in R.module_generating_set())}",
                      flush=True)
        break
    except Exception as error:
        import traceback
        print(f"[{label}] raised {type(error).__name__}: {str(error)[:200]}", flush=True)
