"""Does the retained reduction put the radical last, and can a vector be transported into it?"""
from tests.acceptance.consumer import element, lattice
from sage_indefinite_port.indefinite.indefinite_reduction import indefinite_reduced_presentation
from sage_indefinite_port.indefinite.recursive import _gram_in_label_order

for label, gram, coords in (
    ("radical last", [[2, 0], [0, 0]], [0, 1]),
    ("radical FIRST (permuted)", [[0, 0], [0, 2]], [1, 0]),
    ("mixed", [[2, 0, 0], [0, 0, 0], [0, 0, 2]], [0, 1, 0]),
    ("purely radical", [[0, 0], [0, 0]], [1, 0]),
):
    L = lattice(gram)
    v = element(L, coords)
    try:
        reduction = indefinite_reduced_presentation(L)
        reduced = reduction.reduced
        G = _gram_in_label_order(reduced)
        print(f"[{label}] reduced gram {G.rows()} rank {G.rank()}", flush=True)
        inverse = reduction.isometry.inverse()
        image = inverse(v)
        print(f"    vector transported into reduced coordinates: "
              f"{tuple(int(image.to_vector()(l)) for l in reduced.module_generating_set())}", flush=True)
        back = reduction.isometry(image)
        print(f"    transported back equals original: {back == v}", flush=True)
    except Exception as error:
        print(f"[{label}] raised {type(error).__name__}: {str(error)[:140]}", flush=True)
