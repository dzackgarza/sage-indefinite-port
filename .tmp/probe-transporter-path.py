"""Is the BOUND note's "finite-representation transporter still fails" still true at 42d7361c?"""
from sage_extensions.categories.rational_integral_stabilizers import IntegralStructureAction
from sage_indefinite_port.indefinite.recursive import EichlerRecursion
from tests.acceptance.consumer import element, lattice

original = IntegralStructureAction.transporter
calls = []

def counted(self, *args, **kwargs):
    try:
        result = original(self, *args, **kwargs)
        calls.append(("ok", result is not None))
        return result
    except Exception as error:
        calls.append(("raised", type(error).__name__, str(error)[:120]))
        raise

IntegralStructureAction.transporter = counted
RECURSION = EichlerRecursion()

# a divisibility-four isotropic pointed case: the port route that has a target correction step
specimens = (
    ("div-4 isotropic U(4)+U(4)", [[0, 4, 0, 0], [4, 0, 0, 0], [0, 0, 0, 4], [0, 0, 4, 0]], [1, 0, 0, 0], [0, 0, 1, 0]),
)
for label, gram, source_coords, target_coords in specimens:
    L = lattice(gram)
    source, target = element(L, source_coords), element(L, target_coords)
    try:
        witness = RECURSION.pointed_isometry(source, target)
        print(f"[{label}] witness: {witness is not None}", flush=True)
    except Exception as error:
        print(f"[{label}] raised {type(error).__name__}: {str(error)[:150]}", flush=True)
    print(f"    transporter calls: {calls}", flush=True)
    calls.clear()
