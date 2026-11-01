"""How does an owned permutation element take more than one cycle?"""
from sage_extensions.all import OwnedGroups

symmetric = OwnedGroups().S(4)
attempts = [
    ("list of cycles", lambda: symmetric([(1, 3), (2, 4)])),
    ("tuple of cycles", lambda: symmetric(((1, 3), (2, 4)))),
    ("single flat cycle", lambda: symmetric((1, 3, 2, 4))),
]
for label, call in attempts:
    try:
        element = call()
        print(f"{label:20s} -> {element}  order={element.multiplicative_order() if hasattr(element,'multiplicative_order') else '?'}")
    except Exception as error:
        print(f"{label:20s} raised {type(error).__name__}: {str(error)[:100]}")

# does a flat-cycle-only constructor still generate the same group?
try:
    flat = symmetric((1, 3, 2, 4))
    print("flat cycle subgroup order:", symmetric.subgroup((flat,)).cardinality())
except Exception as error:
    print("subgroup of flat cycle raised:", type(error).__name__, str(error)[:100])
