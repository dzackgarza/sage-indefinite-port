"""Control: is the Isom failure about degeneracy, or about building the two markings in two categories?"""
from tests.acceptance.consumer import (
    element, isomorphisms, lattice, pointed, pointed_category, underlying_isometry,
)


def marking(category, lattice_, vector):
    return category.object(category.base_object().hom([vector], codomain=lattice_))


def probe(label, gram, a, b):
    L = lattice(gram)
    source, target = element(L, a), element(L, b)
    category = pointed_category(int(source.q()))
    same = (marking(category, L, source), marking(category, L, target))
    split = (pointed(L, source), pointed(L, target))
    for how, pair in (("one category", same), ("two category()-calls", split)):
        try:
            empty = isomorphisms(category if how == "one category" else pointed_category(int(source.q())),
                                 pair[0], pair[1]).is_empty()
            print(f"[{label} | {how}] is_empty={empty}", flush=True)
        except Exception as error:
            print(f"[{label} | {how}] raised {type(error).__name__}: {str(error)[:110]}", flush=True)


probe("degenerate <2>+0 (1,0)->(1,1)", [[2, 0], [0, 0]], [1, 0], [1, 1])
probe("nondegenerate U (1,0)->(0,1)", [[0, 1], [1, 0]], [1, 0], [0, 1])
