r"""Congruence subgroups of the Picard groups, by Schreier's lemma (#59).

For generators of a group `G` of matrices over an order with `Z`-coordinates and a
subgroup `S` of the finite image modulo `m`, the preimage of `S` has finite index. Right
cosets are keyed by their finite image, and Schreier generators generate the preimage.
"""

from collections.abc import Callable, Sequence



def reduction(modulus: int) -> Callable:
    r"""Return the map sending a matrix to the tuple of its `Z`-coordinates modulo ``modulus``."""

    def reduce(element):
        return tuple(int(part) % modulus for entry in element.list() for part in entry.vector())

    return reduce


def finite_image(generators: Sequence, reduce: Callable, one) -> dict:
    r"""Return the finite image of the group, as a map from reductions to group elements."""
    image = {reduce(one): one}
    frontier = [one]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            product = current * generator
            key = reduce(product)
            if key not in image:
                image[key] = product
                frontier.append(product)
    return image


def congruence_generators(generators: Sequence, reduce: Callable, one, allowed: Sequence) -> tuple:
    r"""Return ``(index, schreier)`` for the preimage of the subgroup with lifts ``allowed``."""

    def coset_key(element):
        return min(reduce(lift * element) for lift in allowed)

    transversal = {coset_key(one): one}
    frontier = [one]
    schreier = {}
    while frontier:
        current = frontier.pop()
        for generator in generators:
            product = current * generator
            key = coset_key(product)
            if key not in transversal:
                transversal[key] = product
                frontier.append(product)
                continue
            element = product * ~transversal[key]
            if element != one:
                schreier[str(element)] = element
    return len(transversal), tuple(schreier.values())


__all__ = ["congruence_generators", "finite_image", "reduction"]
