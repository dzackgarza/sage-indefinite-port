r"""Primitive isotropic lines, planes and flags of an indefinite lattice up to isometry (map §13).

The induction of references/port-design.md, "Isotropic induction": for a line class Zu with
K = u^perp/u and H the image of Aut((L, Zu)) in O(K), the H-orbits of isotropic lines of K are
the double cosets H \ O(K) / Aut((K, Zx)) over the O(K)-classes of isotropic lines Zx, and each
gives the flag (Zu, J) with J the saturation of the lift of the line through u. Those are the
flag classes with first term Zu; planes are the flags' second terms deduplicated by Isom.

Deciding Isom of flags reduces to K: carry the lines onto each other, then move the image of the
plane by an element of H, which lifts to Aut((L, Zu)). Deciding Isom of planes runs over the
lines of the target plane up to its automorphism group: the double cosets
Gamma \ GL(J) / Stab(line), with Gamma the image of Aut((L, J)) in GL(J).

Parabolic subgroups, Levi actions, double cosets and GL are the preamble's; the recursion
supplies O(K) and the isotropic classes of K.
"""

from collections.abc import Iterator

from dzack_research.preamble.all import ZZ, Cat, Lattices
from dzack_research.preamble.categories.abstract_categories.products import FiniteOrdinalCategory
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from sage.structure.parent import Parent

from sage_indefinite_port.indefinite.lorentzian_cells import LorentzianPerfectComplex, MarkedCellOrbitAlgorithm
from sage_indefinite_port.indefinite.recursive import EichlerRecursion
from sage_indefinite_port.invariants import AttackProfile

ARROWS = Cat().Mor(FiniteOrdinalCategory(2), Lattices(ZZ))
"""Ar(Lattices): a sublattice I of L is the marked lattice (L, I -> L)."""

type Vector = Lattices.ElementMethods
type Flag = tuple[Vector, Parent]
"""A flag Zu < J < L: a primitive isotropic u and the totally isotropic plane J containing it, a subobject of L."""


class UnsupportedRankError(ValueError):
    """Totally isotropic sublattices of rank above 2 are outside the induction (map §13.5)."""


class IsotropicInduction:
    def __init__(self, recursion: EichlerRecursion) -> None:
        self.recursion = recursion

    def isotropic_vector_classes(self, lattice: Lattices.ParentMethods) -> tuple[Vector, ...]:
        r"""One primitive isotropic vector per O(L)-orbit: none for definite L, marked cells for positive index 1, the recursion beyond."""
        profile = AttackProfile.from_lattice(lattice)
        match profile.positive_index:
            case 0:
                return ()
            case 1:
                orbits = MarkedCellOrbitAlgorithm(LorentzianPerfectComplex(lattice, "isotropic")).global_orbits()
                return tuple(orbit[0] for orbit in orbits)
            case _:
                return tuple(vector_class.representative for vector_class in self.recursion.primitive_vector_classes(lattice, 0))

    @staticmethod
    def reduction(lattice: Lattices.ParentMethods, *vectors: Vector) -> Parent:
        r"""The isotropic reduction I^perp/I of the primitive totally isotropic I spanned by the vectors."""
        return lattice.primitive_isotropic_subobject(*vectors).isotropic_reduction()

    @staticmethod
    def plane_vectors(plane: Parent) -> tuple[Vector, ...]:
        inclusion = plane.inclusion()
        return tuple(inclusion(generator) for generator in inclusion.domain().module_generators())

    def line_isometry(self, source: Vector, target: Vector) -> LatticeIsometryMethods | None:
        r"""An isometry carrying Zu onto Zw: one carrying u to w or to -w."""
        isometry = self.recursion.pointed_isometry(source, target)
        return isometry if isometry is not None else self.recursion.pointed_isometry(source, -target)

    @staticmethod
    def quotient_line(reduction: Parent, vectors: tuple[Vector, ...]) -> Parent:
        r"""The line of K = u^perp/u that a plane through u projects to, as a subobject of K."""
        perpendicular = reduction.orthogonal_complement().inclusion()
        projection = reduction.projection()
        images = tuple(image for image in (projection(perpendicular.lift(vector_)) for vector_ in vectors) if image != reduction.zero())
        return reduction.primitive_isotropic_subobject(*images)

    def flags_through(self, lattice: Lattices.ParentMethods, line: Vector) -> Iterator[Flag]:
        r"""One flag (Zu, J) per Aut((L, Zu))-orbit of isotropic planes J through u."""
        reduction = self.reduction(lattice, line)
        quotient = reduction.quotient_lattice()
        levi_image = reduction.levi_action().image()
        orthogonal_group = quotient.O().subgroup(self.recursion.orthogonal_group_generators(quotient))
        perpendicular = reduction.orthogonal_complement().inclusion()
        for isotropic in self.isotropic_vector_classes(quotient):
            line_stabilizer = self.reduction(quotient, isotropic).parabolic_subgroup()
            for coset in orthogonal_group.double_cosets(levi_image, line_stabilizer):
                lifted = perpendicular(reduction.projection().lift(coset(isotropic)))
                yield line, lattice.primitive_isotropic_subobject(line, lifted)

    def flag_classes(self, lattice: Lattices.ParentMethods) -> tuple[Flag, ...]:
        return tuple(flag for line in self.isotropic_vector_classes(lattice) for flag in self.flags_through(lattice, line))

    def plane_classes(self, lattice: Lattices.ParentMethods) -> tuple[Parent, ...]:
        kept: list[Parent] = []
        for _line, plane in self.flag_classes(lattice):
            if all(self.plane_isometry(plane, other) is None for other in kept):
                kept.append(plane)
        return tuple(kept)

    def flag_isometry(self, source: Flag, target: Flag) -> LatticeIsometryMethods | None:
        r"""An isometry carrying the line onto the line and the plane onto the plane, or None."""
        (source_line, source_plane), (target_line, target_plane) = source, target
        carry_line = self.line_isometry(source_line, target_line)
        if carry_line is None:
            return None
        target_lattice = target_line.parent()
        reduction = self.reduction(target_lattice, target_line)
        moved = self.quotient_line(reduction, tuple(carry_line(vector_) for vector_ in self.plane_vectors(source_plane)))
        fixed = self.quotient_line(reduction, self.plane_vectors(target_plane))
        levi = reduction.levi_action()
        isomorphism = ARROWS.Core().Mor(ARROWS.object(moved.inclusion()), ARROWS.object(fixed.inclusion())).an_element_in(levi.image())
        if isomorphism is None:
            return None
        return levi.lift(isomorphism.forward().right()) * carry_line

    def plane_lines(self, plane: Parent) -> tuple[Vector, ...]:
        r"""One primitive vector of J per orbit of Aut((L, J)) on the lines of J."""
        reduction = plane.isotropic_reduction()
        action = reduction.isotropic_action()
        domain = plane.inclusion().domain()
        first = next(iter(domain.module_generators()))
        general_linear = domain.general_linear_group()
        line_stabilizer = general_linear.flag_stabilizer((domain.primitive_isotropic_subobject(first),))
        inclusion = plane.inclusion()
        return tuple(inclusion(coset(first)) for coset in general_linear.double_cosets(action.image(), line_stabilizer))

    def plane_isometry(self, source: Parent, target: Parent) -> LatticeIsometryMethods | None:
        r"""An isometry carrying the plane J onto J', or None: a flag isometry from (Zu, J) to some (Zu', J')."""
        line = self.plane_vectors(source)[0]
        return next(
            (isometry for target_line in self.plane_lines(target) if (isometry := self.flag_isometry((line, source), (target_line, target))) is not None),
            None,
        )

    def line_automorphism_generators(self, lattice: Lattices.ParentMethods, line: Vector) -> tuple[LatticeIsometryMethods, ...]:
        return tuple(self.reduction(lattice, line).parabolic_subgroup().group_generators())

    @staticmethod
    def plane_automorphism_generators(plane: Parent) -> tuple[LatticeIsometryMethods, ...]:
        return tuple(plane.isotropic_reduction().parabolic_subgroup().group_generators())

    def flag_automorphism_generators(self, flag: Flag) -> tuple[LatticeIsometryMethods, ...]:
        r"""Aut of the flag: the elements of Aut((L, Zu)) whose Levi image preserves the line J/u of K."""
        line, plane = flag
        reduction = self.reduction(line.parent(), line)
        quotient_line = self.quotient_line(reduction, self.plane_vectors(plane))
        line_stabilizer = quotient_line.isotropic_reduction().parabolic_subgroup()
        return tuple(reduction.levi_action().preimage_subgroup(line_stabilizer).group_generators())


__all__ = ["IsotropicInduction", "UnsupportedRankError"]
