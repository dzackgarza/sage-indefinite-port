r"""Exact `Gamma`-equivalence on the primitive isotropic shell of `L = U + U`, `Gamma = Stab(e1)` pointwise.

Issues #53 and #62. This is the port-owned half: the actual indefinite `Gamma`-action, the exact equivalence
decision where it holds, and the constructive witnesses. The generic quotient/enumerator framework is research's
and is not reimplemented here -- but the equivalence decision itself belongs to this port, and #54 is the
designated arithmetic source for it in general.

## The specimen and the theorem

`L = U + U`, basis `(e1, f1, e2, f2)`, `b(e_i, f_i) = 1`, the two copies orthogonal,
`q(x) = 2(x1 x2 + x3 x4)`. `Gamma = {g in O(L) : g(e1) = e1}` **pointwise**.

`Q = e1-perp / <e1>` is a rank-two nondegenerate lattice isometric to `U`, and

> **Theorem.** `Gamma = N ⋊ Theta`, where `N` is the kernel of `Gamma -> O(Q)` and `Theta` maps isomorphically
> onto `O(Q)`. Explicitly `N = {u(rho, gam) : rho, gam in Z}` and `Theta = {1, sigma, tau, sigma tau}`, with

    u(rho, gam): e1 |-> e1
                 f1 |-> f1 - gam*rho*e1 - gam*e2 - rho*f2
                 e2 |-> e2 + rho*e1
                 f2 |-> f2 + gam*e1
    sigma:       e1 |-> e1,  f1 |-> f1,  e2 |-> -e2,  f2 |-> -f2
    tau:         e1 |-> e1,  f1 |-> f1,  e2 <-> f2

On `v = a e1 + m f1 + c e2 + d f2` they act by

    u(rho, gam):  (c, d) |-> (c - m*gam, d - m*rho),  a |-> a + c*rho + d*gam - m*gam*rho,  m fixed
    sigma:        (c, d) |-> (-c, -d),                a |-> a
    tau:          (c, d) |-> (d, c),                  a |-> a

*Proof.* `Gamma` fixes `e1`, hence preserves `e1-perp` and the line `<e1>` inside it, so it acts on `Q`. For
`g` in the kernel, `g(e1) = e1`, `g(e2) = e2 + rho e1` and `g(f2) = f2 + gam e1`; writing `g(f1) = lam f1 + z`
and imposing the three conditions `b(g e1, g f1) = 1`, `b(g e2, g f1) = b(e2, f1) = 0`, `b(g f2, g f1) = 0` and
`q(g f1) = 0` gives `lam = 1`, `z = -gam rho e1 - gam e2 - rho f2`: exactly `u(rho, gam)`. Conversely each
`u(rho, gam)` is an isometry fixing `e1` and acting trivially on `Q`. So `N = {u(rho, gam)}`.
`O(U) = {I, -I, swap, -swap}` (from `M [[0,1],[1,0]] M^T = [[0,1],[1,0]]` for integral `M`), and `sigma`, `tau`
realize `-I` and `swap`, so `Theta -> O(Q)` is onto. ∎

The orbit of a **primitive** isotropic `v` (`q(v) = 0`, content one, so `a m + c d = 0`) is therefore:

* `m != 0`: `a = -cd/m`, and the orbit is the class of `(c, d)` in `(Z/m)^2` under `m Z^2` together with the four
  symmetries `(c,d), (-c,-d), (d,c), (-d,-c)`;
* `m = 0`: `c d = 0`, and `a` can be moved by `rho c + gam d`, i.e. by the ideal `gcd(c, d) Z`; the orbit is the
  class of `(c,d)` under the same four symmetries together with `a mod gcd(c, d)` (and for `c = d = 0`, `a`
  itself -- `e1` and `-e1` are different orbits).

**Content.** Every nonzero isotropic `v` is `t w` with `w` primitive isotropic and `t = content(v) >= 1`. An
isometry of `L` is a unimodular change of the basis, so it preserves `content`; it is also `Z`-linear, so
`g(t w) = t g(w)`. Hence `t w` and `t' w'` are `Gamma`-equivalent exactly when `t = t'` and `w, w'` are, and the
`Gamma`-orbits on the **whole isotropic cone** are `{(t, [w]) : t >= 1, [w] a primitive Gamma-orbit}`. The class
key therefore carries `t`, and a witness for a scaled pair is the primitive witness -- certified against the actual
endpoints, so the scaled statement is checked rather than inferred.

This module decides that relation exactly on the whole isotropic cone and returns an explicit element of `Gamma`.
Non-isotropic vectors are outside it and are refused: `q` is the only invariant supplied there and it does not
classify.
"""

from __future__ import annotations

from math import gcd

from sage.all import ZZ, matrix, vector

from sage_indefinite_port.readiness import UnfinishedCapability
from sage_indefinite_port.indefinite.recursive import _isometry_from_row_matrix

__all__ = [
    "U_PLUS_U_GRAM",
    "gamma_action",
    "gamma_class",
    "gamma_decompose",
    "gamma_element",
    "gamma_equivalent",
    "gamma_witness",
    "isotropic_content",
    "pointwise_shell_action",
    "pointwise_shell_group",
    "pointwise_shell_quotient",
    "shell_coordinates",
]

#: `U + U` in the basis `(e1, f1, e2, f2)` with `b(e_i, f_i) = 1`.
U_PLUS_U_GRAM = ((0, 1, 0, 0), (1, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 0))

_SYMMETRIES = (
    lambda c, d: (c, d),
    lambda c, d: (-c, -d),
    lambda c, d: (d, c),
    lambda c, d: (-d, -c),
)


def shell_coordinates(lattice, vector) -> tuple[int, int, int, int]:
    r"""`(a, m, c, d)` of an isotropic vector in the distinguished basis, with the presentation checked.

    The isotropic condition `q(v) = 0` is checked here because the classification below uses `a m + c d = 0`, which
    is not meaningful for a non-isotropic vector. Content is **not** restricted: imprimitive isotropic vectors are
    classified by splitting off content (see :func:`gamma_class` and the module docstring), so the only refusal here
    is the genuine boundary of the decision.
    """
    labels = tuple(lattice.module_generating_set())
    if len(labels) != 4:
        raise ValueError("this decision is for the rank-four specimen L = U + U")
    if tuple(tuple(int(entry) for entry in row) for row in _gram_rows(lattice)) != U_PLUS_U_GRAM:
        raise ValueError("this decision is for the U + U presentation in the basis (e1, f1, e2, f2)")
    raw = vector.to_vector()
    coordinates = tuple(int(raw(label)) for label in labels)
    a, m, c, d = coordinates
    if int(vector.q()) != 0:
        raise UnfinishedCapability(
            f"{coordinates} is not isotropic; this decision covers the isotropic cone only, and `q` alone does not "
            "classify the non-isotropic vectors"
        )
    if a * m + c * d != 0:
        raise ArithmeticError(f"an isotropic vector in this basis must satisfy a*m + c*d = 0, got {coordinates}")
    return coordinates


def isotropic_content(coordinates) -> int:
    r"""`content` of `(a, m, c, d)` -- the gcd of its entries, which is `0` only for the zero vector.

    This is the `t` of the module's content corollary. It is a `Gamma`-invariant: an isometry of `L` is a
    unimodular change of the basis, so it cannot change the gcd of a vector's coordinates.
    """
    a, m, c, d = coordinates
    return gcd(gcd(abs(a), abs(m)), gcd(abs(c), abs(d)))


def _gram_rows(lattice):
    from sage_indefinite_port.indefinite.groupoid_cache import _gram_in_label_order

    gram = _gram_in_label_order(lattice)
    return [[int(gram[row, column]) for column in range(gram.ncols())] for row in range(gram.nrows())]


def _action_coordinates(kind, parameters, coordinates):
    """The action of a named `Gamma` element on `(a, m, c, d)`, in the derived formulas."""
    a, m, c, d = coordinates
    match kind:
        case "u":
            rho, gam = parameters
            return (a + c * rho + d * gam - m * gam * rho, m, c - m * gam, d - m * rho)
        case "sigma":
            return (a, m, -c, -d)
        case "tau":
            return (a, m, d, c)
    raise ValueError(f"unknown Gamma element kind {kind}")


def gamma_action(lattice, kind: str, *parameters):
    r"""The `Gamma` element named by `kind`, as a validated isometry of `lattice`.

    Built from its images of the four basis vectors and passed through the port's checked row-matrix boundary, so
    the Gram identity is verified rather than trusted; every element produced here also fixes `e1`, which is
    asserted.
    """
    labels = tuple(lattice.module_generating_set())
    images = []
    for index in range(4):
        basis_coordinates = tuple(1 if position == index else 0 for position in range(4))
        moved = _action_coordinates(kind, parameters, basis_coordinates)
        images.append([moved[position] for position in range(4)])
    isometry = _isometry_from_row_matrix(lattice, lattice, matrix(ZZ, images))
    e1 = lattice.module_generator(labels[0])
    if isometry(e1) != e1:
        raise ArithmeticError(f"the constructed Gamma element {kind}{parameters} does not fix e1")
    return isometry


#: The proved generating set of `Gamma`, in the order `(sigma, tau, u(1,0), u(0,1))`.
_PROVED_GENERATORS = (("sigma", ()), ("tau", ()), ("u", (1, 0)), ("u", (0, 1)))


def pointwise_shell_group(lattice):
    r"""The pointwise stabilizer `Gamma = {g in O(L) : g(e1) = e1}`, as the provider's subgroup parent.

    Generated by the four elements `sigma`, `tau`, `u(1,0)`, `u(0,1)`, and this set is **proved** to equal `Gamma`,
    not merely sampled:

    * `N = {u(rho, gam)}` is closed and abelian-by-the-identity `u(rho1,gam1) u(rho2,gam2) = u(rho1+rho2, gam1+gam2)`
      -- checked on the four basis images, where the `e2` and `f2` coordinates add and the `f1` coordinate gives
      `-(gam1+gam2)(rho1+rho2) = -gam1rho1 - gam1rho2 - gam2rho1 - gam2rho2`, which is exactly the `f1`
      coordinate of `u(rho1+rho2, gam1+gam2)`. So `N = <u(1,0), u(0,1)>`, isomorphic to `Z^2`;
    * `Theta = {1, sigma, tau, sigma tau} = <sigma, tau>` is the Klein four group, and maps onto `O(Q)`;
    * `Gamma = N ⋊ Theta` is the module theorem, so `<sigma, tau, u(1,0), u(0,1)> = Gamma`.

    Using this set rather than the six generators returned by `pointed_automorphism_generators` is deliberate: those
    are the source's claim, whose completeness this port does not re-derive, whereas this set's completeness follows
    from the decomposition above.
    """
    elements = tuple(gamma_action(lattice, kind, *parameters) for kind, parameters in _PROVED_GENERATORS)
    e1 = lattice.module_generator(tuple(lattice.module_generating_set())[0])
    for element in elements:
        if element(e1) != e1:
            raise ArithmeticError("a proved generator of Gamma does not fix e1")
    return lattice.O().subgroup(elements)


def gamma_element(lattice, symmetry_index: int, rho: int, gam: int):
    r"""The element `u(rho, gam) . S` with `S` the symmetry of that index, applied **after** nothing -- `S` first.

    Public so callers and tests can build the mixed elements that decompose to `(symmetry_index, rho, gam)` without
    reaching into the private composition helper. Validated the same way: an isometry through the checked
    row-matrix boundary, fixing `e1`.
    """
    index = int(symmetry_index)
    if not 0 <= index < len(_SYMMETRIES):
        raise ValueError(f"symmetry index must be in 0..{len(_SYMMETRIES) - 1}, got {symmetry_index}")
    return _compose_symmetry_then_unipotent(lattice, _SYMMETRIES[index], int(rho), int(gam))


def gamma_decompose(lattice, isometry):
    r"""Decompose an isometry of `L` as `u(rho, gam) . S`, or `None` if it is not in `Gamma`.

    This is the exact membership decision the decomposition gives, and it replaces comparing sampled images. Given
    `g`, read `g(e2)` and `g(f2)`: their images in `Q = e1-perp/<e1>` must be `S(1,0)` and `S(0,1)` for one of the
    four symmetries, which determines `S`; the `e1` coefficients give `rho` and `gam`; then `g` is in `Gamma` exactly
    when it agrees with the constructed `u(rho, gam) . S` on **all four** basis vectors. Every retained generator is
    decided this way rather than spot-checked.

    Returns `(symmetry_index, rho, gam)`.
    """
    labels = tuple(lattice.module_generating_set())
    e1 = lattice.module_generator(labels[0])
    if isometry(e1) != e1:
        return None

    def coordinates(element):
        raw = element.to_vector()
        return tuple(int(raw(label)) for label in labels)

    second = coordinates(isometry(lattice.module_generator(labels[2])))   # image of e2
    fourth = coordinates(isometry(lattice.module_generator(labels[3])))   # image of f2
    q_image_of_e2 = (second[2], second[3])
    q_image_of_f2 = (fourth[2], fourth[3])
    for index, symmetry in enumerate(_SYMMETRIES):
        if (symmetry(1, 0), symmetry(0, 1)) != (q_image_of_e2, q_image_of_f2):
            continue
        # The element is `u(rho, gam) . S` with `S` applied FIRST, so `g(e2) = u(S(e2))`: the `e1` coefficient of
        # `g(e2)` is `<S(e2), (rho, gam)>`, not `rho`. For `S = sigma` it is `-rho`, which an extraction reading
        # `rho` off that coefficient gets wrong -- it rebuilds `u(-1,0) sigma` for `u(1,0) sigma` and reports the
        # element as outside `Gamma`. Solve the actual 2x2 system instead: with `A` the matrix whose rows are
        # `S(1,0)` and `S(0,1)`, `A . (rho, gam)^T = (coefficient of e1 in g(e2), coefficient of e1 in g(f2))`.
        system = matrix(ZZ, [[q_image_of_e2[0], q_image_of_e2[1]], [q_image_of_f2[0], q_image_of_f2[1]]])
        if abs(int(system.det())) != 1:
            return None
        solution = system.inverse().change_ring(ZZ) * vector(ZZ, [second[0], fourth[0]])
        rho, gam = int(solution[0]), int(solution[1])
        candidate = _compose_symmetry_then_unipotent(lattice, symmetry, rho, gam)
        for position in range(4):
            if candidate(lattice.module_generator(labels[position])) != isometry(
                lattice.module_generator(labels[position])
            ):
                return None
        return (index, rho, gam)
    return None


def pointwise_shell_action(lattice):
    r"""The `Gamma`-set on the whole isotropic cone `L[0]`, through the provider's own action constructor.

    This is the maintained machinery, not a second framework: `GObjects(Gamma, Sets()).on_set(points, action,
    orbit_relation)` is what :meth:`EichlerRecursion.marked_orbit_action` uses for `O(L)`, and here the group is the
    pointwise stabilizer and the relation is :func:`gamma_equivalent`, which is exact on the whole cone -- imprimitive
    vectors included, by the content corollary. The point set is a predicate cut, so no point is ever materialised;
    `orbits()` keeps returning the provider's lazy orbit set, and its `is_finite()` stays a proposition because the
    quotient here is infinite.

    The zero vector is a member and a fixed point; it has its own class, so including it costs nothing and the point
    set is then exactly the cone `q = 0` rather than the cone minus a point.
    """
    from sage_extensions.all import GObjects, Sets

    cone = Sets().condition_set(lattice, lambda vector: int(vector.q()) == 0)
    return GObjects(pointwise_shell_group(lattice), Sets()).on_set(
        cone,
        lambda isometry, point: isometry(point),
        orbit_relation=lambda left, right: gamma_equivalent(lattice, left, right),
    )


def pointwise_shell_quotient(lattice):
    r"""`L[0] / Gamma` as the provider's orbit set, built on :func:`pointwise_shell_action`."""
    return pointwise_shell_action(lattice).orbits()


def gamma_class(lattice, vector) -> tuple:
    r"""The `Gamma`-orbit class of an isotropic `vector`, as a hashable key -- imprimitive vectors included.

    Two isotropic vectors are `Gamma`-equivalent exactly when their keys agree. For content `t > 1` the key is
    `("content", t, primitive_key)`, by the module's content corollary: an isometry is a unimodular change of basis
    and is `Z`-linear, so it preserves `content` and sends `t w` to `t g(w)`. For a primitive vector the key is the
    class of `(c, d)` described in the module docstring. The zero vector has content `0`, gets the unique key
    `("linear", (0, 0), 0)`, and is its own orbit -- every isometry fixes it.
    """
    coordinates = shell_coordinates(lattice, vector)
    a, m, c, d = coordinates
    t = isotropic_content(coordinates)
    if t > 1:
        return ("content", t, _primitive_class(a // t, m // t, c // t, d // t))
    return _primitive_class(a, m, c, d)


def _primitive_class(a, m, c, d) -> tuple:
    r"""The class key of a **primitive** isotropic `(a, m, c, d)`."""
    if m != 0:
        # SIGNED m is the key. `b(e1, .)` is invariant under every pointwise stabilizer element, so `m` and `-m`
        # are different classes; keying on `abs(m)` merged `f1` with `-f1`, which the invariant separates and no
        # Gamma element can do. Only the residue lattice uses `|m|`, because `u(rho, gam)` shifts `(c, d)` by
        # `m Z^2` and `m Z^2 = |m| Z^2`.
        residue = abs(m)
        candidates = [tuple(value % residue for value in symmetry(c, d)) for symmetry in _SYMMETRIES]
        return ("affine", m, min(candidates))
    g = gcd(abs(c), abs(d))
    candidates = [tuple(symmetry(c, d)) for symmetry in _SYMMETRIES]
    return ("linear", min(candidates), a % g if g else a)


def gamma_witness(lattice, source, target):
    r"""An explicit element of `Gamma` carrying `source` onto `target`, or `None` if there is none.

    `None` is a **decision**, not a failure: it means the two classes differ, which the module's theorem proves
    is the same as being in different `Gamma`-orbits. The witness is constructed -- a symmetry followed by a
    unipotent element solving a linear divisibility condition -- then verified to fix `e1` and to carry the
    *actual* `source` onto the *actual* `target`. Factoring content off first is therefore not a shortcut that
    escapes checking: an imprimitive pair is certified at the scaled endpoints.
    """
    if gamma_class(lattice, source) != gamma_class(lattice, target):
        return None
    source_coordinates = shell_coordinates(lattice, source)
    target_coordinates = shell_coordinates(lattice, target)
    t = isotropic_content(source_coordinates)
    if t > 1:
        source_coordinates = tuple(value // t for value in source_coordinates)
        target_coordinates = tuple(value // t for value in target_coordinates)
    symmetry, rho, gam = _primitive_witness_parameters(source_coordinates, target_coordinates)
    return _certified(lattice, source, target, symmetry, rho, gam)


def _primitive_witness_parameters(source_coordinates, target_coordinates):
    r"""`(symmetry, rho, gam)` carrying one **primitive** class onto another of the same class.

    Split out of :func:`gamma_witness` so content is factored off once: the answer for `t w -> t w'` is the answer
    for `w -> w'`, and the caller certifies it against the scaled endpoints.
    """
    a, m, c, d = source_coordinates
    target_a, _, target_c, target_d = target_coordinates
    if m != 0:
        for symmetry in _SYMMETRIES:
            x, y = symmetry(c, d)
            if (x - target_c) % m == 0 and (y - target_d) % m == 0:
                gam = (x - target_c) // m
                rho = (y - target_d) // m
                return symmetry, rho, gam
        raise ArithmeticError("the classes agree but no symmetry matches modulo m; the classification is wrong")
    g = gcd(abs(c), abs(d))
    for symmetry in _SYMMETRIES:
        x, y = symmetry(c, d)
        if (x, y) != (target_c, target_d):
            continue
        delta = target_a - a
        if g == 0:
            if delta == 0:
                return symmetry, 0, 0
            continue
        rho, gam = _solve(x, y, delta)
        return symmetry, rho, gam
    raise ArithmeticError(
        "the classes agree but no symmetry carries (c, d) onto the target's; the classification is wrong"
    )


def gamma_equivalent(lattice, source, target) -> bool:
    r"""Whether two isotropic vectors lie in one `Gamma`-orbit, exactly."""
    return gamma_class(lattice, source) == gamma_class(lattice, target)


def _solve(x, y, delta):
    r"""`(rho, gam)` with `rho x + gam y = delta`, for `gcd(x, y) | delta`.

    Through Smith normal form of the row `(x0, y0)` of the primitive part: `(x0, y0) V = e_1` gives
    `(rho, gam) = delta0 * (first column of V)`.
    """
    g = gcd(abs(x), abs(y))
    if delta % g:
        raise ArithmeticError(f"{delta} is not in the ideal generated by {x} and {y}")
    x0, y0, delta0 = x // g, y // g, delta // g
    _, _, change = matrix(ZZ, [[x0, y0]]).smith_form()
    if abs(int(change.det())) != 1:
        raise ArithmeticError("Smith form must return a unimodular change")
    first_column = tuple(int(change[row, 0]) for row in range(2))
    rho, gam = delta0 * first_column[0], delta0 * first_column[1]
    if rho * x + gam * y != delta:
        raise ArithmeticError("the Smith-based solution does not satisfy the divisibility equation")
    return rho, gam


def _certified(lattice, source, target, symmetry, rho, gam):
    r"""Build `u(rho, gam) . S` and **check** that it lies in `Gamma` and carries `source` onto `target`.

    Every success path goes through here, so a branch cannot return a constructed element whose endpoint was never
    compared with the actual target. That check is what the previous revision of the `m != 0` branch lacked: it
    returned a composite for a pair whose classes only agreed because the key used `abs(m)`, and the composite was
    the identity.
    """
    candidate = _compose_symmetry_then_unipotent(lattice, symmetry, rho, gam)
    labels = tuple(lattice.module_generating_set())
    e1 = lattice.module_generator(labels[0])
    if candidate(e1) != e1:
        raise ArithmeticError("the constructed Gamma element does not fix e1")
    if candidate(source) != target:
        raise ArithmeticError(
            "the constructed Gamma element does not carry the source onto the target; the classification or the "
            "witness construction is wrong"
        )
    return candidate


def _compose_symmetry_then_unipotent(lattice, symmetry, rho, gam):
    """The composite `u(rho, gam) . S`, built as an isometry and verified on `e1`."""
    rows = []
    for index in range(4):
        basis = tuple(1 if position == index else 0 for position in range(4))
        rows.append(list(_action_coordinates("u", (rho, gam), _action_coordinates_symmetry(symmetry, basis))))
    isometry = _isometry_from_row_matrix(lattice, lattice, matrix(ZZ, rows))
    labels = tuple(lattice.module_generating_set())
    if isometry(lattice.module_generator(labels[0])) != lattice.module_generator(labels[0]):
        raise ArithmeticError("the composite does not fix e1")
    return isometry


def _action_coordinates_symmetry(symmetry, coordinates):
    a, m, c, d = coordinates
    x, y = symmetry(c, d)
    return (a, m, x, y)
