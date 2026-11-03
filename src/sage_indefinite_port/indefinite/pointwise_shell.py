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

The orbit of `v` in the **primitive isotropic shell** (`q(v) = 0`, content one, so `a m + c d = 0`) is therefore:

* `m != 0`: `a = -cd/m`, and the orbit is the class of `(c, d)` in `(Z/m)^2` under `m Z^2` together with the four
  symmetries `(c,d), (-c,-d), (d,c), (-d,-c)`;
* `m = 0`: `c d = 0`, and `a` can be moved by `rho c + gam d`, i.e. by the ideal `gcd(c, d) Z`; the orbit is the
  class of `(c,d)` under the same four symmetries together with `a mod gcd(c, d)` (and for `c = d = 0`, `a`
  itself -- `e1` and `-e1` are different orbits).

This module decides that relation exactly and returns an explicit element of `Gamma`.
"""

from __future__ import annotations

from math import gcd

from sage.all import ZZ, matrix

from sage_indefinite_port.readiness import UnfinishedCapability
from sage_indefinite_port.indefinite.recursive import _isometry_from_row_matrix

__all__ = [
    "U_PLUS_U_GRAM",
    "gamma_action",
    "gamma_class",
    "gamma_equivalent",
    "gamma_witness",
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
    r"""`(a, m, c, d)` of `vector` in the distinguished basis, with the shell hypotheses checked.

    The shell is the set of **primitive isotropic** vectors: `q(v) = 0` and content one. Both are checked here
    because the classification below uses `m != 0`, `a = -cd/m` and `gcd` of the radical-free coordinates, none of
    which are meaningful for an imprimitive or non-isotropic vector.
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
        raise UnfinishedCapability(f"{coordinates} is not isotropic; this decision covers the isotropic shell only")
    if int(vector.content()) != 1:
        raise UnfinishedCapability(
            f"{coordinates} is not primitive (content {int(vector.content())}); the classification below is for "
            "the primitive shell"
        )
    if a * m + c * d != 0:
        raise ArithmeticError(f"an isotropic vector in this basis must satisfy a*m + c*d = 0, got {coordinates}")
    return coordinates


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


def gamma_class(lattice, vector) -> tuple:
    r"""The `Gamma`-orbit class of a primitive isotropic `vector`, as a hashable key.

    Two vectors are `Gamma`-equivalent exactly when their keys agree (see the module docstring for the proof).
    """
    a, m, c, d = shell_coordinates(lattice, vector)
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
    unipotent element solving a linear divisibility condition -- then verified to fix `e1` and to carry `source`
    onto `target`.
    """
    source_class = gamma_class(lattice, source)
    if source_class != gamma_class(lattice, target):
        return None
    a, m, c, d = shell_coordinates(lattice, source)
    target_coordinates = shell_coordinates(lattice, target)
    if m != 0:
        for symmetry in _SYMMETRIES:
            x, y = symmetry(c, d)
            if (x - target_coordinates[2]) % m == 0 and (y - target_coordinates[3]) % m == 0:
                gam = (x - target_coordinates[2]) // m
                rho = (y - target_coordinates[3]) // m
                return _certified(lattice, source, target, symmetry, rho, gam)
        raise ArithmeticError("the classes agree but no symmetry matches modulo m; the classification is wrong")
    g = gcd(abs(c), abs(d))
    for symmetry in _SYMMETRIES:
        x, y = symmetry(c, d)
        if (x, y) != (target_coordinates[2], target_coordinates[3]):
            continue
        delta = target_coordinates[0] - a
        if g == 0:
            if delta == 0:
                return _certified(lattice, source, target, symmetry, 0, 0)
            continue
        rho, gam = _solve(x, y, delta)
        return _certified(lattice, source, target, symmetry, rho, gam)
    raise ArithmeticError(
        "the classes agree but no symmetry carries (c, d) onto the target's; the classification is wrong"
    )


def gamma_equivalent(lattice, source, target) -> bool:
    """Whether two primitive isotropic vectors lie in one `Gamma`-orbit, exactly."""
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
