r"""`O(M)` of an anisotropic binary `O_K`-lattice `M = O_K^2` with Gram `G` (#54).

Torus. With `S = [[0, 1], [-1, 0]]`, `J = G^{-1} S` satisfies `G J = S` (so `J^t G = -G J`) and `J^2 = -\det(G)^{-1}` (`S X S = -\mathrm{adj}\, X`
for symmetric `X`). Put `\delta = -\det G`, `L = K(w)`, `w^2 = \delta`, and `\iota(p + q w) = p + q \delta J`; then
`\iota(\alpha)^t G \iota(\alpha) = N_{L/K}(\alpha) G`, so `SO(M \otimes K) = \iota(L^1)` (`M` anisotropic means `\delta` is not a square). `SO(M)` is the
subgroup of `\iota(O_L^1)` with integral matrices, of finite index: `O_L^1` is the kernel of `N_{L/K}` on `O_L^\times`, computed from the
logarithms of the relative norms of the generators of `O_L^\times` (an integer kernel), and the integral subgroup by Schreier generators.

Improper coset. `\sigma`, the reflection in the first basis vector, is a rational improper automorph, and the improper elements of
`O(M \otimes K)` are `\sigma \iota(t)`, `t \in L^1`. Two integral ones differ by `SO(M)`, so it suffices to find one or prove there is none. If
`\sigma \iota(t)` is integral then `t = \sigma \cdot` (integral) has `D t \in O_K[w] \subseteq O_L` for an explicit integer `D`, and `N(D t) = D^2`, so
`(D t)` is a principal ideal dividing `D^2 O_L` of relative norm `D^2 O_K`. For each such ideal `\beta O_L`, `t = \beta \varepsilon \eta / D` with
one unit `\varepsilon` making the norm `D^2` and `\eta \in O_L^1` modulo `SO(M)`, a finite transversal; integrality is tested on each.
"""

from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage_extensions.algebras.isometry_commutant import _integer_solution
from sage_extensions.categories.lattice_centralizers import finite_preimage_generators


def _integral(g):
    return all(c.is_integral() for c in g.list())


def _lcm_denominator(g):
    result = ZZ(1)
    for c in g.list():
        result = result.lcm(ZZ(c.denominator()))
    return result


class BinaryOrthogonalGroup:
    def __init__(self, field, gram):
        self.field, self.gram = field, matrix(field, gram)
        assert self.gram.nrows() == 2 and self.gram == self.gram.transpose()
        self.delta = -self.gram.det()
        assert not self.delta.is_square(), "M is isotropic"
        S = matrix(field, [[0, 1], [-1, 0]])
        self.J = self.gram.inverse() * S
        one = identity_matrix(field, 2)
        assert self.J * self.J == (1 / self.delta) * one
        self.extension = field.extension(PolynomialRing(field, "x")([-self.delta, 0, 1]), "w")

    def iota(self, alpha):
        p, q = self.extension(alpha).list()
        return p * identity_matrix(self.field, 2) + q * self.delta * self.J

    def _norm_logs(self, units):
        K_units = self.field.unit_group()
        return [list(K_units.log(self.field(u.relative_norm()))) for u in units], K_units

    def norm_one_units(self):
        r"""Return generators of `O_L^1` (module docstring)."""
        units = list(self.extension.unit_group().gens_values())
        logs, K_units = self._norm_logs(units)
        torsion = K_units.gens_orders()[0]
        width = len(logs[0])
        rows = [[ZZ(x) for x in row] for row in logs] + [[torsion if i == 0 else 0 for i in range(width)]]
        kernel = matrix(ZZ, rows).left_kernel().basis()
        result = []
        for c in kernel:
            element = self.extension(1)
            for exponent, u in zip(c[: len(units)], units):
                element *= u**exponent
            if element != 1:
                result.append(element)
        assert all(self.field(u.relative_norm()) == 1 for u in result)
        return result, units, logs, torsion

    def proper_generators(self):
        if not hasattr(self, "_proper"):
            norm_one, *_ = self.norm_one_units()
            moves = [self.iota(u) for u in norm_one] + [-identity_matrix(self.field, 2)]
            subgroup, _index = finite_preimage_generators(moves, lambda g: g, _integral)
            self._proper = tuple(subgroup) + (-identity_matrix(self.field, 2),)
            self._norm_one = norm_one
        return self._proper

    def _transversal(self):
        r"""Return representatives `\eta` of `O_L^1 / SO(M)` (left cosets, compared by integrality of `\eta_2^{-1} \eta_1`)."""
        self.proper_generators()
        one = self.extension(1)
        reps, position = [one], 0
        while position < len(reps):
            current = reps[position]
            position += 1
            for g in self._norm_one:
                moved = g * current
                if not any(_integral(self.iota(~known * moved)) for known in reps):
                    reps.append(moved)
        return reps

    def improper_element(self):
        r"""Return one integral improper automorph, or ``None`` as a proof that there is none (module docstring)."""
        K = self.field
        G = self.gram
        b = vector(K, [1, 0])
        sigma = identity_matrix(K, 2) - (2 / G[0, 0]) * b.column() * (b.row() * G)
        assert sigma.transpose() * G * sigma == G and sigma.det() == -1
        d = _lcm_denominator(sigma)
        entry = next(c for c in (self.delta * self.J).list() if c != 0)
        D = 2 * d * ZZ((~entry).denominator())
        L = self.extension
        _norm_one, units, logs, torsion = self.norm_one_units()
        K_units = K.unit_group()
        width = len(logs[0])
        transversal = self._transversal()
        target_norm = K.ideal(D**2)
        for ideal in L.ideal(D**2).divisors():
            if ideal.relative_norm() != target_norm or not ideal.is_principal():
                continue
            beta = ideal.gens_reduced()[0]
            ratio = K(D**2 / beta.relative_norm())
            target = vector(ZZ, [ZZ(x) for x in K_units.log(ratio)])
            equations = matrix(ZZ, [*([ZZ(x) for x in row] for row in logs), [torsion if i == 0 else 0 for i in range(width)]]).transpose()
            solution = _integer_solution(equations, target)
            if solution is None:
                continue
            epsilon = L(1)
            for exponent, u in zip(solution[: len(units)], units):
                epsilon *= u**exponent
            assert K((beta * epsilon).relative_norm()) == D**2
            for eta in transversal:
                t = beta * epsilon * eta / D
                g = sigma * self.iota(t)
                if _integral(g):
                    assert g.transpose() * G * g == G
                    return g
        return None

    def generators(self):
        improper = self.improper_element()
        return (*self.proper_generators(), *([improper] if improper is not None else []))
