"""The pinned indefinite-form change of integral presentation.

Port of polyhedral_common@1592b246 src_isotropy/Indefinite_LLL.h:34-945.
The preamble LatticeReduction retains an exact integral isometry.
"""

from __future__ import annotations

from random import choice, randrange

from dzack_research.preamble.categories.definite_lattices import (
    LatticeReduction,
    _reduction_from_backend_rows,
)
from dzack_research.preamble.tensors.tensor import _engine_component_matrix
from sage.matrix.constructor import identity_matrix, matrix
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational_field import QQ as SageQQ


def _gram_l1(gram):
    """The source get_l1_norm (lines 368-377)."""
    return sum(abs(entry) for row in gram.rows() for entry in row)


def _indefinite_gram_schmidt(gram, basis):
    """GramSchmidtOrthonormalization (lines 43-84)."""
    rank = gram.nrows()
    mu = matrix(SageQQ, rank, rank)
    qgram = gram.change_ring(SageQQ)
    orthogonal = []
    norms = []
    for i in range(rank):
        residual = vector(SageQQ, basis.row(i))
        for j in range(i):
            value = residual.dot_product(qgram * orthogonal[j]) / norms[j]
            mu[i, j] = value
            residual -= value * orthogonal[j]
        square = residual.dot_product(qgram * residual)
        if square == 0:
            return None, None, None, residual
        orthogonal.append(residual)
        norms.append(square)
    return orthogonal, norms, mu, None


def _indefinite_lll(gram):
    """The signed LLL loop of Denis Simon (lines 94-188)."""
    rank = gram.nrows()
    basis = identity_matrix(SageZZ, rank)
    if rank <= 1:
        if rank == 0:
            isotropic = vector(SageQQ, [])
        elif gram[0, 0] == 0:
            isotropic = vector(SageQQ, [1])
        else:
            isotropic = None
        return basis, gram, isotropic
    k = 1
    while True:
        _, norms, mu, isotropic = _indefinite_gram_schmidt(gram, basis)
        if isotropic is not None:
            return basis, basis * gram * basis.transpose(), isotropic
        for i in range(rank - 1, -1, -1):
            for j in range(i):
                nearest = SageZZ(mu[i, j].round())
                if nearest:
                    basis.add_multiple_of_row(i, j, -nearest)
        sum1 = abs(norms[k] + mu[k, k - 1] ** 2 * norms[k - 1])
        sum2 = SageQQ(7) / 8 * abs(norms[k - 1])
        if sum1 < sum2:
            basis.swap_rows(k, k - 1)
            k = max(k - 1, 1)
        else:
            k += 1
        if k >= rank:
            return basis, basis * gram * basis.transpose(), None


def _random_signed_permutation(rank):
    """get_random_int_matrix, Indefinite_LLL.h:190-213."""
    permutation = list(range(rank))
    for _ in range(4 * rank):
        i, j = randrange(rank), randrange(rank)
        if i != j:
            permutation[i], permutation[j] = permutation[j], permutation[i]
    result = matrix(SageZZ, rank, rank)
    for i, j in enumerate(permutation):
        result[i, j] = choice((-1, 1))
    return result


def _compute_indefinite_reduction(gram, *, look_for_isotropic=False):
    """ComputeReductionIndefinite, Indefinite_LLL.h:215-319.

    Both source modes are retained. In isotropic-search mode a rational
    isotropic vector is returned in the input Gram coordinates. In
    reduction mode the random signed-permutation restarts terminate
    after 2*rank consecutive nonimprovements, exactly as upstream.
    """
    rank = gram.nrows()
    basis = identity_matrix(SageZZ, rank)
    if not rank:
        return basis, gram, vector(SageQQ, []) if look_for_isotropic else None
    working = gram
    best_norm = _gram_l1(working)
    no_improvement = 0
    while True:
        step, reduced, isotropic = _indefinite_lll(working)
        if isotropic is not None and look_for_isotropic:
            original = basis.transpose() * isotropic
            if original.dot_product(gram.change_ring(SageQQ) * original):
                raise ArithmeticError("indefinite LLL returned a nonisotropic vector")
            return step * basis, reduced, original
        basis = step * basis
        working = reduced
        current_norm = _gram_l1(working)
        if current_norm >= best_norm:
            no_improvement += 1
            if no_improvement == 2 * rank:
                return basis, working, None
        else:
            no_improvement = 0
            best_norm = current_norm
        perturbation = _random_signed_permutation(rank)
        basis = perturbation * basis
        working = perturbation * working * perturbation.transpose()


def _simple_indefinite_reduction(gram):
    """SimpleIndefiniteReduction, Indefinite_LLL.h:379-588.

    Every accepted shear strictly decreases the nonnegative integer
    L1 norm. This is the source termination condition, not a cutoff.
    """
    rank = gram.nrows()
    basis = identity_matrix(SageZZ, rank)
    if rank < 2:
        return basis, gram
    working = gram
    indices = list(range(rank))

    def delta(i, j, c):
        delta_off = sum(abs(working[i, k]) - abs(working[i, k] + c * working[j, k]) for k in range(rank) if k != i)
        delta_diag = abs(working[i, i]) - abs(working[i, i] + 2 * c * working[j, i] + c * c * working[j, j])
        return 2 * delta_off + delta_diag

    while True:
        first, second = randrange(rank), randrange(rank)
        if first != second:
            indices[first], indices[second] = indices[second], indices[first]
        found = None
        for i in indices:
            for j in indices:
                if i == j:
                    continue
                for sign in (-1, 1):
                    improvement = delta(i, j, sign)
                    if improvement <= 0:
                        continue
                    coefficient = sign
                    while True:
                        later = delta(i, j, coefficient + sign)
                        if later <= improvement:
                            break
                        coefficient += sign
                        improvement = later
                    found = i, j, coefficient
                    break
                if found is not None:
                    break
            if found is not None:
                break
        if found is None:
            return basis, working
        i, j, coefficient = found
        elementary = identity_matrix(SageZZ, rank)
        elementary[i, j] = coefficient
        basis = elementary * basis
        working = elementary * working * elementary.transpose()


def _gram_connected_components(gram):
    """GetMatrixGraph/MatrixConnectedComponents, lines 340-366."""
    from sage.graphs.graph import Graph

    rank = gram.nrows()
    graph = Graph(rank)
    graph.add_edges((i, j) for i in range(rank) for j in range(i + 1, rank) if gram[i, j] != 0)
    return tuple(tuple(sorted(component)) for component in graph.connected_components())


def _definite_lll(gram, *, dual):
    """ClassicLLL.h:40-284 via Sage's exact Gram LLL.

    Dual LLL reduces the integral adjugate and applies P^(-T) to the
    original Gram, just as in ClassicLLL.h:246-284.
    """
    positive = gram if gram.is_positive_definite() else -gram
    if dual:
        transform = positive.adjugate().LLL_gram().inverse().transpose()
    else:
        transform = positive.LLL_gram()
    transform = transform.change_ring(SageZZ)
    return transform, transform * gram * transform.transpose()


def _nondegenerate_reduction(gram):
    """IndefiniteReductionNonDegenerate, Indefinite_LLL.h:615-822."""
    rank = gram.nrows()
    basis = identity_matrix(SageZZ, rank)
    if rank == 0:
        return basis, gram
    if gram.rank() != rank:
        raise ValueError("expected the nondegenerate block of a Gram form")
    # The source signature branches only choose whether definite LLL
    # and its dual are admissible.
    definite = gram.is_positive_definite() or (-gram).is_positive_definite()
    choices = (0, 1, 2, 3) if definite else (0, 1)
    working = gram
    best_norm = _gram_l1(working)
    while True:
        accepted = 0
        for which in choices:
            blocks = _gram_connected_components(working)
            if len(blocks) > 1:
                block_basis = matrix(SageZZ, rank, rank)
                block_gram = matrix(SageZZ, rank, rank)
                for block in blocks:
                    component = working.matrix_from_rows_and_columns(block, block)
                    step, reduced = _nondegenerate_reduction(component)
                    for i, original_i in enumerate(block):
                        for j, original_j in enumerate(block):
                            block_basis[original_i, original_j] = step[i, j]
                            block_gram[original_i, original_j] = reduced[i, j]
                return block_basis * basis, block_gram
            if which == 0:
                step, candidate, _ = _compute_indefinite_reduction(working)
            elif which == 1:
                step, candidate = _simple_indefinite_reduction(working)
            else:
                step, candidate = _definite_lll(working, dual=which == 3)
            current_norm = _gram_l1(candidate)
            if current_norm < best_norm:
                accepted += 1
                basis = step * basis
                working = candidate
                best_norm = current_norm
        if accepted < 2:
            return basis, working


def _indefinite_reduction(gram):
    """IndefiniteReduction, Indefinite_LLL.h:823-930, radical included.

    Integral row HNF with its unimodular transformation computes a
    primitive integral nullspace and complementary basis, equivalent to
    the source NullspaceIntMat/SubspaceCompletionInt construction.
    """
    if gram.nrows() != gram.ncols() or gram != gram.transpose():
        raise ValueError("a Gram matrix must be square and symmetric")
    rank = gram.nrows()
    if not rank:
        return identity_matrix(SageZZ, 0), gram
    hnf, unimodular = gram.hermite_form(transformation=True, include_zero_rows=True)
    nondegenerate_rank = gram.rank()
    if unimodular * gram != hnf or abs(unimodular.det()) != 1:
        raise ArithmeticError("integral nullspace splitting is not unimodular")
    separated = unimodular * gram * unimodular.transpose()
    if any(separated[i, j] for i in range(nondegenerate_rank, rank) for j in range(rank)):
        raise ArithmeticError("the HNF radical basis fails to annihilate the form")
    nondegenerate = separated.matrix_from_rows_and_columns(range(nondegenerate_rank), range(nondegenerate_rank))
    step, _ = _nondegenerate_reduction(nondegenerate)
    combined = identity_matrix(SageZZ, rank)
    for i in range(nondegenerate_rank):
        for j in range(nondegenerate_rank):
            combined[i, j] = step[i, j]
    basis = combined * unimodular
    reduced = basis * gram * basis.transpose()
    if abs(basis.det()) != 1:
        raise ArithmeticError("indefinite reduction must retain an integral isometry")
    return basis, reduced


def indefinite_reduced_presentation(lattice) -> LatticeReduction:
    """Isometric pre-reduction for the pinned group/isometry recursion.

    IndefiniteForms.g:1379-1388,1480-1492; the change of presentation
    is raised through the existing preamble LatticeReduction constructor.
    """
    original = _engine_component_matrix(lattice.gram_tensor()).change_ring(SageZZ)
    basis, reduced = _indefinite_reduction(original)
    if basis * original * basis.transpose() != reduced:
        raise ArithmeticError("reduced presentation has an incorrect Gram witness")
    return _reduction_from_backend_rows(lattice, basis)
