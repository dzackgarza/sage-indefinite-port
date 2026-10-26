"""The pinned indefinite-form change of integral presentation.

Port of polyhedral_common@1592b246 src_isotropy/Indefinite_LLL.h:34-945.
The preamble LatticeReduction retains an exact integral isometry.
"""

from __future__ import annotations

from functools import cache
from random import choice, randrange

from dzack_research.preamble.categories.definite_lattices import (
    LatticeReduction,
    _reduction_from_backend_rows,
)
from dzack_research.preamble.tensors.tensor import _engine_component_matrix
from sage.matrix.constructor import identity_matrix, matrix
from sage.libs.gap.libgap import libgap
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational_field import QQ as SageQQ


def _gram_l1(gram):
    """The source get_l1_norm (lines 368-377)."""
    return sum(abs(entry) for row in gram.rows() for entry in row)


def _source_nearest_integer(value):
    r"""Port GAP ``NearestInteger`` with exact half-ties toward +infinity."""
    return SageZZ((SageQQ(value) + SageQQ(1) / 2).floor())


def _indefinite_gram_schmidt(gram, basis):
    """GramSchmidtOrthonormalization (lines 43-84), on the Gram matrix of the basis.

    With P = B G B^T, the source's coefficients and squared norms are
    mu_ij = (P_ij - sum_{k<j} mu_jk mu_ik N_k) / N_j and
    N_i = P_ii - sum_{k<i} mu_ik^2 N_k. When N_i = 0 the source returns the
    isotropic residual b_i - sum_j mu_ij b*_j, whose coordinates in the
    basis are tracked in ``coefficients`` (b*_j = sum_k C_jk b_k).
    """
    rank = gram.nrows()
    pairings = (basis * gram * basis.transpose()).change_ring(SageQQ)
    mu = matrix(SageQQ, rank, rank)
    norms = []
    coefficients = []
    for i in range(rank):
        row = [SageQQ.zero()] * rank
        row[i] = SageQQ.one()
        for j in range(i):
            value = (pairings[i, j] - sum(mu[j, k] * mu[i, k] * norms[k] for k in range(j))) / norms[j]
            mu[i, j] = value
            for k in range(j + 1):
                row[k] -= value * coefficients[j][k]
        square = pairings[i, i] - sum(mu[i, k] ** 2 * norms[k] for k in range(i))
        if square == 0:
            return None, None, vector(SageQQ, row) * basis.change_ring(SageQQ)
        coefficients.append(row)
        norms.append(square)
    return norms, mu, None


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
        norms, mu, isotropic = _indefinite_gram_schmidt(gram, basis)
        if isotropic is not None:
            return basis, basis * gram * basis.transpose(), isotropic
        for i in range(rank - 1, -1, -1):
            for j in range(i):
                nearest = _source_nearest_integer(mu[i, j])
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


def _canonicalize_permutation_signs(gram):
    """The signed graph canonicalization in IndefApproxCanonical.h:22-119.

    Canonicalize the absolute weighted adjacency matrix using Sage's labelled
    graph backend. Diagonal weights are vertex colours; off-diagonal weights
    are edge labels. The reordering and subsequent sign propagation are
    integral unimodular transformations. This is only *approximately*
    canonical: differing presentations never certify nonisometry.
    """
    from sage.graphs.graph import Graph

    rank = gram.nrows()
    if rank <= 1:
        return identity_matrix(SageZZ, rank), gram

    graph = Graph()
    graph.add_vertices(range(rank))
    for i in range(rank):
        for j in range(i + 1, rank):
            if gram[i, j]:
                graph.add_edge(i, j, abs(gram[i, j]))
    colours = sorted({abs(gram[i, i]) for i in range(rank)})
    partition = [[i for i in range(rank) if abs(gram[i, i]) == colour] for colour in colours]
    _, canonical = graph.canonical_label(partition=partition, certificate=True, edge_labels=True, algorithm="sage")
    basis = matrix(SageZZ, rank, rank)
    for new, old in enumerate(sorted(range(rank), key=canonical.__getitem__)):
        basis[new, old] = 1

    reordered = basis * gram * basis.transpose()
    signs = identity_matrix(SageZZ, rank)
    for component in _gram_connected_components(reordered):
        reached = {component[0]}
        pending = [component[0]]
        while pending:
            vertex = pending.pop(0)
            for neighbour in component:
                if neighbour not in reached and reordered[vertex, neighbour]:
                    signs[neighbour, neighbour] = signs[vertex, vertex] if reordered[vertex, neighbour] > 0 else -signs[vertex, vertex]
                    reached.add(neighbour)
                    pending.append(neighbour)
    basis = signs * basis
    return basis, basis * gram * basis.transpose()


def _order_blocks_by_signature(gram):
    """Stable block ordering of IndefApproxCanonical.h:125-276.

    The source sorts by two-dimensional Eichler status, dimension,
    coefficient L1 norm and maximal absolute coefficient, in that order.
    """
    blocks = _gram_connected_components(gram)

    def block_order_key(component):
        block = gram.matrix_from_rows_and_columns(component, component)
        eichler = len(component) == 2 and all(block[i, i] == 0 for i in range(2))
        return (
            not eichler,
            len(component),
            _gram_l1(block),
            max((abs(x) for x in block.list()), default=0),
        )

    basis = matrix(SageZZ, gram.nrows(), gram.ncols())
    position = 0
    for component in sorted(blocks, key=block_order_key):
        for old in component:
            basis[position, old] = 1
            position += 1
    return basis, basis * gram * basis.transpose()


def _approximate_indefinite_reduction(gram):
    """The present exact-presentation part of ApproxCanonicalIndefiniteForm.

    The source IndefiniteReduction is followed by blockwise indefinite
    sign/permutation normalization and stable ordering (lines 374-389).
    For definite blocks of rank greater than one, its separate definite
    ComputeCanonicalForm backend (LatticeStabEquiCan.h) is not yet supplied.
    Those blocks retain their *isometric* reduced presentation; no differing
    approximate presentation is used to decide nonisometry.
    """
    first, reduced = _indefinite_reduction(gram)
    second = identity_matrix(SageZZ, gram.nrows())
    for component in _gram_connected_components(reduced):
        block = reduced.matrix_from_rows_and_columns(component, component)
        if not (block.is_positive_definite() or (-block).is_positive_definite()):
            transform, _ = _canonicalize_permutation_signs(block)
            for i, old_i in enumerate(component):
                for j, old_j in enumerate(component):
                    second[old_i, old_j] = transform[i, j]
    reordered = second * reduced * second.transpose()
    third, _ = _order_blocks_by_signature(reordered)
    basis = third * second * first
    return basis, basis * gram * basis.transpose()


def _definite_lll(gram, *, dual):
    """ClassicLLL.h:40-284 via Sage's exact Gram LLL.

    Sage's LLL_gram returns column coordinates U with U^T G U reduced,
    whereas ClassicLLL uses row coordinates P with P G P^T reduced.
    Dual LLL reduces the integral adjugate and applies P^(-T) to the
    original Gram, as in ClassicLLL.h:246-284.
    """
    positive = gram if gram.is_positive_definite() else -gram
    if dual:
        # P = U^T on the dual form, hence P^(-T) = U^(-1).
        transform = positive.adjugate().LLL_gram().inverse()
    else:
        transform = positive.LLL_gram().transpose()
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


def _source_sublattice_basis_reduction(basis):
    r"""Port ``SublatticeBasisReductionKernel`` from ClassicLLL.h:410-448."""
    reduced_basis = matrix(SageZZ, basis)
    norm = _gram_l1(reduced_basis)
    while reduced_basis.nrows() > 1:
        changed = False
        for dual in (False, True):
            row_gram = reduced_basis * reduced_basis.transpose()
            gram = row_gram.adjugate() if dual else row_gram
            reduction = libgap.LLLReducedGramMat(
                [[int(gram[i, j]) for j in range(gram.ncols())] for i in range(gram.nrows())],
                libgap(99) / libgap(100),
            )
            transform = matrix(SageZZ, [[SageZZ(int(entry)) for entry in row] for row in reduction["transformation"]])
            reduced_gram = matrix(SageZZ, [[SageZZ(int(entry)) for entry in row] for row in reduction["remainder"]])
            if abs(transform.det()) != 1:
                raise ArithmeticError("sublattice LLL returned a nonunimodular transform")
            if dual:
                # ClassicLLL's dual transform is P^{-T} for P reducing adj(G).
                transform = transform.inverse().transpose()
                expected = row_gram.det() * reduced_gram.inverse()
            else:
                expected = reduced_gram
            candidate = transform * reduced_basis
            if candidate * candidate.transpose() != expected:
                raise ArithmeticError("sublattice LLL transform and Gram remainder disagree")
            candidate_norm = _gram_l1(candidate)
            if candidate_norm < norm:
                reduced_basis, norm, changed = candidate, candidate_norm, True
        if not changed:
            break
    return reduced_basis


def _source_integer_nullspace(gram):
    r"""Port ``NullspaceIntMat`` followed by ``SublatticeBasisReduction``."""
    rows = [[int(entry) for entry in row] for row in gram.rows()]
    raw = libgap.NullspaceIntMat(rows)
    kernel = matrix(SageZZ, [[SageZZ(int(entry)) for entry in row] for row in raw]) if len(raw) else matrix(SageZZ, 0, gram.ncols())
    if kernel * gram:
        raise ArithmeticError("GAP integer nullspace does not annihilate the Gram matrix")
    if kernel.nrows() != gram.ncols() - gram.rank():
        raise ArithmeticError("GAP integer nullspace has the wrong rank")
    return _source_sublattice_basis_reduction(kernel)


def _source_subspace_completion_int(kernel_basis, dimension):
    r"""Port ``SubspaceCompletionInt`` from IndefiniteForms.g:471-497.

    Smith transforms U and V satisfy U*K*V=D. For a primitive rank-d row
    sublattice, D has d unit pivots; embedding U in the leading block of an
    n-by-n identity and applying U^{-1}V^{-1} yields the pinned completion.
    """
    if dimension < 0 or kernel_basis.ncols() != dimension:
        raise ValueError("kernel basis and ambient dimension disagree")
    d = int(kernel_basis.nrows())
    if d > dimension or kernel_basis.rank() != d:
        raise ValueError("subspace completion requires independent rows")
    if d == 0:
        return identity_matrix(SageZZ, dimension)
    diagonal, row_transform, column_transform = kernel_basis.smith_form(transformation=True)
    if row_transform * kernel_basis * column_transform != diagonal:
        raise ArithmeticError("Smith transformations do not certify the kernel basis")
    if any(diagonal[i, i] != 1 for i in range(d)) or any(
        diagonal[i, j] for i in range(d) for j in range(dimension) if i != j
    ):
        raise ValueError("the row sublattice is not primitive in the ambient integer lattice")
    block = identity_matrix(SageZZ, dimension)
    for i in range(d):
        for j in range(d):
            block[i, j] = row_transform[i, j]
    full_basis = block.inverse() * column_transform.inverse()
    completion = full_basis.matrix_from_rows(range(d, dimension))
    combined = completion.stack(kernel_basis)
    if abs(combined.det()) != 1:
        raise ArithmeticError("subspace completion did not produce a unimodular basis")
    return completion


def _indefinite_reduction(gram):
    """IndefiniteReduction, Indefinite_LLL.h:823-930, radical included.

    The source NullspaceIntMat basis is reduced by alternating direct/dual
    sublattice LLL, then SubspaceCompletionInt supplies the unimodular
    complement before the nondegenerate block is reduced.
    """
    if gram.nrows() != gram.ncols() or gram != gram.transpose():
        raise ValueError("a Gram matrix must be square and symmetric")
    rank = gram.nrows()
    if not rank:
        return identity_matrix(SageZZ, 0), gram
    nondegenerate_rank = int(gram.rank())
    # Indefinite_LLL.h:841-857: NullspaceIntMat, LLL reduction of its row
    # lattice, and Smith-certified SubspaceCompletionInt in source order.
    radical = _source_integer_nullspace(gram)
    complement = _source_subspace_completion_int(radical, rank)
    full_basis = complement.stack(radical)
    separated = full_basis * gram * full_basis.transpose()
    if any(separated[i, j] for i in range(nondegenerate_rank, rank) for j in range(rank)):
        raise ArithmeticError("the completed integral radical basis fails to annihilate the form")
    nondegenerate = separated.matrix_from_rows_and_columns(range(nondegenerate_rank), range(nondegenerate_rank))
    step, _ = _nondegenerate_reduction(nondegenerate)
    combined = identity_matrix(SageZZ, rank)
    for i in range(nondegenerate_rank):
        for j in range(nondegenerate_rank):
            combined[i, j] = step[i, j]
    basis = combined * full_basis
    reduced = basis * gram * basis.transpose()
    if abs(basis.det()) != 1:
        raise ArithmeticError("indefinite reduction must retain an integral isometry")
    return basis, reduced


@cache
def indefinite_reduced_presentation(lattice) -> LatticeReduction:
    """Isometric pre-reduction for the pinned group/isometry recursion.

    IndefiniteForms.g:1379-1388,1480-1492; the change of presentation
    is raised through the existing preamble LatticeReduction constructor.
    The randomized reduction is memoized per lattice, so every consumer of a
    lattice works in one reduced presentation and shares its perfect complex.
    """
    original = _engine_component_matrix(lattice.gram_tensor()).change_ring(SageZZ)
    # IndefiniteForms.g:1379-1383 and :1480-1485 call IndefiniteReduction,
    # not the separate ApproxCanonicalIndefiniteForm (C++ IndefApproxCanonical.h).
    # Canonicalizing the blocks here changes the source recursion's reduction.
    basis, reduced = _indefinite_reduction(original)
    if basis * original * basis.transpose() != reduced:
        raise ArithmeticError("reduced presentation has an incorrect Gram witness")
    return _reduction_from_backend_rows(lattice, basis)
