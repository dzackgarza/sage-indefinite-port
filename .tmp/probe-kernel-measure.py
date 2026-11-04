"""Measured: span of kernel elements reachable by words, versus the full kernel the provider computes."""
import sys
sys.path.insert(0, "/work/tests")
from sage.all import ZZ, identity_matrix, matrix, vector
from tests.acceptance.consumer import element, lattice, rows_of
from sage_indefinite_port.indefinite.recursive import EichlerRecursion, _column_matrix
from test_degenerate_pointed_markings import _split_blocks, _stabilizer_word_shears, _blocks, _hom_element_matrix

RECURSION = EichlerRecursion()
for c_rank, r_rank, marked in ((2, 1, (1, 0, 0)), (3, 1, (1, 0, 0, 0)), (2, 2, (1, 0, 0, 0)), (3, 2, (1, 0, 0, 0, 0))):
    shears = _stabilizer_word_shears(c_rank, r_rank, marked)
    rows = [[int(e) for row in shear.rows() for e in row] for shear in shears]
    reachable = matrix(ZZ, rows).rank() if rows else 0
    complement, radical = _blocks(c_rank, r_rank)
    x = complement.linear_combination({complement.module_generating_set()[0]: 1})
    kernel = RECURSION.shear_kernel(complement, radical, x)
    full = r_rank * (c_rank - 1)
    print(f"C={c_rank} R={r_rank}: words<=3 give {len(shears)} shears, span rank {reachable} | "
          f"kernel rank {len(kernel)} (expected {full}) | {'GENERATED' if reachable == full else 'lower bound only'}",
          flush=True)
