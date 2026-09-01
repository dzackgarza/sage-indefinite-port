"""
Binary IO protocol documentation for INDEF_FORM_* executables.

Derived from source analysis of polyhedral_common/src_indefinite/ and
Indefinite.jl/src/ExternalCalls.jl. These are the serialization contracts
that any SageMath replacement must reproduce exactly.
"""

PROTOCOLS = {
    "INDEF_FORM_TestEquivalence": {
        "input": {
            "format": "two matrix files on command line",
            "file1": "Qmat1: nrows ncols\\n <space> row_entries...\\n",
            "file2": "Qmat2: nrows ncols\\n <space> row_entries...\\n",
            "cli_args": ["gmp", input_file_1, input_file_2],
        },
        "output": {
            "format": "single matrix file or empty",
            "on_equivalent": "nrows ncols\\n <space> P_entries...\\n  (P * Q1 * P^T = Q2)",
            "on_nonequivalent": "empty file",
        },
        "temp_files": ["input1.txt", "input2.txt", "output.txt"],
        "math_precondition": "Qmat1, Qmat2 symmetric integral indefinite",
    },
    "INDEF_FORM_AutomorphismGroup": {
        "input": {
            "format": "single matrix file",
            "file": "Qmat: nrows ncols\\n <space> row_entries...\\n",
            "cli_args": ["gmp", input_file],
        },
        "output": {
            "format": "list of matrices (AST literal via PYTHON backend)",
            "structure": "[[nrows, ncols, [entries...]], ...]",
            "note": "generator matrices in row-action convention",
        },
        "temp_files": ["input.txt", "output.txt"],
        "math_precondition": "Qmat symmetric integral indefinite",
    },
    "INDEF_FORM_GetOrbitRepresentative": {
        "input": {
            "format": "matrix file + scalar on command line",
            "file": "Qmat: nrows ncols\\n <space> row_entries...\\n",
            "cli_args": ["gmp", input_file, str(eNorm)],
        },
        "output": {
            "format": "list of row vectors (AST literal)",
            "structure": "[[n, [entries...]], ...]",
            "note": "each vector is 1 x n matrix",
        },
        "temp_files": ["input.txt", "output.txt"],
        "math_precondition": "Qmat indefinite, eNorm rational",
    },
    "INDEF_FORM_GetOrbit_IsotropicKplane": {
        "input": {
            "format": "matrix file + k + nature on command line",
            "file": "Qmat: nrows ncols\\n <space> row_entries...\\n",
            "cli_args": ["gmp", input_file, str(k), nature],
            "nature": "'plane' or 'flag'",
        },
        "output": {
            "format": "list of basis matrices (AST literal)",
            "structure": "[[k, n, [entries...]], ...]",
        },
        "temp_files": ["input.txt", "output.txt"],
        "math_precondition": "k <= min(p,q) where (p,q) is signature",
    },
    "INDEF_FORM_TestEquivalenceVector": {
        "input": {
            "format": "matrix file + two vector files",
            "cli_args": ["gmp", Q_file, v1_file, v2_file],
            "vector_format": "n\\n <space> entries...\\n",
        },
        "output": {
            "format": "matrix or None (AST literal)",
            "structure": "None if inequivalent, else [[n,n,[entries]],]",
        },
        "temp_files": ["Q.txt", "v1.txt", "v2.txt", "output.txt"],
    },
    "INDEF_FORM_StabilizerVector": {
        "input": {
            "format": "matrix file + vector file",
            "cli_args": ["gmp", Q_file, v_file],
        },
        "output": {
            "format": "list of generator matrices (AST literal)",
            "note": "generators satisfy M * Q * M^T = Q",
        },
        "temp_files": ["Q.txt", "v.txt", "output.txt"],
    },
    "INDEF_FORM_StabilizerIsotropicPlane": {
        "input": {
            "format": "matrix file + basis file + choice",
            "cli_args": ["gmp", Q_file, P_file, choice],
            "choice": "'plane' or 'flag'",
        },
        "output": {
            "format": "list of generator matrices (AST literal)",
        },
        "temp_files": ["Q.txt", "P.txt", "output.txt"],
    },
}

# Common serialization helpers (from ExternalCalls.jl and py_polyhedral/binaries.py)
MATRIX_FORMAT = """\
{nrows} {ncols}
 {e11} {e12} ... {e1n}
 {e21} {e22} ... {e2n}
 ...
 {en1} {en2} ... {enn}"""

VECTOR_FORMAT = """\
{n}
 {v1} {v2} ... {vn}"""
