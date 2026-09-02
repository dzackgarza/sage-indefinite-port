# TASK-01-2: Binary IO Protocol Documentation

## Overview

All `INDEF_FORM_*` executables share a common CLI pattern and IO format.
This document records the exact protocols discovered from source audit.

## Common CLI Pattern

```
EXECUTABLE [arith] [InputFile...] [OutFormat] [OutFile]
```

- `arith`: Selects arithmetic backend.
  Always first argument.

  - `gmp`: GMP native (recommended)

  - `gmp_boost`: Boost.GMP bindings

  - `multi_boost`: Boost.Multiprecision cpp_int

  - `safe`: SafeInt64 overflow-checked

- `InputFile...`: One or more matrix files in polyhedral_common format.

- `OutFormat`: Output serialization format.

  - `GAP`: GAP language (default for most executables)

  - `PYTHON`: Python list-of-lists

  - `CPP`: Native C++ format

- `OutFile`: Output destination.

  - `stdout`: Standard output

  - `stderr`: Standard error (default)

  - Other: Filename

## Matrix File Format

```
nrows ncols
 entry_0_0 entry_0_1 ... entry_0_(ncols-1)
 entry_1_0 entry_1_1 ... entry_1_(ncols-1)
 ...
 entry_(nrows-1)_0 ... entry_(nrows-1)_(ncols-1)
```

All entries are space-separated integers (for integer matrices) or rationals (for rational matrices, format: `numerator/denominator`).

## Vector File Format

Vectors are stored as single-column matrices:
```
n 1
 v_0
 v_1
 ...
 v_(n-1)
```

## Output Formats

### GAP Format

```gap
return [ [[1,0],[0,1]], [[0,1],[1,0]] ];
```
Or for records:
```gap
return rec(B:=[[1,0],[0,1]], Mred:=[[0,1],[1,0]]);
```
Or for single matrix:
```gap
return [[1,0],[0,1]];
```
Or for failure:
```gap
return fail;
```

### PYTHON Format

```python
[[1,0],[0,1]]
```
Or:
```python
None
```

## Executable-Specific Protocols

### INDEF_FORM_AutomorphismGroup

- **Input**: 1 matrix file (Gram matrix Q)

- **Output**: List of generator matrices (P_i with P_i^T Q P_i = Q)

- **Dispatch**: h=0 → definite leaf, h=1 → Lorentzian perfect domain, h>1 → approximate model

- **Dependencies**: libgap (for finite quotient actions), Normaliz/cddlib (for polyhedral cones)

### INDEF_FORM_TestEquivalence

- **Input**: 2 matrix files (Gram matrices Q1, Q2)

- **Output**: Transporter matrix P (P^T Q1 P = Q2) or fail

- **Dispatch**: Same as AutomorphismGroup

- **Dependencies**: Same as AutomorphismGroup

### INDEF_FORM_GetOrbit_IsotropicKplane

- **Input**: 1 matrix file + integer k + choice ("plane" or "flag")

- **Output**: List of matrices (rows span isotropic sublattice orbit reps)

- **Dependencies**: AutomorphismGroup (for stabilizer computation)

### INDEF_FORM_StabilizerIsotropicPlane

- **Input**: 1 matrix file + 1 plane file + choice

- **Output**: List of generator matrices of stabilizer

- **Dependencies**: AutomorphismGroup (recursive call on quotient)

### INDEF_FORM_InvariantIsotropicPlane

- **Input**: 1 matrix file + 1 plane file + choice

- **Output**: Single integer (hash)

- **Dependencies**: None (pure computation)

### INDEF_FORM_GetOrbitRepresentative

- **Input**: 1 matrix file + rational Xnorm

- **Output**: Matrix of orbit representative vectors

- **Dependencies**: AutomorphismGroup (for orbit traversal)

### INDEF_FORM_TestEquivalenceVector

- **Input**: 1 matrix file + 2 vector files

- **Output**: Transporter matrix or fail

- **Dependencies**: AutomorphismGroup (for orbit traversal)

### INDEF_FORM_ApproxCanonicalForm

- **Input**: 1 matrix file

- **Output**: Record with B (transformation) and Mred (reduced form)

- **Dependencies**: None (LLL-like reduction)

## Internal Data Flow

```
Main entry → ReadMatrixFile → IndefiniteCombinedAlgo
  ↓
INDEF_FORM_GetAttackScheme (computes Witt index h)
  ↓
  h=0 → INDEF_FORM_AutomorphismGroup_PosNeg (definite binary)
  h=1 → LORENTZ_GetGeneratorsAutom (perfect domain traversal)
  h>1 → INDEF_FORM_GetApproximateModel → stabilizer + covering orbit
  ↓
WriteListMatrix (output serialization)
```

## Temporary Files

- The executables read from the input file and write to the output file/stdout/stderr.

- No intermediate temporary files are created by the executables themselves.

- The `ReadMatrixFile` function reads directly from the path argument.
