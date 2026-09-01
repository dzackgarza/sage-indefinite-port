"""
Differential test harness for polyhedral_common upstream binaries.

Runs each wrapper function against a corpus of test lattices,
captures exact outputs, and stores them as golden fixtures.
Later: re-run and compare against SageMath port outputs.
"""

import json
import os
import subprocess
import sys
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import IO, Any

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def find_binary(name: str) -> str | None:
    """Locate a polyhedral_common binary on PATH or in common locations."""
    # Check PATH
    result = subprocess.run(["which", name], capture_output=True, text=True)
    if result.returncode == 0:
        return result.stdout.strip()

    # Common install locations
    common = [
        Path.home() / "polyhedral_common" / "build" / "src_indefinite" / name,
        Path("/usr/local/bin") / name,
        Path("/usr/bin") / name,
    ]
    for p in common:
        if p.exists():
            return str(p)

    return None


def write_matrix_to_file(f: IO[str], matrix: Sequence[Sequence[int]]) -> None:
    """Write a matrix in polyhedral_common format: first line nrows ncols, then rows."""
    nrows = len(matrix)
    ncols = len(matrix[0]) if nrows > 0 else 0
    f.write(f"{nrows} {ncols}\n")
    for row in matrix:
        entries = " ".join(str(x) for x in row)
        f.write(f" {entries}\n")


def write_scalar_to_file(f: IO[str], val: int | str) -> None:
    """Write a scalar value."""
    f.write(f"{val}\n")


def parse_matrix_from_file(filepath: str) -> list[list[int]]:
    """Read a matrix from polyhedral_common output format."""
    with open(filepath) as f:
        lines = f.read().strip().split("\n")

    if not lines or not lines[0].strip():
        return []

    first = lines[0].strip().split()
    if len(first) < 2:
        return []

    nrows, _ncols = int(first[0]), int(first[1])
    matrix: list[list[int]] = []
    for i in range(1, nrows + 1):
        if i >= len(lines):
            break
        row = [int(x) for x in lines[i].strip().split()]
        matrix.append(row)
    return matrix


def parse_list_matrix_from_file(filepath: str) -> list[list[list[int]]]:
    """Read a list of matrices from polyhedral_common output format."""
    with open(filepath) as f:
        lines = f.read().strip().split("\n")

    if not lines or not lines[0].strip():
        return []

    n_mat = int(lines[0].strip())
    matrices: list[list[list[int]]] = []
    idx = 1
    for _ in range(n_mat):
        if idx >= len(lines):
            break
        first = lines[idx].strip().split()
        nrows, _ncols = int(first[0]), int(first[1])
        idx += 1
        matrix: list[list[int]] = []
        for _ in range(nrows):
            if idx >= len(lines):
                break
            row = [int(x) for x in lines[idx].strip().split()]
            matrix.append(row)
            idx += 1
        matrices.append(matrix)
    return matrices


class BinaryRunner:
    """Run a polyhedral_common binary and capture output."""

    def __init__(self, binary_name: str) -> None:
        self.binary_name = binary_name
        self.binary_path = find_binary(binary_name)
        self.available = self.binary_path is not None

    def run_with_matrices(
        self,
        input_matrices: Sequence[Sequence[Sequence[int]]],
        extra_args: Sequence[str] | None = None,
    ) -> dict[str, Any] | None:
        """
        Write matrices to temp files, run binary, read output.
        Returns (output_matrices, stderr_text) or None if binary not found.
        """
        if not self.available or self.binary_path is None:
            return None

        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as fin:
            for mat in input_matrices:
                write_matrix_to_file(fin, mat)
            input_path = fin.name

        output_path = tempfile.mktemp(suffix=".txt")

        try:
            cmd = [self.binary_path, input_path, output_path]
            if extra_args:
                cmd.extend(extra_args)

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

            if result.returncode != 0:
                return {"error": result.stderr, "returncode": result.returncode}

            matrices = parse_list_matrix_from_file(output_path)
            return {"matrices": matrices, "stderr": result.stderr}

        finally:
            for p in [input_path, output_path]:
                if os.path.exists(p):
                    os.unlink(p)


def test_testequivalence() -> dict[str, Any] | None:
    """Test INDEF_FORM_TestEquivalence on known pairs."""
    runner = BinaryRunner("INDEF_FORM_TestEquivalence")
    if not runner.available:
        print("SKIP: INDEF_FORM_TestEquivalence binary not found")
        return None

    results: dict[str, Any] = {}

    # Test: identical forms should be equivalent
    gram = [[0, 1], [1, 0]]
    result = runner.run_with_matrices([gram, gram])
    results["self_equivalence_U"] = result

    # Test: non-equivalent forms
    gram_a = [[-1, 1, 0], [1, 63, 0], [0, 0, 2]]
    gram_b = [[-9, 1, 0], [1, 7, 0], [0, 0, 2]]
    result = runner.run_with_matrices([gram_a, gram_b])
    results["SPLAG_51a_vs_51b"] = result

    return results


def test_automorphism_group() -> dict[str, Any] | None:
    """Test INDEF_FORM_AutomorphismGroup on known lattices."""
    runner = BinaryRunner("INDEF_FORM_AutomorphismGroup")
    if not runner.available:
        print("SKIP: INDEF_FORM_AutomorphismGroup binary not found")
        return None

    results: dict[str, Any] = {}

    # Test: E8 automorphism group (known order 696729600)
    e8_gram = [
        [2, -1, 0, 0, 0, 0, 0, 0],
        [-1, 2, -1, 0, 0, 0, 0, 0],
        [0, -1, 2, -1, 0, 0, 0, -1],
        [0, 0, -1, 2, -1, 0, 0, 0],
        [0, 0, 0, -1, 2, -1, 0, 0],
        [0, 0, 0, 0, -1, 2, -1, 0],
        [0, 0, 0, 0, 0, -1, 2, 0],
        [0, 0, -1, 0, 0, 0, 0, 2],
    ]
    result = runner.run_with_matrices([e8_gram])
    results["E8_automorphisms"] = result

    # Test: hyperbolic plane U (infinite group, should get generators)
    u_gram = [[0, 1], [1, 0]]
    result = runner.run_with_matrices([u_gram])
    results["U_automorphisms"] = result

    return results


def test_orbit_representative() -> dict[str, Any] | None:
    """Test INDEF_FORM_GetOrbitRepresentative on known isotropic vectors."""
    runner = BinaryRunner("INDEF_FORM_GetOrbitRepresentative")
    if not runner.available:
        print("SKIP: INDEF_FORM_GetOrbitRepresentative binary not found")
        return None

    results: dict[str, Any] = {}

    # Test: isotropic vectors of U (should have 1 orbit rep)
    u_gram = [[0, 1], [1, 0]]
    result = runner.run_with_matrices([u_gram], extra_args=["0"])
    results["U_isotropic"] = result

    # Test: diag(1,-1) isotropic vectors
    diag_gram = [[1, 0], [0, -1]]
    result = runner.run_with_matrices([diag_gram], extra_args=["0"])
    results["diag_1_neg1_isotropic"] = result

    return results


def run_all_tests() -> dict[str, Any]:
    """Run all differential tests and collect results."""
    all_results: dict[str, Any] = {}

    print("=== Differential Harness: polyhedral_common binary tests ===\n")

    print("1. INDEF_FORM_TestEquivalence")
    all_results["TestEquivalence"] = test_testequivalence()

    print("2. INDEF_FORM_AutomorphismGroup")
    all_results["AutomorphismGroup"] = test_automorphism_group()

    print("3. INDEF_FORM_GetOrbitRepresentative")
    all_results["OrbitRepresentative"] = test_orbit_representative()

    # Save golden fixtures
    fixtures_path = PROJECT_ROOT / "tests" / "fixtures" / "golden_fixtures.json"
    with open(fixtures_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    print(f"\nGolden fixtures saved to {fixtures_path}")
    return all_results


if __name__ == "__main__":
    run_all_tests()
