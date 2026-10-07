#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Extract tests/fixtures/binary_form_automorphisms.json from Hecke's QuadBin tests.

Source: references/vendor/Hecke.jl@e2ab5716/test/QuadForm/QuadBin.jl, testset
"Automormorphism group". Every form a x^2 + b x y + c y^2 in the testset is indefinite;
its lattice has Gram matrix [[2a, b], [b, 2c]]. The testset records three kinds of
fact, each kept in its own list:

- explicit_generators: four forms with Hecke's exact automorphism-group generators;
- generator_counts: six forms where the test records only how many generators
  Hecke's automorphism_group_generators returns. That count is Hecke's output, not an
  invariant of the group, and is recorded under that name;
- improper_automorphisms: forms the test asserts are ambiguous, i.e. have an
  automorphism of determinant -1 (O(L) != SO(L)).

The script asserts that every form is indefinite (b^2 - 4ac > 0), that the testset
contains exactly these 4 + 6 + 1 cases, and determines the action convention from the
data: every recorded generator T must satisfy T G T^T = G for all cases, or T^T G T = G
for all cases; the convention that holds is recorded.

Run from the repository root: uv run references/extract/hecke_binary_forms.py
"""

import json
import re
from pathlib import Path

SOURCE = Path("references/vendor/Hecke.jl@e2ab5716/test/QuadForm/QuadBin.jl")
TARGET = Path("tests/fixtures/binary_form_automorphisms.json")
FORM = re.compile(r"^\s*([fg]) = binary_quadratic_form\((-?\d+), (-?\d+), (-?\d+)\)$")
GENS = re.compile(r"^\s*@test gens == \[(.*)\]$")
COUNT = re.compile(r"^\s*@test length\(gens\) == (\d+)$")
IMPROPER = re.compile(r"^\s*@assert any\(T -> det\(T\) == -1, gens\) # g is ambiguous$")
MATRIX = re.compile(r"ZZ\[([^\]]*)\]")


def mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


def tr(a):
    return [[a[j][i] for j in range(2)] for i in range(2)]


def main() -> None:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    start = next(i for i, line in enumerate(lines) if '@testset "Automormorphism group" begin' in line)
    end = next(i for i in range(start + 1, len(lines)) if lines[i].strip() == "end")
    explicit, counts, improper = [], [], []
    index = start + 1
    while index < end:
        form = FORM.match(lines[index])
        if not form:
            index += 1
            continue
        a, b, c = map(int, form.groups()[1:])
        assert b * b - 4 * a * c > 0, f"{SOURCE}:{index + 1}: form is not indefinite"
        block_end = next(k for k in range(index + 1, end + 1) if k == end or not lines[k].strip())
        case = {"form": [a, b, c], "gram": [[2 * a, b], [b, 2 * c]], "discriminant": b * b - 4 * a * c, "line": index + 1}
        kinds = 0
        for line in lines[index + 1 : block_end]:
            if match := GENS.match(line):
                generators = [[[int(x) for x in row.split()] for row in m.split(";")] for m in MATRIX.findall(match.group(1))]
                explicit.append(case | {"generators": generators})
                kinds += 1
            elif match := COUNT.match(line):
                counts.append(case | {"count": int(match.group(1))})
                kinds += 1
            elif IMPROPER.match(line):
                improper.append(case)
                kinds += 1
        assert kinds == 1, f"{SOURCE}:{index + 1}: expected one recorded fact, found {kinds}"
        index = block_end
    assert (len(explicit), len(counts), len(improper)) == (4, 6, 1), (len(explicit), len(counts), len(improper))
    row = all(mul(mul(t, case["gram"]), tr(t)) == case["gram"] for case in explicit for t in case["generators"])
    column = all(mul(mul(tr(t), case["gram"]), t) == case["gram"] for case in explicit for t in case["generators"])
    assert row or column, "the recorded generators are not isometries under either convention"
    convention = "T G T^T == G (row vectors)" if row else "T^T G T == G (column vectors)"

    def record(case: dict[str, object], **fields: object) -> dict[str, object]:
        form = case["form"]
        return {
            "id": f"hecke_binary_{'_'.join(map(str, form))}",
            "form": form,
            "gram": case["gram"],
            "discriminant": case["discriminant"],
            **fields,
            "source": {"kind": "independent_implementation_test", "file": str(SOURCE), "line": case["line"]},
        }

    data = {
        "explicit_generators": [record(case, automorphism_group_generators=case["generators"], isometry_convention=convention) for case in explicit],
        "generator_counts": [record(case, hecke_generator_count=case["count"]) for case in counts],
        "improper_automorphisms": [record(case, has_improper_automorphism=True) for case in improper],
    }
    TARGET.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {TARGET}: {len(explicit)} explicit, {len(counts)} counts, {len(improper)} improper; convention {convention}")


if __name__ == "__main__":
    main()
