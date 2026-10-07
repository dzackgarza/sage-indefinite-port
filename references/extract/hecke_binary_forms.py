#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Extract tests/fixtures/binary_form_automorphisms.json from Hecke's QuadBin tests.

Source: references/vendor/Hecke.jl@e2ab5716/test/QuadForm/QuadBin.jl, testset
"Automormorphism group". For four indefinite binary forms a x^2 + b x y + c y^2 the test
records Hecke's exact generators of the automorphism group. The lattice of the form has
Gram matrix [[2a, b], [b, 2c]]. Cases where the test records only how many generators
Hecke returns are not included: that count is an implementation detail, not a property
of the group.

The script asserts each form is indefinite (b^2 - 4ac > 0), and determines the action
convention from the data: every recorded generator T must satisfy T G T^T = G for all
cases, or T^T G T = G for all cases; the convention that holds is recorded.

Run from the repository root: uv run references/extract/hecke_binary_forms.py
"""

import json
import re
from pathlib import Path

SOURCE = Path("references/vendor/Hecke.jl@e2ab5716/test/QuadForm/QuadBin.jl")
TARGET = Path("tests/fixtures/binary_form_automorphisms.json")
FORM = re.compile(r"g = binary_quadratic_form\((-?\d+), (-?\d+), (-?\d+)\)")
GENS = re.compile(r"@test gens == \[(.*)\]")
MATRIX = re.compile(r"ZZ\[([^\]]*)\]")


def mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


def tr(a):
    return [[a[j][i] for j in range(2)] for i in range(2)]


def main() -> None:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    cases = []
    for index, line in enumerate(lines):
        form = FORM.search(line)
        if not form:
            continue
        a, b, c = map(int, form.groups())
        recorded = next((GENS.search(lines[k]) for k in range(index + 1, index + 3) if GENS.search(lines[k])), None)
        if recorded is None:
            continue
        generators = [[[int(x) for x in row.split()] for row in m.split(";")] for m in MATRIX.findall(recorded.group(1))]
        assert b * b - 4 * a * c > 0, f"{SOURCE}:{index + 1}: form is not indefinite"
        cases.append({"form": [a, b, c], "gram": [[2 * a, b], [b, 2 * c]], "generators": generators, "line": index + 1})
    assert len(cases) == 4, len(cases)
    row = all(mul(mul(t, case["gram"]), tr(t)) == case["gram"] for case in cases for t in case["generators"])
    column = all(mul(mul(tr(t), case["gram"]), t) == case["gram"] for case in cases for t in case["generators"])
    assert row or column, "the recorded generators are not isometries under either convention"
    convention = "T G T^T == G (row vectors)" if row else "T^T G T == G (column vectors)"
    records = [
        {
            "id": f"hecke_binary_{'_'.join(map(str, case['form']))}",
            "form": case["form"],
            "gram": case["gram"],
            "discriminant": case["form"][1] ** 2 - 4 * case["form"][0] * case["form"][2],
            "automorphism_group_generators": case["generators"],
            "isometry_convention": convention,
            "source": {"kind": "independent_implementation_test", "file": str(SOURCE), "line": case["line"]},
        }
        for case in cases
    ]
    TARGET.write_text(json.dumps(records, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {len(records)} forms to {TARGET}; convention {convention}")


if __name__ == "__main__":
    main()
