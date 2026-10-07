#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Extract tests/fixtures/allcock_rank3_reflective.json from Allcock's rank-3 table.

Source: references/vendor/GeometryDatabase_Rank3_Lorentzian_lattices@63a9067a/RK3_all,
Allcock's classification of the reflective Lorentzian lattices of rank 3 (arXiv:1010.0486,
1111.1264, vendored). Each line is [Gram, elementary divisors, simple roots as columns,
Weyl-group id, lattice id], in Allcock's convention: signature (2,1), roots of positive
norm. This is the same sign convention as polyhedral_common's ListReflect.

The script asserts, for every lattice, that each simple root has positive norm and
defines an integral reflection; and that wherever a Gram matrix also occurs in
reflective_forms_8821.json (polyhedral_common's ListReflect), the two sources give
the same number of simple roots.

Run from the repository root: uv run references/extract/allcock_rank3.py
"""

import json
import re
from pathlib import Path

SOURCE = Path("references/vendor/GeometryDatabase_Rank3_Lorentzian_lattices@63a9067a/RK3_all")
LIST_REFLECT = Path("tests/fixtures/reflective_forms_8821.json")
TARGET = Path("tests/fixtures/allcock_rank3_reflective.json")
LINE = re.compile(r"^\[\[([^\]]*)\],\[([^\]]*)\],\[([^\]]*)\],(\d+),(\d+)\]$")


def pari_matrix(text: str) -> list[list[int]]:
    return [[int(entry) for entry in row.split(",")] for row in text.split(";")]


def main() -> None:
    reflect_counts = {
        json.dumps(record["gram"]): record["num_simple_roots"]
        for record in json.loads(LIST_REFLECT.read_text())
        if record["dimension"] == 3
    }
    records, agreements = [], 0
    for number, line in enumerate(SOURCE.read_text(encoding="utf-8").splitlines(), start=1):
        match = LINE.match(line.strip())
        assert match, f"{SOURCE}:{number}: unparsed line {line[:80]!r}"
        gram = pari_matrix(match.group(1))
        divisors = [int(entry) for entry in match.group(2).split(",")]
        root_matrix = pari_matrix(match.group(3))
        roots = [[root_matrix[i][j] for i in range(3)] for j in range(len(root_matrix[0]))]
        for root in roots:
            pairings = [sum(gram[i][k] * root[k] for k in range(3)) for i in range(3)]
            square = sum(root[i] * pairings[i] for i in range(3))
            assert square > 0, f"{SOURCE}:{number}: root {root} has norm {square}"
            assert all((2 * p) % square == 0 for p in pairings), f"{SOURCE}:{number}: root {root} gives no integral reflection"
        key = json.dumps(gram)
        if key in reflect_counts:
            assert reflect_counts[key] == len(roots), f"{SOURCE}:{number}: ListReflect gives {reflect_counts[key]} roots, Allcock {len(roots)}"
            agreements += 1
        records.append(
            {
                "id": f"allcock_rk3_{int(match.group(5))}",
                "gram": gram,
                "elementary_divisors": divisors,
                "simple_roots": roots,
                "num_simple_roots": len(roots),
                "weyl_group_id": int(match.group(4)),
                "lattice_id": int(match.group(5)),
                "convention": "signature (2,1), simple roots of positive norm",
                "source": {"kind": "published_classification", "file": str(SOURCE), "line": number},
            }
        )
    assert len(records) == 8595, len(records)
    lines = ",\n".join(json.dumps(record, separators=(",", ":"), ensure_ascii=False) for record in records)
    TARGET.write_text(f"[\n{lines}\n]\n", encoding="utf-8")
    print(f"wrote {len(records)} lattices to {TARGET}; {agreements} Gram matrices also in ListReflect, all with equal root counts")


if __name__ == "__main__":
    main()
