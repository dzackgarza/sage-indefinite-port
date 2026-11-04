#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Regenerate the published Conway--Sloane ternary pair from an equation extract.

The committed source is a short transcription of the original displayed
mathematical data (Ch. 15, Sec. 11, (51a)/(51b)), not a copy of the book.
"""

import json
import re
from pathlib import Path

SOURCE = Path("references/source-data/conway_sloane_ch15_sec11_51ab.tex")
TARGET = Path("tests/fixtures/conway_sloane_cases.json")


def main() -> None:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    factual_note = "\n".join(line.removeprefix("% ") for line in lines if line.startswith("% "))
    assert "two spinor genera and two classes" in factual_note
    assert "determinant -128" in factual_note
    display = "\n".join(lines)
    # Parse the two displayed matrix environments, not the preceding binary matrices.
    pair = re.search(r"\$\$\s*\( a \).*?\$\$", display, re.S)
    assert pair is not None, "missing (51a)/(51b) displayed pair"
    environments = re.findall(r"\\begin\{array\} \{ l l l \}(.*?)\\end\{array\}", pair.group(), re.S)
    assert len(environments) == 2, f"expected two ternary forms, got {len(environments)}"

    def parse_matrix(env: str) -> list[list[int]]:
        rows = []
        for row in re.split(r"\\\\", env):
            entries = [re.sub(r"\\[{}]", "", x).strip(" {}$") for x in row.split("&")]
            assert len(entries) == 3, entries
            rows.append([int(x.replace(" ", "")) for x in entries])
        assert len(rows) == 3
        return rows

    a, b = map(parse_matrix, environments)

    def det3(m: list[list[int]]) -> int:
        return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])

    for gram in (a, b):
        assert gram == [list(row) for row in zip(*gram)]
        assert det3(gram) == -128
        assert gram[2][2] == 2 and gram[0][0] < 0 and gram[1][1] > 0
        assert gram[0][0] * gram[1][1] - gram[0][1] ** 2 == -64
    record = {
        "spinor_genus_pair_determinant_minus_128": {
            "id": "cs99_ch15_sec11_forms_51ab",
            "genus": "I_{2,1}(2 x 64)",
            "determinant": -128,
            "signature": [2, 1],
            "form_a": {"equation": "(51a)", "gram": a},
            "form_b": {"equation": "(51b)", "gram": b},
            "same_genus": True,
            "spinor_genera_in_genus": 2,
            "classes_in_genus": 2,
            "integrally_equivalent": False,
            "source": {
                "kind": "published_text",
                "citation": "J. H. Conway and N. J. A. Sloane, Sphere Packings, Lattices and Groups, 3rd ed., Springer, 1999",
                "location": "Chapter 15, Section 11, equations (51a) and (51b)",
            },
        }
    }
    TARGET.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"extracted {TARGET} from {SOURCE}")


if __name__ == "__main__":
    main()
