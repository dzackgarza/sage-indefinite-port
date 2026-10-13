#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Extract tests/fixtures/enriques_classical_indices.json from Dutour Sikirić--Hulek.

Source: arXiv:2302.01679, subsection "Classical cases" (vendored TeX). For the
classical polarizations of degree 2 (double planes, Case 1), 6 (Enriques sextics,
Case 4) and 10 (Reye congruences, Case 7), the paper derives the index
[O^+(N) : Gamma_h] from Barth--Peters's counts of projective models, and then
[Gamma_h : tilde O^+(N)] = |O^+(F_2^10)| / [O^+(N) : Gamma_h], with
|O^+(F_2^10)| = 2^21 3^5 5^2 7 17 31 (line 700). Table table_subgroups1 records
|bar Gamma_h| = [Gamma_h : tilde O^+(N)] (enriques_87_polarizations.json).

Checks asserted here: for degrees 6 and 10 the displayed quotient, its stated value
and the table's |bar Gamma_h| agree. For degree 2 the text is not consistent: it
states [O^+(N) : Gamma_h] = 2^7 7 31 and [Gamma_h : tilde O^+(N)] = 2^14 3^5 5^2 7 31,
while the quotient of the displayed numbers is 2^14 3^5 5^2 17 and the table gives
2^14 3^5 5^2 7, i.e. index 2^7 17 31. The script asserts exactly this disagreement and
records the table-derived index alongside the stated values.

Run from the repository root: uv run references/extract/dh23_classical.py
"""

import json
import math
import re
from pathlib import Path

SOURCE = Path("references/vendor/arxiv/2302.01679/Enriques_compu_rev.tex")
TABLE = Path("tests/fixtures/enriques_87_polarizations.json")
TARGET = Path("tests/fixtures/enriques_classical_indices.json")
FACTOR = re.compile(r"^(\d+)(?:\^\{?(\d+)\}?)?$")
O_PLUS_F2 = 2**21 * 3**5 * 5**2 * 7 * 17 * 31


def product(text: str) -> int:
    values = []
    for factor in text.replace(" ", "").split(r"\cdot"):
        match = FACTOR.match(factor)
        assert match, f"unparsed factor {factor!r} in {text!r}"
        values.append(int(match.group(1)) ** int(match.group(2) or 1))
    return math.prod(values)


def cite(first: int, last: int, *needles: str) -> dict[str, str]:
    text = "\n".join(SOURCE.read_text(encoding="utf-8").splitlines()[first - 1 : last])
    for needle in needles:
        assert needle in text, f"{SOURCE}:{first}-{last} does not contain {needle!r}"
    return {"kind": "published_text", "file": str(SOURCE), "lines": f"{first}-{last}"}


def display(line: int) -> tuple[int, int, int]:
    """The numerator, divisor and stated value of the display [Gamma_h : tilde O^+(N)] = a / b = c."""
    text = SOURCE.read_text(encoding="utf-8").splitlines()[line - 1]
    match = re.fullmatch(r"\[\\Gamma_h: \\tilde\\Orth\^\+\(N\)\]= (.*?)\s*/\s*(.*?)= (.*?)\.?", text.strip())
    assert match, f"{SOURCE}:{line}: unparsed display {text!r}"
    return product(match.group(1)), product(match.group(2)), product(match.group(3))


def main() -> None:
    table = {row["case"]: row["group_order"] for row in json.loads(TABLE.read_text(encoding="utf-8"))}
    order = cite(700, 700, r"\left\vert \Orth^+(\FF_2^{10}) \right\vert = 2^{21} \cdot 3^5 \cdot 5^2 \cdot 7\cdot 17 \cdot 31")
    records = []
    for degree, case, count_lines, count_needle, display_line, case_needle, index in (
        (6, 4, (1291, 1293), r"$2^{10} \cdot 5 \cdot 17 \cdot 31$", 1295, "agreeing with Case 4", 2**10 * 5 * 17 * 31),
        (10, 7, (1305, 1307), r"$2^{13} \cdot 3 \cdot 17 \cdot 31$", 1309, "agrees exactly with Case 7", 2**13 * 3 * 17 * 31),
    ):
        numerator, divisor, stated = display(display_line)
        assert numerator == O_PLUS_F2 and divisor == index and numerator // divisor == stated == table[case], (degree, stated, table[case])
        records.append(
            {
                "degree": degree,
                "case": case,
                "index_in_O_plus_N": index,
                "gamma_h_over_stable_order": stated,
                "source": {
                    "count": cite(*count_lines, count_needle),
                    "display": cite(display_line, display_line + 2, case_needle),
                    "order_of_O_plus_F2": order,
                },
            }
        )
    numerator, divisor, stated = display(1278)
    cite(1275, 1275, r"admits  $2^7\cdot 7 \cdot 31$ different double plane representations", r"$[\Orth^+(N):\Gamma_h]=2^7\cdot 7 \cdot 31$")
    assert numerator == O_PLUS_F2 and divisor == 2**7 * 7 * 31
    assert stated == 2**14 * 3**5 * 5**2 * 7 * 31 and numerator // divisor == 2**14 * 3**5 * 5**2 * 17
    assert table[1] == 2**14 * 3**5 * 5**2 * 7 and O_PLUS_F2 // table[1] == 2**7 * 17 * 31
    inconsistent = [
        {
            "degree": 2,
            "case": 1,
            "index_in_O_plus_N": O_PLUS_F2 // table[1],
            "gamma_h_over_stable_order": table[1],
            "stated_index_in_O_plus_N": divisor,
            "stated_gamma_h_over_stable_order": stated,
            "quotient_of_displayed_numbers": numerator // divisor,
            "source": {
                "count": cite(1275, 1275, "Theorem 3.9"),
                "display": cite(1278, 1280, "This is Case 1"),
                "order_of_O_plus_F2": order,
            },
        }
    ]
    TARGET.write_text(json.dumps({"consistent_with_table": records, "inconsistent_with_table": inconsistent}, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {TARGET}: degrees 6 and 10 agree with the table; degree 2 text disagrees, table-derived index {O_PLUS_F2 // table[1]}")


if __name__ == "__main__":
    main()
