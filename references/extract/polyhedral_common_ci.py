"""Extract the polyhedral_common CI fixtures from the vendored upstream data.

Every value here is the reference implementation's own recorded output (or input),
read from references/vendor/polyhedral_common@1592b246/CI_tests. For a port,
agreement with it is the acceptance condition. GAP literals are evaluated with
Sage's libgap; nothing is built or run from the C++ sources.

Each record carries the vendored file and its index in that file. The script
asserts the invariants that the data itself must satisfy: recorded isometries
and generators preserve the forms, root vectors define integral reflections, and
the paired files describe the same lattices in the same order.

Run from the repository root under Sage's Python:
    "$(dirname $(sage -c 'import sys; print(sys.executable)'))/python3" references/extract/polyhedral_common_ci.py
"""

import json
import tarfile
from pathlib import Path

# Importing from sage.all runs Sage's session startup, which libgap needs.
from sage.all import ZZ, libgap, matrix

CI = Path("references/vendor/polyhedral_common@1592b246/CI_tests")
FIXTURES = Path("tests/fixtures")


def gap_file(path: Path):
    text = path.read_text(encoding="utf-8").strip()
    assert text.startswith("return"), f"{path} does not start with 'return'"
    return libgap.eval(text[len("return") :].strip().rstrip(";"))


def int_matrix(value) -> list[list[int]]:
    return [[int(entry) for entry in row] for row in value.sage()]


def source(path: Path, index: int | None = None) -> dict[str, object]:
    record: dict[str, object] = {"kind": "reference_implementation_output", "file": str(path)}
    if index is not None:
        record["index"] = index
    return record


def write(name: str, data: object) -> None:
    (FIXTURES / name).write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {FIXTURES / name}")


def reflective_and_isotropic() -> list[list[list[int]]]:
    reflect_path = CI / "20_Reflective" / "ListReflect"
    isotropic_path = CI / "DATA" / "IsotropicCases"
    reflect = gap_file(reflect_path)
    isotropic = gap_file(isotropic_path)
    assert reflect.Length() == isotropic.Length() == 8821
    reflective, cases, grams = [], [], []
    for index in range(int(reflect.Length())):
        gram = int_matrix(reflect[index]["LorMat"])
        assert int_matrix(isotropic[index]["M"]) == gram, f"IsotropicCases[{index}] is not ListReflect[{index}]"
        grams.append(gram)
        reflective.append(
            {
                "id": f"ListReflect_{index}",
                "dimension": len(gram),
                "gram": gram,
                "num_simple_roots": int(reflect[index]["n_simple"]),
                "source": source(reflect_path, index),
            }
        )
        cases.append(
            {
                "id": f"IsotropicCases_{index}",
                "dimension": len(gram),
                "gram": gram,
                "has_isotropic": bool(isotropic[index]["has_isotropic"]),
                "source": source(isotropic_path, index),
            }
        )
    write("reflective_forms_8821.json", reflective)
    write("isotropic_cases_8821.json", cases)
    return grams


def root_systems(grams: list[list[list[int]]], reflective: list[dict[str, object]]) -> None:
    path = CI / "01_RatIntAutomorphy" / "ListSimpleRootSystem_4_56_X_5_47"
    systems = gap_file(path)
    paired = [index for index, gram in enumerate(grams) if len(gram) in (4, 5)]
    assert systems.Length() == len(paired) == 103
    records = []
    for position, index in enumerate(paired):
        gram = matrix(ZZ, grams[index])
        roots = int_matrix(systems[position])
        assert len(roots) == reflective[index]["num_simple_roots"], f"root system {position} vs ListReflect[{index}]"
        for root in roots:
            r = matrix(ZZ, [root])
            square = (r * gram * r.transpose())[0, 0]
            assert square != 0 and all((2 * entry) % square == 0 for entry in (gram * r.transpose()).list()), (
                f"root system {position}: {root} does not define an integral reflection of ListReflect[{index}]"
            )
        records.append(
            {
                "id": f"ListSimpleRootSystem_{position}",
                "dimension": len(grams[index]),
                "gram": grams[index],
                "gram_source": source(CI / "20_Reflective" / "ListReflect", index),
                "num_roots": len(roots),
                "roots": roots,
                "source": source(path, position),
            }
        )
    write("root_systems_103.json", records)


def perfect_domains() -> None:
    path = CI / "28B_LorentzianPerfStabEqui" / "Result_Enumeration"
    entries = gap_file(path)
    by_gram: dict[str, dict[str, object]] = {}
    for index in range(int(entries.Length())):
        data, mode, count = entries[index][0], str(entries[index][1]), int(entries[index][2])
        gram = int_matrix(data["M"])
        key = json.dumps(gram)
        record = by_gram.setdefault(
            key,
            {
                "id": f"Result_Enumeration_{len(by_gram)}",
                "dimension": len(gram),
                "gram": gram,
                "has_isotropic": bool(data["has_isotropic"]),
                "source": source(path),
            },
        )
        assert mode in ("isotropic", "total") and f"{mode}_count" not in record
        record[f"{mode}_count"] = count
    records = list(by_gram.values())
    assert len(records) == 40 and all("isotropic_count" in r and "total_count" in r for r in records)
    write("lorentzian_perfect_domains.json", records)


def lorentzian_equivalences_and_stabilizers() -> None:
    equi_tar = CI / "28B_LorentzianPerfStabEqui" / "TestCasesEqui.tar.gz"
    stab_tar = CI / "28B_LorentzianPerfStabEqui" / "TestCasesStab.tar.gz"
    extracted = CI / "28B_LorentzianPerfStabEqui" / "x"
    for archive, folder in ((equi_tar, "TestCasesEqui"), (stab_tar, "TestCasesStab")):
        with tarfile.open(archive) as tar:
            members = sorted(m.name.split("/")[-1] for m in tar.getmembers() if m.isfile())
        on_disk = sorted(p.name for p in (extracted / folder).iterdir())
        assert members == on_disk, f"{extracted / folder} does not match {archive}"
    equivalences = []
    for path in sorted((extracted / "TestCasesEqui").iterdir(), key=lambda p: int(p.name.removeprefix("LORENTZ_TestIsomorphism"))):
        data = gap_file(path)
        mat1, mat2, witness = (int_matrix(data[key]) for key in ("mat1", "mat2", "test"))
        u, a, b = matrix(ZZ, witness), matrix(ZZ, mat1), matrix(ZZ, mat2)
        assert u * a * u.transpose() == b, f"{path}: test * mat1 * test^T != mat2"
        equivalences.append(
            {
                "id": path.name,
                "dimension": len(mat1),
                "mat1": mat1,
                "mat2": mat2,
                "transporter_witness": witness,
                "witness_convention": "witness * mat1 * witness^T == mat2",
                "source": source(path),
            }
        )
    write("lorentzian_equivalence_145.json", equivalences)
    stabilizers = []
    for path in sorted((extracted / "TestCasesStab").iterdir()):
        data = gap_file(path)
        gram = int_matrix(data["mat"])
        generators = [int_matrix(g) for g in data["GRPlor"].GeneratorsOfGroup()]
        g = matrix(ZZ, gram)
        for generator in generators:
            m = matrix(ZZ, generator)
            assert m * g * m.transpose() == g, f"{path}: a recorded generator does not preserve the form"
        stabilizers.append({"id": path.name, "dimension": len(gram), "gram": gram, "generators": generators, "num_generators": len(generators), "source": source(path)})
    write("lorentzian_stabilizers_cases.json", stabilizers)


def main() -> None:
    grams = reflective_and_isotropic()
    reflective = json.loads((FIXTURES / "reflective_forms_8821.json").read_text())
    root_systems(grams, reflective)
    perfect_domains()
    lorentzian_equivalences_and_stabilizers()


if __name__ == "__main__":
    main()
