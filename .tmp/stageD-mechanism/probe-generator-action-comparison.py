"""Locate where the -e sign generator is lost between the rational group and the integral stabilizer.

Specimens: U (control: sign delivered), U(2) (control: sign delivered), U+U(2) (sign lost, measured earlier).
At every layer print the actual tensor so semantics is read, not inferred:
  G[i]   = ambient rational tensor of rational-group generator i
  S[i]   = ambient rational tensor of stabilizer generator i, and its integral_restriction rows
Also: membership of the sign element in the stabilizer, and the INDUCTION end-perception.
"""
import time
import traceback

START = time.monotonic()


def stage(label):
    print(f"[t={time.monotonic()-START:8.1f}s] {label}", flush=True)


stage("imports")
from tests.acceptance.consumer import element, lattice, rows_of, is_isometry
from sage_indefinite_port.indefinite.recursive import EichlerRecursion
from sage_indefinite_port.indefinite.marked_lattices import INDUCTION


def srows(map_, space):
    labels = tuple(space.module_generating_set())
    rows = []
    for label in labels:
        vec = map_(space.module_generator(label)).to_vector()
        rows.append([str(vec(target)) for target in labels])
    return rows


def sign_of(rows):
    diag = tuple(rows[i][i] for i in range(len(rows)))
    offdiag_zero = all(rows[i][j] == "0" for i in range(len(rows)) for j in range(len(rows)) if i != j)
    return diag if offdiag_zero else None


def examine(name, gram, e_coords):
    stage(f"=== {name} ===")
    try:
        L = lattice(gram)
        e = element(L, e_coords)
        REC = EichlerRecursion()
        t0 = time.monotonic()
        action = REC.isotropic_stabilizer_action((e,), ((-e,),))
        stage(f"{name}: isotropic_stabilizer_action {time.monotonic()-t0:.1f}s")
        G = action.rational_group()
        stage(f"{name}: G has_selected_group_resolution={G.has_selected_group_resolution()} "
              f"has_selected_finite_group_generating_set={G.has_selected_finite_group_generating_set()}")
        inc = action._ambient_group_inclusion()
        V = action.lattice_inclusion().codomain().module_over_extension()
        eQ = action.lattice_inclusion()(e).underlying_element()
        sup = tuple(G.group_generators())
        stage(f"{name}: rational group: {len(sup)} generators (ambient tensors):")
        for i, g in enumerate(sup):
            rows = srows(inc(g), V)
            e_image = inc(g)(eQ)
            e_coords_out = [str(e_image.to_vector()(label)) for label in tuple(V.module_generating_set())]
            print(f"  G[{i}] diag={sign_of(rows)} rows={rows} e|->{e_coords_out}", flush=True)
        t0 = time.monotonic()
        stab = action.stabilizer()
        sgens = tuple(stab.group_generators())
        stage(f"{name}: action.stabilizer() {time.monotonic()-t0:.1f}s: {len(sgens)} generators")
        for i, s in enumerate(sgens):
            arows = srows(inc(s), V)
            try:
                integ = L.integral_isometry(s)
                irows = rows_of(integ)
                e_out = [int(x) for x in integ(e).to_vector()]
            except Exception as exc:
                irows = f"<{type(exc).__name__}: {exc}>"
                e_out = "?"
            print(f"  S[{i}] ambient_diag={sign_of(arows)} rows={arows}", flush=True)
            print(f"       integral={irows} e|->{e_out}", flush=True)
        for i, g in enumerate(sup):
            e_image = inc(g)(eQ)
            e_coords_out = [str(e_image.to_vector()(label)) for label in tuple(V.module_generating_set())]
            if e_coords_out and e_coords_out[0].startswith("-"):
                try:
                    member = stab.is_parent_of(g)
                except Exception as exc:
                    member = f"<{type(exc).__name__}: {exc}>"
                try:
                    integ = L.integral_isometry(g)
                    irows = rows_of(integ)
                except Exception as exc:
                    irows = f"<{type(exc).__name__}: {exc}>"
                print(f"  SIGN-CHECK G[{i}]: in_stabilizer={member} integral_restriction_rows={irows}", flush=True)
        t0 = time.monotonic()
        gens = tuple(INDUCTION.line_automorphism_generators(L, e))
        signs = []
        for generator in gens:
            image = generator(e)
            signs.append(1 if image == e else (-1 if image == -e else None))
        stage(f"{name}: INDUCTION.line_automorphism_generators {time.monotonic()-t0:.1f}s: {len(gens)} gens, line signs {signs}")
        stage(f"{name}: all preserve form: {all(is_isometry(rows_of(g), gram) for g in gens)}")
    except Exception:
        traceback.print_exc()


examine("U", [[0, 1], [1, 0]], [1, 0])
examine("U2", [[0, 2], [2, 0]], [1, 0])
examine("U+U2", [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 2], [0, 0, 2, 0]], [1, 0, 0, 0])
stage("PROBE_DONE")
