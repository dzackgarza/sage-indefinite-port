# sage-extensions: a `description` f-string in `_symmetric_group` materialises a whole coset space

Pinned `sage-extensions` **`e83c4d69db6fba9cadc560d6fd8a518bfe6c5dae`**, line numbers its own.
Found from an instrumented run of `test_boundaries.sage::test_enriques_modular_groups_match_dutour_sikiric_hulek[case_1]`,
which reached its 1800 s cap still inside one non-returning call. No library code was changed here.

## The call that does not return

At the cap the test is at `test_boundaries.sage:226` —
`target.subgroup_on(tuple(g.discriminant_morphism() for g in group.group_generators()))` — and the
innermost frames of the last complete dump (76 frames, t≈1740 s) are, innermost first:

```
sage/groups/abelian_gps/abelian_group_gap.py:152   exponents
sage/groups/abelian_gps/abelian_aut.py:203         matrix
sage/groups/fqf_orthogonal.py:484                  _act_
sage/groups/fqf_orthogonal.py:288                  _preserves_form   (<genexpr>)
sage/groups/fqf_orthogonal.py:267                  _element_constructor_
sage_extensions/.../group/groups.py:305            _element_from_engine
sage_extensions/.../group/groups.py:525            representative
sage_extensions/.../group/groups.py:472            representative
sage_extensions/.../group/groups.py:486            _TransversalCoset._repr_
sage_extensions/.../sets/finite_ordered_sets.py:587 <genexpr>
sage_extensions/.../sets/finite_ordered_sets.py:587 _repr_
sage_extensions/.../group/groups.py:2133           _symmetric_group
```

## Root cause

`finite_ordered_sets.py:587`:

```python
def _repr_(self) -> str:
    return "{" + ", ".join(repr(element) for element in self) + "}"
```

`groups.py:486` (`_TransversalCoset`):

```python
def _repr_(self) -> str:
    return f"{self.representative()} * {self.parent().subgroup}"
```

`groups.py:2124-2134`:

```python
    return _own_group(
        engine,
        refinements=refinements,
        description=(
            f"Symmetric group S_{size}"
            if standard
            else f"Symmetric group on {points}"      # <- formats the point set, eagerly
        ),
        symmetric_points=points,
        ...
    )
```

The `description` argument is an eager f-string. On the non-standard branch it formats `points`, a
finite ordered set of transversal cosets. Its `_repr_` iterates every point, each coset's `_repr_`
calls `representative()`, and that materialises a real group element through `_element_from_engine`
into Sage's `FqfOrthogonalGroup._element_constructor_`, which runs `_preserves_form` → `_act_` →
`AbelianGroupAutomorphism.matrix` → `abelian_group_gap.exponents` for every point.

So **building a symmetric group on a non-standard point set costs a full materialisation and
form-validation of that point set, purely to compose a description string that is usually never read.**

## Measured

Phase counters from the same run, wrapped around the provider's own methods:

| phase | calls | total | note |
| --- | --- | --- | --- |
| `group_generators` | 21 | 177.52 s | frozen after t≈380 s |
| `discriminant_morphism` | 21 | 13.75 s | frozen after t≈380 s |
| `cardinality` | **1 570 390** | 9.95 s | still growing at the cap |
| `subgroup_on` (a nested one) | 1 | 0.03 s | **not** the line-226 call, which never returned |

Of 1800 s only ~192 s is inside those wrapped phases; the remainder is the repr chain above. The
1.57 M `cardinality` calls at ~6 µs each are the signature of iterating the point set.

## Expected contract

Composing a group's `description` must not materialise the group's points, and `repr` of a finite
ordered set of cosets must not construct representatives. Either defer the description (a thunk
evaluated only when asked for) or format a bounded, non-materialising summary such as the cardinality
and the parent, and give `_TransversalCoset._repr_` a positional form that does not call
`representative()`.

## Minimal reproduction

Any construction of a symmetric group on a **non-standard** point set whose points are transversal
cosets of a finite quadratic form's orthogonal group reproduces it. The self-contained hot half is
the repr:

```python
# provider-only: repr of a coset set materialises every representative
from sage_extensions.categories.group.groups import _symmetric_group
# points = finite ordered set of _TransversalCoset over a TorsionFormOrthogonalGroup subgroup
_symmetric_group(points)     # never needs the description, yet formats it
repr(points)                 # time this directly: it materialises every point
```

A faithful instance is `.tmp/probe-gamma-h-stabilizer.py` in the port repository, which rebuilds the
Enriques case-1 group; `repr()` on the coset space reached from `predicate_subgroups.coset_space`
shows the same chain without the 1800 s wait.
