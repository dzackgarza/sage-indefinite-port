# Fractional-span integral lifts

Provider `16c798dbea898555a1373686bec1bb893edb973a` retains the selected finite coefficient family when lifting through a restricted fraction-field span. The original O_plus_N proof at a170e2e/0580f98d reaches this lift through `FiniteCommensurabilityQuotient.restricted_automorphism`, during its inverse-preservation assertion. That active proof predates the change.

The owner's `_selected_module_coefficients` contract states that coefficients outside the finite family are zero; its keys contain the support without asserting that every retained coefficient is nonzero. The lift checks that all retained coefficients outside the finite span window vanish. A false result returns no preimage. An undecided result raises. This is the same coordinate test already used by `lift_from_finite_support` in the same file. The native Sage span membership and coordinate solver are unchanged.

The independent example is L=Z(e/2)+Z((e+f)/3). A vector x*e+y*f belongs to L exactly when 3*y and 2*(x-y) are integers, and every other coordinate vanishes. Each accepted lift must map back to its actual owned input.

| Owner assertion | Result | Call seconds |
| --- | --- | ---: |
| Fractional span in a finite named basis | pass | 0.39 |
| Fractional span in an infinite integer-labelled basis | pass | 0.15 |
| Least stable rational orbit span, including bounds and seed factorization | pass | 2.07 |

The tests retain the original owner limits. The span examples include zero, fractional exclusions and a nonzero coordinate outside the support window. Their baseline at cab1d797 also passes; the change preserves the existing mathematical operation.

An isolated measurement prepares fresh owned vectors before each timed interval and asserts every public membership answer. Each repetition tests 200 vectors `(i/2+1/3)*e+(1/3)*f` in the same fractional span.

| Lift implementation | Repetition seconds |
| --- | --- |
| Categorical support construction | 2.4895, 2.1291, 2.2183 |
| Retained coefficient family | 0.7764, 0.6293, 0.6518 |

This measures the membership operation, not an entire lattice consumer. The port pin at `583b36d625bcbadcf1ef397b4140230ce3fa300e` selects this provider and passes `uv lock --check` in the shared Sage container. The full A7, O_plus_N and stable Dawes assertions at the final pair remain required. Provider #294 comment 6718329418 records the input contract and results. Local timing/test logs are optional diagnostics under `.tmp/fractional-span-*`.


At port 583b36d625bcbadcf1ef397b4140230ce3fa300e/provider 16c798dbea898555a1373686bec1bb893edb973a, the unchanged full spanning-cell transporter assertion passes with JUnit case time 33.137 seconds and total time 38.63 seconds. The supplied integral unimodular witness and complete configuration equality are checked. This is additional complete small-consumer evidence for the coefficient-family lift. Its timing does not isolate this change from other wall-time variation. The full A7 assertion remains required.
