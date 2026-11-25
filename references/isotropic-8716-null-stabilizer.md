# Pointed null stabilizer in IsotropicCases_8716

Let L have the frozen Gram matrix of IsotropicCases_8716. Write h=e1+e4 and f=(e1-e4)/2. Then h and f are rational null vectors with b(h,f)=1. The orthogonal complement K of their span has the integral basis e2,e3,e5,e6 and Gram diag(1,1) plus [[10,5],[5,10]]. Let K_even be its even-norm sublattice.

## Full group

The complete integral stabilizer fixing h is

    Stab_O(L)(h) = K_even semidirect O(K).

Every isometry of K extends integrally by fixing e1 and e4. To determine the kernel of the action on h-perp/Zh, suppose g fixes h and induces the identity on K. The equation b(h,g(f))=1 gives

    g(f) = f + t + a h.

The norm equation gives a=-b(t,t)/2. Orthogonality then gives

    g(k) = k - b(k,t) h     for every k in K.

Thus g is the unique null-line translation n_t. Since e1=h/2+f, integrality of g(e1) forces t in K and b(t,t) even. Conversely these conditions make n_t and n_-t integral on every displayed basis vector. The transformations n_t add their parameters. This proves both inclusions and identifies the complete integral kernel, without discarding elements invisible in a finite image.

The even-norm condition is t1+t2=0 modulo two. Its Hermite basis is

    [1 1 0 0]
    [0 2 0 0]
    [0 0 1 0]
    [0 0 0 1].

The norm-one vectors of K are exactly ±e2 and ±e3. Their span is invariant under every isometry. Its orthogonal complement is the scaled positive A2 lattice. The signed permutation group on the two unit axes and the dihedral action on the six minimal A2 vectors give the entire orthogonal group of K.

| Finite factor | Order |
| --- | ---: |
| Orthogonal group of the two unit axes | 8 |
| Orthogonal group of the scaled positive A2 summand | 12 |
| O(K) | 96 |

## Executable complete certificate

`tests/test_isotropic_8716_pointed_stabilizer.py` obtains the public coslice automorphism generators of the pointed lattice (L,h). It checks each generator's integral inverse, form and marking, then their full quotient image on K. The finite image has exactly the independently determined order above.

For every finite image element, GAP supplies a word in the original generator family. The assertion verifies the image of each word before using it as a lift. The complete finite Schreier family then generates the kernel. Every resulting integral action is checked against the explicit translation formula, and the integral span of all its parameters must equal the displayed Hermite module. Surjectivity onto O(K), together with equality of the entire kernel, proves equality with the full pointed stabilizer.

The complete public assertion passes with provider cfee1892 in 43.33 seconds.

Provider cfee1892 retains the original generator frame through its explicit free-group word map. The owner word assertion fails before this repair and passes afterward. The finite word map supplies untrusted candidate words here; the exact quotient and integral action checks certify each use.

This is the stabilizer of the supplied primitive null vector. Complete orbit classification and stabilizers for any other primitive null classes of IsotropicCases_8716 remain required.
