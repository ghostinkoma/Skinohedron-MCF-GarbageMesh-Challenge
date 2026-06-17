"""
s3d_sparameter.py
=================
3D extension, STEP 3 (the idea you proposed): an "as-equally-spaced-as-possible
S-parameter tetrahedron."

The most equally-spaced placement of 4 ports is the regular tetrahedron: the 4
vertices of a regular tetrahedron are the unique 4-point arrangement on a sphere
that maximises the minimum pairwise angle -- every pair subtends the SAME angle
arccos(-1/3) = 109.47 deg. We verify that (all 6 pairwise dot products = -1/3),
then build the scattering matrix of a 4-port junction carrying the full
tetrahedral symmetry (the group T_d), and check the physical invariants that
Part II's S-operator must satisfy:

    unitarity  (lossless)   : S^H S = I
    reciprocity             : S = S^T
    passivity  (|gain| = 1) : every eigenvalue lies on the unit circle

The Td-symmetric lossless reciprocal 4-port.
------------------------------------------
Full tetrahedral symmetry forces every self-term equal and every cross-term
equal, so

        S = e^{i a} P0  +  e^{i b} (I - P0) ,     P0 = (1/4) J  (J = all-ones),

i.e. S_ii = (e^{ia} + 3 e^{ib})/4 and S_ij = (e^{ia} - e^{ib})/4 for i != j.
This is automatically unitary, symmetric (reciprocal), and passive for any real
phases a, b. Its spectrum is

        e^{i a}  (the symmetric mode (1,1,1,1)/2,  multiplicity 1)
        e^{i b}  (the 3-D orthogonal complement,    multiplicity 3),

which is exactly the trivial + standard 3-D irreducible representations of T_d.
The 3-fold degeneracy IS the algebraic fingerprint of the maximal tetrahedral
symmetry -- the scattering analogue of the l=1 triplet seen in the ball
Laplacian (s3d_laplacian, family A).

A clean corollary (exact, not numerical): no choice of (a,b) makes all four
ports simultaneously matched (S_ii = 0), because that needs e^{ia} = -3 e^{ib},
impossible for unit-modulus phases. The symmetric lossless tetrahedral junction
is therefore never reflectionless -- a small, exact, classical impossibility,
recovered here from the symmetry alone.

Honest scope: this is a finite, exact, symmetry object -- a legitimate and tidy
"S-parameter tetrahedron". It does NOT by itself constitute a 3D field solver.
The bridge to the mesh work is described at the end (a TLM-like lattice of these
nodes on the uniform Kuhn/BCC tetrahedra), and is left as the next construction.
"""
from __future__ import annotations
import numpy as np


def regular_tetrahedron():
    """4 unit vertices of a regular tetrahedron (maximally equal-spaced ports)."""
    V = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
    return V / np.linalg.norm(V, axis=1, keepdims=True)


def equal_spacing_residual(V):
    """max | <vi,vj> - (-1/3) | over the 6 port pairs (0 => perfectly equal)."""
    worst = 0.0
    for i in range(4):
        for j in range(i + 1, 4):
            worst = max(worst, abs(float(V[i] @ V[j]) - (-1.0 / 3.0)))
    return worst


def td_scattering_matrix(a, b):
    """Td-symmetric lossless reciprocal 4-port:  S = e^{ia} P0 + e^{ib}(I - P0)."""
    J = np.ones((4, 4))
    P0 = J / 4.0
    return np.exp(1j * a) * P0 + np.exp(1j * b) * (np.eye(4) - P0)


def invariants(S):
    unit = np.linalg.norm(S.conj().T @ S - np.eye(4))
    recip = np.linalg.norm(S - S.T)
    gains = np.abs(np.linalg.eigvals(S))
    return float(unit), float(recip), float(gains.max()), float(gains.min())


def main():
    print("=== 3D-S3: an equally-spaced S-parameter tetrahedron ===\n")

    V = regular_tetrahedron()
    res = equal_spacing_residual(V)
    print("PORT GEOMETRY -- 4 ports at regular-tetrahedron vertices")
    print(f"   all 6 pairwise angles equal? |<vi,vj> + 1/3| <= {res:.2e}")
    print(f"   common angle = arccos(-1/3) = "
          f"{np.degrees(np.arccos(-1/3)):.2f} deg  (maximally equal-spaced)\n")

    print("SCATTERING -- Td-symmetric lossless reciprocal junction S(a,b)")
    print(f'   {"(a,b)":>16} {"unitarity":>11} {"reciprocity":>12} '
          f'{"max|gain|":>10} {"min|gain|":>10}')
    for a, b in [(0.0, 0.0), (0.7, -1.3), (1.1, 2.0), (np.pi, 0.0)]:
        S = td_scattering_matrix(a, b)
        u, r, gmax, gmin = invariants(S)
        print(f'   ({a:5.2f},{b:5.2f}) {u:11.2e} {r:12.2e} {gmax:10.6f} {gmin:10.6f}')
    print()

    # spectrum / degeneracy fingerprint
    S = td_scattering_matrix(0.7, -1.3)
    ev = np.linalg.eigvals(S)
    ang = np.sort(np.angle(ev))
    print("EIGEN-PHASES of S(0.7,-1.3) (group-theory fingerprint):")
    print(f"   phases (rad) = {np.round(ang, 4)}")
    print(f"   => one phase = a = 0.70 (symmetric mode, mult 1),")
    print(f"      three phases = b = -1.30 (mult 3) : trivial + 3-D irrep of Td\n")

    # matched-port impossibility (exact)
    diag = np.abs(np.diag(td_scattering_matrix(np.pi, 0.0)))
    print("MATCHED-PORT CHECK -- can all four ports be reflectionless (S_ii=0)?")
    print(f"   best symmetric attempt (a=pi,b=0): |S_ii| = {diag[0]:.4f} (> 0)")
    print("   exact reason: S_ii=0 needs e^{ia} = -3 e^{ib}, impossible for unit")
    print("   phases. The symmetric lossless tetrahedral junction is never")
    print("   reflectionless -- recovered from symmetry alone.\n")

    print("[結句] A concrete, exact 'S-parameter tetrahedron': 4 maximally")
    print("   equal-spaced ports, a lossless reciprocal junction whose unitarity")
    print("   and reciprocity hold to machine precision for any phases, and whose")
    print("   1+3 eigen-degeneracy is the algebraic signature of full tetrahedral")
    print("   symmetry -- the scattering twin of the Laplacian's l=1 triplet.")
    print("   This is a finite symmetry object, not yet a field solver. The")
    print("   natural bridge to the mesh work: put one such node on every cell of")
    print("   the UNIFORM Kuhn/BCC tetrahedral lattice (s3d_laplacian, family C),")
    print("   ports on the 4 shared faces, to get an isotropic TLM-like 3D")
    print("   scattering network. Because the cells are congruent and equally")
    print("   spaced, every node is identical and the lattice is isotropic by")
    print("   construction. Building and verifying that assembled network")
    print("   (energy conservation across the whole lattice, dispersion vs the")
    print("   exact wave speed) is the next step -- and is where any genuine")
    print("   content beyond standard TLM / FE-exterior-calculus would have to be")
    print("   demonstrated, not assumed.")


if __name__ == "__main__":
    main()
