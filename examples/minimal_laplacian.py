"""
minimal_laplacian.py
====================
Shortest possible example: compute the discrete Laplacian of a spherical
harmonic on the KSF icosphere and check the eigenvalue.

Run from the repo root:
    python3 examples/minimal_laplacian.py

Expected output (approximately):
    Level 3 — 642 vertices, h = 0.1507
    Exact eigenvalue:    6.000000
    Discrete Rayleigh:   6.000082
    Relative error:      0.0014%
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from ksf.mesh import icosphere, mean_edge_length
from ksf.dec  import cotangent_laplacian, mass_matrix_sparse
from ksf.sph  import harmonic

# ── build mesh ────────────────────────────────────────────────────────────────
V, F = icosphere(level=3)          # 642 vertices, h ≈ 0.15
h    = mean_edge_length(V, F)

# ── discrete operators ────────────────────────────────────────────────────────
# L  : sparse (N×N) stiffness, approximates −Δ (positive semidefinite)
# M  : sparse (N×N) diagonal lumped-mass matrix
L, M = cotangent_laplacian(V, F)
M    = mass_matrix_sparse(M)       # convert to sparse for @ operator

# ── test function: Y₂ₓᵧ, exact eigenvalue λ = 6 ─────────────────────────────
f, lam_exact = harmonic(V, 'Y2xy')

# ── discrete Rayleigh quotient: λ_h = (f·Lf) / (f·Mf) ───────────────────────
lam_h = float(f @ L @ f) / float(f @ M @ f)

print(f"Level 3 — {len(V)} vertices, h = {h:.4f}")
print(f"Exact eigenvalue:    {lam_exact:.6f}")
print(f"Discrete Rayleigh:   {lam_h:.6f}")
print(f"Relative error:      {abs(lam_h - lam_exact)/lam_exact * 100:.4f}%")
