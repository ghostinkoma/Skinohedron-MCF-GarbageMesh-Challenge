"""
test_ksf.py
===========
Basic sanity tests for the ksf library.

Run from the repo root:
    python3 tests/test_ksf.py
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import numpy as np
from ksf.mesh       import icosphere, mean_edge_length
from ksf.dec        import cotangent_laplacian, mass_matrix_sparse, coboundary, complex_residual
from ksf.trace_free import check_trace_free, check_frame_invariance
from ksf.sph        import harmonic

PASS = 0; FAIL = 0

def test(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        print(f"  OK   {name}")
        PASS += 1
    else:
        print(f"  FAIL {name}  {detail}")
        FAIL += 1

print("=== ksf library tests ===\n")

# ── mesh ──────────────────────────────────────────────────────────────────────
print("[ mesh ]")
V, F = icosphere(level=3)
test("vertex count level 3",   len(V) == 642)
test("face count level 3",     len(F) == 1280)
test("vertices on unit sphere", np.allclose(np.linalg.norm(V, axis=1), 1.0, atol=1e-10))
test("mesh size h > 0",        mean_edge_length(V, F) > 0)

# ── DEC complex ───────────────────────────────────────────────────────────────
print("\n[ DEC ]")
d0, d1, _, _ = coboundary(V, F)
test("d1 ∘ d0 = 0",            complex_residual(V, F) == 0.0)
V2, E2, F2 = len(V), d0.shape[0], len(F)
test("Euler characteristic = 2", V2 - E2 + F2 == 2)

L, M_arr = cotangent_laplacian(V, F)
M = mass_matrix_sparse(M_arr)
test("L is square (N×N)",      L.shape == (len(V), len(V)))
test("M diagonal entries > 0", np.all(M.diagonal() > 0))

# ── trace-free projection ─────────────────────────────────────────────────────
print("\n[ trace-free ]")
tr_res  = check_trace_free(V, F)
fr_res  = check_frame_invariance(V, F)
test("trace residual < 1e-15",            tr_res < 1e-15, f"got {tr_res:.2e}")
test("frame-invariance residual < 1e-14", fr_res < 1e-14, f"got {fr_res:.2e}")

# ── spherical harmonic eigenvalues ────────────────────────────────────────────
print("\n[ eigenvalues via Rayleigh quotient ]")
for name, lam_exact in [("Y10", 2.0), ("Y2xy", 6.0), ("Y3xyz", 12.0)]:
    f, lam_e = harmonic(V, name)
    lam_h    = float(f @ L @ f) / float(f @ M @ f)
    rel_err  = abs(lam_h - lam_e) / lam_e
    test(f"{name}: relative error < 1%", rel_err < 0.02, f"got {rel_err*100:.4f}%")

# ── summary ───────────────────────────────────────────────────────────────────
print(f"\n{'='*35}")
print(f"  {PASS} passed,  {FAIL} failed")
if FAIL == 0:
    print("  ALL TESTS PASSED ✓")
else:
    print("  SOME TESTS FAILED")
    sys.exit(1)
