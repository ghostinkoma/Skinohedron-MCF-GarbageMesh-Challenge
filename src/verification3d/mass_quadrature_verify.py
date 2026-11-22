"""
mass_quadrature_verify.py  --  when is the lumped mass "exact"? (scope of 06f §1)

06f reports that the lumped-mass form u^T M u equals the continuous integral of u^2
"at every resolution (0.00%)" for the Taylor-Green mode, and earlier versions read this
as "the mass domain carries NO discretisation error". This script measures the lumped
mass M (fem_laplacian on the Kuhn cube) against exact integrals of u^2 in three cases:

  CHECK A  Whole-period sinusoid on the UNIFORM mesh, u = sin(2 pi x) cos(2 pi y)
           (int u^2 = 1/4, the 06f case). Lumped mass on a uniform grid is the
           trapezoidal rule, which is exact for such trigonometric integrands:
           error at machine precision. This reproduces 06f's 0.00%.
  CHECK B  A smooth NON-trigonometric field on the same uniform mesh,
           u = exp(x) * y  (int u^2 = (e^2 - 1)/6). The trapezoidal-rule exactness
           no longer applies: the error is O(h^2).
  CHECK C  The 06f sinusoid again, but on an IRREGULAR mesh (interior vertices of the
           Kuhn cube jittered by 20% of h, fixed seed). The error is no longer zero;
           it decreases fast with refinement (the random, signed local errors largely
           cancel in the sum), so no clean order is claimed for this case.

Conclusion: the 0.00% of 06f is a property of that TEST (whole-period trigonometric
field on a uniform mesh), not of lumped mass in general. For a generic smooth field
(CHECK B) lumped mass is O(h^2), like the stiffness; on an irregular mesh (CHECK C)
it is not exact either.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from ksf3d.mesh3d_uniform import kuhn_cube
from ksf3d.fem3d import fem_laplacian

NS = (8, 12, 16, 24)


def lumped_error(n, field, exact, jitter=0.0, seed=0):
    V, T = kuhn_cube(n)
    if jitter > 0:
        rng = np.random.default_rng(seed)
        interior = ~((np.isclose(V, 0.0) | np.isclose(V, 1.0)).any(axis=1))
        V = V.copy()
        V[interior] += jitter * (1.0 / n) * rng.uniform(-1, 1, size=(interior.sum(), 3))
    _, M = fem_laplacian(V, T)
    u = field(V)
    return abs(np.sum(M * u * u) - exact) / exact


def order(errs):
    h = 1.0 / np.array(NS, float)
    return float(np.polyfit(np.log(h), np.log(np.array(errs)), 1)[0])


def run_case(tag, field, exact, jitter=0.0):
    errs = [lumped_error(n, field, exact, jitter) for n in NS]
    print("    " + "  ".join(f"n={n}: {e:.2e}" for n, e in zip(NS, errs)))
    return errs


def main():
    print("Lumped mass M (fem_laplacian, Kuhn cube): relative error of u^T M u vs int u^2")

    print("\n[A] whole-period sinusoid sin(2 pi x) cos(2 pi y), uniform mesh  (06f case)")
    fA = lambda V: np.sin(2 * np.pi * V[:, 0]) * np.cos(2 * np.pi * V[:, 1])
    eA = run_case("A", fA, 0.25)
    assert max(eA) < 1e-12, "06f's exactness not reproduced"
    print("    PASS  -- exact to machine precision (trapezoidal rule is exact for this field)")

    print("\n[B] non-trigonometric field exp(x) * y, uniform mesh")
    fB = lambda V: np.exp(V[:, 0]) * V[:, 1]
    eB = run_case("B", fB, (np.e**2 - 1.0) / 6.0)
    oB = order(eB)
    print(f"    order ~ h^{oB:.2f}")
    assert min(eB) > 1e-8 and 1.8 < oB < 2.2, "expected a non-zero O(h^2) error"
    print("    PASS  -- not exact: O(h^2)")

    print("\n[C] the same sinusoid on an irregular mesh (interior vertices jittered 20% of h)")
    eC = run_case("C", fA, 0.25, jitter=0.2)
    oC = order(eC)
    print(f"    least-squares fit ~ h^{oC:.2f}  (random signed local errors partly cancel;")
    print("    not a clean order)")
    assert min(eC) > 1e-8 and eC[-1] < eC[0], "expected a non-zero error that decreases"
    print("    PASS  -- not exact on an irregular mesh; decreases with refinement")

    print("\nConclusion: the 0.00% lumped-mass error of 06f holds for a whole-period")
    print("trigonometric field on a uniform mesh (trapezoidal-rule exactness). A generic")
    print("smooth field gives O(h^2) [B], and an irregular mesh is not exact either [C]:")
    print("'the mass domain has no discretisation error' is a property of that test only.")


if __name__ == "__main__":
    main()
