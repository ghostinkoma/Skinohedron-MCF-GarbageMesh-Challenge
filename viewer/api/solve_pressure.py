#!/usr/bin/env python3
"""
solve_pressure.py  --  CLI pressure solver behind the PHP API.

Reads a JSON request on stdin, returns a JSON response on stdout. Uses the SAME
operator and solver as src/verification3d/pressure_field_verify.py, so an API
result is exactly a verified value.

Request  (stdin JSON):
  { "mode": "hydrostatic" | "atmosphere" | "acoustic",
    "rho": 1.0, "g": [0,0,-1.0], "c2": 0.15, "steps": 120 }

Response (stdout JSON):
  { "ok": true, "mode": "...", "nV": 729, "signed": false, "p": [ ... ] }
  values normalised to max|p| = 1 for display.
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve, eigsh
from ksf3d.mesh3d_uniform import kuhn_cube
from ksf3d.fem3d import fem_laplacian

FACE = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]


def operators(V, T):
    tets = []
    for tet in T:
        ix = [int(i) for i in tet]
        P = V[ix]
        vol = abs(np.linalg.det(np.array([P[1]-P[0], P[2]-P[0], P[3]-P[0]]))) / 6.0
        if vol <= 1e-15:
            continue
        C = np.linalg.inv(np.column_stack([np.ones(4), P]))
        tets.append((ix, vol, C[1:4, :]))
    K, M = fem_laplacian(V, T)
    return K, M, tets


def div_source(f_tet, tets, n):
    b = np.zeros(n)
    for t, (ix, vol, G) in enumerate(tets):
        c = vol * (G.T @ f_tet[t])
        for a in range(4):
            b[ix[a]] += c[a]
    return b


def solve_neumann(K, b, pin=0):
    n = K.shape[0]
    keep = np.ones(n, bool); keep[pin] = False
    p = np.zeros(n)
    p[keep] = spsolve(K[keep][:, keep].tocsc(), b[keep])
    return p


def compute(req):
    n = int(req.get("n", 8))
    mode = req.get("mode", "hydrostatic")
    V, T = kuhn_cube(n); V = V - 0.5
    K, M, tets = operators(V, T)
    nV = len(V)

    if mode == "hydrostatic":
        rho = float(req.get("rho", 1.0))
        g = np.array(req.get("g", [0, 0, -1.0]), float)
        f = np.tile(rho * g, (len(tets), 1))
        p = solve_neumann(K, div_source(f, tets, nV))
        signed = False
    elif mode == "atmosphere":
        p0 = float(req.get("rho", 1.0))
        g = float(np.linalg.norm(req.get("g", [0, 0, -1.0])))
        c2 = max(float(req.get("c2", 0.15)), 1e-6)
        z = V[:, 2] - V[:, 2].min()
        p = p0 * np.exp(-g * z / c2)
        signed = False
    elif mode == "acoustic":
        c = float(np.sqrt(max(float(req.get("c2", 1.0)), 1e-9)))
        steps = int(req.get("steps", 120))
        Minv = 1.0 / M
        lam, vec = eigsh(K, k=2, M=diags(M), sigma=1e-8, which="LM")
        mode_v = vec[:, np.argsort(lam)[1]]
        lam_max = 2.2 * float(np.max(K.diagonal() * Minv)) * c * c
        dt = 0.4 * np.sqrt(4.0 / lam_max)
        p0v = mode_v.copy(); p1 = mode_v.copy()
        for _ in range(steps):
            p2 = 2*p1 - p0v - dt*dt*Minv*(c*c*(K @ p1))
            p0v, p1 = p1, p2
        p = p1; signed = True
    else:
        raise ValueError("unknown mode: %s" % mode)

    a = float(np.max(np.abs(p)))
    pl = (p / a) if a > 0 else p
    return {"ok": True, "mode": mode, "nV": int(nV), "signed": bool(signed),
            "p": [round(float(x), 5) for x in pl]}


def main():
    try:
        req = json.load(sys.stdin) if not sys.stdin.isatty() else {}
    except Exception:
        req = {}
    try:
        print(json.dumps(compute(req)))
    except Exception as e:
        print(json.dumps({"ok": False, "error": str(e)}))


if __name__ == "__main__":
    main()
