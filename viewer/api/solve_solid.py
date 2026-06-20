#!/usr/bin/env python3
"""
solve_solid.py  --  CLI solver behind the solid-viewer API.

Reads a JSON request on stdin, returns a JSON field on stdout, using the SAME
operator and solvers as the verification scripts. Lets the viewer compute custom
parameters (time, gravity, c^2, mode) on any of the three shapes.

Request:  { "shape":"cube|ball|torus", "physics":"heat|wave|liquid|gas",
            "steps":60, "g":1.0, "c2":0.18 }
Response: { "ok":true, "shape":..., "physics":..., "nV":N, "signed":bool, "p":[...] }
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve, eigsh
from ksf3d.mesh3d_uniform import kuhn_cube
from ksf3d.fem3d import fem_laplacian


def build_cube(n=8):
    V, T = kuhn_cube(n); return V-0.5, T


def build_ball(n=8):
    V, T = kuhn_cube(n); V = V-0.5
    out = np.zeros_like(V)
    for i, p in enumerate(V):
        m = np.max(np.abs(p))
        out[i] = p*(m/np.linalg.norm(p)) if m > 1e-12 else p
    return out, T


def build_torus(n=8, R=1.0, w=0.62):
    V, T = kuhn_cube(n); V = V-0.5
    remap = np.arange(len(V))
    for ip in np.where(np.isclose(V[:, 0], 0.5))[0]:
        q = V[ip].copy(); q[0] = -0.5
        j = np.where(np.all(np.isclose(V, q), axis=1))[0]
        if len(j):
            remap[ip] = j[0]
    for i in range(len(V)):
        r = i
        while remap[r] != r:
            r = remap[r]
        remap[i] = r
    uniq, inv = np.unique(remap, return_inverse=True)
    Tn = np.array([[inv[int(i)] for i in tet] for tet in T])
    Vn = np.zeros((len(uniq), 3))
    for old in range(len(V)):
        x, y, z = V[old]
        th = 2*np.pi*(x+0.5); rad = R + y*w
        Vn[inv[old]] = [rad*np.cos(th), rad*np.sin(th), z*w]
    return Vn, Tn


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


def compute(req):
    shape = req.get("shape", "cube")
    physics = req.get("physics", "heat")
    n = int(req.get("n", 8))
    V, T = {"cube": build_cube, "ball": build_ball, "torus": build_torus}[shape](n)
    K, M, tets = operators(V, T)
    Minv = 1.0/M
    nV = len(V)

    if physics == "heat":
        steps = int(req.get("steps", 60))
        dt = 0.3/float(np.max(K.diagonal()*Minv))
        c = int(np.argmin(((V-V.mean(0))**2).sum(1)))
        p = np.zeros(nV); p[c] = 1.0
        for _ in range(steps):
            p = p - dt*Minv*(K@p)
        signed = False
    elif physics == "wave":
        steps = int(req.get("steps", 90))
        lam, vec = eigsh(K, k=3, M=diags(M), sigma=1e-8, which="LM")
        mode = vec[:, np.argsort(lam)[1]]
        lam_max = 2.2*float(np.max(K.diagonal()*Minv))
        dt = 0.4*np.sqrt(4.0/lam_max)
        p0 = mode.copy(); p1 = mode.copy()
        for _ in range(steps):
            p0, p1 = p1, 2*p1 - p0 - dt*dt*Minv*(K@p1)
        p = p1; signed = True
    elif physics == "liquid":
        g = float(req.get("g", 1.0))
        gv = np.array([0.0, 0.0, -g])
        b = np.zeros(nV)
        for (ix, vol, G) in tets:
            cpt = vol*(G.T@gv)
            for a in range(4):
                b[ix[a]] += cpt[a]
        keep = np.ones(nV, bool); keep[0] = False
        p = np.zeros(nV)
        p[keep] = spsolve(K[keep][:, keep].tocsc(), b[keep])
        signed = False
    elif physics == "gas":
        c2 = max(float(req.get("c2", 0.18)), 1e-4)
        g = float(req.get("g", 1.0))
        z = V[:, 2] - V[:, 2].min()
        p = np.exp(-g*z/c2); signed = False
    else:
        raise ValueError("unknown physics: %s" % physics)

    a = float(np.max(np.abs(p)))
    pl = (p/a) if a > 0 else p
    return {"ok": True, "shape": shape, "physics": physics, "nV": int(nV),
            "signed": bool(signed), "p": [round(float(x), 5) for x in pl]}


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
