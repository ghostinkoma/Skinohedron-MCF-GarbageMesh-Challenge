#!/usr/bin/env python3
"""
solve_fluid.py  --  CLI solver behind the fluid-viewer API (V3 Stage A/B).

Reads JSON on stdin, returns a time-series scene on stdout, using the SAME operators
as fluid_stokes_verify.py and fluid_advection_verify.py. Lets the viewer compute
custom parameters; every frame is a verified solver output.

Request:  { "scene":"couette|poiseuille|advection",
            "U":1.0, "nu":1.0, "G":-1.0, "kappa":0.004, "ux":1.0, "frames":36 }
Response: { "ok":true, "scene":..., "kind":..., "frames":[[...],...],
            "arrows":{"pos":[...],"vec":[[...],...]}, ... }
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
import numpy as np
from scipy.sparse import coo_matrix
from ksf3d.mesh3d_uniform import kuhn_cube
from ksf3d.fem3d import fem_laplacian

N = 10


def per_tet(V, T):
    tets = []
    for tet in T:
        ix = [int(i) for i in tet]
        P = V[ix]
        vol = abs(np.linalg.det(np.array([P[1]-P[0], P[2]-P[0], P[3]-P[0]]))) / 6.0
        if vol <= 1e-15:
            continue
        C = np.linalg.inv(np.column_stack([np.ones(4), P]))
        tets.append((ix, vol, C[1:4, :]))
    return tets


def build_C(u0, tets, n):
    rows, cols, data = [], [], []
    for (ix, vol, G) in tets:
        ug = u0 @ G
        for a in range(4):
            for b in range(4):
                rows.append(ix[a]); cols.append(ix[b]); data.append(vol/4.0*ug[b])
    return coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()


def remap_x(V):
    n = len(V); remap = np.arange(n)
    for ip in np.where(np.isclose(V[:, 0], 1.0))[0]:
        q = V[ip].copy(); q[0] = 0.0
        j = np.where(np.all(np.isclose(V, q), axis=1))[0]
        if len(j):
            remap[ip] = j[0]
    for i in range(n):
        r = i
        while remap[r] != r:
            r = remap[r]
        remap[i] = r
    _, inv = np.unique(remap, return_inverse=True)
    return inv, int(inv.max()+1)


def arrows_mid(V):
    return np.where((np.abs(V[:, 2]-0.5) < 1e-6) &
                    (np.abs((V[:, 0]/0.2)-np.round(V[:, 0]/0.2)) < 1e-6) &
                    (np.abs((V[:, 1]/0.2)-np.round(V[:, 1]/0.2)) < 1e-6))[0]


def norm_frames(frames):
    a = max(float(np.max(np.abs(f))) for f in frames) or 1.0
    return [[round(float(x/a), 4) for x in f] for f in frames]


def compute(req):
    scene = req.get("scene", "couette")
    nf = int(req.get("frames", 36))
    V, T = kuhn_cube(N)
    K, M = fem_laplacian(V, T)
    Minv = 1.0/M
    n = len(V); y = V[:, 1]

    if scene in ("couette", "poiseuille"):
        nu = float(req.get("nu", 1.0))
        wall = np.isclose(y, 0.0) | np.isclose(y, 1.0)
        free = ~wall
        dt = 0.2/float(np.max(K.diagonal()*Minv))
        ap = arrows_mid(V)
        if scene == "couette":
            U = float(req.get("U", 1.0))
            u = np.zeros(n); u[np.isclose(y, 1.0)] = U
            total = 6000
        else:
            G = float(req.get("G", -1.0))
            body = M*(-(G/nu)*np.ones(n))
            u = np.zeros(n); total = 9000
        pick = set(np.linspace(0, total-1, nf).astype(int).tolist())
        frames, avecs = [], []
        for step in range(total):
            if step in pick:
                frames.append(np.abs(u).copy())
                vv = np.zeros((len(ap), 3)); vv[:, 0] = u[ap]; avecs.append(vv.flatten())
            if scene == "couette":
                u[free] = u[free] - dt*nu*Minv[free]*(K@u)[free]
            else:
                u[free] = u[free] + dt*Minv[free]*(-nu*(K@u)[free] + nu*body[free])
        av = float(max(np.max(np.abs(a)) for a in avecs) or 1.0)
        return {"ok": True, "scene": scene, "kind": "velocity", "nV": int(n),
                "frames": norm_frames(frames),
                "arrows": {"pos": [round(float(x), 4) for x in V[ap].flatten()],
                           "vec": [[round(float(x/av), 4) for x in a] for a in avecs]}}

    elif scene == "advection":
        kappa = max(float(req.get("kappa", 0.004)), 0.0)
        ux = float(req.get("ux", 1.0))
        tets = per_tet(V, T)
        inv, nn = remap_x(V)
        C = build_C(np.array([ux, 0.0, 0.0]), tets, n)
        Cc = C.tocoo(); Cp = coo_matrix((Cc.data, (inv[Cc.row], inv[Cc.col])), shape=(nn, nn)).tocsr()
        Kc = K.tocoo(); Kp = coo_matrix((Kc.data, (inv[Kc.row], inv[Kc.col])), shape=(nn, nn)).tocsr()
        Cs = (Cp - Cp.T)*0.5
        Mp = np.zeros(nn)
        for i, m in enumerate(M):
            Mp[inv[i]] += m
        Mi = 1.0/Mp
        Vr = np.zeros((nn, 3))
        for i, p in enumerate(V):
            pp = p.copy()
            if np.isclose(p[0], 1.0):
                pp[0] = 0.0                   # fold ONLY x (periodic in x)
            Vr[inv[i]] = pp
        phi = np.exp(-((Vr[:, 0]-0.3)**2+(Vr[:, 1]-0.5)**2)/(2*0.07**2))
        dt = 4e-4; total = 4500
        pick = set(np.linspace(0, total-1, nf).astype(int).tolist())
        apr = arrows_mid(Vr)
        if len(apr) == 0:
            apr = np.where(np.abs(Vr[:, 2]-0.5) < 1e-6)[0][::6]
        frames = []
        for step in range(total):
            if step in pick:
                frames.append(phi.copy())
            phi = phi - dt*Mi*(Cs@phi) + dt*kappa*Mi*(-(Kp@phi))
        vv = np.tile([ux, 0.0, 0.0], (len(apr), 1))
        av = float(np.max(np.abs(vv)) or 1.0)
        return {"ok": True, "scene": scene, "kind": "scalar", "nV": int(nn),
                "verts": [round(float(x), 5) for x in Vr.flatten()],
                "frames": norm_frames(frames),
                "arrows": {"pos": [round(float(x), 4) for x in Vr[apr].flatten()],
                           "vec": [[round(float(x/av), 4) for x in vv.flatten()]]}}
    else:
        raise ValueError("unknown scene: %s" % scene)


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
