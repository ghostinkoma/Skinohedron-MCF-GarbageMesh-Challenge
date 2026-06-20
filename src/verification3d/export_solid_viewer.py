"""
export_solid_viewer.py  --  grand finale data builder.

Three solid models (cube, ball, torus) carrying the same operator L = M^{-1}K,
each driven by the verified physics (heat diffusion, scalar wave, liquid/gas
pressure). The point of the viewer this feeds: the SAME framework simulates the
SAME physics on different shapes -- the deformation result of V2.5 made visible.

Emits viewer/data_solid.js with, per shape: mesh (verts, boundary tris, all tris,
wireframe edges) and presets {heat, wave, liquid, gas}, each a verified field
normalised for display. Custom parameters are handled at run time by the server
(api/solve_solid.py); this file is the instant, offline preset set.
"""
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from scipy.sparse import diags, coo_matrix
from scipy.sparse.linalg import spsolve, eigsh
from ksf3d.mesh3d_uniform import kuhn_cube
from ksf3d.fem3d import fem_laplacian

FACE = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]


# --- mesh builders ----------------------------------------------------------
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


def faces_edges(T, nV):
    cnt = {}
    for tet in T:
        for fc in FACE:
            k = tuple(sorted(int(tet[i]) for i in fc))
            cnt[k] = cnt.get(k, 0) + 1
    bnd = [list(k) for k, c in cnt.items() if c == 1]
    allf = [list(k) for k in cnt.keys()]
    es = set()
    for tet in T:
        for a in range(4):
            for b in range(a+1, 4):
                es.add(tuple(sorted((int(tet[a]), int(tet[b])))))
    return bnd, allf, [list(e) for e in es]


def build_cube(n=8):
    V, T = kuhn_cube(n); V = V - 0.5
    return V, T


def build_ball(n=8):
    V, T = kuhn_cube(n); V = V - 0.5
    out = np.zeros_like(V)
    for i, p in enumerate(V):
        m = np.max(np.abs(p))
        out[i] = p*(m/np.linalg.norm(p)) if m > 1e-12 else p   # cube face -> sphere
    return out, T


def build_torus(n=8, R=1.0, w=0.62):
    """Solid square-section ring: identify x=+0.5 with x=-0.5 (periodic major
    angle), place (y,z) as the cross-section. Wave wraps around the ring."""
    V, T = kuhn_cube(n); V = V - 0.5
    L = 1.0
    # identify x = +0.5 -> x = -0.5 (same y,z)
    remap = np.arange(len(V))
    plus = np.where(np.isclose(V[:, 0], 0.5))[0]
    for ip in plus:
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
    nnew = len(uniq)
    Tn = np.array([[inv[int(i)] for i in tet] for tet in T])
    # ring coordinates for each representative
    Vn = np.zeros((nnew, 3))
    for old in range(len(V)):
        x, y, z = V[old]
        th = 2*np.pi*(x + 0.5)                  # x in [-0.5,0.5] -> angle
        rad = R + y*w
        Vn[inv[old]] = [rad*np.cos(th), rad*np.sin(th), z*w]
    return Vn, Tn


# --- physics (verified solvers) ---------------------------------------------
def heat_snapshot(V, K, M, steps=60):
    Minv = 1.0/M
    dt = 0.3/float(np.max(K.diagonal()*Minv))
    c = int(np.argmin(((V-V.mean(0))**2).sum(1)))
    T0 = np.zeros(len(V)); T0[c] = 1.0
    Tf = T0.copy()
    for _ in range(steps):
        Tf = Tf - dt*Minv*(K@Tf)
    return Tf, False


def wave_snapshot(V, K, M, steps=90):
    Minv = 1.0/M
    lam, vec = eigsh(K, k=3, M=diags(M), sigma=1e-8, which="LM")
    order = np.argsort(lam)
    mode = vec[:, order[1]]
    lam_max = 2.2*float(np.max(K.diagonal()*Minv))
    dt = 0.4*np.sqrt(4.0/lam_max)
    p0 = mode.copy(); p1 = mode.copy()
    for _ in range(steps):
        p2 = 2*p1 - p0 - dt*dt*Minv*(K@p1)
        p0, p1 = p1, p2
    return p1, True


def liquid_pressure(V, K, tets):
    g = np.array([0.0, 0.0, -1.0])
    n = len(V)
    b = np.zeros(n)
    for (ix, vol, G) in tets:
        cpt = vol*(G.T@g)
        for a in range(4):
            b[ix[a]] += cpt[a]
    keep = np.ones(n, bool); keep[0] = False
    p = np.zeros(n)
    p[keep] = spsolve(K[keep][:, keep].tocsc(), b[keep])
    return p, False


def gas_atmosphere(V, c2=0.18):
    z = V[:, 2] - V[:, 2].min()
    return np.exp(-1.0*z/c2), False


def norm(p, signed):
    a = float(np.max(np.abs(p)))
    if a <= 0:
        return [0.0]*len(p)
    return [round(float(x/a), 5) for x in p]


def shape_block(V, T):
    K, M, tets = operators(V, T)
    bnd, allf, wed = faces_edges(T, len(V))
    presets = {}
    for name, (p, signed) in {
        "heat":   heat_snapshot(V, K, M),
        "wave":   wave_snapshot(V, K, M),
        "liquid": liquid_pressure(V, K, tets),
        "gas":    gas_atmosphere(V),
    }.items():
        presets[name] = {"signed": signed, "p": norm(p, signed)}
    return {
        "nV": int(len(V)),
        "verts": [round(float(x), 5) for x in V.flatten()],
        "btris": [int(x) for t in bnd for x in t],
        "atris": [int(x) for t in allf for x in t],
        "wedges": [int(x) for e in wed for x in e],
        "presets": presets,
    }


def main():
    shapes = {
        "cube":  build_cube(8),
        "ball":  build_ball(8),
        "torus": build_torus(8),
    }
    data = {"shapes": {}}
    for name, (V, T) in shapes.items():
        print("building", name, "...")
        data["shapes"][name] = shape_block(V, T)
        b = data["shapes"][name]
        print(f"  nV={b['nV']} btris={len(b['btris'])//3} wedges={len(b['wedges'])//2} presets={list(b['presets'])}")
    out = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "viewer", "data_solid.js"))
    with open(out, "w") as f:
        f.write("// Auto-generated by src/verification3d/export_solid_viewer.py\n")
        f.write("// Three solids (cube/ball/torus) x physics (heat/wave/liquid/gas),\n")
        f.write("// all from the verified L = M^-1 K. Same framework, different shapes.\n")
        f.write("window.KSF_SOLID = ")
        json.dump(data, f, separators=(",", ":"))
        f.write(";\n")
    print("wrote", out, "(%d bytes)" % os.path.getsize(out))


if __name__ == "__main__":
    main()
