"""
export_fluid_shapes_viewer.py  --  fluid-on-shapes viewer data builder.

Generates advection animations on the three solids using the SAME operators and the
SAME natural flows as fluid_shapes_verify.py:
  cube  : uniform advection (x-periodic)        -> blob translates and wraps
  torus : uniform around the ring (x-periodic)  -> blob travels around the ring
  ball  : solid-body rotation u = omega x r      -> blob rotates (degraded mesh)

Every frame is a verified-solver snapshot; the viewer only plays them back. The ball
is the geometric (degraded) mesh: its accuracy is FE-degraded (see note), though the
incompressibility/conservation it relies on stay exact (fluid_shapes_verify.py).

Emits viewer/data_fluid_shapes.js.
"""
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from scipy.sparse import coo_matrix
from ksf3d.mesh3d_uniform import kuhn_cube
from ksf3d.fem3d import fem_laplacian

FACE = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
NFRAMES = 30


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


def build_C(ufield, tets, n):
    rows, cols, data = [], [], []
    for (ix, vol, G) in tets:
        ue = ufield[ix].mean(0)
        ug = ue @ G
        for a in range(4):
            for b in range(4):
                rows.append(ix[a]); cols.append(ix[b]); data.append(vol/4.0*ug[b])
    return coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()


def faces_edges(T):
    cnt = {}
    for tet in T:
        for fc in FACE:
            k = tuple(sorted(int(tet[i]) for i in fc))
            cnt[k] = cnt.get(k, 0)+1
    bnd = [list(k) for k, c in cnt.items() if c == 1]
    es = set()
    for tet in T:
        for a in range(4):
            for b in range(a+1, 4):
                es.add(tuple(sorted((int(tet[a]), int(tet[b])))))
    return bnd, [list(e) for e in es]


def norm_frames(frames):
    a = max(float(np.max(np.abs(f))) for f in frames) or 1.0
    return [[round(float(x/a), 4) for x in f] for f in frames]


def arrow_idx(V, mid_axis=2, mid=0.0):
    sel = np.where(np.abs(V[:, mid_axis]-mid) < 1e-6)[0]
    return sel[::max(1, len(sel)//40)]


# --------------------------------------------------------------------------- #
def scene_cube(n=8, kappa=0.004):
    V, T = kuhn_cube(n); V = V-0.5
    K, M = fem_laplacian(V, T); tets = per_tet(V, T); nV = len(V)
    inv, nn = remap_x(V)
    uf = np.zeros((nV, 3)); uf[:, 0] = 1.0
    C = build_C(uf, tets, nV)
    Cp = fold(C, inv, nn); Kp = fold(K, inv, nn)
    Mp = np.zeros(nn)
    for i, m in enumerate(M):
        Mp[inv[i]] += m
    Mi = 1.0/Mp; Cs = (Cp - Cp.T)*0.5
    # periodic representative coords (for the initial blob)
    Vrep = np.zeros((nn, 3))
    for i, p in enumerate(V):
        pp = p.copy()
        if np.isclose(p[0], 0.5):
            pp[0] = -0.5
        Vrep[inv[i]] = pp
    phi = np.exp(-((Vrep[:, 0]+0.2)**2+(Vrep[:, 1])**2)/(2*0.10**2))
    fr_p = integrate(phi, Cs, Kp, Mi, kappa, total=3000)
    # display on the FULL cube (clean cube look): scatter periodic -> full vertices
    frames = [[row[inv[i]] for i in range(nV)] for row in fr_p]
    bnd, wed = faces_edges(T)
    ai = arrow_idx(V, 2, 0.0)
    vec = np.zeros((len(ai), 3)); vec[:, 0] = 1.0
    return pack(V, bnd, wed, frames, V[ai], vec, note="flat cube · uniform flow (FE-accurate)")


def scene_torus(n=8, R=1.0, w=0.62, kappa=0.004):
    V, T = kuhn_cube(n); V = V-0.5
    K, M = fem_laplacian(V, T); tets = per_tet(V, T); nV = len(V)
    inv, nn = remap_x(V)
    uf = np.zeros((nV, 3)); uf[:, 0] = 1.0          # advect around major angle
    C = build_C(uf, tets, nV)
    Cp = fold(C, inv, nn); Kp = fold(K, inv, nn)
    Mp = np.zeros(nn)
    for i, m in enumerate(M):
        Mp[inv[i]] += m
    Mi = 1.0/Mp; Cs = (Cp - Cp.T)*0.5
    # ring geometry for each representative
    Vr = np.zeros((nn, 3))
    foldc = np.zeros((nn, 3))
    for old in range(nV):
        x, y, z = V[old]
        th = 2*np.pi*(x+0.5); rad = R + y*w
        Vr[inv[old]] = [rad*np.cos(th), rad*np.sin(th), z*w]
        pp = V[old].copy()
        if np.isclose(x, 0.5):
            pp[0] = -0.5
        foldc[inv[old]] = pp
    phi = np.exp(-((foldc[:, 0]+0.2)**2)/(2*0.09**2))   # blob localized in major angle
    frames = integrate(phi, Cs, Kp, Mi, kappa, total=3200)
    bnd, wed = faces_edges(fold_T(T, inv))
    # arrows tangent to the ring on the outer rim (y mid, z mid)
    ai = np.where((np.abs(foldc[:, 2]) < 1e-6))[0]
    ai = ai[::max(1, len(ai)//44)]
    vec = np.zeros((len(ai), 3))
    for k, idx in enumerate(ai):
        th = np.arctan2(Vr[idx, 1], Vr[idx, 0])
        vec[k] = [-np.sin(th), np.cos(th), 0.0]          # tangent to ring
    return pack(Vr, bnd, wed, frames, Vr[ai], vec, note="flat torus · ring flow (FE-accurate, wraps)")


def scene_ball(n=8, kappa=0.0):
    V0, T = kuhn_cube(n); V0 = V0-0.5
    Vb = np.zeros_like(V0)
    for i, p in enumerate(V0):
        m = np.max(np.abs(p)); Vb[i] = p*(m/np.linalg.norm(p)) if m > 1e-12 else p
    K, M = fem_laplacian(Vb, T); tets = per_tet(Vb, T); nV = len(Vb)
    omega = np.array([0.0, 0.0, 1.2])
    uf = np.cross(np.tile(omega, (nV, 1)), Vb)
    C = build_C(uf, tets, nV); Cs = (C - C.T)*0.5
    Mi = 1.0/M
    phi = np.exp(-(((Vb[:, 0]-0.28)**2+Vb[:, 1]**2+(Vb[:, 2])**2))/(2*0.14**2))
    frames = integrate(phi, Cs, None, Mi, 0.0, total=2200)
    bnd, wed = faces_edges(T)
    ai = arrow_idx(Vb, 2, 0.0)
    vec = np.cross(np.tile(omega, (len(ai), 1)), Vb[ai])
    flips = count_flips(K)
    return pack(Vb, bnd, wed, frames, Vb[ai], vec,
                note=f"geometric ball · rotation · DEGRADED ({flips} sign-flips, FE accuracy reduced)")


# --- helpers ---------------------------------------------------------------
def remap_x(V):
    n = len(V); remap = np.arange(n)
    for ip in np.where(np.isclose(V[:, 0], 0.5))[0]:
        q = V[ip].copy(); q[0] = -0.5
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


def fold(A, inv, nn):
    Ac = A.tocoo()
    return coo_matrix((Ac.data, (inv[Ac.row], inv[Ac.col])), shape=(nn, nn)).tocsr()


def fold_T(T, inv):
    return np.array([[inv[int(i)] for i in tet] for tet in T])


def ring_or_fold_coords(V, inv, nn, kind="fold"):
    Vr = np.zeros((nn, 3))
    for i, p in enumerate(V):
        pp = p.copy()
        if np.isclose(p[0], 0.5):
            pp[0] = -0.5
        Vr[inv[i]] = pp
    return Vr


def integrate(phi, Cs, Kp, Mi, kappa, total):
    dt = 4e-4
    pick = set(np.linspace(0, total-1, NFRAMES).astype(int).tolist())
    frames = []
    for step in range(total):
        if step in pick:
            frames.append(phi.copy())
        adv = Cs @ phi
        phi = phi - dt*Mi*adv
        if Kp is not None and kappa > 0:
            phi = phi + dt*kappa*Mi*(-(Kp @ phi))
    return norm_frames(frames)


def count_flips(K):
    Kc = K.tocoo()
    return int(((Kc.row != Kc.col) & (Kc.data > 1e-12)).sum())


def pack(V, bnd, wed, frames, apos, avec, note=""):
    av = float(np.max(np.abs(avec)) or 1.0)
    return {"nV": int(len(V)),
            "verts": [round(float(x), 5) for x in V.flatten()],
            "btris": [int(x) for t in bnd for x in t],
            "wedges": [int(x) for e in wed for x in e],
            "frames": frames,
            "note": note,
            "arrows": {"pos": [round(float(x), 4) for x in np.asarray(apos).flatten()],
                       "vec": [round(float(x/av), 4) for x in np.asarray(avec).flatten()]}}


def main():
    data = {"frames": NFRAMES, "shapes": {}}
    print("building cube ...");  data["shapes"]["cube"] = scene_cube()
    print("building torus ..."); data["shapes"]["torus"] = scene_torus()
    print("building ball ...");  data["shapes"]["ball"] = scene_ball()
    out = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "viewer", "data_fluid_shapes.js"))
    with open(out, "w") as f:
        f.write("// Auto-generated by src/verification3d/export_fluid_shapes_viewer.py\n")
        f.write("// Advection on cube/torus/sphere, verified operators + natural flows. Playback only.\n")
        f.write("window.KSF_FLUID_SHAPES = ")
        json.dump(data, f, separators=(",", ":"))
        f.write(";\n")
    print("wrote", out, "(%d KB)" % (os.path.getsize(out)//1024))
    for k, s in data["shapes"].items():
        print(f"  {k}: nV={s['nV']} frames={len(s['frames'])} arrows={len(s['arrows']['pos'])//3}  [{s['note']}]")


if __name__ == "__main__":
    main()
