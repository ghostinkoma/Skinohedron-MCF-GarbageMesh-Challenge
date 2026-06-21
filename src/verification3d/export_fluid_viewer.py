"""
export_fluid_viewer.py  --  V3 fluid viewer data builder.

Precomputes TIME-SERIES snapshots for viewer_fluid.html using the SAME operators as
fluid_stokes_verify.py (Stage A) and fluid_advection_verify.py (Stage B). The viewer
only *plays back* these frames, so every frame the user sees is a verified solver
output -- no in-browser physics.

Scenes:
  couette     : Stokes startup, walls drive flow -> linear profile u = U y.
                heatmap = speed |u|, arrows = velocity, animated startup.
  poiseuille  : Stokes startup, pressure gradient -> parabolic profile.
  advection   : Stage B Gaussian scalar carried by a uniform flow and diffusing.
                heatmap = scalar phi (moves + spreads), arrows = the steady velocity.

Emits viewer/data_fluid.js.
"""
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from scipy.sparse import coo_matrix
from ksf3d.mesh3d_uniform import kuhn_cube
from ksf3d.fem3d import fem_laplacian

FACE = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
NFRAMES = 36


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


def remap_periodic_x(V):
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


def arrow_points(V, every=2):
    """Pick a coarse set of vertices near the z-mid slice for arrows."""
    zmid = 0.5
    near = np.where(np.abs(V[:, 2]-zmid) < 1e-6)[0]
    # subsample on a lattice in x,y
    pts = [i for i in near if (round(V[i, 0]/0.25*1)*1) is not None]
    # thin: keep those whose x,y are multiples of ~0.2
    sel = [i for i in near if (abs((V[i, 0]/0.2) - round(V[i, 0]/0.2)) < 1e-6
                               and abs((V[i, 1]/0.2) - round(V[i, 1]/0.2)) < 1e-6)]
    return np.array(sel if sel else near[::6])


def norm_frames(frames):
    a = max(float(np.max(np.abs(f))) for f in frames) or 1.0
    return [[round(float(x/a), 4) for x in f] for f in frames], a


# --------------------------------------------------------------------------- #
def scene_couette(V, T, K, M, U=1.0, nu=1.0):
    n = len(V); y = V[:, 1]
    wall = np.isclose(y, 0.0) | np.isclose(y, 1.0)
    free = ~wall
    Minv = 1.0/M
    bcv = np.where(np.isclose(y, 1.0), U, 0.0)
    u = np.zeros(n); u[wall] = bcv[wall]
    dt = 0.2/float(np.max(K.diagonal()*Minv))
    total = 6000
    pick = np.linspace(0, total-1, NFRAMES).astype(int)
    ap = arrow_points(V)
    frames, avecs = [], []
    s = 0
    for step in range(total):
        if step in pick:
            speed = np.abs(u)                  # |u_x|
            frames.append(speed.copy())
            vv = np.zeros((len(ap), 3)); vv[:, 0] = u[ap]
            avecs.append(vv.flatten())
        u[free] = u[free] - dt*nu*Minv[free]*(K@u)[free]
    fr, amax = norm_frames(frames)
    av = float(max(np.max(np.abs(a)) for a in avecs) or 1.0)
    return {"kind": "velocity", "signed": False, "frames": fr,
            "arrows": {"pos": [round(float(x), 4) for x in V[ap].flatten()],
                       "vec": [[round(float(x/av), 4) for x in a] for a in avecs]}}


def scene_poiseuille(V, T, K, M, nu=1.0, G=-1.0):
    n = len(V); y = V[:, 1]
    wall = np.isclose(y, 0.0) | np.isclose(y, 1.0)
    free = ~wall
    Minv = 1.0/M
    u = np.zeros(n)
    body = M*(-(G/nu)*np.ones(n))
    dt = 0.2/float(np.max(K.diagonal()*Minv))
    total = 9000
    pick = np.linspace(0, total-1, NFRAMES).astype(int)
    ap = arrow_points(V)
    frames, avecs = [], []
    for step in range(total):
        if step in pick:
            frames.append(np.abs(u).copy())
            vv = np.zeros((len(ap), 3)); vv[:, 0] = u[ap]
            avecs.append(vv.flatten())
        u[free] = u[free] + dt*Minv[free]*(-nu*(K@u)[free] + nu*body[free])
    fr, amax = norm_frames(frames)
    av = float(max(np.max(np.abs(a)) for a in avecs) or 1.0)
    return {"kind": "velocity", "signed": False, "frames": fr,
            "arrows": {"pos": [round(float(x), 4) for x in V[ap].flatten()],
                       "vec": [[round(float(x/av), 4) for x in a] for a in avecs]}}


def scene_advection(V, T, K, M, tets, u0=(1.0, 0.0, 0.0), kappa=0.004):
    n = len(V)
    inv, nn = remap_periodic_x(V)
    C = build_C(np.array(u0), tets, n)
    Cc = C.tocoo(); Cp = coo_matrix((Cc.data, (inv[Cc.row], inv[Cc.col])), shape=(nn, nn)).tocsr()
    Kc = K.tocoo(); Kp = coo_matrix((Kc.data, (inv[Kc.row], inv[Kc.col])), shape=(nn, nn)).tocsr()
    Cs = (Cp - Cp.T)*0.5
    Mp = np.zeros(nn)
    for i, m in enumerate(M):
        Mp[inv[i]] += m
    Minv = 1.0/Mp
    Vr = np.zeros((nn, 3))
    for i, p in enumerate(V):
        pp = p.copy()
        if np.isclose(p[0], 1.0):
            pp[0] = 0.0                       # fold ONLY x (periodic in x)
        Vr[inv[i]] = pp
    phi = np.exp(-((Vr[:, 0]-0.3)**2+(Vr[:, 1]-0.5)**2)/(2*0.07**2))
    dt = 4e-4
    total = 4500
    pick = np.linspace(0, total-1, NFRAMES).astype(int)
    # arrow points on the periodic mesh near z-mid
    apr = np.where((np.abs(Vr[:, 2]-0.5) < 1e-6) &
                   (np.abs((Vr[:, 0]/0.2)-np.round(Vr[:, 0]/0.2)) < 1e-6) &
                   (np.abs((Vr[:, 1]/0.2)-np.round(Vr[:, 1]/0.2)) < 1e-6))[0]
    if len(apr) == 0:
        apr = np.where(np.abs(Vr[:, 2]-0.5) < 1e-6)[0][::6]
    frames = []
    for step in range(total):
        if step in pick:
            frames.append(phi.copy())
        phi = phi - dt*Minv*(Cs@phi) + dt*kappa*Minv*(-(Kp@phi))
    fr, _ = norm_frames(frames)
    # steady uniform velocity arrows (same each frame)
    vv = np.tile(np.array(u0), (len(apr), 1))
    av = float(np.max(np.abs(vv)) or 1.0)
    return {"kind": "scalar", "signed": False, "frames": fr,
            "nV": int(nn),
            "verts": [round(float(x), 5) for x in Vr.flatten()],
            "arrows": {"pos": [round(float(x), 4) for x in Vr[apr].flatten()],
                       "vec": [[round(float(x/av), 4) for x in vv.flatten()]]}}


def mesh_block(V, T):
    bnd, wed = faces_edges(T)
    return {"nV": int(len(V)),
            "verts": [round(float(x), 5) for x in V.flatten()],
            "btris": [int(x) for t in bnd for x in t],
            "wedges": [int(x) for e in wed for x in e]}


def main():
    V, T = kuhn_cube(10)                       # [0,1]^3
    K, M = fem_laplacian(V, T)
    tets = per_tet(V, T)
    print("building couette ...");   cou = scene_couette(V, T, K, M)
    print("building poiseuille ..."); poi = scene_poiseuille(V, T, K, M)
    print("building advection ...");  adv = scene_advection(V, T, K, M, tets)
    data = {
        "mesh": mesh_block(V, T),                       # for velocity scenes
        "frames": NFRAMES,
        "scenes": {"couette": cou, "poiseuille": poi, "advection": adv},
    }
    out = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "viewer", "data_fluid.js"))
    with open(out, "w") as f:
        f.write("// Auto-generated by src/verification3d/export_fluid_viewer.py\n")
        f.write("// Time-series snapshots from the verified Stage A/B solvers. Viewer plays back only.\n")
        f.write("window.KSF_FLUID = ")
        json.dump(data, f, separators=(",", ":"))
        f.write(";\n")
    print("wrote", out, "(%d KB)" % (os.path.getsize(out)//1024))
    for k, s in data["scenes"].items():
        print(f"  {k}: kind={s['kind']} frames={len(s['frames'])} arrows={len(s['arrows']['pos'])//3}")


if __name__ == "__main__":
    main()
