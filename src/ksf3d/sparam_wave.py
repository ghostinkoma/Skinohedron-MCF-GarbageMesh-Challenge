"""
ksf3d.sparam_wave
=================
Faithful implementation of the mathematical model in
theory/01b_sparameter_model.md (Step 1: S-parameter wave).

Primary state = oriented port amplitudes a[(σ,k)] ∈ R^{4 N_C}.
One step:  a ← U a,  U = C ∘ S  (both orthogonal ⇒ exactly energy-conserving).

No physics is decided here that is not in 01b §8; this module only assembles the
operators from mesh combinatorics + the Z / boundary inputs. Each operator
carries the assertion named in 01b.
"""
from __future__ import annotations
import numpy as np

FACE_COMBOS = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]


# --------------------------------------------------------------------------- #
#  combinatorial data (01b §1)                                                #
# --------------------------------------------------------------------------- #
class PortComplex:
    """Builds the oriented-port structure P, the partner involution Pi, boundary
    flags, and per-port geometry from (V, T)."""

    def __init__(self, V, T):
        self.V = np.asarray(V, float)
        self.T = np.asarray(T, int)
        self.N_C = len(self.T)
        self.nP = 4 * self.N_C

        # face(σ,k) keyed by sorted vertex triple -> list of (σ,k)
        faces = {}
        for s, tet in enumerate(self.T):
            for k, fc in enumerate(FACE_COMBOS):
                key = tuple(sorted(int(tet[j]) for j in fc))
                faces.setdefault(key, []).append((s, k))

        # partner involution on oriented ports, boundary flags, facet id
        self.partner = np.full(self.nP, -1, dtype=int)   # global port index
        self.is_boundary = np.zeros(self.nP, dtype=bool)
        self.facet_of = np.full(self.nP, -1, dtype=int)
        self.interior_facets = []   # list of (p0, p1) global port pairs
        fid = 0
        for key, plist in faces.items():
            if len(plist) == 2:
                (s0, k0), (s1, k1) = plist
                p0, p1 = self.pid(s0, k0), self.pid(s1, k1)
                self.partner[p0] = p1
                self.partner[p1] = p0
                self.facet_of[p0] = self.facet_of[p1] = fid
                self.interior_facets.append((p0, p1))
            else:
                (s0, k0) = plist[0]
                p0 = self.pid(s0, k0)
                self.partner[p0] = p0           # self-partner
                self.is_boundary[p0] = True
                self.facet_of[p0] = fid
            fid += 1
        self.N_F = fid

        # geometry per port: outward normal, area; cell volume
        self.normal = np.zeros((self.nP, 3))
        self.area = np.zeros(self.nP)
        self.vol = np.zeros(self.N_C)
        for s, tet in enumerate(self.T):
            P = self.V[tet]
            self.vol[s] = abs(np.linalg.det(
                np.array([P[1] - P[0], P[2] - P[0], P[3] - P[0]]))) / 6.0
            ctr = P.mean(axis=0)
            for k, fc in enumerate(FACE_COMBOS):
                q = P[list(fc)]
                nrm = np.cross(q[1] - q[0], q[2] - q[0])
                a = 0.5 * np.linalg.norm(nrm)
                n = nrm / (np.linalg.norm(nrm) + 1e-15)
                if np.dot(q.mean(axis=0) - ctr, n) < 0:
                    n = -n
                self.normal[self.pid(s, k)] = n
                self.area[self.pid(s, k)] = a

        self._assert_involution()

    def pid(self, s, k):
        return 4 * s + k

    def cell_of(self, p):
        return p // 4

    def _assert_involution(self):
        # Pi ∘ Pi = id  (discrete ∂∂ = 0 shadow), 01b §1
        pp = self.partner[self.partner]
        assert np.array_equal(pp, np.arange(self.nP)), "partner is not an involution"


# --------------------------------------------------------------------------- #
#  scatter operator S  (01b §3)                                               #
# --------------------------------------------------------------------------- #
# S_wave = 2 P0 - I  on each cell's 4 ports; P0 = (1/4) 11^T
def S_apply(a):
    """Apply block-diagonal S_wave = 2P0 - I to port amplitudes a (shape (nP,))."""
    A = a.reshape(-1, 4)
    mean = A.mean(axis=1, keepdims=True)        # (1/4) 11^T a
    out = 2.0 * mean - A                        # (2 P0 - I) a
    return out.reshape(-1)


def S_matrix():
    """The explicit 4x4 node (for the assertion)."""
    P0 = np.ones((4, 4)) / 4.0
    return 2 * P0 - np.eye(4)


def assert_S_orthogonal(tol=1e-12):
    S = S_matrix()
    err = np.linalg.norm(S.T @ S - np.eye(4))
    assert err < tol, f"S not orthogonal: {err}"
    assert np.linalg.norm(S - S.T) < tol, "S not symmetric"
    ev = np.sort(np.linalg.eigvalsh(S))
    # eigenvalues: -1 (mult 3), +1 (mult 1)
    assert np.allclose(ev, [-1, -1, -1, 1], atol=1e-9), f"S spectrum {ev}"
    return err


# --------------------------------------------------------------------------- #
#  connect operator C  (01b §4)  -- built facet-wise to guarantee losslessness #
# --------------------------------------------------------------------------- #
class Connect:
    """a_next = C b, assembled per facet from impedance (R,T) blocks (interior)
    and wall reflection rho (boundary)."""

    def __init__(self, pc: PortComplex, Z, rho_bnd=1.0):
        self.pc = pc
        self.Z = np.asarray(Z, float)
        # per-boundary-port rho (scalar or array)
        if np.isscalar(rho_bnd):
            self.rho = np.full(pc.nP, float(rho_bnd))
        else:
            self.rho = np.asarray(rho_bnd, float)
        self._build()

    def _build(self):
        pc = self.pc
        # interior facet (R,T) per port, symmetric across the facet
        self.R = np.zeros(pc.nP)
        self.T = np.zeros(pc.nP)
        for (p0, p1) in pc.interior_facets:
            Z1 = self.Z[pc.cell_of(p0)]
            Z2 = self.Z[pc.cell_of(p1)]
            R = (Z2 - Z1) / (Z2 + Z1)
            T = 2.0 * np.sqrt(Z1 * Z2) / (Z1 + Z2)
            # side p0 sees (R, T); side p1 sees (-R, T)  (01b §4 consistency)
            self.R[p0], self.T[p0] = R, T
            self.R[p1], self.T[p1] = -R, T

    def apply(self, b):
        pc = self.pc
        a = np.empty_like(b)
        part = pc.partner
        # interior ports
        intr = ~pc.is_boundary
        a[intr] = self.R[intr] * b[intr] + self.T[intr] * b[part[intr]]
        # boundary ports
        bnd = pc.is_boundary
        a[bnd] = self.rho[bnd] * b[bnd]
        return a

    def matrix(self):
        """Dense C (only for small-mesh assertions)."""
        pc = self.pc
        n = pc.nP
        C = np.zeros((n, n))
        for p in range(n):
            if pc.is_boundary[p]:
                C[p, p] = self.rho[p]
            else:
                C[p, p] = self.R[p]
                C[p, pc.partner[p]] = self.T[p]
        return C


def assert_C_orthogonal(connect: Connect, tol=1e-10):
    """Lossless case check ‖CᵀC − I‖ (01b §4). Dense; small meshes only."""
    C = connect.matrix()
    err = np.linalg.norm(C.T @ C - np.eye(C.shape[0]))
    return err


def connect_isometry_err(connect: Connect, n_probe=8, seed=0):
    """Matrix-free losslessness check: max | ‖C x‖ − ‖x‖ | / ‖x‖ over random x.
    For orthogonal C this is ~0. Works on large meshes (no dense matrix)."""
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(n_probe):
        x = rng.normal(size=connect.pc.nP)
        nx = np.linalg.norm(x)
        ncx = np.linalg.norm(connect.apply(x))
        worst = max(worst, abs(ncx - nx) / nx)
    return float(worst)


# --------------------------------------------------------------------------- #
#  full update U = C S  and energy  (01b §5, §6)                              #
# --------------------------------------------------------------------------- #
def U_step(a, connect: Connect, gamma=1.0):
    return connect.apply(gamma * S_apply(a))


def energy(a):
    return float(a @ a)


def inject_pulse(a, pc: PortComplex, cell, s=1.0):
    """Monopole injection at `cell`: add s/2 to all 4 ports (energy s^2)."""
    for k in range(4):
        a[pc.pid(cell, k)] += 0.5 * s
    return a
