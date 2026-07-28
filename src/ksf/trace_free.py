"""
ksf.trace_free
==============
The paper's central constructive object (Part I, Def. 4.2): the trace-free dual
tensor attached to each primal edge.

We implement a concrete, faithful interpretation of the (deliberately
under-specified) text:

  raw dual tensor   ~E(e) = sum_{t > e}  A_t * (n_t (x) n_t)      [Sym^2 R^3]
  trace-free part    E(e) =  ~E(e) - (1/2) tr_g(~E(e)) g|_e

where n_t is the outward unit face normal of a face t incident to edge e, A_t a
dual-area weight, and g|_e the induced metric restricted to the local tangent
plane at the edge midpoint. The projection is done *in the 2-D tangent plane*,
which is what makes the 1/2 factor correct (tr_g g = 2 there).

This module exists to verify two precise, checkable properties the paper asserts:
  (P1)  tr_g E(e) = 0                       (trace-free)
  (P2)  E is coordinate invariant           (frame independence)
It does NOT attempt to show the downstream O(h^2) claim -- that is tested
empirically and honestly in s6 against the exact spectrum.
"""
from __future__ import annotations
import numpy as np
from .dec import _cot  # noqa: F401  (kept for parity with dec angle conventions)


def _tangent_frame(normal: np.ndarray) -> np.ndarray:
    """Orthonormal 3x2 basis [t1 t2] of the plane orthogonal to `normal`."""
    a = np.array([1.0, 0.0, 0.0]) if abs(normal[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    t1 = a - np.dot(a, normal) * normal
    t1 /= np.linalg.norm(t1)
    t2 = np.cross(normal, t1)
    return np.column_stack([t1, t2])          # (3,2)


def _face_data(V, F):
    """Per-face outward unit normal and area."""
    normals = np.zeros((len(F), 3))
    areas = np.zeros(len(F))
    for f, (i, j, k) in enumerate(F):
        cr = np.cross(V[j] - V[i], V[k] - V[i])
        a = np.linalg.norm(cr)
        areas[f] = 0.5 * a
        nrm = cr / a
        # orient outward (sphere centred at origin): align with face centroid
        cen = (V[i] + V[j] + V[k]) / 3.0
        if np.dot(nrm, cen) < 0:
            nrm = -nrm
        normals[f] = nrm
    return normals, areas


def _edge_faces(F):
    """Map each undirected edge to the list of incident face indices."""
    m: dict[tuple[int, int], list[int]] = {}
    for f, (a, b, c) in enumerate(F):
        for u, v in ((a, b), (b, c), (c, a)):
            m.setdefault((min(u, v), max(u, v)), []).append(f)
    return m


def trace_free_dual_tensors(V, F):
    """Return a dict edge -> (E, sphere_normal) where E is the 3x3 trace-free
    dual tensor expressed in ambient coordinates (but trace-free w.r.t. the
    *tangent plane* metric)."""
    normals, areas = _face_data(V, F)
    ef = _edge_faces(F)
    out = {}
    for e, faces in ef.items():
        a, b = e
        mid = 0.5 * (V[a] + V[b])
        nsph = mid / np.linalg.norm(mid)          # surface normal at edge midpoint
        P = np.eye(3) - np.outer(nsph, nsph)       # tangent-plane projector

        raw = np.zeros((3, 3))
        for f in faces:
            n_t = normals[f]
            raw += areas[f] * np.outer(n_t, n_t)
        raw = P @ raw @ P                          # confine to tangent plane

        # trace w.r.t. the 2-D tangent metric == ordinary trace of P-confined tensor
        tr = np.trace(raw)
        E = raw - 0.5 * tr * P                     # (1/2) because dim(tangent)=2
        out[e] = (E, nsph)
    return out


# --------------------------------------------------------------------------- #
#  property checks                                                            #
# --------------------------------------------------------------------------- #
def check_trace_free(V, F) -> float:
    """max |tr_g E(e)| over all edges -- should be ~0 (P1)."""
    tensors = trace_free_dual_tensors(V, F)
    worst = 0.0
    for E, nsph in tensors.values():
        P = np.eye(3) - np.outer(nsph, nsph)
        tr = np.trace(P @ E)                       # tangential trace
        worst = max(worst, abs(tr))
    return float(worst)


def check_frame_invariance(V, F, seed: int = 0) -> float:
    """Rotate the whole mesh by a random R in SO(3); the trace-free tensors must
    transform as R E R^T. Return the max Frobenius mismatch (should be ~0, P2)."""
    rng = np.random.default_rng(seed)
    A = rng.normal(size=(3, 3))
    Q, _ = np.linalg.qr(A)
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]                          # ensure proper rotation

    base = trace_free_dual_tensors(V, F)
    rot = trace_free_dual_tensors((Q @ V.T).T, F)

    worst = 0.0
    for e, (E0, _) in base.items():
        E1 = rot[e][0]
        worst = max(worst, np.linalg.norm(E1 - Q @ E0 @ Q.T))
    return float(worst)
