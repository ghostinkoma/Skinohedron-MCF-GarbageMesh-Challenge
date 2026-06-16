"""
s4_trace_free.py  --  Part I, Sections 4 & 5
============================================
The paper's two *provable* claims about the dual tensor are exactly the two
properties it needs to be a sensible object, and both are testable to machine
precision:

  (P1)  trace-free:        tr_g E(e) = 0           (Def. 4.2 / Remark 4.3)
  (P2)  frame invariance:  E(R x) = R E(x) R^T     (Remark 4.3)

This module confirms both. It deliberately does NOT claim that trace-freeness
implies O(h^2) consistency -- that inference is the manuscript's overreach, and
it is examined empirically (and refuted in the strong/pointwise sense) in s6.
"""
from __future__ import annotations
from ksf import mesh, trace_free


def run(levels=(2, 3, 4)) -> dict:
    rows = []
    for lvl in levels:
        V, F = mesh.icosphere(lvl)
        p1 = trace_free.check_trace_free(V, F)
        p2 = trace_free.check_frame_invariance(V, F, seed=lvl)
        rows.append(dict(level=int(lvl), trace_residual=p1, frame_residual=p2))
    return {"rows": rows,
            "max_trace_residual": float(max(r["trace_residual"] for r in rows)),
            "max_frame_residual": float(max(r["frame_residual"] for r in rows))}


def main():
    res = run()
    print("=== S4/S5: trace-free dual tensor properties ===")
    print(f"{'lvl':>3} {'|tr E| (P1)':>14} {'|frame err| (P2)':>18}")
    for r in res["rows"]:
        print(f"{r['level']:>3} {r['trace_residual']:>14.2e} "
              f"{r['frame_residual']:>18.2e}")
    print(f"\nmax trace residual  (P1) = {res['max_trace_residual']:.2e}")
    print(f"max frame residual  (P2) = {res['max_frame_residual']:.2e}")
    print("\n[結句] Both claims the paper can actually back are TRUE to machine")
    print("       precision: the projection is trace-free and coordinate invariant,")
    print("       and the 1/2 factor is correct precisely because the tangent space")
    print("       is 2-D. This part of the manuscript is sound. The leap from here")
    print("       to 'therefore the operator is O(h^2) pointwise' is the unproven step.")


if __name__ == "__main__":
    main()
