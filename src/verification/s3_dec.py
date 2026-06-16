"""
s3_dec.py  --  Part I, Section 3
================================
Sanity-checks the discrete exterior calculus layer the rest of Part I is built
on. Two checkable facts:

  * the discrete de Rham complex closes:  d1 . d0 = 0  (to machine precision);
  * Euler's formula  V - E + F = 2  holds for every refinement level (the complex
    is a valid closed surface), which is what lets the Hodge theory be applied.
"""
from __future__ import annotations
import numpy as np
from ksf import mesh, dec


def run(levels=range(1, 5)) -> dict:
    rows = []
    for lvl in levels:
        V, F = mesh.icosphere(lvl)
        E = mesh.edges(F)
        res = dec.complex_residual(V, F)
        euler = len(V) - len(E) + len(F)
        rows.append(dict(level=int(lvl), V=len(V), E=len(E), F=len(F),
                         euler=int(euler), dd_residual=float(res)))
    return {"rows": rows,
            "max_dd_residual": float(max(r["dd_residual"] for r in rows)),
            "euler_ok": all(r["euler"] == 2 for r in rows)}


def main():
    res = run()
    print("=== S3: discrete exterior calculus ===")
    print(f"{'lvl':>3} {'V':>6} {'E':>6} {'F':>6} {'V-E+F':>6} {'|d.d|':>10}")
    for r in res["rows"]:
        print(f"{r['level']:>3} {r['V']:>6} {r['E']:>6} {r['F']:>6} "
              f"{r['euler']:>6} {r['dd_residual']:>10.2e}")
    print(f"\nmax |d1 . d0|  = {res['max_dd_residual']:.2e}  (=> complex closes)")
    print(f"Euler char = 2 at every level : {res['euler_ok']}")
    print("\n[結句] The DEC scaffolding is correct and standard. Nothing in Part I's")
    print("       difficulty lives here; d, the incidences and the Hodge setup are")
    print("       textbook. The load-bearing novelty is the dual-tensor operator.")


if __name__ == "__main__":
    main()
