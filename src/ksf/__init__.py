"""
Kosaka Skin-o-hedron Model (KSF) -- reference implementation & verification suite.

A small, dependency-light library (numpy + scipy) that turns the prose claims of
the KSF manuscript into runnable, falsifiable numerical experiments on the
sphere, where the Laplace-Beltrami spectrum is known exactly.

Submodules
----------
mesh        : the S_k icosphere refinement sequence (Part I, Def. 2.1)
dec         : discrete exterior calculus operators + cotangent Laplacian (Part I, §3)
trace_free  : trace-free dual tensors and their invariance checks (Part I, §4-5)
sph         : exact spherical-harmonic ground truth
"""
from . import mesh, dec, trace_free, sph   # noqa: F401

__all__ = ["mesh", "dec", "trace_free", "sph"]
__version__ = "2.0.0"
