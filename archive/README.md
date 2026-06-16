# archive/original_cpp/

This directory contains the original C++ prototype files from the first version
of this repository.

## Files

- `Skinohedron.hpp` — First implementation of the trace-free dual tensor operator,
  written in C++ using Eigen. This was a proof-of-concept sketch.
- `main.cop` — Original driver / notes file.

## Status

These files are **archived for historical reference** and are **not maintained**.

The active implementation is in Python (`src/ksf/`) and is fully verified.
The C++ prototype served its purpose: it demonstrated the core idea works
before the full mathematical analysis was done.

## What changed between the prototype and the current version

The prototype had a subtle issue: the 3D→2D basis projection was ad-hoc
(`block<2,2>(0,0)` of the 3×3 outer product), which is not coordinate-invariant.
The current implementation works directly with the metric tensor in the local
tangent frame, making it provably frame-invariant.

See `src/verification/s4_trace_free.py` for the verification.
