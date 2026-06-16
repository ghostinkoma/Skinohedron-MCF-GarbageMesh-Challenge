# Interactive Viewer

`viewer.html` is a fully self-contained, offline WebGL viewer.

## How to open

```bash
open viewer/viewer.html        # macOS
xdg-open viewer/viewer.html   # Linux
start viewer\viewer.html       # Windows
```

No web server required. No npm. No build step. Works from `file://`.

## What you see

**Left panel — 3D sphere (WebGL)**
- The unit sphere coloured by a spherical-harmonic field (coolwarm: blue=negative, red=positive)
- Drag to rotate
- Auto-spins slowly
- Toggle wireframe with the W button
- Switch harmonic field: Y₁₀ (dipole), Y₂ₓᵧ (quadrupole), Y₃ₓᵧᵤ (octupole)

**Right panel — convergence chart (Canvas 2D)**
- Log-log plot of error vs mesh size h
- Amber lines = regular (icosphere) mesh
- Coral lines = irregular (persistently jittered) mesh
- Solid = spectral (eigenvalue) error
- Dashed = pointwise error
- Grey reference guides for h¹ and h² slopes

**Bottom — verdict cards**
- One card per paper section
- Green chip = "sound" (claims verified)
- Orange chip = "revised" (claims corrected in v2)

## Data source

The viewer reads from `data.js` (in the same folder), which is automatically
regenerated when you run:

```bash
python3 src/verification/run_all.py
```

`data.js` is committed to the repo so the viewer works without running Python.
