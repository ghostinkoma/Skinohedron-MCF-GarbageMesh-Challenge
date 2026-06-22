# Fluid-on-shapes viewer (V3) — cube · torus · sphere

`viewer_fluid_shapes.html` plays back advection animations on the three solids, each
with its natural flow, using snapshots precomputed by
`src/verification3d/export_fluid_shapes_viewer.py` from the SAME operators as
`fluid_shapes_verify.py`. Every frame is a verified-solver output; the browser only
renders and animates.

- **cube**  — uniform flow (x-periodic): the scalar translates and wraps.
- **torus** — flow around the ring (flat, exact): the scalar travels around and returns.
- **sphere**— solid-body rotation `u = ω×r` on the **geometric (degraded) mesh**.
  The on-screen note flags it: ~2844 sign-flipped cotangents, FE accuracy reduced.
  Incompressibility and mass conservation stay exact regardless (see
  `fluid_shapes_verify.py`); only the transported-field accuracy degrades.

Overlay: velocity arrows + scalar heatmap (toggle each). View: wireframe,
transparency, slice, auto-rotate. Transport bar: play/pause/scrub.

This viewer is **presets only** (static hosting is enough) — the animations are the
verified snapshots. No in-browser physics.

Regenerate data:
```sh
python3 src/verification3d/export_fluid_shapes_viewer.py   # -> viewer/data_fluid_shapes.js
```
