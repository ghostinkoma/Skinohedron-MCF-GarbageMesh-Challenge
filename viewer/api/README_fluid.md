# Fluid-viewer API (V3 Stage A/B)

`viewer_fluid.html` plays back **time-series snapshots** of the verified Stage A
(Stokes) and Stage B (advection–diffusion) solvers. Every frame is a verified
solver output — the browser only renders and animates, it does no physics.

- **Preset** (default): instant, from `viewer/data_fluid.js`
  (precomputed by `src/verification3d/export_fluid_viewer.py`). Static hosting is
  enough.
- **Solve on server** (custom params): the viewer POSTs to `api/solve_fluid.php`,
  which pipes to `api/solve_fluid.py` (same operators as the verification scripts)
  and returns a fresh time-series.

## Scenes
- `couette`    — Stokes startup, wall-driven → linear profile `u = U y`.
- `poiseuille` — Stokes startup, pressure-driven → parabolic profile.
- `advection`  — Stage B Gaussian scalar carried by uniform flow + diffusing.

## Files
- `solve_fluid.php` — whitelists scene, clamps params, pipes JSON to Python.
- `solve_fluid.py`  — CLI solver (stdin→stdout), reuses `fem_laplacian`, per-tet `G`.

## Requirements
PHP `proc_open`; `python3` with numpy/scipy; repo layout so `../../src/ksf3d` is importable.

## Test
```sh
echo '{"scene":"advection","ux":1.0,"kappa":0.004,"frames":36}' | python3 solve_fluid.py
echo '{"scene":"poiseuille","G":-2.0,"nu":1.0}'                 | python3 solve_fluid.py
```
Each prints `{"ok":true,...,"frames":[...]}`.

If the API is unreachable the viewer stays on presets — it never breaks.
