# Solid-viewer API (grand finale)

`viewer_solid.html` shows the same operator `L = M⁻¹K` driving heat, scalar wave,
and fluid (liquid/gas) physics on three solids — cube, sphere, torus. Values come
two ways, both verified:

- **Preset** (default): instant, from `viewer/data_solid.js`
  (precomputed by `src/verification3d/export_solid_viewer.py`). Static hosting is
  enough — no server logic needed.
- **Solve on server** (custom params): the viewer POSTs to `api/solve_solid.php`,
  which pipes the request to `api/solve_solid.py` (the same solvers as the
  verification scripts) and returns the vertex field as JSON.

## Files
- `solve_solid.php` — whitelists shape∈{cube,ball,torus} and physics∈{heat,wave,
  liquid,gas}, clamps numeric params, pipes JSON to Python via `proc_open`.
- `solve_solid.py` — CLI solver. Reads JSON on stdin, writes JSON on stdout.

## Requirements
- PHP with `proc_open`; `python3` on PATH with numpy/scipy.
- Deployed inside the repo so `solve_solid.py` reaches `../../src/ksf3d`.

## Quick test
```sh
echo '{"shape":"torus","physics":"wave","steps":120}' | python3 solve_solid.py
echo '{"shape":"ball","physics":"liquid","g":2.0}'    | python3 solve_solid.py
```
Each prints `{"ok":true,...,"p":[...]}` (normalised to max|p|=1).

## Request / response
Request: `{ "shape":"cube|ball|torus", "physics":"heat|wave|liquid|gas",
  "steps":60, "g":1.0, "c2":0.18 }`
Response: `{ "ok":true, "shape":..., "physics":..., "nV":N, "signed":bool, "p":[...] }`

If the API is unreachable, the viewer stays on presets and shows a small note — it
never breaks.
