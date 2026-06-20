# Pressure API (server-side solve)

`viewer_pressure.html` shows **verified** pressure fields two ways:

- **Preset** (default): instant, read from `viewer/data_pressure.js`
  (precomputed by `src/verification3d/export_pressure_viewer.py`). Works with no
  server logic — static hosting is enough.
- **Solve on server** (custom parameters): the viewer POSTs to `api/solve.php`,
  which pipes the request to `api/solve_pressure.py` (the *same* solver as
  `pressure_field_verify.py`) and returns the vertex pressure array as JSON.

Both paths produce values that pass the Step 2 verification — that is the point.

## Files

- `solve.php` — thin, safe wrapper: whitelists the mode, clamps numeric ranges,
  pipes JSON to the Python solver via `proc_open`, returns its stdout.
- `solve_pressure.py` — CLI solver. Reads JSON on stdin, writes JSON on stdout.
  Adds the repo's `src/` to `sys.path` so `ksf3d` is importable.

## Requirements on the server

- PHP with `proc_open` enabled.
- `python3` on PATH with `numpy` and `scipy`.
- The repo checked out so that `solve_pressure.py` can reach
  `../../src/ksf3d/` (it is, when deployed inside the repo's `viewer/api/`).

## Quick test (shell)

```sh
echo '{"mode":"hydrostatic","rho":1.0,"g":[0,0,-1.0]}' | python3 solve_pressure.py
echo '{"mode":"atmosphere","g":[0,0,-1.0],"c2":0.15}'  | python3 solve_pressure.py
echo '{"mode":"acoustic","c2":1.0,"steps":120}'        | python3 solve_pressure.py
```

Each prints `{"ok":true,...,"p":[...]}`.

## Request / response

Request (POST JSON):
```json
{ "mode":"hydrostatic|atmosphere|acoustic",
  "rho":1.0, "g":[0,0,-1.0], "c2":0.15, "steps":120 }
```
Response:
```json
{ "ok":true, "mode":"hydrostatic", "nV":729, "signed":false, "p":[ ... ] }
```
`p` is normalised to `max|p| = 1` for display.

## Notes

- If the API is unreachable (e.g. static-only hosting), the viewer stays on
  presets and shows a small "API unreachable (use preset)" note — it never breaks.
- The solver rebuilds the operator each call (n=8, 729 nodes — cheap). For heavier
  meshes, cache `K`/`M` server-side or precompute more presets.
