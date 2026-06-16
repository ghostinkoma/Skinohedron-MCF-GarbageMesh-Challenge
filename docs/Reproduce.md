# How to Reproduce

This guide explains how to clone this repository, set up the environment,
and run the full verification suite on Windows (PowerShell), Linux/Mac,
and shared hosting servers.

---

## Requirements

- Git
- Python 3.7 or later (3.9+ recommended)
- numpy >= 1.21
- scipy >= 1.7

---

## Windows (PowerShell)

### 1. Clone

```powershell
cd C:\Users\YourName\Desktop
git clone https://github.com/ghostinkoma/Skinohedron-MCF-GarbageMesh-Challenge.git
cd Skinohedron-MCF-GarbageMesh-Challenge
```

### 2. Install dependencies

```powershell
pip install numpy scipy
```

### 3. Run verification

```powershell
cd src\verification
$env:PYTHONPATH = ".."
python run_all.py
```

### 4. View results

`results\results.json` will be generated.
Open `viewer\viewer.html` in any browser to see interactive charts — no server needed.

---

## Linux / Mac

### 1. Clone

```bash
git clone https://github.com/ghostinkoma/Skinohedron-MCF-GarbageMesh-Challenge.git
cd Skinohedron-MCF-GarbageMesh-Challenge
```

### 2. Install dependencies

```bash
# if pip3 is available
pip3 install --user numpy scipy

# if pip3 is not found (common on shared hosting)
python3 -m pip install --user numpy scipy

# verify
python3 -c "import numpy, scipy; print('OK', scipy.__version__)"
```

### 3. Run verification

```bash
cd src/verification
mkdir -p ../results
PYTHONPATH=../ python3 run_all.py
```

### 4. Expected output

```
######################################################################
#  Kosaka Skin-o-hedron Model -- full verification run
######################################################################

=== S2: Skin-o-hedron sequence diagnostics ===
...
surface-area defect order  ~ h^2.00

=== S3: discrete exterior calculus ===
...
max |d1 . d0|  = 0.00e+00  (=> complex closes)

=== S4/S5: trace-free dual tensor properties ===
...
max trace residual  (P1) = 1.08e-19

=== S6: regular icosphere  (harmonic Y2xy) ===
...
  pointwise consistency ~ h^1.06
  spectral convergence  ~ h^3.68

=== S7: non-Delaunay behaviour ===
...

=== S10-12: S-parameter operator ===
...

######################################################################
#  wrote .../results/results.json  (535 KB)
######################################################################
```

All sections complete without error, and `results.json` is written at the end.

---

## Shared hosting (e.g. Lolipop)

Verified on a shared Linux server running Python 3.7 at under $2/month.

```bash
# after SSH login
python3 --version         # confirm 3.7 or later
python3 -m pip --version  # confirm pip is available

python3 -m pip install --user scipy   # numpy is usually pre-installed

git clone https://github.com/ghostinkoma/Skinohedron-MCF-GarbageMesh-Challenge.git
cd Skinohedron-MCF-GarbageMesh-Challenge/src/verification
mkdir -p ../results
PYTHONPATH=../ python3 run_all.py
```

Estimated runtime: 5–10 minutes on a low-spec shared server.

---

## Troubleshooting

### `No module named 'ksf'`

The `PYTHONPATH` environment variable must point to `src/`:

```bash
# when running from src/verification/
PYTHONPATH=../ python3 run_all.py

# when running from the repo root
PYTHONPATH=src/ python3 src/verification/run_all.py
```

### `No such file or directory: '.../results/results.json'`

Create the `results/` directory first:

```bash
mkdir -p ../results   # when inside src/verification/
```

### `pip3: command not found`

Use the module invocation instead:

```bash
python3 -m pip install --user scipy
```

### GitHub push authentication error

GitHub no longer accepts passwords for `git push`.
A Personal Access Token (PAT) is required:

1. Go to https://github.com/settings/tokens/new
2. Check the `repo` scope
3. Click `Generate token` and copy the result
4. Paste the token into the `Password` prompt when pushing

**Never share your token publicly.**

---

## Run individual sections

```bash
cd src/verification
PYTHONPATH=../ python3 s2_layered_complex.py        # §2  mesh quality
PYTHONPATH=../ python3 s3_dec.py                    # §3  DEC complex
PYTHONPATH=../ python3 s4_trace_free.py             # §4-5 trace-free projection
PYTHONPATH=../ python3 s6_laplacian_consistency.py  # §6  convergence (key result)
PYTHONPATH=../ python3 s7_nondelaunay.py            # §7  non-Delaunay meshes
PYTHONPATH=../ python3 s10_sparameter.py            # Part II  S-parameter kernel
```

---

## Verified environments

| Environment | Python | numpy | scipy | Result |
|---|---|---|---|---|
| Windows 11 PowerShell | 3.12 | 1.26 | 1.11 | ✅ all sections pass |
| Lolipop shared server (Linux) | 3.7 | 1.21 | 1.7.3 | ✅ all sections pass |

---

## Author

Shin-Ichiro Kosaka ([@ghostinkoma](https://github.com/ghostinkoma))

License: MIT
