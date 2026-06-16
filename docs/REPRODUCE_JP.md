# How to Reproduce — 再現手順

This guide explains how to clone this repository, set up the environment,
and run the full verification suite — on both Windows (PowerShell) and
Linux/Mac (including shared hosting servers like Lolipop).

この手順書は、Windows（PowerShell）と Linux サーバー（ロリポップ等）の両方で
検証スイートを再現する方法を説明します。

---

## Requirements / 必要なもの

- Git
- Python 3.7 以上（3.9 以上推奨）
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

`results\results.json` が生成されます。
`viewer\viewer.html` をブラウザで開くとグラフが表示されます（サーバー不要）。

---

## Linux / Mac

### 1. Clone

```bash
git clone https://github.com/ghostinkoma/Skinohedron-MCF-GarbageMesh-Challenge.git
cd Skinohedron-MCF-GarbageMesh-Challenge
```

### 2. Install dependencies

```bash
# pip3 が使える場合
pip3 install --user numpy scipy

# pip3 がない場合（ロリポップ等）
python3 -m pip install --user numpy scipy

# すでに入っているか確認
python3 -c "import numpy, scipy; print('OK', scipy.__version__)"
```

### 3. Run verification

```bash
cd src/verification
mkdir -p ../results
PYTHONPATH=../ python3 run_all.py
```

### 4. Expected output / 期待される出力

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

全セクションがエラーなく完走し、最後に `results.json` が生成されれば成功です。

---

## Shared hosting (Lolipop etc.) / レンタルサーバー

ロリポップのような共有サーバーでも動作確認済みです（Python 3.7 で確認）。

```bash
# SSH接続後
python3 --version        # 3.7以上であることを確認
python3 -m pip --version # pip が使えることを確認

python3 -m pip install --user scipy   # numpy は通常インストール済み

git clone https://github.com/ghostinkoma/Skinohedron-MCF-GarbageMesh-Challenge.git
cd Skinohedron-MCF-GarbageMesh-Challenge/src/verification
mkdir -p ../results
PYTHONPATH=../ python3 run_all.py
```

所要時間：低スペックサーバーで約5〜10分。

---

## Troubleshooting / トラブルシューティング

### `No module named 'ksf'`

`PYTHONPATH` の設定が必要です：

```bash
# src/verification/ から実行する場合
PYTHONPATH=../ python3 run_all.py

# リポルートから実行する場合
PYTHONPATH=src/ python3 src/verification/run_all.py
```

### `No such file or directory: '.../results/results.json'`

`results/` フォルダを先に作る必要があります：

```bash
mkdir -p ../results   # src/verification/ にいる場合
```

### `pip3: コマンドが見つかりません`

```bash
python3 -m pip install --user scipy
```

### GitHub push で認証エラー

GitHub はパスワード認証を廃止しています。
Personal Access Token が必要です：

1. https://github.com/settings/tokens/new を開く
2. `repo` にチェックを入れる
3. `Generate token` でトークンを生成
4. push 時の Password 欄にトークンを貼り付ける

**トークンは公開しないでください。**

---

## Verified environments / 動作確認済み環境

| 環境 | Python | numpy | scipy | 結果 |
|---|---|---|---|---|
| Windows 11 PowerShell | 3.12 | 1.26 | 1.11 | ✅ 全セクション通過 |
| Lolipop レンタルサーバー (Linux) | 3.7 | 1.21 | 1.7.3 | ✅ 全セクション通過 |

---

## Run individual sections / セクション単独実行

```bash
cd src/verification
PYTHONPATH=../ python3 s2_layered_complex.py   # §2 メッシュ品質
PYTHONPATH=../ python3 s3_dec.py               # §3 DEC複体
PYTHONPATH=../ python3 s4_trace_free.py        # §4-5 trace-free射影
PYTHONPATH=../ python3 s6_laplacian_consistency.py  # §6 収束性（核心）
PYTHONPATH=../ python3 s7_nondelaunay.py       # §7 非Delaunayメッシュ
PYTHONPATH=../ python3 s10_sparameter.py       # Part II S行列
```

---

## Author / 著者

Shin-Ichiro Kosaka ([@ghostinkoma](https://github.com/ghostinkoma))

License: MIT
