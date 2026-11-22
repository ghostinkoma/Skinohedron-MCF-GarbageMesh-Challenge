// data_precision.js -- VERIFIED accuracy data for the precision-improvement viewer.
// All numbers are outputs of verification scripts; the viewer renders presets only,
// no physics in the browser. Sources cited per series. (See theory/06f, 08b.)
// Some producing scripts were reconstructed in 2026-11 (the originals were never committed);
// see theory/README.md, "Reconstructed and restored items".

const PRECISION = {
  // (3) + high-order enrichment: 1-D eigenvalue error vs local polynomial order.
  // From fluid_highorder_accuracy_verify.py (06f). The route to machine precision
  // for a sinusoidal (field-intrinsic) mode: enrich K.
  enrichment: {
    title: "高次化 (P1→P8, 1次元のみ): 剛性誤差→機械精度",
    source: "06f / fluid_highorder_accuracy_verify.py",
    xlabel: "局所次数", ylabel: "固有値誤差 (log)",
    points: [
      {x:"P1", err:2.30e-2}, {x:"P2", err:1.58e-3}, {x:"P3", err:1.37e-4},
      {x:"P4", err:1.29e-5}, {x:"P5", err:2.64e-7}, {x:"P6", err:3.43e-9},
      {x:"P8", err:2.35e-13}
    ]
  },

  // (1) refinement: stiffness O(h^2) error vs mesh n (the approximation part).
  // From fluid_highorder_accuracy_verify.py / fluid_ns_error_decomposition_verify.py (06f/06d).
  refinement: {
    title: "細分 (h→): 剛性誤差の O(h²) 収束",
    source: "06d / 06f",
    xlabel: "メッシュ n", ylabel: "剛性誤差 %",
    points: [ {x:"8",err:5.04}, {x:"10",err:3.25}, {x:"12",err:2.26}, {x:"16",err:1.28} ]
  },

  // mass domain: exact at every n FOR THIS TEST FIELD (whole-period sinusoid on a uniform
  // mesh -> trapezoidal-rule exactness); O(h^2) in general (mass_quadrature_verify.py).
  massExact: { label:"質量 M (集中・このテスト関数では厳密)", values:[
      {x:"8",err:0.0},{x:"10",err:0.0},{x:"12",err:0.0},{x:"16",err:0.0} ], source:"06f" },

  // (1)+(3) combined on a real benchmark: cavity Ghia error vs resolution at Re=1000.
  // PSPG n64 and TH n32 from cavity_highre_verify.py (08b); PSPG n96/n128 and TH n24/n40
  // from rev_origin_verify.py (10). TH (high-order, =enrichment) is more accurate per node.
  cavity: {
    title: "実ベンチマーク (Ghia Re=1000): 細分・高次化で誤差減",
    source: "08b / 10 / rev_origin_verify.py",
    xlabel: "解像度", ylabel: "Ghia中心線 RMS",
    pspg: [ {x:"n64",err:0.0555}, {x:"n96",err:0.0404}, {x:"n128",err:0.0320} ],
    th:   [ {x:"n24",err:0.0338}, {x:"n32",err:0.0265}, {x:"n40",err:0.0211} ]
  },

  // The rev separation: which part shrinks (1) vs the physical corner eddy.
  // CORRECTED (issue #3, A3): the physical series is now the bottom-right eddy's stream
  // function at Re=1000 (PSPG, cavity_highre_verify.py check C). The earlier series
  // (0.0989 / 0.1413 / 0.1414) was a u_min proxy that measured the primary return flow.
  revSeparation: {
    numerical: { label:"数値のレブ (両法不一致)", source:"10 / rev_origin_verify.py",
      note:"細分で減少 → (1)近似誤差",
      // PSPG-vs-TH mutual RMS at Re=1000 for (n96,TH32)/(n112,TH36)/(n128,TH40); the middle
      // value was 0.0140 from an unrecorded pair, re-measured as 0.0142 (rev_origin_verify.py).
      points:[ {x:"粗",err:0.0161}, {x:"中",err:0.0142}, {x:"密",err:0.0130} ] },
    physical: { label:"物理のレブ (右下二次渦 ψ, Re=1000)", source:"08b",
      note:"n≤128 では増加中 (Ghia の 55→67→74%)。延長計算 (10 §3b) で PSPG・TH ともゆっくり収束し、外挿値 ≈1.9e-3 で一致",
      points:[ {x:"n64",err:9.585e-4}, {x:"n96",err:1.181e-3}, {x:"n128",err:1.304e-3} ] }
  },

  summary: [
    "(1) 近似誤差: 細分で O(h²) 収束 (5.0→1.3%)",
    "(3) 場固有: 高次化で機械精度 (P1 2e-2 → P8 2e-13、1次元のみ)",
    "質量 M はこのテスト関数（一様メッシュ上の周期正弦波）では厳密。一般の場では O(h²)",
    "実ベンチマークでも細分・高次化で Ghia 誤差が減る"
  ],

  // (2) Coulomb/EHD and (4) reaction residual fingerprints, and (1)+(2) compose (12).
  residuals: {
    title: "予想(2)(4)の残差: 機械NSが再現できない構造を項が埋める",
    source: "12 / residual_domains_verify.py",
    coulomb: { label:"(2) クーロン/EHD", off:"純機械NS 最大速度 0 (駆動なし)",
      on:"クーロンON 0.007 (電荷に線形: 0.003/0.007/0.013)",
      note:"Kφ=Mρ_q (06eの同じK) が機械NSの残差を埋める" },
    reaction: { label:"(4) 反応 (Fisher-KPP)", off:"純移流拡散 総量保存 0.060 (薄まるだけ)",
      on:"反応ON 総量 0.125 (増殖・フロント伝播, c→1)",
      note:"反応項が移流拡散の残差(フロント)を埋める, 残差 0.142" },
    compose: { label:"(1)+(2) 合成", points:[
      {x:"n32",err:0.00660},{x:"n48",err:0.00657},{x:"n64",err:0.00658} ],
      note:"EHD流れが細分で収束(連続差 3e-5→1e-5) = 高次化(1)がクーロン駆動流(2)の精度を上げる" }
  }
};
