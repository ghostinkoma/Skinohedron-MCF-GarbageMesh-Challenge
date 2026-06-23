// data_precision.js -- VERIFIED accuracy data for the precision-improvement viewer.
// All numbers are outputs of verification scripts; the viewer renders presets only,
// no physics in the browser. Sources cited per series. (See theory/06f, 08b, 10.)

const PRECISION = {
  // (3) + high-order enrichment: 1-D eigenvalue error vs local polynomial order.
  // From fluid_highorder_accuracy_verify.py (06f). The route to machine precision
  // for a sinusoidal (field-intrinsic) mode: enrich K.
  enrichment: {
    title: "高次化 (P1→P8): 剤性誤差→機械精度",
    source: "06f / fluid_highorder_accuracy_verify.py",
    xlabel: "局所次数", ylabel: "固有値誤差 (log)",
    points: [
      {x:"P1", err:2.30e-2}, {x:"P2", err:1.58e-3}, {x:"P3", err:1.37e-4},
      {x:"P4", err:1.29e-5}, {x:"P5", err:2.64e-7}, {x:"P6", err:3.43e-9},
      {x:"P8", err:2.35e-13}
    ]
  },

  // (1) refinement: stiffness O(h^2) error vs mesh n (the approximation part).
  // From fluid_highorder_accuracy_verify.py / fluid_ns_error_decomposition (06d/06f).
  refinement: {
    title: "細分 (h→): 剤性誤差の O(h²) 収束",
    source: "06d/06f",
    xlabel: "メッシュ n", ylabel: "剤性誤差 %",
    points: [ {x:"8",err:5.04}, {x:"10",err:3.25}, {x:"12",err:2.26}, {x:"16",err:1.28} ]
  },

  // mass domain: exact at every n (for contrast: enrich K, keep M).
  massExact: { label:"質量 M (集中)", values:[
      {x:"8",err:0.0},{x:"10",err:0.0},{x:"12",err:0.0},{x:"16",err:0.0} ], source:"06f" },

  // (1)+(3) combined on a real benchmark: cavity Ghia error vs resolution at Re=1000.
  // From cavity_highre_verify / rev_origin_verify (08b/10). Both elements converge;
  // TH (high-order, =enrichment) is more accurate per node.
  cavity: {
    title: "実ベンチマーク (Ghia Re=1000): 細分・高次化で誤差減",
    source: "08b/10 / cavity_highre_verify.py",
    xlabel: "解像度", ylabel: "Ghia中心線 RMS",
    pspg: [ {x:"n64",err:0.0555}, {x:"n96",err:0.0404}, {x:"n128",err:0.0320} ],
    th:   [ {x:"n24",err:0.0338}, {x:"n32",err:0.0265}, {x:"n40",err:0.0211} ]
  },

  // The rev separation (10): which part shrinks (1) vs converges (3).
  revSeparation: {
    numerical: { label:"数値のレブ (両法不一致)", source:"10",
      note:"細分で減少 → (1)近似誤差",
      points:[ {x:"粗",err:0.0161}, {x:"中",err:0.0140}, {x:"密",err:0.0130} ] },
    physical: { label:"物理のレブ (コーナー渦)", source:"10",
      note:"一定値に収束 → (3)場固有",
      points:[ {x:"n64",err:0.0989}, {x:"n96",err:0.1413}, {x:"n128",err:0.1414} ] }
  },

  summary: [
    "(1) 近似誤差: 細分で O(h²) 収束 (5.0→1.3%)",
    "(3) 場固有: 高次化で機械精度 (P1 2e-2 → P8 2e-13)",
    "質量 M は常に厳密 — K を豊かに、M はそのまま",
    "実ベンチマークでも細分・高次化で Ghia 誤差が減る"
  ]
};
