// data_precision.js -- VERIFIED accuracy data for the precision-improvement viewer.
// All numbers are outputs of verification scripts; the viewer renders presets only,
// no physics in the browser. Sources cited per series. (See theory/06f, 08b, 10.)

const PRECISION = {
  // (3) + high-order enrichment: 1-D eigenvalue error vs local polynomial order.
  // From fluid_highorder_accuracy_verify.py (06f). The route to machine precision
  // for a sinusoidal (field-intrinsic) mode: enrich K.
  enrichment: {
    title: "\u9AD8\u6B21\u5316 (P1\u2192P8): \u5264\u6027\u8AA4\u5DEE\u2192\u6A5F\u68B0\u7CBE\u5EA6",
    source: "06f / fluid_highorder_accuracy_verify.py",
    xlabel: "\u5C40\u6240\u6B21\u6570", ylabel: "\u56FA\u6709\u5024\u8AA4\u5DEE (log)",
    points: [
      {x:"P1", err:2.30e-2}, {x:"P2", err:1.58e-3}, {x:"P3", err:1.37e-4},
      {x:"P4", err:1.29e-5}, {x:"P5", err:2.64e-7}, {x:"P6", err:3.43e-9},
      {x:"P8", err:2.35e-13}
    ]
  },

  // (1) refinement: stiffness O(h^2) error vs mesh n (the approximation part).
  // From fluid_highorder_accuracy_verify.py / fluid_ns_error_decomposition (06d/06f).
  refinement: {
    title: "\u7D30\u5206 (h\u2192): \u5264\u6027\u8AA4\u5DEE\u306E O(h\u00B2) \u53CE\u675F",
    source: "06d/06f",
    xlabel: "\u30E1\u30C3\u30B7\u30E5 n", ylabel: "\u5264\u6027\u8AA4\u5DEE %",
    points: [ {x:"8",err:5.04}, {x:"10",err:3.25}, {x:"12",err:2.26}, {x:"16",err:1.28} ]
  },

  // mass domain: exact at every n (for contrast: enrich K, keep M).
  massExact: { label:"\u8CEA\u91CF M (\u96C6\u4E2D)", values:[
      {x:"8",err:0.0},{x:"10",err:0.0},{x:"12",err:0.0},{x:"16",err:0.0} ], source:"06f" },

  // (1)+(3) combined on a real benchmark: cavity Ghia error vs resolution at Re=1000.
  // From cavity_highre_verify / rev_origin_verify (08b/10). Both elements converge;
  // TH (high-order, =enrichment) is more accurate per node.
  cavity: {
    title: "\u5B9F\u30D9\u30F3\u30C1\u30DE\u30FC\u30AF (Ghia Re=1000): \u7D30\u5206\u30FB\u9AD8\u6B21\u5316\u3067\u8AA4\u5DEE\u6E1B",
    source: "08b/10 / cavity_highre_verify.py",
    xlabel: "\u89E3\u50CF\u5EA6", ylabel: "Ghia\u4E2D\u5FC3\u7DDA RMS",
    pspg: [ {x:"n64",err:0.0555}, {x:"n96",err:0.0404}, {x:"n128",err:0.0320} ],
    th:   [ {x:"n24",err:0.0338}, {x:"n32",err:0.0265}, {x:"n40",err:0.0211} ]
  },

  // The rev separation (10): which part shrinks (1) vs converges (3).
  revSeparation: {
    numerical: { label:"\u6570\u5024\u306E\u30EC\u30D6 (\u4E21\u6CD5\u4E0D\u4E00\u81F4)", source:"10",
      note:"\u7D30\u5206\u3067\u6E1B\u5C11 \u2192 (1)\u8FD1\u4F3C\u8AA4\u5DEE",
      points:[ {x:"\u7C97",err:0.0161}, {x:"\u4E2D",err:0.0140}, {x:"\u5BC6",err:0.0130} ] },
    physical: { label:"\u7269\u7406\u306E\u30EC\u30D6 (\u30B3\u30FC\u30CA\u30FC\u6E26)", source:"10",
      note:"\u4E00\u5B9A\u5024\u306B\u53CE\u675F \u2192 (3)\u5834\u56FA\u6709",
      points:[ {x:"n64",err:0.0989}, {x:"n96",err:0.1413}, {x:"n128",err:0.1414} ] }
  },

  summary: [
    "(1) \u8FD1\u4F3C\u8AA4\u5DEE: \u7D30\u5206\u3067 O(h\u00B2) \u53CE\u675F (5.0\u21921.3%)",
    "(3) \u5834\u56FA\u6709: \u9AD8\u6B21\u5316\u3067\u6A5F\u68B0\u7CBE\u5EA6 (P1 2e-2 \u2192 P8 2e-13)",
    "\u8CEA\u91CF M \u306F\u5E38\u306B\u53B3\u5BC6 \u2014 K \u3092\u8C4A\u304B\u306B\u3001M \u306F\u305D\u306E\u307E\u307E",
    "\u5B9F\u30D9\u30F3\u30C1\u30DE\u30FC\u30AF\u3067\u3082\u7D30\u5206\u30FB\u9AD8\u6B21\u5316\u3067 Ghia \u8AA4\u5DEE\u304C\u6E1B\u308B"
  ]
};
