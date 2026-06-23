// data_model.js  --  VERIFIED solver/script output for the V3 model viewer.
// Every number here was produced by a verification script in src/verification3d/
// and is cited to that script. NO physics runs in the browser; this is a map of
// verified facts, presets only. (See theory/07_v3_synthesis.md and 06z audit.)

const MODEL = {
  operator: { name: "L = M\u207B\u00B9K", note: "one cotangent/FE operator (Kuhn-cube P1): lumped mass M, stiffness K from one per-tet gradient G" },

  // The four domains and their measured error character
  domains: [
    { id: "mass", label: "Mass  M", role: "inertia",
      precision: "exact", color: "#2dd4bf",
      fact: "u\u1D40M u = 0.2500 at every n (0.00% error)",
      detail: "The lumped-mass form equals the continuous integral exactly, at all resolutions. The inertial domain carries NO discretisation error.",
      script: "fluid_highorder_accuracy_verify.py", doc: "06f" },

    { id: "coulomb", label: "Coulomb  K", role: "stiffness / viscosity / electrostatics",
      precision: "O(h\u00B2)", color: "#f59e0b",
      fact: "u\u1D40K u low by 5.0\u21921.3% (n=8\u219216);  1/r law 0.8%",
      detail: "All smooth-mode error lives in the stiffness. The same K gives electrostatics (K\u03C6=M\u03C1, 1/r to 0.8%, superposition 4e-16). Route to machine precision: local enrichment (see rule).",
      script: "coulomb_operator_verify.py / fluid_highorder_accuracy_verify.py", doc: "06e / 06f" },

    { id: "temperature", label: "Temperature  T", role: "parallel domain + glue",
      precision: "fuse + mediate", color: "#a78bfa",
      fact: "fusion total-E drift 0.0001%;  glue interface u=U\u00B5\u2082/(\u00B5\u2081+\u00B5\u2082) to 7.5e-16",
      detail: "Additive (06g): velocity\u2295temperature Boussinesq conserves total energy (no secular drift), incompressibility 1.9e-17. Multiplicative (06h): \u03BD(T) weights the stiffness K_\u03BD; two-viscosity interface is the mechanical twin of 02's T=k\u2082/(k\u2081+k\u2082), machine precision.",
      script: "thermal_fusion_verify.py / thermal_glue_verify.py", doc: "06g / 06h" },

    { id: "velocity", label: "Velocity  u", role: "momentum / incompressibility",
      precision: "div-free (machine) \u00B7 pressure-decoupled", color: "#60a5fa",
      fact: "\u2016Bu\u2016~1e-16 (algebraically trivial);  PSPG moves u only ~1e-4",
      detail: "Divergence-free to machine precision \u2014 but that alone is an algebraic property of any projection. The EARNED result: the velocity is pressure-decoupled (fixing the checkerboard pressure, corr 0.0\u21920.94, changes u by only ~1e-4). Pressure itself is NOT trustworthy (inf-sup).",
      script: "fluid_navier_stokes_verify.py / fluid_ns_velocity_quality_verify.py", doc: "06c" }
  ],

  // Readings of the one operator (operator identities)
  readings: [
    { eq: "\u1E6A = \u2212\u03BA L T", name: "heat / diffusion", doc: "02, 03" },
    { eq: "p\u0308 = \u2212c\u00B2 L p", name: "scalar wave", doc: "01, 03" },
    { eq: "K p = (\u03C1/\u0394t)\u2207\u00B7u*", name: "pressure / incompressibility", doc: "04, 06c" },
    { eq: "K \u03C6 = M \u03C1", name: "electrostatics / Coulomb", doc: "06e" },
    { eq: "u\u0307 = \u2212\u03BD L u + \u2026", name: "momentum (viscous)", doc: "06a\u201306c" }
  ],

  // The machine-precision rule (audit-corrected)
  rule: [
    { kind: "linear \u2192 exact everywhere (in P1)", ex: "Couette; T=k\u2082/(k\u2081+k\u2082); u=\u00B5\u2082/(\u00B5\u2081+\u00B5\u2082)", val: "~1e-15", cls: "exact" },
    { kind: "polynomial-source \u2192 nodally exact (superconvergence)", ex: "parabolic Poiseuille", val: "1.9e-15", cls: "exact" },
    { kind: "eigenmode / sinusoid \u2192 O(h\u00B2)", ex: "Taylor\u2013Green decay; thermal sine; \u03BB_h vs 2k\u00B2", val: "5.0\u21921.3%", cls: "oh2" }
  ],

  // Enrichment route to machine precision (1-D, verified)
  enrichment: [
    { order: "P1", err: 2.30e-2 }, { order: "P2", err: 1.58e-3 },
    { order: "P3", err: 1.37e-4 }, { order: "P4", err: 1.29e-5 },
    { order: "P5", err: 2.64e-7 }, { order: "P6", err: 3.43e-9 },
    { order: "P8", err: 2.35e-13 }
  ],

  // Honest scope (from 06z audit)
  earned: [
    "operator identities (heat/wave/pressure/Coulomb) \u2014 exact, pressure-independent",
    "domain error split: mass exact, Coulomb O(h\u00B2), enrichment \u2192 2e-13",
    "velocity divergence-free AND pressure-decoupled (~1e-4)",
    "temperature fusion invariant (0.0001%) and glue interface (7.5e-16)"
  ],
  notEarned: [
    "\u2016Bu\u2016~1e-16 alone is an algebraic triviality \u2014 not a quality result",
    "equal-order P1 pressure is inf-sup unstable / checkerboard \u2014 NOT trustworthy",
    "battery is idealised: periodic Taylor\u2013Green + analytic interfaces; no cavity/cylinder/turbulence/experiment",
    "the unifying narrative is interpretation on idealised tests \u2014 wider than the proof",
    "no first-principles \u03BD from molecular Coulomb (Green\u2013Kubo/MD out of scope)"
  ],

  // vs state-of-the-art (so the position is not overstated)
  sota: [
    "SOTA has inf-sup-stable pressure (Taylor\u2013Hood / PSPG); here pressure is unreliable",
    "SOTA reaches high order in 3-D; here P1 (high order shown in 1-D only)",
    "SOTA does turbulence (DNS/LES/RANS); here moderate-Re laminar",
    "SOTA uses adaptive unstructured meshes; here structured Kuhn cube n=8\u201316",
    "SOTA runs at 10\u2076\u201310\u2079 DOF on HPC/GPU; here small, direct solvers",
    "SOTA is validated against experiments/benchmarks; here analytic solutions only"
  ]
};
