// data_model.js  --  VERIFIED solver/script output for the V3 model viewer.
// Every number here was produced by a verification script in src/verification3d/
// and is cited to that script. NO physics runs in the browser; this is a map of
// verified facts, presets only. (See theory/07_v3_synthesis.md and 06z audit.)

const MODEL = {
  operator: { name: "L = M⁻¹K", note: "one cotangent/FE operator (Kuhn-cube P1): lumped mass M, stiffness K from one per-tet gradient G" },

  // The four domains and their measured error character
  domains: [
    { id: "mass", label: "Mass  M", role: "inertia",
      precision: "exact", color: "#2dd4bf",
      fact: "uᵀM u = 0.2500 at every n (0.00% error)",
      detail: "The lumped-mass form equals the continuous integral exactly, at all resolutions. The inertial domain carries NO discretisation error.",
      script: "fluid_highorder_accuracy_verify.py", doc: "06f" },

    { id: "coulomb", label: "Coulomb  K", role: "stiffness / viscosity / electrostatics",
      precision: "O(h²)", color: "#f59e0b",
      fact: "uᵀK u low by 5.0→1.3% (n=8→16);  1/r law 0.8%",
      detail: "All smooth-mode error lives in the stiffness. The same K gives electrostatics (Kφ=Mρ, 1/r to 0.8%, superposition 4e-16). Route to machine precision: local enrichment (see rule).",
      script: "coulomb_operator_verify.py / fluid_highorder_accuracy_verify.py", doc: "06e / 06f" },

    { id: "temperature", label: "Temperature  T", role: "parallel domain + glue",
      precision: "fuse + mediate", color: "#a78bfa",
      fact: "fusion total-E drift 0.0001%;  glue interface u=Uµ₂/(µ₁+µ₂) to 7.5e-16",
      detail: "Additive (06g): velocity⊕temperature Boussinesq conserves total energy (no secular drift), incompressibility 1.9e-17. Multiplicative (06h): ν(T) weights the stiffness K_ν; two-viscosity interface is the mechanical twin of 02's T=k₂/(k₁+k₂), machine precision.",
      script: "thermal_fusion_verify.py / thermal_glue_verify.py", doc: "06g / 06h" },

    { id: "velocity", label: "Velocity  u", role: "momentum / incompressibility",
      precision: "div-free (machine) · pressure-decoupled", color: "#60a5fa",
      fact: "‖Bu‖~1e-16 (algebraically trivial);  PSPG moves u only ~1e-4",
      detail: "Divergence-free to machine precision — but that alone is an algebraic property of any projection. The EARNED result: the velocity is pressure-decoupled (fixing the checkerboard pressure, corr 0.0→0.94, changes u by only ~1e-4).",
      script: "fluid_navier_stokes_verify.py / fluid_ns_velocity_quality_verify.py", doc: "06c" },

    { id: "pressure", label: "Pressure  p", role: "incompressibility constraint (P0)",
      precision: "stabilised → benchmark-valid", color: "#34d399",
      fact: "PSPG cavity vs Ghia 1982 Re=100: u,v RMS 0.0095→0.0069 (n=48→64); pressure smooth",
      detail: "Raw equal-order P1 pressure is inf-sup-unstable (06z). PSPG-stabilised, the lid-driven cavity matches the published Ghia benchmark on both centrelines (converging) with a SMOOTH pressure — validated against an EXTERNAL standard. Cross-checked: Taylor-Hood P2-P1 (inf-sup-stable, NO parameter) independently matches Ghia, and agrees with PSPG to RMS 0.0015 — so the match is correct physics, not a τ-tuning artifact. The stabilisation is −τK, the SAME operator K. Pressure is a correlation domain on L=M⁻¹K, in 2-D and 3-D. (Only Re=100 so far.)",
      script: "cavity_pspg_verify.py / cavity_taylorhood_verify.py / cavity3d_demo.py", doc: "08 / 08a" }
  ],

  // Readings of the one operator (operator identities)
  readings: [
    { eq: "Ṫ = −κ L T", name: "heat / diffusion", doc: "02, 03" },
    { eq: "p̈ = −c² L p", name: "scalar wave", doc: "01, 03" },
    { eq: "K p = (ρ/Δt)∇·u*", name: "pressure / incompressibility", doc: "04, 06c" },
    { eq: "K φ = M ρ", name: "electrostatics / Coulomb", doc: "06e" },
    { eq: "u̇ = −ν L u + …", name: "momentum (viscous)", doc: "06a–06c" }
  ],

  // The machine-precision rule (audit-corrected)
  rule: [
    { kind: "linear → exact everywhere (in P1)", ex: "Couette; T=k₂/(k₁+k₂); u=µ₂/(µ₁+µ₂)", val: "~1e-15", cls: "exact" },
    { kind: "polynomial-source → nodally exact (superconvergence)", ex: "parabolic Poiseuille", val: "1.9e-15", cls: "exact" },
    { kind: "eigenmode / sinusoid → O(h²)", ex: "Taylor–Green decay; thermal sine; λ_h vs 2k²", val: "5.0→1.3%", cls: "oh2" }
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
    "operator identities (heat/wave/pressure/Coulomb) — exact, pressure-independent",
    "domain error split: mass exact, Coulomb O(h²), enrichment → 2e-13",
    "velocity divergence-free AND pressure-decoupled (~1e-4)",
    "temperature fusion invariant (0.0001%) and glue interface (7.5e-16)"
  ],
  notEarned: [
    "‖Bu‖~1e-16 alone is an algebraic triviality — not a quality result",
    "RAW equal-order P1 pressure is checkerboard; STABILISED (PSPG) it matches Ghia Re=100 — but only Re=100, and PSPG is not a true inf-sup element",
    "cavity now validates the foundation; cylinder / turbulence / time-dependence still untested",
    "the unifying narrative is interpretation on idealised tests — wider than the proof",
    "no first-principles ν from molecular Coulomb (Green–Kubo/MD out of scope)"
  ],

  // vs state-of-the-art (so the position is not overstated)
  sota: [
    "SOTA has inf-sup-stable pressure (Taylor–Hood / PSPG); here pressure is unreliable",
    "SOTA reaches high order in 3-D; here P1 (high order shown in 1-D only)",
    "SOTA does turbulence (DNS/LES/RANS); here moderate-Re laminar",
    "SOTA uses adaptive unstructured meshes; here structured Kuhn cube n=8–16",
    "SOTA runs at 10⁶–10⁹ DOF on HPC/GPU; here small, direct solvers",
    "SOTA is validated against experiments/benchmarks; here analytic solutions only"
  ]
};
