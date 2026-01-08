# alpha-locked-clockwork — Yukawa ladder test (public PDG data)

This repository implements a *public-data* test inspired by the **Alpha-Locked Clockwork (ALC)** draft
(`ClockworkEMPredictions.pdf`), specifically the **integer-exponent texture** idea:

> small couplings can be expressible as  \( y \sim c\,\alpha^n \) with **integer** \(n\) and \(c = O(1\text{–}10)\).  
> (See Eq. (4) in the draft.)

We test that ansatz on **Standard Model Yukawa eigenvalues** (i.e. singular values of the Yukawa matrices),
computed from PDG fermion masses.

## Why this test

- Yukawas are **dimensionless** couplings, so this avoids the “cross-section → coupling” mediator-mass degeneracy.
- The data are public and widely used (PDG summary tables).
- The test is *not* a proof of clockwork layers; it is a clean “numbers-as-written” stress test.

## What this code does

1. Loads a pinned snapshot of PDG masses (quarks + charged leptons) from `data/pdg_masses_2024_2025.json`.
2. Converts masses to naive Yukawas:
   \[
   y_f = \sqrt{2}\, m_f / v
   \]
   using a fixed electroweak vev \(v\) (configurable).
3. Fits each \(y_f\) to \(c\,\alpha^n\) with integer \(n\) and bounded \(c\) using a pre-registered rule:
   - scan \(n\in[n_{\min}, n_{\max}]\)
   - compute \(c = y/\alpha^n\)
   - choose the \(n\) that makes \(|\log c|\) smallest (i.e. \(c\) closest to 1), optionally requiring \(c\) within bounds
4. Propagates mass uncertainties via Monte Carlo and produces:
   - best-fit \(n\) and \(c\) per fermion
   - pass/fail counts under the chosen coefficient bounds
   - stability summaries (how often the same \(n\) is selected across draws)

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

python -m alc_yukawa --config config/preregistered.yaml fit

# or a stricter interpretation of c ∈ [1,10]
python -m alc_yukawa --config config/preregistered_strict.yaml fit
```

Outputs (tables + JSON artifacts) are written to `artifacts/` by default.

## Results

### Summary (preregistered config, c ∈ [0.1, 10])

| Particle | Sector | Yukawa | n | c | In Bounds |
|----------|--------|--------|---|------|-----------|
| t | up | 0.991 | 0 | 0.991 | Yes |
| c | up | 0.00731 | 1 | 1.002 | Yes |
| u | up | 1.24e-5 | 2 | 0.233 | Yes |
| b | down | 0.0240 | 1 | 3.292 | Yes |
| s | down | 5.37e-4 | 2 | 10.08 | **No** |
| d | down | 2.70e-5 | 2 | 0.507 | Yes |
| tau | lepton | 0.0102 | 1 | 1.399 | Yes |
| mu | lepton | 6.07e-4 | 2 | 11.40 | **No** |
| e | lepton | 2.94e-6 | 3 | 7.553 | Yes |

**Pass count: 7/9** (muon and strange quark fail with c slightly above 10)

### Null-hypothesis testing (100,000 simulations)

To assess whether 7/9 is meaningful, we simulated random log-uniform Yukawas in [10⁻⁶, 1]:

| Test | Observed | Null expectation | p-value | Significant? |
|------|----------|------------------|---------|--------------|
| Absolute fits pass rate | 7/9 | 8.4/9 (mean) | 0.98 | **No** (worse than chance) |
| Ratio fits pass rate | 9/9 | 8.6/9 (mean) | 0.67 | **No** |
| Charm |log(c)| | 0.002 | — | 0.006 | **Yes** (even after look-elsewhere) |

### Interpretation

**The overall α-ladder pattern does not beat random chance.** The observed 7/9 pass rate is actually *below* what random numbers achieve on average. The generous coefficient bounds ([0.1, 10] spans a factor of 100, comparable to the α ≈ 1/137 step size) make most numbers pass by construction.

**The one statistically significant result is charm:** y_c ≈ α¹ × 1.002, sitting almost exactly on the first ladder rung. The probability of any of the 9 fermions being this close by chance is ~0.6% (p = 0.006).

See `report.md` for full analysis including Monte Carlo stability, ratio fits, and method critique.

### Charm Robustness Study: Is the charm coincidence real?

The charm result (c = 1.002) looked striking, so we tested whether it survives **scale-consistent analysis**:

```bash
python -m alc_yukawa --config config/charm_robustness.yaml charm-study
```

**TEST A** (vary α, keep y_c fixed):
| Scale | α(μ) | c = y_c/α | d = \|ln(c)\| |
|-------|------|-----------|---------------|
| α(0) | 1/137.04 | 1.002 | 0.002 |
| α(m_τ) | 1/133.5 | 0.976 | 0.024 |
| α(m_Z) | 1/128.0 | 0.936 | 0.067 |

**TEST B** (run m_c and α to common scale via QCD):
| Scale | m_c(μ) | y_c(μ) | α(μ) | c | d |
|-------|--------|--------|------|---|---|
| Baseline | 1.273 GeV | 0.00731 | 1/137 | **1.002** | 0.002 |
| μ = m_τ | 1.182 GeV | 0.00679 | 1/133.5 | 0.906 | 0.099 |
| μ = m_Z | 0.692 GeV | 0.00397 | 1/128 | 0.509 | 0.676 |

**Key result:**
| Analysis | d = \|ln(c)\| | p-value | Significant? |
|----------|---------------|---------|--------------|
| Baseline α(0) only | 0.002 | 0.0003 | Yes |
| Test A (look-elsewhere) | 0.002 | 0.001 | Yes |
| **Test B (common scale)** | **0.099** | **0.143** | **No** |

**Conclusion: The charm coincidence is a scale/scheme artifact.**

The baseline comparison used mismatched scales (y_c at μ ~ 1.3 GeV vs α at q² = 0). At compatible renormalization scales, charm's coefficient degrades to c = 0.91 with p = 0.14—not statistically significant.

See `report_charm_robustness.md` for full details.

### Global Scale-Consistent Test (Final Analysis)

After finding that the charm coincidence was a scale artifact, we implemented a **global** test: Are ALL SM Yukawas simultaneously close to y = c·α^n when both y and α are evaluated at the **same** renormalization scale μ?

```bash
python -m alc_yukawa --config config/common_scale_global.yaml global-study
```

**Method:**
- Evaluate all 9 fermion Yukawas and α at common scales (μ = m_τ, μ = m_Z)
- Use QCD running for quark masses; leptons are pole masses
- Test statistic: T(μ) = Σ|ln(c_i)| (sum over all fermions)
- 200,000 null simulations with log-uniform random Yukawas

**Results at μ = m_τ (1.78 GeV):**
| Particle | Yukawa | n | c | ln(c) | In Bounds |
|----------|--------|---|-------|-------|-----------|
| t | 0.9914 | 0 | 0.991 | -0.009 | Yes |
| b | 0.0275 | 1 | 3.67 | 1.30 | Yes |
| c | 0.0068 | 1 | 0.906 | -0.099 | Yes |
| s | 6.21e-4 | 2 | 11.1 | 2.40 | **No** |
| d | 3.13e-5 | 2 | 0.557 | -0.585 | Yes |
| u | 1.46e-5 | 2 | 0.260 | -1.35 | Yes |
| tau | 0.0102 | 1 | 1.36 | 0.31 | Yes |
| mu | 6.07e-4 | 2 | 10.8 | 2.38 | **No** |
| e | 2.94e-6 | 3 | 7.01 | 1.95 | Yes |

**Pass count: 7/9** at both scales

**Null hypothesis testing:**
| Metric | μ = m_τ | μ = m_Z | Global |
|--------|---------|---------|--------|
| T(μ) = Σ\|ln(c)\| | 10.59 | 10.97 | min = 10.59 |
| p-value | 0.34 | 0.43 | **0.58** |
| Pass [0.1,10] | 7/9 (p=0.97) | 7/9 (p=0.99) | — |
| Pass [0.5,2] | 3/9 (p=0.84) | 2/9 (p=0.99) | — |

**Excluding top quark:** p_global = 0.83 (even worse)

**Conclusion: The α-ladder hypothesis finds NO statistically significant support in the SM Yukawa spectrum.** When both Yukawas and α are evaluated at consistent renormalization scales, the observed pattern is no better than random numbers. The previously noted charm coincidence (c ≈ 1) was an artifact of comparing quantities at mismatched scales.

See `report_global_common_scale.md` for full details.

## Notes / caveats

- **Quark masses are scheme/scale dependent.** This repo uses the PDG *summary table* masses:
  - \(m_{u,d,s}(\mu=2\,\mathrm{GeV})\) in \(\overline{\mathrm{MS}}\)
  - \(m_c(m_c)\), \(m_b(m_b)\) in \(\overline{\mathrm{MS}}\)
  - \(m_t\) from direct measurements
  This is adequate for a first-pass “ladder proximity” test, but a more careful version would run all Yukawas to a common \(\mu_*\).

- PDG quark masses in the summary table are quoted at **90% CL**. The config lets you decide whether to:
  - treat the quoted \(\pm\) as 1σ (not recommended), or
  - convert to 1σ by dividing by 1.64485 (Gaussian assumption).

## Repo layout

- `src/alc_yukawa/` — library + CLI
- `config/` — pre-registered analysis knobs (alpha, bounds, seeds, etc.)
- `data/` — pinned PDG snapshot used for reproducibility
- `artifacts/` — generated outputs (gitignored)

## License
MIT (see `LICENSE`).

