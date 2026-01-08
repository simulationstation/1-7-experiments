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

