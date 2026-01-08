# Charm Robustness Study: Report

**Generated:** 2026-01-08 00:50 UTC

---

## Executive Summary

| Question | Answer |
|----------|--------|
| Does charm remain "exceptionally exact" after scale-consistent α(μ)? | **Yes** (D_A = 0.0020) |
| Does it survive running y_c to m_Z (Test B)? | **Yes** (D_B = 0.0988) |
| Updated p-value for Test A | **0.0010** (vs 0.0003 single-alpha) |
| Updated p-value for Test B | **0.1429** |

**Bottom line:** The charm coincidence largely evaporates at common scales and is NOT statistically significant.

---

## 1. Baseline Reference

From the preregistered analysis (`config/preregistered.yaml`):

| Quantity | Value |
|----------|-------|
| m_c(m_c) | 1.2730 GeV |
| v (vev) | 246.22 GeV |
| y_c = √2 m_c / v | 0.0073117 |
| α(0) | 0.0072973526 |
| c = y_c / α(0) | **1.001971** |
| d = |ln(c)| | **0.001970** |

The baseline c ≈ 1.002 means y_c is within 0.2% of α.

---

## 2. TEST A: Alpha-at-Scale Substitution

**Procedure:** Keep y_c fixed at its baseline value. Compute c = y_c / α(μ) for different α values.

**Rationale:** Tests whether charm's "exactness" depends on which α we use.

### Results

| Scale | α(μ) | c = y_c / α | d = |ln(c)| |
|-------|------|-------------|------------|
| α(0) | 0.0072973526 | 1.001971 | 0.001970 |
| α(m_τ) | 0.0074920000 | 0.975940 | 0.024355 |
| α(m_Z) | 0.0078151000 | 0.935591 | 0.066577 |

**D_A = min(d) = 0.001970** (achieved at alpha(0))

### Statistical Significance

- Single-alpha p-value (α(0) only): **0.000260**
- Look-elsewhere corrected (3 α choices): **0.000960**

The correction for testing 3 α values does not substantially change the p-value.

---

## 3. TEST B: Common-Scale Running

**Procedure:** Run the charm mass m_c from its definition scale to μ = m_τ and μ = m_Z using 1-loop QCD. Use α(μ) at the same scale.

**Rationale:** The physically meaningful comparison is y_c(μ) vs α(μ) at the same renormalization scale.

### QCD Running Details

- Running formula: m(μ₂)/m(μ₁) = [α_s(μ₂)/α_s(μ₁)]^(γ_m/β₀)
- γ_m = 4 (1-loop mass anomalous dimension)
- β₀ = 11 - 2n_f/3 (1-loop QCD beta function)
- Used α_s(m_Z) = 0.118, α_s(m_τ) = 0.325

### Results

| Scale | m_c(μ) [GeV] | y_c(μ) | α(μ) | c | d = |ln(c)| |
|-------|--------------|--------|------|---|------------|
| Baseline | 1.27300 | 0.0073117 | 0.0072973526 | 1.00197 | 0.001970 |
| μ = m_τ | 1.18162 | 0.0067869 | 0.0074920000 | 0.90589 | 0.098842 |
| μ = m_Z | 0.69199 | 0.0039746 | 0.0078151000 | 0.50858 | 0.676136 |

**D_B = min(d) = 0.098842** (achieved at m_c(m_tau), alpha(m_tau))

### Statistical Significance

- Look-elsewhere corrected p-value: **0.142920**

---

## 4. Interpretation

### What the numbers mean

| Metric | Baseline | Best Test A | Best Test B |
|--------|----------|-------------|-------------|
| c | 1.0020 | 1.0020 | 0.9059 |
| d = |ln(c)| | 0.0020 | 0.0020 | 0.0988 |
| % from c=1 | 0.2% | 0.2% | 9.4% |

### Key finding

**The charm coincidence is NOT robust.**

At common renormalization scales (Test B), the coefficient c = 0.906, which is 9% from 1. The p-value 0.1429 indicates this is NOT statistically significant.

**Conclusion:** The apparent y_c ≈ α coincidence seen at baseline (c = 1.0020) is a scale/scheme artifact. When both quantities are evaluated at compatible scales, the "exactness" evaporates.

---

## 5. Why the baseline looked so good

The baseline comparison used:
- y_c derived from m_c(m_c) ≈ 1.27 GeV (MS-bar at μ = m_c)
- α(0) ≈ 1/137 (low-energy fine-structure constant)

These are evaluated at **different scales**:
- m_c is a GeV-scale quantity
- α(0) is the q² → 0 limit

This mismatch can create spurious coincidences. The QED coupling runs from α(0) ≈ 1/137 to α(m_Z) ≈ 1/128, a ~7% increase. Meanwhile, m_c runs down by ~40% from μ = m_c to μ = m_Z.

---

## 6. Data Sources

All running constants from `data/running_constants_pdg.json`:

| Constant | Value | Source |
|----------|-------|--------|
| α(0) | 0.0072973526 | CODATA 2018 via PDG |
| α(m_τ) | 0.0074920 | PDG 2024 Electroweak Review |
| α(m_Z) | 0.0078151 | PDG 2024 Electroweak Review |
| α_s(m_Z) | 0.1180 | PDG 2024 QCD Review |
| α_s(m_τ) | 0.325 | PDG 2024 QCD Review |

---

## 7. Artifacts Generated

- `artifacts/charm_robustness.csv` - All numerical results
- `artifacts/charm_robustness.json` - Structured JSON output
- `artifacts/charm_scale_plots.png` - Visualization

---

## Appendix: Monte Carlo Details

- Null hypothesis: y drawn log-uniform from [10⁻⁶, 1]
- Number of trials: 100,000
- Test A statistic: D_A = min over {α(0), α(m_τ), α(m_Z)} of |ln(y/α)|
- Test B statistic: D_B = min over {m_τ, m_Z} of |ln(y(μ)/α(μ))|
