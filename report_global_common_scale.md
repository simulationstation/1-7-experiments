# Global Common-Scale Alpha-Ladder Test: Report

**Generated:** 2026-01-08 01:16 UTC

---

## Executive Summary

| Metric | m_τ (1.78 GeV) | m_Z (91.2 GeV) |
|--------|----------------|----------------|
| α(μ) | 0.0074920 | 0.0078151 |
| T(μ) = Σ\|ln c\| | 10.5865 | 10.9684 |
| p-value | 0.339035 | 0.430280 |

**Best scale:** m_tau with T = 10.5865

**Global p-value (look-elsewhere corrected):** 0.578510

**Interpretation:** SM Yukawas show NO global α-ladder structure beyond chance.

---

## 1. Methodology

### Test Statistic
- Primary: T(μ) = Σᵢ |ln(cᵢ)| where cᵢ = yᵢ(μ) / α(μ)^nᵢ
- nᵢ chosen to minimize |ln(cᵢ)| (same rule as baseline)

### Scale Consistency
- All Yukawas y(μ) and α(μ) evaluated at the **same** renormalization scale μ
- Quarks: MS-bar masses run via 1-loop QCD from their definition scales
- Leptons: Pole masses (constant)
- Top: Kept as pole mass (no running)

### Reference Scales (Pre-registered)
- μ = m_τ ≈ 1.78 GeV
- μ = m_Z ≈ 91.2 GeV

### Null Hypothesis
- 9 Yukawas drawn log-uniform from [10⁻⁶, 1]
- Same fitting procedure applied
- 200,000 trials

---

## 2. Results: All 9 Particles

### At μ = m_τ

| Particle | Sector | m(μ) [GeV] | y(μ) | n | c | \|ln c\| | W1 | W2 |
|----------|--------|------------|------|---|------|---------|----|----|
| e | lepton | 0.000511 | 2.9350e-06 | 3 | 6.9794 | 1.9430 | ✓ | ✗ |
| mu | lepton | 0.105658 | 6.0687e-04 | 2 | 10.8119 | 2.3806 | ✗ | ✗ |
| tau | lepton | 1.776930 | 1.0206e-02 | 1 | 1.3623 | 0.3092 | ✓ | ✓ |
| u | up | 0.002212 | 1.2706e-05 | 2 | 0.2264 | 1.4856 | ✓ | ✗ |
| d | down | 0.004814 | 2.7648e-05 | 2 | 0.4926 | 0.7081 | ✓ | ✗ |
| s | down | 0.095759 | 5.5001e-04 | 2 | 9.7989 | 2.2823 | ✓ | ✗ |
| c | up | 1.181623 | 6.7869e-03 | 1 | 0.9059 | 0.0988 | ✓ | ✓ |
| b | down | 5.133610 | 2.9486e-02 | 1 | 3.9357 | 1.3701 | ✓ | ✗ |
| t | up | 172.570000 | 9.9119e-01 | 0 | 0.9912 | 0.0088 | ✓ | ✓ |

**T(m_τ) = 10.5865**

### At μ = m_Z

| Particle | Sector | m(μ) [GeV] | y(μ) | n | c | \|ln c\| | W1 | W2 |
|----------|--------|------------|------|---|------|---------|----|----|
| e | lepton | 0.000511 | 2.9350e-06 | 3 | 6.1491 | 1.8163 | ✓ | ✗ |
| mu | lepton | 0.105658 | 6.0687e-04 | 2 | 9.9364 | 2.2962 | ✓ | ✗ |
| tau | lepton | 1.776930 | 1.0206e-02 | 1 | 1.3060 | 0.2669 | ✓ | ✓ |
| u | up | 0.001307 | 7.5050e-06 | 2 | 0.1229 | 2.0965 | ✓ | ✗ |
| d | down | 0.002843 | 1.6330e-05 | 2 | 0.2674 | 1.3191 | ✓ | ✗ |
| s | down | 0.056561 | 3.2487e-04 | 2 | 5.3191 | 1.6713 | ✓ | ✗ |
| c | up | 0.691990 | 3.9746e-03 | 1 | 0.5086 | 0.6761 | ✓ | ✓ |
| b | down | 3.080282 | 1.7692e-02 | 1 | 2.2639 | 0.8171 | ✓ | ✗ |
| t | up | 172.570000 | 9.9119e-01 | 0 | 0.9912 | 0.0088 | ✓ | ✓ |

**T(m_Z) = 10.9684**

---

## 3. Summary Statistics

| Metric | μ = m_τ | μ = m_Z |
|--------|---------|---------|
| T = Σ\|ln c\| | 10.5865 | 10.9684 |
| T² = Σ(ln c)² | 19.3423 | 18.6964 |
| Pass W1 [0.1,10] | 8/9 | 9/9 |
| Pass W2 [0.5,2] | 3/9 | 3/9 |

---

## 4. Null Hypothesis Testing

### Null Model
- Draw 9 random Yukawas log-uniform in [10⁻⁶, 1]
- Apply same fitting procedure
- Trials: 200,000

### Results

| Scale | T_obs | T_null (mean) | p-value |
|-------|-------|---------------|---------|
| m_τ | 10.5865 | 11.4548 | 0.339035 |
| m_Z | 10.9684 | 11.3286 | 0.430280 |

**Global (look-elsewhere corrected):**
- T_best_obs = 10.5865
- T_best_null (mean) = 10.2146
- **p_global = 0.578510**

---

## 5. Robustness Check: Excluding Top

| Metric | m_τ | m_Z |
|--------|-----|-----|
| T (8 particles) | 10.5777 | 10.9596 |
| p-value | 0.574025 | 0.672320 |
| p_global | 0.829345 | - |

The exclusion of top (which has scheme ambiguities) changes the conclusion.

---

## 6. Interpretation

### Is the SM Yukawa set globally more α-ladder-like than chance?

**NO.**

The global p-value is 0.579, meaning 58% of random Yukawa sets achieve T ≤ 10.59 at their best scale. This is entirely consistent with chance.

The α-ladder hypothesis receives no support from this scale-consistent analysis.

### Key Observations

1. **Scale matters:** The choice of μ affects T significantly. At m_tau, T = 10.59.

2. **Pass counts are not discriminative:** 8/9 pass W1 at m_τ, but random sets also achieve similar pass rates.

3. **The tight window W2 [0.5, 2]:** Only 3/9 pass at m_τ, 3/9 at m_Z. This is a more stringent test.

---

## 7. Conclusion

**The scale-consistent global test finds NO evidence for α-ladder structure in SM Yukawas.**

The test statistic T (sum of |ln c|) is not unusually small compared to random Yukawa sets under the null hypothesis.

---

## Artifacts

- `artifacts/global_common_scale.csv`
- `artifacts/global_common_scale.json`
- `artifacts/global_common_scale_plots.png`
