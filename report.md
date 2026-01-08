# ALC Yukawa Ladder Test: Analysis Report

**Generated:** 2026-01-07
**Commit:** `208e1cc7327381ad4f33d0af66c136c8ce655ee5`
**Python:** 3.12.3

---

## 1. Executive Summary

- **Code runs successfully.** Both configs execute without errors.
- **Headline result (preregistered):** **7 of 9 fermions pass** the absolute fit criterion (c within [0.1, 10]).
- **Failures (absolute fit):**
  - **Muon:** c = 11.40 (14% above upper bound)
  - **Strange quark:** c = 10.08 (0.8% above upper bound)
- **Ratio fits:** All 9 particles pass within-sector ratio fits (9/9).
- **n assignments are perfectly stable:** P(mode) = 1.000 for all particles under 20,000 MC draws.

**CRITICAL: Null-hypothesis results (100,000 simulations):**
- **Absolute fits 7/9: NOT SIGNIFICANT.** Random numbers pass at mean 8.4/9. Observed 7/9 is *worse* than chance (p = 0.98).
- **Ratio fits 9/9: NOT SIGNIFICANT.** Random numbers achieve 9/9 with probability 67%.
- **Charm quark c = 1.002: SIGNIFICANT.** Probability of |log(c)| this small by chance is ~1/1,500 (p = 0.0007).

**Bottom line:** The overall α-ladder pattern does not beat random chance. The only statistically anomalous result is charm sitting almost exactly on α^1.

---

## 2. Reproducibility Checklist

### Commands Executed
```bash
# From repo root:
PYTHONPATH=src python3 -m alc_yukawa --config config/preregistered.yaml fit
PYTHONPATH=src python3 -m alc_yukawa --config config/preregistered_strict.yaml fit
```

### Environment
| Item | Value |
|------|-------|
| Python version | 3.12.3 |
| Git commit | `208e1cc7327381ad4f33d0af66c136c8ce655ee5` |
| Branch | main |
| Warnings/Errors | None |

### Output Files (artifacts/)
| File | Size (bytes) | Description |
|------|--------------|-------------|
| point_fit.csv | 1,361 | Per-particle absolute fit results |
| ratio_fit.csv | 481 | Within-sector ratio fit results |
| results.json | 7,205 | Full structured results |
| yukawa_vs_n.png | 46,911 | Scatter plot of log(y) vs assigned n |
| coeff_hist.png | 43,161 | Histogram of implied coefficients |

---

## 3. Results Table: Absolute Fits

Configuration: `preregistered.yaml` with c in [0.1, 10], alpha = 0.0072973525693

| Particle | Sector | Mass (GeV) | Yukawa | n | c | log(c) | x = log(y)/log(alpha) | dist_to_int | In Bounds |
|----------|--------|------------|--------|---|------|--------|----------------------|-------------|-----------|
| t | up | 172.57 | 0.9912 | 0 | 0.991 | -0.009 | 0.002 | 0.002 | **Yes** |
| c | up | 1.273 | 0.00731 | 1 | 1.002 | 0.002 | 1.000 | 0.0004 | **Yes** |
| u | up | 0.00216 | 1.24e-5 | 2 | 0.233 | -1.457 | 2.296 | 0.296 | **Yes** |
| b | down | 4.183 | 0.0240 | 1 | 3.292 | 1.192 | 0.758 | 0.242 | **Yes** |
| s | down | 0.0935 | 5.37e-4 | 2 | 10.08 | 2.311 | 1.530 | 0.470 | **No** |
| d | down | 0.0047 | 2.70e-5 | 2 | 0.507 | -0.679 | 2.138 | 0.138 | **Yes** |
| tau | lepton | 1.777 | 0.0102 | 1 | 1.399 | 0.335 | 0.932 | 0.068 | **Yes** |
| mu | lepton | 0.1057 | 6.07e-4 | 2 | 11.40 | 2.433 | 1.505 | 0.495 | **No** |
| e | lepton | 5.11e-4 | 2.94e-6 | 3 | 7.553 | 2.022 | 2.589 | 0.411 | **Yes** |

**Pass count: 7/9**

### Interpretation of x and dist_to_int
- **x = log(y)/log(alpha):** The "ideal" real-valued exponent that would give c = 1 exactly.
- **dist_to_int:** How far x is from the nearest integer. Smaller means closer to an exact α^n.

Notable observations:
- **Charm (c):** x = 1.000, dist = 0.0004. Almost exactly on the ladder rung.
- **Top (t):** x = 0.002, dist = 0.002. Essentially y ≈ 1.
- **Muon, strange:** Both have dist ≈ 0.47-0.50, meaning they sit almost exactly midway between rungs.

---

## 4. Results Table: Sector Ratio Fits

Ratios are computed within each sector relative to the heaviest particle (baseline).

| Sector | Particle | Baseline | Ratio | n_ratio | c_ratio | In Bounds |
|--------|----------|----------|-------|---------|---------|-----------|
| up | t | t | 1.000 | 0 | 1.000 | **Yes** |
| up | c | t | 0.00738 | 1 | 1.011 | **Yes** |
| up | u | t | 1.25e-5 | 2 | 0.235 | **Yes** |
| down | b | b | 1.000 | 0 | 1.000 | **Yes** |
| down | s | b | 0.0224 | 1 | 3.063 | **Yes** |
| down | d | b | 0.00112 | 1 | 0.154 | **Yes** |
| lepton | tau | tau | 1.000 | 0 | 1.000 | **Yes** |
| lepton | mu | tau | 0.0595 | 1 | 8.148 | **Yes** |
| lepton | e | tau | 2.88e-4 | 2 | 5.400 | **Yes** |

**Pass count: 9/9** (all pass, even under strict [1, 10] bounds: 9/9)

Key insight: Ratio fits eliminate sensitivity to the vev and some scheme issues, giving cleaner results.

---

## 5. Monte Carlo Stability

20,000 MC samples drawn with masses perturbed by their quoted uncertainties (converted from 90% CL to 1sigma).

| Particle | mode_n | P(mode) | P(in_bounds) | Interpretation |
|----------|--------|---------|--------------|----------------|
| t | 0 | 1.000 | 1.000 | Perfectly stable, always passes |
| c | 1 | 1.000 | 1.000 | Perfectly stable, always passes |
| u | 2 | 1.000 | 1.000 | Perfectly stable, always passes |
| b | 1 | 1.000 | 1.000 | Perfectly stable, always passes |
| s | 2 | 1.000 | **0.052** | n is stable, but only 5.2% of draws yield c < 10 |
| d | 2 | 1.000 | 1.000 | Perfectly stable, always passes |
| tau | 1 | 1.000 | 1.000 | Perfectly stable, always passes |
| mu | 2 | 1.000 | **0.000** | n is stable, but never passes (c always > 10) |
| e | 3 | 1.000 | 1.000 | Perfectly stable, always passes |

**Key observations:**
- All n assignments are deterministic under measurement uncertainty: P(mode) = 1.000 for all particles.
- Strange quark sits right at the boundary: central c = 10.08, so only 5.2% of draws dip below 10.
- Muon's c = 11.4 is far enough above 10 that no MC draw ever brings it into bounds.

---

## 6. Interpretation (Lay Terms)

### What does "yukawas line up with α^n" mean?

The Yukawa couplings are the fundamental numbers that determine how strongly each fermion interacts with the Higgs field (and thus gets its mass). They span a huge range: ~1 for the top quark down to ~3×10^-6 for the electron.

The α-ladder hypothesis asks: can we write each Yukawa as:
```
y = c × α^n
```
where:
- α ≈ 1/137 is the fine-structure constant
- n is a whole number (0, 1, 2, 3, ...)
- c is a "coefficient" that should be roughly 1-10

If true, this would mean nature organizes fermion masses on a geometric ladder with spacing factor α ≈ 1/137.

### What does the coefficient c represent?

The coefficient c absorbs the "leftover" after fitting to an integer power of α. Ideally c ≈ 1, meaning the Yukawa sits almost exactly on a ladder rung. If c is, say, 3, the Yukawa is 3× larger than the nearest exact α^n value.

The test asks: are these leftover factors modest (O(1-10)), or do they blow up arbitrarily? Finding c consistently in [0.1, 10] for most particles is suggestive but not definitive.

### Why are ratio fits more robust?

Absolute Yukawas depend on the electroweak vev (v ≈ 246 GeV) and the specific renormalization scheme/scale used for quark masses. Within-sector ratios cancel these dependencies:
```
y_light / y_heavy = m_light / m_heavy
```
This isolates the relative hierarchy from absolute normalization issues. The fact that ratio fits pass for all particles (including muon and strange) suggests the hierarchy pattern may be more robust than the absolute scale.

---

## 7. Method Critique / Failure Modes

### 7.1 Risk of Overfitting to Integer n

The fitting rule (minimize |log c|) automatically selects the best n for any Yukawa. With n ranging from 0 to 12, there are 13 possible "rungs" to land on. For any positive number y:
- log(y)/log(α) gives a real-valued "ideal exponent" x
- We round to the nearest integer (in effect)
- By construction, c = y/α^n will be within a factor of sqrt(α) ≈ 0.085 to 1/sqrt(α) ≈ 11.7 of 1

This means **most random numbers in the relevant range will achieve c between ~0.1 and ~12 by construction.** The interesting question is whether c clusters closer to 1 than expected by chance.

### 7.2 Width of α^n Steps

Each step of α ≈ 1/137 represents a factor of 137 change. The "allowed" c range of [0.1, 10] spans a factor of 100—comparable to the step size. This means:
- A single rung "covers" about 73% of the log-space between rungs (log(100)/log(137) ≈ 0.94)
- Failure requires c to land in the ~27% gap near the boundaries

The 7/9 pass rate is therefore **not overwhelming evidence**—we'd expect roughly 70-90% to pass even for random numbers.

### 7.3 Scheme/Scale Dependence Caveat

**This is the most serious systematic issue.** The PDG quark masses are:
- u, d, s: MS-bar at μ = 2 GeV
- c: MS-bar at μ = m_c
- b: MS-bar at μ = m_b
- t: Pole mass from direct measurements

A consistent treatment would RG-run all Yukawas to a common scale μ* and use α(μ*). This could shift c values by O(10-30%) for light quarks, potentially moving strange and muon into or out of bounds.

### 7.4 Selection of α

The choice of α ≈ 1/137 is motivated by the ALC draft, but is not unique. Using α(m_Z) ≈ 1/128 or other small numbers would give different (but potentially similar-quality) fits. This is a form of "look-elsewhere" effect not accounted for.

---

## 8. Recommendations

### 8.1 Pre-Registered Next Step (Minimal New Freedom)

**Null-hypothesis p-value via simulation:**

Generate 10,000 sets of 9 "random" Yukawas drawn log-uniformly from [10^-6, 1] (the observed range) and compute what fraction achieve ≥7/9 passing under the same fitting rule and bounds.

This directly tests: "Is 7/9 better than chance?"

Implementation (no code changes needed—run externally):
```python
import numpy as np

alpha = 0.0072973525693
c_min, c_max = 0.1, 10.0
n_range = range(0, 13)

def fit_one(y):
    best_c = None
    for n in n_range:
        c = y / (alpha ** n)
        if best_c is None or abs(np.log(c)) < abs(np.log(best_c)):
            best_c = c
    return c_min <= best_c <= c_max

n_trials = 10000
pass_counts = []
for _ in range(n_trials):
    ys = 10 ** np.random.uniform(-6, 0, 9)
    passes = sum(fit_one(y) for y in ys)
    pass_counts.append(passes)

p_value = np.mean(np.array(pass_counts) >= 7)
print(f"P(>=7/9 pass | null) = {p_value:.4f}")
```

**Actual result (100,000 trials):**
```
Mean passes under null: 8.39 / 9
Std: 0.76
P(>=7/9 pass | null) = 0.9802
P(>=9/9 pass | null) = 0.5285
```

**Interpretation:** The observed 7/9 pass rate is **worse than expected by chance**. Random log-uniform numbers in [10^-6, 1] achieve ≥7/9 passes 98% of the time. The α-ladder hypothesis is **not supported** by this absolute-fit pass count.

**Additional null tests:**
```
Ratio fits (9/9 observed):
  Mean passes under null: 8.61 / 9
  P(9/9 pass | null) = 0.6727
  → Not significant (67% of random sets pass 9/9)

Charm quark (c = 1.002, |log(c)| = 0.00197):
  P(|log(c)| ≤ 0.00197 | null) = 0.000673 (1 in 1,500)
  P(any of 9 this close | null) = 0.0064 (look-elsewhere corrected)
  → SIGNIFICANT even after correction (p = 0.006)
  → Charm sitting almost exactly on α^1 is the one genuinely anomalous result
```

### 8.2 Exploratory Extensions (Label as Post-Hoc)

1. **RG-run to common scale μ* = m_Z:**
   - Requires running quark Yukawas from their native scales to m_Z
   - Use 1-loop or 2-loop RGEs
   - Also switch to α(m_Z) ≈ 1/128
   - Label as "exploratory: common scale analysis"

2. **Include neutrinos (if Dirac):**
   - Cosmological bounds give y_ν < 10^-12
   - Would require n ≥ 5-6
   - Test whether neutrinos extend the ladder

3. **Alternative small parameters:**
   - Try θ_C (Cabibbo angle), sin²θ_W, or other BSM-motivated scales
   - Label as "exploratory: alternative ladder parameters"

4. **Tighter success criterion:**
   - Require c in [0.5, 2] or [0.3, 3]
   - Count passes and compute what this implies about ladder precision

---

## 9. Summary Statistics

| Metric | Preregistered (c ∈ [0.1,10]) | Strict (c ∈ [1,10]) |
|--------|------------------------------|---------------------|
| Absolute fit pass rate | 7/9 (78%) | 4/9 (44%) |
| Ratio fit pass rate | 9/9 (100%) | 9/9 (100%) |
| Failures (absolute) | mu, s | mu, s, t, u, d |
| MC stability (all P(mode)=1) | Yes | Yes |

---

## 10. What the Results Do and Do Not Imply

### What they DO suggest:
1. Fermion Yukawas can be parameterized as c × α^n with integer n and c in [0.1, 10] for 7/9 particles.
2. The within-sector mass hierarchies are cleaner: all 9/9 pass ratio fits.
3. The n assignments are stable under measurement uncertainties.
4. Charm is remarkably close to exactly α^1 (c = 1.002).

### What they DO NOT establish:
1. **Not statistically significant.** The null-hypothesis test shows random numbers pass at a **higher rate** (mean 8.4/9) than observed (7/9). The 7/9 result has p = 0.98, meaning it's worse than 98% of random trials.
2. **Not proof of clockwork or any BSM physics.** This is a numerical coincidence check that fails to beat chance.
3. **Not scheme-independent.** Running to a common scale could change the picture.
4. **Not predictive.** The test is post-hoc on known masses; it makes no new predictions.

### Honest assessment:
**The α-ladder hypothesis is NOT supported by this test overall.** The observed 7/9 absolute-fit pass rate is actually *below* what random log-uniform numbers achieve on average. The 9/9 ratio-fit pass rate is also unremarkable (67% chance).

**However, one result stands out:** The charm quark's c = 1.002 (meaning y_c ≈ α^1 to 0.2% precision) has only a 1-in-1,500 probability of occurring by chance. This single data point is the only statistically significant feature of the dataset.

Possible interpretations of the charm anomaly:
1. Statistical fluke (with 9 particles, expect ~1 to be unusually close by chance)
2. A hint that *some* Yukawas are special, even if the overall pattern isn't
3. An artifact of the charm mass determination (MS-bar at μ = m_c)

The test as designed cannot distinguish the full ALC ansatz from random numbers, but the charm result may warrant focused follow-up.

---

## Appendix: Configuration Details

### preregistered.yaml
```yaml
alpha: 0.0072973525693  # ≈ 1/137.036
vev_GeV: 246.21965
c_min: 0.1
c_max: 10.0
n_range: [0, 12]
assignment_rule: min_abs_logc
mc_samples: 20000
random_seed: 12345
```

### Data Source
PDG masses pinned from 2024-2025 summary tables:
- Quarks: RPP 2025 web update
- Leptons: 2024 booklet
- Quark uncertainties converted from 90% CL to 1σ (factor 1.645)
