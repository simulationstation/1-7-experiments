"""
Charm Robustness Study Pipeline.

This module implements the charm robustness analysis to determine whether
the y_c ≈ α coincidence is robust under scale-consistent comparisons.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .running import (
    RunningConstants,
    CharmScaleResult,
    load_running_constants,
    run_test_A,
    run_test_B,
    compute_null_pvalue_test_A,
    compute_null_pvalue_test_B,
    compute_null_pvalue_simple,
    yukawa_from_mass,
)


def run_charm_study(cfg_path: Path, output_dir_override: Path | None = None) -> None:
    """
    Run the charm robustness study.

    This implements:
    - TEST A: Alpha-at-scale substitution (fixed y_c, varying alpha)
    - TEST B: Common-scale running (both y_c and alpha at same scale)
    - Null-calibrated p-values with look-elsewhere control
    """
    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))

    out_dir = Path(output_dir_override) if output_dir_override else Path(cfg["analysis"]["output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load constants
    consts_path = Path(cfg["data"]["running_constants"])
    consts = load_running_constants(consts_path)

    vev_GeV = float(cfg["constants"]["vev_GeV"])
    n_trials = int(cfg.get("monte_carlo", {}).get("n_trials", 100000))
    seed = int(cfg["analysis"]["random_seed"])

    # Compute baseline charm Yukawa
    y_c_baseline = yukawa_from_mass(consts.m_c_GeV, vev_GeV)
    alpha_0 = consts.alpha_0
    c_baseline = y_c_baseline / alpha_0
    d_baseline = abs(math.log(c_baseline))

    print("=" * 60)
    print("CHARM ROBUSTNESS STUDY")
    print("=" * 60)
    print(f"\nBaseline (from preregistered analysis):")
    print(f"  m_c(m_c) = {consts.m_c_GeV:.4f} GeV")
    print(f"  y_c = {y_c_baseline:.7f}")
    print(f"  alpha(0) = {alpha_0:.10f}")
    print(f"  c = y_c / alpha(0) = {c_baseline:.6f}")
    print(f"  d = |ln(c)| = {d_baseline:.6f}")

    # =====================
    # TEST A: Alpha-at-scale substitution
    # =====================
    print("\n" + "=" * 60)
    print("TEST A: Alpha-at-scale substitution")
    print("(Keep y_c fixed, vary alpha)")
    print("=" * 60)

    results_A = run_test_A(y_c_baseline, consts)

    print("\nResults:")
    print(f"{'Scale':<25} {'alpha':>12} {'c = y_c/alpha':>14} {'d = |ln(c)|':>12}")
    print("-" * 65)
    for r in results_A:
        print(f"{r.scale_name:<25} {r.alpha_value:>12.10f} {r.c_coeff:>14.6f} {r.d_metric:>12.6f}")

    # Compute D_A = min d across scales
    D_A = min(r.d_metric for r in results_A)
    best_A = min(results_A, key=lambda r: r.d_metric)
    print(f"\nD_A = min(d) = {D_A:.6f} (at {best_A.scale_name})")

    # Null p-value for TEST A
    print(f"\nComputing null p-value (n={n_trials} trials)...")
    p_A, mean_D_A = compute_null_pvalue_test_A(D_A, consts, n_trials, seed)
    print(f"  P(D_A <= {D_A:.6f} | null) = {p_A:.6f}")
    print(f"  Mean D_A under null = {mean_D_A:.4f}")

    # Also compute single-alpha p-values for comparison
    print("\nSingle-alpha p-values (no look-elsewhere):")
    p_alpha0, _ = compute_null_pvalue_simple(d_baseline, consts.alpha_0, n_trials, seed)
    p_alpha_tau, _ = compute_null_pvalue_simple(
        results_A[1].d_metric, consts.alpha_m_tau, n_trials, seed
    )
    p_alpha_Z, _ = compute_null_pvalue_simple(
        results_A[2].d_metric, consts.alpha_m_Z, n_trials, seed
    )
    print(f"  alpha(0):     p = {p_alpha0:.6f}")
    print(f"  alpha(m_tau): p = {p_alpha_tau:.6f}")
    print(f"  alpha(m_Z):   p = {p_alpha_Z:.6f}")

    # =====================
    # TEST B: Common-scale running
    # =====================
    print("\n" + "=" * 60)
    print("TEST B: Common-scale running")
    print("(Run m_c to each scale, use alpha at that scale)")
    print("=" * 60)

    results_B = run_test_B(vev_GeV, consts)

    print("\nResults:")
    print(f"{'Scale':<30} {'m_c (GeV)':>10} {'y_c':>12} {'alpha':>12} {'c':>10} {'d':>10}")
    print("-" * 90)
    for r in results_B:
        print(f"{r.scale_name:<30} {r.m_c_GeV:>10.5f} {r.y_c:>12.7f} {r.alpha_value:>12.10f} {r.c_coeff:>10.5f} {r.d_metric:>10.6f}")

    # Compute D_B = min d across non-baseline scales
    results_B_nonbaseline = [r for r in results_B if "baseline" not in r.scale_name]
    D_B = min(r.d_metric for r in results_B_nonbaseline)
    best_B = min(results_B_nonbaseline, key=lambda r: r.d_metric)
    print(f"\nD_B = min(d) = {D_B:.6f} (at {best_B.scale_name})")

    # Null p-value for TEST B
    print(f"\nComputing null p-value (n={n_trials} trials)...")
    p_B, mean_D_B = compute_null_pvalue_test_B(D_B, consts, vev_GeV, n_trials, seed)
    print(f"  P(D_B <= {D_B:.6f} | null) = {p_B:.6f}")
    print(f"  Mean D_B under null = {mean_D_B:.4f}")

    # =====================
    # Summary comparison
    # =====================
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    print("\nComparison of charm's 'specialness' across analyses:")
    print(f"{'Analysis':<40} {'d = |ln(c)|':>12} {'p-value':>12}")
    print("-" * 66)
    print(f"{'Baseline: alpha(0) only':<40} {d_baseline:>12.6f} {p_alpha0:>12.6f}")
    print(f"{'TEST A: best alpha (look-elsewhere)':<40} {D_A:>12.6f} {p_A:>12.6f}")
    print(f"{'TEST B: common-scale (look-elsewhere)':<40} {D_B:>12.6f} {p_B:>12.6f}")

    # Interpretation
    print("\n" + "=" * 60)
    print("INTERPRETATION")
    print("=" * 60)

    if D_B < 0.1:
        print("\n>>> Charm remains close to an alpha ladder rung even at common scales.")
        print(f"    At best common scale, c = {best_B.c_coeff:.4f}, only {100*abs(best_B.c_coeff - 1):.1f}% from 1.")
    elif D_B < 0.5:
        print("\n>>> Charm is moderately close at common scales.")
        print(f"    At best common scale, c = {best_B.c_coeff:.4f}")
    else:
        print("\n>>> Charm's 'exactness' largely evaporates at common scales.")
        print(f"    At best common scale, c = {best_B.c_coeff:.4f}, which is {100*abs(best_B.c_coeff - 1):.0f}% from 1.")

    if p_B < 0.01:
        print(f"\n    The p-value {p_B:.4f} suggests this is still statistically unusual.")
    elif p_B < 0.05:
        print(f"\n    The p-value {p_B:.4f} is marginally significant (p < 0.05).")
    else:
        print(f"\n    The p-value {p_B:.4f} indicates this is NOT statistically significant.")
        print("    The coincidence is likely a scale/scheme artifact.")

    # =====================
    # Write artifacts
    # =====================

    # CSV
    df_A = pd.DataFrame([{
        "test": "A",
        "scale_name": r.scale_name,
        "mu_GeV": r.mu_GeV,
        "alpha": r.alpha_value,
        "m_c_GeV": r.m_c_GeV,
        "y_c": r.y_c,
        "c": r.c_coeff,
        "d": r.d_metric,
    } for r in results_A])

    df_B = pd.DataFrame([{
        "test": "B",
        "scale_name": r.scale_name,
        "mu_GeV": r.mu_GeV,
        "alpha": r.alpha_value,
        "m_c_GeV": r.m_c_GeV,
        "y_c": r.y_c,
        "c": r.c_coeff,
        "d": r.d_metric,
    } for r in results_B])

    df_all = pd.concat([df_A, df_B], ignore_index=True)
    df_all.to_csv(out_dir / "charm_robustness.csv", index=False)
    print(f"\nWrote: {out_dir / 'charm_robustness.csv'}")

    # JSON
    results_json = {
        "metadata": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "config_path": str(cfg_path),
        },
        "baseline": {
            "m_c_GeV": consts.m_c_GeV,
            "y_c": y_c_baseline,
            "alpha_0": alpha_0,
            "c": c_baseline,
            "d": d_baseline,
        },
        "test_A": {
            "description": "Alpha-at-scale substitution (fixed y_c)",
            "results": [asdict(r) for r in results_A],
            "D_A": D_A,
            "best_scale": best_A.scale_name,
            "p_value": p_A,
            "p_value_single_alpha0": p_alpha0,
        },
        "test_B": {
            "description": "Common-scale running (y_c and alpha at same mu)",
            "results": [asdict(r) for r in results_B],
            "D_B": D_B,
            "best_scale": best_B.scale_name,
            "p_value": p_B,
        },
        "monte_carlo": {
            "n_trials": n_trials,
            "seed": seed,
        },
    }
    (out_dir / "charm_robustness.json").write_text(
        json.dumps(results_json, indent=2, default=str),
        encoding="utf-8"
    )
    print(f"Wrote: {out_dir / 'charm_robustness.json'}")

    # Plots
    _make_charm_plots(results_A, results_B, d_baseline, out_dir)
    print(f"Wrote: {out_dir / 'charm_scale_plots.png'}")

    # Markdown report
    _write_markdown_report(
        results_A, results_B,
        d_baseline, D_A, D_B,
        p_alpha0, p_A, p_B,
        consts, vev_GeV, n_trials,
        out_dir,
    )
    print(f"Wrote: {out_dir.parent / 'report_charm_robustness.md'}")


def _make_charm_plots(
    results_A: List[CharmScaleResult],
    results_B: List[CharmScaleResult],
    d_baseline: float,
    out_dir: Path,
) -> None:
    """Generate plots for charm robustness study."""

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Plot 1: c coefficient vs scale (Test A)
    ax1 = axes[0, 0]
    scales_A = [r.scale_name.replace("alpha(", "").replace(")", "") for r in results_A]
    c_vals_A = [r.c_coeff for r in results_A]
    ax1.bar(scales_A, c_vals_A, color=['blue', 'orange', 'green'])
    ax1.axhline(y=1.0, color='red', linestyle='--', label='c = 1 (exact)')
    ax1.set_ylabel('c = y_c / alpha')
    ax1.set_title('TEST A: c vs alpha scale (fixed y_c)')
    ax1.legend()
    ax1.set_ylim(0.9, max(c_vals_A) * 1.1)

    # Plot 2: d metric vs scale (Test A)
    ax2 = axes[0, 1]
    d_vals_A = [r.d_metric for r in results_A]
    ax2.bar(scales_A, d_vals_A, color=['blue', 'orange', 'green'])
    ax2.set_ylabel('d = |ln(c)|')
    ax2.set_title('TEST A: d vs alpha scale (lower is closer to c=1)')

    # Plot 3: c coefficient vs scale (Test B)
    ax3 = axes[1, 0]
    scales_B = ['baseline\nm_c(m_c),α(0)', 'm_tau\nm_c(m_τ),α(m_τ)', 'm_Z\nm_c(m_Z),α(m_Z)']
    c_vals_B = [r.c_coeff for r in results_B]
    colors_B = ['gray', 'orange', 'green']
    ax3.bar(scales_B, c_vals_B, color=colors_B)
    ax3.axhline(y=1.0, color='red', linestyle='--', label='c = 1 (exact)')
    ax3.set_ylabel('c = y_c(μ) / alpha(μ)')
    ax3.set_title('TEST B: c vs common scale')
    ax3.legend()

    # Plot 4: d metric vs scale (Test B)
    ax4 = axes[1, 1]
    d_vals_B = [r.d_metric for r in results_B]
    ax4.bar(scales_B, d_vals_B, color=colors_B)
    ax4.set_ylabel('d = |ln(c)|')
    ax4.set_title('TEST B: d vs common scale')

    plt.tight_layout()
    plt.savefig(out_dir / "charm_scale_plots.png", dpi=150)
    plt.close()


def _write_markdown_report(
    results_A: List[CharmScaleResult],
    results_B: List[CharmScaleResult],
    d_baseline: float,
    D_A: float,
    D_B: float,
    p_alpha0: float,
    p_A: float,
    p_B: float,
    consts: RunningConstants,
    vev_GeV: float,
    n_trials: int,
    out_dir: Path,
) -> None:
    """Write the charm robustness markdown report."""

    # Find best results
    best_A = min(results_A, key=lambda r: r.d_metric)
    results_B_nonbaseline = [r for r in results_B if "baseline" not in r.scale_name]
    best_B = min(results_B_nonbaseline, key=lambda r: r.d_metric)

    y_c_baseline = results_A[0].y_c
    c_baseline = results_A[0].c_coeff

    report = f"""# Charm Robustness Study: Report

**Generated:** {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}

---

## Executive Summary

| Question | Answer |
|----------|--------|
| Does charm remain "exceptionally exact" after scale-consistent α(μ)? | **{"Yes" if D_A < 0.05 else "Partially" if D_A < 0.2 else "No"}** (D_A = {D_A:.4f}) |
| Does it survive running y_c to m_Z (Test B)? | **{"Yes" if D_B < 0.1 else "Partially" if D_B < 0.5 else "No"}** (D_B = {D_B:.4f}) |
| Updated p-value for Test A | **{p_A:.4f}** (vs {p_alpha0:.4f} single-alpha) |
| Updated p-value for Test B | **{p_B:.4f}** |

**Bottom line:** {"The charm coincidence is robust under scale-consistent analysis." if p_B < 0.05 else "The charm coincidence largely evaporates at common scales and is NOT statistically significant."}

---

## 1. Baseline Reference

From the preregistered analysis (`config/preregistered.yaml`):

| Quantity | Value |
|----------|-------|
| m_c(m_c) | {consts.m_c_GeV:.4f} GeV |
| v (vev) | {vev_GeV:.2f} GeV |
| y_c = √2 m_c / v | {y_c_baseline:.7f} |
| α(0) | {consts.alpha_0:.10f} |
| c = y_c / α(0) | **{c_baseline:.6f}** |
| d = |ln(c)| | **{d_baseline:.6f}** |

The baseline c ≈ 1.002 means y_c is within 0.2% of α.

---

## 2. TEST A: Alpha-at-Scale Substitution

**Procedure:** Keep y_c fixed at its baseline value. Compute c = y_c / α(μ) for different α values.

**Rationale:** Tests whether charm's "exactness" depends on which α we use.

### Results

| Scale | α(μ) | c = y_c / α | d = |ln(c)| |
|-------|------|-------------|------------|
| α(0) | {consts.alpha_0:.10f} | {results_A[0].c_coeff:.6f} | {results_A[0].d_metric:.6f} |
| α(m_τ) | {consts.alpha_m_tau:.10f} | {results_A[1].c_coeff:.6f} | {results_A[1].d_metric:.6f} |
| α(m_Z) | {consts.alpha_m_Z:.10f} | {results_A[2].c_coeff:.6f} | {results_A[2].d_metric:.6f} |

**D_A = min(d) = {D_A:.6f}** (achieved at {best_A.scale_name})

### Statistical Significance

- Single-alpha p-value (α(0) only): **{p_alpha0:.6f}**
- Look-elsewhere corrected (3 α choices): **{p_A:.6f}**

The correction for testing 3 α values {"does not substantially change" if abs(p_A - p_alpha0) < 0.01 else "increases"} the p-value.

---

## 3. TEST B: Common-Scale Running

**Procedure:** Run the charm mass m_c from its definition scale to μ = m_τ and μ = m_Z using 1-loop QCD. Use α(μ) at the same scale.

**Rationale:** The physically meaningful comparison is y_c(μ) vs α(μ) at the same renormalization scale.

### QCD Running Details

- Running formula: m(μ₂)/m(μ₁) = [α_s(μ₂)/α_s(μ₁)]^(γ_m/β₀)
- γ_m = 4 (1-loop mass anomalous dimension)
- β₀ = 11 - 2n_f/3 (1-loop QCD beta function)
- Used α_s(m_Z) = {consts.alpha_s_m_Z}, α_s(m_τ) = {consts.alpha_s_m_tau}

### Results

| Scale | m_c(μ) [GeV] | y_c(μ) | α(μ) | c | d = |ln(c)| |
|-------|--------------|--------|------|---|------------|
| Baseline | {results_B[0].m_c_GeV:.5f} | {results_B[0].y_c:.7f} | {results_B[0].alpha_value:.10f} | {results_B[0].c_coeff:.5f} | {results_B[0].d_metric:.6f} |
| μ = m_τ | {results_B[1].m_c_GeV:.5f} | {results_B[1].y_c:.7f} | {results_B[1].alpha_value:.10f} | {results_B[1].c_coeff:.5f} | {results_B[1].d_metric:.6f} |
| μ = m_Z | {results_B[2].m_c_GeV:.5f} | {results_B[2].y_c:.7f} | {results_B[2].alpha_value:.10f} | {results_B[2].c_coeff:.5f} | {results_B[2].d_metric:.6f} |

**D_B = min(d) = {D_B:.6f}** (achieved at {best_B.scale_name})

### Statistical Significance

- Look-elsewhere corrected p-value: **{p_B:.6f}**

---

## 4. Interpretation

### What the numbers mean

| Metric | Baseline | Best Test A | Best Test B |
|--------|----------|-------------|-------------|
| c | {c_baseline:.4f} | {best_A.c_coeff:.4f} | {best_B.c_coeff:.4f} |
| d = |ln(c)| | {d_baseline:.4f} | {D_A:.4f} | {D_B:.4f} |
| % from c=1 | {100*abs(c_baseline - 1):.1f}% | {100*abs(best_A.c_coeff - 1):.1f}% | {100*abs(best_B.c_coeff - 1):.1f}% |

### Key finding

"""

    if p_B < 0.01:
        report += f"""**The charm coincidence IS robust.**

Even at common renormalization scales (Test B), charm's coefficient c = {best_B.c_coeff:.3f} remains close to 1, with p = {p_B:.4f}. This is statistically significant and unlikely to be a scale/scheme artifact.

The y_c ≈ α relationship appears to be a genuine numerical feature of the Standard Model."""

    elif p_B < 0.05:
        report += f"""**The charm coincidence is marginally robust.**

At common scales, c = {best_B.c_coeff:.3f} with p = {p_B:.4f}. This is marginally significant (p < 0.05) but not compelling. The coincidence may be partially due to scale/scheme choices."""

    else:
        report += f"""**The charm coincidence is NOT robust.**

At common renormalization scales (Test B), the coefficient c = {best_B.c_coeff:.3f}, which is {100*abs(best_B.c_coeff - 1):.0f}% from 1. The p-value {p_B:.4f} indicates this is NOT statistically significant.

**Conclusion:** The apparent y_c ≈ α coincidence seen at baseline (c = {c_baseline:.4f}) is a scale/scheme artifact. When both quantities are evaluated at compatible scales, the "exactness" evaporates."""

    report += f"""

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
| α(0) | {consts.alpha_0:.10f} | CODATA 2018 via PDG |
| α(m_τ) | {consts.alpha_m_tau:.7f} | PDG 2024 Electroweak Review |
| α(m_Z) | {consts.alpha_m_Z:.7f} | PDG 2024 Electroweak Review |
| α_s(m_Z) | {consts.alpha_s_m_Z:.4f} | PDG 2024 QCD Review |
| α_s(m_τ) | {consts.alpha_s_m_tau:.3f} | PDG 2024 QCD Review |

---

## 7. Artifacts Generated

- `artifacts/charm_robustness.csv` - All numerical results
- `artifacts/charm_robustness.json` - Structured JSON output
- `artifacts/charm_scale_plots.png` - Visualization

---

## Appendix: Monte Carlo Details

- Null hypothesis: y drawn log-uniform from [10⁻⁶, 1]
- Number of trials: {n_trials:,}
- Test A statistic: D_A = min over {{α(0), α(m_τ), α(m_Z)}} of |ln(y/α)|
- Test B statistic: D_B = min over {{m_τ, m_Z}} of |ln(y(μ)/α(μ))|
"""

    (out_dir.parent / "report_charm_robustness.md").write_text(report, encoding="utf-8")
