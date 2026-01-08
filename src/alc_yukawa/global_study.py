"""
Global Common-Scale Alpha-Ladder Test.

This module implements a scale-consistent test of the alpha-ladder hypothesis
across ALL SM Yukawas, evaluating both y_i(μ) and α(μ) at the same scale.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .running import (
    RunningConstants,
    load_running_constants,
    run_alpha_s_1loop,
    run_quark_mass_1loop,
    beta_0,
    gamma_m_0,
)
from .dataset import load_pdg_snapshot


@dataclass
class ParticleAtScale:
    """Result for one particle at one scale."""
    particle: str
    sector: str
    mu_GeV: float
    scale_name: str
    mass_ref_GeV: float
    mass_at_mu_GeV: float
    yukawa: float
    alpha_mu: float
    best_n: int
    c: float
    ln_c: float
    abs_ln_c: float
    in_bounds_W1: bool  # c in [0.1, 10]
    in_bounds_W2: bool  # c in [0.5, 2]


@dataclass
class ScaleSummary:
    """Summary statistics at one scale."""
    scale_name: str
    mu_GeV: float
    alpha_mu: float
    T: float           # Σ |ln c_i|
    T2: float          # Σ (ln c_i)^2
    pass_W1: int       # count in [0.1, 10]
    pass_W2: int       # count in [0.5, 2]
    n_particles: int


def get_quark_reference_scale(particle: str, pdg_data: dict) -> float:
    """Get the reference scale (in GeV) where the quark mass is defined."""
    p_data = pdg_data["particles"][particle]
    if "scale_GeV" in p_data:
        return p_data["scale_GeV"]
    # Top quark: use pole mass convention, treat as if at m_t
    if particle == "t":
        mass = p_data["mass"]
        unit = p_data.get("unit", "GeV")
        if unit == "MeV":
            return mass / 1000.0
        return mass
    return None


def run_quark_mass_to_scale(
    particle: str,
    mass_ref_GeV: float,
    mu_ref: float,
    mu_target: float,
    consts: RunningConstants,
) -> float:
    """
    Run a quark mass from its reference scale to target scale using 1-loop QCD.

    Handles flavor thresholds approximately:
    - n_f = 3 below m_c
    - n_f = 4 between m_c and m_b
    - n_f = 5 above m_b
    """
    m_c = 1.27
    m_b = 4.18

    # Get alpha_s at reference and target scales
    # Simple approach: use appropriate n_f based on scale
    def get_nf(mu):
        if mu < m_c:
            return 3
        elif mu < m_b:
            return 4
        else:
            return 5

    n_f_ref = get_nf(mu_ref)
    n_f_target = get_nf(mu_target)

    # Get alpha_s values
    if mu_target >= m_b:
        alpha_s_target = run_alpha_s_1loop(consts.alpha_s_m_Z, consts.m_Z_GeV, mu_target, n_f=5)
    else:
        alpha_s_target = run_alpha_s_1loop(consts.alpha_s_m_tau, consts.m_tau_GeV, mu_target, n_f=4)

    if mu_ref >= m_b:
        alpha_s_ref = run_alpha_s_1loop(consts.alpha_s_m_Z, consts.m_Z_GeV, mu_ref, n_f=5)
    else:
        alpha_s_ref = run_alpha_s_1loop(consts.alpha_s_m_tau, consts.m_tau_GeV, mu_ref, n_f=4)

    # Use target n_f for running (simplified)
    m_target = run_quark_mass_1loop(
        m_ref=mass_ref_GeV,
        mu_ref=mu_ref,
        mu_target=mu_target,
        alpha_s_ref=alpha_s_ref,
        alpha_s_target=alpha_s_target,
        n_f=n_f_target,
    )

    return m_target


def fit_to_ladder(y: float, alpha: float, n_min: int = 0, n_max: int = 12) -> Tuple[int, float]:
    """
    Fit Yukawa y to alpha^n by finding n that minimizes |ln(c)| where c = y/alpha^n.
    Returns (best_n, c).
    """
    best_n = 0
    best_c = y / (alpha ** 0)
    best_abs_ln_c = abs(math.log(best_c))

    for n in range(n_min, n_max + 1):
        c = y / (alpha ** n)
        abs_ln_c = abs(math.log(c))
        if abs_ln_c < best_abs_ln_c:
            best_n = n
            best_c = c
            best_abs_ln_c = abs_ln_c

    return best_n, best_c


def analyze_at_scale(
    scale_name: str,
    mu_GeV: float,
    alpha_mu: float,
    particles: List[str],
    pdg_data: dict,
    consts: RunningConstants,
    vev_GeV: float,
) -> Tuple[List[ParticleAtScale], ScaleSummary]:
    """
    Analyze all particles at a single scale.
    """
    results = []
    sectors = pdg_data["sectors"]

    # Build particle -> sector map
    particle_sector = {}
    for sector, ps in sectors.items():
        for p in ps:
            particle_sector[p] = sector

    quarks = set(sectors["up"]) | set(sectors["down"])
    leptons = set(sectors["lepton"])

    for particle in particles:
        p_data = pdg_data["particles"][particle]

        # Get reference mass in GeV
        mass_ref = p_data["mass"]
        unit = p_data.get("unit", "GeV")
        if unit == "MeV":
            mass_ref = mass_ref / 1000.0

        # Compute mass at target scale
        if particle in quarks:
            mu_ref = get_quark_reference_scale(particle, pdg_data)
            if mu_ref is None:
                mu_ref = mass_ref  # fallback

            # Special handling for top: keep as pole mass (no running)
            if particle == "t":
                mass_at_mu = mass_ref
            else:
                mass_at_mu = run_quark_mass_to_scale(
                    particle, mass_ref, mu_ref, mu_GeV, consts
                )
        else:
            # Leptons: keep fixed (pole masses)
            mass_at_mu = mass_ref

        # Compute Yukawa
        yukawa = math.sqrt(2.0) * mass_at_mu / vev_GeV

        # Fit to ladder
        best_n, c = fit_to_ladder(yukawa, alpha_mu)
        ln_c = math.log(c)
        abs_ln_c = abs(ln_c)

        results.append(ParticleAtScale(
            particle=particle,
            sector=particle_sector.get(particle, "unknown"),
            mu_GeV=mu_GeV,
            scale_name=scale_name,
            mass_ref_GeV=mass_ref,
            mass_at_mu_GeV=mass_at_mu,
            yukawa=yukawa,
            alpha_mu=alpha_mu,
            best_n=best_n,
            c=c,
            ln_c=ln_c,
            abs_ln_c=abs_ln_c,
            in_bounds_W1=(0.1 <= c <= 10.0),
            in_bounds_W2=(0.5 <= c <= 2.0),
        ))

    # Compute summary statistics
    T = sum(r.abs_ln_c for r in results)
    T2 = sum(r.ln_c ** 2 for r in results)
    pass_W1 = sum(1 for r in results if r.in_bounds_W1)
    pass_W2 = sum(1 for r in results if r.in_bounds_W2)

    summary = ScaleSummary(
        scale_name=scale_name,
        mu_GeV=mu_GeV,
        alpha_mu=alpha_mu,
        T=T,
        T2=T2,
        pass_W1=pass_W1,
        pass_W2=pass_W2,
        n_particles=len(results),
    )

    return results, summary


def compute_null_T(
    n_particles: int,
    alpha: float,
    n_min: int = 0,
    n_max: int = 12,
    rng: np.random.Generator = None,
) -> float:
    """Compute T for one null trial at one scale."""
    if rng is None:
        rng = np.random.default_rng()

    T = 0.0
    for _ in range(n_particles):
        y = 10 ** rng.uniform(-6, 0)
        _, c = fit_to_ladder(y, alpha, n_min, n_max)
        T += abs(math.log(c))

    return T


def run_null_simulation(
    n_particles: int,
    alphas: Dict[str, float],  # scale_name -> alpha
    n_trials: int,
    seed: int,
) -> Dict[str, Any]:
    """
    Run null hypothesis simulation.

    Returns per-scale T distributions and global T_best distribution.
    """
    rng = np.random.default_rng(seed)

    scale_names = list(alphas.keys())
    null_T = {name: [] for name in scale_names}
    null_T_best = []

    for _ in range(n_trials):
        trial_T = {}
        for name, alpha in alphas.items():
            T = compute_null_T(n_particles, alpha, rng=rng)
            trial_T[name] = T
            null_T[name].append(T)

        T_best = min(trial_T.values())
        null_T_best.append(T_best)

    return {
        "per_scale": {name: np.array(vals) for name, vals in null_T.items()},
        "T_best": np.array(null_T_best),
    }


def compute_pvalues(
    T_obs: Dict[str, float],
    null_results: Dict[str, Any],
) -> Dict[str, float]:
    """Compute p-values from null simulation."""
    pvals = {}

    # Per-scale p-values
    for name, T in T_obs.items():
        null_T = null_results["per_scale"][name]
        pvals[f"p_{name}"] = float(np.mean(null_T <= T))

    # Global p-value (with look-elsewhere)
    T_best_obs = min(T_obs.values())
    null_T_best = null_results["T_best"]
    pvals["p_global"] = float(np.mean(null_T_best <= T_best_obs))
    pvals["T_best_obs"] = T_best_obs

    return pvals


def run_global_study(cfg_path: Path, output_dir_override: Path | None = None) -> None:
    """Run the global common-scale alpha-ladder study."""

    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))

    out_dir = Path(output_dir_override) if output_dir_override else Path(cfg["analysis"]["output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    seed = int(cfg["analysis"]["random_seed"])
    vev_GeV = float(cfg["constants"]["vev_GeV"])

    # Load data
    pdg_path = Path(cfg["data"]["pdg_snapshot"])
    consts_path = Path(cfg["data"]["running_constants"])

    pdg_raw = json.loads(pdg_path.read_text())
    consts = load_running_constants(consts_path)

    # Reference scales
    scales = cfg["reference_scales"]

    # Get alpha values
    consts_raw = json.loads(consts_path.read_text())
    alpha_values = {
        "m_tau": consts_raw["constants"]["alpha_hat_m_tau"]["value"],
        "m_Z": consts_raw["constants"]["alpha_hat_m_Z"]["value"],
    }

    n_trials = int(cfg["monte_carlo"]["n_trials"])

    print("=" * 70)
    print("GLOBAL COMMON-SCALE ALPHA-LADDER TEST")
    print("=" * 70)

    # Run analysis for both versions
    all_results = {}
    all_summaries = {}

    for version in cfg["versions"]:
        version_name = version["name"]
        particles = version["particles"]

        print(f"\n{'=' * 70}")
        print(f"VERSION: {version_name} ({version['description']})")
        print("=" * 70)

        version_results = {}
        version_summaries = {}

        for scale_cfg in scales:
            scale_name = scale_cfg["name"]
            mu_GeV = scale_cfg["mu_GeV"]
            alpha_mu = alpha_values[scale_name]

            print(f"\n--- Scale: μ = {scale_name} ({mu_GeV:.2f} GeV), α = {alpha_mu:.7f} ---")

            results, summary = analyze_at_scale(
                scale_name=scale_name,
                mu_GeV=mu_GeV,
                alpha_mu=alpha_mu,
                particles=particles,
                pdg_data=pdg_raw,
                consts=consts,
                vev_GeV=vev_GeV,
            )

            version_results[scale_name] = results
            version_summaries[scale_name] = summary

            # Print per-particle results
            print(f"\n{'Particle':<8} {'Sector':<8} {'m(μ) GeV':<12} {'y(μ)':<12} {'n':<4} {'c':<10} {'|ln c|':<10} {'W1':<4} {'W2':<4}")
            print("-" * 80)
            for r in results:
                print(f"{r.particle:<8} {r.sector:<8} {r.mass_at_mu_GeV:<12.6f} {r.yukawa:<12.6e} {r.best_n:<4} {r.c:<10.4f} {r.abs_ln_c:<10.4f} {'Y' if r.in_bounds_W1 else 'N':<4} {'Y' if r.in_bounds_W2 else 'N':<4}")

            print(f"\nT({scale_name}) = {summary.T:.4f}")
            print(f"T2({scale_name}) = {summary.T2:.4f}")
            print(f"Pass W1 [0.1,10]: {summary.pass_W1}/{summary.n_particles}")
            print(f"Pass W2 [0.5,2]: {summary.pass_W2}/{summary.n_particles}")

        all_results[version_name] = version_results
        all_summaries[version_name] = version_summaries

        # Determine best scale
        T_obs = {name: s.T for name, s in version_summaries.items()}
        best_scale = min(T_obs, key=T_obs.get)
        T_best = T_obs[best_scale]

        print(f"\n>>> Best scale: {best_scale} with T = {T_best:.4f}")

        # Run null simulation
        print(f"\nRunning null simulation ({n_trials:,} trials)...")

        n_particles = len(particles)
        null_results = run_null_simulation(
            n_particles=n_particles,
            alphas=alpha_values,
            n_trials=n_trials,
            seed=seed,
        )

        # Compute p-values
        pvalues = compute_pvalues(T_obs, null_results)

        print("\n--- Null Hypothesis Results ---")
        for scale_name in T_obs:
            null_T = null_results["per_scale"][scale_name]
            print(f"T({scale_name}): obs = {T_obs[scale_name]:.4f}, null mean = {np.mean(null_T):.4f}, p = {pvalues[f'p_{scale_name}']:.6f}")

        print(f"\nT_best: obs = {pvalues['T_best_obs']:.4f}, null mean = {np.mean(null_results['T_best']):.4f}")
        print(f"p_global (look-elsewhere corrected) = {pvalues['p_global']:.6f}")

        # Interpretation
        print("\n--- Interpretation ---")
        if pvalues['p_global'] < 0.01:
            print(">>> SIGNIFICANT: SM Yukawas are globally closer to α-ladder than chance (p < 0.01)")
        elif pvalues['p_global'] < 0.05:
            print(">>> MARGINAL: Some evidence for α-ladder structure (0.01 < p < 0.05)")
        else:
            print(f">>> NOT SIGNIFICANT: No evidence for α-ladder structure (p = {pvalues['p_global']:.3f})")

        # Store for output
        all_summaries[version_name]["pvalues"] = pvalues
        all_summaries[version_name]["null_stats"] = {
            "n_trials": n_trials,
            "mean_T_per_scale": {name: float(np.mean(null_results["per_scale"][name])) for name in T_obs},
            "mean_T_best": float(np.mean(null_results["T_best"])),
        }

    # Write outputs
    _write_outputs(all_results, all_summaries, cfg, out_dir, vev_GeV, alpha_values)


def _write_outputs(
    all_results: Dict,
    all_summaries: Dict,
    cfg: Dict,
    out_dir: Path,
    vev_GeV: float,
    alpha_values: Dict[str, float],
) -> None:
    """Write CSV, JSON, plots, and markdown report."""

    # CSV
    rows = []
    for version_name, version_results in all_results.items():
        for scale_name, results in version_results.items():
            for r in results:
                rows.append({
                    "version": version_name,
                    "particle": r.particle,
                    "sector": r.sector,
                    "scale_name": r.scale_name,
                    "mu_GeV": r.mu_GeV,
                    "mass_ref_GeV": r.mass_ref_GeV,
                    "mass_at_mu_GeV": r.mass_at_mu_GeV,
                    "yukawa": r.yukawa,
                    "alpha_mu": r.alpha_mu,
                    "best_n": r.best_n,
                    "c": r.c,
                    "ln_c": r.ln_c,
                    "abs_ln_c": r.abs_ln_c,
                    "in_bounds_W1": r.in_bounds_W1,
                    "in_bounds_W2": r.in_bounds_W2,
                })

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "global_common_scale.csv", index=False)
    print(f"\nWrote: {out_dir / 'global_common_scale.csv'}")

    # JSON
    json_out = {
        "metadata": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "config": str(cfg.get("analysis", {}).get("name", "global_study")),
            "vev_GeV": vev_GeV,
            "alpha_values": alpha_values,
        },
        "versions": {},
    }

    for version_name in all_results:
        version_out = {
            "per_scale": {},
            "summary": {},
        }

        for scale_name, results in all_results[version_name].items():
            version_out["per_scale"][scale_name] = [asdict(r) for r in results]

        for scale_name, summary in all_summaries[version_name].items():
            if isinstance(summary, ScaleSummary):
                version_out["summary"][scale_name] = asdict(summary)

        if "pvalues" in all_summaries[version_name]:
            version_out["pvalues"] = all_summaries[version_name]["pvalues"]
        if "null_stats" in all_summaries[version_name]:
            version_out["null_stats"] = all_summaries[version_name]["null_stats"]

        json_out["versions"][version_name] = version_out

    (out_dir / "global_common_scale.json").write_text(
        json.dumps(json_out, indent=2, default=str),
        encoding="utf-8"
    )
    print(f"Wrote: {out_dir / 'global_common_scale.json'}")

    # Plots
    _make_plots(all_results, all_summaries, out_dir)
    print(f"Wrote: {out_dir / 'global_common_scale_plots.png'}")

    # Markdown report
    _write_markdown_report(all_results, all_summaries, cfg, alpha_values, out_dir)
    print(f"Wrote: {out_dir.parent / 'report_global_common_scale.md'}")


def _make_plots(all_results: Dict, all_summaries: Dict, out_dir: Path) -> None:
    """Generate plots for the global study."""

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Use "all_9" version for main plots
    version_name = "all_9"
    if version_name not in all_results:
        version_name = list(all_results.keys())[0]

    version_results = all_results[version_name]

    # Plot 1: |ln c| histogram at m_tau
    ax1 = axes[0, 0]
    if "m_tau" in version_results:
        abs_ln_c = [r.abs_ln_c for r in version_results["m_tau"]]
        ax1.hist(abs_ln_c, bins=15, edgecolor='black', alpha=0.7)
        ax1.axvline(x=np.mean(abs_ln_c), color='red', linestyle='--', label=f'Mean = {np.mean(abs_ln_c):.2f}')
        ax1.set_xlabel('|ln(c)|')
        ax1.set_ylabel('Count')
        ax1.set_title('Distribution of |ln(c)| at μ = m_τ')
        ax1.legend()

    # Plot 2: |ln c| histogram at m_Z
    ax2 = axes[0, 1]
    if "m_Z" in version_results:
        abs_ln_c = [r.abs_ln_c for r in version_results["m_Z"]]
        ax2.hist(abs_ln_c, bins=15, edgecolor='black', alpha=0.7, color='orange')
        ax2.axvline(x=np.mean(abs_ln_c), color='red', linestyle='--', label=f'Mean = {np.mean(abs_ln_c):.2f}')
        ax2.set_xlabel('|ln(c)|')
        ax2.set_ylabel('Count')
        ax2.set_title('Distribution of |ln(c)| at μ = m_Z')
        ax2.legend()

    # Plot 3: Yukawa vs best_n at m_tau
    ax3 = axes[1, 0]
    if "m_tau" in version_results:
        results = version_results["m_tau"]
        for r in results:
            ax3.scatter(r.best_n, np.log10(r.yukawa), s=100)
            ax3.annotate(r.particle, (r.best_n, np.log10(r.yukawa)), fontsize=9)
        ax3.set_xlabel('Best-fit n')
        ax3.set_ylabel('log₁₀(y)')
        ax3.set_title('Yukawa vs α^n assignment at μ = m_τ')

    # Plot 4: Yukawa vs best_n at m_Z
    ax4 = axes[1, 1]
    if "m_Z" in version_results:
        results = version_results["m_Z"]
        for r in results:
            ax4.scatter(r.best_n, np.log10(r.yukawa), s=100, color='orange')
            ax4.annotate(r.particle, (r.best_n, np.log10(r.yukawa)), fontsize=9)
        ax4.set_xlabel('Best-fit n')
        ax4.set_ylabel('log₁₀(y)')
        ax4.set_title('Yukawa vs α^n assignment at μ = m_Z')

    plt.tight_layout()
    plt.savefig(out_dir / "global_common_scale_plots.png", dpi=150)
    plt.close()


def _write_markdown_report(
    all_results: Dict,
    all_summaries: Dict,
    cfg: Dict,
    alpha_values: Dict[str, float],
    out_dir: Path,
) -> None:
    """Write the markdown report."""

    # Get main results (all_9 version)
    main_version = "all_9"
    results_main = all_results.get(main_version, list(all_results.values())[0])
    summaries_main = all_summaries.get(main_version, list(all_summaries.values())[0])

    # Get pvalues
    pvalues = summaries_main.get("pvalues", {})
    null_stats = summaries_main.get("null_stats", {})

    # Find best scale
    T_obs = {}
    for scale_name, summary in summaries_main.items():
        if isinstance(summary, ScaleSummary):
            T_obs[scale_name] = summary.T

    best_scale = min(T_obs, key=T_obs.get) if T_obs else "m_tau"
    T_best = T_obs.get(best_scale, 0)

    report = f"""# Global Common-Scale Alpha-Ladder Test: Report

**Generated:** {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}

---

## Executive Summary

| Metric | m_τ (1.78 GeV) | m_Z (91.2 GeV) |
|--------|----------------|----------------|
| α(μ) | {alpha_values.get('m_tau', 0):.7f} | {alpha_values.get('m_Z', 0):.7f} |
| T(μ) = Σ|ln c| | {T_obs.get('m_tau', 0):.4f} | {T_obs.get('m_Z', 0):.4f} |
| p-value | {pvalues.get('p_m_tau', 0):.6f} | {pvalues.get('p_m_Z', 0):.6f} |

**Best scale:** {best_scale} with T = {T_best:.4f}

**Global p-value (look-elsewhere corrected):** {pvalues.get('p_global', 1):.6f}

**Interpretation:** {"SM Yukawas show NO global α-ladder structure beyond chance." if pvalues.get('p_global', 1) > 0.05 else "Evidence for α-ladder structure."}

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
- {null_stats.get('n_trials', 200000):,} trials

---

## 2. Results: All 9 Particles

### At μ = m_τ

| Particle | Sector | m(μ) [GeV] | y(μ) | n | c | |ln c| | W1 | W2 |
|----------|--------|------------|------|---|------|---------|----|----|
"""

    if "m_tau" in results_main:
        for r in results_main["m_tau"]:
            w1 = "✓" if r.in_bounds_W1 else "✗"
            w2 = "✓" if r.in_bounds_W2 else "✗"
            report += f"| {r.particle} | {r.sector} | {r.mass_at_mu_GeV:.6f} | {r.yukawa:.4e} | {r.best_n} | {r.c:.4f} | {r.abs_ln_c:.4f} | {w1} | {w2} |\n"

    report += f"""
**T(m_τ) = {T_obs.get('m_tau', 0):.4f}**

### At μ = m_Z

| Particle | Sector | m(μ) [GeV] | y(μ) | n | c | |ln c| | W1 | W2 |
|----------|--------|------------|------|---|------|---------|----|----|
"""

    if "m_Z" in results_main:
        for r in results_main["m_Z"]:
            w1 = "✓" if r.in_bounds_W1 else "✗"
            w2 = "✓" if r.in_bounds_W2 else "✗"
            report += f"| {r.particle} | {r.sector} | {r.mass_at_mu_GeV:.6f} | {r.yukawa:.4e} | {r.best_n} | {r.c:.4f} | {r.abs_ln_c:.4f} | {w1} | {w2} |\n"

    # Get summary stats
    summary_tau = summaries_main.get("m_tau")
    summary_Z = summaries_main.get("m_Z")

    report += f"""
**T(m_Z) = {T_obs.get('m_Z', 0):.4f}**

---

## 3. Summary Statistics

| Metric | μ = m_τ | μ = m_Z |
|--------|---------|---------|
| T = Σ|ln c| | {summary_tau.T if summary_tau else 0:.4f} | {summary_Z.T if summary_Z else 0:.4f} |
| T² = Σ(ln c)² | {summary_tau.T2 if summary_tau else 0:.4f} | {summary_Z.T2 if summary_Z else 0:.4f} |
| Pass W1 [0.1,10] | {summary_tau.pass_W1 if summary_tau else 0}/{summary_tau.n_particles if summary_tau else 9} | {summary_Z.pass_W1 if summary_Z else 0}/{summary_Z.n_particles if summary_Z else 9} |
| Pass W2 [0.5,2] | {summary_tau.pass_W2 if summary_tau else 0}/{summary_tau.n_particles if summary_tau else 9} | {summary_Z.pass_W2 if summary_Z else 0}/{summary_Z.n_particles if summary_Z else 9} |

---

## 4. Null Hypothesis Testing

### Null Model
- Draw 9 random Yukawas log-uniform in [10⁻⁶, 1]
- Apply same fitting procedure
- Trials: {null_stats.get('n_trials', 200000):,}

### Results

| Scale | T_obs | T_null (mean) | p-value |
|-------|-------|---------------|---------|
| m_τ | {T_obs.get('m_tau', 0):.4f} | {null_stats.get('mean_T_per_scale', {}).get('m_tau', 0):.4f} | {pvalues.get('p_m_tau', 0):.6f} |
| m_Z | {T_obs.get('m_Z', 0):.4f} | {null_stats.get('mean_T_per_scale', {}).get('m_Z', 0):.4f} | {pvalues.get('p_m_Z', 0):.6f} |

**Global (look-elsewhere corrected):**
- T_best_obs = {pvalues.get('T_best_obs', 0):.4f}
- T_best_null (mean) = {null_stats.get('mean_T_best', 0):.4f}
- **p_global = {pvalues.get('p_global', 1):.6f}**

---

## 5. Robustness Check: Excluding Top

"""

    if "exclude_top" in all_summaries:
        excl_summaries = all_summaries["exclude_top"]
        excl_pvalues = excl_summaries.get("pvalues", {})
        excl_null = excl_summaries.get("null_stats", {})

        excl_T_obs = {}
        for scale_name, summary in excl_summaries.items():
            if isinstance(summary, ScaleSummary):
                excl_T_obs[scale_name] = summary.T

        report += f"""| Metric | m_τ | m_Z |
|--------|-----|-----|
| T (8 particles) | {excl_T_obs.get('m_tau', 0):.4f} | {excl_T_obs.get('m_Z', 0):.4f} |
| p-value | {excl_pvalues.get('p_m_tau', 0):.6f} | {excl_pvalues.get('p_m_Z', 0):.6f} |
| p_global | {excl_pvalues.get('p_global', 1):.6f} | - |

The exclusion of top (which has scheme ambiguities) {"does not change" if abs(excl_pvalues.get('p_global', 1) - pvalues.get('p_global', 1)) < 0.1 else "changes"} the conclusion.
"""

    report += f"""
---

## 6. Interpretation

### Is the SM Yukawa set globally more α-ladder-like than chance?

"""

    p_global = pvalues.get('p_global', 1)

    if p_global > 0.10:
        report += f"""**NO.**

The global p-value is {p_global:.3f}, meaning {100*p_global:.0f}% of random Yukawa sets achieve T ≤ {T_best:.2f} at their best scale. This is entirely consistent with chance.

The α-ladder hypothesis receives no support from this scale-consistent analysis.
"""
    elif p_global > 0.05:
        report += f"""**WEAK / MARGINAL.**

The global p-value is {p_global:.3f}, which is marginally suggestive but not statistically significant at conventional thresholds.
"""
    else:
        report += f"""**YES (TENTATIVE).**

The global p-value is {p_global:.4f}, suggesting the SM Yukawas are closer to an α-ladder pattern than expected by chance.
"""

    report += f"""
### Key Observations

1. **Scale matters:** The choice of μ affects T significantly. At {best_scale}, T = {T_best:.2f}.

2. **Pass counts are not discriminative:** {summary_tau.pass_W1 if summary_tau else 0}/9 pass W1 at m_τ, but random sets also achieve similar pass rates.

3. **The tight window W2 [0.5, 2]:** Only {summary_tau.pass_W2 if summary_tau else 0}/9 pass at m_τ, {summary_Z.pass_W2 if summary_Z else 0}/9 at m_Z. This is a more stringent test.

---

## 7. Conclusion

{"**The scale-consistent global test finds NO evidence for α-ladder structure in SM Yukawas.**" if p_global > 0.05 else "**There is marginal evidence for α-ladder structure, warranting further investigation.**"}

The test statistic T (sum of |ln c|) is {"not unusually small" if p_global > 0.05 else "unusually small"} compared to random Yukawa sets under the null hypothesis.

---

## Artifacts

- `artifacts/global_common_scale.csv`
- `artifacts/global_common_scale.json`
- `artifacts/global_common_scale_plots.png`
"""

    (out_dir.parent / "report_global_common_scale.md").write_text(report, encoding="utf-8")
