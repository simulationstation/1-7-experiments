"""
QCD/QED running and charm robustness analysis.

This module implements:
1. Loading of PDG running coupling constants
2. QCD running of alpha_s and quark masses
3. Charm robustness tests (Test A and Test B as defined in the study)
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np


@dataclass(frozen=True)
class RunningConstants:
    """Container for PDG running coupling constants."""
    alpha_0: float           # alpha(q^2=0)
    alpha_m_tau: float       # alpha_hat(m_tau)
    alpha_m_Z: float         # alpha_hat(m_Z)
    alpha_s_m_Z: float       # alpha_s(m_Z)
    alpha_s_m_tau: float     # alpha_s(m_tau)
    m_tau_GeV: float
    m_Z_GeV: float
    m_c_GeV: float           # m_c(m_c) MS-bar


def load_running_constants(path: Path) -> RunningConstants:
    """Load PDG running constants from JSON file."""
    data = json.loads(path.read_text(encoding="utf-8"))
    c = data["constants"]
    s = data["reference_scales_GeV"]
    return RunningConstants(
        alpha_0=c["alpha_0"]["value"],
        alpha_m_tau=c["alpha_hat_m_tau"]["value"],
        alpha_m_Z=c["alpha_hat_m_Z"]["value"],
        alpha_s_m_Z=c["alpha_s_m_Z"]["value"],
        alpha_s_m_tau=c["alpha_s_m_tau"]["value"],
        m_tau_GeV=s["m_tau"]["value"],
        m_Z_GeV=s["m_Z"]["value"],
        m_c_GeV=s["m_c_at_m_c"]["value"],
    )


def beta_0(n_f: int) -> float:
    """1-loop QCD beta function coefficient: beta_0 = 11 - 2*n_f/3."""
    return 11.0 - 2.0 * n_f / 3.0


def gamma_m_0() -> float:
    """1-loop quark mass anomalous dimension coefficient."""
    return 4.0


def run_alpha_s_1loop(alpha_s_ref: float, mu_ref: float, mu_target: float, n_f: int) -> float:
    """
    Run alpha_s from mu_ref to mu_target using 1-loop RGE.

    alpha_s(mu) = alpha_s(mu_ref) / (1 + alpha_s(mu_ref) * beta_0 / (2*pi) * ln(mu/mu_ref))
    """
    if mu_target <= 0 or mu_ref <= 0:
        raise ValueError("Scales must be positive")

    b0 = beta_0(n_f)
    L = math.log(mu_target / mu_ref)
    denom = 1.0 + alpha_s_ref * b0 / (2.0 * math.pi) * L

    if denom <= 0:
        raise ValueError(f"alpha_s running hit Landau pole: denom={denom}")

    return alpha_s_ref / denom


def run_quark_mass_1loop(
    m_ref: float,
    mu_ref: float,
    mu_target: float,
    alpha_s_ref: float,
    alpha_s_target: float,
    n_f: int,
) -> float:
    """
    Run MS-bar quark mass from mu_ref to mu_target using 1-loop RGE.

    m(mu2) / m(mu1) = [alpha_s(mu2) / alpha_s(mu1)]^(gamma_m_0 / beta_0)

    This is the leading-order formula for quark mass running.
    """
    b0 = beta_0(n_f)
    gm0 = gamma_m_0()

    if alpha_s_ref <= 0 or alpha_s_target <= 0:
        raise ValueError("alpha_s must be positive")

    ratio = alpha_s_target / alpha_s_ref
    exponent = gm0 / b0

    return m_ref * (ratio ** exponent)


def compute_alpha_s_at_scale(
    mu_target: float,
    consts: RunningConstants,
) -> Tuple[float, int]:
    """
    Compute alpha_s at mu_target using 1-loop running from appropriate reference.

    Returns (alpha_s, n_f) where n_f is the number of active flavors.

    Strategy:
    - If mu_target >= m_b (~4.18 GeV): run down from m_Z with n_f=5
    - If m_c <= mu_target < m_b: run from m_tau with n_f=4
    - If mu_target < m_c: run from m_tau with n_f=4 (simplified)

    Note: This is a simplified approach. Full treatment would match at thresholds.
    """
    m_b = 4.18  # approximate b threshold
    m_c = consts.m_c_GeV

    if mu_target >= m_b:
        # Run from m_Z with n_f=5
        alpha_s = run_alpha_s_1loop(
            consts.alpha_s_m_Z,
            consts.m_Z_GeV,
            mu_target,
            n_f=5
        )
        return alpha_s, 5
    else:
        # Run from m_tau with n_f=4
        alpha_s = run_alpha_s_1loop(
            consts.alpha_s_m_tau,
            consts.m_tau_GeV,
            mu_target,
            n_f=4
        )
        return alpha_s, 4


def run_charm_mass_to_scale(
    mu_target: float,
    consts: RunningConstants,
) -> float:
    """
    Run the charm quark MS-bar mass from m_c(m_c) to mu_target.

    Uses 1-loop QCD running. The number of active flavors depends on the scale.
    """
    m_c_ref = consts.m_c_GeV
    mu_ref = m_c_ref  # m_c is defined at mu = m_c

    # Get alpha_s at both scales
    alpha_s_ref, n_f_ref = compute_alpha_s_at_scale(mu_ref, consts)
    alpha_s_target, n_f_target = compute_alpha_s_at_scale(mu_target, consts)

    # For simplicity, use the target n_f for running
    # (A more sophisticated treatment would match at thresholds)
    m_target = run_quark_mass_1loop(
        m_ref=m_c_ref,
        mu_ref=mu_ref,
        mu_target=mu_target,
        alpha_s_ref=alpha_s_ref,
        alpha_s_target=alpha_s_target,
        n_f=n_f_target,
    )

    return m_target


def yukawa_from_mass(mass_GeV: float, vev_GeV: float) -> float:
    """Compute Yukawa coupling y = sqrt(2) * m / v."""
    return math.sqrt(2.0) * mass_GeV / vev_GeV


@dataclass(frozen=True)
class CharmScaleResult:
    """Result of charm analysis at a single scale."""
    scale_name: str
    mu_GeV: float
    alpha_value: float
    m_c_GeV: float           # charm mass at this scale (or reference)
    y_c: float               # Yukawa at this scale
    c_coeff: float           # c = y_c / alpha
    d_metric: float          # d = |ln(c)|
    test_type: str           # "A" (alpha-only) or "B" (full running)


def run_test_A(
    y_c_ref: float,
    consts: RunningConstants,
) -> List[CharmScaleResult]:
    """
    TEST A: Alpha-at-scale substitution.

    Keep y_c fixed (from baseline), compute c = y_c / alpha(mu) for each scale.
    This tests: "Is charm special because of the alpha value used?"
    """
    results = []

    # Baseline: alpha(0)
    c_0 = y_c_ref / consts.alpha_0
    results.append(CharmScaleResult(
        scale_name="alpha(0)",
        mu_GeV=0.0,  # q^2 = 0
        alpha_value=consts.alpha_0,
        m_c_GeV=consts.m_c_GeV,
        y_c=y_c_ref,
        c_coeff=c_0,
        d_metric=abs(math.log(c_0)),
        test_type="A",
    ))

    # alpha(m_tau)
    c_tau = y_c_ref / consts.alpha_m_tau
    results.append(CharmScaleResult(
        scale_name="alpha(m_tau)",
        mu_GeV=consts.m_tau_GeV,
        alpha_value=consts.alpha_m_tau,
        m_c_GeV=consts.m_c_GeV,
        y_c=y_c_ref,
        c_coeff=c_tau,
        d_metric=abs(math.log(c_tau)),
        test_type="A",
    ))

    # alpha(m_Z)
    c_Z = y_c_ref / consts.alpha_m_Z
    results.append(CharmScaleResult(
        scale_name="alpha(m_Z)",
        mu_GeV=consts.m_Z_GeV,
        alpha_value=consts.alpha_m_Z,
        m_c_GeV=consts.m_c_GeV,
        y_c=y_c_ref,
        c_coeff=c_Z,
        d_metric=abs(math.log(c_Z)),
        test_type="A",
    ))

    return results


def run_test_B(
    vev_GeV: float,
    consts: RunningConstants,
) -> List[CharmScaleResult]:
    """
    TEST B: Common-scale running.

    Run both m_c and use alpha to common scales (m_tau, m_Z).
    This tests: "Is charm special when y_c and alpha are at the same scale?"
    """
    results = []

    # Baseline: m_c(m_c), alpha(0) - for comparison
    y_c_ref = yukawa_from_mass(consts.m_c_GeV, vev_GeV)
    c_0 = y_c_ref / consts.alpha_0
    results.append(CharmScaleResult(
        scale_name="baseline: m_c(m_c), alpha(0)",
        mu_GeV=consts.m_c_GeV,
        alpha_value=consts.alpha_0,
        m_c_GeV=consts.m_c_GeV,
        y_c=y_c_ref,
        c_coeff=c_0,
        d_metric=abs(math.log(c_0)),
        test_type="B",
    ))

    # At m_tau: run m_c to m_tau, use alpha(m_tau)
    m_c_at_tau = run_charm_mass_to_scale(consts.m_tau_GeV, consts)
    y_c_at_tau = yukawa_from_mass(m_c_at_tau, vev_GeV)
    c_tau = y_c_at_tau / consts.alpha_m_tau
    results.append(CharmScaleResult(
        scale_name="m_c(m_tau), alpha(m_tau)",
        mu_GeV=consts.m_tau_GeV,
        alpha_value=consts.alpha_m_tau,
        m_c_GeV=m_c_at_tau,
        y_c=y_c_at_tau,
        c_coeff=c_tau,
        d_metric=abs(math.log(c_tau)),
        test_type="B",
    ))

    # At m_Z: run m_c to m_Z, use alpha(m_Z)
    m_c_at_Z = run_charm_mass_to_scale(consts.m_Z_GeV, consts)
    y_c_at_Z = yukawa_from_mass(m_c_at_Z, vev_GeV)
    c_Z = y_c_at_Z / consts.alpha_m_Z
    results.append(CharmScaleResult(
        scale_name="m_c(m_Z), alpha(m_Z)",
        mu_GeV=consts.m_Z_GeV,
        alpha_value=consts.alpha_m_Z,
        m_c_GeV=m_c_at_Z,
        y_c=y_c_at_Z,
        c_coeff=c_Z,
        d_metric=abs(math.log(c_Z)),
        test_type="B",
    ))

    return results


def compute_null_pvalue_test_A(
    d_observed: float,
    consts: RunningConstants,
    n_trials: int = 100000,
    seed: int = 42,
) -> Tuple[float, float]:
    """
    Compute p-value for TEST A under null hypothesis.

    Null: y is log-uniform in [1e-6, 1].
    Statistic: D_A = min over {alpha(0), alpha(m_tau), alpha(m_Z)} of |ln(y/alpha)|

    Returns (p_value, mean_D_under_null).
    """
    rng = np.random.default_rng(seed)

    alphas = [consts.alpha_0, consts.alpha_m_tau, consts.alpha_m_Z]

    null_D = []
    for _ in range(n_trials):
        y = 10 ** rng.uniform(-6, 0)
        # Compute d = |ln(c)| = |ln(y/alpha)| for each alpha
        ds = [abs(math.log(y / a)) for a in alphas]
        D = min(ds)
        null_D.append(D)

    null_D = np.array(null_D)
    p_value = float(np.mean(null_D <= d_observed))
    mean_D = float(np.mean(null_D))

    return p_value, mean_D


def compute_null_pvalue_test_B(
    d_observed: float,
    consts: RunningConstants,
    vev_GeV: float,
    n_trials: int = 100000,
    seed: int = 42,
) -> Tuple[float, float]:
    """
    Compute p-value for TEST B under null hypothesis.

    This is more complex because we need to simulate running a random "charm mass"
    to different scales. We approximate by:
    1. Generate a random mass m in [0.1, 10] GeV (charm-like range)
    2. Run it to m_tau and m_Z
    3. Compute D_B = min over {m_tau, m_Z} of |ln(y(mu)/alpha(mu))|

    Returns (p_value, mean_D_under_null).
    """
    rng = np.random.default_rng(seed)

    null_D = []
    for _ in range(n_trials):
        # Generate a random "charm-like" mass
        # Use log-uniform in a range that produces Yukawas in [1e-6, 1]
        # y = sqrt(2) * m / v, so m = y * v / sqrt(2)
        # For y in [1e-6, 1], m in [~1.7e-4, ~174] GeV
        # But for quark running to make sense, use [0.5, 5] GeV
        m_ref = 10 ** rng.uniform(-0.3, 0.7)  # ~0.5 to 5 GeV

        # Compute Yukawa at reference scale
        y_ref = yukawa_from_mass(m_ref, vev_GeV)

        # For running, we need alpha_s at the reference scale
        # Use a simplified model: assume the ratio m(mu2)/m(mu1) ~ (mu1/mu2)^0.2
        # This is roughly what QCD running gives for light quarks

        # At m_tau
        m_at_tau = m_ref * (m_ref / consts.m_tau_GeV) ** 0.2 if m_ref < consts.m_tau_GeV else m_ref * (consts.m_tau_GeV / m_ref) ** 0.2
        y_at_tau = yukawa_from_mass(m_at_tau, vev_GeV)
        c_tau = y_at_tau / consts.alpha_m_tau
        d_tau = abs(math.log(c_tau)) if c_tau > 0 else float('inf')

        # At m_Z
        m_at_Z = m_ref * (m_ref / consts.m_Z_GeV) ** 0.2 if m_ref < consts.m_Z_GeV else m_ref * (consts.m_Z_GeV / m_ref) ** 0.2
        y_at_Z = yukawa_from_mass(m_at_Z, vev_GeV)
        c_Z = y_at_Z / consts.alpha_m_Z
        d_Z = abs(math.log(c_Z)) if c_Z > 0 else float('inf')

        D = min(d_tau, d_Z)
        null_D.append(D)

    null_D = np.array(null_D)
    p_value = float(np.mean(null_D <= d_observed))
    mean_D = float(np.mean(null_D))

    return p_value, mean_D


def compute_null_pvalue_simple(
    d_observed: float,
    alpha: float,
    n_trials: int = 100000,
    seed: int = 42,
) -> Tuple[float, float]:
    """
    Compute p-value for a single alpha value (no look-elsewhere).

    Null: y is log-uniform in [1e-6, 1].
    Statistic: d = |ln(y/alpha)|

    Returns (p_value, mean_d_under_null).
    """
    rng = np.random.default_rng(seed)

    null_d = []
    for _ in range(n_trials):
        y = 10 ** rng.uniform(-6, 0)
        c = y / alpha
        d = abs(math.log(c))
        null_d.append(d)

    null_d = np.array(null_d)
    p_value = float(np.mean(null_d <= d_observed))
    mean_d = float(np.mean(null_d))

    return p_value, mean_d
