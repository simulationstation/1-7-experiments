"""Tests for the running module (QCD/QED running and charm robustness)."""

import math
import pytest
from pathlib import Path

from alc_yukawa.running import (
    beta_0,
    gamma_m_0,
    run_alpha_s_1loop,
    run_quark_mass_1loop,
    yukawa_from_mass,
    RunningConstants,
    load_running_constants,
    run_test_A,
    run_test_B,
)


def test_beta_0():
    """Test QCD beta function coefficient."""
    # beta_0 = 11 - 2*n_f/3
    assert beta_0(3) == 11 - 2  # = 9
    assert beta_0(4) == pytest.approx(11 - 8/3)  # ≈ 8.33
    assert beta_0(5) == pytest.approx(11 - 10/3)  # ≈ 7.67
    assert beta_0(6) == 11 - 4  # = 7


def test_gamma_m_0():
    """Test mass anomalous dimension."""
    assert gamma_m_0() == 4.0


def test_run_alpha_s_1loop_identity():
    """Running to same scale should return same alpha_s."""
    alpha_s = 0.12
    result = run_alpha_s_1loop(alpha_s, mu_ref=91.0, mu_target=91.0, n_f=5)
    assert result == pytest.approx(alpha_s, rel=1e-10)


def test_run_alpha_s_1loop_direction():
    """alpha_s should decrease at higher scales (asymptotic freedom)."""
    alpha_s_low = 0.3
    mu_low = 2.0
    mu_high = 91.0
    n_f = 4

    alpha_s_high = run_alpha_s_1loop(alpha_s_low, mu_low, mu_high, n_f)
    assert alpha_s_high < alpha_s_low  # asymptotic freedom


def test_run_quark_mass_1loop_direction():
    """Quark mass should decrease at higher scales."""
    m_ref = 1.27  # GeV
    mu_ref = 1.27
    mu_high = 91.0
    alpha_s_ref = 0.35
    alpha_s_high = 0.12

    m_high = run_quark_mass_1loop(
        m_ref, mu_ref, mu_high,
        alpha_s_ref, alpha_s_high, n_f=5
    )
    assert m_high < m_ref  # mass runs down at higher scales


def test_yukawa_from_mass():
    """Test Yukawa calculation."""
    # y = sqrt(2) * m / v
    m = 1.0  # GeV
    v = 246.0  # GeV
    y = yukawa_from_mass(m, v)
    expected = math.sqrt(2) * 1.0 / 246.0
    assert y == pytest.approx(expected)


def test_load_running_constants():
    """Test loading PDG constants from JSON."""
    path = Path("data/running_constants_pdg.json")
    if not path.exists():
        pytest.skip("PDG constants file not found")

    consts = load_running_constants(path)

    # Check reasonable values
    assert 0.007 < consts.alpha_0 < 0.008
    assert 0.007 < consts.alpha_m_tau < 0.008
    assert 0.007 < consts.alpha_m_Z < 0.008
    assert consts.alpha_m_Z > consts.alpha_0  # alpha runs up
    assert 0.1 < consts.alpha_s_m_Z < 0.13
    assert 1.7 < consts.m_tau_GeV < 1.8
    assert 90 < consts.m_Z_GeV < 92
    assert 1.2 < consts.m_c_GeV < 1.4


def test_run_test_A_baseline():
    """Test that TEST A includes baseline alpha(0) result."""
    consts = RunningConstants(
        alpha_0=0.0073,
        alpha_m_tau=0.0075,
        alpha_m_Z=0.0078,
        alpha_s_m_Z=0.118,
        alpha_s_m_tau=0.32,
        m_tau_GeV=1.78,
        m_Z_GeV=91.2,
        m_c_GeV=1.27,
    )
    y_c = 0.00731

    results = run_test_A(y_c, consts)

    assert len(results) == 3
    assert results[0].scale_name == "alpha(0)"
    assert results[0].y_c == y_c
    assert results[0].c_coeff == pytest.approx(y_c / consts.alpha_0)


def test_run_test_B_has_running():
    """Test that TEST B actually runs the charm mass."""
    consts = RunningConstants(
        alpha_0=0.0073,
        alpha_m_tau=0.0075,
        alpha_m_Z=0.0078,
        alpha_s_m_Z=0.118,
        alpha_s_m_tau=0.32,
        m_tau_GeV=1.78,
        m_Z_GeV=91.2,
        m_c_GeV=1.27,
    )
    vev = 246.0

    results = run_test_B(vev, consts)

    assert len(results) == 3

    # Baseline should have original mass
    assert results[0].m_c_GeV == consts.m_c_GeV

    # At m_Z, mass should be lower than at m_c
    m_c_at_Z = results[2].m_c_GeV
    assert m_c_at_Z < consts.m_c_GeV


def test_d_metric_zero_for_exact_match():
    """If c = 1 exactly, d = |ln(c)| = 0."""
    # d = |ln(c)| = 0 when c = 1
    c = 1.0
    d = abs(math.log(c))
    assert d == 0.0
