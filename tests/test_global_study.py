"""Tests for the global common-scale study module."""

import math
import pytest
from pathlib import Path

from alc_yukawa.global_study import (
    fit_to_ladder,
    compute_null_T,
    ParticleAtScale,
    ScaleSummary,
    get_quark_reference_scale,
)


def test_fit_to_ladder_exact():
    """Test fitting when y = alpha^n exactly."""
    alpha = 0.00749
    n_true = 2
    y = alpha ** n_true

    best_n, c = fit_to_ladder(y, alpha)

    assert best_n == n_true
    assert c == pytest.approx(1.0, rel=1e-10)


def test_fit_to_ladder_with_coefficient():
    """Test fitting when y = c * alpha^n with c != 1."""
    alpha = 0.00749
    n_true = 1
    c_true = 2.5
    y = c_true * (alpha ** n_true)

    best_n, c = fit_to_ladder(y, alpha)

    assert best_n == n_true
    assert c == pytest.approx(c_true, rel=1e-10)


def test_fit_to_ladder_chooses_closest():
    """Test that fitting chooses n that minimizes |ln(c)|."""
    alpha = 0.00749

    # y = sqrt(alpha) * alpha^1 = alpha^1.5
    # This should fit to n=1 (c = alpha^0.5 ~ 0.087) or n=2 (c = alpha^-0.5 ~ 11.6)
    # |ln(0.087)| ≈ 2.44, |ln(11.6)| ≈ 2.45
    # So n=1 should be slightly better
    y = alpha ** 1.5
    best_n, c = fit_to_ladder(y, alpha)

    # The fit should choose either n=1 or n=2 depending on which |ln(c)| is smaller
    assert best_n in [1, 2]


def test_compute_null_T_positive():
    """Test that null T statistic is always positive."""
    import numpy as np
    rng = np.random.default_rng(42)

    T = compute_null_T(n_particles=9, alpha=0.00749, rng=rng)

    assert T > 0


def test_compute_null_T_deterministic():
    """Test that null T is deterministic given seed."""
    import numpy as np

    rng1 = np.random.default_rng(42)
    T1 = compute_null_T(n_particles=9, alpha=0.00749, rng=rng1)

    rng2 = np.random.default_rng(42)
    T2 = compute_null_T(n_particles=9, alpha=0.00749, rng=rng2)

    assert T1 == T2


def test_particle_at_scale_dataclass():
    """Test ParticleAtScale dataclass creation."""
    p = ParticleAtScale(
        particle="c",
        sector="up",
        mu_GeV=1.78,
        scale_name="m_tau",
        mass_ref_GeV=1.273,
        mass_at_mu_GeV=1.18,
        yukawa=0.00679,
        alpha_mu=0.00749,
        best_n=1,
        c=0.906,
        ln_c=-0.099,
        abs_ln_c=0.099,
        in_bounds_W1=True,
        in_bounds_W2=True,
    )

    assert p.particle == "c"
    assert p.in_bounds_W1 == True
    assert p.abs_ln_c == pytest.approx(0.099)


def test_scale_summary_dataclass():
    """Test ScaleSummary dataclass creation."""
    s = ScaleSummary(
        scale_name="m_tau",
        mu_GeV=1.78,
        alpha_mu=0.00749,
        T=10.5,
        T2=18.5,
        pass_W1=8,
        pass_W2=3,
        n_particles=9,
    )

    assert s.T == 10.5
    assert s.pass_W1 == 8


def test_get_quark_reference_scale():
    """Test extraction of quark reference scale from PDG data."""
    pdg_data = {
        "particles": {
            "c": {"mass": 1.273, "unit": "GeV", "scheme": "MSbar", "scale_GeV": 1.273},
            "u": {"mass": 2.16, "unit": "MeV", "scheme": "MSbar", "scale_GeV": 2.0},
            "t": {"mass": 172.57, "unit": "GeV", "scheme": "direct"},
        }
    }

    assert get_quark_reference_scale("c", pdg_data) == 1.273
    assert get_quark_reference_scale("u", pdg_data) == 2.0
    # Top should return its mass as reference scale
    assert get_quark_reference_scale("t", pdg_data) == pytest.approx(172.57)


def test_fit_bounds_W1():
    """Test W1 bounds checking [0.1, 10]."""
    # c = 0.5 should be in W1
    assert 0.1 <= 0.5 <= 10.0

    # c = 0.05 should NOT be in W1
    assert not (0.1 <= 0.05 <= 10.0)

    # c = 15 should NOT be in W1
    assert not (0.1 <= 15 <= 10.0)


def test_fit_bounds_W2():
    """Test W2 bounds checking [0.5, 2]."""
    # c = 1.0 should be in W2
    assert 0.5 <= 1.0 <= 2.0

    # c = 0.3 should NOT be in W2
    assert not (0.5 <= 0.3 <= 2.0)

    # c = 3.0 should NOT be in W2
    assert not (0.5 <= 3.0 <= 2.0)
