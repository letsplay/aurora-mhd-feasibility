# tests/test_v11.py
"""
AURORA-MHD v11 Physics Alignment Tests
========================================
These tests verify that the v11 corrections bring the repo
into alignment with the v11 feasibility paper. They run
alongside the 58 v1 regression tests — both must pass.

Test naming convention:
  test_d5X_description → maps to v3 task D5.X
"""
import numpy as np
import pytest


# ═══════════════════════════════════════════════════════════
# D5.1 — Finite Solenoid (Fabry Factor)
# ═══════════════════════════════════════════════════════════

class TestD51_FiniteSolenoid:

    def test_fabry_factor_paper_geometry(self):
        """F(R=0.5, L=0.6) = 0.514 ± 0.01."""
        from src.physics_v11 import fabry_factor
        F = fabry_factor(0.5, 0.6)
        assert abs(F - 0.514) < 0.01, f"F = {F:.4f}"

    def test_fabry_factor_infinite_limit(self):
        """F → 1 as L/R → ∞."""
        from src.physics_v11 import fabry_factor
        F = fabry_factor(0.5, 100.0)
        assert F > 0.999

    def test_fabry_factor_pancake_limit(self):
        """F → 0 as L/R → 0."""
        from src.physics_v11 import fabry_factor
        F = fabry_factor(10.0, 0.01)
        assert F < 0.001

    def test_v11_more_turns_than_v1(self):
        """v11 finite solenoid needs more turns than v1
        infinite solenoid at same B."""
        from src.magnets import CoilDesigner
        cd = CoilDesigner()
        d1 = cd.design_solenoid(2.0, 0.5, T_op=20, length=0.6)
        d11 = cd.design_solenoid_v11(2.0, 0.5, T_op=20, length=0.6)
        assert d1 is not None and d11 is not None
        assert d11['N_turns'] > d1['N_turns'], \
            f"v11={d11['N_turns']}, v1={d1['N_turns']}"

    def test_v11_heavier_than_v1(self):
        """v11 coil mass > v1 coil mass at same B."""
        from src.magnets import CoilDesigner
        cd = CoilDesigner()
        d1 = cd.design_solenoid(2.0, 0.5, T_op=20, length=0.6)
        d11 = cd.design_solenoid_v11(2.0, 0.5, T_op=20, length=0.6)
        assert d11['m_coil'] > d1['m_total'], \
            f"v11={d11['m_coil']:.0f}, v1={d1['m_total']:.0f}"

    def test_v11_2T_coil_mass_range(self):
        """2T coil assembly: 250–700 kg.
        Paper: 498 kg. Tolerance wide due to Ic table differences.
        """
        from src.magnets import CoilDesigner
        cd = CoilDesigner()
        d = cd.design_solenoid_v11(2.0, 0.5, T_op=20, length=0.6)
        assert d is not None
        assert 250 < d['m_coil'] < 700, \
            f"m_coil = {d['m_coil']:.0f} kg"

    def test_v11_2T_flight_system_range(self):
        """2T flight system (coil + auxiliaries): 400–950 kg.
        Paper: 782 kg.
        """
        from src.magnets import CoilDesigner
        cd = CoilDesigner()
        d = cd.design_solenoid_v11(2.0, 0.5, T_op=20, length=0.6)
        assert d is not None
        assert 400 < d['m_flight'] < 950, \
            f"m_flight = {d['m_flight']:.0f} kg"

    def test_v11_returns_model_tag(self):
        """v11 design dict contains 'model': 'v11'."""
        from src.magnets import CoilDesigner
        cd = CoilDesigner()
        d = cd.design_solenoid_v11(2.0, 0.5)
        assert d['model'] == 'v11'


# ═══════════════════════════════════════════════════════════
# D5.2 — Post-Shock Velocity
# ═══════════════════════════════════════════════════════════

class TestD52_PostShockVelocity:

    def test_leo_v_ps(self):
        """v_ps(7800) = 1560 ± 10 m/s."""
        from src.physics_v11 import post_shock_velocity
        v_ps = post_shock_velocity(7800)
        assert abs(v_ps - 1560) < 10

    def test_linear_scaling(self):
        """v_ps ∝ v_free."""
        from src.physics_v11 import post_shock_velocity
        v1 = post_shock_velocity(5000)
        v2 = post_shock_velocity(10000)
        assert abs(v2 / v1 - 2.0) < 0.001

    def test_custom_ratio(self):
        """Adjustable compression ratio."""
        from src.physics_v11 import post_shock_velocity
        v = post_shock_velocity(7800, ratio=3.0)
        assert abs(v - 2600) < 1


# ═══════════════════════════════════════════════════════════
# D5.3 — Hall Parameter & Effective Area
# ═══════════════════════════════════════════════════════════

class TestD53_HallAndArea:

    def test_hall_correction_baseline(self):
        """Hall correction at β=0.82: 0.598 ± 0.01."""
        from src.physics_v11 import hall_correction
        c = hall_correction(0.82)
        assert abs(c - 0.598) < 0.01, f"c = {c:.4f}"

    def test_hall_no_effect(self):
        """β=0 → correction = 1.0."""
        from src.physics_v11 import hall_correction
        assert hall_correction(0.0) == 1.0

    def test_hall_strong_effect(self):
        """β >> 1 → correction → 0."""
        from src.physics_v11 import hall_correction
        assert hall_correction(10.0) < 0.01

    def test_B2dA_scales_B_squared(self):
        """∫B²dA ∝ B²."""
        from src.physics_v11 import compute_B2dA
        b1 = compute_B2dA(1.0, 0.5)
        b2 = compute_B2dA(2.0, 0.5)
        assert abs(b2 / b1 - 4.0) < 0.01

    def test_effective_area_2T(self):
        """A_eff(2T, R=0.5m) in 0.2–1.5 m²."""
        from src.physics_v11 import compute_B2dA, effective_area
        B2dA = compute_B2dA(2.0, 0.5)
        A = effective_area(B2dA, 2.0)
        assert 0.2 < A < 1.5, f"A_eff = {A:.2f} m²"

    def test_effective_area_much_less_than_frontal(self):
        """A_eff << frontal area of 4.5m body (~63 m²)."""
        from src.physics_v11 import compute_B2dA, effective_area
        B2dA = compute_B2dA(2.0, 0.5)
        A = effective_area(B2dA, 2.0)
        assert A < 5.0  # well below 63 m²


# ═══════════════════════════════════════════════════════════
# D5.4 — Kinetic Ceiling & Stuart Number
# ═══════════════════════════════════════════════════════════

class TestD54_KineticCeiling:

    def test_stuart_number_leo(self):
        """S(LEO, 2T) >> 100."""
        from src.physics_v11 import stuart_number
        # σ=50 S/m, B=2T, v_ps=1560, L=0.6, ρ_ps=5×3e-5
        S = stuart_number(50, 2.0, 1560, 0.6, 5 * 3e-5)
        assert S > 100, f"S = {S:.0f}"

    def test_stuart_scales_B_squared(self):
        """S ∝ B²."""
        from src.physics_v11 import stuart_number
        S1 = stuart_number(50, 1.0, 1560, 0.6, 1.5e-4)
        S2 = stuart_number(50, 2.0, 1560, 0.6, 1.5e-4)
        assert abs(S2 / S1 - 4.0) < 0.01

    def test_kinetic_ceiling_leo(self):
        """P_max(LEO) in 200–800 kW range.
        Paper: 456 kW at K=0.5, ṁ=0.75 kg/s.
        """
        from src.physics_v11 import kinetic_ceiling
        P = kinetic_ceiling(0.75, 1560, K=0.5)
        assert 200e3 < P < 800e3, \
            f"P_max = {P/1e3:.0f} kW"

    def test_ceiling_scales_v_squared(self):
        """P_max ∝ v²_ps."""
        from src.physics_v11 import kinetic_ceiling
        P1 = kinetic_ceiling(1.0, 1000, K=0.5)
        P2 = kinetic_ceiling(1.0, 2000, K=0.5)
        assert abs(P2 / P1 - 4.0) < 0.01

    def test_ceiling_scales_linearly_with_K(self):
        """P_max ∝ K."""
        from src.physics_v11 import kinetic_ceiling
        P1 = kinetic_ceiling(1.0, 1000, K=0.3)
        P2 = kinetic_ceiling(1.0, 1000, K=0.9)
        assert abs(P2 / P1 - 3.0) < 0.01


# ═══════════════════════════════════════════════════════════
# D5.4 — Full v11 Pipeline
# ═══════════════════════════════════════════════════════════

class TestD54_V11Pipeline:

    @pytest.fixture
    def atmo(self):
        from src.atmosphere import EarthAtmosphere
        return EarthAtmosphere()

    def test_v11_extraction_positive(self, atmo):
        """v11 extraction > 0 at LEO conditions."""
        from src.physics_v11 import faraday_power_v11
        from src.plasma import sigma_at
        sig = sigma_at(60e3, 7800, atmo)
        P = faraday_power_v11(sig, 7800, 2.0, 60e3, atmo)
        assert P > 0, f"P = {P:.0f} W"

    def test_v11_much_less_than_v1(self, atmo):
        """v11 extraction << v1 extraction (>50× reduction).
        v1 uses v_free and 20 m² area → MW range.
        v11 uses v_ps, A_eff, Hall, ceiling → kW range.
        """
        from src.physics_v11 import faraday_power_v11
        from src.mhd import faraday_power
        from src.plasma import sigma_at

        sig = sigma_at(60e3, 7800, atmo)
        P_v1 = faraday_power(sig, 7800, 2.0)
        P_v11 = faraday_power_v11(sig, 7800, 2.0, 60e3, atmo)
        ratio = P_v1 / P_v11 if P_v11 > 0 else float('inf')
        assert ratio > 50, \
            f"v1={P_v1/1e3:.0f} kW, v11={P_v11/1e3:.0f} kW, ratio={ratio:.0f}×"

    def test_v11_zero_at_zero(self, atmo):
        """P = 0 when any input is 0."""
        from src.physics_v11 import faraday_power_v11
        assert faraday_power_v11(0, 7800, 2.0, 60e3, atmo) == 0
        assert faraday_power_v11(50, 0, 2.0, 60e3, atmo) == 0
        assert faraday_power_v11(50, 7800, 0, 60e3, atmo) == 0


# ═══════════════════════════════════════════════════════════
# D5.5 — Stuart Number Grid
# ═══════════════════════════════════════════════════════════

class TestD55_StuartGrid:

    def test_grid_shape(self):
        from src.atmosphere import EarthAtmosphere
        from src.envelope import compute_stuart_grid
        atmo = EarthAtmosphere()
        h, v, S = compute_stuart_grid(atmo, B=2.0, n_h=10, n_v=10)
        assert S.shape == (10, 10)

    def test_effective_zone_exists(self):
        """At B=2T, S > 10 exists somewhere."""
        from src.atmosphere import EarthAtmosphere
        from src.envelope import compute_stuart_grid, classify_stuart_zones
        atmo = EarthAtmosphere()
        _, _, S = compute_stuart_grid(atmo, B=2.0, n_h=20, n_v=20)
        zones = classify_stuart_zones(S)
        assert np.any(zones == 2), "No S > 10 zone found"

    def test_S_scales_with_B_squared(self):
        """S(2T) ≈ 4 × S(1T) at same conditions."""
        from src.atmosphere import EarthAtmosphere
        from src.envelope import compute_stuart_grid
        atmo = EarthAtmosphere()
        _, _, S1 = compute_stuart_grid(atmo, B=1.0, n_h=5, n_v=5)
        _, _, S2 = compute_stuart_grid(atmo, B=2.0, n_h=5, n_v=5)
        # Check ratio at a point where both are nonzero
        mask = (S1 > 0.01) & (S2 > 0.01)
        if mask.any():
            ratios = S2[mask] / S1[mask]
            assert np.mean(ratios) > 3.0  # should be ~4


# ═══════════════════════════════════════════════════════════
# D5.7 — COP & Ride-Through Demand
# ═══════════════════════════════════════════════════════════

class TestD57_Demand:

    def test_ride_through_under_15kW(self):
        """Ride-through demand < 15 kW."""
        from src.mhd import PowerDemand
        coil = {'I_op': 336, 'radius': 0.5}
        pd = PowerDemand(coil)
        P = pd.total_v11(mode='ride_through')
        assert P < 15e3, f"P = {P/1e3:.1f} kW"

    def test_conservative_higher_than_ride_through(self):
        """Conservative demand > ride-through during reentry
        (when aero heating drives cryocooler load)."""
        from src.mhd import PowerDemand
        coil = {'I_op': 336, 'radius': 0.5}
        pd = PowerDemand(coil)
        P_rt = pd.total_v11(q_aero=200e3, mode='ride_through')
        P_con = pd.total_v11(q_aero=200e3, mode='conservative')
        assert P_con > P_rt

    def test_v11_demand_higher_than_v1(self):
        """v11 conservative demand > v1 demand
        (because v1 used COP=0.015 vs v11 COP=0.002).
        """
        from src.mhd import PowerDemand
        coil = {'I_op': 336, 'radius': 0.5}
        pd = PowerDemand(coil)
        P_v1 = pd.total(q_aero=0)
        P_v11 = pd.total_v11(q_aero=0, mode='conservative')
        assert P_v11 > P_v1

    def test_invalid_mode_raises(self):
        from src.mhd import PowerDemand
        coil = {'I_op': 336, 'radius': 0.5}
        pd = PowerDemand(coil)
        with pytest.raises(ValueError):
            pd.total_v11(mode='nonsense')


# ═══════════════════════════════════════════════════════════
# D5.8 — v11 Energy Balance
# ═══════════════════════════════════════════════════════════

class TestD58_EnergyBalance:

    @pytest.fixture
    def leo_setup(self):
        from src.atmosphere import EarthAtmosphere
        from src.trajectory import ReentryTrajectory
        from src.magnets import CoilDesigner
        atmo = EarthAtmosphere()
        traj = ReentryTrajectory(atmo)
        r = traj.run(v0=7800, gamma0_deg=-2)
        cd = CoilDesigner()
        coil = cd.design_solenoid_v11(2.0, 0.5, T_op=20, length=0.6)
        return r, atmo, coil

    def test_v11_balance_runs(self, leo_setup):
        """energy_balance_v11 completes without error."""
        from src.mhd import energy_balance_v11
        r, atmo, coil = leo_setup
        eb = energy_balance_v11(r, 2.0, atmo, coil)
        assert 'P_extract' in eb
        assert 'peak_margin' in eb
        assert eb['model'] == 'v11'

    def test_v11_margin_positive(self, leo_setup):
        """2T ride-through margin > 1.0 (energy-positive)."""
        from src.mhd import energy_balance_v11
        r, atmo, coil = leo_setup
        eb = energy_balance_v11(r, 2.0, atmo, coil)
        assert eb['peak_margin'] > 1.0, \
            f"margin = {eb['peak_margin']:.1f}×"

    def test_v11_margin_less_than_v1(self, leo_setup):
        """v11 margin << v1 margin (v1 was 24×+)."""
        from src.mhd import energy_balance, energy_balance_v11
        r, atmo, coil = leo_setup
        eb_v1 = energy_balance(r, 2.0, atmo, coil)
        eb_v11 = energy_balance_v11(r, 2.0, atmo, coil)

        peak_v1 = eb_v1['P_extract'].max()
        peak_v11 = eb_v11['peak_extraction_W']
        assert peak_v1 > peak_v11 * 10, \
            f"v1 peak={peak_v1/1e3:.0f} kW, v11 peak={peak_v11/1e3:.0f} kW"

    def test_v11_battery_tracked(self, leo_setup):
        """Battery state is tracked and has correct shape."""
        from src.mhd import energy_balance_v11
        r, atmo, coil = leo_setup
        eb = energy_balance_v11(r, 2.0, atmo, coil)
        assert 'E_battery' in eb
        assert len(eb['E_battery']) == len(eb['t'])
        # Starts at full charge
        assert eb['E_battery'][0] == 50e3 * 3600

    def test_higher_B_better_margin(self, leo_setup):
        """Higher B → better margin in v11 too."""
        from src.mhd import energy_balance_v11
        from src.magnets import CoilDesigner
        r, atmo, _ = leo_setup
        cd = CoilDesigner()

        margins = []
        for B in [1.0, 2.0, 3.0]:
            coil = cd.design_solenoid_v11(B, 0.5, T_op=20, length=0.6)
            if coil is None:
                continue
            eb = energy_balance_v11(r, B, atmo, coil)
            margins.append(eb['peak_margin'])

        # Margins should generally increase with B
        assert margins[-1] >= margins[0], \
            f"margins = {margins}"


# ═══════════════════════════════════════════════════════════
# v1 Regression Sanity
# ═══════════════════════════════════════════════════════════

class TestV1Regression:
    """Verify v1 functions still work unchanged."""

    def test_v1_faraday_unchanged(self):
        """v1 faraday_power still returns MW-range at LEO."""
        from src.mhd import faraday_power
        P = faraday_power(50, 7800, 2.0)
        assert P > 1e6, f"P = {P/1e3:.0f} kW"

    def test_v1_design_unchanged(self):
        """v1 design_solenoid still works and has model tag."""
        from src.magnets import CoilDesigner
        cd = CoilDesigner()
        d = cd.design_solenoid(2.0, 0.5, T_op=20)
        assert d is not None
        assert d['model'] == 'v1'

    def test_v1_energy_balance_unchanged(self):
        """v1 energy_balance returns model tag."""
        from src.atmosphere import EarthAtmosphere
        from src.trajectory import ReentryTrajectory
        from src.magnets import CoilDesigner
        from src.mhd import energy_balance
        atmo = EarthAtmosphere()
        traj = ReentryTrajectory(atmo)
        r = traj.run(v0=7800)
        cd = CoilDesigner()
        coil = cd.design_solenoid(2.0, 2.0, T_op=20)
        eb = energy_balance(r, 2.0, atmo, coil)
        assert eb['model'] == 'v1'
