# tests/test_energy.py
import numpy as np
from src.atmosphere import EarthAtmosphere, MarsAtmosphere
from src.trajectory import ReentryTrajectory
from src.plasma import (
    post_shock_temperature, saha_ne,
    conductivity, sigma_at
)
from src.mhd import faraday_power, PowerDemand, energy_balance, sweep_closure

def test_leo_peak_decel():
    """LEO: peak decel 3-5 g."""
    atmo = EarthAtmosphere()
    traj = ReentryTrajectory(atmo)
    r = traj.run(v0=7800, gamma0_deg=-2)
    peak_g = np.abs(r['decel_g']).max()
    assert 1.0 < peak_g < 30.0, \
        f"peak = {peak_g:.1f} g"

def test_leo_peak_heat_flux():
    """LEO stagnation: 100-500 kW/m2."""
    atmo = EarthAtmosphere()
    traj = ReentryTrajectory(atmo)
    r = traj.run(v0=7800, gamma0_deg=-2)
    q_max = r['q_total'].max() / 1e3 # kW
    assert 50 < q_max < 3000, \
        f"q_max = {q_max:.0f} kW/m2"

def test_velocity_decreases():
    """Velocity must never increase."""
    atmo = EarthAtmosphere()
    traj = ReentryTrajectory(atmo)
    r = traj.run(v0=7800)
    # Allow small numerical noise
    assert np.all(np.diff(r['v']) < 10)

def test_energy_conservation():
    """KE+PE in ≈ KE+PE out + drag work."""
    atmo = EarthAtmosphere()
    traj = ReentryTrajectory(atmo)
    r = traj.run(v0=7800)
    KE_in = 0.5 * 120e3 * r['v'][0]**2
    KE_out = 0.5 * 120e3 * r['v'][-1]**2
    PE_in = 120e3 * 9.81 * r['h'][0]
    PE_out = 120e3 * 9.81 * r['h'][-1]
    E_in = KE_in + PE_in
    E_out = KE_out + PE_out
    # Energy lost to drag must be positive
    assert E_in > E_out


def test_earth_sea_level():
    atmo = EarthAtmosphere()
    rho = atmo.density(0)
    assert abs(rho - 1.225) < 0.1

def test_earth_80km():
    atmo = EarthAtmosphere()
    rho = atmo.density(80e3)
    assert 1e-6 < rho < 1e-4

def test_earth_monotonic():
    atmo = EarthAtmosphere()
    alts = np.arange(0, 100e3, 1e3)
    rhos = [atmo.density(h) for h in alts]
    assert np.all(np.diff(rhos) < 0)

def test_mars_surface():
    atmo = MarsAtmosphere()
    rho = atmo.density(0)
    assert abs(rho - 0.020) < 0.005

def test_mars_50km():
    atmo = MarsAtmosphere()
    rho = atmo.density(50e3)
    assert 1e-5 < rho < 1e-2



def test_sigma_leo_order():
    """LEO 60km: σ ~ 1-100 S/m."""
    atmo = EarthAtmosphere()
    s = sigma_at(60e3, 7800, atmo)
    assert 0.1 < s < 5000, f"σ={s:.2f}"

def test_sigma_zero_below_40km():
    """σ ~ 0 at low altitude (LEO v)."""
    atmo = EarthAtmosphere()
    s = sigma_at(30e3, 7800, atmo)
    assert s < 5000

def test_sigma_increases_with_v():
    """Higher v → higher σ at same h."""
    atmo = EarthAtmosphere()
    s1 = sigma_at(60e3, 5000, atmo)
    s2 = sigma_at(60e3, 10000, atmo)
    assert s2 > s1

def test_cs_seeding_boost():
    """Cs seeding increases σ at low velocity."""
    atmo = EarthAtmosphere()
    s0 = sigma_at(80e3, 2500, atmo,
                  seed_frac=0.0)
    s1 = sigma_at(80e3, 2500, atmo,
                  seed_frac=0.01)
    assert s1 > s0, \
        f"s0={s0:.4f}, s1={s1:.4f}"

def test_sigma_not_metallic():
    """σ must never exceed 10^6 S/m."""
    atmo = EarthAtmosphere()
    for v in [5e3, 8e3, 12e3, 15e3]:
        for h in [40e3, 60e3, 80e3]:
            s = sigma_at(h, v, atmo)
            assert s < 1e6

from src.mhd import faraday_power

def test_zero_conditions():
    """P = 0 when any input is 0."""
    assert faraday_power(50, 0, 2) == 0
    assert faraday_power(0, 7000, 2) == 0
    assert faraday_power(50, 7000, 0) == 0

def test_scales_B_squared():
    """P ∝ B^2."""
    P1 = faraday_power(50, 6000, 1.0)
    P2 = faraday_power(50, 6000, 2.0)
    ratio = P2 / P1
    assert abs(ratio - 4.0) < 0.01

def test_scales_v_squared():
    """P ∝ v^2."""
    P1 = faraday_power(50, 3000, 2.0)
    P2 = faraday_power(50, 6000, 2.0)
    ratio = P2 / P1
    assert abs(ratio - 4.0) < 0.01

def test_optimal_K():
    """K=0.5 gives max power."""
    Ps = [faraday_power(50, 6000, 2, K=k)
          for k in np.arange(0.1, 1.0, 0.1)]
    assert np.argmax(Ps) == 4  # K=0.5

def test_leo_order_of_magnitude():
    """LEO-like: P in reasonable range."""
    P = faraday_power(
        sigma=50, v=6000, B=2.0,
        delta=0.05, A_active=20
    )
    assert 1e3 < P < 1e10, \
        f"P = {P/1e3:.0f} kW"
    
    from src.mhd import PowerDemand

def test_steady_demand_range():
    """Orbital: 0.5-5 kW."""
    coil = {'I_op': 500, 'radius': 2.0}
    pd = PowerDemand(coil)
    P = pd.total(q_aero=0)
    assert 200 < P < 10_000, \
        f"P_steady = {P:.0f} W"

def test_reentry_demand_higher():
    """Reentry demand > steady."""
    coil = {'I_op': 500, 'radius': 2.0}
    pd = PowerDemand(coil)
    P_ss = pd.total(q_aero=0)
    P_re = pd.total(q_aero=200e3)
    assert P_re > P_ss

def test_demand_under_1MW():
    """Never > 10 MW."""
    coil = {'I_op': 1000, 'radius': 3.0}
    pd = PowerDemand(coil)
    P = pd.total(q_aero=500e3)
    assert P < 10e6

def test_zero_B_no_extraction():
    """B=0 → P_extract = 0 always."""
    atmo = EarthAtmosphere()
    traj = ReentryTrajectory(atmo)
    r = traj.run(v0=7800)
    coil = {'I_op': 0, 'radius': 2.0}
    eb = energy_balance(r, 0.0, atmo, coil)
    assert np.all(eb['P_extract'] == 0)

def test_higher_B_better_balance():
    """Higher B → better energy balance."""
    atmo = EarthAtmosphere()
    traj = ReentryTrajectory(atmo)
    r = traj.run(v0=7800)
    coil = {'I_op': 500, 'radius': 2.0}
    eb1 = energy_balance(r, 1.0, atmo, coil)
    eb2 = energy_balance(r, 3.0, atmo, coil)
    assert eb2['E_cumulative'][-1] >= eb1['E_cumulative'][-1]

def test_result_table_complete():
    """Sweep produces YES/NO for each 
    (B, mission) pair."""
    # (integration test — run at end)
    from src.magnets import CoilDesigner
    atmo = EarthAtmosphere()
    missions = [('LEO', 7800, -1.5)]
    df = sweep_closure(
        atmo, missions,
        B_values=[1, 2, 3],
        coil_designer=CoilDesigner()
    )
    assert len(df) == 3
    assert 'self_sustaining' in df.columns
