# src/envelope.py
import numpy as np
from src.plasma import sigma_at
from src.atmosphere import (
    EarthAtmosphere, MarsAtmosphere
)

def compute_sigma_grid(atmo, species='air',
                       seed_frac=0.0,
                       h_range=(30e3, 120e3),
                       v_range=(2e3, 15e3),
                       n_h=200, n_v=200):
    """Compute σ over (altitude, velocity) 
    grid.
    
    Returns: h_arr, v_arr, sigma_2d
    """
    h_arr = np.linspace(*h_range, n_h)
    v_arr = np.linspace(*v_range, n_v)
    sigma_2d = np.zeros((n_h, n_v))
    
    for i, h in enumerate(h_arr):
        for j, v in enumerate(v_arr):
            sigma_2d[i, j] = sigma_at(
                h, v, atmo,
                species=species,
                seed_frac=seed_frac
            )
    
    return h_arr, v_arr, sigma_2d

def classify_zones(sigma_2d):
    """Classify into MHD effectiveness zones.
    Returns: zone_2d (0=ineffective, 
    1=marginal, 2=effective)
    """
    zones = np.zeros_like(sigma_2d, dtype=int)
    zones[sigma_2d >= 0.1] = 1   # marginal
    zones[sigma_2d >= 10.0] = 2  # effective
    return zones

# src/envelope.py (continued)

def time_in_mhd_zone(traj_result, atmo,
                     species='air',
                     seed_frac=0.0,
                     sigma_threshold=10.0):
    """Compute total seconds where σ 
    exceeds threshold along trajectory."""
    t = traj_result['t']
    h = traj_result['h']
    v = traj_result['v']
    
    in_zone = np.zeros(len(t), dtype=bool)
    for i in range(len(t)):
        sig = sigma_at(
            h[i], v[i], atmo,
            species=species,
            seed_frac=seed_frac
        )
        in_zone[i] = sig >= sigma_threshold
    
    # Integrate time where in_zone is True
    dt = np.diff(t)
    time_s = np.sum(dt[in_zone[:-1]])
    return float(time_s)


# ─────────────────────────────────────────────────────────
# v11: Stuart Number Grid (D5.5)
# ─────────────────────────────────────────────────────────

def compute_stuart_grid(atmo, B=2.0,
                        species='air',
                        seed_frac=0.0,
                        h_range=(30e3, 120e3),
                        v_range=(2e3, 15e3),
                        n_h=200, n_v=200,
                        L_char=0.6,
                        v_ps_ratio=5.0):
    """Compute Stuart number S(h, v) over grid.

    S = σ B² L / (ρ_ps × v_ps)

    Zone classification:
      S > 10:  MHD effective (strong deflection)
      S > 1:   MHD significant (partial interaction)
      S < 0.1: MHD negligible

    Args:
        atmo: atmosphere object
        B: magnetic field [T]
        species: 'air' or 'co2'
        seed_frac: Cs seeding fraction
        h_range, v_range: grid bounds
        n_h, n_v: grid resolution
        L_char: interaction length [m]
        v_ps_ratio: v_free/v_ps compression ratio

    Returns:
        h_arr, v_arr, S_2d
    """
    from src.physics_v11 import (
        stuart_number, post_shock_velocity,
        COMPRESSION_RATIO
    )

    h_arr = np.linspace(*h_range, n_h)
    v_arr = np.linspace(*v_range, n_v)
    S_2d = np.zeros((n_h, n_v))

    for i, h in enumerate(h_arr):
        rho = atmo.density(h)
        rho_ps = COMPRESSION_RATIO * rho
        for j, v in enumerate(v_arr):
            sig = sigma_at(
                h, v, atmo,
                species=species,
                seed_frac=seed_frac
            )
            v_ps = post_shock_velocity(
                v, ratio=v_ps_ratio)
            S_2d[i, j] = stuart_number(
                sig, B, v_ps, L_char, rho_ps)

    return h_arr, v_arr, S_2d


def classify_stuart_zones(S_2d):
    """Classify Stuart number grid into MHD zones.

    Returns:
        zones: 0=negligible, 1=significant, 2=effective
    """
    zones = np.zeros_like(S_2d, dtype=int)
    zones[S_2d >= 1.0] = 1    # significant
    zones[S_2d >= 10.0] = 2   # effective
    return zones


