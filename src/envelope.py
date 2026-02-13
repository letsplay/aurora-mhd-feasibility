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


