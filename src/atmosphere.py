# src/atmosphere.py
import numpy as np
from scipy.interpolate import interp1d
from datetime import datetime

class EarthAtmosphere:
    """NRLMSISE-00 atmospheric model with
    fast interpolation lookup."""
    
    def __init__(self, f107=150, ap=4):
        from nrlmsise00 import msise_model
        self._alts = np.arange(0, 121e3, 500)
        rhos, temps = [], []
        dt = datetime(2025, 6, 21, 12)
        
        for h in self._alts:
            out = msise_model(
                dt, h/1e3, 0, 0,
                f107, f107, ap
            )
            rhos.append(out[0][5]*1000)   # g/cm3 -> kg/m3
            temps.append(out[1][1])  # K
        
        self._rho_interp = interp1d(
            self._alts, rhos,
            kind='linear',
            fill_value=(rhos[0], rhos[-1]),
            bounds_error=False
        )
        self._T_interp = interp1d(
            self._alts, temps,
            kind='linear',
            fill_value=(temps[0], temps[-1]),
            bounds_error=False
        )
    
    def density(self, h):
        return float(self._rho_interp(h))
    
    def temperature(self, h):
        return float(self._T_interp(h))

class MarsAtmosphere:
    """Simplified Mars atmospheric model.
    Exponential fit from Mars-GRAM."""
    
    RHO_0 = 0.020     # kg/m3 at surface
    SCALE_H = 11.1e3   # m
    T_SURFACE = 210     # K
    
    def density(self, h):
        return self.RHO_0 * np.exp(
            -h / self.SCALE_H
        )
    
    def temperature(self, h):
        # Simple lapse rate
        return max(self.T_SURFACE - 2.5e-3*h,
                   130)
    