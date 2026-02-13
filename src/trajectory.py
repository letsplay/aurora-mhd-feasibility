# src/trajectory.py
import numpy as np
from scipy.integrate import solve_ivp

R_EARTH = 6371e3       # m
G_EARTH = 9.81         # m/s2

class ReentryTrajectory:
    """3-DOF reentry trajectory model."""
    
    def __init__(self, atmo, mass=120e3,
                 Cd=1.2, Cl=0.3, A_ref=70.0,
                 R_nose=4.5):
        """
        atmo: atmosphere object with rho(h)
        mass: vehicle mass [kg]
        Cd, Cl: drag, lift coefficients
        A_ref: reference area [m2]
        R_nose: nose radius [m] for heat flux
        """
        self.atmo = atmo
        self.mass = mass
        self.Cd = Cd
        self.Cl = Cl
        self.A_ref = A_ref
        self.R_nose = R_nose
    
    def _odes(self, t, y):
        v, gamma, h = y
        if h < 0 or v < 100:
            return [0, 0, 0]
        
        rho = self.atmo.density(h)
        q = 0.5 * rho * v**2
        D = q * self.Cd * self.A_ref
        L = q * self.Cl * self.A_ref
        m = self.mass
        g = G_EARTH * (R_EARTH/(R_EARTH+h))**2
        
        dv = -D/m - g * np.sin(gamma)
        dgamma = (L/m - g*np.cos(gamma) 
                  + v*np.cos(gamma)
                  /(R_EARTH + h)) / v
        dh = v * np.sin(gamma)
        return [dv, dgamma, dh]
    
    def _stop_ground(self, t, y):
        return y[2] - 10e3  # stop at 10 km
    _stop_ground.terminal = True
    
    def run(self, v0, gamma0_deg=-2.0,
            h0=120e3, t_max=800):
        """Integrate trajectory.
        Returns dict with time series.
        """
        gamma0 = np.radians(gamma0_deg)
        sol = solve_ivp(
            self._odes,
            [0, t_max],
            [v0, gamma0, h0],
            method='RK45',
            max_step=0.5,
            events=[self._stop_ground],
            dense_output=True
        )
        
        t = sol.t
        v, gamma, h = sol.y
        
        # Heat flux (Sutton-Graves convective)
        rho = np.array([self.atmo.density(hi) 
                        for hi in h])
        q_conv = (1.7415e-4 
                  * np.sqrt(rho / self.R_nose)
                  * v**3)  # W/m2
        
        # Tauber-Sutton radiative (simplified)
        q_rad = np.where(
            v > 4000,
            4.736e4 * self.R_nose**0.2 
            * rho**1.22 * (v/1e4)**8.5,
            0.0
        )
        
        return {
            't': t, 'v': v, 'gamma': gamma,
            'h': h, 'rho': rho,
            'q_conv': q_conv, 'q_rad': q_rad,
            'q_total': q_conv + q_rad,
            'decel_g': np.gradient(v, t) / G_EARTH
        }
    