# src/mhd.py
import numpy as np
import pandas as pd
from src.trajectory import ReentryTrajectory

def faraday_power(sigma, v, B,
                  delta=0.05,
                  A_active=20.0,
                  K=0.5):
    """Faraday MHD power extraction [W].
    
    sigma: conductivity [S/m]
    v: flow velocity [m/s]
    B: magnetic field [T]
    delta: interaction layer thickness [m]
    A_active: electrode area [m2]
    K: loading factor (0-1, optimal=0.5)
    """
    sigma = np.asarray(sigma)
    v = np.asarray(v)
    B = np.asarray(B)
    
    # EMF per unit length = v * B
    # Power density = sigma * (vB)^2 * K(1-K)
    p_density = (sigma * (v * B)**2 
                 * K * (1 - K))
    
    # Total power = density * volume
    P = p_density * delta * A_active
    return np.maximum(P, 0.0)

def optimal_loading_factor():
    """Optimal K = 0.5 maximizes K(1-K)."""
    return 0.5

# src/mhd.py (continued)

class PowerDemand:
    """MHD TPS power demand model."""
    
    def __init__(self, coil_design,
                 n_joints=30,
                 R_joint=30e-9,
                 cryo_cop_at_20K=0.015,
                 Q_leak_per_m2=3.0,
                 P_controls=300):
        self.coil = coil_design
        self.n_joints = n_joints
        self.R_joint = R_joint
        self.cop = cryo_cop_at_20K
        self.Q_leak = Q_leak_per_m2
        self.P_controls = P_controls
    
    def joint_losses(self):
        """Resistive loss at tape splices [W]."""
        I = self.coil['I_op']
        return self.n_joints * self.R_joint * I**2
    
    def cryo_steady(self):
        """Steady-state cryocooler power [W].
        (orbital / pre-reentry)"""
        # Cold surface area estimate
        r = self.coil['radius']
        A_cold = 2 * np.pi * r * 0.5  # cylinder
        Q_total = (A_cold * self.Q_leak 
                   + self.joint_losses())
        return Q_total / self.cop
    
    def cryo_reentry(self, q_aero,
                     penetration=0.01):
        """Cryocooler power during reentry.
        q_aero: local heat flux [W/m2]
        penetration: fraction reaching cryostat
        """
        r = self.coil['radius']
        A_cold = 2 * np.pi * r * 0.5
        Q_aero = q_aero * penetration * A_cold
        Q_total = (A_cold * self.Q_leak 
                   + self.joint_losses()
                   + Q_aero)
        return Q_total / self.cop
    
    def total(self, q_aero=0):
        """Total demand [W]."""
        if q_aero > 0:
            P_cryo = self.cryo_reentry(q_aero)
        else:
            P_cryo = self.cryo_steady()
        return P_cryo + self.P_controls
    
    # src/mhd.py (continued)

def energy_balance(traj_result, B_field,
                   atmo, coil_design,
                   seed_frac=0.0):
    """Compute energy balance over trajectory.
    
    Returns dict with time series of
    P_extract, P_demand, P_net, E_cum.
    """
    from src.plasma import sigma_at
    
    t = traj_result['t']
    v = traj_result['v']
    h = traj_result['h']
    q = traj_result['q_total']
    
    pd = PowerDemand(coil_design)
    
    P_ext = np.zeros_like(t)
    P_dem = np.zeros_like(t)
    
    for i in range(len(t)):
        sig = sigma_at(
            h[i], v[i], atmo,
            seed_frac=seed_frac
        )
        P_ext[i] = faraday_power(
            sig, v[i], B_field
        )
        P_dem[i] = pd.total(q_aero=q[i])
    
    P_net = P_ext - P_dem
    # Cumulative energy [J]
    E_cum = np.cumsum(
        P_net[:-1] * np.diff(t)
    )
    E_cum = np.insert(E_cum, 0, 0)
    
    return {
        't': t, 'P_extract': P_ext,
        'P_demand': P_dem, 'P_net': P_net,
        'E_cumulative': E_cum,
        'self_sustaining': np.all(E_cum >= 
            -50e3 * 3600),  # 50 kWh battery
        'min_E_cum': E_cum.min(),
    }

def sweep_closure(atmo, missions, 
                  B_values, coil_designer):
    """Sweep B × missions → result table."""
    results = []
    for name, v0, gamma in missions:
        traj = ReentryTrajectory(atmo)
        r = traj.run(v0=v0, gamma0_deg=gamma)
        for B in B_values:
            cd = coil_designer.design_solenoid(
                B, radius=2.0, T_op=20
            )
            if cd is None:
                results.append({
                    'mission': name, 'B': B,
                    'feasible': False
                })
                continue
            eb = energy_balance(r, B, atmo, cd)
            results.append({
                'mission': name, 'B': B,
                'self_sustaining': 
                    eb['self_sustaining'],
                'min_energy_J': eb['min_E_cum'],
                'mass_kg': cd['m_total'],
                'feasible': True,
            })
    return pd.DataFrame(results)