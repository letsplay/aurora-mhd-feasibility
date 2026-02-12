# src/magnets.py
import numpy as np
from scipy.special import ellipk, ellipe

MU_0 = 4 * np.pi * 1e-7  # T·m/A

def b_field_loop(a, I, rho, z):
    """B-field of a single circular loop 
    of radius a, current I.
    
    Args:
        a: loop radius [m]
        I: current [A]
        rho: cylindrical radial distance [m]
            (array OK)
        z: axial distance from loop plane [m]
            (array OK)
    
    Returns:
        B_rho, B_z in Tesla (same shape 
        as rho, z)
    """
    rho = np.asarray(rho, dtype=float)
    z = np.asarray(z, dtype=float)
    
    alpha2 = (a + rho)**2 + z**2
    beta2 = (a - rho)**2 + z**2
    k2 = 1 - beta2 / alpha2
    # Clamp k2 to avoid numerical issues
    k2 = np.clip(k2, 0, 1 - 1e-15)
    
    K = ellipk(k2)
    E = ellipe(k2)
    
    C = MU_0 * I / (2 * np.pi)
    sqrt_alpha2 = np.sqrt(alpha2)
    
    B_z = C / sqrt_alpha2 * (
        K + (a**2 - rho**2 - z**2) 
        / beta2 * E
    )
    
    # B_rho: handle rho=0 (on axis)
    B_rho = np.where(
        rho < 1e-12,
        0.0,
        C * z / (rho * sqrt_alpha2) * (
            -K + (a**2 + rho**2 + z**2) 
            / beta2 * E
        )
    )
    return B_rho, B_z

def b_field_coil(loops, eval_points):
    """Total B-field from multiple coaxial
    loops at evaluation points.
    
    Args:
        loops: list of (z_pos, radius, current)
        eval_points: (N,3) array of (x,y,z)
    
    Returns:
        (N,3) array of (Bx, By, Bz)
    """
    pts = np.asarray(eval_points)
    rho = np.sqrt(pts[:,0]**2 + pts[:,1]**2)
    z_pts = pts[:, 2]
    phi = np.arctan2(pts[:,1], pts[:,0])
    
    B_rho_total = np.zeros(len(pts))
    B_z_total = np.zeros(len(pts))
    
    for z_loop, a, I in loops:
        Br, Bz = b_field_loop(
            a, I, rho, z_pts - z_loop
        )
        B_rho_total += Br
        B_z_total += Bz
    
    # Convert cylindrical to Cartesian
    Bx = B_rho_total * np.cos(phi)
    By = B_rho_total * np.sin(phi)
    
    return np.stack([Bx, By, B_z_total], 
                    axis=-1)

def b_on_axis(a, I, z):
    """Analytic on-axis field for validation.
    B_z = mu0*I*a^2 / (2*(a^2+z^2)^1.5)
    """
    return (MU_0 * I * a**2 
            / (2 * (a**2 + z**2)**1.5))

# src/magnets.py (continued)
from scipy.interpolate import (
    RegularGridInterpolator
)

class REBCOTape:
    """REBCO tape Jc(B,T) model based on 
    SuperPower SCS4050-AP published data.
    
    Tape specs:
      width: 4mm or 12mm
      thickness: ~100 um total
      substrate: Hastelloy (density 8890 kg/m3)
    """
    
    # Digitized Ic(B) at various T [A/cm-width]
    # B in Tesla, Ic per cm of tape width
    _T_data = np.array([20, 40, 65, 77])  # K
    _B_data = np.array(
        [0, 1, 2, 3, 5, 8, 10, 15]
    )  # T
    # Ic [A/cm-w] — worst case (B perp)
    _Ic_table = np.array([
        [550, 420, 350, 300, 230, 170, 145, 100],
        [400, 290, 230, 190, 140, 95, 78, 48],
        [230, 140, 100, 78, 50, 30, 22, 10],
        [150, 80, 52, 38, 20, 10, 6, 2],
    ])  # shape (4_T, 8_B)
    
    TAPE_THICKNESS = 100e-6   # m
    TAPE_DENSITY = 8890       # kg/m3 (Hastelloy)
    
    def __init__(self, width_mm=12):
        self.width = width_mm * 1e-3  # m
        self.cross_section = (
            self.width * self.TAPE_THICKNESS
        )
        
        # Build interpolator for Jc [A/m2]
        Jc_table = (
            self._Ic_table * 1e2  
            # A/cm-w → A/m-width
            / self.TAPE_THICKNESS  # → A/m2
        )
        self._interp = RegularGridInterpolator(
            (self._T_data, self._B_data),
            Jc_table,
            method='linear',
            bounds_error=False,
            fill_value=None
        )
    
    def Jc(self, B, T):
        """Critical current density [A/m2].
        B: magnetic field [T]
        T: temperature [K]
        """
        B = np.atleast_1d(np.clip(np.asarray(B, dtype=float), 0, 15))
        T = np.atleast_1d(np.clip(np.asarray(T, dtype=float), 20, 77))
        B, T = np.broadcast_arrays(B, T)
        pts = np.column_stack([T.ravel(), B.ravel()])
    
        result = self._interp(pts)
        return np.maximum(result, 0.0)
    
    def Ic(self, B, T):
        """Critical current [A] for one tape.
        """
        return self.Jc(B, T) * self.cross_section
    
    def mass_per_meter(self):
        """Tape mass per meter [kg/m]."""
        return (self.cross_section 
                * self.TAPE_DENSITY)
    
    # src/magnets.py (continued)
from scipy.optimize import brentq

class CoilDesigner:
    """Design an HTS coil system and 
    compute total mass."""
    
    K_STRUCT = 2.0     # structural multiplier
    CRYO_SURFACE = 5.0 # kg/m2 cryostat
    CRYO_HEAD = 15.0   # kg cryocooler
    
    def __init__(self, tape=None):
        self.tape = tape or REBCOTape()
    
    def design_solenoid(self, B_target, radius,
                        T_op=20.0, length=0.5):
        """Design a solenoid coil that produces
        B_target [T] at center.
        
        Returns dict with mass breakdown 
        or None if infeasible.
        """
        tape = self.tape
        
        # Solenoid: B_center ≈ mu0 * n * I
        # where n = N/length, I = current
        # Tape carries I <= Ic(B_peak, T_op)
        # B_peak ≈ B_target for a solenoid
        
        Ic_at_target = float(
            np.squeeze(tape.Ic(B_target, T_op))
        )
        if Ic_at_target < 10.0:
            return None  # infeasible
        
        I_op = 0.8 * Ic_at_target  # 80% margin
        n = B_target / (MU_0 * I_op)  # turns/m
        N_turns = int(np.ceil(n * length))
        
        # Tape length and mass
        circumference = 2 * np.pi * radius
        tape_length = N_turns * circumference
        m_tape = tape_length * tape.mass_per_meter()
        
        # Structural mass
        m_struct = m_tape * self.K_STRUCT
        
        # Cryostat: cylinder around coil
        cryo_area = (2 * np.pi * (radius + 0.05) 
                     * (length + 0.1))
        m_cryo = cryo_area * self.CRYO_SURFACE
        
        m_total = (m_tape + m_struct 
                   + m_cryo + self.CRYO_HEAD)
        
        return {
            'B_target': B_target,
            'radius': radius,
            'T_op': T_op,
            'N_turns': N_turns,
            'I_op': I_op,
            'm_tape': m_tape,
            'm_struct': m_struct,
            'm_cryo': m_cryo,
            'm_total': m_total,
            'Ic_margin': I_op / Ic_at_target,
            'feasible': True,
        }
        # src/magnets.py (continued)
import pandas as pd
from itertools import product

def run_trade_sweep(
    B_range=np.arange(0.5, 5.5, 0.5),
    R_range=np.arange(0.5, 4.5, 0.5),
    T_range=np.array([20, 30, 40, 50, 65]),
    save_path='results/trade_space.parquet'
):
    """Run parametric sweep over coil 
    design space.
    
    Returns: pd.DataFrame with all 
    configurations.
    """
    designer = CoilDesigner()
    records = []
    
    for B, R, T in product(
        B_range, R_range, T_range
    ):
        result = designer.design_solenoid(
            B_target=float(B),
            radius=float(R),
            T_op=float(T)
        )
        if result is not None:
            records.append(result)
        else:
            records.append({
                'B_target': float(B),
                'radius': float(R),
                'T_op': float(T),
                'm_total': np.nan,
                'feasible': False,
            })
    
    df = pd.DataFrame(records)
    if save_path:
        df.to_parquet(save_path, index=False)
    return df

def pareto_front(df):
    """Extract Pareto-optimal configs:
    max B for given mass."""
    feasible = df[df['feasible']].copy()
    feasible = feasible.sort_values('m_total')
    
    pareto = []
    best_B = -1
    for _, row in feasible.iterrows():
        if row['B_target'] > best_B:
            pareto.append(row)
            best_B = row['B_target']
    
    return pd.DataFrame(pareto)
