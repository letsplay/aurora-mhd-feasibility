# src/plasma.py
import numpy as np

K_B = 1.381e-23     # J/K
E_CHARGE = 1.602e-19 # C
M_E = 9.109e-31     # kg
H_PLANCK = 6.626e-34 # J·s
EV_TO_J = 1.602e-19

# Ionization energies [eV]
E_ION = {
    'N': 14.53, 'O': 13.62,
    'CO2': 13.77, 'Cs': 3.89
}

def post_shock_temperature(v, gamma=1.3,
                           M_mol=0.029):
    """Simplified post-shock T from 
    Rankine-Hugoniot. 
    v: freestream velocity [m/s]
    M_mol: mean molecular mass [kg/mol]
    """
    R_gas = 8.314  # J/(mol·K)
    return (2 * (gamma-1) * M_mol * v**2 
            / ((gamma+1)**2 * R_gas))

def saha_ne(n_neutral, T, E_ion_eV,
            g_ratio=1.0):
    """Saha equation: solve for n_e.
    n_e^2 / n_n = f(T).
    Returns electron density [m-3].
    """
    T = np.maximum(T, 300)  # floor
    E_ion = E_ion_eV * EV_TO_J
    
    thermal = (2*np.pi*M_E*K_B*T 
               / H_PLANCK**2)**1.5
    saha_rhs = 2 * g_ratio * thermal * np.exp(
        -E_ion / (K_B * T)
    )
    
    # n_e^2 = saha_rhs * (n_neutral - n_e)
    # ≈ saha_rhs * n_neutral for weak ioniz.
    # Full quadratic: n_e^2 + saha*n_e 
    #                 - saha*n_n = 0
    disc = saha_rhs**2 + 4*saha_rhs*n_neutral
    n_e = (-saha_rhs + np.sqrt(disc)) / 2
    return np.maximum(n_e, 0)

def conductivity(n_e, T, n_neutral,
                 sigma_en=1e-19):
    """Electrical conductivity [S/m].
    Handles both weakly and partially 
    ionized regimes."""
    if np.any(n_e < 1):
        return np.where(n_e > 1,
            conductivity(
                np.maximum(n_e, 1), T,
                n_neutral, sigma_en
            ), 0.0)
    
    v_thermal = np.sqrt(8*K_B*T / (np.pi*M_E))
    nu_en = n_neutral * sigma_en * v_thermal
    sigma_weakly = n_e * E_CHARGE**2 / (M_E * nu_en)
    
    # Spitzer limit for ionized plasma
    ln_lambda = 15.0
    sigma_spitzer = 1.53e-2 * T**1.5 / ln_lambda
    
    return np.minimum(sigma_weakly, sigma_spitzer)

def sigma_at(h, v, atmo, species='air',
             seed_frac=0.0):
    """Full pipeline: h, v → σ.
    atmo: atmosphere object
    species: 'air' or 'co2'
    seed_frac: Cs mass fraction (0-0.01)
    """
    rho = atmo.density(h)
    T = atmo.temperature(h)
    
    if species == 'air':
        M_mol = 0.029
        E_ion_primary = (0.78 * E_ION['N'] 
                        + 0.22 * E_ION['O'])
    else:  # CO2
        M_mol = 0.044
        E_ion_primary = E_ION['CO2']
    
    T_post = post_shock_temperature(
        v, M_mol=M_mol
    )
    n_total = rho / (M_mol / 6.022e23)
    
    # Primary species ionization
    n_e = saha_ne(n_total, T_post,
                  E_ion_primary)
    
    # Cs seeding
    if seed_frac > 0:
        n_cs = seed_frac * n_total
        n_e_cs = saha_ne(n_cs, T_post,
                         E_ION['Cs'])
        n_e = np.sqrt(n_e**2 + n_e_cs**2)
    
    n_neutral = np.maximum(n_total - n_e, 1.0)
    return conductivity(n_e, T_post, n_neutral)
