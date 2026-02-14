# src/physics_v11.py
"""
AURORA-MHD v11 Physics Corrections
===================================
This module contains all corrections that align the repo
with the v11 feasibility paper. Each function wraps or
extends a v1 calculation with physically motivated fixes.

Key corrections:
  D5.1 — Finite solenoid (Fabry factor)
  D5.2 — Post-shock velocity (v_ps = v_free/5)
  D5.3 — Hall parameter + effective area
  D5.4 — Kinetic ceiling (P_max = K × ½ṁv²_ps)

The v1 functions remain untouched as regression baselines.
"""
import numpy as np

# ─────────────────────────────────────────────────────────
# D5.1 — Finite Solenoid Correction
# ─────────────────────────────────────────────────────────

MU_0 = 4 * np.pi * 1e-7  # T·m/A


def fabry_factor(R, L):
    """Fabry correction for finite solenoid.

    B_center(finite) = B_center(infinite) × F
    where F = L / sqrt(L² + 4R²)

    For R=0.5m, L=0.6m: F ≈ 0.514 (need ~2× more turns).

    Args:
        R: coil radius [m]
        L: coil length [m]

    Returns:
        F: dimensionless correction factor (0, 1)
    """
    return L / np.sqrt(L**2 + 4 * R**2)


def solenoid_turns_v11(B_target, R, L, I_op):
    """Number of turns needed for B_target at center
    of a finite solenoid.

    B = mu0 × (N/L) × I × F  →  N = B×L / (mu0 × I × F)

    Args:
        B_target: desired center field [T]
        R: coil radius [m]
        L: coil length [m]
        I_op: operating current [A]

    Returns:
        N_turns: integer turn count
    """
    F = fabry_factor(R, L)
    n = B_target / (MU_0 * I_op * F)  # turns per meter
    return int(np.ceil(n * L))


# ─────────────────────────────────────────────────────────
# D5.2 — Post-Shock Velocity Model
# ─────────────────────────────────────────────────────────

# Physics note:
# Normal shock (γ=1.2): v_ps/v_free = (γ-1)/(γ+1) = 0.091 → v_free/11
# But MHD extraction occurs at the shoulder, not stagnation point.
# Tangential (cross-flow) velocity ~ v_free × sin(θ) where
# θ ~ 30°–60° over the shoulder region.
# Angular average: v_cross ~ v_free/5 (v11 baseline).
# Table 5 in paper shows sensitivity: v_free/3 to v_free/8.
COMPRESSION_RATIO = 5.0  # post-shock density / freestream


def post_shock_velocity(v_free, ratio=5.0):
    """Post-shock cross-flow velocity for MHD extraction.

    This is the effective velocity driving Faraday EMF in the
    shoulder region of the shock layer, NOT the normal shock
    velocity at stagnation.

    Args:
        v_free: freestream velocity [m/s]
        ratio: v_free/v_ps ratio (default 5.0, paper baseline)

    Returns:
        v_ps: post-shock cross-flow velocity [m/s]
    """
    return v_free / ratio


def mass_flux_through_channel(v_free, h, atmo, R_body=4.5):
    """Mass flux through the MHD extraction region.

    ṁ = ρ_free × v_free × A_capture

    A_capture = π × R_body² (frontal area captures flow
    that enters the shock layer)

    Args:
        v_free: freestream velocity [m/s]
        h: altitude [m]
        atmo: atmosphere object with .density(h) method
        R_body: vehicle nose radius [m]

    Returns:
        m_dot: mass flow rate [kg/s]
    """
    rho = atmo.density(h)
    A_cap = np.pi * R_body**2
    return rho * v_free * A_cap


# ─────────────────────────────────────────────────────────
# D5.3 — Hall Parameter & Effective Area
# ─────────────────────────────────────────────────────────

def hall_correction(beta_eff=0.82):
    """Effective extraction reduction due to Hall effect.

    At β_eff = 0.82 (paper baseline for N=40 electrode pairs),
    correction factor = 1/(1 + β²) ≈ 0.598

    Segmented electrodes partially compensate the Hall effect;
    β_eff accounts for this (β_bulk would be higher).

    Args:
        beta_eff: effective Hall parameter (0 = no Hall, >1 = strong)

    Returns:
        correction factor in (0, 1]
    """
    return 1.0 / (1.0 + beta_eff**2)


def compute_B2dA(B_center, R, fill_factor=0.5):
    """Estimate ∫B²dA over the extraction region.

    Conservative hand estimate:
    ∫B²dA ≈ B²_center × A_frontal × fill_factor

    Biot-Savart integration gives ~2.0 T²·m² at 2T;
    this estimate gives ~1.35 T²·m² (deliberately conservative).

    Args:
        B_center: on-axis field [T]
        R: body radius [m]
        fill_factor: fraction of frontal area with effective B
                     (0.5 = conservative, accounts for field
                     decay away from axis)

    Returns:
        B2dA: field-area integral [T²·m²]
    """
    A_frontal = np.pi * R**2
    return B_center**2 * A_frontal * fill_factor


def effective_area(B2dA, B_center):
    """Effective MHD extraction area.

    A_eff = ∫B²dA / B²_center

    This is the area-weighted average where the field is
    strong enough for extraction. Typically 0.4–1.0 m² for
    R=0.5m, much less than the full frontal area (~20 m²
    for R=4.5m body).

    The extraction zone is an annular region around the
    shoulder where v ⊥ B (not the stagnation point where
    v → 0 and B ∥ v).

    Args:
        B2dA: field-area integral [T²·m²]
        B_center: on-axis field [T]

    Returns:
        A_eff: effective area [m²]
    """
    if B_center <= 0:
        return 0.0
    return B2dA / B_center**2


# ─────────────────────────────────────────────────────────
# D5.4 — Kinetic Ceiling & Stuart Number
# ─────────────────────────────────────────────────────────

def stuart_number(sigma, B, v_ps, L_char, rho_ps):
    """Stuart number (magnetic interaction parameter).

    S = σ B² L / (ρ v)

    S >> 1: MHD dominates flow (deflection robust)
    S ~ 1:  MHD significant
    S << 1: MHD negligible

    Args:
        sigma: conductivity [S/m]
        B: magnetic field [T]
        v_ps: post-shock velocity [m/s]
        L_char: interaction length [m] (≈ coil length)
        rho_ps: post-shock density [kg/m³]

    Returns:
        S: dimensionless Stuart number
    """
    if rho_ps <= 0 or v_ps <= 0:
        return 0.0
    return sigma * B**2 * L_char / (rho_ps * v_ps)


def kinetic_ceiling(m_dot, v_ps, K=0.5):
    """Maximum extractable power bounded by post-shock
    kinetic energy flux.

    P_max = K × ½ṁv²_ps

    At S >> 1, the MHD force decelerates the flow. You cannot
    extract more energy than the kinetic energy carried by the
    flow through the extraction volume.

    K = 0.5 (conservative): extract half the kinetic flux.
    K = 0.9 (achievable): still satisfies full-deceleration
        condition S(1-K)δ/L ≈ 13.

    Args:
        m_dot: mass flow rate through channel [kg/s]
        v_ps: post-shock cross-flow velocity [m/s]
        K: extraction fraction (0 to ~0.95)

    Returns:
        P_max: kinetic ceiling [W]
    """
    return K * 0.5 * m_dot * v_ps**2


# ─────────────────────────────────────────────────────────
# Complete v11 Extraction Pipeline
# ─────────────────────────────────────────────────────────

def faraday_power_v11(sigma, v_free, B, h, atmo,
                      B2dA=None, R_coil=0.5, L_coil=0.6,
                      beta_eff=0.82, channel_eff=3.0,
                      K_ceiling=0.5, R_body=4.5,
                      v_ps_ratio=5.0):
    """Complete v11 Faraday power extraction pipeline.

    Steps:
      1. v_ps = v_free / v_ps_ratio
      2. σ already provided (from Saha model)
      3. P_F_raw = σ × (v_ps × B)² × 0.25 × δ × A_eff
         (with Hall correction and effective area)
      4. P_max = K × ½ṁv²_ps (kinetic ceiling)
      5. P_raw = min(P_F, P_max)
      6. P_extract = P_raw / channel_eff

    Args:
        sigma: conductivity [S/m]
        v_free: freestream velocity [m/s]
        B: magnetic field [T]
        h: altitude [m]
        atmo: atmosphere object
        B2dA: pre-computed ∫B²dA [T²·m²] (or None to compute)
        R_coil: coil radius [m]
        L_coil: coil length [m]
        beta_eff: effective Hall parameter
        channel_eff: total channel loss factor (1.5–6×)
        K_ceiling: kinetic extraction fraction
        R_body: vehicle nose radius [m]
        v_ps_ratio: v_free/v_ps

    Returns:
        P_extract: net extractable power [W]
    """
    if sigma <= 0 or v_free <= 0 or B <= 0:
        return 0.0

    # Step 1: Post-shock velocity
    v_ps = post_shock_velocity(v_free, ratio=v_ps_ratio)

    # Step 2: Field-area integral
    if B2dA is None:
        B2dA = compute_B2dA(B, R_coil)
    A_eff = effective_area(B2dA, B)

    # Step 3: Raw Faraday power with Hall correction
    delta = 0.05  # shock layer thickness [m]
    K_load = 0.5  # optimal loading factor
    p_density = sigma * (v_ps * B)**2 * K_load * (1 - K_load)
    P_F = p_density * delta * A_eff * hall_correction(beta_eff)

    # Step 4: Kinetic ceiling
    m_dot = mass_flux_through_channel(
        v_free, h, atmo, R_body=R_body)
    P_max = kinetic_ceiling(m_dot, v_ps, K=K_ceiling)

    # Step 5: Binding constraint
    P_raw = min(P_F, P_max)

    # Step 6: Channel losses
    P_extract = P_raw / channel_eff

    return max(P_extract, 0.0)
