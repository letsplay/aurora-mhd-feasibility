# tests/test_geometry.py
import os
import numpy as np
from src.geometry import BluntBody

def test_surface_area():
    """Hemisphere area = 2*pi*R^2."""
    body = BluntBody(R_nose=4.5, 
                     theta_max=np.pi/2)
    expected = 2 * np.pi * 4.5**2  # ~127.2 m²
    assert abs(body.surface_area() - expected) \
        / expected < 0.001

def test_normals_outward():
    """All normals point outward: 
    dot(normal, position) > 0."""
    body = BluntBody()
    dots = np.sum(
        body.normals * body.points, axis=1
    )
    assert np.all(dots > 0)

def test_rotational_symmetry():
    """Points at same theta but different phi
    have same radius."""
    body = BluntBody(n_theta=20, n_phi=20)
    pts = body.points.reshape(20, 20, 3)
    # radius at each theta (first ring)
    r_ring = np.linalg.norm(pts[5, :, :], 
                            axis=1)
    assert np.std(r_ring) < 1e-10

def test_point_count():
    body = BluntBody(n_theta=50, n_phi=40)
    assert body.points.shape == (2000, 3)
    assert body.normals.shape == (2000, 3)

    # tests/test_magnet.py
import numpy as np
import pytest
from src.magnets import (
    b_field_loop, b_on_axis, 
    b_field_coil, MU_0
)

def test_on_axis_single_loop():
    """On-axis B matches analytic formula."""
    a, I = 1.0, 1000.0
    z_vals = np.linspace(0.01, 3.0, 50)
    _, Bz_numeric = b_field_loop(
        a, I, rho=np.zeros_like(z_vals), z=z_vals
    )
    Bz_analytic = b_on_axis(a, I, z_vals)
    rel_err = np.abs(
        (Bz_numeric - Bz_analytic) / Bz_analytic
    )
    assert np.all(rel_err < 1e-4), \
        f"Max rel error: {rel_err.max():.2e}"

def test_center_field():
    """B at center = mu0*I/(2a)."""
    a, I = 0.5, 500.0
    _, Bz = b_field_loop(a, I, rho=0.0, z=0.0)
    expected = MU_0 * I / (2 * a)
    assert abs(Bz - expected)/expected < 1e-6

def test_helmholtz_uniformity():
    """Helmholtz pair: uniform within 1%
    over central 50% of gap."""
    a, I = 1.0, 1000.0
    d = a  # Helmholtz spacing
    loops = [(-d/2, a, I), (d/2, a, I)]
    z_eval = np.linspace(-d/4, d/4, 20)
    pts = np.column_stack([
        np.zeros(20), np.zeros(20), z_eval
    ])
    B = b_field_coil(loops, pts)
    Bz = B[:, 2]
    uniformity = (Bz.max()-Bz.min()) / Bz.mean()
    assert uniformity < 0.01

def test_b_rho_zero_on_axis():
    """B_rho must be 0 on the axis."""
    a, I = 1.0, 1000.0
    Br, _ = b_field_loop(a, I, rho=0.0, z=0.5)
    assert abs(Br) < 1e-15

def test_computation_speed():
    """1000 points from 100 loops < 2s."""
    import time
    loops = [(i*0.01, 1.0, 100) 
             for i in range(100)]
    pts = np.random.randn(1000, 3)
    t0 = time.time()
    b_field_coil(loops, pts)
    assert time.time() - t0 < 2.0

    # tests/test_magnet.py (continued)
from src.magnets import REBCOTape, CoilDesigner, run_trade_sweep, pareto_front

def test_Ic_77K_self_field():
    """77K self-field Ic ~ 100-180 A/cm-w."""
    tape = REBCOTape(width_mm=12)
    Ic = tape.Ic(B=0.0, T=77.0)
    Ic_per_cm = Ic / (tape.width * 100)
    assert 100 < Ic_per_cm < 200, \
        f"Ic/cm-w = {Ic_per_cm:.0f}"

def test_Ic_20K_5T():
    """20K, 5T perp: Ic ~ 200-300 A/cm-w."""
    tape = REBCOTape(width_mm=12)
    Ic = tape.Ic(B=5.0, T=20.0)
    Ic_per_cm = Ic / (tape.width * 100)
    assert 150 < Ic_per_cm < 350

def test_Jc_decreases_with_B():
    """Jc must decrease monotonically 
    with B at fixed T."""
    tape = REBCOTape()
    B_vals = np.linspace(0, 10, 20)
    Jc_vals = tape.Jc(B_vals, T=40.0)
    assert np.all(np.diff(Jc_vals) <= 0)

def test_Jc_decreases_with_T():
    """Jc must decrease with T at fixed B."""
    tape = REBCOTape()
    T_vals = np.linspace(20, 77, 20)
    Jc_vals = tape.Jc(B=3.0, T=T_vals)
    assert np.all(np.diff(Jc_vals) <= 0)

def test_Jc_never_negative():
    tape = REBCOTape()
    B = np.linspace(0, 15, 50)
    T = np.linspace(20, 77, 50)
    BB, TT = np.meshgrid(B, T)
    Jc = tape.Jc(BB.ravel(), TT.ravel())
    assert np.all(Jc >= 0)

def test_mass_sanity_1T():
    """1T, 1m radius: mass 20-300 kg."""
    cd = CoilDesigner()
    d = cd.design_solenoid(
        B_target=1.0, radius=1.0, T_op=20
    )
    assert d is not None
    assert 20 < d['m_total'] < 300, \
        f"mass={d['m_total']:.0f} kg"

def test_mass_scales_with_B():
    """Doubling B should roughly 4x mass
    (tape ∝ N ∝ B, length ∝ N → m ∝ B²)."""
    cd = CoilDesigner()
    d1 = cd.design_solenoid(1.0, 1.0, 20)
    d2 = cd.design_solenoid(2.0, 1.0, 20)
    ratio = d2['m_total'] / d1['m_total']
    assert 2.0 < ratio < 8.0, \
        f"Mass ratio: {ratio:.1f}"

def test_infeasible_detection():
    """Very high B at 77K should be 
    infeasible (Jc too low)."""
    cd = CoilDesigner()
    d = cd.design_solenoid(
        B_target=10.0, radius=1.0, T_op=77
    )
    assert d is None

def test_cryo_less_than_tape():
    """Cryostat mass < tape+struct mass."""
    cd = CoilDesigner()
    d = cd.design_solenoid(2.0, 1.5, 20)
    assert d['m_cryo'] < (
        d['m_tape'] + d['m_struct']
    ) 
def test_sweep_completes():
    """Small sweep runs in < 30s."""
    import time
    t0 = time.time()
    df = run_trade_sweep(
        B_range=np.arange(1, 4, 1),
        R_range=np.arange(1, 3, 1),
        T_range=np.array([20, 40]),
        save_path=None
    )
    assert time.time() - t0 < 30
    assert len(df) > 0
    assert 'feasible' in df.columns

def test_no_nans_in_required_cols():
    """B_target, radius, T_op never NaN."""
    df = run_trade_sweep(
        B_range=[1,2], R_range=[1,2],
        T_range=[20], save_path=None
    )
    for col in ['B_target','radius','T_op']:
        assert df[col].isna().sum() == 0

def test_feasible_fraction():
    """> 20% configs should be feasible."""
    df = run_trade_sweep(
        B_range=np.arange(0.5, 4, 0.5),
        R_range=np.arange(1, 4, 1),
        T_range=[20, 40], save_path=None
    )
    frac = df['feasible'].mean()
    assert frac > 0.2, f"Only {frac:.0%}"

def test_pareto_monotonic():
    """Pareto front: increasing B requires 
    increasing mass."""
    df = run_trade_sweep(
        B_range=np.arange(0.5, 5, 0.5),
        R_range=[1.5, 2.5],
        T_range=[20], save_path=None
    )
    pf = pareto_front(df)
    if len(pf) > 1:
        assert np.all(
            np.diff(pf['m_total']) >= 0
        )

def test_plots_generated():
    """Verify output files exist and 
    are reasonable size."""
    # (run after notebook execution)
    png = 'results/figures/d1_trade.png'
    pdf = 'results/figures/d1_trade.pdf'
    assert os.path.exists(png)
    assert os.path.exists(pdf)
    # PNG should be 100KB - 5MB
    sz = os.path.getsize(png)
    assert 100_000 < sz < 5_000_000
    # PDF should be vector (usually < 1MB)
    assert os.path.getsize(pdf) < 3_000_000

def test_budget_line_intersects():
    """At least one config under 200 kg 
    at >= 1T."""
    df = run_trade_sweep(
        B_range=np.arange(0.5, 5, 0.5),
        R_range=np.arange(1, 4, 1),
        T_range=[20, 30], save_path=None
    )
    under = df[
        (df['feasible']) & 
        (df['m_total'] <= 200) &
        (df['B_target'] >= 1.0)
    ]
    assert len(under) > 0, \
        "No configs under 200kg at 1T+"
    
def test_plots_generated():
    """Verify output files exist and are reasonable size."""
    import os
    png = 'results/figures/d1_trade.png'
    pdf = 'results/figures/d1_trade.pdf'
    assert os.path.exists(png)
    assert os.path.exists(pdf)
    sz = os.path.getsize(png)
    assert 100_000 < sz < 5_000_000

def test_budget_line_intersects():
    """At least one config under 300 kg at >= 1T."""
    df = run_trade_sweep(
        B_range=np.arange(0.5, 4, 0.5),
        R_range=np.arange(1, 4, 1),
        T_range=[20, 30], save_path=None
    )
    under = df[
        (df['feasible']) &
        (df['m_total'] <= 300) &
        (df['B_target'] >= 1.0)
    ]
    assert len(under) > 0, "No configs under 300kg at 1T+"

def test_notebook_execution():
    """Notebook runs end-to-end."""
    import subprocess
    result = subprocess.run([
        'jupyter', 'nbconvert',
        '--execute', '--to', 'notebook',
        '--output', '/tmp/d1_test.ipynb',
        'notebooks/01_magnet_trade.ipynb'
    ], capture_output=True, timeout=1200)
    assert result.returncode == 0, \
        result.stderr.decode()[-500:]