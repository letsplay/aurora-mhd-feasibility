# AURORA-MHD Feasibility Study

**Can magnetohydrodynamic thermal protection be self-sustaining for Starship reentry? I built the models to find out.**

![Operating Envelope](results/figures/d3_envelope.png)

## Key Results

| Result | Value |
|--------|-------|
| Minimum magnet mass for 2T field | 264 kg (REBCO @ 20K, R=0.5m) |
| Energy closure | ✓ Self-sustaining at all B ≥ 0.5T |
| Safety margin (with 1000× loss factor) | 24× (LEO), 58× (Lunar), 77× (Mars) |
| MHD effective zone | v > 3 km/s, h > 40 km |

## What This Is

A physics-based feasibility study for replacing Starship's ~18,000 ceramic heat shield tiles with superconducting magnets that deflect reentry plasma via magnetohydrodynamic (MHD) interaction. The study covers:

- **D1**: REBCO magnet mass-field trade space (Biot-Savart + SuperPower tape model)
- **D2**: Time-resolved energy balance — can MHD extraction power the system?
- **D3**: Operating envelope — where does MHD TPS work (Earth vs Mars)?

## How to Run
```bash
git clone https://github.com/letsplay/aurora-mhd-feasibility.git
cd aurora-mhd-feasibility
pip install -r requirements.txt
pytest tests/ -v                    # 58 tests
jupyter notebook notebooks/         # open any notebook
```

## Test Suite

58 tests covering geometry, magnetics, atmosphere, trajectory, plasma physics, MHD power extraction, energy balance, and operating envelope.

## Author

Luc — PhD Applied Physics (CEA Saclay), 25+ years tech program management.
Exploring AI-accelerated MHD thermal protection for next-generation spacecraft.

## License

MIT