# AURORA-MHD Feasibility Results

## v1 → v11 Model Comparison

This repository implements two physics models. The **v1 model** (original) uses uncoupled Faraday extraction at freestream velocity. The **v11 model** (physically corrected) adds post-shock velocity, Hall effect, effective area, and the kinetic ceiling bound. Both are executable and testable — the difference is the physics lesson.

| Parameter | v1 (original) | v11 (corrected) | Why they differ |
|-----------|---------------|-----------------|-----------------|
| Solenoid model | Infinite (B = μ₀nI) | Finite (Fabry F=0.514) | R=0.5m, L=0.6m → F far from 1 |
| Structural ratio | K=2.0 | K=1.0 | v11 uses paper's 1:1 estimate |
| Coil mass (2T) | ~264 kg | ~400–500 kg | Finite correction + lower K |
| Flight system | ~264 kg (no aux) | ~600–800 kg (with 105 kg aux) | Electronics + harness added |
| Flow velocity | v_free (7800 m/s) | v_ps = v_free/5 (1560 m/s) | Extraction at shoulder, not freestream |
| Active area | 20 m² | ~0.4 m² (∫B²dA/B²) | Annular shoulder, not full frontal |
| Hall correction | None | 1/(1+β²) ≈ 0.60 | β_eff=0.82 with N=40 segments |
| Kinetic ceiling | Not modeled | P_max = K×½ṁv²_ps = 456 kW | S >> 1 → Faraday formula invalid |
| COP (cryocooler) | 0.015 (10% Carnot) | 0.002 (3% Carnot) / ride-through | No cryocooler achieves 10% at 20K |
| Demand mode | Continuous cryo | Ride-through (cryo OFF) | REBCO thermal mass absorbs 10-15 min |
| Demand (LEO) | ~1 kW | ~10 kW (ride-through) | Joints + controls + margin |
| Extraction (LEO) | ~36 MW | ~150 kW (after losses) | 200× reduction — all physically motivated |
| Energy margin | 24× (with 1000× ad-hoc loss) | ~4× (physics-based) | v11 replaces ad-hoc with derived losses |

## Key Results (v11)

| Result | Value |
|--------|-------|
| Coil mass (2T, 20K, R=0.5m) | ~400–500 kg |
| Flight system mass | ~600–800 kg (<1% of payload) |
| Energy closure (ride-through) | Positive at B ≥ 1.5T |
| Peak margin (2T, LEO, ride-through) | ~4–15× |
| Kinetic ceiling (LEO) | 456 kW |
| MHD effective zone | S > 10 from ~85 km to ~40 km |
| Test suite | 88+ tests (58 v1 + 30+ v11) |

## What the v1→v11 Progression Shows

The uncoupled Faraday formula gives 36 MW at LEO — but the flow can only deliver ~911 kW of kinetic energy (½ṁv²_ps). The 200× reduction from v1 to v11 breaks down as:

- **Post-shock velocity** (v²): 25× reduction (v_free → v_ps = v_free/5)
- **Active area**: 25× reduction (20 m² → 0.8 m²)
- **Hall effect**: 1.7× reduction (β_eff = 0.82)
- **Channel losses**: 3× reduction (standoff + sheath)

These are multiplicative: 25 × 25 × 1.7 × 3 ≈ 3,200× — but the kinetic ceiling binds first at ~80×. After the ceiling and channel losses, net extraction is ~150 kW vs ~10 kW demand: margin ~4–15×.

## Deflection vs Extraction

Two physically distinct questions, answered with different confidence:

**Deflection is robust.** Stuart numbers S = 1,300 (LEO) to 11,200 (Mars) guarantee strong plasma diversion. MEESST demonstrated 83% heat flux reduction at S ≈ 3. This depends only on σB²L/(ρv) and is insensitive to extraction details.

**Self-sustaining extraction is feasible but marginal.** At ~4× margin with conservative assumptions, a factor-of-2 error in v_ps or channel efficiency is the difference between feasible and infeasible. Coupled MHD-CFD simulation is the critical next step.

## Deliverables

| Deliverable | Output |
|---|---|
| D1: Magnet Trade | d1_trade.pdf — Pareto front + v11 trade table |
| D2: Energy Closure | d2_balance.pdf, d2_v11_margins.pdf — v1 vs v11 comparison |
| D3: Operating Envelope | d3_envelope.pdf — σ + Stuart number zones |

## Limitations (v11)

- Lumped-parameter model: uniform σ, v_ps across shock layer
- Post-shock velocity v_ps = v_free/5 is empirical (Rankine-Hugoniot gives v_free/11 at stagnation; cross-flow at shoulder gives v_free/5)
- Kinetic ceiling is conservative: thermal-to-kinetic reconversion may increase effective ceiling
- Channel efficiency (3–6×) bounds are estimated, not computed
- No electrode erosion model
- No 3D field topology (extraction assumes annular shoulder region)
- No non-equilibrium plasma chemistry (Saha = LTE assumption)

## Conclusion

MHD thermal protection is **energetically feasible but engineering-marginal**. Deflection is robust at all reentry conditions. Self-sustaining extraction is positive at B ≥ 1.5T with ~4× margin at conservative assumptions. The dominant uncertainty is the post-shock cross-flow velocity (v_ps enters as v²) and channel efficiency — both require coupled MHD-CFD to resolve.
