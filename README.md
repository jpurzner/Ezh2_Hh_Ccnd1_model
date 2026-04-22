# Ezh2_Hh_Ccnd1_model

ODE model of MYCN–EZH2–Cyclin D1–Hedgehog regulation of the mammalian cell cycle, applied to cerebellar granule neuron progenitors (GNPs) and SHH-subgroup medulloblastoma (SHH-MB). This is the v42 model accompanying Chahin et al. on Ezh2/Cyclin D1 negative feedback in GNPs and medulloblastoma.

## Model overview

The model integrates four modules into the Gérard & Goldbeter (2009) mammalian cell cycle engine:

1. **Hedgehog signaling** — SHH → Ptch1 → Smo → Gli_act/Gli_rep → Gli1 → CycD1 transcription
2. **EZH2 module** — E2F- and pRB-hyperphosphorylation–driven transcription; represses CycD1 via `K / (K + EZH2·(1 − EZH2i))`
3. **MYCN module** — autonomous basal synthesis (amplified 2.8× in MB) plus a small Gli-dependent term; drives CycD1 through a Hill function (n = 3), producing a threshold effect
4. **Cell cycle engine (Gérard 2009)** — Rb–E2F bistable switch, p27/Skp2, cyclin D/E/A/B–CDK complexes, Cdc25 phosphatases, APC/C (Cdc20, Cdh1), time-scaled by ε = 150 to produce ~22 h periods

CycD1 mRNA is the integration node:

```
d[Cd_mRNA]/dt = [ k_basal
                  + k_Gli · (Gli_act + Gli1)^2 / (K_Gli^2 + (Gli_act + Gli1)^2)
                                                · K_Gli_rep^2 / (K_Gli_rep^2 + Gli_rep^2)
                  + k_MYCN · MYCN^3 / (K_MYCN^3 + MYCN^3) ]
                · K_EZH2 / (K_EZH2 + EZH2 · (1 − EZH2i))
                − k_deg · Cd_mRNA
```

## Key predictions

- **GNPs**: CycD1 is ~86% Gli-driven; HHi causes near-complete CycD1 loss and cell cycle arrest.
- **MYCN-amplified MB**: MYCN contributes ~41% of CycD1, rising to ~82% under HHi, producing partial HHi resistance.
- **EZH2 brake**: EZH2 feedback lengthens the GNP cycle by ~4.3 h but MB by only ~0.9 h (MYCN saturates the Rb–E2F switch).
- **Rescue specificity**: EZH2i rescues HHi-arrested MB (upstream blockade) but *not* CDK4/6i-arrested MB (downstream blockade), because EZH2i acts by de-repressing CycD1 transcription.

Sixteen of 18 validation targets from GNP/MB RNA-seq and functional assays are recovered. See `simulations/V42_MODEL_SUMMARY.md` for the full calibration/validation table.

## Repository layout

```
src/
  build_model_v42_mycn.py          # Antimony model builder
simulations/
  validate_v42.py                  # 10-condition validation (drug panel, rescue, CDK4/6i)
  sim_ezh2_feedback_impact.py      # EZH2 feedback on/off comparison (168 h)
  sim_ezh2_compensation.py         # EZH2-boost sweep
  sim_ezh2_data_constrained.py     # RNA-seq-anchored EZH2 analysis
  sweep_mycn_hill.py               # MYCN Hill parameter sweep (960 combos)
  fig_v42_architecture.py          # Model wiring diagram
  V42_MODEL_SUMMARY.md             # Construction, calibration, results
  V42_MODEL_EQUATIONS.md           # Full ODE listing
  PAPER_MODEL_SECTIONS.md/.txt     # Results + Methods sections for the paper
  fig_v42_*.png/pdf                # Architecture, validation, rescue, CycB traces
  fig_ezh2_*.png/pdf               # EZH2 feedback impact, periods, compensation
```

## Getting started

```bash
pip install -r requirements.txt
```

Use system `python3` (tellurium can be fragile inside some virtualenvs).

## Reproducing the figures

All scripts must be run **from the repository root** (they use relative paths for the output figures and resolve the `src/` module with `os.path.dirname(__file__)`). Each script regenerates its figures in-place under `simulations/`.

| Figure(s) | Script | Approx. runtime |
|-----------|--------|-----------------|
| `fig_v42_architecture.png/pdf` | `python simulations/fig_v42_architecture.py` | <10 s (schematic, no simulation) |
| `fig_v42_validation.png/pdf`, `fig_v42_cycb_traces.png`, `fig_v42_ezh2i_rescue.png`, `validation_v42_results.json` | `python simulations/validate_v42.py` | ~3–5 min (10 conditions × 96 h + 5-condition 168 h rescue panel) |
| `fig_ezh2_feedback_impact.png/pdf`, `fig_ezh2_feedback_periods.png` | `python simulations/sim_ezh2_feedback_impact.py` | ~2 min (4 × 168 h simulations) |
| `fig_ezh2_compensation.png/pdf` | `python simulations/sim_ezh2_compensation.py` | ~2 min |
| `fig_ezh2_data_constrained.png/pdf` | `python simulations/sim_ezh2_data_constrained.py` | ~5 min |
| `sweep_mycn_hill_results.json` (parameter calibration) | `python simulations/sweep_mycn_hill.py` | ~30–60 min (960 parameter combinations) |

The model itself can be printed/inspected with:

```bash
python src/build_model_v42_mycn.py
```

All simulations use the CVODE stiff integrator via tellurium/roadrunner (default). Because the ODEs are deterministic and CVODE is stable, re-running should produce results numerically identical to the checked-in figures modulo tellurium/roadrunner version differences.

## Drug implementations

| Drug | Target | Model action |
|------|--------|--------------|
| HHi (vismodegib) | SMO | `Smo_activation · (1 − GDC0449)` |
| EZH2i (tazemetostat) | EZH2 methyltransferase | Removes EZH2 repression term |
| CDK4/6i (palbociclib) | CDK4 kinase | Sets `V1 = 0` (blocks pRB phosphorylation by CycD/CDK4) |

## Citation

If you use this model, please cite Gérard & Goldbeter (2009) *PNAS* for the cell cycle engine and Chahin et al. (in preparation) for the EZH2/CycD1/MYCN/HH integration.
