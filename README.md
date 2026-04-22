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
python src/build_model_v42_mycn.py
python simulations/validate_v42.py
python simulations/sim_ezh2_feedback_impact.py
```

Use system `python3` (tellurium can be fragile inside some virtualenvs).

## Drug implementations

| Drug | Target | Model action |
|------|--------|--------------|
| HHi (vismodegib) | SMO | `Smo_activation · (1 − GDC0449)` |
| EZH2i (tazemetostat) | EZH2 methyltransferase | Removes EZH2 repression term |
| CDK4/6i (palbociclib) | CDK4 kinase | Sets `V1 = 0` (blocks pRB phosphorylation by CycD/CDK4) |

## Citation

If you use this model, please cite Gérard & Goldbeter (2009) *PNAS* for the cell cycle engine and Chahin et al. (in preparation) for the EZH2/CycD1/MYCN/HH integration.
