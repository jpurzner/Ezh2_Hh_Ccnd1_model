# V42 Model Equations — MYCN-EZH2-CyclinD1-Hedgehog Cell Cycle

## Module 1: Hedgehog Signaling

### Ptch1 receptor dynamics
```
d[Ptch1_mRNA]/dt = k_Ptch1_tx * Ptch1_copy_number - k_Ptch1_mRNA_deg * Ptch1_mRNA

d[Ptch1_free]/dt = k_Ptch1_translation * Ptch1_mRNA
                   - k_SHH_Ptch_bind * SHH * Ptch1_free
                   + k_SHH_Ptch_release * SHH_Ptch
                   - k_Ptch1_deg * Ptch1_free

d[SHH_Ptch]/dt = k_SHH_Ptch_bind * SHH * Ptch1_free
                 - k_SHH_Ptch_release * SHH_Ptch
                 - k_SHH_Ptch_deg * SHH_Ptch
```

### Smoothened activation (GDC0449 target)
```
d[Smo_active]/dt = k_Smo_act / (1 + Ptch1_free / K_Ptch_Smo) * (1 - GDC0449)
                   - k_Smo_inact * Smo_active
```

### Gli processing (repressor/activator switch)
```
d[Gli_act]/dt = k_Gli_rep_to_act * Smo^2 / (K_Smo_Gli^2 + Smo^2) * Gli_rep
               - k_Gli_act_to_rep * (1 - Smo^2 / (K_Smo_Gli^2 + Smo^2)) * Gli_act

d[Gli_rep]/dt = -d[Gli_act]/dt  (interconversion)
```

### Gli1 amplification (positive feedback)
```
d[Gli1_mRNA]/dt = Vmax_Gli1_tx * Gli_act^2 / (K_Gli_act^2 + Gli_act^2)
                  * K_Gli_rep^2 / (K_Gli_rep^2 + Gli_rep^2)
                  - k_Gli1_mRNA_deg * Gli1_mRNA

d[Gli1]/dt = k_Gli1_translation * Gli1_mRNA - k_Gli1_deg * Gli1
```

---

## Module 1.5: EZH2 Regulation (E2F/pRBpp-dependent)

```
d[EZH2_mRNA]/dt = k_EZH2_basal + k_EZH2_E2F * E2F / (K_E2F + E2F) * pRBpp / (K_pRBpp + pRBpp)
                  - k_EZH2_mRNA_deg * EZH2_mRNA

d[EZH2]/dt = k_EZH2_translation * EZH2_mRNA - k_EZH2_deg * EZH2
```

EZH2 repression of CycD1 (acts on all CycD1 transcription):
```
EZH2_repression_factor = K_EZH2 / (K_EZH2 + EZH2 * (1 - EZH2i))
```
where EZH2i = 0 (no drug) or 1 (full inhibition).

---

## Module 1.75: MYCN (GDC-resistant CycD1 driver)

```
d[MYCN]/dt = k_MYCN_basal * MYCN_amplification
             + k_MYCN_Gli * (Gli_act + Gli1) / (K_Gli_MYCN + Gli_act + Gli1)
             - k_MYCN_deg * MYCN
```
- `MYCN_amplification` = 1.0 (WT) or 2.8 (MB)
- Basal component is amplified in MB (autonomous, GDC-resistant)
- Gli component is NOT amplified (~34% of basal in WT)

---

## CycD1 mRNA Integration (central node)

CycD1 transcription integrates three inputs, all subject to EZH2 repression:

```
d[Cd_mRNA]/dt = [ k_basal
                  + k_Gli * (Gli_act + Gli1)^n / (K_Gli^n + (Gli_act + Gli1)^n)
                         * K_Gli_rep^n / (K_Gli_rep^n + Gli_rep^n)
                  + k_MYCN * MYCN^3 / (K_MYCN^3 + MYCN^3) ]        ← Hill n=3
                * K_EZH2 / (K_EZH2 + EZH2 * (1 - EZH2i))          ← EZH2 repression
                - k_Cd_mRNA_deg * Cd_mRNA
```

| Term | Description | Parameters |
|------|-------------|------------|
| Basal | Constitutive low-level | k_basal = 0.3 |
| Gli-driven | SHH pathway-dependent | k_Gli = 4.0, K_Gli = 0.4, n = 2 |
| MYCN-driven | Hill function, threshold | k_MYCN = 15.0, K_MYCN = 1.5, n = 3 |
| EZH2 repression | H3K27me3 at CycD1 promoter | K_EZH2 = 0.5 |

---

## Module 2: Cell Cycle Engine (Gerard 2009)

All cell cycle reactions are scaled by `eps = 150` (time-scaling factor).

### Rb-E2F Bistable Switch

```
d[pRB]/dt = vsprb - kdprb*pRB
            - V1*(pRB/(K1+pRB))*(Md+Mdp27) + V2*(pRBp/(K2+pRBp))
            - kpc1*pRB*E2F + kpc2*pRBc1

d[pRBp]/dt = V1*(pRB/(K1+pRB))*(Md+Mdp27) - V2*(pRBp/(K2+pRBp))
             - V3*(pRBp/(K3+pRBp))*Me + V4*(pRBpp/(K4+pRBpp))
             - kdprbp*pRBp - kpc3*pRBp*E2F + kpc4*pRBc2

d[pRBpp]/dt = V3*(pRBp/(K3+pRBp))*Me - V4*(pRBpp/(K4+pRBpp)) - kdprbpp*pRBpp

d[E2F]/dt = vse2f - kde2f*E2F
            - V1e2f*(E2F/(K1e2f+E2F))*Ma + V2e2f*(E2Fp/(K2e2f+E2Fp))
            - kpc1*pRB*E2F + kpc2*pRBc1 - kpc3*pRBp*E2F + kpc4*pRBc2

d[E2Fp]/dt = V1e2f*(E2F/(K1e2f+E2F))*Ma - V2e2f*(E2Fp/(K2e2f+E2Fp)) - kde2fp*E2Fp
```

### Cyclin D / CDK4 (G1 phase)

```
d[Cd]/dt = k_Cd_translation * Cd_mRNA - Vdd*(Cd/(Kdd+Cd)) - kddd*Cd
           - kcom1*Cd*(Cdk4_tot - Mdi - Md - Mdp27) + kdecom1*Mdi

d[Mdi]/dt = kcom1*Cd*(Cdk4_free) - kdecom1*Mdi
            - Vm1d*(Mdi/(K1d+Mdi)) + Vm2d*(Md/(K2d+Md))

d[Md]/dt = Vm1d*(Mdi/(K1d+Mdi)) - Vm2d*(Md/(K2d+Md))
           - kc1*Md*p27 + kc2*Mdp27

d[Mdp27]/dt = kc1*Md*p27 - kc2*Mdp27
```

### Cyclin E / CDK2 (G1/S transition)

```
d[Ce]/dt = kce*E2F*(Ki9/(Ki9+pRB))*(Ki10/(Ki10+pRBp))
           - Vde*(Skp2/(Kdceskp2+Skp2))*(Ce/(Kde+Ce)) - kdde*Ce
           - kcom2*Ce*Cdk2_free + kdecom2*Mei

d[Me]/dt = Vm1e*(Mei/(K1e+Mei))*Pe - Vm2e*ib1*(Me/(K2e+Me))
           - kc3*Me*p27 + kc4*Mep27
```

### Cyclin A / CDK2 (S/G2 phase)

```
d[Ca]/dt = kca*E2F*(Ki11/(Ki11+pRB))*(Ki12/(Ki12+pRBp))
           - Vda*(Ca/(Kda+Ca))*(Cdc20a/(Kacdc20+Cdc20a)) - kdda*Ca

d[Ma]/dt = Vm1a*(Mai/(K1a+Mai))*Pa - Vm2a*ib2*(Ma/(K2a+Ma))
           - kc5*Ma*p27 + kc6*Map27
```

### Cyclin B / CDK1 (M phase)

```
d[Cb]/dt = vcb - Vdb*(Cb/(Kdb+Cb))*(Cdc20a/(Kdbcdc20+Cdc20a))
               - Vdb*(Cb/(Kdb+Cb))*(Cdh1a/(Kdbcdh1+Cdh1a)) - kddb*Cb

d[Mb]/dt = Vm1b*(Mbi/(K1b+Mbi))*Pb - Vm2b*ib3*(Mb/(K2b+Mb))
           - kc7*Mb*p27 + kc8*Mbp27
```

### p27/Kip1 (CDK inhibitor — G0 arrest barrier)

```
d[p27]/dt = vs1p27 + vs2p27*E2F*(Ki13/(Ki13+pRB))*(Ki14/(Ki14+pRBp))
            - V1p27*(p27/(K1p27+p27))*Me + V2p27*(p27p/(K2p27+p27p))
            - kddp27*p27
            - kc1*Md*p27 + kc2*Mdp27 - kc3*Me*p27 + kc4*Mep27
            - kc5*Ma*p27 + kc6*Map27 - kc7*Mb*p27 + kc8*Mbp27

d[p27p]/dt = V1p27*(p27/(K1p27+p27))*Me - V2p27*(p27p/(K2p27+p27p))
             - Vdp27p*(Skp2/(Kdp27skp2+Skp2))*(p27p/(Kdp27p+p27p)) - kddp27p*p27p
```

### Skp2 (E2F target, p27p degradation)

```
d[Skp2]/dt = vsskp2 - Vdskp2*(Skp2/(Kdskp2+Skp2))*(Cdh1a/(Kcdh1+Cdh1a)) - kddskp2*Skp2
```

### APC/C co-activators

```
d[Cdh1a]/dt = vscdh1a + V1cdh1*(Cdh1i/(K1cdh1+Cdh1i))
              - V2cdh1*(Cdh1a/(K2cdh1+Cdh1a))*(Ma+Mb) - kdcdh1a*Cdh1a

d[Cdc20a]/dt = Vm3b*(Cdc20i/(K3b+Cdc20i))*Mb - Vm4b*(Cdc20a/(K4b+Cdc20a)) - kdcdc20a*Cdc20a
```

### Cdc25 Phosphatases (CDK activators)

```
d[Pe]/dt = Vm5e*(Me+ae)*(Pei/(K5e+Pei)) - V6e*xe1*(Pe/(K6e+Pe)) - kdpe*Pe
d[Pa]/dt = Vm5a*(Ma+aa)*(Pai/(K5a+Pai)) - V6a*xa1*(Pa/(K6a+Pa)) - kdpa*Pa
d[Pb]/dt = Vm5b*(Mb+ab)*(Pbi/(K5b+Pbi)) - V6b*xb1*(Pb/(K6b+Pb)) - kdpb*Pb
```

---

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| eps | 150 | Cell cycle time scaling |
| k_Cd_tx_basal | 0.3 | Basal CycD1 transcription |
| k_Cd_tx_Gli_max | 4.0 | Gli-driven CycD1 Vmax |
| k_Cd_tx_MYCN | 15.0 | MYCN-driven CycD1 Vmax |
| K_MYCN_Cd | 1.5 | MYCN half-activation for CycD1 |
| n_MYCN_Cd | 3 | Hill coefficient (cooperativity) |
| K_EZH2_repression | 0.5 | EZH2 half-repression for CycD1 |
| k_MYCN_synth_basal | 0.3 | MYCN autonomous synthesis |
| k_MYCN_synth_Gli | 0.102 | MYCN Gli-dependent synthesis |
| MYCN_amplification | 1.0 / 2.8 | WT / MB amplification |
| Ptch1_copy_number | 1.0 / 0.0 | WT / MB (Ptch1-/-) |
| vs1p27 | 4.0 | p27 basal synthesis (G0 barrier) |
| k_EZH2_mRNA_synth_basal | 0.06 | EZH2 basal transcription |
| k_EZH2_mRNA_synth_E2F | 2.0 | E2F-driven EZH2 transcription |
| K_pRBpp_EZH2 | 0.03 | pRBpp activation of EZH2 |
| k_EZH2_deg | 0.02 | EZH2 protein degradation (slow) |

## Drug Targets

| Drug | Target | Mechanism |
|------|--------|-----------|
| GDC0449 | Smo_active | Blocks Smo activation (SMO inhibitor) |
| EZH2i | EZH2 | Blocks EZH2 catalytic activity (PRC2 inhibitor) |

## Key Feedback Loops

1. **Negative (Ptch1)**: Gli → Ptch1 up → Smo down (homeostatic)
2. **Positive (Gli1)**: Gli_act → Gli1 up → amplifies Gli signal
3. **Negative (EZH2)**: Cycling → E2F → EZH2 up → CycD1 down → slows cycle
4. **Bistable (Rb-E2F)**: CycD/CDK4 → pRB phosph → E2F release → CycE → more pRB phosph
5. **Inhibitory (p27)**: p27 inhibits CDKs (G0 lock); Me → p27 phosph → Skp2 degradation (S-phase entry)
6. **Mitotic exit (APC/C)**: Mb → Cdc20 → CycB degradation; low CDK → Cdh1 → reset
