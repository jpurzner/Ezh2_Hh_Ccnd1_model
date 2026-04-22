"""
Model v42 - MYCN Module Addition to EZH2-CyclinD1-HH Model

Changes from v41:
1. Added MYCN module: MYCN drives CycD1 independently of Gli/SHH
2. MYCN synthesis: largely autonomous + small Gli-dependent component (~22% in WT)
3. MYCN_amplification parameter: 1.0 for WT, ~2.8 for MB (MYCN amplified in SHH-MB)
4. CycD1 transcription restructured: (basal + Gli*Gli_rep + MYCN_Hill) * EZH2_repression
5. MYCN->CycD1 uses Hill function (n_MYCN_Cd=3): threshold effect where low MYCN (WT)
   has negligible CycD1 contribution but high MYCN (amplified MB) drives CycD1 strongly
6. GDC0449 no longer fully ablates CycD1 in MB (MYCN maintains ~40%)

Key experimental constraints (from RNA-seq data):
- MYCN 2.8x higher in MB vs WT (amplified)
- MYCN only 14% decrease with GDC in MB (GDC-resistant)
- MYCN only 22% decrease with GDC in WT
- CycD1: WT+GDC = 86% reduction, MB+GDC = only 60% reduction
- Gli1: WT+GDC = 99% reduction, MB+GDC = 77% reduction
- EZH2: GDC reduces 25-47% (via cell cycle arrest)
- MYCN drives CycD1 independently of SHH/Gli, explaining GDC resistance in MB
"""


def build_model_v42():
    """Build v42 model with MYCN module."""

    antimony_model = r"""
// Model v42 - MYCN-EZH2-CyclinD1-HH Cell Cycle Model
// ==================================================================

model ezh2_cyclind1_v42

// ===================================================================
// MODULE 1: HEDGEHOG SIGNALING -> CycD1_mRNA
// ===================================================================

species $SHH = 0.5;
species SHH_Ptch = 0.0;
species Ptch1_free = 0.5;
species Ptch1_mRNA = 0.5;
species Smo_active = 0.1;
species Gli_rep = 0.5;
species Gli_act = 0.1;
species Gli1_mRNA = 0.1;
species Gli1 = 0.1;

// OUTPUT: CycD1 mRNA level
species Cd_mRNA = 0.5;

// HH pathway parameters
k_Ptch1_tx = 0.5;
k_Ptch1_mRNA_deg = 0.8;
k_Ptch1_translation = 1.0;
k_Ptch1_deg = 0.5;
k_SHH_Ptch_bind = 5.0;
k_SHH_Ptch_release = 0.1;
k_SHH_Ptch_deg = 0.8;
k_Smo_act = 1.5;
k_Smo_inact = 1.2;
K_Ptch_Smo = 0.4;
k_Gli_rep_to_act = 3.0;
k_Gli_act_to_rep = 1.5;
K_Smo_Gli_switch = 0.6;
Vmax_Gli1_tx = 1.2;
K_Gli_act_Gli1 = 0.3;
n_Gli_act = 2;
K_Gli_rep_Gli1 = 0.4;
n_Gli_rep = 2;
k_Gli1_mRNA_deg = 0.8;
k_Gli1_translation = 1.2;
k_Gli1_deg = 0.8;
k_Cd_tx_basal = 0.3;
k_Cd_tx_Gli_max = 4.0;
K_Gli_act_CycD = 0.4;
K_Gli_rep_CycD = 0.3;
k_Cd_mRNA_deg = 0.8;

// Ptch1 gene copy number (1.0=WT, 0.5=het, 0.0=hom loss/MB)
Ptch1_copy_number = 1.0;

// GDC0449 (Vismodegib) SMO inhibitor (0=no drug, 1=full inhibition)
species $GDC0449 = 0.0;

// HH pathway reactions
Ptch1_transcription: -> Ptch1_mRNA; k_Ptch1_tx * Ptch1_copy_number;
Ptch1_mRNA_degradation: Ptch1_mRNA -> ; k_Ptch1_mRNA_deg * Ptch1_mRNA;
Ptch1_translation: Ptch1_mRNA -> Ptch1_mRNA + Ptch1_free; k_Ptch1_translation * Ptch1_mRNA;
Ptch1_degradation: Ptch1_free -> ; k_Ptch1_deg * Ptch1_free;
SHH_Ptch_binding: SHH + Ptch1_free -> SHH_Ptch; k_SHH_Ptch_bind * SHH * Ptch1_free;
SHH_Ptch_release: SHH_Ptch -> SHH + Ptch1_free; k_SHH_Ptch_release * SHH_Ptch;
SHH_Ptch_degradation: SHH_Ptch -> ; k_SHH_Ptch_deg * SHH_Ptch;
Smo_activation: -> Smo_active; k_Smo_act / (1 + Ptch1_free / K_Ptch_Smo) * (1 - GDC0449);
Smo_inactivation: Smo_active -> ; k_Smo_inact * Smo_active;
Gli_rep_to_act: Gli_rep -> Gli_act; k_Gli_rep_to_act * Smo_active^2 / (K_Smo_Gli_switch^2 + Smo_active^2) * Gli_rep;
Gli_act_to_rep: Gli_act -> Gli_rep; k_Gli_act_to_rep * (1 - Smo_active^2 / (K_Smo_Gli_switch^2 + Smo_active^2)) * Gli_act;
Gli1_transcription: -> Gli1_mRNA; Vmax_Gli1_tx * Gli_act^n_Gli_act / (K_Gli_act_Gli1^n_Gli_act + Gli_act^n_Gli_act) * (1 - Gli_rep^n_Gli_rep / (K_Gli_rep_Gli1^n_Gli_rep + Gli_rep^n_Gli_rep));
Gli1_mRNA_degradation: Gli1_mRNA -> ; k_Gli1_mRNA_deg * Gli1_mRNA;
Gli1_translation: Gli1_mRNA -> Gli1_mRNA + Gli1; k_Gli1_translation * Gli1_mRNA;
Gli1_degradation: Gli1 -> ; k_Gli1_deg * Gli1;

// CycD1 transcription: basal + Gli-driven + MYCN-driven, all subject to EZH2 repression
// Gli_rep repression only applies to Gli-driven term (it acts on Gli binding sites)
// EZH2 repression applies to all (H3K27me3 at CycD1 promoter)
CycD1_transcription: -> Cd_mRNA; (k_Cd_tx_basal + k_Cd_tx_Gli_max * (Gli_act + Gli1)^n_Gli_act / (K_Gli_act_CycD^n_Gli_act + (Gli_act + Gli1)^n_Gli_act) * (K_Gli_rep_CycD^n_Gli_rep / (K_Gli_rep_CycD^n_Gli_rep + Gli_rep^n_Gli_rep)) + k_Cd_tx_MYCN * MYCN^n_MYCN_Cd / (K_MYCN_Cd^n_MYCN_Cd + MYCN^n_MYCN_Cd)) * (K_EZH2_repression / (K_EZH2_repression + EZH2 * (1 - EZH2i)));
CycD1_mRNA_degradation: Cd_mRNA -> ; k_Cd_mRNA_deg * Cd_mRNA;

// ===================================================================
// MODULE 1.5: EZH2 REGULATION (Rb-E2F dependent, via pRBpp)
// ===================================================================

species EZH2_mRNA = 0.03;
species EZH2 = 0.5;

k_EZH2_mRNA_synth_basal = 0.06;
k_EZH2_mRNA_synth_E2F = 2.0;
K_E2F_EZH2 = 5.0;
K_pRBpp_EZH2 = 0.03;
k_EZH2_mRNA_deg = 1.5;
k_EZH2_translation = 0.12;
k_EZH2_deg = 0.02;

K_EZH2_repression = 0.5;

species $EZH2i = 0.0;

EZH2_mRNA_synthesis: -> EZH2_mRNA; k_EZH2_mRNA_synth_basal + k_EZH2_mRNA_synth_E2F * E2F / (K_E2F_EZH2 + E2F) * pRBpp / (K_pRBpp_EZH2 + pRBpp);
EZH2_mRNA_degradation: EZH2_mRNA ->; k_EZH2_mRNA_deg * EZH2_mRNA;
EZH2_translation: EZH2_mRNA -> EZH2_mRNA + EZH2; k_EZH2_translation * EZH2_mRNA;
EZH2_degradation: EZH2 ->; k_EZH2_deg * EZH2;

// ===================================================================
// MODULE 1.75: MYCN (GDC-resistant CycD1 driver)
// ===================================================================
// MYCN is amplified in SHH-MB, drives CycD1 independently of Gli.
// Data: MYCN 2.8x higher in MB vs WT, only 14% decrease with GDC in MB.
// MYCN amplification applies to basal (autonomous) component only,
// so MB MYCN is more GDC-resistant than WT MYCN (matches data).
// The Gli-dependent component is a fixed amount (not amplified).

species MYCN = 0.5;

// MYCN parameters
k_MYCN_synth_basal = 0.3;     // Autonomous synthesis (dominant, GDC-resistant)
k_MYCN_synth_Gli = 0.102;    // Small Gli-dependent component (~34% of basal)
K_Gli_MYCN = 0.5;             // Half-activation for Gli on MYCN
k_MYCN_deg = 1.0;             // Degradation rate

// MYCN amplification: 1.0 for WT, ~2.8 for MB
// Amplification applies only to basal (autonomous) component
// This makes MB MYCN more GDC-resistant than WT MYCN (data: 86% vs 78%)
MYCN_amplification = 1.0;

// MYCN-driven CycD1 transcription parameters
k_Cd_tx_MYCN = 15.0;          // MYCN drives CycD1 (high Vmax, offset by high K)
K_MYCN_Cd = 1.5;              // Half-activation for MYCN on CycD1
n_MYCN_Cd = 3;                // Hill coefficient for MYCN on CycD1 (cooperativity)

// MYCN reactions
// Basal*amplification is the autonomous (amplified in MB) component
// Gli-dependent is a small fixed component (not amplified)
MYCN_synthesis: -> MYCN; k_MYCN_synth_basal * MYCN_amplification + k_MYCN_synth_Gli * (Gli_act + Gli1) / (K_Gli_MYCN + Gli_act + Gli1);
MYCN_degradation: MYCN -> ; k_MYCN_deg * MYCN;

// ===================================================================
// MODULE 2: GERARD 2009 CELL CYCLE (Simplified)
// ===================================================================

eps = 150;

species pRB = 1.0;
species pRBp = 0.25;
species pRBpp = 0.1;
species pRBc1 = 0.1;
species pRBc2 = 0.05;

species E2F = 0.1;
species E2Fp = 0.001;

species Cd = 0.01;
species Mdi = 0.01;
species Md = 0.01;
species Mdp27 = 0.01;

species Ce = 0.01;
species Mei = 0.01;
species Me = 0.01;
species Mep27 = 0.01;

species Ca = 0.01;
species Mai = 0.01;
species Ma = 0.01;
species Map27 = 0.01;

species Cb = 0.01;
species Mbi = 0.01;
species Mb = 0.01;
species Mbp27 = 0.01;

species p27 = 0.25;
species p27p = 0.01;

species Skp2 = 0.01;

species Cdh1i = 0.01;
species Cdh1a = 0.01;

species Cdc20i = 0.01;
species Cdc20a = 0.01;

species Pei = 0.01;
species Pe = 0.01;
species Pai = 0.01;
species Pa = 0.01;
species Pbi = 0.01;
species Pb = 0.1;

// Gerard parameters
vsprb = 0.8;
kdprb = 0.01;
kdprbp = 0.06;
kdprbpp = 0.04;
kpc1 = 0.05;
kpc2 = 0.5;
kpc3 = 0.025;
kpc4 = 0.5;
V1 = 2.2;
V2 = 2.0;
V3 = 1.0;
V4 = 2.0;
K1 = 0.1;
K2 = 0.1;
K3 = 0.1;
K4 = 0.1;

vse2f = 0.15;
kde2f = 0.002;
kde2fp = 1.1;
V1e2f = 4.0;
V2e2f = 0.75;
K1e2f = 5.0;
K2e2f = 5.0;

k_Cd_translation = 0.25;
Cdk4_tot = 1.5;
kcom1 = 0.175;
kdecom1 = 0.1;
Vm1d = 1.0;
Vm2d = 0.2;
K1d = 0.1;
K2d = 0.1;
kc1 = 0.15;
kc2 = 0.05;
Vdd = 5.0;
Kdd = 0.1;
kddd = 0.005;
Ki7 = 0.1;
Ki8 = 2.0;

Cdk2_tot = 2.0;
kce = 0.25;
kcom2 = 0.2;
kdecom2 = 0.1;
Vm1e = 2.0;
Vm2e = 1.4;
K1e = 0.1;
K2e = 0.1;
kc3 = 0.2;
kc4 = 0.1;
Vde = 3.0;
Kde = 0.1;
kdde = 0.005;
Ki9 = 0.1;
Ki10 = 2.0;
ae = 0.25;
ib1 = 0.5;

kca = 0.0375;
kcom3 = 0.2;
kdecom3 = 0.1;
Vm1a = 2.0;
Vm2a = 1.85;
K1a = 0.1;
K2a = 0.1;
kc5 = 0.15;
kc6 = 0.125;
Vda = 2.5;
Kda = 1.1;
kdda = 0.005;
Ki11 = 0.1;
Ki12 = 2.0;
aa = 0.2;
ib2 = 0.5;

Cdk1_tot = 0.5;
vcb = 0.05;
kcom4 = 0.2;
kdecom4 = 0.1;
Vm1b = 1.9;
Vm2b = 1.75;
K1b = 0.1;
K2b = 0.1;
kc7 = 0.175;
kc8 = 0.125;
Vdb = 2.5;
Kdb = 0.1;
kddb = 0.005;
Kdbcdh1 = 0.1;
Kdbcdc20 = 0.2;
ab = 0.11;
ib3 = 0.5;

vs1p27 = 4.0;
vs2p27 = 0.1;
kddp27 = 0.06;
kddp27p = 0.01;
V1p27 = 100.0;
V2p27 = 0.1;
K1p27 = 0.5;
K2p27 = 0.5;
Vdp27p = 5.0;
Kdp27p = 0.1;
Kdp27skp2 = 0.1;
Ki13 = 0.1;
Ki14 = 2.0;

vsskp2 = 0.15;
Vdskp2 = 1.1;
Kdskp2 = 0.5;
Kdceskp2 = 2.0;
kddskp2 = 0.005;

vscdh1a = 0.11;
V1cdh1 = 1.25;
V2cdh1 = 8.0;
K1cdh1 = 0.01;
K2cdh1 = 0.01;
Kcdh1 = 0.4;
kdcdh1i = 0.2;
kdcdh1a = 0.1;

vscdc20i = 0.1;
Vm3b = 5.0;
Vm4b = 2.5;
K3b = 0.1;
K4b = 0.1;
Kacdc20 = 2.0;
kdcdc20i = 0.05;
kdcdc20a = 0.05;

vspei = 0.13;
Vm5e = 5.0;
V6e = 0.8;
K5e = 0.1;
K6e = 0.1;
kdpei = 0.15;
kdpe = 0.075;
xe1 = 1.0;
xe2 = 0.0;

vspai = 0.105;
Vm5a = 4.0;
V6a = 1.0;
K5a = 0.1;
K6a = 0.1;
kdpai = 0.15;
kdpa = 0.075;
xa1 = 1.0;
xa2 = 0.0;

vspbi = 0.07;
Vm5b = 4.5;
V6b = 0.9;
K5b = 0.1;
K6b = 0.1;
kdpbi = 0.15;
kdpb = 0.075;
xb1 = 1.0;
xb2 = 0.0;

// ===================================================================
// CELL CYCLE REACTIONS
// ===================================================================

pRB_synthesis: -> pRB; vsprb * eps;
pRB_decay: pRB -> ; kdprb * pRB * eps;
pRBp_decay: pRBp -> ; kdprbp * pRBp * eps;
pRBpp_decay: pRBpp -> ; kdprbpp * pRBpp * eps;

pRB_E2F_complex_formation: pRB + E2F -> pRBc1; kpc1 * pRB * E2F * eps;
pRB_E2F_complex_dissociation: pRBc1 -> pRB + E2F; kpc2 * pRBc1 * eps;
pRBp_E2F_complex_formation: pRBp + E2F -> pRBc2; kpc3 * pRBp * E2F * eps;
pRBp_E2F_complex_dissociation: pRBc2 -> pRBp + E2F; kpc4 * pRBc2 * eps;

pRB_phosphorylation: pRB -> pRBp; V1 * (pRB / (K1 + pRB)) * (Md + Mdp27) * eps;
pRBp_dephosphorylation: pRBp -> pRB; V2 * (pRBp / (K2 + pRBp)) * eps;
pRBp_phosphorylation: pRBp -> pRBpp; V3 * (pRBp / (K3 + pRBp)) * Me * eps;
pRBpp_dephosphorylation: pRBpp -> pRBp; V4 * (pRBpp / (K4 + pRBpp)) * eps;

E2F_synthesis: -> E2F; vse2f * eps;
E2F_decay: E2F -> ; kde2f * E2F * eps;
E2F_phosphorylation: E2F -> E2Fp; V1e2f * (E2F / (K1e2f + E2F)) * Ma * eps;
E2Fp_dephosphorylation: E2Fp -> E2F; V2e2f * (E2Fp / (K2e2f + E2Fp)) * eps;
E2Fp_decay: E2Fp -> ; kde2fp * E2Fp * eps;

Cd_translation: Cd_mRNA -> Cd_mRNA + Cd; k_Cd_translation * Cd_mRNA * eps;
Cd_degradation: Cd -> ; Vdd * (Cd / (Kdd + Cd)) * eps;
Cd_decay: Cd -> ; kddd * Cd * eps;

Cd_CDK4_formation: Cd -> Mdi; kcom1 * Cd * (Cdk4_tot - (Mdi + Md + Mdp27)) * eps;
Cd_CDK4_dissociation: Mdi -> Cd; kdecom1 * Mdi * eps;
Md_activation: Mdi -> Md; Vm1d * (Mdi / (K1d + Mdi)) * eps;
Md_inactivation: Md -> Mdi; Vm2d * (Md / (K2d + Md)) * eps;
Md_p27_binding: Md + p27 -> Mdp27; kc1 * Md * p27 * eps;
Md_p27_release: Mdp27 -> Md + p27; kc2 * Mdp27 * eps;

Ce_synthesis: -> Ce; kce * E2F * (Ki9 / (Ki9 + pRB)) * (Ki10 / (Ki10 + pRBp)) * eps;
Ce_degradation: Ce -> ; Vde * (Skp2 / (Kdceskp2 + Skp2)) * (Ce / (Kde + Ce)) * eps;
Ce_decay: Ce -> ; kdde * Ce * eps;

Ce_CDK2_formation: Ce -> Mei; kcom2 * Ce * (Cdk2_tot - (Mei + Me + Mep27 + Mai + Ma + Map27)) * eps;
Ce_CDK2_dissociation: Mei -> Ce; kdecom2 * Mei * eps;
Me_activation: Mei -> Me; Vm1e * (Mei / (K1e + Mei)) * Pe * eps;
Me_inactivation: Me -> Mei; Vm2e * (ib1) * (Me / (K2e + Me)) * eps;
Me_p27_binding: Me + p27 -> Mep27; kc3 * Me * p27 * eps;
Me_p27_release: Mep27 -> Me + p27; kc4 * Mep27 * eps;

Ca_synthesis: -> Ca; kca * E2F * (Ki11 / (Ki11 + pRB)) * (Ki12 / (Ki12 + pRBp)) * eps;
Ca_degradation: Ca -> ; Vda * (Ca / (Kda + Ca)) * (Cdc20a / (Kacdc20 + Cdc20a)) * eps;
Ca_decay: Ca -> ; kdda * Ca * eps;

Ca_CDK2_formation: Ca -> Mai; kcom3 * Ca * (Cdk2_tot - (Mei + Me + Mep27 + Mai + Ma + Map27)) * eps;
Ca_CDK2_dissociation: Mai -> Ca; kdecom3 * Mai * eps;
Ma_activation: Mai -> Ma; Vm1a * (Mai / (K1a + Mai)) * Pa * eps;
Ma_inactivation: Ma -> Mai; Vm2a * (ib2) * (Ma / (K2a + Ma)) * eps;
Ma_p27_binding: Ma + p27 -> Map27; kc5 * Ma * p27 * eps;
Ma_p27_release: Map27 -> Ma + p27; kc6 * Map27 * eps;

Cb_synthesis: -> Cb; vcb * eps;
Cb_degradation_Cdc20: Cb -> ; Vdb * (Cb / (Kdb + Cb)) * (Cdc20a / (Kdbcdc20 + Cdc20a)) * eps;
Cb_degradation_Cdh1: Cb -> ; Vdb * (Cb / (Kdb + Cb)) * (Cdh1a / (Kdbcdh1 + Cdh1a)) * eps;
Cb_decay: Cb -> ; kddb * Cb * eps;

Cb_CDK1_formation: Cb -> Mbi; kcom4 * Cb * (Cdk1_tot - (Mbi + Mb + Mbp27)) * eps;
Cb_CDK1_dissociation: Mbi -> Cb; kdecom4 * Mbi * eps;
Mb_activation: Mbi -> Mb; Vm1b * (Mbi / (K1b + Mbi)) * Pb * eps;
Mb_inactivation: Mb -> Mbi; Vm2b * (ib3) * (Mb / (K2b + Mb)) * eps;
Mb_p27_binding: Mb + p27 -> Mbp27; kc7 * Mb * p27 * eps;
Mb_p27_release: Mbp27 -> Mb + p27; kc8 * Mbp27 * eps;

p27_synthesis: -> p27; vs1p27 * eps;
p27_synthesis_E2F: -> p27; vs2p27 * E2F * (Ki13 / (Ki13 + pRB)) * (Ki14 / (Ki14 + pRBp)) * eps;
p27_phosphorylation: p27 -> p27p; V1p27 * (p27 / (K1p27 + p27)) * Me * eps;
p27p_dephosphorylation: p27p -> p27; V2p27 * (p27p / (K2p27 + p27p)) * eps;
p27_decay: p27 -> ; kddp27 * p27 * eps;
p27p_degradation: p27p -> ; Vdp27p * (Skp2 / (Kdp27skp2 + Skp2)) * (p27p / (Kdp27p + p27p)) * eps;
p27p_decay: p27p -> ; kddp27p * p27p * eps;

Skp2_synthesis: -> Skp2; vsskp2 * eps;
Skp2_degradation: Skp2 -> ; Vdskp2 * (Skp2 / (Kdskp2 + Skp2)) * (Cdh1a / (Kcdh1 + Cdh1a)) * eps;
Skp2_decay: Skp2 -> ; kddskp2 * Skp2 * eps;

Cdh1_synthesis: -> Cdh1a; vscdh1a * eps;
Cdh1_activation: Cdh1i -> Cdh1a; V1cdh1 * (Cdh1i / (K1cdh1 + Cdh1i)) * eps;
Cdh1_inactivation: Cdh1a -> Cdh1i; V2cdh1 * (Cdh1a / (K2cdh1 + Cdh1a)) * (Ma + Mb) * eps;
Cdh1i_decay: Cdh1i -> ; kdcdh1i * Cdh1i * eps;
Cdh1a_decay: Cdh1a -> ; kdcdh1a * Cdh1a * eps;

Cdc20i_synthesis: -> Cdc20i; vscdc20i * eps;
Cdc20_activation: Cdc20i -> Cdc20a; Vm3b * (Cdc20i / (K3b + Cdc20i)) * Mb * eps;
Cdc20_inactivation: Cdc20a -> Cdc20i; Vm4b * (Cdc20a / (K4b + Cdc20a)) * eps;
Cdc20i_decay: Cdc20i -> ; kdcdc20i * Cdc20i * eps;
Cdc20a_decay: Cdc20a -> ; kdcdc20a * Cdc20a * eps;

Cdc25E_synthesis: -> Pei; vspei * eps;
Cdc25E_activation: Pei -> Pe; Vm5e * (Me + ae) * (Pei / (K5e + Pei)) * eps;
Cdc25E_inactivation: Pe -> Pei; V6e * (xe1) * (Pe / (K6e + Pe)) * eps;
Pei_decay: Pei -> ; kdpei * Pei * eps;
Pe_decay: Pe -> ; kdpe * Pe * eps;

Cdc25A_synthesis: -> Pai; vspai * eps;
Cdc25A_activation: Pai -> Pa; Vm5a * (Ma + aa) * (Pai / (K5a + Pai)) * eps;
Cdc25A_inactivation: Pa -> Pai; V6a * (xa1) * (Pa / (K6a + Pa)) * eps;
Pai_decay: Pai -> ; kdpai * Pai * eps;
Pa_decay: Pa -> ; kdpa * Pa * eps;

Cdc25B_synthesis: -> Pbi; vspbi * eps;
Cdc25B_activation: Pbi -> Pb; Vm5b * (Mb + ab) * (Pbi / (K5b + Pbi)) * eps;
Cdc25B_inactivation: Pb -> Pbi; V6b * (xb1) * (Pb / (K6b + Pb)) * eps;
Pbi_decay: Pbi -> ; kdpbi * Pbi * eps;
Pb_decay: Pb -> ; kdpb * Pb * eps;

// ===================================================================
// INITIAL CONDITIONS
// ===================================================================

SHH = 0.5;
Ptch1_free = 0.5;
Smo_active = 0.1;
Gli_rep = 0.5;
Gli_act = 0.1;
Gli1 = 0.1;
Cd_mRNA = 0.5;

MYCN = 0.5;

pRB = 1.0;
pRBp = 0.25;
pRBpp = 0.1;
E2F = 0.1;
E2Fp = 0.001;
pRBc1 = 0.1;
pRBc2 = 0.05;
Cd = 0.01;
Mdi = 0.01;
Md = 0.01;
Mdp27 = 0.01;
Ce = 0.01;
Mei = 0.01;
Me = 0.01;
Mep27 = 0.01;
Pei = 0.01;
Pe = 0.01;
Ca = 0.01;
Mai = 0.01;
Ma = 0.01;
Map27 = 0.01;
p27 = 0.25;
p27p = 0.01;
Cdh1i = 0.01;
Cdh1a = 0.01;
Pai = 0.01;
Pa = 0.01;
Cb = 0.01;
Mbi = 0.01;
Mb = 0.01;
Mbp27 = 0.01;
Cdc20i = 0.01;
Cdc20a = 0.01;
Pbi = 0.01;
Pb = 0.1;
Skp2 = 0.01;

end
"""
    return antimony_model


if __name__ == "__main__":
    print("v42 Model: MYCN-EZH2-CyclinD1-HH Cell Cycle")
    print("=" * 70)
    print("\nChanges from v41:")
    print("  1. Added MYCN module (GDC-resistant CycD1 driver)")
    print("  2. MYCN_amplification parameter (1.0 WT, ~2.8 MB)")
    print("  3. CycD1 = basal + Gli*Gli_rep + MYCN, all x EZH2 repression")
    print("  4. MYCN basal is amplified in MB, Gli component is not")

    model_str = build_model_v42()
    print(f"\nModel built: {len(model_str)} characters")
