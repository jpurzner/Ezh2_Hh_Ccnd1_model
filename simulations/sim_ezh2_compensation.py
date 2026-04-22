"""
Does higher EZH2 in MB compensate for MYCN-driven CycD1?

Tests:
1. Current model: GNP vs MB EZH2 levels and CycD1
2. MB with boosted EZH2 (simulating MYCN→EZH2 drive): how much CycD1 is suppressed?
3. MB with no EZH2 feedback: what would CycD1 be without any compensation?

This quantifies how much the EZH2→CycD1 repression buffers MYCN's proliferative drive.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import tellurium as te
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.signal import find_peaks
from src.build_model_v42_mycn import build_model_v42

model_str = build_model_v42()
MYCN_AMP_MB = 2.8


def run_sim(context='wt_shh', ezh2_boost=1.0, no_ezh2_fb=False,
            t_end=96, n_pts=9600):
    rr = te.loada(model_str)
    if context == 'wt_shh':
        rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 1.0
        rr['MYCN_amplification'] = 1.0
    elif context == 'mb':
        rr['SHH'] = 0.0; rr['Ptch1_copy_number'] = 0.0
        rr['MYCN_amplification'] = MYCN_AMP_MB
        rr['Ptch1_mRNA'] = 0.0; rr['Ptch1_free'] = 0.0
        rr['SHH_Ptch'] = 0.0; rr['Smo_active'] = 1.0
    rr['GDC0449'] = 0.0; rr['EZH2i'] = 0.0

    # Boost EZH2 synthesis (simulating MYCN→EZH2 transcription)
    rr['k_EZH2_mRNA_synth_basal'] *= ezh2_boost

    if no_ezh2_fb:
        rr['K_EZH2_repression'] = 1e6

    return rr.simulate(0, t_end, n_pts)


def analyze(result, t_start=10):
    t = result['time']; cb = result['[Cb]']
    mask = t >= t_start
    peaks, _ = find_peaks(cb[mask], prominence=0.05, distance=10)
    peak_times = t[mask][peaks]
    periods = np.diff(peak_times) if len(peak_times) > 1 else np.array([])
    idx = -1500
    return {
        'n_div': len(peaks),
        'period': np.mean(periods) if len(periods) > 0 else None,
        'cd_mean': np.mean(result['[Cd_mRNA]'][idx:]),
        'ezh2_mean': np.mean(result['[EZH2]'][idx:]),
        'ezh2_mRNA_mean': np.mean(result['[EZH2_mRNA]'][idx:]),
        'e2f_mean': np.mean(result['[E2F]'][idx:]),
        'p27_mean': np.mean(result['[p27]'][idx:]),
        'mycn_mean': np.mean(result['[MYCN]'][idx:]),
        'prbpp_mean': np.mean(result['[pRBpp]'][idx:]),
    }


# =====================================================================
# PART 1: Current model — GNP vs MB EZH2 and CycD1
# =====================================================================

print("=" * 90)
print("PART 1: Current model — GNP vs MB")
print("=" * 90)

gnp = analyze(run_sim('wt_shh'))
mb = analyze(run_sim('mb'))
mb_no_fb = analyze(run_sim('mb', no_ezh2_fb=True))
gnp_no_fb = analyze(run_sim('wt_shh', no_ezh2_fb=True))

print(f"\n  {'':25s} {'GNP+SHH':>10s} {'MB':>10s} {'MB/GNP':>10s}")
print(f"  {'-'*60}")
print(f"  {'EZH2 protein (mean)':25s} {gnp['ezh2_mean']:10.3f} {mb['ezh2_mean']:10.3f} {mb['ezh2_mean']/gnp['ezh2_mean']:10.2f}x")
print(f"  {'EZH2 mRNA (mean)':25s} {gnp['ezh2_mRNA_mean']:10.4f} {mb['ezh2_mRNA_mean']:10.4f} {mb['ezh2_mRNA_mean']/gnp['ezh2_mRNA_mean']:10.2f}x")
print(f"  {'CycD1 mRNA (mean)':25s} {gnp['cd_mean']:10.3f} {mb['cd_mean']:10.3f} {mb['cd_mean']/gnp['cd_mean']:10.2f}x")
print(f"  {'MYCN (mean)':25s} {gnp['mycn_mean']:10.3f} {mb['mycn_mean']:10.3f} {mb['mycn_mean']/gnp['mycn_mean']:10.2f}x")
print(f"  {'Period (h)':25s} {gnp['period']:10.1f} {mb['period']:10.1f}")
print(f"  {'E2F (mean)':25s} {gnp['e2f_mean']:10.1f} {mb['e2f_mean']:10.1f} {mb['e2f_mean']/gnp['e2f_mean']:10.2f}x")

# Compute theoretical CycD1 without any EZH2 repression
print(f"\n  CycD1 without EZH2 fb:  GNP={gnp_no_fb['cd_mean']:.2f}, MB={mb_no_fb['cd_mean']:.2f}")
print(f"  Compensation by EZH2:   GNP: {(1 - gnp['cd_mean']/gnp_no_fb['cd_mean'])*100:.0f}% of CycD1 suppressed")
print(f"                          MB:  {(1 - mb['cd_mean']/mb_no_fb['cd_mean'])*100:.0f}% of CycD1 suppressed")

# The EZH2 repression factor = K/(K + EZH2)
K = 0.5
rep_gnp = K / (K + gnp['ezh2_mean'])
rep_mb = K / (K + mb['ezh2_mean'])
print(f"\n  EZH2 repression factor: GNP={rep_gnp:.3f}, MB={rep_mb:.3f}")
print(f"  → MB CycD1 is {(1-rep_mb)/(1-rep_gnp)*100 - 100:+.0f}% more repressed by EZH2 than GNP")

# =====================================================================
# PART 2: Sweep EZH2 boost in MB — how much does extra EZH2 help?
# =====================================================================

print(f"\n{'=' * 90}")
print("PART 2: EZH2 boost sweep in MB (simulating MYCN→EZH2 transcription)")
print("=" * 90)

boosts = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 8.0]
boost_results = []

print(f"\n  {'EZH2 boost':>10s} {'EZH2 prot':>10s} {'EZH2 mRNA':>10s} {'CycD1':>8s} {'Period':>8s} {'Div':>5s} {'Rep factor':>12s}")
print(f"  {'-'*70}")

for b in boosts:
    r = run_sim('mb', ezh2_boost=b)
    m = analyze(r)
    rep = K / (K + m['ezh2_mean'])
    boost_results.append(m)
    p_str = f"{m['period']:.1f}h" if m['period'] else "N/A"
    print(f"  {b:10.1f}x {m['ezh2_mean']:10.3f} {m['ezh2_mRNA_mean']:10.4f} {m['cd_mean']:8.2f} {p_str:>8s} {m['n_div']:5d} {rep:12.3f}")

# Reference: what would CycD1 be with no EZH2 at all?
print(f"\n  No feedback:  CycD1={mb_no_fb['cd_mean']:.2f} (100% of unrepressed level)")
print(f"  Current (1x): CycD1={mb['cd_mean']:.2f} ({mb['cd_mean']/mb_no_fb['cd_mean']*100:.0f}% of unrepressed)")

# Find how much boost needed to bring CycD1 down to GNP level
print(f"\n  GNP CycD1 = {gnp['cd_mean']:.2f}")
for i, (b, m) in enumerate(zip(boosts, boost_results)):
    if m['cd_mean'] <= gnp['cd_mean']:
        print(f"  → Need ~{b:.0f}x EZH2 boost to bring MB CycD1 to GNP level")
        break
else:
    print(f"  → Even {boosts[-1]}x EZH2 boost can't bring MB CycD1 to GNP level")
    # Check the last one
    print(f"    (at {boosts[-1]}x: CycD1={boost_results[-1]['cd_mean']:.2f} vs GNP={gnp['cd_mean']:.2f})")

# =====================================================================
# PART 3: Decompose CycD1 sources in MB
# =====================================================================

print(f"\n{'=' * 90}")
print("PART 3: CycD1 transcription decomposition")
print("=" * 90)

# Run MB to get steady-state species values
r_mb = run_sim('mb')
idx = -1500
gli_act = np.mean(r_mb['[Gli_act]'][idx:])
gli1 = np.mean(r_mb['[Gli1]'][idx:])
gli_rep = np.mean(r_mb['[Gli_rep]'][idx:])
mycn = np.mean(r_mb['[MYCN]'][idx:])
ezh2 = np.mean(r_mb['[EZH2]'][idx:])

# CycD1 transcription components
k_basal = 0.3
k_gli = 4.0
K_gli_act = 0.4; n_gli = 2
K_gli_rep = 0.3; n_gli_rep = 2
k_mycn = 15.0
K_mycn = 1.5; n_mycn = 3
K_ezh2 = 0.5

gli_total = gli_act + gli1
gli_drive = k_gli * gli_total**n_gli / (K_gli_act**n_gli + gli_total**n_gli) * \
            K_gli_rep**n_gli_rep / (K_gli_rep**n_gli_rep + gli_rep**n_gli_rep)
mycn_drive = k_mycn * mycn**n_mycn / (K_mycn**n_mycn + mycn**n_mycn)
total_unrepressed = k_basal + gli_drive + mycn_drive
ezh2_factor = K_ezh2 / (K_ezh2 + ezh2)
total_repressed = total_unrepressed * ezh2_factor

print(f"\n  MB steady-state species:")
print(f"    Gli_act={gli_act:.3f}, Gli1={gli1:.3f}, Gli_rep={gli_rep:.3f}")
print(f"    MYCN={mycn:.3f}, EZH2={ezh2:.3f}")

print(f"\n  CycD1 transcription breakdown (MB):")
print(f"    Basal:          {k_basal:6.2f}  ({k_basal/total_unrepressed*100:4.1f}%)")
print(f"    Gli-driven:     {gli_drive:6.2f}  ({gli_drive/total_unrepressed*100:4.1f}%)")
print(f"    MYCN-driven:    {mycn_drive:6.2f}  ({mycn_drive/total_unrepressed*100:4.1f}%)")
print(f"    ─────────────────────────")
print(f"    Total (unrepr): {total_unrepressed:6.2f}")
print(f"    EZH2 factor:    {ezh2_factor:6.3f}  (→ {(1-ezh2_factor)*100:.0f}% repressed)")
print(f"    Total (repr):   {total_repressed:6.2f}")

# Same for GNP
r_gnp = run_sim('wt_shh')
gli_act_g = np.mean(r_gnp['[Gli_act]'][idx:])
gli1_g = np.mean(r_gnp['[Gli1]'][idx:])
gli_rep_g = np.mean(r_gnp['[Gli_rep]'][idx:])
mycn_g = np.mean(r_gnp['[MYCN]'][idx:])
ezh2_g = np.mean(r_gnp['[EZH2]'][idx:])

gli_total_g = gli_act_g + gli1_g
gli_drive_g = k_gli * gli_total_g**n_gli / (K_gli_act**n_gli + gli_total_g**n_gli) * \
              K_gli_rep**n_gli_rep / (K_gli_rep**n_gli_rep + gli_rep_g**n_gli_rep)
mycn_drive_g = k_mycn * mycn_g**n_mycn / (K_mycn**n_mycn + mycn_g**n_mycn)
total_unrep_g = k_basal + gli_drive_g + mycn_drive_g
ezh2_factor_g = K_ezh2 / (K_ezh2 + ezh2_g)
total_rep_g = total_unrep_g * ezh2_factor_g

print(f"\n  CycD1 transcription breakdown (GNP + SHH):")
print(f"    Basal:          {k_basal:6.2f}  ({k_basal/total_unrep_g*100:4.1f}%)")
print(f"    Gli-driven:     {gli_drive_g:6.2f}  ({gli_drive_g/total_unrep_g*100:4.1f}%)")
print(f"    MYCN-driven:    {mycn_drive_g:6.2f}  ({mycn_drive_g/total_unrep_g*100:4.1f}%)")
print(f"    ─────────────────────────")
print(f"    Total (unrepr): {total_unrep_g:6.2f}")
print(f"    EZH2 factor:    {ezh2_factor_g:6.3f}  (→ {(1-ezh2_factor_g)*100:.0f}% repressed)")
print(f"    Total (repr):   {total_rep_g:6.2f}")

print(f"\n  KEY: MYCN contributes {mycn_drive/total_unrepressed*100:.0f}% of MB CycD1 vs {mycn_drive_g/total_unrep_g*100:.0f}% of GNP CycD1")
print(f"  KEY: EZH2 represses {(1-ezh2_factor)*100:.0f}% of MB CycD1 vs {(1-ezh2_factor_g)*100:.0f}% of GNP CycD1")

# How much of the MYCN-driven excess is compensated by the extra EZH2?
# Counterfactual: MB CycD1 if EZH2 were at GNP level
total_mb_with_gnp_ezh2 = total_unrepressed * ezh2_factor_g
compensation = total_mb_with_gnp_ezh2 - total_repressed
print(f"\n  If MB had GNP-level EZH2: CycD1 would be {total_mb_with_gnp_ezh2:.2f} (vs {total_repressed:.2f})")
print(f"  Extra EZH2 in MB suppresses {compensation:.2f} units of CycD1 ({compensation/total_mb_with_gnp_ezh2*100:.1f}% additional repression)")


# =====================================================================
# FIGURE
# =====================================================================

fig = plt.figure(figsize=(18, 14))
gs = GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.35)

# --- Panel A: CycD1 decomposition bar chart ---
ax = fig.add_subplot(gs[0, 0])
categories = ['GNP\n+SHH', 'MB\n(Ptch1-/-)']
basal_vals = [k_basal * ezh2_factor_g, k_basal * ezh2_factor]
gli_vals = [gli_drive_g * ezh2_factor_g, gli_drive * ezh2_factor]
mycn_vals = [mycn_drive_g * ezh2_factor_g, mycn_drive * ezh2_factor]

x = np.arange(2)
w = 0.5
b1 = ax.bar(x, basal_vals, w, color='#bdc3c7', edgecolor='black', label='Basal')
b2 = ax.bar(x, gli_vals, w, bottom=basal_vals, color='#27ae60', edgecolor='black', label='Gli-driven')
b3 = ax.bar(x, mycn_vals, w, bottom=[b+g for b,g in zip(basal_vals, gli_vals)],
            color='#e67e22', edgecolor='black', label='MYCN-driven')
# Show EZH2-repressed portion as faded
unrep_vals = [total_unrep_g, total_unrepressed]
rep_vals = [total_rep_g, total_repressed]
for i in range(2):
    ax.bar(x[i], unrep_vals[i] - rep_vals[i], w, bottom=rep_vals[i],
           color='#e8daef', edgecolor='#8e44ad', alpha=0.5, hatch='///',
           label='EZH2-repressed' if i == 0 else None)
ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=10)
ax.set_ylabel('CycD1 transcription rate', fontsize=10)
ax.set_title('A. CycD1 Source Decomposition', fontsize=11, fontweight='bold', loc='left')
ax.legend(fontsize=8, loc='upper left')
ax.grid(True, alpha=0.2, axis='y')

# --- Panel B: EZH2 levels comparison ---
ax = fig.add_subplot(gs[0, 1])
ezh2_vals = [gnp['ezh2_mean'], mb['ezh2_mean']]
ezh2_mRNA_vals = [gnp['ezh2_mRNA_mean'], mb['ezh2_mRNA_mean']]
x = np.arange(2)
ax.bar(x - 0.15, ezh2_vals, 0.3, color=['#1b9e77', '#762A83'], edgecolor='black',
       alpha=0.8, label='Protein')
ax2 = ax.twinx()
ax2.bar(x + 0.15, ezh2_mRNA_vals, 0.3, color=['#1b9e77', '#762A83'], edgecolor='black',
        alpha=0.4, hatch='...', label='mRNA')
ax.set_xticks(x)
ax.set_xticklabels(['GNP+SHH', 'MB'], fontsize=10)
ax.set_ylabel('EZH2 Protein', fontsize=10, color='black')
ax2.set_ylabel('EZH2 mRNA', fontsize=10, color='#7f8c8d')
ax.set_title('B. EZH2 Levels: GNP vs MB', fontsize=11, fontweight='bold', loc='left')
# Manual legend
from matplotlib.patches import Patch
ax.legend([Patch(facecolor='gray', edgecolor='black', alpha=0.8),
           Patch(facecolor='gray', edgecolor='black', alpha=0.4, hatch='...')],
          ['Protein', 'mRNA'], fontsize=8, loc='upper left')
ax.grid(True, alpha=0.2, axis='y')

# --- Panel C: EZH2 boost effect on CycD1 ---
ax = fig.add_subplot(gs[0, 2])
ezh2_prot_vals = [m['ezh2_mean'] for m in boost_results]
cd_vals = [m['cd_mean'] for m in boost_results]
ax.plot(ezh2_prot_vals, cd_vals, 'o-', color='#e67e22', linewidth=2, markersize=6)
# Mark current MB and GNP levels
ax.axhline(gnp['cd_mean'], color='#1b9e77', linestyle='--', alpha=0.6, label=f'GNP CycD1={gnp["cd_mean"]:.1f}')
ax.axhline(mb_no_fb['cd_mean'], color='#e74c3c', linestyle=':', alpha=0.6, label=f'No fb CycD1={mb_no_fb["cd_mean"]:.1f}')
ax.axvline(mb['ezh2_mean'], color='#762A83', linestyle='--', alpha=0.5, label=f'Current MB EZH2')
for b, m in zip(boosts, boost_results):
    if b in [1.0, 2.0, 4.0, 8.0]:
        ax.annotate(f'{b:.0f}x', (m['ezh2_mean'], m['cd_mean']),
                    textcoords="offset points", xytext=(8, 5), fontsize=7.5)
ax.set_xlabel('EZH2 Protein Level', fontsize=10)
ax.set_ylabel('CycD1 mRNA', fontsize=10)
ax.set_title('C. EZH2 Boost → CycD1 in MB', fontsize=11, fontweight='bold', loc='left')
ax.legend(fontsize=7.5)
ax.grid(True, alpha=0.2)

# --- Panel D: EZH2 boost effect on period ---
ax = fig.add_subplot(gs[1, 0])
periods = [m['period'] if m['period'] else 0 for m in boost_results]
ax.plot(boosts, periods, 's-', color='#2980b9', linewidth=2, markersize=6)
ax.axhline(gnp['period'], color='#1b9e77', linestyle='--', alpha=0.6, label=f'GNP period={gnp["period"]:.1f}h')
ax.set_xlabel('EZH2 Synthesis Boost (fold)', fontsize=10)
ax.set_ylabel('Cell Cycle Period (h)', fontsize=10)
ax.set_title('D. EZH2 Boost → Period in MB', fontsize=11, fontweight='bold', loc='left')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.2)

# --- Panel E: EZH2 boost effect on divisions ---
ax = fig.add_subplot(gs[1, 1])
divs = [m['n_div'] for m in boost_results]
ax.bar(range(len(boosts)), divs, color='#2980b9', edgecolor='black', alpha=0.7)
ax.set_xticks(range(len(boosts)))
ax.set_xticklabels([f'{b:.1f}x' for b in boosts], fontsize=8)
ax.set_xlabel('EZH2 Synthesis Boost', fontsize=10)
ax.set_ylabel('Divisions (96h)', fontsize=10)
ax.set_title('E. EZH2 Boost → Divisions in MB', fontsize=11, fontweight='bold', loc='left')
for i, v in enumerate(divs):
    ax.text(i, v + 0.3, str(v), ha='center', fontweight='bold', fontsize=9)
ax.grid(True, alpha=0.2, axis='y')

# --- Panel F: Repression factor diagram ---
ax = fig.add_subplot(gs[1, 2])
ezh2_range = np.linspace(0, 3, 200)
rep_factor = K_ezh2 / (K_ezh2 + ezh2_range)
ax.plot(ezh2_range, rep_factor, 'k-', linewidth=2)
ax.fill_between(ezh2_range, rep_factor, 1.0, alpha=0.15, color='#8e44ad', label='Repressed fraction')
# Mark GNP and MB operating points
ax.plot(gnp['ezh2_mean'], K_ezh2/(K_ezh2 + gnp['ezh2_mean']), 'o', color='#1b9e77',
        markersize=10, zorder=5, label=f'GNP (EZH2={gnp["ezh2_mean"]:.2f})')
ax.plot(mb['ezh2_mean'], K_ezh2/(K_ezh2 + mb['ezh2_mean']), 's', color='#762A83',
        markersize=10, zorder=5, label=f'MB (EZH2={mb["ezh2_mean"]:.2f})')
ax.set_xlabel('EZH2 Protein', fontsize=10)
ax.set_ylabel('CycD1 repression factor\nK/(K+EZH2)', fontsize=10)
ax.set_title('F. EZH2 Repression Curve', fontsize=11, fontweight='bold', loc='left')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.2)
ax.set_ylim(0, 1.05)

# --- Row 3: Time courses comparing MB with different EZH2 levels ---
boost_show = [0.5, 1.0, 2.0, 4.0]
boost_colors = ['#3498db', '#2c3e50', '#e67e22', '#c0392b']

ax = fig.add_subplot(gs[2, 0])
for b, c in zip(boost_show, boost_colors):
    r = run_sim('mb', ezh2_boost=b)
    m = analyze(r)
    p_str = f"{m['period']:.1f}h" if m['period'] else "arr"
    ax.plot(r['time'], r['[Cb]'], color=c, linewidth=1, alpha=0.8,
            label=f'EZH2 {b}x ({m["n_div"]}div, {p_str})')
ax.set_ylabel('Cyclin B', fontsize=10)
ax.set_xlabel('Time (h)', fontsize=10)
ax.set_title('G. MB CycB: EZH2 Boost Comparison', fontsize=11, fontweight='bold', loc='left')
ax.legend(fontsize=7.5, loc='upper right')
ax.grid(True, alpha=0.2)

ax = fig.add_subplot(gs[2, 1])
for b, c in zip(boost_show, boost_colors):
    r = run_sim('mb', ezh2_boost=b)
    ax.plot(r['time'], r['[Cd_mRNA]'], color=c, linewidth=1, alpha=0.8,
            label=f'EZH2 {b}x')
ax.set_ylabel('CycD1 mRNA', fontsize=10)
ax.set_xlabel('Time (h)', fontsize=10)
ax.set_title('H. MB CycD1: EZH2 Boost Effect', fontsize=11, fontweight='bold', loc='left')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.2)

ax = fig.add_subplot(gs[2, 2])
for b, c in zip(boost_show, boost_colors):
    r = run_sim('mb', ezh2_boost=b)
    ax.plot(r['time'], r['[EZH2]'], color=c, linewidth=1, alpha=0.8,
            label=f'EZH2 synth {b}x')
ax.set_ylabel('EZH2 Protein', fontsize=10)
ax.set_xlabel('Time (h)', fontsize=10)
ax.set_title('I. MB EZH2 Protein Levels', fontsize=11, fontweight='bold', loc='left')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.2)

fig.suptitle('Does Higher EZH2 in MB Compensate for MYCN-Driven CycD1?',
             fontsize=14, fontweight='bold', y=1.01)
plt.savefig('simulations/fig_ezh2_compensation.png', dpi=200, bbox_inches='tight')
plt.savefig('simulations/fig_ezh2_compensation.pdf', bbox_inches='tight')
plt.close()
print("\nSaved: fig_ezh2_compensation.png/pdf")
