"""
EZH2 compensation with data-constrained levels.

RNA-seq data:
  P7_wt (GNP):  5,346
  MB:           9,450  → 1.77x higher

Current model only gives 1.13x (from faster cycling/E2F).
Gap: 1.77/1.13 = 1.57x unexplained.

Test: What's the impact of data-level EZH2 on CycD1 and cell cycle in MB?
Should we add MYCN→EZH2 transcription to close the gap?
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
K_EZH2 = 0.5


def run_sim(context='wt_shh', ezh2_boost=1.0, no_ezh2_fb=False,
            gdc=0.0, ezh2i=0.0, t_end=800, n_pts=8000):
    rr = te.loada(model_str)
    if context == 'wt_shh':
        rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 1.0
        rr['MYCN_amplification'] = 1.0
    elif context == 'mb':
        rr['SHH'] = 0.0; rr['Ptch1_copy_number'] = 0.0
        rr['MYCN_amplification'] = MYCN_AMP_MB
        rr['Ptch1_mRNA'] = 0.0; rr['Ptch1_free'] = 0.0
        rr['SHH_Ptch'] = 0.0; rr['Smo_active'] = 1.0
    rr['GDC0449'] = gdc; rr['EZH2i'] = ezh2i
    rr['k_EZH2_mRNA_synth_basal'] *= ezh2_boost
    if no_ezh2_fb:
        rr['K_EZH2_repression'] = 1e6
    return rr.simulate(0, t_end, n_pts)


def analyze(result, t_start=50):
    t = result['time']; cb = result['[Cb]']
    mask = t >= t_start
    peaks, _ = find_peaks(cb[mask], prominence=0.05, distance=10)
    peak_times = t[mask][peaks]
    periods = np.diff(peak_times) if len(peak_times) > 1 else np.array([])
    idx = -2000
    return {
        'n_div': len(peaks),
        'period': np.mean(periods) if len(periods) > 0 else None,
        'cd_mean': np.mean(result['[Cd_mRNA]'][idx:]),
        'ezh2_mean': np.mean(result['[EZH2]'][idx:]),
        'ezh2_mRNA_mean': np.mean(result['[EZH2_mRNA]'][idx:]),
        'e2f_mean': np.mean(result['[E2F]'][idx:]),
        'p27_mean': np.mean(result['[p27]'][idx:]),
        'mycn_mean': np.mean(result['[MYCN]'][idx:]),
        'gli1_mean': np.mean(result['[Gli1]'][idx:]),
    }


# =====================================================================
# Find the EZH2 boost that matches data ratio
# =====================================================================

print("=" * 90)
print("FINDING EZH2 BOOST TO MATCH DATA (MB/GNP EZH2 mRNA = 1.77x)")
print("=" * 90)

gnp_base = analyze(run_sim('wt_shh'))
gnp_ezh2_mRNA = gnp_base['ezh2_mRNA_mean']

# Binary search for the boost that gives MB EZH2 mRNA = 1.77x GNP
target_ratio = 9450 / 5346  # 1.767x
print(f"\nTarget: MB EZH2 mRNA / GNP EZH2 mRNA = {target_ratio:.3f}x")
print(f"GNP EZH2 mRNA = {gnp_ezh2_mRNA:.4f}")
print(f"Target MB EZH2 mRNA = {gnp_ezh2_mRNA * target_ratio:.4f}")

# Scan boost values
boosts_scan = np.arange(1.0, 3.5, 0.1)
ratios = []
for b in boosts_scan:
    m = analyze(run_sim('mb', ezh2_boost=b))
    ratio = m['ezh2_mRNA_mean'] / gnp_ezh2_mRNA
    ratios.append(ratio)

# Interpolate to find exact boost
from scipy.interpolate import interp1d
f_interp = interp1d(ratios, boosts_scan)
try:
    optimal_boost = float(f_interp(target_ratio))
    print(f"\nOptimal EZH2 synthesis boost: {optimal_boost:.2f}x")
except:
    # If target is beyond range, use closest
    optimal_boost = boosts_scan[np.argmin(np.abs(np.array(ratios) - target_ratio))]
    print(f"\nClosest EZH2 synthesis boost: {optimal_boost:.2f}x")

# Verify
mb_matched = analyze(run_sim('mb', ezh2_boost=optimal_boost))
print(f"Achieved MB EZH2 mRNA = {mb_matched['ezh2_mRNA_mean']:.4f}")
print(f"Ratio = {mb_matched['ezh2_mRNA_mean']/gnp_ezh2_mRNA:.3f}x (target: {target_ratio:.3f}x)")


# =====================================================================
# Compare: current model vs data-matched EZH2
# =====================================================================

print(f"\n{'=' * 90}")
print("COMPARISON: Current Model vs Data-Matched EZH2 in MB")
print("=" * 90)

mb_current = analyze(run_sim('mb', ezh2_boost=1.0))
mb_data = mb_matched  # already computed
mb_no_fb = analyze(run_sim('mb', no_ezh2_fb=True))

print(f"\n  {'Metric':30s} {'GNP+SHH':>10s} {'MB (curr)':>10s} {'MB (data)':>10s} {'MB (no fb)':>10s}")
print(f"  {'-'*75}")

def row(name, vg, v1, v2, v3, fmt='.2f'):
    vals = [vg, v1, v2, v3]
    strs = [f"{v:{fmt}}" if v is not None else "N/A" for v in vals]
    print(f"  {name:30s} {''.join(f'{s:>10s}' for s in strs)}")

row("EZH2 mRNA", gnp_base['ezh2_mRNA_mean'], mb_current['ezh2_mRNA_mean'],
    mb_data['ezh2_mRNA_mean'], mb_no_fb['ezh2_mRNA_mean'], fmt='.4f')
row("EZH2 mRNA ratio (vs GNP)",
    1.0, mb_current['ezh2_mRNA_mean']/gnp_ezh2_mRNA,
    mb_data['ezh2_mRNA_mean']/gnp_ezh2_mRNA, mb_no_fb['ezh2_mRNA_mean']/gnp_ezh2_mRNA)
row("EZH2 protein", gnp_base['ezh2_mean'], mb_current['ezh2_mean'],
    mb_data['ezh2_mean'], mb_no_fb['ezh2_mean'], fmt='.3f')
row("CycD1 mRNA", gnp_base['cd_mean'], mb_current['cd_mean'],
    mb_data['cd_mean'], mb_no_fb['cd_mean'])
row("Period (h)", gnp_base['period'], mb_current['period'],
    mb_data['period'], mb_no_fb['period'], fmt='.1f')
row("Divisions", gnp_base['n_div'], mb_current['n_div'],
    mb_data['n_div'], mb_no_fb['n_div'], fmt='d')
row("E2F (mean)", gnp_base['e2f_mean'], mb_current['e2f_mean'],
    mb_data['e2f_mean'], mb_no_fb['e2f_mean'], fmt='.1f')
row("MYCN (mean)", gnp_base['mycn_mean'], mb_current['mycn_mean'],
    mb_data['mycn_mean'], mb_no_fb['mycn_mean'], fmt='.3f')

# Quantify compensation
print(f"\n--- Compensation Analysis ---")
cd_no_fb = mb_no_fb['cd_mean']  # CycD1 with no EZH2 feedback at all
cd_curr = mb_current['cd_mean']  # CycD1 with current (model-predicted) EZH2
cd_data = mb_data['cd_mean']  # CycD1 with data-level EZH2

total_mycn_excess = cd_no_fb - gnp_base['cd_mean']  # How much MYCN adds over GNP
curr_compensation = cd_no_fb - cd_curr  # How much current EZH2 compensates
data_compensation = cd_no_fb - cd_data  # How much data-level EZH2 compensates
extra_compensation = cd_curr - cd_data  # Additional compensation from higher EZH2

print(f"  MYCN-driven CycD1 excess over GNP (no EZH2 fb):  {total_mycn_excess:.2f}")
print(f"  Current EZH2 compensates:                         {curr_compensation:.2f} ({curr_compensation/total_mycn_excess*100:.0f}%)")
print(f"  Data-level EZH2 compensates:                      {data_compensation:.2f} ({data_compensation/total_mycn_excess*100:.0f}%)")
print(f"  Extra compensation from data-level EZH2:          {extra_compensation:.2f} ({extra_compensation/total_mycn_excess*100:.0f}%)")
print(f"  Remaining MYCN excess (data-level):               {cd_data - gnp_base['cd_mean']:.2f}")

# =====================================================================
# How does this affect drug responses?
# =====================================================================

print(f"\n{'=' * 90}")
print("DRUG RESPONSES: Current vs Data-Matched EZH2")
print("=" * 90)

# MB + GDC
mb_gdc_curr = analyze(run_sim('mb', gdc=1.0))
mb_gdc_data = analyze(run_sim('mb', gdc=1.0, ezh2_boost=optimal_boost))

# MB + EZH2i
mb_ezh2i_curr = analyze(run_sim('mb', ezh2i=1.0))
mb_ezh2i_data = analyze(run_sim('mb', ezh2i=1.0, ezh2_boost=optimal_boost))

# MB + combo
mb_combo_curr = analyze(run_sim('mb', gdc=1.0, ezh2i=1.0))
mb_combo_data = analyze(run_sim('mb', gdc=1.0, ezh2i=1.0, ezh2_boost=optimal_boost))

print(f"\n  {'Condition':25s} {'CycD1 (curr)':>12s} {'CycD1 (data)':>12s} {'Div (curr)':>10s} {'Div (data)':>10s}")
print(f"  {'-'*75}")
print(f"  {'MB untreated':25s} {mb_current['cd_mean']:12.2f} {mb_data['cd_mean']:12.2f} {mb_current['n_div']:10d} {mb_data['n_div']:10d}")
print(f"  {'MB + GDC':25s} {mb_gdc_curr['cd_mean']:12.2f} {mb_gdc_data['cd_mean']:12.2f} {mb_gdc_curr['n_div']:10d} {mb_gdc_data['n_div']:10d}")
print(f"  {'MB + EZH2i':25s} {mb_ezh2i_curr['cd_mean']:12.2f} {mb_ezh2i_data['cd_mean']:12.2f} {mb_ezh2i_curr['n_div']:10d} {mb_ezh2i_data['n_div']:10d}")
print(f"  {'MB + GDC + EZH2i':25s} {mb_combo_curr['cd_mean']:12.2f} {mb_combo_data['cd_mean']:12.2f} {mb_combo_curr['n_div']:10d} {mb_combo_data['n_div']:10d}")

# GDC CycD1 ratio
cd_ratio_curr = mb_gdc_curr['cd_mean'] / mb_current['cd_mean']
cd_ratio_data = mb_gdc_data['cd_mean'] / mb_data['cd_mean']
print(f"\n  MB+GDC CycD1 ratio: current={cd_ratio_curr:.3f}, data-EZH2={cd_ratio_data:.3f} (target: 0.40)")

# EZH2i fold change
ezh2i_fc_curr = mb_ezh2i_curr['cd_mean'] / mb_current['cd_mean']
ezh2i_fc_data = mb_ezh2i_data['cd_mean'] / mb_data['cd_mean']
print(f"  MB EZH2i CycD1 FC:  current={ezh2i_fc_curr:.2f}x, data-EZH2={ezh2i_fc_data:.2f}x")

# Combo rescue
print(f"  MB combo divisions:  current={mb_combo_curr['n_div']}, data-EZH2={mb_combo_data['n_div']}")


# =====================================================================
# FIGURE
# =====================================================================

fig = plt.figure(figsize=(16, 12))
gs = GridSpec(3, 3, figure=fig, hspace=0.5, wspace=0.35)

# RNA-seq data
data_groups = ['P7_wt\n(GNP)', 'P7_Ptch\n(het)', 'P28_Ptch', 'MB']
data_vals = [5346, 8794, 7119, 9450]
data_colors = ['#1b9e77', '#66a61e', '#7570b3', '#762A83']

# Panel A: RNA-seq EZH2 data
ax = fig.add_subplot(gs[0, 0])
bars = ax.bar(range(4), data_vals, color=data_colors, edgecolor='black', alpha=0.8)
for bar, val in zip(bars, data_vals):
    ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 100,
            f'{val:,}', ha='center', va='bottom', fontsize=9, fontweight='bold')
ax.set_xticks(range(4))
ax.set_xticklabels(data_groups, fontsize=9)
ax.set_ylabel('EZH2 mRNA (counts)', fontsize=10)
ax.set_title('A. RNA-seq: EZH2 Expression', fontsize=11, fontweight='bold', loc='left')
ax.grid(True, alpha=0.2, axis='y')
# Add fold change annotations
ax.annotate('1.77x', xy=(3, 9450), xytext=(3.3, 7000),
            fontsize=9, fontweight='bold', color='#762A83',
            arrowprops=dict(arrowstyle='->', color='#762A83'))

# Panel B: Model EZH2 mRNA — current vs data-matched
ax = fig.add_subplot(gs[0, 1])
labels = ['GNP', 'MB\n(model)', 'MB\n(data-matched)']
vals = [gnp_base['ezh2_mRNA_mean'], mb_current['ezh2_mRNA_mean'], mb_data['ezh2_mRNA_mean']]
cols = ['#1b9e77', '#762A83', '#c0392b']
bars = ax.bar(range(3), vals, color=cols, edgecolor='black', alpha=0.8)
for i, (bar, val) in enumerate(zip(bars, vals)):
    ratio = val / gnp_base['ezh2_mRNA_mean']
    ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.001,
            f'{ratio:.2f}x', ha='center', va='bottom', fontsize=9, fontweight='bold')
ax.set_xticks(range(3))
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel('EZH2 mRNA (model units)', fontsize=10)
ax.set_title('B. Model EZH2 mRNA Levels', fontsize=11, fontweight='bold', loc='left')
ax.grid(True, alpha=0.2, axis='y')

# Panel C: CycD1 compensation waterfall
ax = fig.add_subplot(gs[0, 2])
# Waterfall: GNP level → MYCN excess → current EZH2 comp → extra EZH2 comp → final
waterfall_labels = ['GNP\nbaseline', 'MYCN\nexcess', 'EZH2 fb\n(model)', 'Extra EZH2\n(data)', 'MB\nfinal']
waterfall_vals = [gnp_base['cd_mean'], total_mycn_excess, -curr_compensation,
                  -extra_compensation, cd_data]
waterfall_bottoms = [0, gnp_base['cd_mean'], cd_no_fb, cd_curr, 0]
waterfall_colors = ['#1b9e77', '#e74c3c', '#8e44ad', '#c0392b', '#762A83']

for i in range(5):
    if i == 4:  # final bar from 0
        ax.bar(i, waterfall_vals[i], color=waterfall_colors[i], edgecolor='black',
               alpha=0.8, width=0.6)
    else:
        ax.bar(i, waterfall_vals[i], bottom=waterfall_bottoms[i],
               color=waterfall_colors[i], edgecolor='black', alpha=0.7, width=0.6)
    # Value label
    top = waterfall_bottoms[i] + waterfall_vals[i] if i < 4 else waterfall_vals[i]
    ax.text(i, top + 0.15, f'{waterfall_vals[i]:+.1f}' if i in [1,2,3] else f'{waterfall_vals[i]:.1f}',
            ha='center', fontsize=8, fontweight='bold')

ax.set_xticks(range(5))
ax.set_xticklabels(waterfall_labels, fontsize=8)
ax.set_ylabel('CycD1 mRNA', fontsize=10)
ax.set_title('C. CycD1 Compensation Waterfall', fontsize=11, fontweight='bold', loc='left')
ax.axhline(gnp_base['cd_mean'], color='#1b9e77', linestyle='--', alpha=0.4)
ax.grid(True, alpha=0.2, axis='y')

# Panel D: CycB traces — current vs data-matched
ax = fig.add_subplot(gs[1, 0])
r1 = run_sim('mb', ezh2_boost=1.0)
r2 = run_sim('mb', ezh2_boost=optimal_boost)
m1 = analyze(r1); m2 = analyze(r2)
ax.plot(r1['time'], r1['[Cb]'], color='#762A83', linewidth=1,
        label=f'Model EZH2 ({m1["n_div"]}div, {m1["period"]:.1f}h)')
ax.plot(r2['time'], r2['[Cb]'], color='#c0392b', linewidth=1, linestyle='--',
        label=f'Data EZH2 ({m2["n_div"]}div, {m2["period"]:.1f}h)')
ax.set_ylabel('Cyclin B', fontsize=10)
ax.set_xlabel('Time (h)', fontsize=10)
ax.set_title('D. MB CycB: Model vs Data-Level EZH2', fontsize=11, fontweight='bold', loc='left')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.2)

# Panel E: CycD1 traces
ax = fig.add_subplot(gs[1, 1])
ax.plot(r1['time'], r1['[Cd_mRNA]'], color='#762A83', linewidth=1,
        label=f'Model EZH2 (CycD1={m1["cd_mean"]:.1f})')
ax.plot(r2['time'], r2['[Cd_mRNA]'], color='#c0392b', linewidth=1, linestyle='--',
        label=f'Data EZH2 (CycD1={m2["cd_mean"]:.1f})')
ax.set_ylabel('CycD1 mRNA', fontsize=10)
ax.set_xlabel('Time (h)', fontsize=10)
ax.set_title('E. MB CycD1 mRNA', fontsize=11, fontweight='bold', loc='left')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.2)

# Panel F: EZH2 traces
ax = fig.add_subplot(gs[1, 2])
ax.plot(r1['time'], r1['[EZH2]'], color='#762A83', linewidth=1,
        label=f'Model (EZH2={m1["ezh2_mean"]:.2f})')
ax.plot(r2['time'], r2['[EZH2]'], color='#c0392b', linewidth=1, linestyle='--',
        label=f'Data-matched (EZH2={m2["ezh2_mean"]:.2f})')
ax.set_ylabel('EZH2 Protein', fontsize=10)
ax.set_xlabel('Time (h)', fontsize=10)
ax.set_title('F. MB EZH2 Protein', fontsize=11, fontweight='bold', loc='left')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.2)

# Panel G: Drug response comparison
ax = fig.add_subplot(gs[2, 0])
drug_labels = ['Untx', '+GDC', '+EZH2i', '+Both']
cd_curr_vals = [mb_current['cd_mean'], mb_gdc_curr['cd_mean'],
                mb_ezh2i_curr['cd_mean'], mb_combo_curr['cd_mean']]
cd_data_vals = [mb_data['cd_mean'], mb_gdc_data['cd_mean'],
                mb_ezh2i_data['cd_mean'], mb_combo_data['cd_mean']]
x = np.arange(4)
ax.bar(x - 0.15, cd_curr_vals, 0.3, color='#762A83', edgecolor='black',
       alpha=0.7, label='Model EZH2')
ax.bar(x + 0.15, cd_data_vals, 0.3, color='#c0392b', edgecolor='black',
       alpha=0.7, label='Data EZH2')
ax.set_xticks(x)
ax.set_xticklabels(drug_labels, fontsize=10)
ax.set_ylabel('CycD1 mRNA', fontsize=10)
ax.set_title('G. Drug Responses: CycD1', fontsize=11, fontweight='bold', loc='left')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.2, axis='y')

# Panel H: Divisions
ax = fig.add_subplot(gs[2, 1])
div_curr = [mb_current['n_div'], mb_gdc_curr['n_div'],
            mb_ezh2i_curr['n_div'], mb_combo_curr['n_div']]
div_data = [mb_data['n_div'], mb_gdc_data['n_div'],
            mb_ezh2i_data['n_div'], mb_combo_data['n_div']]
ax.bar(x - 0.15, div_curr, 0.3, color='#762A83', edgecolor='black',
       alpha=0.7, label='Model EZH2')
ax.bar(x + 0.15, div_data, 0.3, color='#c0392b', edgecolor='black',
       alpha=0.7, label='Data EZH2')
for i in range(4):
    ax.text(x[i]-0.15, div_curr[i]+0.3, str(div_curr[i]), ha='center', fontsize=8, fontweight='bold')
    ax.text(x[i]+0.15, div_data[i]+0.3, str(div_data[i]), ha='center', fontsize=8, fontweight='bold',
            color='#c0392b')
ax.set_xticks(x)
ax.set_xticklabels(drug_labels, fontsize=10)
ax.set_ylabel('Divisions (800h)', fontsize=10)
ax.set_title('H. Drug Responses: Divisions', fontsize=11, fontweight='bold', loc='left')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.2, axis='y')

# Panel I: Summary — what higher EZH2 means
ax = fig.add_subplot(gs[2, 2])
ax.axis('off')
summary_text = (
    f"Summary: EZH2 Compensation in MB\n"
    f"{'─' * 40}\n\n"
    f"RNA-seq: MB EZH2 = {target_ratio:.2f}x GNP\n"
    f"Model predicts: {mb_current['ezh2_mRNA_mean']/gnp_ezh2_mRNA:.2f}x (cycling alone)\n"
    f"Gap: {target_ratio/(mb_current['ezh2_mRNA_mean']/gnp_ezh2_mRNA):.2f}x unexplained\n"
    f"→ suggests MYCN or other drive on EZH2\n\n"
    f"EZH2 boost needed: {optimal_boost:.1f}x basal synth\n\n"
    f"Impact on CycD1 (data-level EZH2):\n"
    f"  CycD1: {mb_current['cd_mean']:.1f} → {mb_data['cd_mean']:.1f} "
    f"({(1-mb_data['cd_mean']/mb_current['cd_mean'])*100:.0f}% reduction)\n"
    f"  Period: {mb_current['period']:.1f}h → {mb_data['period']:.1f}h\n"
    f"  Divisions: {mb_current['n_div']} → {mb_data['n_div']}\n\n"
    f"Compensates {extra_compensation/total_mycn_excess*100:.0f}% of MYCN\n"
    f"excess CycD1 (on top of {curr_compensation/total_mycn_excess*100:.0f}%\n"
    f"from cycling-level EZH2)"
)
ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=9.5,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#fef9e7', alpha=0.8))
ax.set_title('I. Summary', fontsize=11, fontweight='bold', loc='left')

fig.suptitle('EZH2 Compensation for MYCN: Model vs RNA-seq Data',
             fontsize=14, fontweight='bold', y=1.01)
plt.savefig('simulations/fig_ezh2_data_constrained.png', dpi=200, bbox_inches='tight')
plt.savefig('simulations/fig_ezh2_data_constrained.pdf', bbox_inches='tight')
plt.close()
print("\nSaved: fig_ezh2_data_constrained.png/pdf")
