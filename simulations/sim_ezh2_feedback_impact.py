"""
Impact of EZH2 feedback on CycD1 in GNPs (WT+SHH) and MB (Ptch1-/-).

Compares normal model (EZH2 represses CycD1) vs no-feedback model
(EZH2 repression removed by setting K_EZH2_repression very high).
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


def run_sim(context='wt_shh', ezh2_feedback=True, t_end=168, n_pts=16800):
    """Run simulation with or without EZH2 feedback on CycD1."""
    rr = te.loada(model_str)

    # Context setup
    if context == 'wt_shh':
        rr['SHH'] = 0.5
        rr['Ptch1_copy_number'] = 1.0
        rr['MYCN_amplification'] = 1.0
    elif context == 'mb':
        rr['SHH'] = 0.0
        rr['Ptch1_copy_number'] = 0.0
        rr['MYCN_amplification'] = MYCN_AMP_MB
        rr['Ptch1_mRNA'] = 0.0
        rr['Ptch1_free'] = 0.0
        rr['SHH_Ptch'] = 0.0
        rr['Smo_active'] = 1.0

    rr['GDC0449'] = 0.0
    rr['EZH2i'] = 0.0

    # Remove EZH2 feedback by making K_EZH2_repression very large
    # This makes the repression factor ≈ 1.0 regardless of EZH2 level
    if not ezh2_feedback:
        rr['K_EZH2_repression'] = 1e6

    return rr.simulate(0, t_end, n_pts)


def analyze(result, t_start=10):
    """Extract key metrics from simulation."""
    t = result['time']
    cb = result['[Cb]']
    mask = t >= t_start

    peaks, _ = find_peaks(cb[mask], prominence=0.05, distance=10)
    peak_times = t[mask][peaks]
    periods = np.diff(peak_times) if len(peak_times) > 1 else np.array([])

    idx = -1500
    return {
        'n_div': len(peaks),
        'period': np.mean(periods) if len(periods) > 0 else None,
        'peak_times': peak_times,
        'periods': periods,
        'cd_mRNA_mean': np.mean(result['[Cd_mRNA]'][idx:]),
        'cd_mRNA_min': np.min(result['[Cd_mRNA]'][idx:]),
        'cd_mRNA_max': np.max(result['[Cd_mRNA]'][idx:]),
        'ezh2_mean': np.mean(result['[EZH2]'][idx:]),
        'ezh2_min': np.min(result['[EZH2]'][idx:]),
        'ezh2_max': np.max(result['[EZH2]'][idx:]),
        'cb_mean': np.mean(result['[Cb]'][idx:]),
        'cb_max': np.max(result['[Cb]'][idx:]),
        'e2f_mean': np.mean(result['[E2F]'][idx:]),
        'p27_mean': np.mean(result['[p27]'][idx:]),
        'prbpp_mean': np.mean(result['[pRBpp]'][idx:]),
        'mycn_mean': np.mean(result['[MYCN]'][idx:]),
    }


# =====================================================================
# RUN SIMULATIONS
# =====================================================================

contexts = {
    'GNP + SHH': 'wt_shh',
    'MB (Ptch1-/-)': 'mb',
}

results = {}
metrics = {}
for label, ctx in contexts.items():
    for fb in [True, False]:
        key = (label, fb)
        results[key] = run_sim(context=ctx, ezh2_feedback=fb)
        metrics[key] = analyze(results[key])

# =====================================================================
# PRINT COMPARISON TABLE
# =====================================================================

print("=" * 100)
print("IMPACT OF EZH2 FEEDBACK ON CycD1 — GNPs vs MB")
print("=" * 100)

for label in contexts:
    print(f"\n--- {label} ---")
    with_fb = metrics[(label, True)]
    no_fb = metrics[(label, False)]

    print(f"  {'Metric':28s} | {'With EZH2→CycD1':>16s} | {'No EZH2→CycD1':>16s} | {'Change':>10s}")
    print(f"  {'-'*80}")

    def row(name, v1, v2, fmt='.2f', pct=False):
        s1 = f"{v1:{fmt}}" if v1 is not None else "N/A"
        s2 = f"{v2:{fmt}}" if v2 is not None else "N/A"
        if v1 is not None and v2 is not None and v1 != 0:
            if pct:
                ch = f"{(v2/v1 - 1)*100:+.0f}%"
            else:
                ch = f"{v2/v1:.2f}x"
        else:
            ch = "—"
        print(f"  {name:28s} | {s1:>16s} | {s2:>16s} | {ch:>10s}")

    row("Divisions (168h)", with_fb['n_div'], no_fb['n_div'], fmt='d')
    row("Period (h)", with_fb['period'], no_fb['period'], fmt='.1f')
    row("CycD1 mRNA (mean)", with_fb['cd_mRNA_mean'], no_fb['cd_mRNA_mean'])
    row("CycD1 mRNA (min)", with_fb['cd_mRNA_min'], no_fb['cd_mRNA_min'])
    row("CycD1 mRNA (max)", with_fb['cd_mRNA_max'], no_fb['cd_mRNA_max'])
    row("CycD1 amplitude", with_fb['cd_mRNA_max'] - with_fb['cd_mRNA_min'],
        no_fb['cd_mRNA_max'] - no_fb['cd_mRNA_min'])
    row("EZH2 (mean)", with_fb['ezh2_mean'], no_fb['ezh2_mean'])
    row("EZH2 (min-max)", with_fb['ezh2_max'] - with_fb['ezh2_min'],
        no_fb['ezh2_max'] - no_fb['ezh2_min'])
    row("CycB peak", with_fb['cb_max'], no_fb['cb_max'])
    row("E2F (mean)", with_fb['e2f_mean'], no_fb['e2f_mean'])
    row("p27 (mean)", with_fb['p27_mean'], no_fb['p27_mean'])
    row("pRBpp (mean)", with_fb['prbpp_mean'], no_fb['prbpp_mean'])

# =====================================================================
# INTERPRETATION
# =====================================================================

print("\n" + "=" * 100)
print("INTERPRETATION")
print("=" * 100)

for label in contexts:
    wfb = metrics[(label, True)]
    nfb = metrics[(label, False)]
    print(f"\n--- {label} ---")

    if wfb['period'] and nfb['period']:
        dp = nfb['period'] - wfb['period']
        print(f"  Period change: {dp:+.1f}h ({wfb['period']:.1f}h → {nfb['period']:.1f}h)")
    elif wfb['period'] and not nfb['period']:
        print(f"  Cell cycle LOST without EZH2 feedback (was {wfb['period']:.1f}h)")
    elif not wfb['period'] and nfb['period']:
        print(f"  Cell cycle GAINED without EZH2 feedback ({nfb['period']:.1f}h)")

    cd_ratio = nfb['cd_mRNA_mean'] / wfb['cd_mRNA_mean'] if wfb['cd_mRNA_mean'] > 0 else 0
    print(f"  CycD1 mRNA: {cd_ratio:.2f}x without feedback")

    if wfb['cd_mRNA_max'] - wfb['cd_mRNA_min'] > 0:
        amp_ratio = (nfb['cd_mRNA_max'] - nfb['cd_mRNA_min']) / (wfb['cd_mRNA_max'] - wfb['cd_mRNA_min'])
        print(f"  CycD1 oscillation amplitude: {amp_ratio:.2f}x without feedback")

    # EZH2 repression factor in cycling vs G1
    # With feedback: repression = K/(K + EZH2) where K=0.5
    # When EZH2 is high (cycling): stronger repression
    # When EZH2 is low (G1): weaker repression
    K = 0.5
    rep_high = K / (K + wfb['ezh2_max'])
    rep_low = K / (K + wfb['ezh2_min'])
    print(f"  EZH2 repression factor range: {rep_high:.3f} (peak EZH2) to {rep_low:.3f} (trough EZH2)")
    print(f"  → CycD1 modulated {rep_low/rep_high:.1f}x by EZH2 oscillation within each cycle")


# =====================================================================
# FIGURE
# =====================================================================

fig = plt.figure(figsize=(18, 16))
gs = GridSpec(4, 2, figure=fig, hspace=0.45, wspace=0.3)

ctx_labels = list(contexts.keys())
colors = {'with': '#2c3e50', 'without': '#e74c3c'}
lstyles = {'with': '-', 'without': '--'}

# Compute shared CycD1 y-limit across both panels
cd_max_all = max(
    max(results[(l, fb)]['[Cd_mRNA]'].max() for fb in [True, False])
    for l in ctx_labels
)
cd_ylim = (0, cd_max_all * 1.08)

for col, label in enumerate(ctx_labels):
    r_with = results[(label, True)]
    r_without = results[(label, False)]
    m_with = metrics[(label, True)]
    m_without = metrics[(label, False)]

    t_w = r_with['time']
    t_wo = r_without['time']

    # Row 0: CycB (cell cycle oscillations)
    ax = fig.add_subplot(gs[0, col])
    ax.plot(t_w, r_with['[Cb]'], color=colors['with'], linewidth=1.2,
            label=f"With EZH2 fb ({m_with['n_div']} div, "
                  f"{m_with['period']:.1f}h)" if m_with['period'] else
                  f"With EZH2 fb ({m_with['n_div']} div)")
    ax.plot(t_wo, r_without['[Cb]'], color=colors['without'], linewidth=1.2,
            linestyle='--',
            label=f"No EZH2 fb ({m_without['n_div']} div, "
                  f"{m_without['period']:.1f}h)" if m_without['period'] else
                  f"No EZH2 fb ({m_without['n_div']} div)")
    ax.set_title(f'{label}', fontsize=13, fontweight='bold')
    ax.set_ylabel('Cyclin B (CycB)', fontsize=10)
    ax.set_xlabel('Time (h)', fontsize=10)
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.2)

    # Row 1: CycD1 mRNA (shared y-axis across panels)
    ax = fig.add_subplot(gs[1, col])
    ax.plot(t_w, r_with['[Cd_mRNA]'], color=colors['with'], linewidth=1.2,
            label='With EZH2→CycD1')
    ax.plot(t_wo, r_without['[Cd_mRNA]'], color=colors['without'], linewidth=1.2,
            linestyle='--', label='No EZH2→CycD1')
    ax.set_ylabel('CycD1 mRNA', fontsize=10)
    ax.set_xlabel('Time (h)', fontsize=10)
    ax.set_title(f'{label}: CycD1 mRNA', fontsize=11, fontweight='bold')
    ax.set_ylim(cd_ylim)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.2)

    # Row 2: EZH2 protein
    ax = fig.add_subplot(gs[2, col])
    ax.plot(t_w, r_with['[EZH2]'], color=colors['with'], linewidth=1.2,
            label='With EZH2→CycD1')
    ax.plot(t_wo, r_without['[EZH2]'], color=colors['without'], linewidth=1.2,
            linestyle='--', label='No EZH2→CycD1')
    ax.set_ylabel('EZH2 Protein', fontsize=10)
    ax.set_xlabel('Time (h)', fontsize=10)
    ax.set_title(f'{label}: EZH2 Dynamics', fontsize=11, fontweight='bold')
    ax.legend(fontsize=8)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.2)

    # Row 3: Overlay CycD1 and EZH2 (normalized) to show phase relationship
    ax = fig.add_subplot(gs[3, col])
    # Zoom into last 70h for clarity
    t_zoom = 70
    mask_w = t_w >= (t_w[-1] - t_zoom)
    mask_wo = t_wo >= (t_wo[-1] - t_zoom)

    # Normalize for overlay
    cd_w = r_with['[Cd_mRNA]'][mask_w]
    ezh2_w = r_with['[EZH2]'][mask_w]
    cb_w = r_with['[Cb]'][mask_w]
    t_plot = t_w[mask_w]

    cd_norm = (cd_w - cd_w.min()) / (cd_w.max() - cd_w.min() + 1e-10)
    ezh2_norm = (ezh2_w - ezh2_w.min()) / (ezh2_w.max() - ezh2_w.min() + 1e-10)
    cb_norm = (cb_w - cb_w.min()) / (cb_w.max() - cb_w.min() + 1e-10)

    ax.plot(t_plot, cd_norm, color='#e67e22', linewidth=1.5, label='CycD1 mRNA')
    ax.plot(t_plot, ezh2_norm, color='#8e44ad', linewidth=1.5, label='EZH2')
    ax.plot(t_plot, cb_norm, color='#2980b9', linewidth=1.2, alpha=0.6, label='CycB')
    ax.set_ylabel('Normalized (0–1)', fontsize=10)
    ax.set_xlabel('Time (h)', fontsize=10)
    ax.set_title(f'{label}: Phase Relationship (with feedback, last {t_zoom}h)',
                 fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.2)

fig.suptitle('Impact of EZH2→CycD1 Feedback on Cell Cycle Dynamics',
             fontsize=15, fontweight='bold', y=1.01)
plt.savefig('simulations/fig_ezh2_feedback_impact.png', dpi=200, bbox_inches='tight')
plt.savefig('simulations/fig_ezh2_feedback_impact.pdf', bbox_inches='tight')
plt.close()
print("\nSaved: fig_ezh2_feedback_impact.png/pdf")


# =====================================================================
# ADDITIONAL: Period and CycD1 bar comparison
# =====================================================================

fig2, axes = plt.subplots(1, 3, figsize=(14, 4.5))
fig2.suptitle('EZH2 Feedback Impact: Key Metrics', fontsize=13, fontweight='bold')

# Panel 1: Period comparison
ax = axes[0]
x = np.arange(len(ctx_labels))
per_w = [metrics[(l, True)]['period'] or 0 for l in ctx_labels]
per_wo = [metrics[(l, False)]['period'] or 0 for l in ctx_labels]
ax.bar(x - 0.15, per_w, 0.3, color=colors['with'], edgecolor='black', alpha=0.8, label='With EZH2 fb')
ax.bar(x + 0.15, per_wo, 0.3, color=colors['without'], edgecolor='black', alpha=0.8, label='No EZH2 fb')
for i in range(len(ctx_labels)):
    ax.text(x[i]-0.15, per_w[i]+0.3, f'{per_w[i]:.1f}h', ha='center', fontsize=8, fontweight='bold')
    ax.text(x[i]+0.15, per_wo[i]+0.3, f'{per_wo[i]:.1f}h', ha='center', fontsize=8, fontweight='bold', color=colors['without'])
ax.set_xticks(x)
ax.set_xticklabels(ctx_labels, fontsize=10)
ax.set_ylabel('Cell Cycle Period (h)', fontsize=10)
ax.set_title('Period', fontsize=11, fontweight='bold')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.2, axis='y')

# Panel 2: CycD1 comparison
ax = axes[1]
cd_w = [metrics[(l, True)]['cd_mRNA_mean'] for l in ctx_labels]
cd_wo = [metrics[(l, False)]['cd_mRNA_mean'] for l in ctx_labels]
ax.bar(x - 0.15, cd_w, 0.3, color=colors['with'], edgecolor='black', alpha=0.8, label='With EZH2 fb')
ax.bar(x + 0.15, cd_wo, 0.3, color=colors['without'], edgecolor='black', alpha=0.8, label='No EZH2 fb')
for i in range(len(ctx_labels)):
    ax.text(x[i]-0.15, cd_w[i]+0.1, f'{cd_w[i]:.1f}', ha='center', fontsize=8, fontweight='bold')
    ax.text(x[i]+0.15, cd_wo[i]+0.1, f'{cd_wo[i]:.1f}', ha='center', fontsize=8, fontweight='bold', color=colors['without'])
ax.set_xticks(x)
ax.set_xticklabels(ctx_labels, fontsize=10)
ax.set_ylabel('CycD1 mRNA (mean)', fontsize=10)
ax.set_title('CycD1 mRNA', fontsize=11, fontweight='bold')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.2, axis='y')

# Panel 3: Divisions
ax = axes[2]
div_w = [metrics[(l, True)]['n_div'] for l in ctx_labels]
div_wo = [metrics[(l, False)]['n_div'] for l in ctx_labels]
ax.bar(x - 0.15, div_w, 0.3, color=colors['with'], edgecolor='black', alpha=0.8, label='With EZH2 fb')
ax.bar(x + 0.15, div_wo, 0.3, color=colors['without'], edgecolor='black', alpha=0.8, label='No EZH2 fb')
for i in range(len(ctx_labels)):
    ax.text(x[i]-0.15, div_w[i]+0.1, str(div_w[i]), ha='center', fontsize=9, fontweight='bold')
    ax.text(x[i]+0.15, div_wo[i]+0.1, str(div_wo[i]), ha='center', fontsize=9, fontweight='bold', color=colors['without'])
ax.set_xticks(x)
ax.set_xticklabels(ctx_labels, fontsize=10)
ax.set_ylabel('Divisions (168h)', fontsize=10)
ax.set_title('Cell Divisions', fontsize=11, fontweight='bold')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.2, axis='y')

plt.tight_layout()
plt.savefig('simulations/fig_ezh2_feedback_periods.png', dpi=200, bbox_inches='tight')
plt.close()
print("Saved: fig_ezh2_feedback_periods.png")
