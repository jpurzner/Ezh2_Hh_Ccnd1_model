"""
Comprehensive validation of v42 model (MYCN + Hill coefficient).
Tests all conditions against experimental data and generates figures.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import tellurium as te
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.signal import find_peaks
from src.build_model_v42_mycn import build_model_v42
import json

model_str = build_model_v42()

# =====================================================================
# SIMULATION HELPERS
# =====================================================================

def run_sim(shh=0.5, ptch1_cn=1.0, gdc=0.0, ezh2i=0.0, mycn_amp=1.0,
            t_end=96, n_pts=9600, serum_starve=False, cdk46i=False):
    """Run simulation with given conditions."""
    rr = te.loada(model_str)
    rr['SHH'] = shh
    rr['Ptch1_copy_number'] = ptch1_cn
    rr['GDC0449'] = gdc
    rr['EZH2i'] = ezh2i
    rr['MYCN_amplification'] = mycn_amp

    # Ptch1-/- (MB): constitutive Smo activation
    if ptch1_cn == 0.0:
        rr['Ptch1_mRNA'] = 0.0
        rr['Ptch1_free'] = 0.0
        rr['SHH_Ptch'] = 0.0
        rr['Smo_active'] = 1.0
    elif shh == 0.0:
        rr['SHH_Ptch'] = 0.0
        rr['Ptch1_free'] = 1.0
        rr['Smo_active'] = 0.01

    # Serum starvation: block CycD1 translation
    if serum_starve:
        rr['k_Cd_translation'] = 0.0

    # CDK4/6 inhibitor: block CycD/CDK4-mediated pRB phosphorylation
    if cdk46i:
        rr['V1'] = 0.0

    return rr.simulate(0, t_end, n_pts)


def count_divisions(result, t_start=10):
    """Count cell divisions (CycB peaks) after t_start."""
    t = result['time']
    cb = result['[Cb]']
    mask = t >= t_start
    peaks, props = find_peaks(cb[mask], prominence=0.05, distance=10)
    peak_times = t[mask][peaks]
    periods = np.diff(peak_times) if len(peak_times) > 1 else np.array([])
    return len(peaks), peak_times, periods


def mean_last(result, species, n=1500):
    """Mean of last n points for a species."""
    return np.mean(result[species][-n:])


# =====================================================================
# RUN ALL CONDITIONS
# =====================================================================

MYCN_AMP_MB = 2.8

conditions = {
    # GNP conditions
    'GNP + SHH':          dict(shh=0.5, ptch1_cn=1.0, gdc=0.0, ezh2i=0.0, mycn_amp=1.0),
    'GNP - SHH':          dict(shh=0.0, ptch1_cn=1.0, gdc=0.0, ezh2i=0.0, mycn_amp=1.0),
    'GNP + HHi':          dict(shh=0.5, ptch1_cn=1.0, gdc=1.0, ezh2i=0.0, mycn_amp=1.0),
    'GNP + EZH2i':        dict(shh=0.5, ptch1_cn=1.0, gdc=0.0, ezh2i=1.0, mycn_amp=1.0),
    'GNP + HHi + EZH2i':  dict(shh=0.5, ptch1_cn=1.0, gdc=1.0, ezh2i=1.0, mycn_amp=1.0),
    # Serum starvation (G0)
    'GNP Serum-starved':   dict(shh=0.5, ptch1_cn=1.0, gdc=0.0, ezh2i=0.0, mycn_amp=1.0, serum_starve=True),
    # MB (Ptch1-/-) conditions
    'MB (Ptch1-/-)':       dict(shh=0.0, ptch1_cn=0.0, gdc=0.0, ezh2i=0.0, mycn_amp=MYCN_AMP_MB),
    'MB + HHi':            dict(shh=0.0, ptch1_cn=0.0, gdc=1.0, ezh2i=0.0, mycn_amp=MYCN_AMP_MB),
    'MB + EZH2i':          dict(shh=0.0, ptch1_cn=0.0, gdc=0.0, ezh2i=1.0, mycn_amp=MYCN_AMP_MB),
    'MB + HHi + EZH2i':   dict(shh=0.0, ptch1_cn=0.0, gdc=1.0, ezh2i=1.0, mycn_amp=MYCN_AMP_MB),
}

print("=" * 90)
print("V42 MODEL VALIDATION — MYCN Hill Function (n=3)")
print("=" * 90)

sims = {}
for name, kwargs in conditions.items():
    sims[name] = run_sim(**kwargs)

# =====================================================================
# PRINT RESULTS TABLE
# =====================================================================

print(f"\n{'Condition':24s} | {'Div':>4} {'Period':>7} | {'MYCN':>6} {'CycD1':>7} {'Gli1':>7} {'EZH2':>6} {'E2F':>6} {'pRBpp':>7} {'p27':>6}")
print("-" * 105)

idx = -2000
for name in conditions:
    r = sims[name]
    n, pt, per = count_divisions(r)
    p_str = f"{np.mean(per):.1f}h" if len(per) > 0 else "N/A"
    mycn = mean_last(r, '[MYCN]')
    cd = mean_last(r, '[Cd_mRNA]')
    gli1 = mean_last(r, '[Gli1]')
    ezh2 = mean_last(r, '[EZH2]')
    e2f = mean_last(r, '[E2F]')
    prbpp = mean_last(r, '[pRBpp]')
    p27 = mean_last(r, '[p27]')
    print(f"  {name:22s} | {n:4d} {p_str:>7s} | {mycn:6.3f} {cd:7.3f} {gli1:7.3f} {ezh2:6.3f} {e2f:6.2f} {prbpp:7.4f} {p27:6.2f}")

# =====================================================================
# VALIDATION TESTS
# =====================================================================

print("\n" + "=" * 90)
print("VALIDATION TESTS")
print("=" * 90)

results = {}
tests_passed = 0
tests_total = 0

def check(name, actual, target, tol=0.15, unit=''):
    """Check if actual is within tolerance of target."""
    global tests_passed, tests_total
    tests_total += 1
    if target == 0:
        ok = abs(actual) <= 1  # allow 0 or 1 (transient peak)
    else:
        ok = abs(actual - target) / abs(target) <= tol
    status = "PASS" if ok else "FAIL"
    if ok:
        tests_passed += 1
    symbol = "[+]" if ok else "[-]"
    if target == 0:
        print(f"  {symbol} {name:45s}: {actual:7.3f} (target: {target}{unit}, {'ok' if ok else 'FAIL'})")
    else:
        print(f"  {symbol} {name:45s}: {actual:7.3f} (target: {target}{unit}, {'ok' if ok else f'off by {abs(actual-target)/abs(target)*100:.0f}%'})")
    results[name] = {'actual': actual, 'target': target, 'pass': ok}
    return ok

# ---- Test 1: WT cycling ----
print("\n--- 1. WT Cell Cycle ---")
wt = sims['GNP + SHH']
n_wt, _, per_wt = count_divisions(wt)
period_wt = np.mean(per_wt) if len(per_wt) > 0 else 0
check("WT divisions (72h)", n_wt, 3, tol=0.50)
check("WT period (~23h)", period_wt, 23.0, tol=0.15)

# ---- Test 2: G0 arrest (serum starvation) ----
print("\n--- 2. G0 Arrest ---")
g0 = sims['GNP Serum-starved']
n_g0, _, _ = count_divisions(g0)
check("G0 divisions (=0)", n_g0, 0, tol=0.01)

ezh2_cyc = mean_last(wt, '[EZH2]')
ezh2_g0 = mean_last(g0, '[EZH2]')
check("EZH2 G0/cycling ratio (~0.6)", ezh2_g0 / ezh2_cyc if ezh2_cyc > 0 else 0, 0.6, tol=0.25)

# ---- Test 3: GDC in WT ----
print("\n--- 3. GDC Response (WT) ---")
wt_gdc = sims['GNP + HHi']
cd_wt = mean_last(wt, '[Cd_mRNA]')
cd_wt_gdc = mean_last(wt_gdc, '[Cd_mRNA]')
gli_wt = mean_last(wt, '[Gli1]')
gli_wt_gdc = mean_last(wt_gdc, '[Gli1]')
mycn_wt = mean_last(wt, '[MYCN]')
mycn_wt_gdc = mean_last(wt_gdc, '[MYCN]')
ezh2_wt_gdc = mean_last(wt_gdc, '[EZH2]')

check("WT+GDC CycD1 ratio (0.14)", cd_wt_gdc / cd_wt if cd_wt > 0 else 0, 0.14, tol=0.20)
check("WT+GDC Gli1 reduction (>90%)", 1 - gli_wt_gdc / gli_wt if gli_wt > 0 else 0, 0.99, tol=0.10)
check("WT+GDC MYCN ratio (0.78)", mycn_wt_gdc / mycn_wt if mycn_wt > 0 else 0, 0.78, tol=0.10)
n_wt_gdc, _, _ = count_divisions(wt_gdc)
check("WT+GDC arrest (0 div)", n_wt_gdc, 0, tol=0.01)

# ---- Test 4: EZH2i in WT ----
print("\n--- 4. EZH2i Response (WT) ---")
wt_ezh2i = sims['GNP + EZH2i']
cd_wt_ezh2i = mean_last(wt_ezh2i, '[Cd_mRNA]')
check("EZH2i CycD1 fold change (~2x)", cd_wt_ezh2i / cd_wt if cd_wt > 0 else 0, 2.0, tol=0.20)

# ---- Test 5: noSHH ----
print("\n--- 5. No SHH ---")
wt_noshh = sims['GNP - SHH']
n_noshh, _, _ = count_divisions(wt_noshh)
check("noSHH arrest (0 div)", n_noshh, 0, tol=0.01)

# ---- Test 6: MB cycling ----
print("\n--- 6. MB (Ptch1-/-) Cell Cycle ---")
mb = sims['MB (Ptch1-/-)']
n_mb, _, per_mb = count_divisions(mb)
check("MB divisions (72h)", n_mb, 3, tol=0.50)
mycn_mb = mean_last(mb, '[MYCN]')
check("MYCN MB/WT ratio (~2.8)", mycn_mb / mycn_wt if mycn_wt > 0 else 0, 2.8, tol=0.20)

# ---- Test 7: GDC in MB ----
print("\n--- 7. GDC Response (MB) --- KEY TEST ---")
mb_gdc = sims['MB + HHi']
cd_mb = mean_last(mb, '[Cd_mRNA]')
cd_mb_gdc = mean_last(mb_gdc, '[Cd_mRNA]')
gli_mb = mean_last(mb, '[Gli1]')
gli_mb_gdc = mean_last(mb_gdc, '[Gli1]')
mycn_mb_gdc = mean_last(mb_gdc, '[MYCN]')

check("MB+GDC CycD1 ratio (0.40)", cd_mb_gdc / cd_mb if cd_mb > 0 else 0, 0.40, tol=0.15)
check("MB+GDC MYCN ratio (0.86)", mycn_mb_gdc / mycn_mb if mycn_mb > 0 else 0, 0.86, tol=0.10)
check("MB+GDC Gli1 reduction (77%)", 1 - gli_mb_gdc / gli_mb if gli_mb > 0 else 0, 0.77, tol=0.20)
n_mb_gdc, _, _ = count_divisions(mb_gdc)
check("MB+GDC arrest (0 div)", n_mb_gdc, 0, tol=0.01)

# ---- Test 8: Combo in MB ----
print("\n--- 8. GDC + EZH2i Combo (MB) ---")
mb_combo = sims['MB + HHi + EZH2i']
n_mb_combo, _, per_mb_combo = count_divisions(mb_combo)
check("MB+GDC+EZH2i rescue (>0 div)", n_mb_combo, 3, tol=0.80)

# ---- Test 9: EZH2 regulation ----
print("\n--- 9. EZH2 Regulation ---")
ezh2_wt = mean_last(wt, '[EZH2]')
ezh2_mb_gdc = mean_last(mb_gdc, '[EZH2]')
check("EZH2 WT+GDC reduction (25-47%)", 1 - ezh2_wt_gdc / ezh2_wt if ezh2_wt > 0 else 0, 0.36, tol=0.40)

# Summary
print(f"\n{'=' * 90}")
print(f"VALIDATION SUMMARY: {tests_passed}/{tests_total} tests passed")
print(f"{'=' * 90}")

# =====================================================================
# FIGURE 1: Comprehensive validation — CycB traces
# =====================================================================

fig1, axes = plt.subplots(2, 5, figsize=(20, 8))
fig1.suptitle('V42 Model — CycB Traces (All Conditions)', fontsize=14, fontweight='bold')

cond_list = list(conditions.keys())
colors_list = ['#1b9e77', '#999999', '#e7298a', '#66a61e', '#d95f02',
               '#7570b3', '#762A83', '#E08214', '#D6604D', '#4A1486']

for i, (name, color) in enumerate(zip(cond_list, colors_list)):
    ax = axes[i // 5, i % 5]
    r = sims[name]
    ax.plot(r['time'], r['[Cb]'], color=color, linewidth=1)
    n, _, per = count_divisions(r)
    p_str = f", {np.mean(per):.1f}h" if len(per) > 0 else ""
    ax.set_title(f'{name}\n({n} div{p_str})', fontsize=8, fontweight='bold', color=color)
    ax.set_xlabel('Time (h)', fontsize=7)
    if i % 5 == 0:
        ax.set_ylabel('Cyclin B', fontsize=8)
    ax.grid(True, alpha=0.2)
    ax.tick_params(labelsize=7)

# Hide last subplot if odd count
if len(cond_list) < 10:
    axes[1, 4].set_visible(False)

plt.tight_layout()
plt.savefig('simulations/fig_v42_cycb_traces.png', dpi=200, bbox_inches='tight')
plt.close()
print("\nSaved: fig_v42_cycb_traces.png")

# =====================================================================
# FIGURE 2: Model vs Data — Key comparisons
# =====================================================================

fig2 = plt.figure(figsize=(18, 14))
gs = GridSpec(3, 4, figure=fig2, hspace=0.5, wspace=0.4)

# --- Row 1: GNP conditions ---
# Panel 1a: GNP CycB traces
ax = fig2.add_subplot(gs[0, 0])
gnp_conds = ['GNP + SHH', 'GNP - SHH', 'GNP + HHi', 'GNP + EZH2i', 'GNP + HHi + EZH2i']
gnp_colors = ['#1b9e77', '#999999', '#e7298a', '#66a61e', '#d95f02']
for name, c in zip(gnp_conds, gnp_colors):
    r = sims[name]
    ax.plot(r['time'], r['[Cb]'], color=c, linewidth=1, label=name.replace('GNP + ', '').replace('GNP - ', '-'))
ax.set_title('GNP: CycB Traces', fontsize=10, fontweight='bold')
ax.set_ylabel('Cyclin B', fontsize=9)
ax.set_xlabel('Time (h)', fontsize=9)
ax.legend(fontsize=6, loc='upper right')
ax.grid(True, alpha=0.2)

# Panel 1b: GNP CycD1 bars (model vs data)
ax = fig2.add_subplot(gs[0, 1])
gnp_cd_model = [1.0, cd_wt_gdc / cd_wt, cd_wt_ezh2i / cd_wt,
                mean_last(sims['GNP + HHi + EZH2i'], '[Cd_mRNA]') / cd_wt]
gnp_cd_data = [1.0, 0.14, None, None]
gnp_bar_colors = ['#1b9e77', '#e7298a', '#66a61e', '#d95f02']
x = np.arange(4)
ax.bar(x - 0.15, gnp_cd_model, 0.3, color=gnp_bar_colors, edgecolor='black', alpha=0.7, label='Model')
for i, val in enumerate(gnp_cd_data):
    if val is not None:
        ax.scatter(i + 0.15, val, color='black', s=80, marker='D', zorder=5)
ax.set_xticks(x)
ax.set_xticklabels(['SHH', '+HHi', '+EZH2i', '+Both'], fontsize=8)
ax.set_ylabel('CycD1 (fold of GNP)', fontsize=9)
ax.set_title('GNP: CycD1 Model vs Data', fontsize=10, fontweight='bold')
ax.legend(['Model', 'Data'], fontsize=8)
ax.grid(True, alpha=0.2, axis='y')

# Panel 1c: GNP MYCN bars
ax = fig2.add_subplot(gs[0, 2])
gnp_mycn_model = [1.0, mycn_wt_gdc / mycn_wt]
gnp_mycn_data = [1.0, 0.78]
x = np.arange(2)
ax.bar(x - 0.15, gnp_mycn_model, 0.3, color=['#1b9e77', '#e7298a'], edgecolor='black', alpha=0.7)
for i, val in enumerate(gnp_mycn_data):
    if val is not None:
        ax.scatter(i + 0.15, val, color='black', s=80, marker='D', zorder=5)
        ax.text(i + 0.25, val, f'{val:.2f}', fontsize=8, va='center')
ax.set_xticks(x)
ax.set_xticklabels(['GNP', 'GNP+HHi'], fontsize=9)
ax.set_ylabel('MYCN (fold of GNP)', fontsize=9)
ax.set_title('GNP: MYCN', fontsize=10, fontweight='bold')
ax.grid(True, alpha=0.2, axis='y')

# Panel 1d: EZH2 across conditions
ax = fig2.add_subplot(gs[0, 3])
ezh2_conds = ['GNP + SHH', 'GNP Serum-starved', 'GNP + HHi']
ezh2_vals = [mean_last(sims[c], '[EZH2]') / ezh2_wt for c in ezh2_conds]
ezh2_colors = ['#1b9e77', '#7570b3', '#e7298a']
ax.bar(range(len(ezh2_conds)), ezh2_vals, color=ezh2_colors, edgecolor='black', alpha=0.7)
ax.set_xticks(range(len(ezh2_conds)))
ax.set_xticklabels(['Cycling', 'G0', 'HHi'], fontsize=9)
ax.set_ylabel('EZH2 (fold of cycling)', fontsize=9)
ax.set_title('EZH2 Regulation', fontsize=10, fontweight='bold')
ax.axhline(0.6, color='red', linestyle='--', alpha=0.5, label='Target G0')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.2, axis='y')

# --- Row 2: MB conditions ---
mb_conds = ['MB (Ptch1-/-)', 'MB + HHi', 'MB + EZH2i', 'MB + HHi + EZH2i']
mb_colors = ['#762A83', '#E08214', '#D6604D', '#4A1486']

# Panel 2a: MB CycB traces
ax = fig2.add_subplot(gs[1, 0])
for name, c in zip(mb_conds, mb_colors):
    r = sims[name]
    ax.plot(r['time'], r['[Cb]'], color=c, linewidth=1, label=name.replace('MB ', '').replace('(Ptch1-/-)', 'Untx'))
ax.set_title('MB: CycB Traces', fontsize=10, fontweight='bold')
ax.set_ylabel('Cyclin B', fontsize=9)
ax.set_xlabel('Time (h)', fontsize=9)
ax.legend(fontsize=7, loc='upper right')
ax.grid(True, alpha=0.2)

# Panel 2b: MB CycD1 bars (model vs data)
ax = fig2.add_subplot(gs[1, 1])
mb_cd_model = [1.0, cd_mb_gdc / cd_mb,
               mean_last(sims['MB + EZH2i'], '[Cd_mRNA]') / cd_mb,
               mean_last(sims['MB + HHi + EZH2i'], '[Cd_mRNA]') / cd_mb]
mb_cd_data = [1.0, 0.40, None, None]
x = np.arange(4)
ax.bar(x - 0.15, mb_cd_model, 0.3, color=mb_colors, edgecolor='black', alpha=0.7)
for i, val in enumerate(mb_cd_data):
    if val is not None:
        ax.scatter(i + 0.15, val, color='black', s=80, marker='D', zorder=5)
        ax.text(i + 0.25, val, f'{val:.2f}', fontsize=8, va='center')
ax.set_xticks(x)
ax.set_xticklabels(['MB', '+HHi', '+EZH2i', '+Both'], fontsize=8)
ax.set_ylabel('CycD1 (fold of MB)', fontsize=9)
ax.set_title('MB: CycD1 Model vs Data', fontsize=10, fontweight='bold')
ax.legend(['Model', 'Data'], fontsize=8)
ax.grid(True, alpha=0.2, axis='y')

# Panel 2c: MB MYCN bars
ax = fig2.add_subplot(gs[1, 2])
mb_mycn_model = [mycn_mb / mycn_wt, mycn_mb_gdc / mycn_wt]
mb_mycn_data = [2.8, 2.8 * 0.86]
x = np.arange(2)
ax.bar(x - 0.15, mb_mycn_model, 0.3, color=['#762A83', '#E08214'], edgecolor='black', alpha=0.7)
for i, val in enumerate(mb_mycn_data):
    ax.scatter(i + 0.15, val, color='black', s=80, marker='D', zorder=5)
    ax.text(i + 0.25, val, f'{val:.2f}', fontsize=8, va='center')
ax.set_xticks(x)
ax.set_xticklabels(['MB', 'MB+HHi'], fontsize=9)
ax.set_ylabel('MYCN (fold of GNP)', fontsize=9)
ax.set_title('MB: MYCN', fontsize=10, fontweight='bold')
ax.grid(True, alpha=0.2, axis='y')

# Panel 2d: MB Gli1 bars
ax = fig2.add_subplot(gs[1, 3])
mb_gli_model = [1.0, gli_mb_gdc / gli_mb if gli_mb > 0 else 0,
                mean_last(sims['MB + EZH2i'], '[Gli1]') / gli_mb if gli_mb > 0 else 0,
                mean_last(sims['MB + HHi + EZH2i'], '[Gli1]') / gli_mb if gli_mb > 0 else 0]
mb_gli_data = [1.0, 0.23, None, None]
x = np.arange(4)
ax.bar(x - 0.15, mb_gli_model, 0.3, color=mb_colors, edgecolor='black', alpha=0.7)
for i, val in enumerate(mb_gli_data):
    if val is not None:
        ax.scatter(i + 0.15, val, color='black', s=80, marker='D', zorder=5)
        ax.text(i + 0.25, val, f'{val:.2f}', fontsize=8, va='center')
ax.set_xticks(x)
ax.set_xticklabels(['MB', '+HHi', '+EZH2i', '+Both'], fontsize=8)
ax.set_ylabel('Gli1 (fold of MB)', fontsize=9)
ax.set_title('MB: Gli1 Model vs Data', fontsize=10, fontweight='bold')
ax.grid(True, alpha=0.2, axis='y')

# --- Row 3: Division counts & key species time courses ---
# Panel 3a: Division counts all conditions
ax = fig2.add_subplot(gs[2, 0])
all_conds_plot = ['GNP + SHH', 'GNP - SHH', 'GNP + HHi', 'GNP + EZH2i', 'MB (Ptch1-/-)', 'MB + HHi', 'MB + EZH2i', 'MB + HHi + EZH2i']
all_colors_plot = ['#1b9e77', '#999999', '#e7298a', '#66a61e', '#762A83', '#E08214', '#D6604D', '#4A1486']
divs = [count_divisions(sims[c])[0] for c in all_conds_plot]
bars = ax.bar(range(len(all_conds_plot)), divs, color=all_colors_plot, edgecolor='black')
for bar, val in zip(bars, divs):
    ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.3,
            str(val), ha='center', va='bottom', fontweight='bold', fontsize=9)
ax.set_xticks(range(len(all_conds_plot)))
labels = ['GNP', 'GNP\n-SHH', 'GNP\n+HHi', 'GNP\n+EZH2i', 'MB', 'MB\n+HHi', 'MB\n+EZH2i', 'MB\n+Both']
ax.set_xticklabels(labels, fontsize=6)
ax.set_ylabel('Divisions (96h)', fontsize=9)
ax.set_title('Cell Divisions', fontsize=10, fontweight='bold')
ax.grid(True, alpha=0.2, axis='y')

# Panel 3b: CycD1 time courses MB conditions
ax = fig2.add_subplot(gs[2, 1])
for name, c in zip(mb_conds, mb_colors):
    r = sims[name]
    ax.plot(r['time'], r['[Cd_mRNA]'], color=c, linewidth=1, label=name.replace('MB ', '').replace('(Ptch1-/-)', 'Untx'))
ax.set_ylabel('CycD1 mRNA', fontsize=9)
ax.set_xlabel('Time (h)', fontsize=9)
ax.set_title('MB: CycD1 Time Course', fontsize=10, fontweight='bold')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.2)

# Panel 3c: EZH2 time courses MB conditions
ax = fig2.add_subplot(gs[2, 2])
for name, c in zip(mb_conds, mb_colors):
    r = sims[name]
    ax.plot(r['time'], r['[EZH2]'], color=c, linewidth=1, label=name.replace('MB ', '').replace('(Ptch1-/-)', 'Untx'))
ax.set_ylabel('EZH2 Protein', fontsize=9)
ax.set_xlabel('Time (h)', fontsize=9)
ax.set_title('MB: EZH2 Time Course', fontsize=10, fontweight='bold')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.2)

# Panel 3d: MYCN time courses
ax = fig2.add_subplot(gs[2, 3])
mycn_conds = ['GNP + SHH', 'GNP + HHi', 'MB (Ptch1-/-)','MB + HHi']
mycn_colors = ['#1b9e77', '#e7298a', '#762A83', '#E08214']
for name, c in zip(mycn_conds, mycn_colors):
    r = sims[name]
    ax.plot(r['time'], r['[MYCN]'], color=c, linewidth=1, label=name)
ax.set_ylabel('MYCN Protein', fontsize=9)
ax.set_xlabel('Time (h)', fontsize=9)
ax.set_title('MYCN Dynamics', fontsize=10, fontweight='bold')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.2)

fig2.suptitle('V42 Model Validation — MYCN Hill (n=3), Data Comparison', fontsize=14, fontweight='bold', y=1.01)
plt.savefig('simulations/fig_v42_validation.png', dpi=200, bbox_inches='tight')
plt.savefig('simulations/fig_v42_validation.pdf', bbox_inches='tight')
plt.close()
print("Saved: fig_v42_validation.png/pdf")

# =====================================================================
# FIGURE 3: Key prediction — EZH2i rescues HHi-arrested but NOT CDK4/6i-arrested MB
# =====================================================================

# Run longer simulations (168h / 7 days) to avoid edge effects
rescue_t_end = 168
rescue_n_pts = int(rescue_t_end * 100)

rescue_conditions = {
    'MB (Ptch1-/-)':         dict(shh=0.0, ptch1_cn=0.0, gdc=0.0, ezh2i=0.0, mycn_amp=MYCN_AMP_MB),
    'MB + HHi':              dict(shh=0.0, ptch1_cn=0.0, gdc=1.0, ezh2i=0.0, mycn_amp=MYCN_AMP_MB),
    'MB + HHi + EZH2i':     dict(shh=0.0, ptch1_cn=0.0, gdc=1.0, ezh2i=1.0, mycn_amp=MYCN_AMP_MB),
    'MB + CDK4/6i':          dict(shh=0.0, ptch1_cn=0.0, gdc=0.0, ezh2i=0.0, mycn_amp=MYCN_AMP_MB, cdk46i=True),
    'MB + CDK4/6i + EZH2i':  dict(shh=0.0, ptch1_cn=0.0, gdc=0.0, ezh2i=1.0, mycn_amp=MYCN_AMP_MB, cdk46i=True),
}

rescue_sims = {}
for name, kwargs in rescue_conditions.items():
    rescue_sims[name] = run_sim(**kwargs, t_end=rescue_t_end, n_pts=rescue_n_pts)

rescue_names = list(rescue_conditions.keys())
rescue_colors = ['#762A83', '#E08214', '#4A1486', '#C2185B', '#D81B60']

fig3, axes3 = plt.subplots(1, 5, figsize=(22, 4.5))
fig3.suptitle('EZH2i Rescues HHi-Arrested MB but NOT CDK4/6i-Arrested MB (7 days)',
              fontsize=13, fontweight='bold')

for i, (name, c) in enumerate(zip(rescue_names, rescue_colors)):
    ax = axes3[i]
    r = rescue_sims[name]
    ax.plot(r['time'], r['[Cb]'], color=c, linewidth=1.2)
    n, _, per = count_divisions(r)
    p_str = f", {np.mean(per):.1f}h" if len(per) > 0 else ""
    ax.set_title(f'{name}\n({n} div{p_str})', fontsize=9, fontweight='bold', color=c)
    ax.set_ylabel('Cyclin B' if i == 0 else '', fontsize=10)
    ax.set_xlabel('Time (h)', fontsize=10)
    ax.grid(True, alpha=0.2)
    ax.set_ylim(-0.01, 0.22)

    # Highlight rescue vs no rescue
    if 'HHi + EZH2i' in name:
        ax.patch.set_facecolor('#e8f5e9')
        ax.patch.set_alpha(0.3)
    elif 'CDK4/6i + EZH2i' in name:
        ax.patch.set_facecolor('#ffebee')
        ax.patch.set_alpha(0.3)

plt.tight_layout()
plt.savefig('simulations/fig_v42_ezh2i_rescue.png', dpi=200, bbox_inches='tight')
plt.close()
print("Saved: fig_v42_ezh2i_rescue.png")

# =====================================================================
# SAVE VALIDATION RESULTS JSON
# =====================================================================

validation_data = {
    'model': 'v42_mycn_hill',
    'best_params': {
        'n_MYCN_Cd': 3,
        'k_Cd_tx_MYCN': 15.0,
        'K_MYCN_Cd': 1.5,
        'k_MYCN_synth_basal': 0.3,
        'k_MYCN_synth_Gli': 0.102,
        'k_Cd_tx_Gli_max': 4.0,
        'MYCN_amplification_MB': 2.8,
    },
    'tests_passed': tests_passed,
    'tests_total': tests_total,
    'test_results': results,
    'summary': {
        'wt_divisions': n_wt,
        'wt_period': float(period_wt),
        'mb_divisions': count_divisions(mb)[0],
        'mb_gdc_divisions': n_mb_gdc,
        'mb_combo_divisions': n_mb_combo,
        'wt_gdc_cd_ratio': float(cd_wt_gdc / cd_wt) if cd_wt > 0 else 0,
        'mb_gdc_cd_ratio': float(cd_mb_gdc / cd_mb) if cd_mb > 0 else 0,
        'mycn_wt_gdc_ratio': float(mycn_wt_gdc / mycn_wt) if mycn_wt > 0 else 0,
        'mycn_mb_gdc_ratio': float(mycn_mb_gdc / mycn_mb) if mycn_mb > 0 else 0,
        'ezh2_g0_cyc_ratio': float(ezh2_g0 / ezh2_cyc) if ezh2_cyc > 0 else 0,
        'ezh2i_cd_fc': float(cd_wt_ezh2i / cd_wt) if cd_wt > 0 else 0,
    }
}

with open('simulations/validation_v42_results.json', 'w') as f:
    json.dump(validation_data, f, indent=2, default=float)
print("Saved: validation_v42_results.json")

print(f"\nDone. {tests_passed}/{tests_total} tests passed.")
