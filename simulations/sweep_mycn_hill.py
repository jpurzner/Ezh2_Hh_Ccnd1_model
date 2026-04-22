"""
Focused MYCN Hill sweep for v42 model.

The previous sweep showed that simple Michaelis-Menten (n=1) can't differentiate
WT (MYCN~0.3) from MB (MYCN~0.84) enough to match data:
  - WT+GDC CycD1: 14% (86% reduction) — needs MYCN to be negligible in WT
  - MB+GDC CycD1: 40% (60% reduction) — needs MYCN to be significant in MB

Solution: Hill coefficient n>=2 for MYCN->CycD1 creates threshold behavior.

Sweeps: n_MYCN_Cd, k_Cd_tx_MYCN, K_MYCN_Cd, MYCN_amplification_MB, k_Cd_tx_Gli_max
"""
import sys
sys.path.insert(0, '/Users/jpurzner/Dropbox/Q_research/py_projects/ezh2_cyclind1_sym')
import tellurium as te
import numpy as np
from scipy.signal import find_peaks
from src.build_model_v42_mycn import build_model_v42
import time
import json

base_model = build_model_v42()


def run_condition(model_str, shh, cn, gdc, ezh2i, mycn_amp, t_end=600, n_pts=6000):
    """Run a single condition and extract metrics."""
    try:
        rr = te.loada(model_str)
        rr['SHH'] = shh
        rr['Ptch1_copy_number'] = cn
        rr['GDC0449'] = gdc
        rr['EZH2i'] = ezh2i
        rr['MYCN_amplification'] = mycn_amp
        if cn == 0.0:
            rr['Ptch1_mRNA'] = 0.0
            rr['Ptch1_free'] = 0.0
            rr['SHH_Ptch'] = 0.0
            rr['Smo_active'] = 1.0
        elif shh == 0.0:
            rr['SHH_Ptch'] = 0.0
            rr['Ptch1_free'] = 1.0
            rr['Smo_active'] = 0.01
        result = rr.simulate(0, t_end, n_pts)
    except Exception:
        return None

    idx = -1200
    t = result['time']
    cb = result['[Cb]']
    mask = t >= 50
    peaks, _ = find_peaks(cb[mask], prominence=0.05, distance=10)
    peak_times = t[mask][peaks]
    periods = np.diff(peak_times) if len(peak_times) > 1 else np.array([])

    return {
        'n_div': len(peaks),
        'period': float(np.mean(periods)) if len(periods) > 0 else 0.0,
        'EZH2': float(np.mean(result['[EZH2]'][idx:])),
        'Cd_mRNA': float(np.mean(result['[Cd_mRNA]'][idx:])),
        'Gli1': float(np.mean(result['[Gli1]'][idx:])),
        'MYCN': float(np.mean(result['[MYCN]'][idx:])),
        'E2F': float(np.mean(result['[E2F]'][idx:])),
        'Cb_mean': float(np.mean(result['[Cb]'][idx:])),
    }


def test_params(n_MYCN_Cd, k_Cd_tx_MYCN, K_MYCN_Cd, k_Cd_tx_Gli_max,
                k_MYCN_synth_basal, MYCN_amp_MB):
    """Test a parameter set across all conditions."""
    model_str = base_model

    # Replace parameters via string substitution
    model_str = model_str.replace(
        'n_MYCN_Cd = 2',
        f'n_MYCN_Cd = {n_MYCN_Cd}')
    model_str = model_str.replace(
        'k_Cd_tx_Gli_max = 4.0',
        f'k_Cd_tx_Gli_max = {k_Cd_tx_Gli_max}')
    model_str = model_str.replace(
        'k_Cd_tx_MYCN = 1.5',
        f'k_Cd_tx_MYCN = {k_Cd_tx_MYCN}')
    model_str = model_str.replace(
        'K_MYCN_Cd = 1.0',
        f'K_MYCN_Cd = {K_MYCN_Cd}')
    model_str = model_str.replace(
        'k_MYCN_synth_basal = 0.5',
        f'k_MYCN_synth_basal = {k_MYCN_synth_basal}')

    # Derive k_MYCN_synth_Gli for ~22% Gli-dependence in WT
    k_gli = 0.34 * k_MYCN_synth_basal
    model_str = model_str.replace(
        'k_MYCN_synth_Gli = 0.15',
        f'k_MYCN_synth_Gli = {k_gli:.4f}')

    # Conditions
    conditions = [
        ('WT_SHH',     0.5, 1.0, 0.0, 0.0, 1.0),
        ('WT_noSHH',   0.0, 1.0, 0.0, 0.0, 1.0),
        ('WT_GDC',     0.5, 1.0, 1.0, 0.0, 1.0),
        ('MB',         0.0, 0.0, 0.0, 0.0, MYCN_amp_MB),
        ('MB_GDC',     0.0, 0.0, 1.0, 0.0, MYCN_amp_MB),
        ('MB_EZH2i',   0.0, 0.0, 0.0, 1.0, MYCN_amp_MB),
        ('MB_combo',   0.0, 0.0, 1.0, 1.0, MYCN_amp_MB),
    ]

    results = {}
    for name, shh, cn, gdc, ezh2i, mycn_amp in conditions:
        r = run_condition(model_str, shh, cn, gdc, ezh2i, mycn_amp)
        if r is None:
            return None
        results[name] = r
        if name == 'WT_SHH' and r['n_div'] < 10:
            return None
        if name == 'WT_noSHH' and r['n_div'] > 2:
            return None

    return results


def score_results(results):
    """Score how well results match data. Lower = better."""
    wt = results['WT_SHH']
    g0 = results['WT_noSHH']
    wt_gdc = results['WT_GDC']
    mb = results['MB']
    mb_gdc = results['MB_GDC']
    mb_ezh2i = results['MB_EZH2i']
    mb_combo = results['MB_combo']

    # Hard criteria
    if wt['n_div'] < 15:
        return None, "WT too few div"
    if g0['n_div'] > 0:
        return None, "G0 not arrested"
    if mb['n_div'] < 10:
        return None, "MB too few div"

    # CycD1 ratios
    if wt['Cd_mRNA'] <= 0 or mb['Cd_mRNA'] <= 0:
        return None, "CycD1 zero"
    wt_gdc_ratio = wt_gdc['Cd_mRNA'] / wt['Cd_mRNA']
    mb_gdc_ratio = mb_gdc['Cd_mRNA'] / mb['Cd_mRNA']

    wt_gdc_err = abs(wt_gdc_ratio - 0.14) / 0.14
    mb_gdc_err = abs(mb_gdc_ratio - 0.40) / 0.40

    # MYCN GDC resistance
    if wt['MYCN'] <= 0 or mb['MYCN'] <= 0:
        return None, "MYCN zero"
    mycn_wt_ratio = wt_gdc['MYCN'] / wt['MYCN']
    mycn_mb_ratio = mb_gdc['MYCN'] / mb['MYCN']
    mycn_wt_err = abs(mycn_wt_ratio - 0.78) / 0.78
    mycn_mb_err = abs(mycn_mb_ratio - 0.86) / 0.86

    # EZH2 regulation
    if wt['EZH2'] <= 0:
        return None, "EZH2 zero"
    ezh2_ratio = g0['EZH2'] / wt['EZH2']
    ezh2_err = 0.0
    if ezh2_ratio < 0.3:
        ezh2_err = (0.3 - ezh2_ratio) / 0.3
    elif ezh2_ratio > 0.7:
        ezh2_err = (ezh2_ratio - 0.7) / 0.7

    # EZH2 GDC reduction in MB
    mb_ezh2_red = (1 - mb_gdc['EZH2'] / mb['EZH2']) * 100 if mb['EZH2'] > 0 else 0
    ezh2_gdc_err = 0.0
    if mb_ezh2_red < 15:
        ezh2_gdc_err = (15 - mb_ezh2_red) / 15
    elif mb_ezh2_red > 55:
        ezh2_gdc_err = (mb_ezh2_red - 55) / 55

    # EZH2i CycD1 fold change
    mb_fc = mb_ezh2i['Cd_mRNA'] / mb['Cd_mRNA'] if mb['Cd_mRNA'] > 0 else 0
    fc_err = 0.0
    if mb_fc < 1.3:
        fc_err = (1.3 - mb_fc) / 1.3
    elif mb_fc > 3.0:
        fc_err = (mb_fc - 3.0) / 3.0

    # Period
    period_err = abs(wt['period'] - 23.0) / 23.0 if wt['period'] > 0 else 1.0

    # Composite score (strongly prioritize MB+GDC CycD1 ratio)
    score = (
        4.0 * mb_gdc_err +
        2.0 * wt_gdc_err +
        1.5 * mycn_wt_err +
        1.0 * mycn_mb_err +
        1.5 * ezh2_err +
        1.0 * ezh2_gdc_err +
        1.0 * fc_err +
        0.5 * period_err
    )

    details = {
        'wt_gdc_cd_ratio': wt_gdc_ratio,
        'mb_gdc_cd_ratio': mb_gdc_ratio,
        'mycn_wt_gdc_ratio': mycn_wt_ratio,
        'mycn_mb_gdc_ratio': mycn_mb_ratio,
        'ezh2_g0_ratio': ezh2_ratio,
        'mb_ezh2_gdc_red': mb_ezh2_red,
        'mb_ezh2i_fc': mb_fc,
        'wt_period': wt['period'],
        'wt_div': wt['n_div'],
        'mb_div': mb['n_div'],
        'mb_gdc_div': mb_gdc['n_div'],
        'mb_combo_div': mb_combo['n_div'],
        'mb_gdc_cd': mb_gdc['Cd_mRNA'],
        'mb_combo_cd': mb_combo['Cd_mRNA'],
    }

    return score, details


# =====================================================================
# PARAMETER SWEEP
# =====================================================================
print("=" * 130)
print("MYCN HILL SWEEP (v42 with cooperative MYCN->CycD1)")
print("=" * 130)
print(f"\nTargets:")
print(f"  WT+GDC CycD1: 14% of WT (86% reduction)")
print(f"  MB+GDC CycD1: 40% of MB (60% reduction)")
print(f"  MYCN WT+GDC/WT: 0.78, MB+GDC/MB: 0.86")
print(f"  EZH2 G0/cycling: 0.5-0.7, EZH2i FC: ~2x")
print()

# Sweep parameters
n_MYCN_Cd_vals = [2, 3]
k_Cd_tx_MYCN_vals = [1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0]
K_MYCN_Cd_vals = [0.3, 0.5, 0.7, 1.0, 1.5]
k_Cd_tx_Gli_max_vals = [2.5, 3.0, 3.5, 4.0]
MYCN_amp_MB_vals = [2.8, 4.0, 5.0]
# Fix k_MYCN_synth_basal at 0.3 (best from prev sweep)
k_MYCN_synth_basal = 0.3

total = (len(n_MYCN_Cd_vals) * len(k_Cd_tx_MYCN_vals) * len(K_MYCN_Cd_vals) *
         len(k_Cd_tx_Gli_max_vals) * len(MYCN_amp_MB_vals))
print(f"Total combinations: {total}")
print(f"Fixed: k_MYCN_synth_basal = {k_MYCN_synth_basal}")
print()

# Header
hdr = (f"{'n':>2} {'kMYCN':>6} {'K_Cd':>6} {'kGli':>6} {'Amp':>5} | "
       f"{'Score':>6} | "
       f"{'WTgdc%':>7} {'MBgdc%':>7} {'MYCNwt':>7} {'MYCNmb':>7} | "
       f"{'EZH2r':>6} {'EZH2g':>6} {'FC':>5} | "
       f"{'nWT':>4} {'nMB':>4} {'nMBg':>5} {'nComb':>5} {'Per':>5} | "
       f"{'MBg_Cd':>7} {'Cmb_Cd':>7}")
print(hdr)
print("-" * len(hdr))

good = []
n_tested = 0
n_passed = 0
n_skipped = 0
t_start = time.time()

for n_MYCN_Cd in n_MYCN_Cd_vals:
    for MYCN_amp_MB in MYCN_amp_MB_vals:
        for k_Cd_tx_Gli_max in k_Cd_tx_Gli_max_vals:
            for K_MYCN_Cd in K_MYCN_Cd_vals:
                for k_Cd_tx_MYCN in k_Cd_tx_MYCN_vals:
                    n_tested += 1

                    results = test_params(
                        n_MYCN_Cd, k_Cd_tx_MYCN, K_MYCN_Cd,
                        k_Cd_tx_Gli_max, k_MYCN_synth_basal, MYCN_amp_MB)
                    if results is None:
                        n_skipped += 1
                        continue

                    score, details = score_results(results)
                    if score is None:
                        n_skipped += 1
                        continue

                    n_passed += 1
                    params = {
                        'n_MYCN_Cd': n_MYCN_Cd,
                        'k_Cd_tx_MYCN': k_Cd_tx_MYCN,
                        'K_MYCN_Cd': K_MYCN_Cd,
                        'k_Cd_tx_Gli_max': k_Cd_tx_Gli_max,
                        'k_MYCN_synth_basal': k_MYCN_synth_basal,
                        'MYCN_amplification_MB': MYCN_amp_MB,
                    }
                    good.append((score, params, details, results))

                    d = details
                    wt_gdc_pct = (1 - d['wt_gdc_cd_ratio']) * 100
                    mb_gdc_pct = (1 - d['mb_gdc_cd_ratio']) * 100
                    print(f"{n_MYCN_Cd:2d} {k_Cd_tx_MYCN:6.1f} {K_MYCN_Cd:6.2f} "
                          f"{k_Cd_tx_Gli_max:6.2f} {MYCN_amp_MB:5.1f} | "
                          f"{score:6.2f} | "
                          f"{wt_gdc_pct:6.1f}% {mb_gdc_pct:6.1f}% "
                          f"{d['mycn_wt_gdc_ratio']:7.3f} {d['mycn_mb_gdc_ratio']:7.3f} | "
                          f"{d['ezh2_g0_ratio']:6.3f} {d['mb_ezh2_gdc_red']:5.1f}% "
                          f"{d['mb_ezh2i_fc']:5.2f} | "
                          f"{d['wt_div']:4d} {d['mb_div']:4d} {d['mb_gdc_div']:5d} "
                          f"{d['mb_combo_div']:5d} {d['wt_period']:5.1f} | "
                          f"{d['mb_gdc_cd']:7.3f} {d['mb_combo_cd']:7.3f}")

                    if n_tested % 100 == 0:
                        elapsed = time.time() - t_start
                        rate = n_tested / elapsed
                        remaining = (total - n_tested) / rate / 60
                        print(f"  ... [{n_tested}/{total}] {elapsed/60:.1f}min, "
                              f"~{remaining:.1f}min left, {n_passed} passed")

elapsed = time.time() - t_start
print(f"\n{'='*130}")
print(f"SWEEP COMPLETE: {n_tested} tested, {n_passed} passed, {n_skipped} skipped")
print(f"Time: {elapsed/60:.1f} minutes")

if not good:
    print("\nNO PASSING COMBINATIONS FOUND.")
    sys.exit(1)

good.sort(key=lambda x: x[0])

print(f"\n{'='*130}")
print(f"TOP 20 RESULTS (sorted by score, lower=better)")
print(f"{'='*130}")

for i, (score, params, details, results) in enumerate(good[:20]):
    d = details
    p = params
    wt_gdc_pct = (1 - d['wt_gdc_cd_ratio']) * 100
    mb_gdc_pct = (1 - d['mb_gdc_cd_ratio']) * 100
    print(f"\n--- Rank {i+1} (score={score:.3f}) ---")
    print(f"  n_MYCN={p['n_MYCN_Cd']}, k_Cd_tx_MYCN={p['k_Cd_tx_MYCN']:.1f}, "
          f"K_MYCN_Cd={p['K_MYCN_Cd']:.2f}, k_Cd_tx_Gli_max={p['k_Cd_tx_Gli_max']:.1f}, "
          f"MYCN_amp={p['MYCN_amplification_MB']:.1f}")
    print(f"  WT+GDC CycD1: {wt_gdc_pct:.1f}% reduction (target: 86%)")
    print(f"  MB+GDC CycD1: {mb_gdc_pct:.1f}% reduction (target: 60%)")
    print(f"  MYCN WT ratio: {d['mycn_wt_gdc_ratio']:.3f} (target: 0.78)")
    print(f"  MYCN MB ratio: {d['mycn_mb_gdc_ratio']:.3f} (target: 0.86)")
    print(f"  EZH2 G0/cyc: {d['ezh2_g0_ratio']:.3f}, EZH2 GDC red: {d['mb_ezh2_gdc_red']:.1f}%")
    print(f"  EZH2i FC: {d['mb_ezh2i_fc']:.2f}x, Period: {d['wt_period']:.1f}h")
    print(f"  Divisions: WT={d['wt_div']}, MB={d['mb_div']}, MB+GDC={d['mb_gdc_div']}, "
          f"MB+combo={d['mb_combo_div']}")

# Detailed species table for top 5
print(f"\n{'='*130}")
print("DETAILED TOP 5")
print(f"{'='*130}")

for i, (score, params, details, results) in enumerate(good[:5]):
    p = params
    print(f"\n--- Rank {i+1} (score={score:.3f}) ---")
    print(f"  {'Condition':20s} | {'Div':>4} {'Per':>6} | {'CycD1':>7} {'EZH2':>6} "
          f"{'Gli1':>6} {'MYCN':>6} {'E2F':>5}")
    print(f"  {'-'*80}")
    for cname in ['WT_SHH', 'WT_noSHH', 'WT_GDC', 'MB', 'MB_GDC', 'MB_EZH2i', 'MB_combo']:
        r = results[cname]
        per_str = f"{r['period']:.1f}h" if r['period'] > 0 else "N/A"
        print(f"  {cname:20s} | {r['n_div']:4d} {per_str:>6s} | "
              f"{r['Cd_mRNA']:7.3f} {r['EZH2']:6.3f} {r['Gli1']:6.3f} "
              f"{r['MYCN']:6.3f} {r['E2F']:5.1f}")

# Save results
output = {
    'sweep_params': {
        'n_MYCN_Cd_vals': n_MYCN_Cd_vals,
        'k_Cd_tx_MYCN_vals': k_Cd_tx_MYCN_vals,
        'K_MYCN_Cd_vals': K_MYCN_Cd_vals,
        'k_Cd_tx_Gli_max_vals': k_Cd_tx_Gli_max_vals,
        'MYCN_amp_MB_vals': MYCN_amp_MB_vals,
        'k_MYCN_synth_basal': k_MYCN_synth_basal,
    },
    'stats': {
        'total': total,
        'tested': n_tested,
        'passed': n_passed,
        'skipped': n_skipped,
        'time_minutes': elapsed / 60,
    },
    'top_20': [],
}
for score, params, details, results in good[:20]:
    entry = {
        'score': score,
        'params': params,
        'details': details,
    }
    output['top_20'].append(entry)

with open('simulations/sweep_mycn_hill_results.json', 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nResults saved to simulations/sweep_mycn_hill_results.json")
