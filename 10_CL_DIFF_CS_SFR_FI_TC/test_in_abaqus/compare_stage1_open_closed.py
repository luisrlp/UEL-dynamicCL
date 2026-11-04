# -*- coding: utf-8 -*-
#
# Comparison of the stage-1 results (run_stage1.py, single element in simple shear, ramp + hold)
# for an open system (mu = INITMU prescribed at all nodes: bath) and a closed system (no-flux).
# Reads the stage1_output.pkl files written by getoutput_stage1.py in both study folders and
# pairs the runs by <sweep>/<case>.
#
# Usage (from test_in_abaqus/, plain Python with numpy and matplotlib):
#   python3 compare_stage1_open_closed.py
#   python3 compare_stage1_open_closed.py --open SA_stage1_open --closed SA_stage1_closed \
#           --out SA_stage1_compare --gammas 0.45 0.5 0.55
#
# Output in <out>/ (default SA_stage1_compare/):
#   compare_summary.csv        paired metrics of every case present in both folders
#   fig1_reservoir.pdf         theta_f/theta_f(0) and c/c(0) vs time (A_core, selected GAMMA):
#                              the free pool moves in the closed system, the content in the open one
#   fig2_cbtot.pdf             c_b,tot/c_b,tot(0) vs time (A_core, selected GAMMA)
#   fig3_exchange.pdf          net exchange with the bath, c(end)/c(0) - 1, vs GAMMA (A_core) and vs
#                              the catch fraction (B_catch_frac), open system
#   fig4_stress.pdf            sigma_12 vs time (absolute) and (sigma_open - sigma_closed)/sigma_peak,closed
#   fig5_directions.pdf        c_b,i(open)/c_b,i(closed) at the end of the hold vs lambda_i (A_core)
#   fig6_overview.pdf          open vs closed for all paired cases: change of c_b,tot, end stress
#                              relative to the closed sigma_peak, t_half (values below T_HALF_MIN
#                              are limited by the time increment and drawn as open markers)
# Time axis: t - t_ramp (the loading ramp is at negative times, the hold at positive times).
# Stresses are compared in absolute terms or relative to the closed sigma_peak (one common
# reference): the ratio to each run's own sigma_peak exaggerates the differences when the ramps
# differ (e.g. GAMMA = 0.55, where the open system unbinds more during loading).

import os
import sys
import csv
import glob
import pickle
import argparse
import numpy as np

from getoutput_stage1 import fortran_directions, read_global_factor, case_metrics

OUT_PKL = 'stage1_output.pkl'
LAWS = ('catch', 'slip')
T_HALF_MIN = 0.3        # s: shorter t_half values are limited by the time increment


def load_study(sa_dir):
    """{'<sweep>/<case>': output dict} from the stage1_output.pkl files of a study folder."""
    outs = {}
    for path in sorted(glob.glob(os.path.join(sa_dir, '*', '*', OUT_PKL))):
        with open(path, 'rb') as f:
            out = pickle.load(f)
        key = os.path.relpath(os.path.dirname(path), sa_dir).replace('\\', '/')
        outs[key] = out
    return outs


def time_axis(out):
    """t - t_ramp (end of the loading step)."""
    return out['t'] - out['step_ends'][0]


def hold_start(out):
    return int(np.nonzero(out['t'] >= out['step_ends'][0] - 1e-12)[0][0])


def a_core_key(law, gamma):
    return f'A_core/{law}_gamma{float(gamma):g}'


def paired_metrics(o, c, dirs):
    """Metrics of an open/closed pair (stresses relative to the closed sigma_peak)."""
    mo, mc = case_metrics(o, dirs), case_metrics(c, dirs)
    s_ref = mc['sigma_peak']
    row = {}
    for name, out, m in (('open', o, mo), ('closed', c, mc)):
        row[f'sigma_peak_{name}'] = m['sigma_peak']
        row[f'sigma_end_{name}'] = m['sigma_end']
        row[f'sigma_end_rel_{name}'] = m['sigma_end'] / s_ref if s_ref else np.nan
        row[f'R200_{name}'] = m.get('R200', np.nan)
        row[f't_half_{name}'] = m['t_half']
        row[f'dcb_tot_rel_{name}'] = m['dcb_tot_rel']
        row[f'thetaf_drift_{name}'] = m['thetaf_drift']
        row[f'dc_rel_{name}'] = out['c'][-1] / out['c'][0] - 1.0
    row['cb_amplification'] = (row['dcb_tot_rel_open'] / row['dcb_tot_rel_closed']
                               if abs(row['dcb_tot_rel_closed']) > 1e-12 else np.nan)
    row['dsigma_end_rel'] = row['sigma_end_rel_open'] - row['sigma_end_rel_closed']
    return row, mo, mc


def write_summary(pairs, out_dir):
    rows = []
    for key, (o, c, row, mo, mc) in sorted(pairs.items()):
        sweep, case = key.split('/')
        p = o['props']
        rows.append({'sweep': sweep, 'case': case,
                     **{k: p.get(k, np.nan) for k in ('GAMMA', 'KOFF0', 'KCATCH0', 'DX', 'DXC', 'ETAC')},
                     **row})
    path = os.path.join(out_dir, 'compare_summary.csv')
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows:
            writer.writerow({k: (f'{v:.6g}' if isinstance(v, float) else v) for k, v in r.items()})
    print(f"Summary written to {path} ({len(rows)} paired cases)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--open', default='SA_stage1_open', help='open-system study folder')
    parser.add_argument('--closed', default='SA_stage1_closed', help='closed-system study folder')
    parser.add_argument('--out', default='SA_stage1_compare', help='output folder')
    parser.add_argument('--gammas', nargs='+', type=float, default=[0.45, 0.5, 0.55],
                        help='GAMMA values of the time plots (A_core)')
    args = parser.parse_args()

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("matplotlib is needed: run with a Python that has numpy and matplotlib")

    opens, closeds = load_study(args.open), load_study(args.closed)
    if not opens or not closeds:
        sys.exit(f"No {OUT_PKL} files found in {args.open} and/or {args.closed} "
                 "(run getoutput_stage1.py on both studies first)")
    common = sorted(set(opens) & set(closeds))
    only = sorted(set(opens) ^ set(closeds))
    print(f"{len(common)} paired cases ({len(opens)} open, {len(closeds)} closed)")
    if only:
        print(f"Not paired (present in one folder only): {only}")
    os.makedirs(args.out, exist_ok=True)

    dirs = fortran_directions(read_global_factor())
    pairs = {}
    for key in common:
        o, c = opens[key], closeds[key]
        row, mo, mc = paired_metrics(o, c, dirs)
        pairs[key] = (o, c, row, mo, mc)
    write_summary(pairs, args.out)

    gammas = [g for g in args.gammas if all(a_core_key(l, g) in pairs for l in LAWS)]
    gcolor = {g: plt.cm.viridis(v) for g, v in zip(gammas, np.linspace(0.0, 0.85, max(len(gammas), 1)))}
    style = {'open': dict(ls='-'), 'closed': dict(ls='--')}

    def save(fig, name):
        fig.tight_layout()
        path = os.path.join(args.out, f'{name}.pdf')
        fig.savefig(path)
        plt.close(fig)
        print(f"Plot saved: {path}")

    def legend_keys(ax, loc):
        """Colour = GAMMA, line style = open/closed."""
        from matplotlib.lines import Line2D
        handles = [Line2D([], [], color=gcolor[g], label=f'GAMMA = {g:g}') for g in gammas]
        handles += [Line2D([], [], color='k', label=n, **style[n]) for n in ('open', 'closed')]
        ax.legend(handles=handles, fontsize=7, frameon=False, loc=loc)

    # Fig. 1: the reservoir (theta_f and total content)
    fig, axs = plt.subplots(2, 2, figsize=(10.0, 6.8), sharex=True)
    for i, law in enumerate(LAWS):
        for g in gammas:
            o, c = pairs[a_core_key(law, g)][:2]
            for name, out in (('open', o), ('closed', c)):
                t = time_axis(out)
                axs[i, 0].plot(t, out['thetaf'] / out['thetaf'][0], color=gcolor[g], **style[name])
                axs[i, 1].plot(t, out['c'] / out['c'][0], color=gcolor[g], **style[name])
        axs[i, 0].set_ylabel(rf'{law}: $\theta_f/\theta_f(0)$')
        axs[i, 1].set_ylabel(rf'{law}: $c/c(0)$')
    for ax in axs.flat:
        ax.grid(True, alpha=0.3)
    for ax in axs[1]:
        ax.set_xlabel('Time from the end of loading (s)')
    legend_keys(axs[0, 0], 'upper right')
    fig.suptitle('Free pool (closed: $\\theta_f$ moves) vs total content (open: $c$ moves)')
    save(fig, 'fig1_reservoir')

    # Fig. 2: total bound crosslinkers
    fig, axs = plt.subplots(1, 2, figsize=(10.0, 3.8), sharex=True)
    for ax, law in zip(axs, LAWS):
        for g in gammas:
            for name, out in zip(('open', 'closed'), pairs[a_core_key(law, g)][:2]):
                ax.plot(time_axis(out), out['cb_tot'] / out['cb_tot'][0], color=gcolor[g], **style[name])
        ax.set(title=law, xlabel='Time from the end of loading (s)', ylabel=r'$c_{b,tot}/c_{b,tot}(0)$')
        ax.grid(True, alpha=0.3)
    legend_keys(axs[0], 'lower right')
    fig.suptitle('Bound crosslinkers: open vs closed')
    save(fig, 'fig2_cbtot')

    # Fig. 3: net exchange with the bath (open system)
    fig, axs = plt.subplots(1, 2, figsize=(10.0, 3.8))
    for law, marker in zip(LAWS, ('o', 's')):
        pts = sorted((v[0]['props']['GAMMA'], v[2]['dc_rel_open'], v[2]['dc_rel_closed'])
                     for k, v in pairs.items() if k.startswith(f'A_core/{law}_gamma'))
        if pts:
            x, yo, yc = zip(*pts)
            axs[0].plot(x, 100 * np.array(yo), marker=marker, label=f'{law}, open')
            axs[0].plot(x, 100 * np.array(yc), marker=marker, ls='--', mfc='none', label=f'{law}, closed')
    axs[0].axhline(0.0, color='k', lw=0.6)
    axs[0].set(xlabel='GAMMA', ylabel=r'$c(end)/c(0) - 1$ (%)', title='A_core: uptake (> 0) / release (< 0)')
    axs[0].legend(fontsize=7, frameon=False, loc='lower left')
    b = {}
    for k, v in pairs.items():
        if k.startswith('B_catch_frac/'):
            p = v[0]['props']
            frac = p['KCATCH0'] / (p['KOFF0'] + p['KCATCH0'])
            b.setdefault(p['GAMMA'], []).append((frac, v[2]['dc_rel_open']))
    bcolor = {g: plt.cm.plasma(v) for g, v in zip(sorted(b), np.linspace(0.0, 0.85, max(len(b), 1)))}
    for g in sorted(b):
        x, y = zip(*sorted(b[g]))
        axs[1].plot(x, 100 * np.array(y), marker='o', color=bcolor[g], label=f'GAMMA = {g:g}')
    axs[1].axhline(0.0, color='k', lw=0.6)
    axs[1].set(xlabel=r'Catch fraction $k_{c0}/k_{off}(0)$', ylabel=r'$c(end)/c(0) - 1$ (%)',
               title='B_catch_frac: open system')
    if b:
        axs[1].legend(fontsize=7, frameon=False, loc='upper left')
    for ax in axs:
        ax.grid(True, alpha=0.3)
    fig.suptitle('Net crosslinker exchange with the bath at the end of the hold')
    save(fig, 'fig3_exchange')

    # Fig. 4: stress, absolute and difference relative to the closed sigma_peak
    fig, axs = plt.subplots(2, 2, figsize=(10.0, 6.8), sharex=True)
    for j, law in enumerate(LAWS):
        for g in gammas:
            o, c = pairs[a_core_key(law, g)][:2]
            ref = c['sig'][hold_start(c), 3]
            for name, out in (('open', o), ('closed', c)):      # hold only (the ramp starts at 0)
                i0 = hold_start(out)
                axs[0, j].plot(time_axis(out)[i0:], out['sig'][i0:, 3], color=gcolor[g], **style[name])
            i0 = hold_start(c)
            t = time_axis(c)[i0:]
            s_open = np.interp(t, time_axis(o)[hold_start(o):], o['sig'][hold_start(o):, 3])
            axs[1, j].plot(t, 100 * (s_open - c['sig'][i0:, 3]) / ref, color=gcolor[g], label=f'GAMMA = {g:g}')
        axs[0, j].set(title=law, ylabel=r'$\sigma_{12}$ (Pa)', yscale='log')
        axs[1, j].set(xlabel='Hold time (s)',
                      ylabel=r'$(\sigma_{open} - \sigma_{closed})/\sigma_{peak,closed}$ (%)')
        axs[1, j].axhline(0.0, color='k', lw=0.6)
    for ax in axs.flat:
        ax.grid(True, alpha=0.3)
    legend_keys(axs[0, 0], 'center right')
    fig.suptitle('Shear stress: open vs closed (common reference: closed $\\sigma_{12}$ at the end of loading)')
    save(fig, 'fig4_stress')

    # Fig. 5: per-direction bound crosslinkers, open/closed, at the end of the hold
    fig, axs = plt.subplots(1, 2, figsize=(10.0, 3.8), sharey=True)
    a_gammas = sorted({v[0]['props']['GAMMA'] for k, v in pairs.items() if k.startswith('A_core/')})
    acolor = {g: plt.cm.viridis(v) for g, v in zip(a_gammas, np.linspace(0.0, 0.9, max(len(a_gammas), 1)))}
    for ax, law in zip(axs, LAWS):
        for g in a_gammas:
            key = a_core_key(law, g)
            if key not in pairs:
                continue
            o, c, row, mo, mc = pairs[key]
            nd = min(len(dirs), o['cb_dirs'].shape[1])
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = o['cb_dirs'][-1, :nd] / c['cb_dirs'][-1, :nd]
            ok = np.isfinite(ratio) & (c['cb_dirs'][-1, :nd] > 1e-6 * c['cb_dirs'][0, :nd])
            ax.plot(mc['lambda_i'][:nd][ok], ratio[ok], 'o', ms=3, color=acolor[g], label=f'GAMMA = {g:g}')
        ax.axhline(1.0, color='k', lw=0.6)
        ax.set(title=law, xlabel=r'$\lambda_i$')
        ax.grid(True, alpha=0.3)
    axs[0].set_ylabel(r'$c_{b,i}(open)/c_{b,i}(closed)$ at the end of the hold')
    axs[0].legend(fontsize=7, frameon=False, loc='best')
    fig.suptitle('Per-direction bound crosslinkers: open vs closed (depleted directions omitted)')
    save(fig, 'fig5_directions')

    # Fig. 6: overview of all paired cases
    sweeps = sorted({k.split('/')[0] for k in pairs})
    scolor = {s: plt.cm.tab10(i % 10) for i, s in enumerate(sweeps)}
    fig, axs = plt.subplots(1, 3, figsize=(13.0, 4.2))
    panels = [('dcb_tot_rel', 100.0, r'$c_{b,tot}$ change (%)'),
              ('sigma_end_rel', 1.0, r'$\sigma_{12,end}/\sigma_{peak,closed}$'),
              ('t_half', 1.0, r'$t_{1/2}$ (s)')]
    for ax, (name, scale, label) in zip(axs, panels):
        allv = []
        for s in sweeps:
            for key, v in pairs.items():
                if not key.startswith(s + '/'):
                    continue
                row, law = v[2], key.split('/')[1].split('_')[0]
                xo, xc = row[f'{name}_closed'] * scale, row[f'{name}_open'] * scale
                if not (np.isfinite(xo) and np.isfinite(xc)):
                    continue
                limited = name == 't_half' and min(xo, xc) < T_HALF_MIN
                ax.plot(xo, xc, 'o' if law != 'slip' else 's', ms=4, color=scolor[s],
                        mfc='none' if limited else scolor[s])
                allv += [xo, xc]
        if allv:
            lo, hi = min(allv), max(allv)
            ax.plot([lo, hi], [lo, hi], 'k-', lw=0.6)
        ax.set(xlabel=f'closed: {label}', ylabel=f'open: {label}')
        if name == 't_half':
            ax.set(xscale='log', yscale='log')
        ax.grid(True, alpha=0.3)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker='o', ls='none', color=scolor[s], label=s) for s in sweeps]
    handles += [Line2D([], [], marker='s', ls='none', color='gray', label='slip'),
                Line2D([], [], marker='o', ls='none', color='gray', label='catch / other'),
                Line2D([], [], marker='o', ls='none', color='gray', mfc='none',
                       label=f't_half < {T_HALF_MIN:g} s (time-step limited)')]
    axs[0].legend(handles=handles, fontsize=6, frameon=False, loc='upper left')
    fig.suptitle('All paired cases: open vs closed')
    save(fig, 'fig6_overview')

    # Short report of the main differences
    print("\nA_core: change of c_b,tot (open / closed, amplification) and end stress relative to the closed sigma_peak")
    for law in LAWS:
        for g in a_gammas:
            key = a_core_key(law, g)
            if key in pairs:
                r = pairs[key][2]
                print(f"  {law:5s} GAMMA = {g:<5g} dcb_tot {100 * r['dcb_tot_rel_open']:+6.2f}% / {100 * r['dcb_tot_rel_closed']:+6.2f}%"
                      f" (x{r['cb_amplification']:.2f})  sigma_end/sigma_peak,closed {r['sigma_end_rel_open']:.3f} / "
                      f"{r['sigma_end_rel_closed']:.3f}  exchange {100 * r['dc_rel_open']:+.2f}%")


if __name__ == '__main__':
    main()
