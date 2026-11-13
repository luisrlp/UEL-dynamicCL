"""
LAOS loops (Unterberger et al. 2013, Fig. 4a) of saved optimisation runs against the data, with the Fortran engine
(closed-system driver laos_closed, same path as optimise_laos_v2.py --engine fortran).

  python -B plot_laos_runs.py --runs final_F1p                       # last cycle (the fitted one)
  python -B plot_laos_runs.py --runs final_F1p --cycle 1 3 7         # chosen cycles
  python -B plot_laos_runs.py --runs final_F1p --cycle all           # all cycles
  python -B plot_laos_runs.py --runs final_F1p final_Pdxc_0 --r2     # several runs, R^2 in the legend
  python -B plot_laos_runs.py --runs final_F1p --set dxc_nm=0.3      # variant of a saved parameter set

Parameters: optimum and fixed values of optimisation_runs/<run>/result.json; dt, factor and gamma_A from its cfg.
The fitted number of cycles is simulated by default (--cycles N to change; --dt and --factor to change the time step
and the number of directions); with the fitted settings, the last cycle is checked against the run's best_loop.npz. R^2 = 1 - SS_res/SS_tot over the data points of
both branches, with the model interpolated at the data strains as in the fit (unit weights).
Output: figures/validation/laos_<runs>.pdf (or --out).
"""
import argparse, json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import optimise_laos_v2 as opt


def load_run(name):
    """Parameter dict (fixed values and optimum) and fit settings of a saved run."""
    r = json.load(open(os.path.join(opt.REPO, 'optimisation_runs', name, 'result.json')))
    return dict(r['fixed'], **r['optimum']), dict(r['cfg'])


def r2(gm, sm, data):
    """Coefficient of determination of the loop (gm, sm) against the data points of both branches."""
    res = opt.residuals_from_loop(gm, sm, {b: (g, s, np.ones_like(g)) for b, (g, s, w) in data.items()})
    s = np.concatenate([data[b][1] for b in (1, -1)])
    return 1.0 - np.sum(res**2) / np.sum((s - s.mean())**2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--runs', nargs='+', required=True, help='run folders in optimisation_runs/')
    ap.add_argument('--cycle', nargs='+', default=['last'], help="cycles to plot: numbers, 'last' or 'all'")
    ap.add_argument('--cycles', type=int, default=None, help='number of simulated cycles (default: as in the fit)')
    ap.add_argument('--dt', type=float, default=None, help='time step [s] (default: as in the fit)')
    ap.add_argument('--factor', type=int, default=None, help='integration factor, NDIR = 10 factor^2 (default: as in the fit)')
    ap.add_argument('--set', nargs='+', default=[], metavar='NAME=VALUE', help='override parameters of all runs')
    ap.add_argument('--r2', action='store_true', help='show R^2 in the legend')
    ap.add_argument('--out', default=None, help='output file (default figures/validation/laos_<runs>.pdf)')
    a = ap.parse_args()

    data, _ = opt.load_data()
    overrides = {k: float(v) for k, v in (s.split('=') for s in a.set)}
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for b, lab in ((1, 'Unterberger et al. 2013'), (-1, None)):
        ax.plot(data[b][0], data[b][1], 'ko', ms=3, mfc='none', label=lab)

    for j, name in enumerate(a.runs):
        p, cfg = load_run(name)
        p.update(overrides)
        cfg['fortran_engine'] = 1
        n_fit, cfg['cycles'] = cfg['cycles'], a.cycles or cfg['cycles']
        dt_fit, cfg['dt'] = cfg['dt'], a.dt or cfg['dt']
        factor_fit, cfg['factor'] = cfg['factor'], a.factor or cfg['factor']
        opt.check_fortran(cfg['factor'])
        t, g, s = opt.simulate(p, cfg)
        cycles = opt.split_cycles(t, g, s, cfg['cycles'])

        # with the fitted number of cycles, the last cycle of an unmodified run must reproduce the saved fitted loop
        best = os.path.join(opt.REPO, 'optimisation_runs', name, 'best_loop.npz')
        if not overrides and cfg['cycles'] == n_fit and cfg['dt'] == dt_fit and cfg['factor'] == factor_fit and os.path.exists(best):
            print(f'{name}: max |sigma - best_loop| = {np.abs(cycles[-1][1] - np.load(best)["sigma"]).max():.2e} Pa')

        n = cfg['cycles']
        ks = (list(range(1, n + 1)) if 'all' in a.cycle else
              [n if c == 'last' else int(c) for c in a.cycle])
        if max(ks) > n:
            raise SystemExit(f'{name}: cycle {max(ks)} requested but only {n} simulated (use --cycles {max(ks)})')
        for i, k in enumerate(ks):
            gm, sm = cycles[k - 1]
            if len(a.runs) == 1:                         # one run: colour per cycle
                color, alpha = plt.cm.viridis(i / max(len(ks) - 1, 1)), 1.0
            else:                                        # several runs: colour per run, later cycles darker
                color, alpha = f'C{j}', 0.35 + 0.65 * (i + 1) / len(ks)
            lab = f'{name}, cycle {k}' + (f', $R^2$ = {r2(gm, sm, data):.3f}' if a.r2 else '')
            ax.plot(gm, sm, '-', color=color, alpha=alpha, lw=1.5, label=lab)
            print(f'{name} cycle {k}: R^2 = {r2(gm, sm, data):.4f}')

    ax.set_xlabel(r'Shear strain $\gamma$')
    ax.set_ylabel(r'Shear stress $\sigma_{12}$ (Pa)')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    out = a.out or os.path.join(opt.REPO, 'figures', 'validation',
                                'laos_' + '_'.join(a.runs) + ('_set' if overrides else '') + (f'_dt{a.dt:g}' if a.dt else '') + (f'_f{a.factor}' if a.factor else '') + '.pdf')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, bbox_inches='tight')
    print(f'saved {out}')


if __name__ == '__main__':
    main()
