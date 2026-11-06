"""
Amplitude sweep (Unterberger et al. 2013, Fig. 4b) for the optimum of saved LAOS fits in optimisation_runs/.

One LAOS run per amplitude from rest (fresh material), 5 cycles by default (--cycles N), G' and G'' of the last cycle
from main.moduli (Eq. (1) of the paper, positive peaks), as in main.py. Runs fitted with the consistent (C4) kinetics
use the closed-system Fortran driver (optimise_laos_v2.simulate, as plot_laos_runs.py), with the parameters of
result.json (fixed values and optimum). Older runs use the Python material model, with the constants of setup.json
(or of optimise_laos.py v1).

  python -B amplitude_sweep_runs.py --runs final_F1p                     # one run
  python -B amplitude_sweep_runs.py --runs final_F1p final_Pdxc_0        # several runs in one figure
  python -B amplitude_sweep_runs.py --runs final_F1p --set dxc_nm=0.3    # variant of a saved set (C4 runs only)
  python -B amplitude_sweep_runs.py --runs final_F1p --max-gamma 0.5     # all data amplitudes (default up to 0.2)
  python -B amplitude_sweep_runs.py --runs final_F1p --dt 0.025          # finer time step (default: as in the fit)
  python -B amplitude_sweep_runs.py --runs final_F1p --factor 6          # more directions (default: as in the fit)
  python -B amplitude_sweep_runs.py --verify --runs laos0                # reproduce best_loop.npz (mapping check)

Outputs: validation/0_amplitude_sweep_<run>.csv per run, figures/validation/0_amplitude_sweep_<runs>.pdf (or --out)
"""
import argparse, contextlib, io, json, os
from concurrent.futures import ProcessPoolExecutor
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, 'optimisation_runs')
FREQ = 0.05
MAX_GAMMA = 0.2                    # default largest amplitude (--max-gamma)
# Constants of optimise_laos.py (v1) for runs saved before setup.json existed (laos0, laos1)
V1_CONSTANTS = {'FIXED_DIFF': {'CHI': 0.1, 'THETA': 294.0, 'D': 0.1, 'MU0': 0.0, 'VMOL': 0.15,
                               'RGAS': 8.31446261815324, 'CACTIN': 0.0095, 'R': 0.1, 'RBMAX': 0.25, 'RFMAX': 1.0},
                'FIXED_MECH': {'filament': 'ext', 'KBULK': 1000.0, 'MU0STRETCH': 38600.0, 'Lp': 16.0, 'BETA': 0.438,
                               'R0C': 0.014, 'ETA': 1.0, 'b_disp': 1e-06, 'prefdir': [1.0, 0.0, 0.0],
                               'stretch_free_binding': True, 'xi_t': 1.0},
                'KEQ': 5.0, 'THETA': 294.0}
# laos0 was run with ETA = 0.989 (checked with --verify: reproduces best_loop.npz exactly; laos1 used ETA = 1)
V1_MECH_OVERRIDES = {'laos0': {'ETA': 0.989}}


def load_run(name, eta=None):
    """Material properties (diff, mech) and fit settings (cfg) of the optimum of one saved run."""
    res = json.load(open(os.path.join(RUNS, name, 'result.json')))
    setup_file = os.path.join(RUNS, name, 'setup.json')
    const = json.load(open(setup_file))['constants'] if os.path.exists(setup_file) else V1_CONSTANTS
    p = dict(res.get('fixed', {}), **res['optimum'])
    kt = 1.38065e-5 * const['THETA']
    keq = p.get('keq', const.get('KEQ'))
    diff = dict(const['FIXED_DIFF'], k_off_0=(1.0 - p['phi']) * p['koff0'], k_catch_0=p['phi'] * p['koff0'],
                k_on_0=keq * p['koff0'], dx_kT=p['dxs_nm'] * 1e-3 / kt, dx_c_kT=p['dxc_nm'] * 1e-3 / kt)
    mech = dict(const['FIXED_MECH'], LAMBDA0=p['lambda0'], a_ratio=p['s_taut'] * p['lambda0'])
    if not os.path.exists(setup_file):
        mech.update(V1_MECH_OVERRIDES.get(name, {}))
    mech['prefdir'] = np.array(mech['prefdir'])
    if 'beta' in p:
        mech['BETA'] = p['beta']
    if 'mu0stretch' in p:
        mech['MU0STRETCH'] = p['mu0stretch']
    if 'eta' in p:
        mech['ETA'] = p['eta']
    if eta is not None:
        mech['ETA'] = eta
    return diff, mech, res['cfg']


def loop_worker(args):
    """Last cycle of the fitted LAOS run (same protocol as optimise_laos*.simulate) for the mapping check."""
    name, eta = args
    from material import DynamicCrosslinking
    diff, mech, cfg = load_run(name, eta)
    if cfg.get('plate_nodes') or cfg.get('consistent_kinetics'):   # same model as the fit (optimise_laos_v2.simulate)
        import optimise_laos_v2 as O
        p = json.load(open(os.path.join(RUNS, name, 'result.json')))['optimum']
        t, g, s = O.simulate(p, cfg)
    else:
        sim = {'deformation_type': 'shear', 'factor': cfg['factor'], 'max_iter': 50, 'tol': 1.0e-12}
        t = np.arange(0.0, cfg['cycles'] / FREQ + 0.5 * cfg['dt'], cfg['dt'])
        g = cfg['gamma_A'] * np.sin(2 * np.pi * FREQ * t)
        s = np.empty_like(t)
        with contextlib.redirect_stdout(io.StringIO()), np.errstate(all='ignore'):
            m = DynamicCrosslinking(diff, mech, sim)
            for i in range(len(t)):
                F = np.array([[1.0, g[i], 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
                s[i] = m.evaluate_state(F, None, 0.0 if i == 0 else t[i] - t[i - 1])[0, 1]
    last = t >= t[-1] - 1.0 / FREQ - 0.5 * cfg['dt']
    ref = np.load(os.path.join(RUNS, name, 'best_loop.npz'))
    n = min(len(ref['sigma']), int(last.sum()))
    return name, eta, float(np.max(np.abs(s[last][:n] - ref['sigma'][:n])))


def sweep_job(args):
    """One amplitude of one run: (run, gamma_A, G', G'', delta)."""
    import main
    name, gA, cycles, dt, factor, overrides = args
    diff, mech, cfg = load_run(name)
    cfg['cycles'] = cycles
    cfg['dt'] = dt or cfg['dt']
    cfg['factor'] = factor or cfg['factor']
    sim = {'deformation_type': 'shear', 'factor': cfg['factor'], 'max_iter': 50, 'tol': 1.0e-12}
    try:
        if cfg.get('plate_nodes') or cfg.get('consistent_kinetics'):   # same model as the fit (optimise_laos_v2.simulate)
            import optimise_laos_v2 as O
            res = json.load(open(os.path.join(RUNS, name, 'result.json')))
            p = dict(res.get('fixed', {}), **res['optimum'])
            p.update(overrides)
            if cfg.get('consistent_kinetics'):
                cfg['fortran_engine'] = 1                              # closed-system Fortran driver
            t, g, s = O.simulate(p, dict(cfg, gamma_A=gA))
            last = t >= t[-1] - 1.0 / FREQ - 0.5 * cfg['dt']
            G1, G2, delta = main.moduli(t[last], g[last], s[last], FREQ)
        else:
            _, G1, G2, delta = main.sweep_worker((gA, diff, mech, sim, FREQ, cfg['cycles'], cfg['dt'], None))
    except Exception as e:                    # keep the other amplitudes if one fails
        print(f'{name} gamma_A {gA:.4f} failed: {e}', flush=True)
        G1 = G2 = delta = np.nan
    return name, gA, G1, G2, delta


def plot_runs(res, dG1, dG2, cycles, out):
    """G' (solid) and G'' (dashed) of every run (one colour per run) against the data of Fig. 4b."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.figure(figsize=(6, 4.5))
    plt.loglog(dG1[:, 0], dG1[:, 1], 'kx', ms=6, label="Unterberger et al. 2013 $G'$")
    plt.loglog(dG2[:, 0], dG2[:, 1], 'k+', ms=7, label="Unterberger et al. 2013 $G''$")
    for j, (name, r) in enumerate(res.items()):
        plt.loglog(r[:, 0], r[:, 1], '-', color=f'C{j}', lw=1.5, label=f"{name} $G'$")
        plt.loglog(r[:, 0], r[:, 2], '--', color=f'C{j}', lw=1.5, label=f"{name} $G''$")
    plt.xlabel(r'Amplitude of shear strain $\gamma_A$')
    plt.ylabel(r"Moduli $G'$, $G''$ (Pa)")
    plt.title(f'Amplitude sweep (cycle {cycles})')
    plt.legend(fontsize=8)
    plt.grid(True, which='both', alpha=0.3)
    for axis, ticks in ((plt.gca().xaxis, [0.01, 0.02, 0.05, 0.1, 0.2, 0.5]), (plt.gca().yaxis, [1, 2, 5, 10, 20, 50, 100])):
        axis.set_major_locator(matplotlib.ticker.FixedLocator(ticks))      # plain labels at 1-2-5 positions
        axis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%g'))
        axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    plt.savefig(out, bbox_inches='tight')
    plt.close()


def main_cli():
    ap = argparse.ArgumentParser()
    ap.add_argument('--runs', nargs='+', required=True, help='run folders in optimisation_runs/')
    ap.add_argument('--cycles', type=int, default=5, help='number of cycles per amplitude (default 5)')
    ap.add_argument('--dt', type=float, default=None, help='time step [s] (default: as in the fit)')
    ap.add_argument('--factor', type=int, default=None, help='integration factor, NDIR = 10 factor^2 (default: as in the fit)')
    ap.add_argument('--max-gamma', type=float, default=MAX_GAMMA, help='largest amplitude (default 0.2; data up to about 0.5)')
    ap.add_argument('--set', nargs='+', default=[], metavar='NAME=VALUE', help='override parameters (C4 runs only)')
    ap.add_argument('--out', default=None, help='output figure (default figures/validation/0_amplitude_sweep_<runs>.pdf)')
    ap.add_argument('--workers', type=int, default=min(4, os.cpu_count() or 1))
    ap.add_argument('--verify', action='store_true', help='only reproduce best_loop.npz of each run')
    a = ap.parse_args()
    os.chdir(HERE)
    runs = a.runs

    if a.verify:
        jobs = [(r, None) for r in runs] + [(r, 0.989) for r in runs if not os.path.exists(os.path.join(RUNS, r, 'setup.json'))]
        with ProcessPoolExecutor(a.workers) as ex:
            for name, eta, err in ex.map(loop_worker, jobs):
                print(f'{name:10s} ETA {"setup" if eta is None else eta}: max |sigma - best_loop| = {err:.2e} Pa', flush=True)
        return

    overrides = {k: float(v) for k, v in (s.split('=') for s in a.set)}
    for r in runs:                            # stop early if the Fortran driver is missing or out of date
        cfg = load_run(r)[2]
        if cfg.get('consistent_kinetics'):
            import optimise_laos_v2 as O
            O.check_fortran(a.factor or cfg['factor'])

    dG1 = np.loadtxt('validation_data/4b_G_.csv', delimiter=',')
    dG2 = np.loadtxt('validation_data/4b_G__.csv', delimiter=',')
    dG2 = dG2[np.argsort(dG2[:, 0])]
    gammas = dG1[dG1[:, 0] <= 1.02 * a.max_gamma, 0]      # data amplitudes within 2 % of the limit (0.5042 for 0.5)
    jobs = [(r, g, a.cycles, a.dt, a.factor, overrides) for r in runs for g in gammas]
    out = {r: [] for r in runs}
    with ProcessPoolExecutor(a.workers) as ex:
        for name, gA, G1, G2, delta in ex.map(sweep_job, jobs):
            out[name].append((gA, G1, G2, delta))
            if len(out[name]) == len(gammas):
                cfg = load_run(name)[2]
                tag = name + ('_set' if overrides else '') + (f'_dt{a.dt:g}' if a.dt else '') + (f'_f{a.factor}' if a.factor else '')
                np.savetxt(os.path.join('validation', f'0_amplitude_sweep_{tag}.csv'), np.array(out[name])[:, :3],
                           delimiter=',', header=f"run {name} {a.set} (dt {a.dt or cfg['dt']}, factor {a.factor or cfg['factor']}, "
                                                 f"cycle {a.cycles}); gamma_A, G' (Pa), G'' (Pa)")
                print(f'{name}: done', flush=True)
    res = {r: np.array(v) for r, v in out.items()}
    fig = a.out or os.path.join('figures', 'validation',
                                '0_amplitude_sweep_' + '_'.join(runs) + ('_set' if overrides else '') + (f'_dt{a.dt:g}' if a.dt else '') + (f'_f{a.factor}' if a.factor else '') + '.pdf')
    plot_runs(res, dG1, dG2, a.cycles, fig)
    print(f'saved {fig}')


if __name__ == '__main__':
    main_cli()
