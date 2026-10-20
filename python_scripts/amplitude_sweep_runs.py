"""
Amplitude sweep (Unterberger et al. 2013, Fig. 4b) for the optimum of every saved LAOS fit in optimisation_runs/.

For each run, the material properties are rebuilt from result.json (optimum) and setup.json (fixed constants; runs
without setup.json use the constants of optimise_laos.py v1). One LAOS run per amplitude from rest (fresh material),
G' and G'' of the last cycle from main.moduli (Eq. (1) of the paper, positive peaks), as in main.py.

  python -B amplitude_sweep_runs.py --verify          # reproduce each run's best_loop.npz (mapping check)
  python -B amplitude_sweep_runs.py                   # sweeps of all runs
  python -B amplitude_sweep_runs.py --runs laos_v2_3  # selected runs

Outputs: validation/0_amplitude_sweep_<run>.csv, figures/validation/0_amplitude_sweep_<run>.pdf
"""
import argparse, contextlib, io, json, os
from concurrent.futures import ProcessPoolExecutor
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, 'optimisation_runs')
FREQ = 0.05
MAX_GAMMA = 0.2
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
    if eta is not None:
        mech['ETA'] = eta
    return diff, mech, res['cfg']


def loop_worker(args):
    """Last cycle of the fitted LAOS run (same protocol as optimise_laos*.simulate) for the mapping check."""
    name, eta = args
    from material import DynamicCrosslinking
    diff, mech, cfg = load_run(name, eta)
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
    name, gA = args
    diff, mech, cfg = load_run(name)
    sim = {'deformation_type': 'shear', 'factor': cfg['factor'], 'max_iter': 50, 'tol': 1.0e-12}
    try:
        _, G1, G2, delta = main.sweep_worker((gA, diff, mech, sim, FREQ, cfg['cycles'], cfg['dt'], None))
    except Exception as e:                    # keep the other amplitudes if one fails
        print(f'{name} gamma_A {gA:.4f} failed: {e}', flush=True)
        G1 = G2 = delta = np.nan
    return name, gA, G1, G2, delta


def plot_run(name, gA, Gp, Gpp, dG1, dG2, cfg):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.figure(figsize=(8, 6))
    plt.loglog(gA, Gp, 'b-', lw=2, label="Model $G'$")
    plt.loglog(gA, Gpp, 'r-', lw=2, label="Model $G''$")
    plt.loglog(dG1[:, 0], dG1[:, 1], 'bx', ms=6, label="Unterberger et al. 2013 $G'$")
    plt.loglog(dG2[:, 0], dG2[:, 1], 'r+', ms=7, label="Unterberger et al. 2013 $G''$")
    plt.xlabel('Amplitude of shear strain $\\gamma_A$')
    plt.ylabel("Moduli $G'$, $G''$ (Pa)")
    plt.title(f"Amplitude sweep, parameters of optimisation run '{name}' (cycle {cfg['cycles']})".replace('_', '\\_'))
    plt.legend()
    plt.grid(True, which='both')
    plt.savefig(os.path.join(HERE, 'figures', 'validation', f'0_amplitude_sweep_{name}.pdf'), bbox_inches='tight')
    plt.close()


def main_cli():
    ap = argparse.ArgumentParser()
    ap.add_argument('--runs', nargs='*', default=None, help='run folders (default: all with result.json)')
    ap.add_argument('--workers', type=int, default=min(15, os.cpu_count() or 1))
    ap.add_argument('--verify', action='store_true', help='only reproduce best_loop.npz of each run')
    a = ap.parse_args()
    os.chdir(HERE)
    runs = a.runs or sorted(d for d in os.listdir(RUNS) if os.path.exists(os.path.join(RUNS, d, 'result.json')))

    if a.verify:
        jobs = [(r, None) for r in runs] + [(r, 0.989) for r in runs if not os.path.exists(os.path.join(RUNS, r, 'setup.json'))]
        with ProcessPoolExecutor(a.workers) as ex:
            for name, eta, err in ex.map(loop_worker, jobs):
                print(f'{name:10s} ETA {"setup" if eta is None else eta}: max |sigma - best_loop| = {err:.2e} Pa', flush=True)
        return

    dG1 = np.loadtxt('validation_data/4b_G_.csv', delimiter=',')
    dG2 = np.loadtxt('validation_data/4b_G__.csv', delimiter=',')
    dG2 = dG2[np.argsort(dG2[:, 0])]
    gammas = dG1[dG1[:, 0] <= MAX_GAMMA, 0]
    jobs = [(r, g) for r in runs for g in gammas]
    out = {r: [] for r in runs}
    with ProcessPoolExecutor(a.workers) as ex:
        for name, gA, G1, G2, delta in ex.map(sweep_job, jobs):
            out[name].append((gA, G1, G2, delta))
            if len(out[name]) == len(gammas):
                res = np.array(out[name])
                cfg = load_run(name)[2]
                np.savetxt(os.path.join('validation', f'0_amplitude_sweep_{name}.csv'), res[:, :3], delimiter=',',
                           header=f"run {name} (dt {cfg['dt']}, factor {cfg['factor']}, cycle {cfg['cycles']}); "
                                  f"gamma_A, G' (Pa), G'' (Pa)")
                plot_run(name, res[:, 0], res[:, 1], res[:, 2], dG1, dG2, cfg)
                print(f'{name}: done', flush=True)


if __name__ == '__main__':
    main_cli()
