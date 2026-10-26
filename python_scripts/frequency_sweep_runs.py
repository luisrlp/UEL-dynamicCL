"""
Frequency sweep (Unterberger et al. 2013, Fig. 5: gamma_A = 0.5%, f = 0.01-10 Hz) for the optimum of saved LAOS fits.

Material properties of each run are rebuilt with amplitude_sweep_runs.load_run (result.json + setup.json).
One run per frequency from rest (fresh material), with the run's number of cycles and the run's number of time steps
per cycle (dt scaled with the period). G' and G'' of the last cycle from main.moduli, averaged over the positive and
negative peaks (cancels the start-up stress offset, which decays only on the relaxation time and remains at high f).

  python -B frequency_sweep_runs.py                                    # all runs with result.json
  python -B frequency_sweep_runs.py --runs sweep_w15_c5_g0158 c5_w15_guo
  python -B frequency_sweep_runs.py --runs c5_w15_guo --cycles 3       # fewer cycles (cheaper)

Note: Fig. 5 was measured on another rheometer and sample (Anton Paar MCR301) than Fig. 4; its G' at 0.05 Hz
(~20 Pa) is ~1.5x the amplitude-sweep G'_lin (13.5 Pa), while G'' agrees (~4.3 vs 4.4 Pa).

Outputs: validation/0_frequency_sweep_<run>.csv, figures/validation/0_frequency_sweep_<run>.pdf
"""
import argparse, contextlib, os
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from amplitude_sweep_runs import HERE, RUNS, FREQ, load_run

GAMMA_A = 0.005


def sweep_job(args):
    """One frequency of one run: (run, f, G', G'', delta)."""
    import main
    name, f, n_cycles = args
    diff, mech, cfg = load_run(name)
    sim = {'deformation_type': 'shear', 'factor': cfg['factor'], 'max_iter': 50, 'tol': 1.0e-12}
    n_cycles = n_cycles or cfg['cycles']
    dt = cfg['dt'] * FREQ / f                   # same number of steps per cycle as the fit
    try:
        t = np.arange(0.0, n_cycles / f + 0.5 * dt, dt)
        with open(os.devnull, 'w') as devnull, contextlib.redirect_stdout(devnull), np.errstate(all='ignore'):
            mat = main.DynamicCrosslinking(diff, mech, sim)
            prof = lambda tt: main.profile_cyclic(tt, amplitude=GAMMA_A, frequency=f, deformation_type='shear')
            res = main.run_simulation(mat, t, prof, mode='shear', mu=None)
        th = np.array(res['time'])
        gh = np.array(res['def_val'])
        sh = np.array([s[0, 1] for s in res['stress']])
        m = th >= th[-1] - 1.0 / f - 0.5 * dt                                    # last full cycle
        G1, G2, delta = main.moduli(th[m], gh[m], sh[m], f, both_halves=True)
    except Exception as e:                    # keep the other frequencies if one fails
        print(f'{name} f {f:.4g} Hz failed: {e}', flush=True)
        G1 = G2 = delta = np.nan
    return name, f, G1, G2, delta


def plot_run(name, f, Gp, Gpp, dG1, dG2, n_cycles):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.figure(figsize=(8, 6))
    plt.loglog(f, Gp, 'b-', lw=2, label="Model $G'$")
    plt.loglog(f, Gpp, 'r-', lw=2, label="Model $G''$")
    plt.loglog(dG1[:, 0], dG1[:, 1], 'bx', ms=6, label="Unterberger et al. 2013 $G'$")
    plt.loglog(dG2[:, 0], dG2[:, 1], 'r+', ms=7, label="Unterberger et al. 2013 $G''$")
    plt.axvline(FREQ, color='grey', ls=':', lw=1)
    plt.xlabel('Oscillation frequency $f$ (Hz)')
    plt.ylabel("Moduli $G'$, $G''$ (Pa)")
    plt.title(f"Frequency sweep ($\\gamma_A$ = {100 * GAMMA_A:g}%), parameters of optimisation run '{name}' "
              f"(cycle {n_cycles})".replace('_', '\\_'))
    plt.legend()
    plt.grid(True, which='both')
    plt.savefig(os.path.join(HERE, 'figures', 'validation', f'0_frequency_sweep_{name}.pdf'), bbox_inches='tight')
    plt.close()


def main_cli():
    ap = argparse.ArgumentParser()
    ap.add_argument('--runs', nargs='*', default=None, help='run folders (default: all with result.json)')
    ap.add_argument('--workers', type=int, default=min(6, os.cpu_count() or 1))
    ap.add_argument('--cycles', type=int, default=None, help="cycles per frequency (default: the run's cfg)")
    a = ap.parse_args()
    os.chdir(HERE)
    runs = a.runs or sorted(d for d in os.listdir(RUNS) if os.path.exists(os.path.join(RUNS, d, 'result.json')))

    dG1 = np.loadtxt('validation_data/5_G_.csv', delimiter=',')
    dG2 = np.loadtxt('validation_data/5_G__.csv', delimiter=',')
    freqs = dG1[:, 0]
    jobs = [(r, f, a.cycles) for r in runs for f in freqs]
    out = {r: [] for r in runs}
    with ProcessPoolExecutor(a.workers) as ex:
        for name, f, G1, G2, delta in ex.map(sweep_job, jobs):
            out[name].append((f, G1, G2, delta))
            if len(out[name]) == len(freqs):
                res = np.array(out[name])
                cfg = load_run(name)[2]
                n_cycles = a.cycles or cfg['cycles']
                np.savetxt(os.path.join('validation', f'0_frequency_sweep_{name}.csv'), res[:, :3], delimiter=',',
                           header=f"run {name} (gamma_A {GAMMA_A}, {1.0 / (cfg['dt'] * FREQ):g} steps per cycle, factor "
                                  f"{cfg['factor']}, cycle {n_cycles}); f (Hz), G' (Pa), G'' (Pa)")
                plot_run(name, res[:, 0], res[:, 1], res[:, 2], dG1, dG2, n_cycles)
                print(f'{name}: done', flush=True)


if __name__ == '__main__':
    main_cli()
