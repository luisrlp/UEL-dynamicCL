"""Least-squares fit of the stretch-free-binding model to the frequency sweep of Unterberger et al. (2013), Fig. 5
(gamma_A = 0.5 %), alone: no LAOS loop, no amplitude sweep.

Model : python_scripts/material.py (DynamicCrosslinking), simple shear at one material point, closed system,
        gamma(t) = gamma_A sin(2 pi f t) from rest, one run per fitted frequency, --steps time steps per cycle and
        --cycles cycles; G', G'' of the last cycle from main.moduli (both half-cycles, as frequency_sweep_runs.py).
Data  : validation_data/5_G_.csv, 5_G__.csv at the frequencies --freqs (Hz; must be data points, f <= 0.8 Hz
        by default, below the high-frequency G'' rise).
Residuals (2 per frequency): log10(G'_model / G'_data), log10(G''_model / G''_data). With --shape, G' and G'' are
        both divided by G'(--f-ref) of the model and of the data (removes the sample offset; the G' row at f_ref is
        then zero and dropped). A failed run gives a log10 error of FAIL.
Parameters, bounds, material constants and the parameter mapping are those of optimise_laos_v2.py (PARAMS,
DEFAULT_FIXED, material_props). The Bell distances barely act at 0.5 % strain, so fix dxs_nm and dxc_nm.

Usage (results go to --out):
  python optimise_freq.py --fix dxs_nm=1.1325 --fix dxc_nm=1.5 --fix phi=0.907 --workers 16 --out optimisation_runs/freq_F1
  python optimise_freq.py ... --free beta --free mu0stretch          # also fit the force-law constants
  python optimise_freq.py ... --shape                                 # shape only (normalised at G'(f_ref))
  python optimise_freq.py ... --evaluate-only                         # one evaluation at x0
Outputs (in --out): evals.csv, cache.json, result.json, setup.json, setup_<start time>.json. result.json and
setup.json can be read by amplitude_sweep_runs.load_run, so frequency_sweep_runs.py and amplitude_sweep_runs.py
work on the optimum (cfg['dt'] is the 0.05 Hz step with the same number of steps per cycle).
"""
import argparse, contextlib, hashlib, io, json, math, os, sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import scipy
from scipy.optimize import least_squares

REPO = os.path.dirname(os.path.abspath(__file__))     # python_scripts/
sys.path.insert(0, REPO)
import optimise_laos_v2 as O

FREQ_DATA = (os.path.join(REPO, 'validation_data', '5_G_.csv'),    # Fig. 5: G' and G'' vs f
             os.path.join(REPO, 'validation_data', '5_G__.csv'))
GAMMA_A = 0.005
FREQS = (0.01, 0.0207, 0.0428, 0.0886, 0.1833, 0.3793, 0.7848)   # every second data point up to 0.8 Hz
F_REF = 0.0428                                                    # normalisation frequency for --shape
FAIL = 2.0                                                        # log10 error of a failed run


def load_freq_data(freqs):
    """G' and G'' of Fig. 5 at the data frequencies closest to freqs; returns (freqs at data points, (n, 2))."""
    d1, d2 = (np.loadtxt(f, delimiter=',', comments='#') for f in FREQ_DATA)
    idx = [int(np.argmin(np.abs(np.log(d1[:, 0] / f)))) for f in freqs]
    for f, i in zip(freqs, idx):
        if abs(math.log(d1[i, 0] / f)) > 0.02:
            sys.exit(f'frequency {f} Hz is not a data point of Fig. 5 (nearest {d1[i, 0]:.4g} Hz)')
    return d1[idx, 0], np.column_stack([d1[idx, 1], d2[idx, 1]])


def run_freq(args):
    """Worker: one oscillation from rest at frequency f. Returns (G', G'', wall time); NaN if the run fails."""
    from main import moduli
    from material import DynamicCrosslinking
    p, cfg, f = args
    t0 = time.time()
    try:
        diff, mech = O.material_props(p)
        sim = {'deformation_type': 'shear', 'factor': cfg['factor'], 'max_iter': 50, 'tol': 1.0e-12}
        dt = 1.0 / (f * cfg['steps'])
        t = np.arange(0.0, cfg['cycles'] / f + 0.5 * dt, dt)
        g = cfg['gamma_A'] * np.sin(2 * np.pi * f * t)
        s = np.empty_like(t)
        with contextlib.redirect_stdout(io.StringIO()), np.errstate(all='ignore'):
            m = DynamicCrosslinking(diff, mech, sim)
            for i in range(len(t)):
                F = np.array([[1.0, g[i], 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
                s[i] = m.evaluate_state(F, None, 0.0 if i == 0 else t[i] - t[i - 1])[0, 1]
        last = t >= t[-1] - 1.0 / f - 0.5 * dt
        G1, G2, _ = moduli(t[last], g[last], s[last], f, both_halves=True)
        if not (np.isfinite(G1) and np.isfinite(G2)):
            raise FloatingPointError('non-finite moduli')
    except Exception as e:
        print(f'   ! run failed at f {f} Hz for {p}: {e!r}', flush=True)
        G1 = G2 = float('nan')
    return G1, G2, time.time() - t0


def settings_signature(cfg):
    """Hash of the run settings, fixed constants, data, model source and moduli evaluation (part of every cache key)."""
    import inspect
    from main import moduli
    h = hashlib.sha1()
    fixed_mech = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in O.FIXED_MECH.items()}
    h.update(json.dumps({'cfg': cfg, 'diff': O.FIXED_DIFF, 'mech': fixed_mech, 'FAIL': FAIL}, sort_keys=True).encode())
    for f in ('material.py', 'kinematics.py', 'dfd.py'):
        h.update(open(os.path.join(REPO, f), 'rb').read())
    for f in FREQ_DATA:
        h.update(open(f, 'rb').read())
    h.update(inspect.getsource(moduli).encode())
    return h.hexdigest()[:12]


class Problem:
    def __init__(self, free, fixed, cfg, out, workers, h):
        self.free, self.fixed, self.cfg, self.out, self.h = free, fixed, cfg, out, h
        self.freqs = list(cfg['freqs'])
        _, self.Gd = load_freq_data(self.freqs)                       # (n_freqs, 2): G', G'' data
        self.iref = self.freqs.index(cfg['f_ref']) if cfg['shape'] else None
        self.lo = np.array([O.PARAMS[k][1] for k in free])
        self.hi = np.array([O.PARAMS[k][2] for k in free])
        self.pool = ProcessPoolExecutor(max_workers=workers)
        self.cache_file = os.path.join(out, 'cache.json')
        self.cache = json.load(open(self.cache_file)) if os.path.exists(self.cache_file) else {}
        self.signature = settings_signature(cfg)
        self.n_points = 0
        self.log = os.path.join(out, 'evals.csv')
        header = ','.join(['point'] + list(O.PARAMS) + ['rms_log10'] + [f'{q}_f{f:.4g}' for q in ('G1', 'G2') for f in self.freqs] +
                          ['wall_s'])
        if not os.path.exists(self.log):
            with open(self.log, 'w') as fh:
                fh.write(header + '\n')
        elif open(self.log).readline().strip() != header:
            sys.exit(f'{self.log} has different columns (other frequencies): use a new --out')

    def phys(self, z):
        p = dict(self.fixed)
        p.update({k: float(lo + zi * (hi - lo)) for k, zi, lo, hi in zip(self.free, z, self.lo, self.hi)})
        return p

    def key(self, p):
        return self.signature + ' ' + json.dumps({k: round(p[k], 10) for k in O.PARAMS}, sort_keys=True)

    def residuals(self, G):
        """Residual vector from the model moduli G (n_freqs, 2): absolute or shape (normalised at f_ref)."""
        ok = np.isfinite(G) & (G > 0)
        e = np.where(ok, np.log10(np.where(ok, G, 1.0) / self.Gd), FAIL)          # (n, 2)
        if self.iref is None:
            return e.T.ravel()                                                     # G' rows, then G''
        if not ok[self.iref, 0]:
            return np.full(2 * len(self.freqs) - 1, FAIL)
        e = e - e[self.iref, 0]                         # log10[(G/G'_ref)_model / (G/G'_ref)_data]
        return np.concatenate([np.delete(e[:, 0], self.iref), e[:, 1]])

    def evaluate(self, zs):
        """Moduli and residual vectors for a list of normalised points (cached; missing points run in parallel)."""
        ps = [self.phys(z) for z in zs]
        todo = list({self.key(p): p for p in ps if self.key(p) not in self.cache}.values())
        if todo:
            futs = [[self.pool.submit(run_freq, (p, self.cfg, f)) for f in self.freqs] for p in todo]
            for p, fs in zip(todo, futs):
                res = np.array([fu.result() for fu in fs])                       # rows: G', G'', wall
                G = res[:, :2]
                r = self.residuals(G)
                self.cache[self.key(p)] = {'G': G.tolist(), 'r': r.tolist()}
                self.n_points += 1
                with open(self.log, 'a') as fh:
                    fh.write(','.join([str(self.n_points)] + [f'{p[q]:.6g}' for q in O.PARAMS] +
                                      [f'{np.sqrt(np.mean(r**2)):.5f}'] + [f'{x:.5g}' for x in G[:, 0]] +
                                      [f'{x:.5g}' for x in G[:, 1]] + [f'{res[:, 2].max():.0f}']) + '\n')
            json.dump(self.cache, open(self.cache_file, 'w'))
        return [self.cache[self.key(p)] for p in ps]

    def fun(self, z):
        c = self.evaluate([z])[0]
        r = np.array(c['r'])
        print(f'   f: {self.fmt(z)} -> RMS log10 {np.sqrt(np.mean(r**2)):.4f}  (points so far: {self.n_points})', flush=True)
        return r

    def jac(self, z):
        steps = [self.h if zi + self.h <= 1.0 else -self.h for zi in z]
        pts = [z] + [z + s * e for s, e in zip(steps, np.eye(len(z)))]
        rs = [np.array(c['r']) for c in self.evaluate(pts)]
        return np.column_stack([(rs[i + 1] - rs[0]) / steps[i] for i in range(len(z))])

    def fmt(self, z):
        p = self.phys(z)
        return ', '.join(f'{k} {p[k]:.4g}' for k in self.free)


def print_table(prob, G):
    print(f"   f (Hz) |  G' model  G' data | G'' model G'' data   [Pa]")
    for f, (G1, G2), (d1, d2) in zip(prob.freqs, G, prob.Gd):
        print(f'   {f:6.4f} | {G1:9.3f} {d1:8.3f} | {G2:9.3f} {d2:9.3f}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default=os.path.join(REPO, 'optimisation_runs', 'freq'))
    ap.add_argument('--fix', action='append', default=[], help='name=value (repeatable)')
    ap.add_argument('--free', action='append', default=[], help=f'release a default-fixed parameter {O.DEFAULT_FIXED} (repeatable)')
    ap.add_argument('--x0', action='append', default=[], help='name=value to override a starting value')
    ap.add_argument('--freqs', type=float, nargs='+', default=list(FREQS), help='fitted frequencies (Fig. 5 data points, Hz)')
    ap.add_argument('--shape', action='store_true', help="fit the shape only: G', G'' normalised by G'(--f-ref)")
    ap.add_argument('--f-ref', type=float, default=F_REF, help='normalisation frequency of --shape (one of --freqs)')
    ap.add_argument('--steps', type=int, default=200, help='time steps per cycle')
    ap.add_argument('--cycles', type=int, default=4, help='cycles from rest; the last one is evaluated')
    ap.add_argument('--factor', type=int, default=3)
    ap.add_argument('--workers', type=int, default=7)
    ap.add_argument('--max-nfev', type=int, default=25)
    ap.add_argument('--h', type=float, default=0.01, help='finite-difference step in normalised variables')
    ap.add_argument('--evaluate-only', action='store_true')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    freqs, _ = load_freq_data(sorted(a.freqs))
    freqs = [float(f) for f in freqs]
    f_ref = float(load_freq_data([a.f_ref])[0][0])
    if a.shape and f_ref not in freqs:
        sys.exit(f'--f-ref {a.f_ref} must be one of --freqs')
    # cfg: 'dt' and 'cycles' as read by frequency_sweep_runs / amplitude_sweep_runs (dt of the 0.05 Hz cycle)
    cfg = {'gamma_A': GAMMA_A, 'freqs': freqs, 'shape': bool(a.shape), 'f_ref': f_ref if a.shape else None,
           'steps': a.steps, 'cycles': a.cycles, 'factor': a.factor, 'dt': 1.0 / (O.FREQ * a.steps)}
    x0 = {k: O.PARAMS[k][0] for k in O.PARAMS}
    x0.update({k: float(v) for k, v in (s.split('=') for s in a.x0)})
    for k in a.free:
        if k not in O.PARAMS:
            sys.exit(f'--free {k}: unknown parameter')
    user_fixed = {k: float(v) for k, v in (s.split('=') for s in a.fix)}
    if set(a.free) & set(user_fixed):
        sys.exit(f'parameters both freed and fixed: {set(a.free) & set(user_fixed)}')
    fixed = {k: x0[k] for k in O.DEFAULT_FIXED if k not in a.free}
    fixed.update(user_fixed)
    free = [k for k in O.PARAMS if k not in fixed]
    for k in free:
        if not O.PARAMS[k][1] <= x0[k] <= O.PARAMS[k][2]:
            sys.exit(f'x0 of {k} = {x0[k]} outside bounds {O.PARAMS[k][1:]}')
    print(f'cfg = {cfg}\nfree = {free}\nfixed = {fixed}', flush=True)

    prob = Problem(free, {**{k: x0[k] for k in O.PARAMS}, **fixed}, cfg, a.out, a.workers, a.h)
    z0 = np.array([(x0[k] - O.PARAMS[k][1]) / (O.PARAMS[k][2] - O.PARAMS[k][1]) for k in free])

    started = time.strftime('%Y-%m-%d_%H%M%S')
    setup = {'started': started, 'command': ' '.join(sys.argv), 'signature': prob.signature,
             'data_files': list(FREQ_DATA), 'cfg': cfg, 'free': free, 'fixed': fixed, 'default_fixed': list(O.DEFAULT_FIXED),
             'x0': x0, 'bounds': {k: [O.PARAMS[k][1], O.PARAMS[k][2]] for k in O.PARAMS},
             'constants': {'FIXED_DIFF': O.FIXED_DIFF,
                           'FIXED_MECH': {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in O.FIXED_MECH.items()},
                           'FREQ': O.FREQ, 'THETA': O.THETA},
             'FAIL': FAIL,
             'optimiser': {**O.LSQ, 'max_nfev': a.max_nfev, 'h': a.h, 'workers': a.workers, 'evaluate_only': a.evaluate_only},
             'versions': {'numpy': np.__version__, 'scipy': scipy.__version__, 'python': sys.version.split()[0]}}
    for name in ('setup.json', f'setup_{started}.json'):
        json.dump(setup, open(os.path.join(a.out, name), 'w'), indent=2)

    if a.evaluate_only:
        r = prob.fun(z0)
        print_table(prob, np.array(prob.evaluate([z0])[0]['G']))
        prob.pool.shutdown()
        return

    t0 = time.time()
    res = least_squares(prob.fun, z0, jac=prob.jac, bounds=(0.0, 1.0), max_nfev=a.max_nfev, verbose=1, **O.LSQ)
    p_opt = prob.phys(res.x)
    J = res.jac
    dof = max(len(res.fun) - len(free), 1)
    s2 = 2 * res.cost / dof
    try:
        cov = s2 * np.linalg.inv(J.T @ J) * np.outer(prob.hi - prob.lo, prob.hi - prob.lo)
        se = np.sqrt(np.diag(cov))
        corr = cov / np.outer(se, se)
    except np.linalg.LinAlgError:
        se, corr = np.full(len(free), np.nan), np.full((len(free), len(free)), np.nan)
    G = np.array(prob.evaluate([res.x])[0]['G'])
    out = {'status': res.status, 'message': res.message, 'nfev': res.nfev, 'njev': res.njev, 'points': prob.n_points,
           'wall_h': (time.time() - t0) / 3600, 'cfg': cfg, 'free': free, 'fixed': fixed, 'optimum': p_opt,
           'a_ratio': p_opt['s_taut'] * p_opt['lambda0'],
           'k_off_0 (slip)': (1 - p_opt['phi']) * p_opt['koff0'], 'k_catch_0': p_opt['phi'] * p_opt['koff0'],
           'k_on_0': p_opt['keq'] * p_opt['koff0'],
           'rms_log10': float(np.sqrt(np.mean(res.fun**2))),
           'freq': {'f': prob.freqs, 'G1_model': G[:, 0].tolist(), 'G1_data': prob.Gd[:, 0].tolist(),
                    'G2_model': G[:, 1].tolist(), 'G2_data': prob.Gd[:, 1].tolist()},
           'std_error': dict(zip(free, se.tolist())), 'correlation': corr.tolist(),
           'at_bound': [k for k, z in zip(free, res.x) if z < 1e-3 or z > 1 - 1e-3]}
    json.dump(out, open(os.path.join(a.out, 'result.json'), 'w'), indent=2)

    print(f'\n{res.message}  (nfev {res.nfev}, njev {res.njev}, points {prob.n_points}, {out["wall_h"]:.2f} h)')
    print(f'RMS log10 {out["rms_log10"]:.4f}' + (' (shape)' if a.shape else ''))
    for k, e in zip(free, se):
        print(f'  {k:10s} = {p_opt[k]:.5g}  +- {e:.2g}' + ('   (at bound)' if k in out['at_bound'] else ''))
    print('correlation matrix:\n' + '\n'.join(f'  {k:10s} ' + ' '.join(f'{c:+.2f}' for c in row) for k, row in zip(free, corr)))
    print_table(prob, G)
    prob.pool.shutdown()


if __name__ == '__main__':
    main()
