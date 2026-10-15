"""Least-squares fit of the stretch-free-binding model to the LAOS loop of Unterberger et al. (2013), Fig. 4a.

Model : python_scripts/material.py (DynamicCrosslinking), simple shear at one material point, closed system,
        gamma(t) = gamma_A sin(2 pi f t), f = 0.05 Hz; the last simulated cycle is compared with the data.
Data  : python_scripts/validation_data/unterberger2013_fig4a_laos_dataset0.csv (gamma, sigma_xy, branch).
Method: scipy.optimize.least_squares (trust-region reflective, bounds) in normalised parameters z in [0, 1];
        the forward-difference Jacobian is evaluated in parallel (one LAOS run per free parameter).

Free parameters (defaults; any of them can be fixed with --fix name=value):
  lambda0  prestretch LAMBDA0
  s_taut   a_ratio / LAMBDA0 (affine stretch at which a filament becomes taut; decorrelates LAMBDA0 and a_ratio)
  koff0    unloaded off-rate k_s + k_c [1/s]
  phi      catch fraction k_c / (k_s + k_c)
  dxc_nm   catch reactive distance [nm]
  dxs_nm   slip reactive distance [nm]
Held fixed: ETA, KEQ = k_on_0 / koff0 (so c_b0 and the elastic state at rest do not depend on the kinetics), and all
other material constants (see FIXED_DIFF / FIXED_MECH).

Residuals: model - data at the data points of both branches, weighted for the extraction uncertainty
(delta_sigma, delta_gamma): r_w = r * delta_sigma / sqrt(delta_sigma^2 + (dsigma/dgamma * delta_gamma)^2),
which down-weights the near-vertical tips. Units stay Pa (r_w = r where the loop is flat).

Usage (from any directory; results go to --out, default python_scripts/optimisation_runs/laos):
  python optimise_laos.py --evaluate-only                 # one LAOS run at x0, prints the residual summary
  python optimise_laos.py --workers 7 --max-nfev 20        # fit
  python optimise_laos.py --fix dxs_nm=0.5 --workers 6     # fit with dxs fixed
Outputs (in --out): evals.csv (every LAOS run), cache.json (residuals per parameter set, reused on restart
                    only with the same settings, fixed constants, data and model source),
                    result.json (optimum, RMS, correlation matrix, standard errors), best_loop.npz,
                    setup.json (full setup of the latest session, written before the first LAOS run) and
                    setup_<start time>.json (one copy per session, kept when the folder is reused).
"""
import argparse, contextlib, hashlib, io, json, math, os, sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import scipy
from scipy.optimize import least_squares

REPO = os.path.dirname(os.path.abspath(__file__))     # python_scripts/
DATA = os.path.join(REPO, 'validation_data', 'unterberger2013_fig4a_laos_dataset0.csv')
sys.path.insert(0, REPO)

THETA = 294.0
KT = 1.38065e-5 * THETA            # pN um
KEQ = 5.0                          # k_on_0 / koff0, held fixed
FREQ = 0.05                        # Hz

# name: (x0, lower, upper)
PARAMS = {
    'lambda0': (1.0115, 1.000, 1.050), # 1.0467
    's_taut':  (1.0776, 1.050, 1.200), # 1.0861 a_ratio/lambda0
    'koff0':   (0.24991,  0.020, 5.00), # [0.08] k_c + k_s
    'phi':     (0.98578, 0.1, 0.999), # [0.9091] k_c / (k_s + k_c) = k_c / koff0
    'dxc_nm':  (17.209,    0.5,   40.0), # 2.
    'dxs_nm':  (0.5419,    0.2,   40.0), # 0.5
}

FIXED_DIFF = {'CHI': 0.1, 
              'THETA': THETA, 
              'D': 0.1, 
              'MU0': 0.0, 
              'VMOL': 0.15, 
              'RGAS': 8.31446261815324,
              'CACTIN': 0.0095, 
              'R': 0.1, 
              'RBMAX': 0.25, 
              'RFMAX': 1.0}

FIXED_MECH = {'filament': 'ext', 
              'KBULK': 1.e3, 
              'MU0STRETCH': 38600.0,
              'Lp': 16.0, 
              'BETA': 0.438, 
              'R0C': 0.014,
              'ETA': 1., 
              'b_disp': 1e-6, 
              'prefdir': np.array([1.0, 0.0, 0.0]),
              'stretch_free_binding': True, 
              'xi_t': 1.0}

DSIG, DGAM = 0.025, 0.0007         # extraction uncertainty (Pa, -)
LSQ = {'method': 'trf', 'x_scale': 1.0, 'ftol': 1e-4, 'xtol': 1e-4, 'gtol': 1e-6}   # least_squares settings


# ---------------------------------------------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------------------------------------------
def load_data():
    # Reads .txt, .csv, .dat
    d = np.loadtxt(DATA, delimiter=',', comments='#')
    out = {}
    for b in (1, -1): # 1: loading, -1: unloading
        x = d[d[:, 2] == b]
        o = np.argsort(x[:, 0], kind='stable')
        gamma, stress = x[o, 0], x[o, 1]
        # Aux variables for extraction uncertainty
        slope = np.gradient(stress, gamma + 1e-12 * np.arange(len(gamma)))      # tiny offset: two tip points share gamma
        w = DSIG / np.sqrt(DSIG**2 + (slope * DGAM)**2)
        out[b] = (gamma, stress, w)
    # Strain amplitude
    gamma_A = 0.5 * (d[:, 0].max() - d[:, 0].min())
    return out, gamma_A


# ---------------------------------------------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------------------------------------------
def material_props(p):
    koff0, phi = p['koff0'], p['phi']
    diff = dict(FIXED_DIFF, k_off_0=(1.0 - phi) * koff0, k_catch_0=phi * koff0, k_on_0=KEQ * koff0,
                dx_kT=p['dxs_nm'] * 1e-3 / KT, dx_c_kT=p['dxc_nm'] * 1e-3 / KT)
    mech = dict(FIXED_MECH, LAMBDA0=p['lambda0'], a_ratio=p['s_taut'] * p['lambda0'])
    return diff, mech


def simulate(p, cfg):
    """Strain-controlled LAOS; returns gamma and sigma_xy of the last cycle."""
    from material import DynamicCrosslinking
    diff, mech = material_props(p)
    sim = {'deformation_type': 'shear', 'factor': cfg['factor'], 'max_iter': 50, 'tol': 1.0e-12}
    period = 1.0 / FREQ
    t = np.arange(0.0, cfg['cycles'] * period + 0.5 * cfg['dt'], cfg['dt'])
    g = cfg['gamma_A'] * np.sin(2 * np.pi * FREQ * t)
    s = np.empty_like(t)
    with contextlib.redirect_stdout(io.StringIO()), np.errstate(all='ignore'):
        m = DynamicCrosslinking(diff, mech, sim)
        for i in range(len(t)):
            F = np.array([[1.0, g[i], 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
            s[i] = m.evaluate_state(F, None, 0.0 if i == 0 else t[i] - t[i - 1])[0, 1]
    last = t >= (cfg['cycles'] - 1) * period - 1e-9
    return g[last], s[last]


def residuals_from_loop(gm, sm, data):
    gd = np.gradient(gm)
    r = []
    for b in (1, -1):
        msk = gd > 0 if b == 1 else gd < 0
        o = np.argsort(gm[msk])
        g, s, w = data[b]
        # Interpolate the model stress at the data points of branch b, 
        # and weight by the extraction uncertainty
        r.append(w * (np.interp(g, gm[msk][o], sm[msk][o]) - s))
    return np.concatenate(r)


def run_one(args):
    """Worker: one LAOS run for the physical parameter dict p. 
    Returns (residual vector, loop, wall time)."""
    p, cfg, data = args
    t0 = time.time()
    try:
        gm, sm = simulate(p, cfg)
        r = residuals_from_loop(gm, sm, data)
        if not np.all(np.isfinite(r)):
            raise FloatingPointError('non-finite stress')
    except Exception as e:                           # failed run: large, finite penalty
        print(f'   ! run failed at {p}: {e!r}', flush=True)
        gm = sm = np.array([])
        r = np.full(sum(len(data[b][0]) for b in (1, -1)), 10.0)
    return r, (gm, sm), time.time() - t0


# ---------------------------------------------------------------------------------------------------------------
# Optimiser wrapper: normalised variables, cache, parallel forward-difference Jacobian
# ---------------------------------------------------------------------------------------------------------------
def settings_signature(cfg):
    """Hash of everything besides PARAMS that changes the residuals: run settings, fixed constants, weights,
    data file and model source. It is part of every cache key, so residuals computed with other settings
    (or before an edit of material.py) are never reused."""
    h = hashlib.sha1()
    fixed_mech = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in FIXED_MECH.items()}
    h.update(json.dumps({'cfg': {k: float(v) for k, v in cfg.items()}, 'diff': FIXED_DIFF, 'mech': fixed_mech,
                         'KEQ': KEQ, 'FREQ': FREQ, 'DSIG': DSIG, 'DGAM': DGAM}, sort_keys=True).encode())
    for f in ('material.py', 'kinematics.py', 'dfd.py'):
        h.update(open(os.path.join(REPO, f), 'rb').read())
    h.update(open(DATA, 'rb').read())
    return h.hexdigest()[:12]


class Problem:
    def __init__(self, free, fixed, cfg, data, out, workers, h):
        self.free, self.fixed, self.cfg, self.data, self.out, self.h = free, fixed, cfg, data, out, h
        self.lo = np.array([PARAMS[k][1] for k in free])
        self.hi = np.array([PARAMS[k][2] for k in free])
        self.pool = ProcessPoolExecutor(max_workers=workers)
        self.cache_file = os.path.join(out, 'cache.json')
        self.cache = json.load(open(self.cache_file)) if os.path.exists(self.cache_file) else {}
        self.signature = settings_signature(cfg)
        self.loops = {}
        self.n_runs = 0
        self.log = os.path.join(out, 'evals.csv')
        if not os.path.exists(self.log):
            with open(self.log, 'w') as fh:
                fh.write(','.join(['run'] + list(PARAMS) + ['rms_weighted', 'rms_raw', 'wall_s']) + '\n')

    def phys(self, z):
        # Normalised -> physical parameters
        p = dict(self.fixed)
        p.update({k: float(lo + zi * (hi - lo)) for k, zi, lo, hi in zip(self.free, z, self.lo, self.hi)})
        return p

    def key(self, p):
        return self.signature + ' ' + json.dumps({k: round(p[k], 10) for k in PARAMS}, sort_keys=True)

    def evaluate(self, zs):
        """Residual vectors for a list of normalised points 
        (cached; missing ones run in parallel)."""
        # Requested points in physical parameters
        ps = [self.phys(z) for z in zs]
        # Keep only the ones not yet in cache, without duplicates
        todo = [p for p in ps if self.key(p) not in self.cache]
        todo = list({self.key(p): p for p in todo}.values())
        if todo:
            for p, (r, loop, wall) in zip(todo, self.pool.map(run_one, [(p, self.cfg, self.data) for p in todo])):
                k = self.key(p)
                self.cache[k] = r.tolist()
                self.loops[k] = loop
                self.n_runs += 1
                with open(self.log, 'a') as fh:
                    fh.write(','.join([str(self.n_runs)] + [f'{p[q]:.6g}' for q in PARAMS] +
                                      [f'{np.sqrt(np.mean(r**2)):.5f}', f'{self.raw_rms(loop):.5f}', f'{wall:.0f}']) + '\n')
            json.dump(self.cache, open(self.cache_file, 'w'))
        return [np.array(self.cache[self.key(p)]) for p in ps]

    def raw_rms(self, loop):
        gm, sm = loop
        if len(gm) == 0:
            return float('nan')
        # data points with unit weights
        unweighted = {b: (g, s, np.ones_like(g)) for b, (g, s, w) in self.data.items()}
        return float(np.sqrt(np.mean(residuals_from_loop(gm, sm, unweighted)**2)))

    def fun(self, z):
        r = self.evaluate([z])[0]
        print(f'   f: {self.fmt(z)} -> RMS_w {np.sqrt(np.mean(r**2)):.4f} Pa  (LAOS runs so far: {self.n_runs})', flush=True)
        return r

    def jac(self, z):
        # Forward-difference Jacobian: each free parameter j is perturbed alone by h (1% of its bound range, backwards at the upper bound) 
        # while the others stay at z, one LAOS run per parameter in parallel, column j = (r(z + h e_j) - r(z)) / h with r(z) taken from the cache.
        steps = [self.h if zi + self.h <= 1.0 else -self.h for zi in z]
        pts = [z] + [z + s * e for s, e in zip(steps, np.eye(len(z)))]
        rs = self.evaluate(pts)
        return np.column_stack([(rs[i + 1] - rs[0]) / steps[i] for i in range(len(z))])

    def fmt(self, z):
        # Free parameter formatting
        p = self.phys(z)
        return ', '.join(f'{k} {p[k]:.4g}' for k in self.free)


# ---------------------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default=os.path.join(REPO, 'optimisation_runs', 'laos'))
    ap.add_argument('--fix', action='append', default=[], help='name=value, e.g. dxs_nm=0.5 (repeatable)')
    ap.add_argument('--x0', action='append', default=[], help='name=value to override a starting value')
    ap.add_argument('--dt', type=float, default=0.1)
    ap.add_argument('--factor', type=int, default=3)
    ap.add_argument('--cycles', type=int, default=3)
    ap.add_argument('--workers', type=int, default=7)
    ap.add_argument('--max-nfev', type=int, default=25, help='max residual evaluations of least_squares (Jacobians not counted)')
    ap.add_argument('--h', type=float, default=0.01, help='finite-difference step in normalised variables (fraction of the bound range)')
    ap.add_argument('--evaluate-only', action='store_true')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    data, gamma_A = load_data()
    cfg = {'dt': a.dt, 'factor': a.factor, 'cycles': a.cycles, 'gamma_A': gamma_A}
    fixed = {k: float(v) for k, v in (s.split('=') for s in a.fix)}
    x0 = {k: PARAMS[k][0] for k in PARAMS}
    x0.update({k: float(v) for k, v in (s.split('=') for s in a.x0)})
    free = [k for k in PARAMS if k not in fixed]
    for k in free:
        if not PARAMS[k][1] <= x0[k] <= PARAMS[k][2]:
            sys.exit(f'x0 of {k} = {x0[k]} outside bounds {PARAMS[k][1:]}')
    print(f'gamma_A = {gamma_A:.5f}, cfg = {cfg}\nfree = {free}\nfixed = {fixed}', flush=True)

    prob = Problem(free, {**{k: x0[k] for k in PARAMS}, **fixed}, cfg, data, a.out, a.workers, a.h)
    z0 = np.array([(x0[k] - PARAMS[k][1]) / (PARAMS[k][2] - PARAMS[k][1]) for k in free])

    # full setup of this session, written before any LAOS run (so it also exists for interrupted fits)
    started = time.strftime('%Y-%m-%d_%H%M%S')
    setup = {'started': started, 'command': ' '.join(sys.argv), 'signature': prob.signature,
             'data_file': DATA, 'cfg': cfg, 'free': free, 'fixed': fixed,
             'x0': x0, 'bounds': {k: [PARAMS[k][1], PARAMS[k][2]] for k in PARAMS},
             'constants': {'FIXED_DIFF': FIXED_DIFF,
                           'FIXED_MECH': {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in FIXED_MECH.items()},
                           'KEQ': KEQ, 'FREQ': FREQ, 'THETA': THETA},
             'weights': {'DSIG': DSIG, 'DGAM': DGAM},
             'optimiser': {**LSQ, 'max_nfev': a.max_nfev, 'h': a.h, 'workers': a.workers, 'evaluate_only': a.evaluate_only},
             'versions': {'numpy': np.__version__, 'scipy': scipy.__version__, 'python': sys.version.split()[0]}}
    for name in ('setup.json', f'setup_{started}.json'):
        json.dump(setup, open(os.path.join(a.out, name), 'w'), indent=2)

    if a.evaluate_only:
        r = prob.fun(z0)
        k = prob.key(prob.phys(z0))
        print(f'weighted RMS {np.sqrt(np.mean(r**2)):.4f} Pa, raw RMS {prob.raw_rms(prob.loops[k]) if k in prob.loops else float("nan"):.4f} Pa')
        return

    t0 = time.time()
    res = least_squares(prob.fun, z0, jac=prob.jac, bounds=(0.0, 1.0), max_nfev=a.max_nfev, verbose=1, **LSQ)
    p_opt = prob.phys(res.x)

    # linearised uncertainty at the optimum (normalised -> physical), from the last Jacobian
    J = res.jac
    dof = max(len(res.fun) - len(free), 1)
    s2 = 2 * res.cost / dof
    try:
        cov_z = s2 * np.linalg.inv(J.T @ J)
        scale = prob.hi - prob.lo
        cov = cov_z * np.outer(scale, scale)
        se = np.sqrt(np.diag(cov))
        corr = cov / np.outer(se, se)
    except np.linalg.LinAlgError:
        se, corr = np.full(len(free), np.nan), np.full((len(free), len(free)), np.nan)

    key = prob.key(p_opt)
    loop = prob.loops.get(key)
    if loop is None:                                                     # optimum from the cache of an earlier session
        _, loop, _ = run_one((p_opt, cfg, data))
    np.savez(os.path.join(a.out, 'best_loop.npz'), gamma=loop[0], sigma=loop[1])
    out = {'status': res.status, 'message': res.message, 'nfev': res.nfev, 'njev': res.njev, 'laos_runs': prob.n_runs,
           'wall_h': (time.time() - t0) / 3600, 'cfg': cfg, 'free': free, 'fixed': fixed,
           'optimum': p_opt, 'a_ratio': p_opt['s_taut'] * p_opt['lambda0'],
           'k_off_0 (slip)': (1 - p_opt['phi']) * p_opt['koff0'], 'k_catch_0': p_opt['phi'] * p_opt['koff0'],
           'k_on_0': KEQ * p_opt['koff0'],
           'rms_weighted': float(np.sqrt(np.mean(res.fun**2))), 'rms_raw': prob.raw_rms(loop),
           'std_error': dict(zip(free, se.tolist())), 'correlation': corr.tolist(),
           'at_bound': [k for k, z in zip(free, res.x) if z < 1e-3 or z > 1 - 1e-3]}
    json.dump(out, open(os.path.join(a.out, 'result.json'), 'w'), indent=2)

    print(f'\n{res.message}  (nfev {res.nfev}, njev {res.njev}, LAOS runs {prob.n_runs}, {out["wall_h"]:.2f} h)')
    print(f'RMS weighted {out["rms_weighted"]:.4f} Pa, raw {out["rms_raw"]:.4f} Pa')
    for k, e in zip(free, se):
        print(f'  {k:8s} = {p_opt[k]:.5g}  +- {e:.2g}' + ('   (at bound)' if k in out['at_bound'] else ''))
    print(f'  derived: a_ratio {out["a_ratio"]:.5g}, k_off_0 {out["k_off_0 (slip)"]:.5g}, k_catch_0 {out["k_catch_0"]:.5g}, k_on_0 {out["k_on_0"]:.5g}')
    print('correlation matrix:\n' + '\n'.join(f'  {k:8s} ' + ' '.join(f'{c:+.2f}' for c in row) for k, row in zip(free, corr)))
    prob.pool.shutdown()


if __name__ == '__main__':
    main()
