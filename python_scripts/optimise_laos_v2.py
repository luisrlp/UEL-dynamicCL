"""Least-squares fit of the stretch-free-binding model to the LAOS loop of Unterberger et al. (2013), Fig. 4a.
Version 2 of optimise_laos.py: per-cycle diagnostics, candidate parameters keq / beta / mu0stretch, optional
early-cycle overshoot penalty (Step 2a), steady-state fits with --cycles N (Step 2b) and an optional joint fit with
the amplitude sweep of Fig. 4b (--sweep-weight). With the defaults (keq, beta, mu0stretch fixed at x0,
--overshoot-weight 0, --sweep-weight 0, --cycles 3) the residuals equal those of optimise_laos.py.

Model : python_scripts/material.py (DynamicCrosslinking), simple shear at one material point, closed system,
        gamma(t) = gamma_A sin(2 pi f t), f = 0.05 Hz, started from rest; the last simulated cycle is compared
        with the data (--cycles 3: the measured 3rd cycle; --cycles N > 3: a steady-state cycle).
Data  : python_scripts/validation_data/unterberger2013_fig4a_laos_dataset0.csv (gamma, sigma_xy, branch).
Method: scipy.optimize.least_squares (trust-region reflective, bounds) in normalised parameters z in [0, 1];
        the forward-difference Jacobian is evaluated in parallel (one LAOS run per free parameter).

Parameters (PARAMS; free by default unless listed in DEFAULT_FIXED):
  lambda0     prestretch LAMBDA0
  s_taut      a_ratio / LAMBDA0 (affine stretch at which a filament becomes taut; decorrelates LAMBDA0 and a_ratio)
  koff0       unloaded off-rate k_s + k_c [1/s]
  phi         catch fraction k_c / (k_s + k_c)
  dxc_nm      catch reactive distance [nm]
  dxs_nm      slip reactive distance [nm]
  keq         k_on_0 / koff0                                   (default fixed; release with --free keq)
  beta        effective extensional number BETA                (default fixed; release with --free beta)
  mu0stretch  stretch modulus MU0STRETCH [pN]                  (default fixed; release with --free mu0stretch)
Any parameter can be fixed with --fix name=value. All other material constants: FIXED_DIFF / FIXED_MECH.

Residuals: model - data at the data points of both branches, weighted for the extraction uncertainty
(delta_sigma, delta_gamma): r_w = r * delta_sigma / sqrt(delta_sigma^2 + (dsigma/dgamma * delta_gamma)^2),
which down-weights the near-vertical tips. Units stay Pa (r_w = r where the loop is flat).
Optional overshoot penalty (--overshoot-weight w > 0): 2 extra rows per early cycle k = 1, 2,
w_eff * max(0, sigma_max,k - sigma_max,last) and w_eff * max(0, sigma_min,last - sigma_min,k),
w_eff = w * sqrt(N_data / N_pen), so only an overshoot of the early cycles beyond the fitted cycle is penalised.
Optional amplitude sweep (--sweep-weight w > 0): one extra LAOS run from rest per amplitude of --sweep-gammas (same dt,
factor and cycles; G', G'' of the last cycle from main.moduli, Eq. (1) of the paper), and 2 rows per amplitude,
w_eff * log10(G'_model / G'_data) and w_eff * log10(G''_model / G''_data), w_eff = w * sqrt(N_data / N_sweep).
Data: validation_data/4b_G_.csv, 4b_G__.csv, interpolated in log-log at the requested amplitudes. With w = 1 an RMS
log10 error of 0.1 (26 %) costs as much as an RMS loop misfit of 0.1 Pa. A failed sweep run (or G' or G'' <= 0)
gives a log10 error of SWEEP_FAIL.

Usage (from any directory; results go to --out, default python_scripts/optimisation_runs/laos_v2):
  python optimise_laos_v2.py --evaluate-only                  # one LAOS run at x0, residual summary + cycle table
  python optimise_laos_v2.py --workers 7 --max-nfev 20        # fit
  python optimise_laos_v2.py --fix dxs_nm=0.5 --workers 6     # fit with dxs fixed
  python optimise_laos_v2.py --free keq --workers 8           # also fit keq
  python optimise_laos_v2.py --overshoot-weight 1             # Step 2a: penalise the first-cycle overshoot
  python optimise_laos_v2.py --cycles 6                       # Step 2b: fit a (near) steady-state cycle
  python optimise_laos_v2.py --sweep-weight 1 --workers 16    # joint fit: LAOS loop + amplitude sweep (Fig. 4b)
  python optimise_laos_v2.py --sweep-weight 1 --sweep-gammas 0.0126 0.1 0.2   # other sweep amplitudes
Outputs (in --out): evals.csv (every LAOS run, with per-cycle peaks and cycle-to-cycle change, and the sweep moduli),
                    cache.json (residuals per parameter set, reused on restart
                    only with the same settings, fixed constants, data and model source),
                    result.json (optimum, RMS, correlation matrix, standard errors, cycle diagnostics), best_loop.npz,
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
SWEEP_DATA = (os.path.join(REPO, 'validation_data', '4b_G_.csv'),     # Fig. 4b: G' and G'' vs gamma_A
              os.path.join(REPO, 'validation_data', '4b_G__.csv'))
sys.path.insert(0, REPO)

THETA = 294.0
KT = 1.38065e-5 * THETA            # pN um
FREQ = 0.05                        # Hz

# name: (x0, lower, upper)
PARAMS = {
    'lambda0': (1.0172, 1.000, 1.05), # 1.0467
    's_taut':  (1.072, 1.050, 1.200), # 1.0861 a_ratio/lambda0
    'koff0':   (0.279,  0.020, 2.00), # [0.08] k_c + k_s
    'phi':     (0.99876, 0.000, 0.999), # [0.9091] k_c / (k_s + k_c) = k_c / koff0
    'dxc_nm':  (15.0,    0.5,   35.0), # 2.
    'dxs_nm':  (1.66,    0.5,   20.0), # 0.5
    'keq':        (10.0,     3.0,    50.0),    # k_on_0 / koff0
    'beta':       (0.438,   0.30,   0.70),    # BETA (Unterberger 2013, from Liu & Pollack)
    'mu0stretch': (38600.0, 2.0e4,  6.0e4),   # MU0STRETCH [pN]
}
DEFAULT_FIXED = ('keq', 'beta', 'mu0stretch')     # fixed at their x0 unless released with --free

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
              'Lp': 16.0,
              'R0C': 0.014,
              'ETA': 1.,
              'b_disp': 1e-6,
              'prefdir': np.array([1.0, 0.0, 0.0]),
              'stretch_free_binding': True,
              'xi_t': 1.0}

DSIG, DGAM = 0.025, 0. # 0.0007         # extraction uncertainty (Pa, -)
STEADY_TOL = 0.01                  # Pa: RMS change of the last cycle vs the previous one above which it is not steady
N_EARLY = 2                        # early cycles checked by the overshoot penalty
SWEEP_GAMMAS = (0.0126, 0.05, 0.1, 0.158, 0.2)   # default sweep amplitudes (0.126 is the fitted loop itself)
SWEEP_FAIL = 2.0                   # log10 error assigned to a failed sweep run
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


def n_data(data):
    return sum(len(data[b][0]) for b in (1, -1))


def n_penalty(cfg):
    """Number of overshoot-penalty rows (0 when the penalty is off)."""
    return 2 * min(N_EARLY, cfg['cycles'] - 1) if cfg['overshoot_weight'] > 0 else 0


def n_sweep(cfg):
    """Number of amplitude-sweep rows (0 when the sweep is off)."""
    return 2 * len(cfg['sweep_gammas']) if cfg.get('sweep_weight', 0.0) > 0 else 0


def load_sweep_data(gammas):
    """G' and G'' of Fig. 4b at the amplitudes gammas (log-log interpolation); array (n_gammas, 2)."""
    out = []
    for f in SWEEP_DATA:
        d = np.loadtxt(f, delimiter=',', comments='#')
        d = d[np.argsort(d[:, 0])]                    # 4b_G__.csv is not sorted
        out.append(10.0**np.interp(np.log10(gammas), np.log10(d[:, 0]), np.log10(d[:, 1])))
    return np.column_stack(out)


def sweep_range():
    """Amplitude range covered by both sweep data files."""
    g = [np.loadtxt(f, delimiter=',', comments='#')[:, 0] for f in SWEEP_DATA]
    return max(x.min() for x in g), min(x.max() for x in g)


# ---------------------------------------------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------------------------------------------
def material_props(p):
    koff0, phi = p['koff0'], p['phi']
    diff = dict(FIXED_DIFF, k_off_0=(1.0 - phi) * koff0, k_catch_0=phi * koff0, k_on_0=p['keq'] * koff0,
                dx_kT=p['dxs_nm'] * 1e-3 / KT, dx_c_kT=p['dxc_nm'] * 1e-3 / KT)
    mech = dict(FIXED_MECH, LAMBDA0=p['lambda0'], a_ratio=p['s_taut'] * p['lambda0'],
                BETA=p['beta'], MU0STRETCH=p['mu0stretch'])
    return diff, mech


def simulate(p, cfg):
    """Strain-controlled LAOS from rest; returns time, gamma and sigma_xy of all cycles."""
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
    return t, g, s


def split_cycles(t, g, s, ncyc):
    """Per-cycle (gamma, sigma); cycle k covers [(k-1) T, k T] including both end points."""
    period = 1.0 / FREQ
    out = []
    for k in range(1, ncyc + 1):
        c = (t >= (k - 1) * period - 1e-9) & (t <= k * period + 1e-9)
        out.append((g[c], s[c]))
    return out


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


def cycle_diagnostics(cycles, data):
    """Peaks of every cycle and the RMS change of each cycle vs the previous one at the data strains [Pa]."""
    unweighted = {b: (g, s, np.ones_like(g)) for b, (g, s, w) in data.items()}
    on_data = [residuals_from_loop(gm, sm, unweighted) for gm, sm in cycles]     # model - data, unit weights
    dprev = [float('nan')] + [float(np.sqrt(np.mean((on_data[k] - on_data[k - 1])**2))) for k in range(1, len(cycles))]
    return {'smax': [float(sm.max()) for gm, sm in cycles], 'smin': [float(sm.min()) for gm, sm in cycles],
            'dprev': dprev}


def overshoot_rows(diag, cfg, n_dat):
    """Penalty rows: overshoot of the early cycles' peaks beyond the peaks of the fitted (last) cycle."""
    npen = n_penalty(cfg)
    if npen == 0:
        return np.array([])
    w_eff = cfg['overshoot_weight'] * math.sqrt(n_dat / npen)
    smax_l, smin_l = diag['smax'][-1], diag['smin'][-1]
    rows = []
    for k in range(min(N_EARLY, cfg['cycles'] - 1)):
        rows += [max(0.0, diag['smax'][k] - smax_l),          # positive peak above the fitted cycle's
                 max(0.0, smin_l - diag['smin'][k])]          # negative peak below the fitted cycle's
    return w_eff * np.array(rows)


def run_one(args):
    """Worker: one LAOS run for the physical parameter dict p.
    Returns (residual vector, last-cycle loop, cycle diagnostics, wall time)."""
    p, cfg, data = args
    t0 = time.time()
    try:
        t, g, s = simulate(p, cfg)
        cycles = split_cycles(t, g, s, cfg['cycles'])
        gm, sm = cycles[-1]
        diag = cycle_diagnostics(cycles, data)
        r = np.concatenate([residuals_from_loop(gm, sm, data), overshoot_rows(diag, cfg, n_data(data))])
        if not np.all(np.isfinite(r)):
            raise FloatingPointError('non-finite stress')
    except Exception as e:                           # failed run: large, finite penalty
        print(f'   ! run failed at {p}: {e!r}', flush=True)
        gm = sm = np.array([])
        diag = None
        r = np.full(n_data(data) + n_penalty(cfg), 10.0)
    return r, (gm, sm), diag, time.time() - t0


def run_sweep(args):
    """Worker: one LAOS run from rest at the sweep amplitude gA (same protocol as the fitted loop).
    Returns (G', G'' of the last cycle, wall time); NaN if the run fails."""
    from main import moduli
    p, cfg, gA = args
    t0 = time.time()
    try:
        t, g, s = simulate(p, dict(cfg, gamma_A=gA))
        last = t >= t[-1] - 1.0 / FREQ - 0.5 * cfg['dt']
        G1, G2, _ = moduli(t[last], g[last], s[last], FREQ)
    except Exception as e:
        print(f'   ! sweep run failed at gamma_A {gA} for {p}: {e!r}', flush=True)
        G1 = G2 = float('nan')
    return G1, G2, time.time() - t0


def sweep_errors(G, Gd):
    """log10(model / data), G' rows first, then G''; SWEEP_FAIL where the model modulus is NaN or <= 0."""
    G, Gd = G.T.ravel(), Gd.T.ravel()
    ok = np.isfinite(G) & (G > 0)
    return np.where(ok, np.log10(np.where(ok, G, 1.0) / Gd), SWEEP_FAIL)


def sweep_weight_eff(cfg, n_dat):
    return cfg['sweep_weight'] * math.sqrt(n_dat / n_sweep(cfg))


# ---------------------------------------------------------------------------------------------------------------
# Optimiser wrapper: normalised variables, cache, parallel forward-difference Jacobian
# ---------------------------------------------------------------------------------------------------------------
def settings_signature(cfg):
    """Hash of everything besides PARAMS that changes the residuals: run settings, fixed constants, weights,
    data file and model source. It is part of every cache key, so residuals computed with other settings
    (or before an edit of material.py) are never reused."""
    h = hashlib.sha1()
    fixed_mech = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in FIXED_MECH.items()}
    h.update(json.dumps({'cfg': {k: ([float(x) for x in v] if isinstance(v, (list, tuple)) else float(v))
                                 for k, v in cfg.items()},
                         'diff': FIXED_DIFF, 'mech': fixed_mech,
                         'FREQ': FREQ, 'DSIG': DSIG, 'DGAM': DGAM, 'N_EARLY': N_EARLY}, sort_keys=True).encode())
    for f in ('material.py', 'kinematics.py', 'dfd.py'):
        h.update(open(os.path.join(REPO, f), 'rb').read())
    h.update(open(DATA, 'rb').read())
    if n_sweep(cfg):                                  # sweep settings only when it is on (signature unchanged otherwise)
        import inspect
        from main import moduli
        h.update(json.dumps({'SWEEP_FAIL': SWEEP_FAIL}).encode())
        h.update(inspect.getsource(moduli).encode())
        for f in SWEEP_DATA:
            h.update(open(f, 'rb').read())
    return h.hexdigest()[:12]


class Problem:
    def __init__(self, free, fixed, cfg, data, out, workers, h):
        self.free, self.fixed, self.cfg, self.data, self.out, self.h = free, fixed, cfg, data, out, h
        self.n_dat = n_data(data)
        self.n_pen, self.n_sw = n_penalty(cfg), n_sweep(cfg)
        self.sweep_g = list(cfg['sweep_gammas']) if self.n_sw else []
        self.sweep_d = load_sweep_data(self.sweep_g) if self.n_sw else None      # (n_gammas, 2): G', G'' data
        self.lo = np.array([PARAMS[k][1] for k in free])
        self.hi = np.array([PARAMS[k][2] for k in free])
        self.pool = ProcessPoolExecutor(max_workers=workers)
        self.cache_file = os.path.join(out, 'cache.json')
        self.cache = json.load(open(self.cache_file)) if os.path.exists(self.cache_file) else {}
        self.signature = settings_signature(cfg)
        self.loops = {}
        self.diags = {}
        self.n_runs = 0
        self.log = os.path.join(out, 'evals.csv')
        ncyc = cfg['cycles']
        header = ','.join(['run'] + list(PARAMS) + ['rms_weighted', 'rms_penalty', 'rms_raw', 'last_cycle_change'] +
                          [f'{q}_c{k}' for k in range(1, ncyc + 1) for q in ('smax', 'smin', 'dprev')] + ['wall_s'] +
                          (['rms_sweep_log10'] + [f'{q}_g{g:.4g}' for q in ('G1', 'G2') for g in self.sweep_g] +
                           ['sweep_wall_s'] if self.n_sw else []))
        if not os.path.exists(self.log):
            with open(self.log, 'w') as fh:
                fh.write(header + '\n')
        elif open(self.log).readline().strip() != header:
            sys.exit(f'{self.log} has different columns (other parameters or number of cycles): use a new --out')

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
            # all LAOS runs and sweep runs of all points in one pool (a point needs 1 + n_gammas runs)
            laos = [self.pool.submit(run_one, (p, self.cfg, self.data)) for p in todo]
            sweep = [[self.pool.submit(run_sweep, (p, self.cfg, g)) for g in self.sweep_g] for p in todo]
            for p, fl, fs in zip(todo, laos, sweep):
                r, loop, diag, wall = fl.result()
                sw = np.array([f.result() for f in fs]).reshape(-1, 3)          # rows: G', G'', wall per amplitude
                if self.n_sw:
                    e = sweep_errors(sw[:, :2], self.sweep_d)
                    r = np.concatenate([r, sweep_weight_eff(self.cfg, self.n_dat) * e])
                k = self.key(p)
                self.cache[k] = r.tolist()
                self.loops[k] = loop
                self.diags[k] = diag
                self.n_runs += 1
                rd, rp, _ = self.split(r)
                ncyc = self.cfg['cycles']
                cyc = ([f'{diag[q][c]:.5g}' for c in range(ncyc) for q in ('smax', 'smin', 'dprev')] if diag
                       else ['nan'] * 3 * ncyc)
                last = diag['dprev'][-1] if diag else float('nan')
                swcols = ([f'{np.sqrt(np.mean(e**2)):.5f}'] + [f'{x:.5g}' for x in sw[:, 0]] +
                          [f'{x:.5g}' for x in sw[:, 1]] + [f'{sw[:, 2].sum():.0f}'] if self.n_sw else [])
                with open(self.log, 'a') as fh:
                    fh.write(','.join([str(self.n_runs)] + [f'{p[q]:.6g}' for q in PARAMS] +
                                      [f'{np.sqrt(np.mean(rd**2)):.5f}',
                                       f'{np.sqrt(np.mean(rp**2)) if len(rp) else 0.0:.5f}',
                                       f'{self.raw_rms(loop):.5f}', f'{last:.5f}'] + cyc + [f'{wall:.0f}'] +
                                      swcols) + '\n')
            json.dump(self.cache, open(self.cache_file, 'w'))
        return [np.array(self.cache[self.key(p)]) for p in ps]

    def split(self, r):
        """Residual blocks: loop data, overshoot penalty, amplitude sweep."""
        n1 = self.n_dat + self.n_pen
        return r[:self.n_dat], r[self.n_dat:n1], r[n1:]

    def sweep_moduli(self, r):
        """Model G', G'' (n_gammas, 2) and log10 errors recovered from the sweep rows of a residual vector
        (a failed run shows as G = data * 10**SWEEP_FAIL)."""
        e = self.split(r)[2] / sweep_weight_eff(self.cfg, self.n_dat)
        return self.sweep_d * 10.0**e.reshape(2, -1).T, e

    def raw_rms(self, loop):
        gm, sm = loop
        if len(gm) == 0:
            return float('nan')
        # data points with unit weights
        unweighted = {b: (g, s, np.ones_like(g)) for b, (g, s, w) in self.data.items()}
        return float(np.sqrt(np.mean(residuals_from_loop(gm, sm, unweighted)**2)))

    def fun(self, z):
        r = self.evaluate([z])[0]
        rd, rp, rs = self.split(r)
        pen = f', penalty RMS {np.sqrt(np.mean(rp**2)):.4f}' if len(rp) else ''
        pen += f', sweep RMS log10 {np.sqrt(np.mean(self.sweep_moduli(r)[1]**2)):.4f}' if len(rs) else ''
        print(f'   f: {self.fmt(z)} -> RMS_w {np.sqrt(np.mean(rd**2)):.4f} Pa{pen}  (points so far: {self.n_runs})', flush=True)
        diag = self.diags.get(self.key(self.phys(z)))
        if diag and self.cfg['cycles'] > 1 and diag['dprev'][-1] > STEADY_TOL:
            print(f'   ! fitted cycle not steady: RMS change vs previous cycle {diag["dprev"][-1]:.3f} Pa > {STEADY_TOL} Pa', flush=True)
        return r

    def jac(self, z):
        # forward differences in normalised variables; step backwards at the upper bound
        steps = [self.h if zi + self.h <= 1.0 else -self.h for zi in z]
        pts = [z] + [z + s * e for s, e in zip(steps, np.eye(len(z)))]
        rs = self.evaluate(pts)
        return np.column_stack([(rs[i + 1] - rs[0]) / steps[i] for i in range(len(z))])

    def fmt(self, z):
        # Free parameter formatting
        p = self.phys(z)
        return ', '.join(f'{k} {p[k]:.4g}' for k in self.free)


def print_cycles(diag):
    print('   cycle | sigma_max  sigma_min | RMS change vs previous cycle [Pa]')
    for c in range(len(diag['smax'])):
        print(f'   {c + 1:5d} | {diag["smax"][c]:9.3f} {diag["smin"][c]:9.3f} | {diag["dprev"][c]:.4f}')


def print_sweep(prob, r):
    G, e = prob.sweep_moduli(r)
    print(f"   gamma_A |  G' model  G' data | G'' model G'' data   [Pa]   (RMS log10 error {np.sqrt(np.mean(e**2)):.4f})")
    for g, (G1, G2), (d1, d2) in zip(prob.sweep_g, G, prob.sweep_d):
        print(f'   {g:7.4f} | {G1:9.3f} {d1:8.3f} | {G2:9.3f} {d2:9.3f}')


# ---------------------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default=os.path.join(REPO, 'optimisation_runs', 'laos_v2'))
    ap.add_argument('--fix', action='append', default=[], help='name=value, e.g. dxs_nm=0.5 (repeatable)')
    ap.add_argument('--free', action='append', default=[], help=f'release a default-fixed parameter {DEFAULT_FIXED} (repeatable)')
    ap.add_argument('--x0', action='append', default=[], help='name=value to override a starting value')
    ap.add_argument('--dt', type=float, default=0.1)
    ap.add_argument('--factor', type=int, default=3)
    ap.add_argument('--cycles', type=int, default=3, help='simulated cycles from rest; the last one is fitted')
    ap.add_argument('--overshoot-weight', type=float, default=0.0, help='weight of the early-cycle overshoot penalty (0 = off)')
    ap.add_argument('--sweep-weight', type=float, default=0.0, help='weight of the amplitude-sweep rows (Fig. 4b, 0 = off)')
    ap.add_argument('--sweep-gammas', type=float, nargs='+', default=list(SWEEP_GAMMAS),
                    help='sweep amplitudes (one extra LAOS run each)')
    ap.add_argument('--workers', type=int, default=7)
    ap.add_argument('--max-nfev', type=int, default=25, help='max residual evaluations of least_squares (Jacobians not counted)')
    ap.add_argument('--h', type=float, default=0.01, help='finite-difference step in normalised variables (fraction of the bound range)')
    ap.add_argument('--evaluate-only', action='store_true')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    data, gamma_A = load_data()
    cfg = {'dt': a.dt, 'factor': a.factor, 'cycles': a.cycles, 'gamma_A': gamma_A, 'overshoot_weight': a.overshoot_weight}
    if a.sweep_weight > 0:                           # only when on, so the cache signature of loop-only fits is unchanged
        lo, hi = sweep_range()
        if not all(lo <= g <= hi for g in a.sweep_gammas):
            sys.exit(f'--sweep-gammas outside the data range [{lo:.4g}, {hi:.4g}]')
        cfg.update(sweep_weight=a.sweep_weight, sweep_gammas=sorted(a.sweep_gammas))
    x0 = {k: PARAMS[k][0] for k in PARAMS}
    x0.update({k: float(v) for k, v in (s.split('=') for s in a.x0)})
    for k in a.free:
        if k not in PARAMS:
            sys.exit(f'--free {k}: unknown parameter')
    user_fixed = {k: float(v) for k, v in (s.split('=') for s in a.fix)}
    if set(a.free) & set(user_fixed):
        sys.exit(f'parameters both freed and fixed: {set(a.free) & set(user_fixed)}')
    fixed = {k: x0[k] for k in DEFAULT_FIXED if k not in a.free}
    fixed.update(user_fixed)
    free = [k for k in PARAMS if k not in fixed]
    for k in free:
        if not PARAMS[k][1] <= x0[k] <= PARAMS[k][2]:
            sys.exit(f'x0 of {k} = {x0[k]} outside bounds {PARAMS[k][1:]}')
    if a.overshoot_weight > 0 and a.cycles < 2:
        sys.exit('--overshoot-weight needs --cycles >= 2')
    print(f'gamma_A = {gamma_A:.5f}, cfg = {cfg}\nfree = {free}\nfixed = {fixed}', flush=True)

    prob = Problem(free, {**{k: x0[k] for k in PARAMS}, **fixed}, cfg, data, a.out, a.workers, a.h)
    z0 = np.array([(x0[k] - PARAMS[k][1]) / (PARAMS[k][2] - PARAMS[k][1]) for k in free])

    # full setup of this session, written before any LAOS run (so it also exists for interrupted fits)
    started = time.strftime('%Y-%m-%d_%H%M%S')
    setup = {'started': started, 'command': ' '.join(sys.argv), 'signature': prob.signature,
             'data_file': DATA, 'cfg': cfg, 'free': free, 'fixed': fixed, 'default_fixed': list(DEFAULT_FIXED),
             'x0': x0, 'bounds': {k: [PARAMS[k][1], PARAMS[k][2]] for k in PARAMS},
             'constants': {'FIXED_DIFF': FIXED_DIFF,
                           'FIXED_MECH': {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in FIXED_MECH.items()},
                           'FREQ': FREQ, 'THETA': THETA},
             'weights': {'DSIG': DSIG, 'DGAM': DGAM, 'N_EARLY': N_EARLY, 'n_penalty_rows': n_penalty(cfg),
                         'n_sweep_rows': n_sweep(cfg), 'SWEEP_FAIL': SWEEP_FAIL},
             'sweep_data_files': list(SWEEP_DATA) if n_sweep(cfg) else [],
             'steady_tol': STEADY_TOL,
             'optimiser': {**LSQ, 'max_nfev': a.max_nfev, 'h': a.h, 'workers': a.workers, 'evaluate_only': a.evaluate_only},
             'versions': {'numpy': np.__version__, 'scipy': scipy.__version__, 'python': sys.version.split()[0]}}
    for name in ('setup.json', f'setup_{started}.json'):
        json.dump(setup, open(os.path.join(a.out, name), 'w'), indent=2)

    if a.evaluate_only:
        r = prob.fun(z0)
        k = prob.key(prob.phys(z0))
        print(f'weighted RMS (data) {np.sqrt(np.mean(r[:prob.n_dat]**2)):.4f} Pa, '
              f'raw RMS {prob.raw_rms(prob.loops[k]) if k in prob.loops else float("nan"):.4f} Pa')
        if prob.diags.get(k):
            print_cycles(prob.diags[k])
        if prob.n_sw:
            print_sweep(prob, r)
        prob.pool.shutdown()
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
    loop, diag = prob.loops.get(key), prob.diags.get(key)
    if loop is None or diag is None:                                     # optimum from the cache of an earlier session
        _, loop, diag, _ = run_one((p_opt, cfg, data))
    np.savez(os.path.join(a.out, 'best_loop.npz'), gamma=loop[0], sigma=loop[1])
    rd, rp, rs = prob.split(res.fun)
    if len(rs):
        G, e = prob.sweep_moduli(res.fun)
        sweep = {'gamma_A': prob.sweep_g, 'G1_model': G[:, 0].tolist(), 'G1_data': prob.sweep_d[:, 0].tolist(),
                 'G2_model': G[:, 1].tolist(), 'G2_data': prob.sweep_d[:, 1].tolist(),
                 'rms_log10': float(np.sqrt(np.mean(e**2)))}
    out = {'status': res.status, 'message': res.message, 'nfev': res.nfev, 'njev': res.njev, 'laos_runs': prob.n_runs,
           'wall_h': (time.time() - t0) / 3600, 'cfg': cfg, 'free': free, 'fixed': fixed,
           'optimum': p_opt, 'a_ratio': p_opt['s_taut'] * p_opt['lambda0'],
           'k_off_0 (slip)': (1 - p_opt['phi']) * p_opt['koff0'], 'k_catch_0': p_opt['phi'] * p_opt['koff0'],
           'k_on_0': p_opt['keq'] * p_opt['koff0'],
           'rms_weighted': float(np.sqrt(np.mean(rd**2))),
           'rms_penalty': float(np.sqrt(np.mean(rp**2))) if len(rp) else 0.0,
           'rms_raw': prob.raw_rms(loop),
           'cycles': diag,
           **({'sweep_runs': prob.n_runs * len(prob.sweep_g), 'sweep': sweep} if len(rs) else {}),
           'std_error': dict(zip(free, se.tolist())), 'correlation': corr.tolist(),
           'at_bound': [k for k, z in zip(free, res.x) if z < 1e-3 or z > 1 - 1e-3]}
    json.dump(out, open(os.path.join(a.out, 'result.json'), 'w'), indent=2)

    print(f'\n{res.message}  (nfev {res.nfev}, njev {res.njev}, LAOS runs {prob.n_runs}, {out["wall_h"]:.2f} h)')
    print(f'RMS weighted {out["rms_weighted"]:.4f} Pa, raw {out["rms_raw"]:.4f} Pa' +
          (f', penalty {out["rms_penalty"]:.4f}' if len(rp) else ''))
    for k, e in zip(free, se):
        print(f'  {k:10s} = {p_opt[k]:.5g}  +- {e:.2g}' + ('   (at bound)' if k in out['at_bound'] else ''))
    print(f'  derived: a_ratio {out["a_ratio"]:.5g}, k_off_0 {out["k_off_0 (slip)"]:.5g}, k_catch_0 {out["k_catch_0"]:.5g}, k_on_0 {out["k_on_0"]:.5g}')
    print('correlation matrix:\n' + '\n'.join(f'  {k:10s} ' + ' '.join(f'{c:+.2f}' for c in row) for k, row in zip(free, corr)))
    if diag:
        print_cycles(diag)
    if len(rs):
        print_sweep(prob, res.fun)
    prob.pool.shutdown()


if __name__ == '__main__':
    main()
