"""
One-at-a-time sensitivity analysis (Part A, simple shear) with the material-point Fortran driver main_sa.f90.

Models : tc    = 10_CL_DIFF_CS_SFR_FI_TC (consistent kinetics)
         nontc = 9_CL_DIFF_CS_SFR_FI     (non-consistent kinetics), same parameter sets
Base   : optimisation_runs/final_F1p_dt025_f6 (properties through optimise_laos_v2.fortran_props, as the fit)
Cases  : the base and 15 parameters x 2 levels (LEVELS), each in 7 protocols (PROTOCOLS), closed and open system.

  python -B sensitivity_runs.py run                                  # everything (missing runs only)
  python -B sensitivity_runs.py run --models tc --protocols cyc_126  # a subset
  python -B sensitivity_runs.py run --factor 3 --out sensitivity_smoke --params koff0 chi   # smoke test
  python -B sensitivity_runs.py analyse                              # summary.csv, sensitivity.csv, figures

Executables: <model folder>/sa_point_f<factor> (bash build_sa.sh 10 6; bash build_sa.sh 9 6).
Results : <out>/<model>/<system>/<protocol>/<case>/ with run.npz (every increment), snap.npz (all directions at the
          snapshot times), stdout.txt (kinetics warnings of the UMAT) and meta.json (signature, wall time). A run is
          repeated only if its signature (properties, protocol, system, dt, executable) changed.
Analysis: <out>/summary.csv (one row per run, metrics()), <out>/sensitivity.csv (effect y/y0 - 1, or y - y0 for the
          metrics of ABSOLUTE_EFFECT, and elasticity effect/(q/q0 - 1) of each level, q the natural variable of the
          parameter), figures in
          figures/sensitivity/<out name>/ (tornado and overlay plots per model, system and protocol).
"""
import argparse, csv, hashlib, json, os, subprocess, sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import matplotlib
matplotlib.use('Agg')
import optimise_laos_v2 as opt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
BASE_RUN = 'final_F1p_dt025_f6'
MODELS = {'tc': '10_CL_DIFF_CS_SFR_FI_TC', 'nontc': '9_CL_DIFF_CS_SFR_FI'}
SYSTEMS = {'closed': 1, 'open': 0}
FREQ = 0.05
NCYC = 7

# Loading segments (type, v1, v2, v3): 1 ramp (target gamma, duration), 2 hold (duration), 3 sine (amplitude,
# frequency, cycles, about the gamma at its start); snapshots: times at which all directions are written.
# main_sa.f90 steps adaptively between breakpoints (segment ends, snapshot times, sine points every DT_SINE, gamma
# linear in between) under the PNEWDT of the UMAT (Backward-Euler error estimate of the kinetics), with steps of at
# most dt_max. Background (2026-11-07): with fixed dt = 0.025 the TC response in the strain-stiffening range was too
# soft (ramp in 1 s: -18 % stress at gamma 0.2; cyc_158: -1.4 % peak, -5 % for a = 1.1; 0.1 s step ramp: up to 5x).
DT_SINE = 0.025          # breakpoint spacing of the sines (linear interpolation error of gamma ~1e-5 relative)
DT_INIT = 1.0e-4         # first trial step
PROTOCOLS = {
    'ramp_slow': {'segments': [(1, 0.2, 100.0, 0)], 'snaps': [100.0], 'dt_max': 0.5},
    'ramp_mid':  {'segments': [(1, 0.2, 16.0, 0)], 'snaps': [16.0], 'dt_max': 0.08},
    'ramp_fast': {'segments': [(1, 0.2, 1.0, 0)], 'snaps': [1.0], 'dt_max': 0.005},
    'step':      {'segments': [(1, 0.126, 0.1, 0), (2, 200.0, 0, 0), (1, 0.0, 0.1, 0), (2, 200.0, 0, 0)],
                  'snaps': [0.1, 200.1, 200.2, 400.2], 'dt_max': 0.5},
    'cyc_lin':   {'segments': [(3, 0.0126, FREQ, NCYC)], 'snaps': [125.0, 130.0, 135.0, 140.0], 'dt_max': DT_SINE},
    'cyc_126':   {'segments': [(3, 0.126, FREQ, NCYC)], 'snaps': [125.0, 130.0, 135.0, 140.0], 'dt_max': DT_SINE},
    'cyc_158':   {'segments': [(3, 0.158, FREQ, NCYC)], 'snaps': [125.0, 130.0, 135.0, 140.0], 'dt_max': DT_SINE},
}

# Composition and mixing constants (not fitted): index in PROPS (0-based) and base value (optimise_laos_v2.FIXED_DIFF)
PROPS_INDEX = {'cactin': 14, 'R': 15, 'Rfmax': 16, 'Rbmax': 17, 'chi': 18}
PROPS_BASE = {'cactin': 'CACTIN', 'R': 'R', 'Rfmax': 'RFMAX', 'Rbmax': 'RBMAX', 'chi': 'CHI'}

# Levels of each parameter (functions of the base value b) and its natural variable q (for the elasticity)
LEVELS = {
    'lambda0':    (lambda b: (1.03, 1.07),          lambda v: v - 1.0),     # at fixed a
    'a':          (lambda b: (1.10, 1.20),          lambda v: v),           # at fixed lambda0 (s_taut = a/lambda0)
    'koff0':      (lambda b: (0.25, 1.0),           lambda v: v),           # at fixed keq
    'keq':        (lambda b: (2.0, 50.0),           lambda v: v),           # at fixed koff0
    'phi':        (lambda b: (0.5, 0.1),            lambda v: 1.0 - v),     # slip fraction 1 - phi
    'dxs_nm':     (lambda b: (b / 2.0, 2.0 * b),    lambda v: v),
    'dxc_nm':     (lambda b: (0.0, 0.3),            lambda v: v),
    'beta':       (lambda b: (b - 0.1, b + 0.1),    lambda v: v),
    'mu0stretch': (lambda b: (b / 1.5, 1.5 * b),    lambda v: v),
    'eta':        (lambda b: (0.8, 0.6),            lambda v: v),
    'cactin':     (lambda b: (b / 2.0, 2.0 * b),    lambda v: v),
    'R':          (lambda b: (b / 2.0, 2.0 * b),    lambda v: v),
    'Rbmax':      (lambda b: (b / 2.0, 2.0 * b),    lambda v: v),
    'Rfmax':      (lambda b: (b / 2.0, 2.0 * b),    lambda v: v),
    'chi':        (lambda b: (0.0, 1.0),            lambda v: v),
}

# Columns of <base>.dat (main_sa.f90)
COLS = ['t', 'gamma', 's11', 's22', 's33', 's12', 's13', 's23', 'thetaf', 'mu', 'cb_tot', 'content',
        'cb_mean', 'cb_std', 'cb_min', 'cb_max', 'lamn_mean', 'lamn_std', 'lamn_min', 'lamn_max',
        'lame_mean', 'lame_std', 'lame_min', 'lame_max', 'f_max', 'd_remod', 'd_react', 'p_mech', 'exchange',
        'it_mu', 'converged', 'pnewdt']
SNAP_COLS = ['isnap', 't', 'gamma', 'idir', 'm0x', 'm0y', 'm0z', 'w', 'cb', 'lambar', 'lamn', 'lame', 'f', 'g', 'x',
             'affinity', 'rate']
C = {c: i for i, c in enumerate(COLS)}


# ----------------------------------------------------------------------------------------------- cases and running
def base_params():
    r = json.load(open(os.path.join(HERE, 'optimisation_runs', BASE_RUN, 'result.json')))
    return dict(r['fixed'], **r['optimum'])


def base_value(name, p0):
    if name == 'a':
        return p0['s_taut'] * p0['lambda0']
    if name in PROPS_INDEX:
        return opt.FIXED_DIFF[PROPS_BASE[name]]
    return p0[name]


def case_list(params):
    """[(case name, parameter, level 0/1/2, value)]: the base and two levels per parameter."""
    p0 = base_params()
    cases = [('base', 'base', 0, float('nan'))]
    for name in params:
        for k, v in enumerate(LEVELS[name][0](base_value(name, p0)), start=1):
            cases.append((f'{name}_{v:.6g}', name, k, float(v)))
    return cases


def case_props(name, value):
    """PROPS(1:28) of a case: the base with one parameter changed."""
    p = base_params()
    a0 = p['s_taut'] * p['lambda0']
    if name == 'lambda0':                            # at fixed a = s_taut * lambda0
        p['lambda0'], p['s_taut'] = value, a0 / value
    elif name == 'a':                                # at fixed lambda0
        p['s_taut'] = value / p['lambda0']
    elif name in p and name != 'base':
        p[name] = value
    props = [float(x) for x in opt.fortran_props(p)]
    if name in PROPS_INDEX:
        props[PROPS_INDEX[name]] = value
    return props


def exe_path(model, factor):
    return os.path.join(REPO, MODELS[model], f'sa_point_f{factor}')


def check_exe(model, factor):
    """Stop if the driver is missing or older than its sources (driver, global.f90, material routines)."""
    exe = exe_path(model, factor)
    if not os.access(exe, os.X_OK):
        sys.exit(f'{exe} not found: bash build_sa.sh {MODELS[model].split("_")[0]} {factor}')
    fdir = os.path.join(REPO, MODELS[model])
    src = [os.path.join(REPO, MODELS['tc'], 'main_sa.f90'), os.path.join(fdir, 'global.f90')]
    mat = os.path.join(fdir, 'src', 'material_AFFCL')
    src += [os.path.join(mat, f) for f in os.listdir(mat) if f.endswith('.f90')]
    newest = max(src, key=os.path.getmtime)
    if os.path.getmtime(exe) < os.path.getmtime(newest):
        sys.exit(f'{exe} is older than {newest}: rebuild (bash build_sa.sh {MODELS[model].split("_")[0]} {factor})')


def breakpoints(protocol):
    """Breakpoints (t, gamma) of the loading (t = 0, gamma = 0 excluded): segment ends, sine points every DT_SINE,
    snapshot times."""
    pr = PROTOCOLS[protocol]
    t, g, ts, gs = 0.0, 0.0, [], []
    for typ, v1, v2, v3 in pr['segments']:
        if typ == 1:                                 # ramp to v1 in v2
            tau, gseg = np.array([v2]), np.array([v1])
        elif typ == 2:                               # hold for v1
            tau, gseg = np.array([v1]), np.array([g])
        else:                                        # sine: amplitude v1, frequency v2, v3 cycles
            n = int(round(v3 / v2 / DT_SINE))
            tau = (v3 / v2) * np.arange(1, n + 1) / n
            gseg = g + v1 * np.sin(2.0 * np.pi * v2 * tau)
        ts.append(t + tau)
        gs.append(gseg)
        t, g = t + tau[-1], gseg[-1]
    t, g = np.concatenate(ts), np.concatenate(gs)
    extra = [x for x in pr['snaps'] if np.min(np.abs(t - x)) > 1e-9]
    if extra:                                        # snapshot times inside a segment become breakpoints
        tt = np.sort(np.concatenate([t, extra]))
        g, t = np.interp(tt, np.r_[0.0, t], np.r_[0.0, g]), tt
    return t, g


def driver_input(props, system, dt_max, protocol, base):
    """Input of main_sa.f90; dt_max overrides the maximum step of the protocol (None: protocol value)."""
    pr = PROTOCOLS[protocol]
    t, g = breakpoints(protocol)
    lines = [' '.join(repr(v) for v in props),
             f'{SYSTEMS[system]} {len(t)} {len(pr["snaps"])} {dt_max or pr["dt_max"]!r} {DT_INIT!r}']
    lines += [f'{ti!r} {gi!r}' for ti, gi in zip(t.tolist(), g.tolist())]
    if pr['snaps']:
        lines.append(' '.join(repr(float(x)) for x in pr['snaps']))
    lines.append(base)
    return '\n'.join(lines) + '\n'


def signature(props, system, dt, protocol, exe):
    st = os.stat(exe)
    s = json.dumps({'props': props, 'system': system, 'dt': dt, 'protocol': PROTOCOLS[protocol],
                    'exe': [os.path.basename(exe), st.st_size, st.st_mtime]}, sort_keys=True)
    return hashlib.sha1(s.encode()).hexdigest()[:16]


def run_job(job):
    """Worker: one driver run; saves run.npz, snap.npz, stdout.txt, meta.json in the case folder."""
    folder, inp, sig, meta = job
    os.makedirs(folder, exist_ok=True)
    base = os.path.join(folder, 'raw')
    t0 = time.time()
    with open(os.path.join(folder, 'stdout.txt'), 'w') as out:
        proc = subprocess.run([meta['exe']], input=inp, text=True, stdout=out, stderr=subprocess.PIPE)
    wall = time.time() - t0
    if proc.returncode != 0 or not os.path.exists(base + '.dat'):
        return folder, False, wall, proc.stderr.strip()
    d = np.loadtxt(base + '.dat', comments='#')
    s = np.loadtxt(base + '_snap.dat', comments='#', ndmin=2)
    header = [l for l in open(base + '.dat') if l.startswith('# ')][:2]
    np.savez_compressed(os.path.join(folder, 'run.npz'), data=d, header=np.array(header))
    np.savez_compressed(os.path.join(folder, 'snap.npz'), data=s)
    os.remove(base + '.dat')
    os.remove(base + '_snap.dat')
    nwarn = sum('WARNING' in l for l in open(os.path.join(folder, 'stdout.txt')))
    json.dump(dict(meta, signature=sig, wall_s=wall, warnings=nwarn, stderr=proc.stderr.strip()),
              open(os.path.join(folder, 'meta.json'), 'w'), indent=1)
    return folder, True, wall, proc.stderr.strip()


def cmd_run(a):
    params = a.params or list(LEVELS)
    for m in a.models:
        check_exe(m, a.factor)
    cost = {pr: breakpoints(pr)[0][-1] * (5.0 if pr.startswith('ramp_fast') else 1.0) for pr in PROTOCOLS}
    jobs = []
    for m in a.models:
        exe = exe_path(m, a.factor)
        for s in a.systems:
            for pr in a.protocols:
                for case, name, lev, val in case_list(params):
                    folder = os.path.join(a.out, m, s, pr, case)
                    props = case_props(name, val)
                    sig = signature(props, s, a.dt, pr, exe)
                    mf = os.path.join(folder, 'meta.json')
                    if os.path.exists(mf) and json.load(open(mf)).get('signature') == sig and not a.force:
                        continue
                    inp = driver_input(props, s, a.dt, pr, os.path.join(folder, 'raw'))
                    meta = {'model': m, 'system': s, 'protocol': pr, 'case': case, 'param': name, 'level': lev,
                            'value': val, 'props': props, 'dt': a.dt, 'factor': a.factor, 'exe': exe}
                    jobs.append((cost[pr] * (2.5 if s == 'closed' else 1.0), (folder, inp, sig, meta)))
    jobs.sort(key=lambda j: -j[0])                   # longest first
    print(f'{len(jobs)} runs to do ({a.workers} workers)', flush=True)
    t0, done = time.time(), 0
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for folder, ok, wall, err in pool.map(run_job, [j[1] for j in jobs]):
            done += 1
            rel = os.path.relpath(folder, a.out)
            print(f'[{done}/{len(jobs)} {time.time() - t0:7.0f} s] {rel}: {"ok" if ok else "FAILED"} '
                  f'{wall:.0f} s' + (f'  ({err})' if err else ''), flush=True)


# ---------------------------------------------------------------------------------------------------------- metrics
def last_cycle(d, freq=FREQ, ncyc=NCYC, k=0):
    """Rows of cycle ncyc - k (closed interval, both ends)."""
    t, eps = d[:, C['t']], 1e-6 / freq
    t1 = (ncyc - k) / freq
    return d[(t >= t1 - 1.0 / freq - eps) & (t <= t1 + eps)]


def integral(rate, t):
    """Time integral of a per-increment rate (row j: rate over the increment ending at t_j; row 0 is the rest)."""
    return float(np.sum(rate[1:] * np.diff(t)))


# Metrics that are already changes or fractions (near zero, either sign): their effect is the difference y - y0,
# the others use the relative change y/y0 - 1
ABSOLUTE_EFFECT = {'cb_end', 'content_min', 'content_end', 'react_frac', 'd_react_min', 'overshoot', 'last_change',
                   'cb_mean_last', 'react_frac_last', 'thetaf_end', 'lamn_recovered', 'cb_hold', 'energy_ratio'}


def metrics(d, protocol, cb0, c0):
    """Scalar outputs of one run (time series d, variable time steps)."""
    from main import moduli
    t, g, s = d[:, C['t']], d[:, C['gamma']], d[:, C['s12']]
    eps = 1e-9
    dtot = d[:, C['d_remod']] + d[:, C['d_react']]
    diss = integral(dtot, t)
    m = {'cb_end': d[-1, C['cb_tot']] / cb0 - 1.0,
         'spread_max': d[:, C['cb_std']].max() / cb0,
         'lamn_spread_max': d[:, C['lamn_std']].max(),
         'lame_max': d[:, C['lame_max']].max(),
         'f_max': d[:, C['f_max']].max(),
         'diss_total': diss,
         'react_frac': integral(d[:, C['d_react']], t) / max(diss, 1e-300),
         'd_react_min': d[:, C['d_react']].min(),
         'content_min': d[:, C['content']].min() / c0 - 1.0,
         'content_end': d[-1, C['content']] / c0 - 1.0,
         'unconverged': int((d[:, C['converged']] == 0).sum()),
         'pnewdt_fail': int((d[:, C['pnewdt']] < 1.0).sum())}
    if protocol.startswith('cyc'):
        lc, pc = last_cycle(d), last_cycle(d, k=1)
        tl = lc[:, C['t']]
        tu = np.linspace(tl[0], tl[-1], 4001)       # moduli() assumes uniform steps: resample the cycle
        G1, G2, _ = moduli(tu, np.interp(tu, tl, lc[:, C['gamma']]), np.interp(tu, tl, lc[:, C['s12']]), FREQ)
        i = np.argmax(lc[:, C['gamma']])
        work = integral(lc[:, C['p_mech']], tl)
        dcyc = integral(lc[:, C['d_remod']] + lc[:, C['d_react']], tl)
        s1 = s[t <= 1.0 / FREQ + eps]                # first cycle
        smax = np.abs(lc[:, C['s12']]).max()
        m.update(G1=G1, G2=G2, GM=lc[i, C['s12']] / lc[i, C['gamma']], loop_work=work,
                 energy_ratio=dcyc / work if work else float('nan'),
                 overshoot=np.abs(s1).max() / smax - 1.0,            # largest |stress| of cycle 1 vs steady cycle
                 virgin_peak=s1.max() / lc[:, C['s12']].max(),       # first positive peak (virgin loading) vs steady
                 last_change=np.sqrt(np.mean((lc[:, C['s12']] - pc[:, C['s12']])**2)) / smax,
                 cb_mean_last=lc[:, C['cb_tot']].mean() / cb0 - 1.0,
                 react_frac_last=integral(lc[:, C['d_react']], tl) / max(dcyc, 1e-300))
    elif protocol.startswith('ramp'):
        K = np.gradient(s, g)
        m.update(s_005=np.interp(0.05, g, s), s_01=np.interp(0.1, g, s), s_02=s[-1],
                 stiffening=K[-1] / np.interp(0.01, g, K), thetaf_end=d[-1, C['thetaf']] / d[0, C['thetaf']] - 1.0)
    elif protocol == 'step':
        t0, t1 = 0.1, 200.1
        s0 = np.interp(t0, t, s)
        hold = (t >= t0 - eps) & (t <= t1 + eps)
        th, sh = t[hold], s[hold]
        below = np.where(sh <= 0.5 * s0)[0]
        if len(below) and below[0] > 0:              # linear interpolation of the half-relaxation time
            j = below[0]
            half_time = th[j - 1] + (0.5 * s0 - sh[j - 1]) * (th[j] - th[j - 1]) / (sh[j] - sh[j - 1]) - t0
        else:
            half_time = float('nan')
        std_hold = np.interp(t1, t, d[:, C['lamn_std']])
        m.update(s0=s0, s_peak=s.max(), relax_10=np.interp(t0 + 10.0, t, s) / s0,
                 relax_50=np.interp(t0 + 50.0, t, s) / s0, residual=np.interp(t1, t, s) / s0, half_time=half_time,
                 s_unload=np.interp(t1 + 0.1, t, s) / s0,
                 lamn_recovered=1.0 - d[-1, C['lamn_std']] / std_hold if std_hold > 0 else float('nan'),
                 cb_hold=np.interp(t1, t, d[:, C['cb_tot']]) / cb0 - 1.0)
    return m


def header_value(header, key):
    words = ' '.join(header).replace('#', ' ').split()
    return float(words[words.index(key) + 1])


def cmd_analyse(a):
    rows = []
    for root, dirs, files in os.walk(a.out):
        if 'meta.json' not in files or 'run.npz' not in files:
            continue
        meta = json.load(open(os.path.join(root, 'meta.json')))
        z = np.load(os.path.join(root, 'run.npz'))
        d, header = z['data'], list(z['header'])
        cb0 = header_value(header, 'cb0')
        c0 = d[0, C['content']]
        r = {k: meta[k] for k in ('model', 'system', 'protocol', 'case', 'param', 'level', 'value')}
        r.update(warnings=meta['warnings'], wall_s=meta['wall_s'], cb0=cb0, thetaf0=header_value(header, 'thetaf0'),
                 mu_rest=header_value(header, 'mu_rest'))
        r.update(metrics(d, meta['protocol'], cb0, c0))
        rows.append(r)
    if not rows:
        sys.exit(f'no runs in {a.out}')
    rows.sort(key=lambda r: (r['model'], r['system'], r['protocol'], r['param'] != 'base', r['param'], r['level']))
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(a.out, 'summary.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    # sensitivities relative to the base of the same model, system and protocol
    p0 = base_params()
    sens = []
    base = {(r['model'], r['system'], r['protocol']): r for r in rows if r['param'] == 'base'}
    info = ('model', 'system', 'protocol', 'case', 'param', 'level', 'value', 'warnings', 'wall_s', 'cb0',
            'thetaf0', 'mu_rest')
    for r in rows:
        b = base.get((r['model'], r['system'], r['protocol']))
        if r['param'] == 'base' or b is None:
            continue
        qf = LEVELS[r['param']][1]
        q0, q = qf(base_value(r['param'], p0)), qf(r['value'])
        for k, y in r.items():
            if k in info or not isinstance(y, (int, float)) or not isinstance(b.get(k), (int, float)):
                continue
            y0 = b[k]
            if k in ABSOLUTE_EFFECT:
                eff = y - y0
            else:
                eff = y / y0 - 1.0 if y0 not in (0, 0.0) and np.isfinite(y0) else float('nan')
            ela = eff / (q / q0 - 1.0) if q0 != 0 and q != q0 else float('nan')
            sens.append({'model': r['model'], 'system': r['system'], 'protocol': r['protocol'], 'metric': k,
                         'param': r['param'], 'level': r['level'], 'value': r['value'], 'y0': y0, 'y': y,
                         'effect': eff, 'elasticity': ela})
    with open(os.path.join(a.out, 'sensitivity.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(sens[0]))
        w.writeheader()
        w.writerows(sens)
    print(f'{len(rows)} runs -> {a.out}/summary.csv, {len(sens)} rows -> {a.out}/sensitivity.csv')
    if not a.no_plots:
        plots(a.out, rows, sens)


# ----------------------------------------------------------------------------------------------------------- plots
KEY_METRICS = {'cyc': ['G1', 'G2', 'GM', 'overshoot', 'cb_mean_last', 'spread_max'],
               'ramp': ['s_005', 's_02', 'stiffening', 'cb_end', 'spread_max', 'lame_max'],
               'step': ['s0', 'relax_10', 'half_time', 'residual', 's_unload', 'cb_hold']}


def plots(out, rows, sens):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['text.usetex'] = False
    fdir = os.path.join(REPO, 'figures', 'sensitivity', os.path.basename(os.path.normpath(out)))
    os.makedirs(fdir, exist_ok=True)
    groups = sorted({(r['model'], r['system'], r['protocol']) for r in rows})
    for m, s, pr in groups:
        kind = 'cyc' if pr.startswith('cyc') else ('ramp' if pr.startswith('ramp') else 'step')
        # tornado: effect of both levels per parameter, sorted by the largest effect
        sel = [x for x in sens if (x['model'], x['system'], x['protocol']) == (m, s, pr)]
        if sel:
            mets = KEY_METRICS[kind]
            fig, axs = plt.subplots(2, 3, figsize=(15, 9))
            for ax, met in zip(axs.flat, mets):
                xs = [x for x in sel if x['metric'] == met]
                par = sorted({x['param'] for x in xs},
                             key=lambda p: max(abs(x['effect']) if np.isfinite(x['effect']) else 0
                                               for x in xs if x['param'] == p))
                for lev, col in ((1, 'C0'), (2, 'C3')):
                    e = [next((100 * x['effect'] for x in xs if x['param'] == p and x['level'] == lev), np.nan)
                         for p in par]
                    ax.barh(np.arange(len(par)) + (0.2 if lev == 2 else -0.2), e, height=0.4, color=col,
                            label=f'level {lev}')
                ax.set_yticks(np.arange(len(par)))
                ax.set_yticklabels(par, fontsize=8)
                ax.axvline(0, color='k', lw=0.5)
                ax.set_title(met, fontsize=10)
                ax.set_xlabel('difference from base (x100)' if met in ABSOLUTE_EFFECT else 'change from base (%)',
                              fontsize=8)
            axs.flat[0].legend(fontsize=8)
            fig.suptitle(f'{m}, {s}, {pr}')
            fig.tight_layout()
            fig.savefig(os.path.join(fdir, f'tornado_{m}_{s}_{pr}.pdf'))
            plt.close(fig)
        # overlay: stress response of the base and both levels, one panel per parameter
        params = sorted({r['param'] for r in rows if (r['model'], r['system'], r['protocol']) == (m, s, pr)} - {'base'})
        if not params:
            continue
        nc = 5
        fig, axs = plt.subplots(int(np.ceil(len(params) / nc)), nc, figsize=(3.2 * nc, 2.6 * np.ceil(len(params) / nc)),
                                squeeze=False)
        bdir = os.path.join(out, m, s, pr, 'base')
        bd = np.load(os.path.join(bdir, 'run.npz'))['data'] if os.path.exists(os.path.join(bdir, 'run.npz')) else None
        for ax, par in zip(axs.flat, params):
            curves = [('base', bd, 'k')]
            for r in rows:
                if (r['model'], r['system'], r['protocol'], r['param']) == (m, s, pr, par):
                    z = np.load(os.path.join(out, m, s, pr, r['case'], 'run.npz'))['data']
                    curves.append((f'{r["value"]:.4g}', z, 'C0' if r['level'] == 1 else 'C3'))
            for lab, z, col in curves:
                if z is None:
                    continue
                if kind == 'cyc':
                    z = last_cycle(z)
                    ax.plot(z[:, C['gamma']], z[:, C['s12']], color=col, lw=1, label=lab)
                elif kind == 'ramp':
                    ax.plot(z[:, C['gamma']], z[:, C['s12']], color=col, lw=1, label=lab)
                else:
                    ax.plot(z[:, C['t']], z[:, C['s12']], color=col, lw=1, label=lab)
            ax.set_title(par, fontsize=9)
            ax.legend(fontsize=6)
            ax.tick_params(labelsize=7)
        for ax in list(axs.flat)[len(params):]:
            ax.axis('off')
        fig.suptitle(f'{m}, {s}, {pr}: sigma_12 (Pa) vs ' + ('t (s)' if kind == 'step' else 'gamma'))
        fig.tight_layout()
        fig.savefig(os.path.join(fdir, f'overlay_{m}_{s}_{pr}.pdf'))
        plt.close(fig)
    print(f'figures -> {fdir}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('action', choices=['run', 'analyse'])
    ap.add_argument('--models', nargs='+', default=list(MODELS), choices=list(MODELS))
    ap.add_argument('--systems', nargs='+', default=list(SYSTEMS), choices=list(SYSTEMS))
    ap.add_argument('--protocols', nargs='+', default=list(PROTOCOLS), choices=list(PROTOCOLS))
    ap.add_argument('--params', nargs='+', default=None, choices=list(LEVELS), help='default: all')
    ap.add_argument('--factor', type=int, default=6)
    ap.add_argument('--dt', type=float, default=None, help='override the maximum time step of the protocols')
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--out', default='sensitivity_runs', help='results folder (relative to python_scripts/)')
    ap.add_argument('--force', action='store_true', help='rerun even if the signature is unchanged')
    ap.add_argument('--no-plots', action='store_true')
    a = ap.parse_args()
    a.out = a.out if os.path.isabs(a.out) else os.path.join(HERE, a.out)
    cmd_run(a) if a.action == 'run' else cmd_analyse(a)


if __name__ == '__main__':
    main()
