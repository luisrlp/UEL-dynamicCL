"""Step 2 of the chain-relaxation checks (chain_relaxation_rules.md): does the strain stiffening of Fig. 4 survive
when the reference stretch relaxes through the chain mechanism at the rate needed by Fig. 5?

Model : material_chain.py (DynamicCrosslinking), simple shear at one material point, closed system, LAOS from rest at
        f = 0.05 Hz, the protocol of optimise_laos_v2 (same FIXED_DIFF/FIXED_MECH, parameter mapping, data).
Base  : the optimum of optimisation_runs/c5_w15_plate_guo_dxs (not refitted; copied into BASE_PARAMS, since
        optimisation_runs/ is not tracked by git).
Cases : (relaxation, <L_fil> [um], target rest rate k = B/(2 c_b0) [1/s]; None = the base koff0)
          ref       current model (SFR), base koff0
          sfr_k035  control: current model with koff0 scaled to k = 0.35 (single relaxation time at the option-1 speed)
          pops7     6 populations, <L> 7 um,  k 0.35 (option 1: spontaneous polymerisation, k ~ literature)
          pops15    6 populations, <L> 15 um, k 0.87 (option 2: best Fig. 5 fit)
          chain7    exact half-chains, <L> 7 um, k 0.35 (checks the 6-population reduction; expensive)
        The rate is raised by scaling k_off_0, k_catch_0 and k_on_0 together (keq fixed), so the rest state
        (c_b0, r_0f, l_c, G_0) is unchanged and only the rates change. k targets: linear fit of the exact mixture to
        Fig. 5 (f <= 0.8 Hz), see chain_relaxation_rules.md.
Runs  : one LAOS run from rest per case and amplitude (--gammas; 0.125875 is the Fig. 4a loop).
Output (--out): runs/<case>_g<gamma>.npz (t, gamma, sigma of all cycles), summary.csv, setup_<cases>.json.
        Finished runs are skipped on restart if their settings match (delete runs/ to recompute).
        Per run: G', G'' of the last cycle (main.moduli, both halves), G_M (slope at gamma = 0) and G_L (secant at
        gamma_A) of the last cycle, cycle-to-cycle RMS change, peaks, range of c_b/c_b0 over directions and time
        (checks the frozen populations); at gamma_A of Fig. 4a also the loop RMS vs data.

  python -B laos_chain_check.py --workers 8                          # all cases, default amplitudes
  python -B laos_chain_check.py --cases ref pops7 --workers 4
  python -B laos_chain_check.py --summary                            # only print the table of finished runs
Keep --workers <= the number of physical cores; pin BLAS threads (OMP/OPENBLAS/MKL_NUM_THREADS=1).
"""
import argparse, contextlib, io, json, os, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import optimise_laos_v2 as O

BASE_RUN = 'c5_w15_plate_guo_dxs'
BASE_PARAMS = {'lambda0': 1.0270240289881567, 's_taut': 1.0745682572324684, 'koff0': 0.06321273017602523,
               'phi': 0.907, 'dxc_nm': 1.5, 'dxs_nm': 1.132504971376299, 'keq': 14.0, 'beta': 0.438,
               'mu0stretch': 38600.0}
CASES = {'ref':      ('sfr',   None, None),
         'sfr_k035': ('sfr',   None, 0.35),
         'pops7':    ('pops',  7.0,  0.35),
         'pops15':   ('pops',  15.0, 0.87),
         'chain7':   ('chain', 7.0,  0.35)}
GAMMAS = (0.0126, 0.05, 0.1, 0.125875, 0.158)
N_POPS = 6


def base_params():
    return dict(BASE_PARAMS)


def rest_rate(p, factor):
    """Rest rate per crosslink k = B/(2 c_b0) of the current model for the parameter set p."""
    from material import DynamicCrosslinking
    diff, mech = O.material_props(p)
    sim = {'deformation_type': 'shear', 'factor': factor, 'max_iter': 50, 'tol': 1.0e-12}
    with contextlib.redirect_stdout(io.StringIO()):
        m = DynamicCrosslinking(diff, mech, sim)
    B = m._binding_rate(m.cf0 / m.CF_MAX)[1]
    return B / (2.0 * m.cb0)


def case_params(case, cfg):
    """Physical parameters, extra mech properties and koff scale of a case."""
    relax, L, k_target = CASES[case]
    p = base_params()
    scale = 1.0 if k_target is None else k_target / rest_rate(p, cfg['factor'])
    p = dict(p, koff0=p['koff0'] * scale)
    extra = {'relaxation': relax}
    if relax != 'sfr':
        extra['L_mean'] = L
        extra['n_pops'] = N_POPS
    return p, extra, scale


def simulate(p, extra, cfg, gA):
    """LAOS from rest at amplitude gA (optionally plate-averaged); time, gamma and sigma_xy of all cycles, and the
    range of c_b/c_b0 over the directions per step (min, max; over all plate nodes)."""
    if cfg['plate_nodes'] > 0:
        x, w = O.plate_nodes(cfg['plate_nodes'])
        runs = [simulate(p, extra, dict(cfg, plate_nodes=0), xi * gA) for xi in x]
        t = runs[0][0]
        cbr = np.column_stack([np.min([r[3][:, 0] for r in runs], 0), np.max([r[3][:, 1] for r in runs], 0)])
        return t, gA * np.sin(2 * np.pi * O.FREQ * t), sum(wi * r[2] for wi, r in zip(w, runs)), cbr
    from material_chain import DynamicCrosslinking
    diff, mech = O.material_props(p)
    mech.update(extra)
    sim = {'deformation_type': 'shear', 'factor': cfg['factor'], 'max_iter': 50, 'tol': 1.0e-12}
    t = np.arange(0.0, cfg['cycles'] / O.FREQ + 0.5 * cfg['dt'], cfg['dt'])
    g = gA * np.sin(2 * np.pi * O.FREQ * t)
    s = np.empty_like(t)
    cbr = np.empty((len(t), 2))
    with contextlib.redirect_stdout(io.StringIO()), np.errstate(all='ignore'):
        m = DynamicCrosslinking(diff, mech, sim)
        for i in range(len(t)):
            F = np.array([[1.0, g[i], 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
            s[i] = m.evaluate_state(F, None, 0.0 if i == 0 else t[i] - t[i - 1])[0, 1]
            cbr[i] = m.cb_history.min() / m.cb0, m.cb_history.max() / m.cb0
    return t, g, s, cbr


def run_file(out, case, gA):
    return os.path.join(out, 'runs', f'{case}_g{gA:.6g}.npz')


def run_job(args):
    """Worker: one case at one amplitude; saves the npz and returns its wall time (or the error)."""
    case, gA, cfg, out = args
    t0 = time.time()
    try:
        p, extra, scale = case_params(case, cfg)
        t, g, s, cbr = simulate(p, extra, cfg, gA)
        np.savez(run_file(out, case, gA), t=t, g=g, s=s, cb_range=cbr, cfg=json.dumps(cfg), wall=time.time() - t0,
                 scale=scale)
        return case, gA, time.time() - t0, None
    except Exception as e:
        return case, gA, time.time() - t0, repr(e)


def branch_slope(g, s, gA):
    """Slope of one loop branch at gamma = 0: central difference over gamma = +-0.05 gamma_A (interpolated)."""
    o = np.argsort(g, kind='stable')
    d = 0.05 * gA
    return float((np.interp(d, g[o], s[o]) - np.interp(-d, g[o], s[o])) / (2 * d))


def slope_at_zero(gm, sm, gA):
    """G_M: mean slope of both branches of the cycle at gamma = 0."""
    dg = np.gradient(gm)
    return float(np.mean([branch_slope(gm[msk], sm[msk], gA) for msk in (dg > 0, dg < 0)]))


def analyse(out, case, gA, data):
    from main import moduli
    f = run_file(out, case, gA)
    if not os.path.exists(f):
        return None
    z = np.load(f)
    t, g, s = z['t'], z['g'], z['s']
    cyc = O.split_cycles(t, g, s, json.loads(str(z['cfg']))['cycles'])
    gm, sm = cyc[-1]
    last = t >= t[-1] - 1.0 / O.FREQ - 1e-9
    G1, G2, _ = moduli(t[last], g[last], s[last], O.FREQ, both_halves=True)
    i = np.argmax(gm)
    row = {'case': case, 'gamma_A': gA, 'koff_scale': float(z['scale']), 'wall_h': float(z['wall']) / 3600,
           "G'": G1, "G''": G2, 'G_M': slope_at_zero(gm, sm, gA), 'G_L': float(sm[i] / gm[i]),
           'smax': float(sm.max()), 'smin': float(sm.min()), 'smax_c1': float(cyc[0][1].max()),
           'cb_min': float(z['cb_range'][:, 0].min()) if 'cb_range' in z else float('nan'),
           'cb_max': float(z['cb_range'][:, 1].max()) if 'cb_range' in z else float('nan'),
           'dprev_last': float('nan'), 'loop_rms': float('nan')}
    if len(cyc) > 1:
        a, b = cyc[-2], cyc[-1]
        n = min(len(a[1]), len(b[1]))
        row['dprev_last'] = float(np.sqrt(np.mean((a[1][-n:] - b[1][-n:])**2)))
    if abs(gA - data['gamma_A']) < O.REUSE_TOL:
        unweighted = {k: (x, y, np.ones_like(x)) for k, (x, y, w) in data['loop'].items()}
        row['loop_rms'] = float(np.sqrt(np.mean(O.residuals_from_loop(gm, sm, unweighted)**2)))
    return row


def summary(out, cases, gammas):
    loop, gA_data = O.load_data()
    data = {'loop': loop, 'gamma_A': gA_data}
    lo, hi = O.sweep_range()
    rows = [r for c in cases for gA in gammas if (r := analyse(out, c, gA, data)) is not None]
    if not rows:
        print('no finished runs')
        return
    keys = list(rows[0].keys())
    with open(os.path.join(out, 'summary.csv'), 'w') as fh:
        fh.write(','.join(keys) + '\n')
        for r in rows:
            fh.write(','.join(str(r[k]) for k in keys) + '\n')
    # data: G', G'' of Fig. 4b, and G_M, G_L of the Fig. 4a loop
    Gd = {gA: O.load_sweep_data(np.array([gA]))[0] for gA in gammas if lo <= gA <= hi}
    dg, ds = np.concatenate([loop[1][0], loop[-1][0]]), np.concatenate([loop[1][1], loop[-1][1]])
    print(f"\nData Fig. 4a loop (gamma_A {gA_data:.4f}): G_M {slope_at_zero_data(loop, gA_data):.2f} Pa, "
          f"G_L {ds[np.argmax(dg)] / dg.max():.2f} Pa, smax {ds.max():.2f} Pa")
    print(f"{'case':>9} {'gamma_A':>8} {'Gp':>6} {'Gp/data':>7} {'Gpp':>6} {'Gpp/data':>8} {'Gp/Gplin':>8} "
          f"{'data':>5} {'G_M':>6} {'G_L':>6} {'smax':>6} {'smax_c1':>7} {'dprev':>6} {'loopRMS':>7} {'cb/cb0':>11} {'h':>5}")
    for c in cases:
        rc = [r for r in rows if r['case'] == c]
        lin = next((r for r in rc if r['gamma_A'] == min(gammas)), None)
        for r in rc:
            gA = r['gamma_A']
            d = Gd.get(gA)
            q1 = r["G'"] / d[0] if d is not None else float('nan')
            q2 = r["G''"] / d[1] if d is not None else float('nan')
            st = r["G'"] / lin["G'"] if lin else float('nan')
            dst = d[0] / Gd[min(gammas)][0] if d is not None and min(gammas) in Gd else float('nan')
            G1, G2 = r["G'"], r["G''"]
            print(f"{c:>9} {gA:8.4f} {G1:6.2f} {q1:7.2f} {G2:6.2f} {q2:8.2f} "
                  f"{st:8.2f} {dst:5.2f} {r['G_M']:6.2f} {r['G_L']:6.2f} {r['smax']:6.2f} {r['smax_c1']:7.2f} "
                  f"{r['dprev_last']:6.3f} {r['loop_rms']:7.3f} {r['cb_min']:5.3f}-{r['cb_max']:5.3f} {r['wall_h']:5.2f}")
    print(f"\nkoff scale per case: " + ', '.join(f"{c} {next(r['koff_scale'] for r in rows if r['case'] == c):.2f}"
                                               for c in cases if any(r['case'] == c for r in rows)))
    print(f"summary.csv written to {out}")


def slope_at_zero_data(loop, gA):
    return float(np.mean([branch_slope(x, y, gA) for x, y, w in loop.values()]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cases', nargs='+', default=list(CASES), choices=list(CASES))
    ap.add_argument('--gammas', nargs='+', type=float, default=list(GAMMAS))
    ap.add_argument('--dt', type=float, default=0.1)
    ap.add_argument('--factor', type=int, default=3)
    ap.add_argument('--cycles', type=int, default=5)
    ap.add_argument('--chain-cycles', type=int, default=None, help='cycles of chain7 (default: --cycles)')
    ap.add_argument('--plate-nodes', type=int, default=0, help='parallel-plate averaging nodes (0: point stress)')
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--out', default=os.path.join(HERE, 'optimisation_runs', 'laos_chain_check'))
    ap.add_argument('--summary', action='store_true', help='only analyse finished runs')
    a = ap.parse_args()

    os.makedirs(os.path.join(a.out, 'runs'), exist_ok=True)
    cfg = {'dt': a.dt, 'factor': a.factor, 'cycles': a.cycles, 'plate_nodes': a.plate_nodes}
    cfg_case = {c: dict(cfg, cycles=a.chain_cycles) if c == 'chain7' and a.chain_cycles else cfg for c in a.cases}
    if not a.summary:
        setup = {'base_run': BASE_RUN, 'base_params': base_params(), 'cases': {c: CASES[c] for c in a.cases},
                 'n_pops': N_POPS, 'gammas': a.gammas, 'cfg': cfg_case, 'FIXED_DIFF': O.FIXED_DIFF,
                 'FIXED_MECH': {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in O.FIXED_MECH.items()}}
        with open(os.path.join(a.out, f"setup_{'-'.join(a.cases)}.json"), 'w') as fh:
            json.dump(setup, fh, indent=1)
        jobs = []
        for c in a.cases:
            for gA in a.gammas:
                f = run_file(a.out, c, gA)
                if os.path.exists(f) and json.loads(str(np.load(f)['cfg'])) == cfg_case[c]:
                    print(f'skip {c} gamma_A {gA} (done)')
                    continue
                jobs.append((c, gA, cfg_case[c], a.out))
        # expensive jobs first (chain, then pops, larger amplitudes)
        order = {'chain': 0, 'pops': 1, 'sfr': 2}
        jobs.sort(key=lambda j: (order[CASES[j[0]][0]], -j[1]))
        print(f'{len(jobs)} runs, {a.workers} workers, out {a.out}', flush=True)
        t0 = time.time()
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            for fut in as_completed([ex.submit(run_job, j) for j in jobs]):
                case, gA, wall, err = fut.result()
                msg = f'FAILED: {err}' if err else 'done'
                print(f'[{(time.time() - t0) / 3600:6.2f} h] {case:>9} gamma_A {gA:<8g} {wall / 60:7.1f} min  {msg}',
                      flush=True)
    summary(a.out, a.cases, a.gammas)


if __name__ == '__main__':
    main()
