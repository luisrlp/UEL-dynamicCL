# -*- coding: utf-8 -*-
#
# Outputs of the stage-1 study (run_stage1.py): single-element simple shear (ramp + hold),
# catch-slip vs slip kinetics. Results in SA_stage1/<sweep>/<case>/.
#
# Usage (from test_in_abaqus/):
#   abaqus python getoutput_stage1.py                          -> extract all runs in SA_stage1/
#   abaqus python getoutput_stage1.py SA_stage1/A_core/catch_gamma0.5   -> extract one run
#   abaqus python getoutput_stage1.py --plot                   -> extract all, then analyse/plot
#   python3 getoutput_stage1.py --plot-only                    -> only analyse/plot (no Abaqus), from
#                                                                 the saved stage1_output.pkl files
#   (abaqus cae noGUI=getoutput_stage1.py -- <args> also works)
#
# Extraction, per run -> <run folder>/stage1_output.pkl:
#   t (total time), step_names, step_ends, sig (n x 6: S11 S22 S33 S12 S13 S23, UVARM5-10),
#   thetaf (UVARM1), c (UVARM3), cb_tot (UVARM4), cb_dirs (n x NDIR: UVARM17-...), averaged over
#   the 8 integration points of the dummy element; rf1_top (sum of RF1 over the top nodes);
#   rfl (sum of RFL over the gel nodes: crosslinker exchange with the bath); props (properties.inp
#   of the run), label, warnings (counts in the Abaqus output file), completed (.sta).
#
# Analysis -> SA_stage1/stage1_summary.csv and SA_stage1/plots/<sweep>_*.pdf:
#   sigma_peak (S12 at the end of the ramp), R(10, 50, 200 s) = S12(hold time)/sigma_peak,
#   t_half (hold time to lose half of the stress that relaxes: S12 = S_end + (S_peak - S_end)/2),
#   settled (relaxation finished by the end of the hold), relative change of cb_tot, theta_f drift,
#   checks: sum(RF1 top)/area vs S12, integral of sum(RFL) vs change of total content c.
#   Plots per sweep: S12/sigma_peak vs hold time; t_half and R(200) vs GAMMA (or koff(0));
#   cb_tot/cb_tot(0) vs hold time; cb_i/cb0 at the end of the hold vs direction stretch lambda_i.
#   The direction vectors are generated in the same order as the Fortran loop
#   (affclnetfic_discrete: icosahedron faces 1-10, sub-triangles of order FACTOR).

import os
import re
import sys
import csv
import glob
import pickle
import numpy as np

try:
    from odbAccess import openOdb
    from abaqusConstants import INTEGRATION_POINT
    IN_ABAQUS = True
except ImportError:
    IN_ABAQUS = False

JOB_NAME = 'stage1_shear_hold'
sa_dir = os.path.join(os.getcwd(), 'SA_stage1')
file = f'{JOB_NAME}.odb'
OUT_PKL = 'stage1_output.pkl'
ELEM_OFFSET = 100000     # dummy (visualisation) elements: UEL element label + ElemOffset
AREA_TOP = 1.0           # top face area of the unit cube (um^2)
VOLUME = 1.0             # element volume (um^3); J = 1 in simple shear
HOLD_TIMES = (10.0, 50.0, 200.0)


def read_global_factor():
    """FACTOR (icosahedron subdivision) from ../global.f90, default 3."""
    path = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0] or '.')), '..', 'global.f90')
    if os.path.exists(path):
        m = re.search(r'PARAMETER\s*\(\s*FACTOR\s*=\s*(\d+)', open(path).read(), re.IGNORECASE)
        if m:
            return int(m.group(1))
    return 3


def fortran_directions(factor):
    """Unit direction vectors in the order of the Fortran loop (affclnetfic_discrete):
    icosahedron (icos_shape) faces 1..10, two sub-triangle orientations, sphere01_triangle_project."""
    phi = 0.5 * (np.sqrt(5.0) + 1.0)
    a, b, z = phi / np.sqrt(1.0 + phi * phi), 1.0 / np.sqrt(1.0 + phi * phi), 0.0
    pts = np.array([[a, b, z], [a, -b, z], [b, z, a], [b, z, -a], [z, a, b], [z, a, -b],
                    [z, -a, b], [z, -a, -b], [-b, z, a], [-b, z, -a], [-a, b, z], [-a, -b, z]])
    faces = [(1, 2, 4), (1, 3, 2), (1, 4, 6), (1, 5, 3), (1, 6, 5), (2, 3, 7), (2, 7, 8), (2, 8, 4),
             (3, 5, 9), (3, 9, 7), (4, 8, 10), (4, 10, 6), (5, 6, 11), (5, 11, 9), (6, 10, 11),
             (7, 9, 12), (7, 12, 8), (8, 12, 10), (9, 11, 12), (10, 12, 11)]
    f3_start, f3_end, f2_start = (1, 2), (3 * factor - 2, 3 * factor - 4), (1, 2)
    dirs = []
    for face in faces[:10]:
        pa, pb, pc = (pts[i - 1] for i in face)
        for k, ifacedir in enumerate((1, 2)):
            for f3 in range(f3_start[k], f3_end[k] + 1, 3):
                for f2 in range(f2_start[k], 3 * factor - f3 - ifacedir + 1, 3):
                    f1 = 3 * factor - f3 - f2
                    v = (f1 * pa + f2 * pb + f3 * pc) / (f1 + f2 + f3)
                    dirs.append(v / np.linalg.norm(v))
    return np.array(dirs)


def read_run_properties(run_dir):
    """Properties (NAME = value) and label (first '** ' line) of a run's properties.inp."""
    props, label = {}, ''
    path = os.path.join(run_dir, 'properties.inp')
    if os.path.exists(path):
        for line in open(path):
            line = line.strip()
            if line.startswith('** ') and not label:
                label = line[3:]
            elif '=' in line and not line.startswith('*'):
                key, value = (s.strip() for s in line.split('=', 1))
                try:
                    props[key] = float(eval(value, {'__builtins__': {}}))
                except Exception:
                    pass
    return props, label


def count_warnings(run_dir):
    path = os.path.join(run_dir, f'{JOB_NAME}_output.txt')
    text = open(path, errors='ignore').read() if os.path.exists(path) else ''
    return {'accuracy': text.count('accuracy check not met'),
            'kinetics_failed': text.count('kinetics failed'),
            'detF_negative': text.count('detF.lt.zero')}


def job_completed(run_dir):
    sta = os.path.join(run_dir, f'{JOB_NAME}.sta')
    return os.path.exists(sta) and 'COMPLETED SUCCESSFULLY' in open(sta, errors='ignore').read()


# ------------------------------------------------------------------------------------------
# Extraction (Abaqus Python)
# ------------------------------------------------------------------------------------------
def scalar_sum(values):
    return float(sum(np.sum(np.atleast_1d(v.data)) for v in values))


def extract(odb_files):
    print(f"ODB files to process: {len(odb_files)}")
    n_ok = 0
    for i_file, odb_file in enumerate(odb_files):
        run_dir = os.path.dirname(odb_file)
        folder = os.path.relpath(run_dir, sa_dir)
        print(f"Processing file {i_file + 1}/{len(odb_files)}: {folder}")
        if not os.path.exists(odb_file):
            print(f"  ODB not found, skipped: {odb_file}")
            continue
        if os.path.exists(odb_file.replace('.odb', '.lck')):
            print("  Job still running (.lck file present), skipped")
            continue
        try:
            odb = openOdb(path=odb_file, readOnly=True)
        except Exception as e:
            print(f"  Could not open ODB, skipped: {e}")
            continue

        # Flat input: one automatic instance holding the UEL nodes, dummy and extra elements
        inst = list(odb.rootAssembly.instances.values())[0]
        dummy_labels = sorted(e.label for e in inst.elements if e.label > ELEM_OFFSET)
        dummy = inst.getElementFromLabel(dummy_labels[0])
        top = inst.nodeSets['TOP_NODES']
        gel = inst.nodeSets['ALL_NODES']

        steps = list(odb.steps.values())
        names0 = steps[0].frames[-1].fieldOutputs.keys()
        n_uvarm = max(int(k[5:]) for k in names0 if k.startswith('UVARM') and k[5:].isdigit())
        t, uvarm, rf1_top, rfl = [], [], [], []
        for step in steps:
            for frame in step.frames:
                fo = frame.fieldOutputs
                if 'UVARM1' not in fo.keys():
                    continue
                row = []
                for k in range(1, n_uvarm + 1):
                    vals = fo[f'UVARM{k}'].getSubset(region=dummy, position=INTEGRATION_POINT).values
                    row.append(np.mean([v.data for v in vals]))
                uvarm.append(row)
                t.append(step.totalTime + frame.frameValue)
                rf1_top.append(sum(v.data[0] for v in fo['RF'].getSubset(region=top).values))
                rfl.append(scalar_sum(fo['RFL'].getSubset(region=gel).values) if 'RFL' in fo.keys() else np.nan)
        # Copy the step data before closing: ODB objects are invalid after odb.close()
        step_names = [s.name for s in steps]
        step_ends = [s.totalTime + s.timePeriod for s in steps]
        odb.close()

        uvarm = np.array(uvarm)
        props, label = read_run_properties(run_dir)
        output = {'t': np.array(t), 'step_names': step_names, 'step_ends': step_ends,
                  'thetaf': uvarm[:, 0], 'c': uvarm[:, 2], 'cb_tot': uvarm[:, 3],
                  'sig': uvarm[:, 4:10], 'cb_dirs': uvarm[:, 16:],
                  'rf1_top': np.array(rf1_top), 'rfl': np.array(rfl),
                  'props': props, 'label': label, 'folder': folder, 'odb': odb_file,
                  'warnings': count_warnings(run_dir), 'completed': job_completed(run_dir)}
        with open(os.path.join(run_dir, OUT_PKL), 'wb') as f:
            pickle.dump(output, f)
        print(f"  Saved {os.path.join(run_dir, OUT_PKL)} ({len(t)} frames)")
        n_ok += 1
    print(f"Done: {n_ok}/{len(odb_files)} ODB files processed.")


# ------------------------------------------------------------------------------------------
# Analysis and plots (plain Python)
# ------------------------------------------------------------------------------------------
def equilibrium_cb0(p, tol=1.e-14):
    """Equilibrium bound concentration per direction at t = 0 (as compute_initmu / pullchem)."""
    cr, fmax, cmax = p['CACTIN'] * p['R'], p['CACTIN'] * p['RFMAX'], p['CACTIN'] * p['RBMAX']
    def h(cb0):
        cf0 = cr - cb0
        th = cf0 / fmax
        return p['KEQ'] * cf0 * (1.0 - cb0 / cmax) - cb0 * (1.0 - th) * np.exp(-p['CHI'] * (1.0 - 2.0 * th))
    a, b = 0.0, min(cr, cmax)
    while (b - a) / 2.0 > tol:
        m = 0.5 * (a + b)
        if h(a) * h(m) <= 0.0:
            b = m
        else:
            a = m
    return 0.5 * (a + b)


def case_metrics(out, dirs):
    """Relaxation and consistency metrics of one run."""
    t, s12 = out['t'], out['sig'][:, 3]
    t_ramp = out['step_ends'][0]
    hold = t >= t_ramp - 1e-12
    th, sh = t[hold] - t_ramp, s12[hold]
    sig_peak, sig_end = sh[0], sh[-1]
    m = {'sigma_peak': sig_peak, 'sigma_end': sig_end}
    for tt in HOLD_TIMES:
        m[f'R{tt:g}'] = np.interp(tt, th, sh) / sig_peak if th[-1] >= tt and sig_peak != 0 else np.nan
    # t_half: first hold time at which (S - S_end)/(S_peak - S_end) <= 0.5 (relaxation or rise)
    dsig = sig_peak - sig_end
    m['t_half'] = np.nan
    if abs(dsig) > 1e-3 * abs(sig_peak):
        s_norm = (sh - sig_end) / dsig
        idx = np.nonzero(s_norm <= 0.5)[0]
        if len(idx) and idx[0] > 0:
            i = idx[0]
            m['t_half'] = np.interp(0.5, [s_norm[i], s_norm[i - 1]], [th[i], th[i - 1]])
    # settled: change over the last 10% of the hold small compared with the total relaxation
    t10 = th[-1] - 0.1 * th[-1]
    m['settled'] = bool(abs(sig_end - np.interp(t10, th, sh)) <= 0.05 * abs(dsig)) if abs(dsig) > 0 else True
    cb = out['cb_tot']
    m['dcb_tot_rel'] = cb[-1] / cb[0] - 1.0
    m['thetaf_drift'] = np.max(np.abs(out['thetaf'] / out['thetaf'][0] - 1.0))
    # Checks: reaction force vs stress (end of ramp); bath exchange vs change of total content
    i_peak = np.nonzero(hold)[0][0]
    m['rf_check'] = (out['rf1_top'][i_peak] / AREA_TOP) / sig_peak - 1.0 if sig_peak != 0 else np.nan
    trapezoid = getattr(np, 'trapezoid', None) or np.trapz    # np.trapz removed in recent NumPy
    rfl_int = trapezoid(np.nan_to_num(out['rfl']), t)
    dc = (out['c'][-1] - out['c'][0]) * VOLUME
    m['rfl_integral'], m['dc_total'] = rfl_int, dc
    # Direction stretches (simple shear F = I + GAMMA e_x (x) e_y, J = 1) and final cb_i/cb0
    gamma = out['props'].get('GAMMA', np.nan)
    F = np.eye(3)
    F[0, 1] = gamma
    m['lambda_i'] = np.linalg.norm(dirs @ F.T, axis=1)
    m['cb_rel_end'] = out['cb_dirs'][-1, :len(dirs)] / equilibrium_cb0(out['props'])
    return m


def group_and_x(sweep, case, props):
    """Line group and x value of a case: GAMMA sweeps grouped by the folder name without
    '_gamma<value>'; E_rate grouped by law with x = koff(0)."""
    if '_gamma' in case or case.startswith('gamma'):
        return re.sub(r'_?gamma[^_]*$', '', case) or 'all', props.get('GAMMA', np.nan), 'GAMMA'
    if sweep == 'E_rate':
        return re.sub(r'_koff[^_]*$', '', case), props['KOFF0'] + props['KCATCH0'], 'koff(0) (1/s)'
    return case, np.nan, ''


def load_outputs():
    outs = {}
    for path in sorted(glob.glob(os.path.join(sa_dir, '*', '*', OUT_PKL))):
        with open(path, 'rb') as f:
            out = pickle.load(f)
        out['folder'] = os.path.relpath(os.path.dirname(path), sa_dir)
        outs[out['folder']] = out
    return outs


def analyse(outs):
    dirs = fortran_directions(read_global_factor())
    rows, by_sweep = [], {}
    for folder, out in outs.items():
        sweep, case = folder.replace('\\', '/').split('/')[:2]
        m = case_metrics(out, dirs)
        group, x, xlabel = group_and_x(sweep, case, out['props'])
        p = out['props']
        rows.append({'sweep': sweep, 'case': case, 'label': out['label'], 'group': group, 'x': x,
                     **{k: p.get(k, np.nan) for k in ('GAMMA', 'KOFF0', 'KCATCH0', 'DX', 'DXC', 'ETAC')},
                     **{k: v for k, v in m.items() if k not in ('lambda_i', 'cb_rel_end')},
                     **{f'warn_{k}': v for k, v in out['warnings'].items()}, 'completed': out['completed']})
        by_sweep.setdefault(sweep, []).append((group, x, xlabel, out, m))
    csv_path = os.path.join(sa_dir, 'stage1_summary.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in sorted(rows, key=lambda r: (r['sweep'], r['group'], r['x'])):
            writer.writerow({k: (f'{v:.6g}' if isinstance(v, float) else v) for k, v in r.items()})
    print(f"Summary written to {csv_path} ({len(rows)} runs)")
    unsettled = [f"{r['sweep']}/{r['case']}" for r in rows if not r['settled']]
    if unsettled:
        print(f"WARNING: stress still changing at the end of the hold (t_half underestimated): {unsettled}")
    return by_sweep


def plot_all(by_sweep):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not available in this Python. Plot the extracted results with a "
              "Python that has it, e.g. '<venv>/bin/python getoutput_stage1.py --plot-only'")
        return
    plot_dir = os.path.join(sa_dir, 'plots')
    os.makedirs(plot_dir, exist_ok=True)

    for sweep, items in sorted(by_sweep.items()):
        groups = sorted(set(it[0] for it in items))
        xlabel = items[0][2] or 'case'
        xs = sorted(set(it[1] for it in items if np.isfinite(it[1])))
        cmap = plt.cm.viridis
        color = {x: cmap(v) for x, v in zip(xs, np.linspace(0.0, 0.9, max(len(xs), 1)))}
        ncol = min(len(groups), 3)
        nrow = int(np.ceil(len(groups) / ncol))

        def panels(title, ylabel, plot_fn, name, xlab='Hold time (s)'):
            fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.3 * nrow), squeeze=False, sharey=True)
            for ax, grp in zip(axs.flat, groups):
                for group, x, _, out, m in sorted((it for it in items if it[0] == grp), key=lambda it: it[1]):
                    plot_fn(ax, out, m, color.get(x, 'k'), f'{xlabel} = {x:g}' if np.isfinite(x) else group)
                ax.set_title(grp, fontsize=9)
                ax.set_xlabel(xlab)
                ax.grid(True, alpha=0.3)
                ax.legend(fontsize=7, frameon=False)
            for ax in axs[:, 0]:
                ax.set_ylabel(ylabel)
            for ax in list(axs.flat)[len(groups):]:
                ax.set_visible(False)
            fig.suptitle(f'{sweep}: {title}')
            fig.tight_layout()
            path = os.path.join(plot_dir, f'{sweep}_{name}.pdf')
            fig.savefig(path)
            plt.close(fig)
            print(f"Plot saved: {path}")

        def hold_series(out, y):
            t_ramp = out['step_ends'][0]
            hold = out['t'] >= t_ramp - 1e-12
            return out['t'][hold] - t_ramp, y[hold]

        panels('shear stress during the hold', r'$\sigma_{12}/\sigma_{12,peak}$',
               lambda ax, out, m, c, lab: ax.plot(*hold_series(out, out['sig'][:, 3] / m['sigma_peak']), color=c, label=lab),
               'stress_time')
        panels('total bound crosslinkers during the hold', r'$c_{b,tot}/c_{b,tot}(0)$',
               lambda ax, out, m, c, lab: ax.plot(*hold_series(out, out['cb_tot'] / out['cb_tot'][0]), color=c, label=lab),
               'cbtot_time')
        panels('bound crosslinkers per direction at the end of the hold', r'$c_{b,i}/c_{b0}$',
               lambda ax, out, m, c, lab: ax.scatter(m['lambda_i'], m['cb_rel_end'], s=10, color=c, label=lab),
               'cb_dirs', xlab=r'Direction stretch $\lambda_i$')

        # Relaxation metrics vs x (the catch-slip signature: t_half vs GAMMA)
        if xs:
            fig, axs = plt.subplots(1, 2, figsize=(9.0, 3.6))
            for grp, marker in zip(groups, 'osD^v<>ph*'):
                pts = sorted((it[1], it[4]) for it in items if it[0] == grp and np.isfinite(it[1]))
                x = [p[0] for p in pts]
                axs[0].plot(x, [p[1]['t_half'] for p in pts], marker=marker, label=grp)
                axs[1].plot(x, [p[1]['R200'] for p in pts], marker=marker, label=grp)
            axs[0].set_ylabel(r'$t_{1/2}$ (s)')
            axs[1].set_ylabel(r'$R(200\,s) = \sigma_{12}/\sigma_{12,peak}$')
            for ax in axs:
                ax.set_xlabel(xlabel)
                ax.grid(True, alpha=0.3)
                ax.legend(fontsize=7, frameon=False)
            fig.suptitle(f'{sweep}: relaxation vs {xlabel}')
            fig.tight_layout()
            path = os.path.join(plot_dir, f'{sweep}_relaxation_vs_x.pdf')
            fig.savefig(path)
            plt.close(fig)
            print(f"Plot saved: {path}")


# ------------------------------------------------------------------------------------------
if __name__ == '__main__':
    user_args = sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else sys.argv[1:]
    plot_only = '--plot-only' in user_args
    do_plot = plot_only or '--plot' in user_args
    paths = [a for a in user_args if not a.startswith('--')]

    if not plot_only:
        if not IN_ABAQUS:
            sys.exit("Extraction needs Abaqus: 'abaqus python getoutput_stage1.py' "
                     "(or use --plot-only to analyse existing results)")
        if paths:
            target = os.path.abspath(paths[0])
            odb_files = [target] if target.endswith('.odb') else [os.path.join(target, file)]
        else:
            odb_files = sorted(glob.glob(os.path.join(sa_dir, '*', '*', file)))
        extract(odb_files)

    if do_plot:
        outs = load_outputs()
        if not outs:
            sys.exit(f"No {OUT_PKL} files found in {sa_dir}")
        plot_all(analyse(outs))
