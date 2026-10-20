# -*- coding: utf-8 -*-
#
# Outputs of the stage-2 study (run_stage2.py): column along x with a bath on the end face
# x = L, homogeneous simple shear (ramp + hold). Results in SA_stage2/<sweep>/<case>/.
#
# Usage (from test_in_abaqus/):
#   abaqus python getoutput_stage2.py                          -> extract all runs in SA_stage2/
#   abaqus python getoutput_stage2.py SA_stage2/A_core/catch_gamma0.5   -> extract one run
#   abaqus python getoutput_stage2.py --plot                   -> extract all, then analyse/plot
#   python3 getoutput_stage2.py --plot-only                    -> only analyse/plot (no Abaqus), from
#                                                                 the saved stage2_output.pkl files
#   (abaqus cae noGUI=getoutput_stage2.py -- <args> also works)
#
# Extraction, per run -> <run folder>/stage2_output.pkl:
#   t (total time), step_names, step_ends; per element (UEL element = dummy element - offset),
#   sorted by position along the column: x (centroid), vol (volume); per frame and element
#   (average over the 8 integration points): thetaf (UVARM1), c (UVARM3, total content),
#   cb_tot (UVARM4), s12 (UVARM8); per frame: rfl_bath (sum of RFL over the bath nodes:
#   crosslinker exchange with the bath), rf1_top (sum of RF1 over the top nodes); x_bath, area_top,
#   props, label, warnings, completed.
#
# Analysis -> SA_stage2/stage2_summary.csv and SA_stage2/plots/<sweep>_*.pdf:
#   uptake(t) = sum_e (c_e(t) - c_e(0)) V_e (net crosslinker taken up from (> 0) or released to
#   (< 0) the bath), relative to the initial total content; t_half_uptake; penetration depth
#   (distance from the bath face where |dc| falls to 10% of its maximum, at 100 s, 300 s and the
#   end of the hold; the end value saturates once the front reaches the closed end); conservation check
#   (time integral of the bath flux vs uptake); average shear stress sum(RF1)/area: peak, end,
#   R_end = end/peak; theta_f at the bath face and at the closed end.
#   Plots per sweep: profiles along x at several hold times (relative change of c and cb_tot,
#   theta_f, S12); uptake vs hold time; average S12 vs hold time; summary vs the swept parameter.

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

JOB_NAME = 'stage2_column_bath'
sa_dir = os.path.join(os.getcwd(), 'SA_stage2')
file = f'{JOB_NAME}.odb'
OUT_PKL = 'stage2_output.pkl'
ELEM_OFFSET = 100000                                # dummy elements: UEL element label + ElemOffset
PROFILE_TIMES = (0.0, 10.0, 50.0, 100.0, 300.0, 1000.0)   # hold times of the profiles (s)
PENETRATION_FRACTION = 0.1
PENETRATION_TIMES = (100.0, 300.0)                  # hold times of the penetration depth (s), plus the end


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


def element_field(fo, name, region, index):
    """Per-element average over the integration points of a scalar field (UVARM)."""
    sums, counts = np.zeros(len(index)), np.zeros(len(index))
    for v in fo[name].getSubset(region=region, position=INTEGRATION_POINT).values:
        k = index.get(v.elementLabel)
        if k is not None:
            sums[k] += v.data
            counts[k] += 1
    return sums / np.maximum(counts, 1)


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

        # Flat input: one automatic instance with the UEL nodes, dummy and extra elements
        inst = list(odb.rootAssembly.instances.values())[0]
        coords = {n.label: np.array(n.coordinates) for n in inst.nodes}
        dummy = inst.elementSets['DUMMY_MESH']
        elems = sorted(dummy.elements, key=lambda e: np.mean([coords[n][0] for n in e.connectivity]))
        index = {e.label: k for k, e in enumerate(elems)}
        x_el, vol = [], []
        for e in elems:
            p = np.array([coords[n] for n in e.connectivity])
            x_el.append(p[:, 0].mean())
            vol.append(np.prod(p.max(axis=0) - p.min(axis=0)))   # box elements
        bath = inst.nodeSets['OPEN_CHEM_NODES']
        top = inst.nodeSets['TOP_NODES']
        x_bath = float(np.mean([n.coordinates[0] for n in bath.nodes]))
        top_xz = np.array([[n.coordinates[0], n.coordinates[2]] for n in top.nodes])
        area_top = float(np.ptp(top_xz[:, 0]) * np.ptp(top_xz[:, 1]))

        steps = list(odb.steps.values())
        t, thetaf, c, cb_tot, s12, rfl_bath, rf1_top = [], [], [], [], [], [], []
        for step in steps:
            for frame in step.frames:
                fo = frame.fieldOutputs
                if 'UVARM1' not in fo.keys():
                    continue
                t.append(step.totalTime + frame.frameValue)
                thetaf.append(element_field(fo, 'UVARM1', dummy, index))
                c.append(element_field(fo, 'UVARM3', dummy, index))
                cb_tot.append(element_field(fo, 'UVARM4', dummy, index))
                s12.append(element_field(fo, 'UVARM8', dummy, index))
                rfl_bath.append(scalar_sum(fo['RFL'].getSubset(region=bath).values) if 'RFL' in fo.keys() else np.nan)
                rf1_top.append(sum(v.data[0] for v in fo['RF'].getSubset(region=top).values))
        # Copy the step data before closing: ODB objects are invalid after odb.close()
        step_names = [s.name for s in steps]
        step_ends = [s.totalTime + s.timePeriod for s in steps]
        odb.close()

        props, label = read_run_properties(run_dir)
        output = {'t': np.array(t), 'step_names': step_names, 'step_ends': step_ends,
                  'x': np.array(x_el), 'vol': np.array(vol), 'x_bath': x_bath, 'area_top': area_top,
                  'thetaf': np.array(thetaf), 'c': np.array(c), 'cb_tot': np.array(cb_tot),
                  's12': np.array(s12), 'rfl_bath': np.array(rfl_bath), 'rf1_top': np.array(rf1_top),
                  'props': props, 'label': label, 'folder': folder, 'odb': odb_file,
                  'warnings': count_warnings(run_dir), 'completed': job_completed(run_dir)}
        with open(os.path.join(run_dir, OUT_PKL), 'wb') as f:
            pickle.dump(output, f)
        print(f"  Saved {os.path.join(run_dir, OUT_PKL)} ({len(t)} frames, {len(elems)} elements)")
        n_ok += 1
    print(f"Done: {n_ok}/{len(odb_files)} ODB files processed.")


# ------------------------------------------------------------------------------------------
# Analysis and plots (plain Python)
# ------------------------------------------------------------------------------------------
def hold_index(out):
    """Frames of the hold step and the hold time (t - end of the ramp)."""
    t_ramp = out['step_ends'][0]
    hold = out['t'] >= t_ramp - 1e-12
    return hold, out['t'][hold] - t_ramp


def frame_at(out, t_hold):
    """Index of the frame closest to a given hold time."""
    hold, th = hold_index(out)
    return np.nonzero(hold)[0][np.argmin(np.abs(th - t_hold))]


def penetration_depth(out, i):
    """Distance from the bath face where |c - c(0)| falls to PENETRATION_FRACTION of its maximum."""
    dc = np.abs(out['c'][i] - out['c'][0])
    if dc.max() <= 0.0:
        return 0.0
    d = out['x_bath'] - out['x']                        # distance from the bath face
    order = np.argsort(d)
    d, dc = d[order], dc[order] / dc.max()
    below = np.nonzero(dc <= PENETRATION_FRACTION)[0]
    if len(below) == 0:
        return float(d[-1])                             # whole column affected
    k = below[0]
    if k == 0:
        return float(d[0])
    return float(np.interp(PENETRATION_FRACTION, [dc[k], dc[k - 1]], [d[k], d[k - 1]]))


def case_metrics(out):
    trapezoid = getattr(np, 'trapezoid', None) or np.trapz    # np.trapz removed in recent NumPy
    t, vol = out['t'], out['vol']
    hold, th = hold_index(out)
    total0 = float(np.sum(out['c'][0] * vol))
    uptake = (out['c'] - out['c'][0]) @ vol             # net amount taken up from the bath
    m = {'total0': total0, 'uptake_t': uptake, 'uptake_end': float(uptake[-1]),
         'uptake_rel_end': float(uptake[-1] / total0)}
    # t_half of the exchange: hold time to reach half of the final uptake
    uh = uptake[hold]
    m['t_half_uptake'] = np.nan
    if abs(uptake[-1]) > 1e-9 * total0:
        s = uh / uptake[-1]
        idx = np.nonzero(s >= 0.5)[0]
        if len(idx):
            i = idx[0]
            m['t_half_uptake'] = th[i] if i == 0 else float(np.interp(0.5, [s[i - 1], s[i]], [th[i - 1], th[i]]))
    for tp in PENETRATION_TIMES:          # while the front is still inside the column
        m[f'penetration_{tp:g}s'] = penetration_depth(out, frame_at(out, tp)) if th[-1] >= tp else np.nan
    m['penetration_end'] = penetration_depth(out, len(t) - 1)
    # Conservation: time integral of the bath flux vs the change of total content
    m['rfl_integral'] = float(trapezoid(np.nan_to_num(out['rfl_bath']), t))
    m['conservation_ratio'] = m['rfl_integral'] / uptake[-1] if uptake[-1] != 0 else np.nan
    # Average shear stress from the reaction force (homogeneous deformation)
    s_avg = out['rf1_top'] / out['area_top']
    m['sigma_peak'], m['sigma_end'] = float(s_avg[hold][0]), float(s_avg[-1])
    m['R_end'] = m['sigma_end'] / m['sigma_peak'] if m['sigma_peak'] != 0 else np.nan
    d = out['x_bath'] - out['x']
    i_bath, i_far = np.argmin(d), np.argmax(d)
    m['thetaf_bath_end'] = float(out['thetaf'][-1, i_bath] / out['thetaf'][0, i_bath])
    m['thetaf_far_end'] = float(out['thetaf'][-1, i_far] / out['thetaf'][0, i_far])
    return m


SWEEP_VAR = {'gamma': ('GAMMA', 'GAMMA'), 'D': ('D', r'D ($\mu$m$^2$/s)'),
             'keq': ('KEQ', 'KEQ'), 'rfmax': ('RFMAX', 'RFMAX')}


def group_and_x(case, props):
    """Law (group) and swept value of a case folder '<law>_<var><value>'."""
    m = re.match(r'^([A-Za-z]+)_([A-Za-z]+)([-+.\deE]+)$', case)
    if m and m.group(2) in SWEEP_VAR:
        key, label = SWEEP_VAR[m.group(2)]
        return m.group(1), props.get(key, float(m.group(3))), label
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
    rows, by_sweep = [], {}
    for folder, out in outs.items():
        sweep, case = folder.replace('\\', '/').split('/')[:2]
        m = case_metrics(out)
        group, x, xlabel = group_and_x(case, out['props'])
        p = out['props']
        rows.append({'sweep': sweep, 'case': case, 'label': out['label'], 'group': group, 'x': x,
                     **{k: p.get(k, np.nan) for k in ('GAMMA', 'KOFF0', 'KCATCH0', 'D', 'KEQ', 'RFMAX')},
                     **{k: v for k, v in m.items() if k != 'uptake_t'},
                     **{f'warn_{k}': v for k, v in out['warnings'].items()}, 'completed': out['completed']})
        by_sweep.setdefault(sweep, []).append((group, x, xlabel, out, m))
    csv_path = os.path.join(sa_dir, 'stage2_summary.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in sorted(rows, key=lambda r: (r['sweep'], r['group'], r['x'])):
            writer.writerow({k: (f'{v:.6g}' if isinstance(v, float) else v) for k, v in r.items()})
    print(f"Summary written to {csv_path} ({len(rows)} runs)")
    bad = [f"{r['sweep']}/{r['case']}" for r in rows
           if np.isfinite(r['conservation_ratio']) and abs(abs(r['conservation_ratio']) - 1.0) > 0.05]
    if bad:
        print(f"WARNING: bath flux integral and change of content differ by > 5%: {bad}")
    return by_sweep


def plot_all(by_sweep):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not available in this Python. Plot the extracted results with a "
              "Python that has it, e.g. '<venv>/bin/python getoutput_stage2.py --plot-only'")
        return
    print(f"Plotting with matplotlib {matplotlib.__version__} (Python {sys.version.split()[0]})")
    try:                                                # watchdog: stack trace if a plot hangs
        import faulthandler
        faulthandler.dump_traceback_later(120, repeat=True)
    except Exception:
        faulthandler = None
    plot_dir = os.path.join(sa_dir, 'plots')
    os.makedirs(plot_dir, exist_ok=True)
    lstyle = {'slip': '--', 'catch': '-'}

    def save(fig, name):
        fig.tight_layout()
        path = os.path.join(plot_dir, f'{name}.pdf')
        fig.savefig(path)
        plt.close(fig)
        print(f"Plot saved: {path}")

    for sweep, items in sorted(by_sweep.items()):
        items = sorted(items, key=lambda it: (it[0], it[1]))
        xlabel = items[0][2] or 'case'
        xs = sorted(set(it[1] for it in items if np.isfinite(it[1])))
        color = {x: plt.cm.viridis(v) for x, v in zip(xs, np.linspace(0.0, 0.9, max(len(xs), 1)))}

        # Profiles along the column at several hold times: one panel per case
        tcolors = plt.cm.plasma(np.linspace(0.0, 0.85, len(PROFILE_TIMES)))
        profiles = [('c', r'$c/c(0) - 1$', lambda out, i: out['c'][i] / out['c'][0] - 1.0),
                    ('cb', r'$c_{b,tot}/c_{b,tot}(0) - 1$', lambda out, i: out['cb_tot'][i] / out['cb_tot'][0] - 1.0),
                    ('thetaf', r'$\theta_f/\theta_f(0)$', lambda out, i: out['thetaf'][i] / out['thetaf'][0]),
                    ('stress', r'$\sigma_{12}$ (Pa)', lambda out, i: out['s12'][i])]
        ncol = min(len(items), 4)
        nrow = int(np.ceil(len(items) / ncol))
        for name, ylabel, fn in profiles:
            print(f"Plotting {sweep}_profiles_{name} ...", flush=True)
            fig, axs = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 3.0 * nrow), squeeze=False, sharey=True)
            for ax, (group, x, _, out, m) in zip(axs.flat, items):
                for th, tc in zip(PROFILE_TIMES, tcolors):
                    i = frame_at(out, th)
                    ax.plot(out['x_bath'] - out['x'], fn(out, i), color=tc, marker='.', ms=3,
                            label=f"t = {out['t'][i] - out['step_ends'][0]:.0f} s")
                ax.set_title(f'{group}, {xlabel} = {x:g}' if np.isfinite(x) else group, fontsize=8)
                ax.set_xlabel(r'Distance from bath face ($\mu$m)')
                ax.grid(True, alpha=0.3)
            axs.flat[0].legend(fontsize=6, frameon=False, loc='upper right')
            for ax in axs[:, 0]:
                ax.set_ylabel(ylabel)
            for ax in list(axs.flat)[len(items):]:
                ax.set_visible(False)
            fig.suptitle(f'{sweep}: profiles along the column')
            save(fig, f'{sweep}_profiles_{name}')

        # Uptake and average stress vs hold time: one panel, lines coloured by x, styled by law
        print(f"Plotting {sweep}_uptake_stress_time ...", flush=True)
        fig, axs = plt.subplots(1, 2, figsize=(10.0, 3.8))
        for group, x, _, out, m in items:
            hold, th = hold_index(out)
            lab = f'{group}, {xlabel} = {x:g}' if np.isfinite(x) else group
            axs[0].plot(th, m['uptake_t'][hold] / m['total0'], ls=lstyle.get(group, '-'), color=color.get(x, 'k'), label=lab)
            s_avg = out['rf1_top'] / out['area_top']
            axs[1].plot(th, s_avg[hold] / m['sigma_peak'], ls=lstyle.get(group, '-'), color=color.get(x, 'k'), label=lab)
        axs[0].set_ylabel('Net uptake / initial content')
        axs[1].set_ylabel(r'$\bar\sigma_{12}/\bar\sigma_{12,peak}$')
        for ax in axs:
            ax.set_xlabel('Hold time (s)')
            ax.grid(True, alpha=0.3)
        axs[0].legend(fontsize=6, frameon=False, loc='upper left', ncol=2)
        fig.suptitle(f'{sweep}: exchange with the bath and average shear stress')
        save(fig, f'{sweep}_uptake_stress_time')

        # Summary vs the swept parameter, one line per law
        if xs:
            print(f"Plotting {sweep}_summary_vs_x ...", flush=True)
            fig, axs = plt.subplots(1, 3, figsize=(12.0, 3.6))
            for group in sorted(set(it[0] for it in items)):
                pts = sorted((it[1], it[4]) for it in items if it[0] == group and np.isfinite(it[1]))
                xv = [p[0] for p in pts]
                kw = dict(ls=lstyle.get(group, '-'), marker='o', label=group)
                axs[0].plot(xv, [p[1]['uptake_rel_end'] for p in pts], **kw)
                axs[1].plot(xv, [p[1][f'penetration_{PENETRATION_TIMES[0]:g}s'] for p in pts], **kw)
                axs[2].plot(xv, [p[1]['R_end'] for p in pts], **kw)
            axs[0].set_ylabel('Net uptake at end / initial content')
            axs[1].set_ylabel(f'Penetration depth at {PENETRATION_TIMES[0]:g} s' + r' ($\mu$m)')
            axs[2].set_ylabel(r'$\bar\sigma_{12,end}/\bar\sigma_{12,peak}$')
            for ax in axs:
                ax.set_xlabel(xlabel)
                ax.grid(True, alpha=0.3)
                ax.legend(fontsize=7, frameon=False, loc='upper right')
                if sweep == 'B_D' or xlabel.startswith('KEQ'):
                    ax.set_xscale('log')
            fig.suptitle(f'{sweep}: summary vs {xlabel}')
            save(fig, f'{sweep}_summary_vs_x')
    if faulthandler is not None:
        faulthandler.cancel_dump_traceback_later()


# ------------------------------------------------------------------------------------------
if __name__ == '__main__':
    user_args = sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else sys.argv[1:]
    plot_only = '--plot-only' in user_args
    do_plot = plot_only or '--plot' in user_args
    paths = [a for a in user_args if not a.startswith('--')]

    if not plot_only:
        if not IN_ABAQUS:
            sys.exit("Extraction needs Abaqus: 'abaqus python getoutput_stage2.py' "
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
