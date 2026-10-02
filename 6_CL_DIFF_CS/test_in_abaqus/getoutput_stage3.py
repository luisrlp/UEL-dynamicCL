# -*- coding: utf-8 -*-
#
# Outputs of the stage-3 study (run_stage3.py): AFM indentation (spherical indenter, quarter model),
# Equilibrate (indenter fixed) -> Indentation (ramp) -> Hold. Results in
# SA_stage3/<scale>_<chem>/<sweep>/<case>/ (job stage3_afm_<scale>_<chem>).
#
# Usage (from test_in_abaqus/):
#   abaqus python getoutput_stage3.py                          -> extract all runs in SA_stage3/
#   abaqus python getoutput_stage3.py SA_stage3/large_closed/A_core/catch_d0.15   -> extract one run
#   abaqus python getoutput_stage3.py --plot                   -> extract all, then analyse/plot
#   python3 getoutput_stage3.py --plot-only                    -> only analyse/plot (no Abaqus), from
#                                                                 the saved stage3_output.pkl files
#   (abaqus cae noGUI=getoutput_stage3.py -- <args> also works)
#
# Extraction, per run -> <run folder>/stage3_output.pkl:
#   indenter (history output, every increment): th (total time), u2, rf2 of the reference node;
#   per field frame: t; per element (average over the 8 integration points) thetaf (UVARM1),
#   c (UVARM3), cb_tot (UVARM4), s22 (UVARM6); top-surface nodes: u2 (for the contact radius);
#   probe elements (column along the indenter axis): per-direction cb (UVARM17-...) and LE;
#   rfl_bath (sum of RFL over the bath nodes, open faces); element centroids and volumes
#   (reference), indenter radius and gap (from the '** INDENTER:' line of the input), props.
#
# Analysis -> SA_stage3/stage3_summary.csv and SA_stage3/plots/<scale>_<chem>_<sweep>_*.pdf:
#   F = 4 RF2 (full sphere, quarter model; pN with um and Pa), depth delta from the contact onset
#   (first F > 1e-4 F_peak), Hertz shear modulus G_H = 3 F / (16 sqrt(R) delta^1.5) (incompressible
#   half-space) at the end of the indentation; relaxation in the hold: R(10, 50, 200 s) and
#   R_end = F/F_peak, t_half, settled; contact radius a (top nodes on the sphere) vs sqrt(R delta);
#   under the indenter (top element of the probe column): cb_tot, c and thetaf relative to the
#   start of the indentation; Equilibrate: mean cb_tot change (prestress) and max |S22| (should be
#   ~0: the network prestress is isotropic); conservation of the total content (closed faces).
#   Plots per sweep: F vs delta and G_H vs delta; F/F_peak vs hold time; cb_tot and c under the
#   indenter vs hold time; profiles along the axis (end of the hold); per-direction cb_i (top probe
#   element) vs the direction stretch lambda_i at the end of the indentation and of the hold
#   (lambda_i = |V m_i|, V = exp(LE): exact on the symmetry axis, where the rotation vanishes).
#   Scaling check: F/R^2 vs hold time, large scale vs small scale with D (R_small/R_large)^2.

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

from getoutput_stage1 import fortran_directions, read_global_factor, read_run_properties

sa_dir = os.path.join(os.getcwd(), 'SA_stage3')
OUT_PKL = 'stage3_output.pkl'
REF_NODE = 100000                  # indenter reference node (2_uel_afm.py ref_node_id)
N_QUARTER = 4.0                    # quarter model: total force = 4 RF2
CONTACT_FRACTION = 1e-4            # contact onset: F > CONTACT_FRACTION * F_peak
CONTACT_TOL = 2e-3                 # top node in contact: within CONTACT_TOL * R of the sphere
HOLD_TIMES = (10.0, 50.0, 200.0)
PROFILE_TIMES = (0.0, 10.0, 100.0, 600.0)   # hold times of the axis profiles (s)


def odb_job_name(run_dir):
    """Job name of a run folder: stage3_afm_<scale>_<chem> from SA_stage3/<scale>_<chem>/..."""
    group = os.path.relpath(run_dir, sa_dir).replace('\\', '/').split('/')[0]
    return f'stage3_afm_{group}'


def read_indenter(inp_path):
    with open(inp_path) as f:
        for line in f:
            if line.startswith('** INDENTER:'):
                vals = dict(re.findall(r'(\w+) = ([^,\s]+)', line))
                return float(vals['radius']), float(vals['initial_gap'])
    return np.nan, np.nan


def count_warnings(run_dir, job):
    path = os.path.join(run_dir, f'{job}_output.txt')
    text = open(path, errors='ignore').read() if os.path.exists(path) else ''
    return {'accuracy': text.count('accuracy check not met'),
            'kinetics_failed': text.count('kinetics failed'),
            'detF_negative': text.count('detF.lt.zero')}


def job_completed(run_dir, job):
    sta = os.path.join(run_dir, f'{job}.sta')
    return os.path.exists(sta) and 'COMPLETED SUCCESSFULLY' in open(sta, errors='ignore').read()


# ------------------------------------------------------------------------------------------
# Extraction (Abaqus Python)
# ------------------------------------------------------------------------------------------
def element_average(field, region, lookup, n):
    """Per-element average over the integration points (bulk data: fast for large meshes).
    lookup: array element label -> row (-1 if not in the region). Returns (n, ncomp)."""
    sums, counts = None, np.zeros(n)
    for blk in field.getSubset(region=region, position=INTEGRATION_POINT).bulkDataBlocks:
        rows = lookup[np.asarray(blk.elementLabels, dtype=int)]
        data = np.asarray(blk.data, dtype=float).reshape(len(rows), -1)
        keep = rows >= 0
        if sums is None:
            sums = np.zeros((n, data.shape[1]))
        np.add.at(sums, rows[keep], data[keep])
        np.add.at(counts, rows[keep], 1.0)
    return sums / np.maximum(counts, 1.0)[:, None]


def node_values(field, region, lookup, n, comp):
    out = np.full(n, np.nan)
    for blk in field.getSubset(region=region).bulkDataBlocks:
        rows = lookup[np.asarray(blk.nodeLabels, dtype=int)]
        data = np.asarray(blk.data, dtype=float).reshape(len(rows), -1)
        keep = rows >= 0
        out[rows[keep]] = data[keep, comp]
    return out


def extract(odb_files):
    print(f"ODB files to process: {len(odb_files)}")
    n_ok = 0
    for i_file, odb_file in enumerate(odb_files):
        run_dir = os.path.dirname(odb_file)
        job = os.path.splitext(os.path.basename(odb_file))[0]
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
        radius, gap = read_indenter(os.path.join(run_dir, f'{job}.inp'))

        # Flat input: one automatic instance with the UEL nodes, dummy and extra elements
        inst = list(odb.rootAssembly.instances.values())[0]
        coords = {n.label: np.array(n.coordinates) for n in inst.nodes}
        dummy = inst.elementSets['DUMMY_MESH']
        elems = list(dummy.elements)
        lookup = -np.ones(max(e.label for e in elems) + 1, dtype=int)
        xc, vol = np.zeros((len(elems), 3)), np.zeros(len(elems))
        for k, e in enumerate(elems):
            lookup[e.label] = k
            p = np.array([coords[n] for n in e.connectivity])
            xc[k] = p.mean(axis=0)
            vol[k] = np.prod(p.max(axis=0) - p.min(axis=0))       # box elements
        probe_set = inst.elementSets['PROBE_MESH'] if 'PROBE_MESH' in inst.elementSets.keys() else None
        probe_labels = sorted((e.label for e in probe_set.elements), key=lambda l: -xc[lookup[l], 1]) if probe_set else []
        probe_lookup = -np.ones(len(lookup), dtype=int)
        for k, l in enumerate(probe_labels):                       # sorted from the top surface down
            probe_lookup[l] = k
        top = inst.nodeSets['TOP_NODES']
        top_labels = [n.label for n in top.nodes]
        top_lookup = -np.ones(max(coords) + 1, dtype=int)
        for k, l in enumerate(top_labels):
            top_lookup[l] = k
        top_xyz = np.array([coords[l] for l in top_labels])
        bath = inst.nodeSets['OPEN_CHEM_NODES'] if 'OPEN_CHEM_NODES' in inst.nodeSets.keys() else None

        steps = list(odb.steps.values())
        # Indenter history (every increment)
        th, u2, rf2 = [], [], []
        for step in steps:
            for key, region in step.historyRegions.items():
                if key.endswith(f'.{REF_NODE}') and 'RF2' in region.historyOutputs.keys():
                    ho = region.historyOutputs
                    th += [step.totalTime + d[0] for d in ho['RF2'].data]
                    rf2 += [d[1] for d in ho['RF2'].data]
                    u2 += [d[1] for d in ho['U2'].data]
        # Field frames
        t, thetaf, c, cb_tot, s22, top_u2, rfl_bath, probe_cb, probe_le = [], [], [], [], [], [], [], [], []
        ndir = None
        for step in steps:
            for frame in step.frames:
                fo = frame.fieldOutputs
                if 'UVARM1' not in fo.keys():
                    continue
                t.append(step.totalTime + frame.frameValue)
                n = len(elems)
                thetaf.append(element_average(fo['UVARM1'], dummy, lookup, n)[:, 0])
                c.append(element_average(fo['UVARM3'], dummy, lookup, n)[:, 0])
                cb_tot.append(element_average(fo['UVARM4'], dummy, lookup, n)[:, 0])
                s22.append(element_average(fo['UVARM6'], dummy, lookup, n)[:, 0])
                top_u2.append(node_values(fo['U'], top, top_lookup, len(top_labels), 1))
                rfl_bath.append(float(sum(np.sum(np.asarray(b.data)) for b in fo['RFL'].getSubset(region=bath).bulkDataBlocks))
                                if bath is not None and 'RFL' in fo.keys() else np.nan)
                if probe_labels:
                    if ndir is None:
                        ndir = len([k for k in fo.keys() if re.match(r'UVARM\d+$', k) and int(k[5:]) >= 17])
                    np_ = len(probe_labels)
                    probe_cb.append(np.hstack([element_average(fo[f'UVARM{17 + j}'], probe_set, probe_lookup, np_)
                                               for j in range(ndir)]))
                    probe_le.append(element_average(fo['LE'], probe_set, probe_lookup, np_))
        # Copy the step data before closing: ODB objects are invalid after odb.close()
        step_names = [s.name for s in steps]
        step_ends = [s.totalTime + s.timePeriod for s in steps]
        odb.close()

        props, label = read_run_properties(run_dir)
        output = {'th': np.array(th), 'u2': np.array(u2), 'rf2': np.array(rf2),
                  't': np.array(t), 'step_names': step_names, 'step_ends': step_ends,
                  'xc': xc, 'vol': vol, 'thetaf': np.array(thetaf), 'c': np.array(c),
                  'cb_tot': np.array(cb_tot), 's22': np.array(s22),
                  'top_xyz': top_xyz, 'top_u2': np.array(top_u2), 'rfl_bath': np.array(rfl_bath),
                  'xc_axis': (coords[REF_NODE][[0, 2]] if REF_NODE in coords
                              else top_xyz[np.argmax(top_xyz[:, 0] + top_xyz[:, 2]), [0, 2]]),
                  'probe_y': np.array([xc[lookup[l], 1] for l in probe_labels]),
                  'probe_cb': np.array(probe_cb), 'probe_le': np.array(probe_le),
                  'probe_top_row': int(lookup[probe_labels[0]]) if probe_labels else -1,
                  'radius': radius, 'gap': gap, 'props': props, 'label': label, 'folder': folder,
                  'odb': odb_file, 'warnings': count_warnings(run_dir, job),
                  'completed': job_completed(run_dir, job)}
        with open(os.path.join(run_dir, OUT_PKL), 'wb') as f:
            pickle.dump(output, f)
        print(f"  Saved {os.path.join(run_dir, OUT_PKL)} ({len(th)} increments, {len(t)} frames, "
              f"{len(elems)} elements, {len(probe_labels)} probe elements)")
        n_ok += 1
    print(f"Done: {n_ok}/{len(odb_files)} ODB files processed.")


# ------------------------------------------------------------------------------------------
# Analysis and plots (plain Python)
# ------------------------------------------------------------------------------------------
def step_end(out, name, default=0.0):
    return out['step_ends'][out['step_names'].index(name)] if name in out['step_names'] else default


def indenter_force(out):
    """Total force F = 4 RF2 (sign: > 0 when pushing into the gel) and the indenter displacement."""
    f = N_QUARTER * out['rf2']
    t_ind = step_end(out, 'Indentation')
    i_end = np.argmin(np.abs(out['th'] - t_ind))
    if f[i_end] < 0:
        f = -f
    return f, -out['u2'], i_end


def frame_index(out, t_abs):
    return int(np.argmin(np.abs(out['t'] - t_abs)))


def contact_radius(out, i, u_ind):
    """Largest distance from the axis of the top nodes lying on the sphere (frame i)."""
    R = out['radius']
    xyz = out['top_xyz']
    axis = out['xc_axis']                                               # indenter axis (x, z)
    r = np.hypot(xyz[:, 0] - axis[0], xyz[:, 2] - axis[1])
    y_top = xyz[:, 1].max()
    y_tip = y_top + out['gap'] - u_ind                                  # current sphere tip
    y_sphere = y_tip + R - np.sqrt(np.maximum(R ** 2 - r ** 2, 0.0))   # sphere surface above r
    y_node = xyz[:, 1] + out['top_u2'][i]
    on = (r < R) & (y_node >= y_sphere - CONTACT_TOL * R)
    return float(r[on].max()) if on.any() else 0.0


def direction_stretch(le, dirs):
    """lambda_i = |V m_i| with V = exp(LE) (LE components 11 22 33 12 13 23)."""
    L = np.array([[le[0], le[3], le[4]], [le[3], le[1], le[5]], [le[4], le[5], le[2]]])
    w, Q = np.linalg.eigh(L)
    V = Q @ np.diag(np.exp(w)) @ Q.T
    return np.linalg.norm(dirs @ V.T, axis=1)


def case_metrics(out, dirs):
    m = {}
    R, gap = out['radius'], out['gap']
    F, uind, i_end = indenter_force(out)
    th = out['th']
    t_eq, t_ind = step_end(out, 'Equilibrate'), step_end(out, 'Indentation')
    load = (th > t_eq + 1e-12) & (th <= t_ind + 1e-9)
    hold = th >= t_ind - 1e-9
    F_peak = F[i_end]
    m['F_peak'] = F_peak
    # Contact onset and depth from contact
    # (F^(2/3) is linear in the depth for Hertz contact: extrapolated to zero from the first two
    # increments in contact, so that the onset is not biased by the increment size)
    idx = np.nonzero(load & (F > CONTACT_FRACTION * abs(F_peak)))[0]
    u_c = uind[idx[0]] if len(idx) else gap
    if len(idx) > 1:
        i1, i2 = idx[0], idx[1]
        q1, q2 = F[i1] ** (2.0 / 3.0), F[i2] ** (2.0 / 3.0)
        if q2 > q1 and uind[i2] > uind[i1]:
            u_c = max(uind[i1] - q1 * (uind[i2] - uind[i1]) / (q2 - q1), uind[i1 - 1] if i1 > 0 else -np.inf)
    m['u_contact'], m['gap'] = u_c, gap
    delta = uind - u_c
    m['delta_end'] = delta[i_end]
    m['delta_R'] = delta[i_end] / R
    with np.errstate(divide='ignore', invalid='ignore'):
        m['G_H'] = 3.0 * F_peak / (16.0 * np.sqrt(R) * delta[i_end] ** 1.5) if delta[i_end] > 0 else np.nan
    m['F_over_R2'] = F_peak / R ** 2
    # Relaxation in the hold (history output, every increment)
    t_h, F_h = th[hold] - t_ind, F[hold]
    F_end = F_h[-1]
    m['F_end'], m['R_end'] = F_end, F_end / F_peak if F_peak else np.nan
    for tt in HOLD_TIMES:
        m[f'R{tt:g}'] = np.interp(tt, t_h, F_h) / F_peak if t_h[-1] >= tt and F_peak else np.nan
    dF = F_peak - F_end
    m['t_half'] = np.nan
    if abs(dF) > 1e-3 * abs(F_peak):
        s = (F_h - F_end) / dF
        k = np.nonzero(s <= 0.5)[0]
        if len(k) and k[0] > 0:
            m['t_half'] = np.interp(0.5, [s[k[0]], s[k[0] - 1]], [t_h[k[0]], t_h[k[0] - 1]])
    t90 = 0.9 * t_h[-1]
    m['settled'] = bool(abs(F_end - np.interp(t90, t_h, F_h)) <= 0.05 * abs(dF)) if dF else True
    # Field frames: start of the indentation (end of Equilibrate), end of indentation, end of hold
    i0, i_eq, i_in, i_ho = 0, frame_index(out, t_eq), frame_index(out, t_ind), len(out['t']) - 1
    cb, c, tf, V = out['cb_tot'], out['c'], out['thetaf'], out['vol']
    m['equil_dcb_mean'] = float(np.sum(cb[i_eq] * V) / np.sum(cb[i0] * V) - 1.0)
    m['equil_s22_max'] = float(np.max(np.abs(out['s22'][i_eq])))
    k = out['probe_top_row']
    if k >= 0:
        m['dcb_under_end_indent'] = cb[i_in, k] / cb[i_eq, k] - 1.0
        m['dcb_under_end_hold'] = cb[i_ho, k] / cb[i_eq, k] - 1.0
        m['dc_under_end_hold'] = c[i_ho, k] / c[i_eq, k] - 1.0
        m['thetaf_under_end_hold'] = tf[i_ho, k] / tf[i_eq, k]
    total = c @ V
    m['content_drift_max'] = float(np.max(np.abs(total / total[0] - 1.0)))
    m['uptake_rel_end'] = float(total[-1] / total[0] - 1.0)
    # Contact radius at the end of the indentation vs Hertz sqrt(R delta)
    u_in = np.interp(out['t'][i_in], th, uind)
    m['a_end'] = contact_radius(out, i_in, u_in)
    m['a_hertz'] = np.sqrt(R * m['delta_end']) if m['delta_end'] > 0 else np.nan
    # Per-direction state of the top probe element
    if len(out['probe_cb']):
        nd = min(len(dirs), out['probe_cb'].shape[2])
        m['lambda_end_indent'] = direction_stretch(out['probe_le'][i_in, 0], dirs[:nd])
        m['lambda_end_hold'] = direction_stretch(out['probe_le'][i_ho, 0], dirs[:nd])
        ref = out['probe_cb'][i_eq, 0, :nd]
        m['cbi_rel_end_indent'] = out['probe_cb'][i_in, 0, :nd] / ref
        m['cbi_rel_end_hold'] = out['probe_cb'][i_ho, 0, :nd] / ref
    # Series for the plots
    m['delta_t'], m['F_t'], m['load'], m['hold_t'], m['hold_F'] = delta, F, load, t_h, F_h
    return m


SWEEP_X = {'d': ('delta/R', None), 'lambda0': ('LAMBDA0', 'LAMBDA0'), 'D': (r'D ($\mu$m$^2$/s)', 'D')}


def group_and_x(case, props, radius):
    """Law and swept value of a case folder '<law>_<var><value>' (or '<law>_nomatrix')."""
    law, _, rest = case.partition('_')
    for prefix in sorted(SWEEP_X, key=len, reverse=True):     # 'lambda0' before shorter prefixes
        value = rest[len(prefix):]
        if rest.startswith(prefix) and re.fullmatch(r'\d[-+.\deE]*', value):
            xlabel, key = SWEEP_X[prefix]
            return law, (float(value) if key is None else props.get(key, float(value))), xlabel
    return law, np.nan, ''


def load_outputs():
    outs = {}
    for path in sorted(glob.glob(os.path.join(sa_dir, '*', '*', '*', OUT_PKL))):
        with open(path, 'rb') as f:
            out = pickle.load(f)
        out['folder'] = os.path.relpath(os.path.dirname(path), sa_dir).replace('\\', '/')
        outs[out['folder']] = out
    return outs


SKIP_CSV = ('delta_t', 'F_t', 'load', 'hold_t', 'hold_F', 'lambda_end_indent', 'lambda_end_hold',
            'cbi_rel_end_indent', 'cbi_rel_end_hold')


def analyse(outs, dirs):
    rows, by_sweep, metrics = [], {}, {}
    for folder, out in outs.items():
        group, sweep, case = folder.split('/')[:3]
        try:
            m = case_metrics(out, dirs)
        except Exception as e:                      # e.g. an incomplete run without a hold
            print(f"WARNING: metrics of {folder} failed ({e}); skipped")
            continue
        law, x, xlabel = group_and_x(case, out['props'], out['radius'])
        p = out['props']
        rows.append({'group': group, 'sweep': sweep, 'case': case, 'label': out['label'], 'law': law, 'x': x,
                     'R': out['radius'], **{k: p.get(k, np.nan) for k in ('KOFF0', 'KCATCH0', 'LAMBDA0', 'D', 'PHINET', 'C10')},
                     **{k: v for k, v in m.items() if k not in SKIP_CSV},
                     **{f'warn_{k}': v for k, v in out['warnings'].items()}, 'completed': out['completed']})
        by_sweep.setdefault((group, sweep), []).append((law, x, xlabel, out, m))
        metrics[folder] = (out, m)
    if not rows:
        return by_sweep, metrics
    csv_path = os.path.join(sa_dir, 'stage3_summary.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in sorted(rows, key=lambda r: (r['group'], r['sweep'], r['law'], r['x'])):
            writer.writerow({k: (f'{v:.6g}' if isinstance(v, float) else v) for k, v in r.items()})
    print(f"Summary written to {csv_path} ({len(rows)} runs)")
    drift = [f"{r['group']}/{r['sweep']}/{r['case']}" for r in rows
             if r['group'].endswith('closed') and r['content_drift_max'] > 1e-3]
    if drift:
        print(f"WARNING: total crosslinker content not conserved (> 0.1%) in closed runs: {drift}")
    eq = [f"{r['group']}/{r['sweep']}/{r['case']}" for r in rows if r['equil_s22_max'] > 1e-3 * abs(r['F_over_R2'])]
    if eq:
        print(f"NOTE: stress at the end of Equilibrate > 0.1% of F/R^2 (prestress not isotropic?): {eq}")
    return by_sweep, metrics


def plot_all(by_sweep, metrics):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not available in this Python. Plot the extracted results with a "
              "Python that has it, e.g. '<venv>/bin/python getoutput_stage3.py --plot-only'")
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

    for (group, sweep), items in sorted(by_sweep.items()):
        items = sorted(items, key=lambda it: (it[0], it[1]))
        xlabel = items[0][2] or 'case'
        xs = sorted(set(it[1] for it in items if np.isfinite(it[1])))
        color = {x: plt.cm.viridis(v) for x, v in zip(xs, np.linspace(0.0, 0.9, max(len(xs), 1)))}
        name = f'{group}_{sweep}'

        def lab(law, x):
            return f'{law}, {xlabel} = {x:g}' if np.isfinite(x) else law

        print(f"Plotting {name} ...", flush=True)
        # Force-depth (indentation), Hertz modulus vs depth, relaxation, cb under the indenter
        fig, axs = plt.subplots(2, 2, figsize=(10.0, 7.2))
        for law, x, _, out, m in items:
            kw = dict(ls=lstyle.get(law, '-'), color=color.get(x, 'k'), label=lab(law, x))
            ld = m['load'] & (m['delta_t'] > 0)
            d, F, R = m['delta_t'][ld], m['F_t'][ld], out['radius']
            axs[0, 0].plot(d / R, F / R ** 2, **kw)
            axs[0, 1].plot(d / R, 3.0 * F / (16.0 * np.sqrt(R) * d ** 1.5), **kw)
            axs[1, 0].plot(m['hold_t'], m['hold_F'] / m['F_peak'], **kw)
            k = out['probe_top_row']
            if k >= 0:
                t_eq, t_ind = step_end(out, 'Equilibrate'), step_end(out, 'Indentation')
                i_eq = frame_index(out, t_eq)
                sel = out['t'] >= t_eq - 1e-9
                axs[1, 1].plot(out['t'][sel] - t_ind, out['cb_tot'][sel, k] / out['cb_tot'][i_eq, k] - 1.0, **kw)
        axs[0, 0].set(xlabel=r'$\delta/R$', ylabel=r'$F/R^2$ (Pa)')
        axs[0, 1].set(xlabel=r'$\delta/R$', ylabel=r'Hertz $G_H = 3F/(16\sqrt{R}\delta^{3/2})$ (Pa)')
        axs[1, 0].set(xlabel='Hold time (s)', ylabel=r'$F/F_{peak}$', xscale='symlog')
        axs[1, 1].set(xlabel='Time from the end of the indentation (s)',
                      ylabel=r'$c_{b,tot}/c_{b,tot}(t_0) - 1$ under the indenter', xscale='symlog')
        for ax in axs.flat:
            ax.grid(True, alpha=0.3)
        axs[0, 0].legend(fontsize=6, frameon=False, loc='upper left')
        fig.suptitle(f'{group} / {sweep}: force, Hertz modulus, relaxation, binding under the indenter')
        save(fig, f'{name}_force_relax')

        # Profiles along the indenter axis (top probe column = elements on the axis) at hold times
        prof = [it for it in items if it[3]['probe_top_row'] >= 0]
        if prof:
            ncol = min(len(prof), 4)
            nrow = int(np.ceil(len(prof) / ncol))
            tcolors = plt.cm.plasma(np.linspace(0.0, 0.85, len(PROFILE_TIMES)))
            fig, axs = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 3.0 * nrow), squeeze=False, sharey=True)
            for ax, (law, x, _, out, m) in zip(axs.flat, prof):
                R = out['radius']
                xc = out['xc']
                axis_xz = xc[out['probe_top_row'], [0, 2]]
                col = np.nonzero(np.hypot(xc[:, 0] - axis_xz[0], xc[:, 2] - axis_xz[1]) < 1e-6 * R + 1e-12)[0]
                col = col[np.argsort(-xc[col, 1])]
                depth = (xc[:, 1].max() - xc[col, 1]) / R
                t_eq, t_ind = step_end(out, 'Equilibrate'), step_end(out, 'Indentation')
                i_eq = frame_index(out, t_eq)
                for tt, tc in zip(PROFILE_TIMES, tcolors):
                    i = frame_index(out, t_ind + tt)
                    ax.plot(depth, out['cb_tot'][i, col] / out['cb_tot'][i_eq, col] - 1.0, color=tc, marker='.',
                            ms=3, label=f"t = {out['t'][i] - t_ind:.0f} s")
                ax.set_title(lab(law, x), fontsize=8)
                ax.set_xlabel('Depth below the surface / R')
                ax.grid(True, alpha=0.3)
            axs.flat[0].legend(fontsize=6, frameon=False, loc='lower right')
            for ax in axs[:, 0]:
                ax.set_ylabel(r'$c_{b,tot}/c_{b,tot}(t_0) - 1$')
            for ax in list(axs.flat)[len(prof):]:
                ax.set_visible(False)
            fig.suptitle(f'{group} / {sweep}: bound crosslinkers along the indenter axis')
            save(fig, f'{name}_axis_profiles')

            # Per-direction binding vs direction stretch (top probe element)
            fig, axs = plt.subplots(1, 2, figsize=(10.0, 3.8), sharey=True)
            for law, x, _, out, m in prof:
                if 'lambda_end_hold' not in m:
                    continue
                kw = dict(marker='o' if law == 'catch' else 'x', ls='none', ms=3, color=color.get(x, 'k'), label=lab(law, x))
                axs[0].plot(m['lambda_end_indent'], m['cbi_rel_end_indent'], **kw)
                axs[1].plot(m['lambda_end_hold'], m['cbi_rel_end_hold'], **kw)
            axs[0].set(xlabel=r'$\lambda_i$', ylabel=r'$c_{b,i}/c_{b,i}(t_0)$', title='End of the indentation')
            axs[1].set(xlabel=r'$\lambda_i$', title='End of the hold')
            for ax in axs:
                ax.grid(True, alpha=0.3)
            axs[0].legend(fontsize=6, frameon=False, loc='best')
            fig.suptitle(f'{group} / {sweep}: per-direction binding under the indenter (top axis element)')
            save(fig, f'{name}_directions')

        # Summary vs the swept parameter
        if len(xs) > 1:
            fig, axs = plt.subplots(1, 3, figsize=(12.0, 3.6))
            for law in sorted(set(it[0] for it in items)):
                pts = sorted((it[1], it[4]) for it in items if it[0] == law and np.isfinite(it[1]))
                xv = [p[0] for p in pts]
                kw = dict(ls=lstyle.get(law, '-'), marker='o', label=law)
                axs[0].plot(xv, [p[1]['G_H'] for p in pts], **kw)
                axs[1].plot(xv, [p[1]['R_end'] for p in pts], **kw)
                axs[2].plot(xv, [p[1]['t_half'] for p in pts], **kw)
            axs[0].set_ylabel(r'$G_H$ at the end of the indentation (Pa)')
            axs[1].set_ylabel(r'$F_{end}/F_{peak}$')
            axs[2].set_ylabel(r'$t_{1/2}$ (s)')
            for ax in axs:
                ax.set_xlabel(xlabel)
                ax.grid(True, alpha=0.3)
                ax.legend(fontsize=7, frameon=False, loc='best')
            if sweep.startswith('B'):
                axs[0].set_yscale('log')
            fig.suptitle(f'{group} / {sweep}: summary vs {xlabel}')
            save(fig, f'{name}_summary')

    # Scaling check: large (D) vs small (D (R_s/R_l)^2) and small (D), same delta/R
    pairs = [('large_closed/A_core/{law}_d0.15', 'Large, D'), ('small_closed/D_scaling/{law}_D*', 'Small, scaled D'),
             ('small_closed/A_core/{law}_d0.15', 'Small, D')]
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    found = False
    for law in ('slip', 'catch'):
        ref = None
        for pattern, name in pairs:
            keys = [k for k in metrics if re.fullmatch(pattern.format(law=law).replace('.', r'\.').replace('*', '.*'), k)]
            if not keys:
                continue
            out, m = metrics[keys[0]]
            found = True
            y = m['hold_F'] / out['radius'] ** 2
            ax.plot(m['hold_t'], y, ls=lstyle[law], label=f'{law}: {name}')
            if ref is None:
                ref = (m['hold_t'], y)
            elif name == 'Small, scaled D':
                diff = np.max(np.abs(np.interp(ref[0], m['hold_t'], y) - ref[1])) / np.max(np.abs(ref[1]))
                print(f"Scaling check ({law}): max |F/R^2 (small, scaled D) - F/R^2 (large)| / max = {diff:.3g}")
    if found:
        ax.set(xlabel='Hold time (s)', ylabel=r'$F/R^2$ (Pa)', xscale='symlog')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, frameon=False)
        fig.suptitle(r'Scaling check: $\delta/R = 0.15$')
        save(fig, 'scaling_check')
    else:
        plt.close(fig)
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
            sys.exit("Extraction needs Abaqus: 'abaqus python getoutput_stage3.py' "
                     "(or use --plot-only to analyse existing results)")
        if paths:
            target = os.path.abspath(paths[0])
            odb_files = [target] if target.endswith('.odb') else [os.path.join(target, odb_job_name(target) + '.odb')]
        else:
            odb_files = sorted(glob.glob(os.path.join(sa_dir, '*', '*', '*', 'stage3_afm_*.odb')))
        extract(odb_files)

    if do_plot:
        outs = load_outputs()
        if not outs:
            sys.exit(f"No {OUT_PKL} files found in {sa_dir}")
        plot_all(*analyse(outs, fortran_directions(read_global_factor())))
