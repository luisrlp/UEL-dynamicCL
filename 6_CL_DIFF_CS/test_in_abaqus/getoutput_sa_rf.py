# -*- coding: utf-8 -*-
# (the old replay file used 'mbcs', a Windows-only encoding that fails on Linux)
#
# Reaction force at the indenter reference point (node set IndenterRef) for the
# sensitivity analysis runs created by run_sa.py (based on getoutput_gel_final.py),
# and plots of reaction force vs indentation for each parameter studied.
#
# Usage (from test_in_abaqus/):
#   abaqus cae noGUI=getoutput_sa_rf.py                           -> extract all runs in SA/
#   abaqus cae noGUI=getoutput_sa_rf.py -- SA/DX/DX_0.001          -> extract one run folder
#   abaqus cae noGUI=getoutput_sa_rf.py -- SA/DX/DX_0.001/cube_indent_uel_relax_auto.odb
#   abaqus cae noGUI=getoutput_sa_rf.py -- --plot                  -> extract all, then plot
#   python3 getoutput_sa_rf.py --plot-only                        -> only plot (no Abaqus),
#                                                                    from the saved rf_indenter.pkl
#
# For each run, the indenter time history is saved to <run folder>/rf_indenter.pkl:
#   {'t': total time, 'u2': indenter displacement U2, 'rf2': reaction force RF2,
#    'folder': run folder, 'odb': odb path}
# and all runs together to SA/rf_indenter_all.pkl as {run folder: dict above}.
# Plots: SA/plots/RF_vs_indentation_<VAR>.pdf, one line per value of <VAR>.
# Note: RF2 is the raw reaction at the reference point (sign as in Abaqus), and the model is a
# quarter of the sample (XSYMM and ZSYMM planes through the indenter axis): see FORCE_FACTOR.

import os
import sys
import glob
import pickle
import numpy as np

try:
    from abaqus import *
    from abaqusConstants import *
    IN_ABAQUS = True
except ImportError:
    IN_ABAQUS = False

if IN_ABAQUS:
    session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=243.416656494141,
        height=165.277770996094)
    session.viewports['Viewport: 1'].makeCurrent()
    session.viewports['Viewport: 1'].maximize()
    from viewerModules import *
    from driverUtils import executeOnCaeStartup
    executeOnCaeStartup()

# ODB file name in every run folder
file = 'cube_indent_uel_relax_auto.odb'
# Sensitivity analysis folder and node set of the indenter reference point
sa_dir = os.path.join(os.getcwd(), 'SA')
ref_set = 'INDENTERREF'

# Plot settings: force plotted = FORCE_FACTOR * RF2 (e.g. -1 to flip the sign,
# 4 or -4 for the full sample instead of the quarter model); indentation depth = -U2
FORCE_FACTOR = 1.0
FORCE_LABEL = r'Reaction force at indenter, $RF_2$ (pN)'
DEPTH_LABEL = r'Indentation depth, $-U_2$ ($\mu$m)'


def extract(odb_files):
    """Extract t, U2 and RF2 at the indenter reference point from each ODB (Abaqus only)."""
    print(f"ODB files to process: {len(odb_files)}")
    all_outputs = {}
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
            o1 = session.openOdb(name=odb_file, readOnly=True)
        except Exception as e:
            print(f"  Could not open ODB, skipped: {e}")
            continue
        session.viewports['Viewport: 1'].setValues(displayedObject=o1)
        odb = session.odbs[odb_file]

        ###########################     GET NODES IN SET     ############################
        # Flat input file: the set belongs to the automatic part instance (e.g. PART-1-1)
        nodeLabelsTuple = ()
        for inst_name, instance in odb.rootAssembly.instances.items():
            if ref_set in instance.nodeSets.keys():
                for node in instance.nodeSets[ref_set].nodes:
                    nodeLabelsTuple += ((inst_name, (str(node.label), )), )
        if not nodeLabelsTuple:
            print(f"  Node set {ref_set} not found, skipped")
            o1.close()
            continue
        ref_label = nodeLabelsTuple[0][1][0]
        # print('Node labels tuple:', nodeLabelsTuple)
        #################################################################################

        # Without the next step the cycle does not work.
        # It always gets the output from the first file
        # **Clear xyDataObjects before creating new ones**
        for key in list(session.xyDataObjects.keys()):
            del session.xyDataObjects[key]

        session.xyDataListFromField(odb=odb,
                                    outputPosition=NODAL,
                                    variable=(
                                        ('U', NODAL, ((COMPONENT, 'U2'), )),
                                        ('RF', NODAL, ((COMPONENT, 'RF2'), )), ),
                                    nodeLabels=nodeLabelsTuple,
                                    )

        # Data organized as: t, U2, RF2 (over all steps: Indentation, Hold, Withdraw, Relax)
        def find_xy(variable):
            for key in session.xyDataObjects.keys():
                if key.startswith(variable) and key.endswith(f'N: {ref_label}'):
                    return session.xyDataObjects[key]
            raise KeyError(f"No xyData for {variable} at node {ref_label}: "
                           f"{list(session.xyDataObjects.keys())}")

        u2 = find_xy('U:U2')
        rf2 = find_xy('RF:RF2')
        output = {'t': np.array([pt[0] for pt in u2.data]),
                  'u2': np.array([pt[1] for pt in u2.data]),
                  'rf2': np.array([pt[1] for pt in rf2.data]),
                  'folder': folder,
                  'odb': odb_file}

        output_path = os.path.join(run_dir, 'rf_indenter.pkl')
        print(f"  Saving output to: {output_path}")
        with open(output_path, 'wb') as f:
            pickle.dump(output, f)
        all_outputs[folder] = output

        # **Close ODB after extraction**
        o1.close()
    print(f"Done: {len(all_outputs)}/{len(odb_files)} ODB files processed.")
    return all_outputs


def load_outputs():
    """Load the rf_indenter.pkl of every run folder in SA/ -> {run folder: output dict}."""
    outputs = {}
    for path in sorted(glob.glob(os.path.join(sa_dir, '*', 'rf_indenter.pkl'))
                       + glob.glob(os.path.join(sa_dir, '*', '*', 'rf_indenter.pkl'))):
        with open(path, 'rb') as f:
            outputs[os.path.relpath(os.path.dirname(path), sa_dir)] = pickle.load(f)
    return outputs


def plot_sa(outputs):
    """One PDF per studied parameter: reaction force vs indentation depth, one line per value."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not available in this Python. Plot the extracted results with a "
              "Python that has it, e.g. '<venv>/bin/python getoutput_sa_rf.py --plot-only'")
        return

    # Group the runs by parameter, from the folder names <VAR>/<VAR>_<value> (base excluded)
    groups = {}
    for folder, out in outputs.items():
        parts = folder.replace('\\', '/').split('/')
        if len(parts) != 2 or not parts[1].startswith(parts[0] + '_'):
            continue
        var, value = parts[0], float(parts[1][len(parts[0]) + 1:])
        groups.setdefault(var, []).append((value, out))

    plot_dir = os.path.join(sa_dir, 'plots')
    os.makedirs(plot_dir, exist_ok=True)
    for var, runs in sorted(groups.items()):
        runs.sort(key=lambda r: r[0])
        colors = plt.cm.viridis(np.linspace(0.0, 0.9, len(runs)))
        fig, ax = plt.subplots(figsize=(5.0, 3.8))
        for (value, out), color in zip(runs, colors):
            ax.plot(-out['u2'], FORCE_FACTOR * out['rf2'], color=color, lw=1.5,
                    label=f'{var} = {value:g}')
        ax.set_xlabel(DEPTH_LABEL)
        ax.set_ylabel(FORCE_LABEL)
        ax.set_title(f'Sensitivity to {var}')
        ax.grid(True, alpha=0.3)
        ax.legend(frameon=False)
        fig.tight_layout()
        pdf_path = os.path.join(plot_dir, f'RF_vs_indentation_{var}.pdf')
        fig.savefig(pdf_path)
        plt.close(fig)
        print(f"Plot saved: {pdf_path} ({len(runs)} lines)")


###########################     ODB FILES TO PROCESS     #########################
# Inside Abaqus, script arguments come after '--'; with plain Python, after the script name
user_args = (sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else []) if IN_ABAQUS else sys.argv[1:]
plot_only = '--plot-only' in user_args
do_plot = plot_only or '--plot' in user_args
paths = [a for a in user_args if not a.startswith('--')]

if not plot_only:
    if not IN_ABAQUS:
        sys.exit("Extraction needs Abaqus: 'abaqus cae noGUI=getoutput_sa_rf.py' "
                 "(or use --plot-only to plot existing results)")
    if paths:
        target = os.path.abspath(paths[0])
        odb_files = [target] if target.endswith('.odb') else [os.path.join(target, file)]
    else:
        odb_files = sorted(glob.glob(os.path.join(sa_dir, '*', file))
                           + glob.glob(os.path.join(sa_dir, '*', '*', file)))
    all_outputs = extract(odb_files)
    # All runs together (only when processing the whole sensitivity analysis)
    if not paths and all_outputs:
        all_path = os.path.join(sa_dir, 'rf_indenter_all.pkl')
        print(f"Saving all outputs to: {all_path}")
        with open(all_path, 'wb') as f:
            pickle.dump(all_outputs, f)
#################################################################################

if do_plot:
    plot_sa(load_outputs())
