"""One-at-a-time sensitivity analysis of the UEL model in Abaqus (indentation study, stage 3).

The machinery (properties, run folders, skip logic, parallel runs, log, summary) is in
run_study(), which other study scripts (e.g. run_stage1.py) import with their own job, files,
results folder and list of cases.

For every variable in STUDY_PROPS_INFO and every value swept for it, a run folder
SA/<VAR>/<VAR>_<value>/ is created with:
  - properties.inp : base properties with only <VAR> changed, and INITMU recomputed for
                     that set of properties (initial chemical equilibrium, as get_initmu.py)
  - the files needed to run the job: JOB_NAME.inp, sec_uel_cube.inp, prefdir.inp (read by
    uexternaldb at start-up), uel.f90, abaqus_v6.env, aba_param.inc
All jobs are then run in parallel (MAX_WORKERS at a time), each in its own folder.

Base properties are read from ./properties.inp (single source of truth), with optional
overrides in BASE_OVERRIDES. A summary of all cases is written to SA/sa_cases.csv, and every
message (cases skipped/prepared, start and end of each job, jobs running at that moment) is
printed and appended with a timestamp to SA/output.txt.

Usage (from test_in_abaqus/):
  python3 run_sa.py                  # create all folders and run all jobs
  python3 run_sa.py --dry-run        # only create the folders (no Abaqus)
  python3 run_sa.py --only DX KEQ    # restrict to some variables
  python3 run_sa.py --workers 8      # number of simultaneous Abaqus jobs
  python3 run_sa.py --force          # rerun jobs that already completed
Completed jobs (.sta containing "COMPLETED SUCCESSFULLY") whose inputs (properties.inp and the
copied files) are identical to the current ones are skipped and their folders left untouched,
so the script can be relaunched to finish or retry a study. A completed job whose inputs
changed (e.g. new base properties or uel.f90) is rerun.
"""
import argparse
import ast
import concurrent.futures
import csv
import filecmp
import math
import operator
import os
import shutil
import subprocess
import threading
import time
from datetime import datetime

# ------------------------------------------------------------------------------------------
# Study definition
# ------------------------------------------------------------------------------------------
JOB_NAME = 'cube_indent_uel_relax_auto'
SA_DIR_NAME = 'SA'
MAX_WORKERS = 4        # simultaneous Abaqus jobs (limited by CPUs and license tokens)
CPUS_PER_JOB = 4
RUN_BASE_CASE = True   # also run the unmodified base properties in SA/base/

# Changes to the base properties read from properties.inp (applied to every case)
BASE_OVERRIDES = {
    # 'KCATCH0': 0.04,
}

# Variables and values to sweep (one variable at a time)
STUDY_PROPS_INFO = {
    'ETAC': [1. / 3., 2. / 3., 1.0],
    'DX': [1.e-3, 5.e-3, 1.e-2, 5.e-2],
    'DXC': [1.e-3, 5.e-3, 1.e-2, 5.e-2],
    'CACTIN': [9.5e-3 * f for f in (1, 2, 4)],
    'R': [0.1 * f for f in (0.2, 0.5, 1, 2)],
    'RFMAX': [0.25 * f for f in (1, 2, 4)],
    'CHI': [0.1 * f for f in (1, 2, 4)],
    'D': [1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1],
    'KOFF0': [0.01, 0.1, 1.0],
    'KCATCH0': [0.01, 0.1, 1.0],
    'KEQ': [0.25 * f for f in (1 / 4, 1, 4, 20)],
}

# Log of the study (console + SA/output.txt), shared by all worker threads
_LOG_FILE = None
_LOG_LOCK = threading.Lock()
_RUN_LOCK = threading.Lock()
_RUNNING = set()        # folders of the jobs currently running


def log(message):
    """Print a message and append it, with a timestamp, to SA/output.txt (thread-safe)."""
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}"
    with _LOG_LOCK:
        print(line, flush=True)
        if _LOG_FILE is not None:
            with open(_LOG_FILE, 'a') as f:
                f.write(line + '\n')


# Files copied into every run folder besides the job input (properties.inp is written, not copied)
COMMON_FILES = ['sec_uel_cube.inp', 'prefdir.inp', 'uel.f90', 'abaqus_v6.env', 'aba_param.inc']
FILES_TO_COPY = [f'{JOB_NAME}.inp'] + COMMON_FILES

# ------------------------------------------------------------------------------------------
# Properties
# ------------------------------------------------------------------------------------------
_OPS = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
        ast.Div: operator.truediv, ast.Pow: operator.pow, ast.USub: operator.neg}


def _safe_eval(expr):
    """Evaluate a numeric *parameter expression such as '25.0 + 273.0'."""
    def ev(node):
        if isinstance(node, ast.Expression):
            return ev(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return float(node.value)
        if isinstance(node, ast.BinOp) and type(node.op) in _OPS:
            return _OPS[type(node.op)](ev(node.left), ev(node.right))
        if isinstance(node, ast.UnaryOp) and type(node.op) in _OPS:
            return _OPS[type(node.op)](ev(node.operand))
        raise ValueError(f"Unsupported expression in properties.inp: {expr}")
    return ev(ast.parse(expr, mode='eval'))


def read_properties(path):
    """Read the active 'NAME = value' lines of an Abaqus *parameter file (INITMU excluded)."""
    props = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('*') or '=' not in line:
                continue
            key, value = (s.strip() for s in line.split('=', 1))
            if key != 'INITMU':
                props[key] = _safe_eval(value)
    return props


def compute_initmu(p, rgas=8.31446261815324, tol=1.e-12):
    """Initial chemical potential at chemical equilibrium (same equations as get_initmu.py)."""
    cr = p['CACTIN'] * p['R']
    fmax = p['CACTIN'] * p['RFMAX']
    cmax = p['CACTIN'] * p['RBMAX']

    def evalh(cb0):
        cf0 = cr - cb0
        thetaf0 = cf0 / fmax
        return p['KEQ'] * cf0 * (1.0 - cb0 / cmax) - cb0 * (1.0 - thetaf0) * math.exp(-p['CHI'] * (1.0 - 2.0 * thetaf0))

    a, b = 0.0, min(cr, cmax)
    while (b - a) / 2.0 > tol:
        m = 0.5 * (a + b)
        if evalh(m) == 0.0:
            a = b = m
            break
        if evalh(a) * evalh(m) < 0.0:
            b = m
        else:
            a = m
    cb0 = 0.5 * (a + b)
    thetaf0 = (cr - cb0) / fmax
    if not 0.0 < thetaf0 < 1.0:
        raise ValueError(f"Initial free fraction thetaf0 = {thetaf0:.4g} outside (0, 1): "
                         f"free crosslinkers exceed RFMAX*CACTIN for these properties.")
    jc = 1.0    # J^c = 1 + VMOL*(c - c0) is measured from the initial content c0 = cr (stress-free reference)
    return p['MU0'] + rgas * p['THETA'] * (
        math.log(thetaf0 / (1.0 - thetaf0)) + p['CHI'] * (1.0 - 2.0 * thetaf0)
        - (p['K'] * p['VMOL'] / (rgas * p['THETA'])) * (math.log(1.0 / jc) / jc))


def render_properties(props, initmu, header):
    lines = ['*parameter', f'** {header}']
    lines += [f'{key} = {value:.12g}' for key, value in props.items()]
    lines += ['*' * 79, f'INITMU = {initmu:.6f}']
    return '\n'.join(lines) + '\n'


# ------------------------------------------------------------------------------------------
# Cases and folders
# ------------------------------------------------------------------------------------------
def base_properties(overrides=None):
    """Base properties: ./properties.inp (INITMU excluded) with the given overrides applied."""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    props = read_properties(os.path.join(base_dir, 'properties.inp'))
    props.update({k: float(v) for k, v in (overrides or {}).items()})
    return props


def make_case(folder, label, base_props, overrides):
    """A case: run folder (relative to the study folder), label, and its properties, i.e. the
    base properties with the overrides applied (overrides may add parameters, e.g. GAMMA)."""
    props = dict(base_props)
    props.update({k: float(v) for k, v in overrides.items()})
    return dict(folder=folder, label=label, overrides=dict(overrides), props=props)


def build_cases(base_props, study, only=None):
    """One-at-a-time cases: <VAR>/<VAR>_<value>, plus base/ (base properties)."""
    cases = []
    if RUN_BASE_CASE and not only:
        cases.append(make_case('base', 'base properties', base_props, {}))
    for var, values in study.items():
        if only and var not in only:
            continue
        if var not in base_props:
            raise KeyError(f"{var} is not a property in properties.inp")
        if var == 'DXC' and base_props.get('KCATCH0', 0.0) == 0.0:
            log("WARNING: sweeping DXC with KCATCH0 = 0 has no effect on the results "
                "(no catch pathway). Set KCATCH0 in BASE_OVERRIDES.")
        for value in values:
            cases.append(make_case(os.path.join(var, f'{var}_{float(value):g}'),
                                   f'{var} = {float(value):g}', base_props, {var: value}))
    return cases


def set_case_inputs(case, sa_dir):
    """Run folder, INITMU and properties.inp content of a case (nothing written yet)."""
    case['run_dir'] = os.path.join(sa_dir, case['folder'])
    case['initmu'] = compute_initmu(case['props'])
    case['properties'] = render_properties(case['props'], case['initmu'], case['label'])


def inputs_unchanged(case, base_dir, files):
    """True if the run folder already holds exactly the inputs that would be written now."""
    run_dir = case['run_dir']
    path = os.path.join(run_dir, 'properties.inp')
    if not os.path.exists(path):
        return False
    with open(path) as f:
        if f.read() != case['properties']:
            return False
    for name in files:
        dest = os.path.join(run_dir, name)
        if not (os.path.exists(dest) and filecmp.cmp(os.path.join(base_dir, name), dest, shallow=False)):
            return False
    return True


def prepare_case(case, base_dir, files):
    """Write properties.inp and copy the job files into the run folder."""
    run_dir = case['run_dir']
    os.makedirs(run_dir, exist_ok=True)
    with open(os.path.join(run_dir, 'properties.inp'), 'w') as f:
        f.write(case['properties'])
    for name in files:
        shutil.copy(os.path.join(base_dir, name), os.path.join(run_dir, name))


def job_completed(run_dir, job_name=JOB_NAME):
    sta = os.path.join(run_dir, f'{job_name}.sta')
    if not os.path.exists(sta):
        return False
    with open(sta, errors='ignore') as f:
        return 'COMPLETED SUCCESSFULLY' in f.read()


def clean_previous_attempt(run_dir, job_name=JOB_NAME):
    """Remove the files of an unfinished attempt (the .lck file blocks a new run)."""
    for ext in ('.lck', '.odb', '.sta', '.msg', '.dat', '.com', '.prt', '.sim', '.log', '.exception'):
        path = os.path.join(run_dir, f'{job_name}{ext}')
        if os.path.exists(path):
            os.remove(path)


def run_abaqus(run_dir, job_name=JOB_NAME, cpus=CPUS_PER_JOB, poll_seconds=10):
    """Run one job in the foreground ('interactive'), with its output written to a file.

    The worker only returns once the job has finished, so that at most `workers` jobs run
    at the same time: 'interactive' keeps the launcher in the foreground, and as a safeguard
    (in case the launcher ever returns early, e.g. queued or background submission) the
    worker also waits while the job lock file <job>.lck exists (Abaqus removes it at the end).
    stdin is closed so that no Abaqus prompt can block a worker.
    """
    cmd = ['abaqus', f'job={job_name}', 'user=uel.f90', f'cpus={cpus}',
           'ask_delete=off', 'interactive']
    with open(os.path.join(run_dir, f'{job_name}_output.txt'), 'w') as out:
        result = subprocess.run(cmd, cwd=run_dir, stdin=subprocess.DEVNULL,
                                stdout=out, stderr=subprocess.STDOUT)
    lock = os.path.join(run_dir, f'{job_name}.lck')
    while os.path.exists(lock):
        time.sleep(poll_seconds)
    return result.returncode, job_completed(run_dir, job_name)


def run_case(case, sa_dir, job_name=JOB_NAME, cpus=CPUS_PER_JOB):
    """Run one case, logging its start and end and the jobs running at that moment."""
    folder = case['folder']
    # The set update and its log line are done together, so the logged lines are in order
    with _RUN_LOCK:
        _RUNNING.add(folder)
        log(f"START  {folder} ({case['label']}): {job_name}.inp with user=uel.f90 in "
            f"{os.path.join(os.path.basename(sa_dir), folder)} | running now ({len(_RUNNING)}): "
            f"{', '.join(sorted(_RUNNING))}")
    t_start = time.time()
    code, completed = run_abaqus(case['run_dir'], job_name, cpus)
    status = 'completed' if completed else f'FAILED (return code {code})'
    with _RUN_LOCK:
        _RUNNING.discard(folder)
        log(f"END    {folder}: {status} after {(time.time() - t_start) / 60:.1f} min | "
            f"still running ({len(_RUNNING)}): {', '.join(sorted(_RUNNING)) if _RUNNING else '-'}")
    return code, completed


# ------------------------------------------------------------------------------------------
# Study runner (shared by all study scripts)
# ------------------------------------------------------------------------------------------
def run_study(cases, job_name, sa_dir_name, files_to_copy=None, workers=MAX_WORKERS,
              cpus_per_job=CPUS_PER_JOB, dry_run=False, force=False, title=''):
    """Prepare the run folder of every case in <sa_dir_name>/, skip completed and unchanged
    cases, run the others in parallel (`workers` at a time), and write the log
    (<sa_dir_name>/output.txt) and the summary (<sa_dir_name>/sa_cases.csv)."""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    files = files_to_copy if files_to_copy is not None else [f'{job_name}.inp'] + COMMON_FILES
    missing = [f for f in files if not os.path.exists(os.path.join(base_dir, f))]
    if missing:
        raise FileNotFoundError(f"Missing files in {base_dir}: {missing}")
    sa_dir = os.path.join(base_dir, sa_dir_name)
    os.makedirs(sa_dir, exist_ok=True)

    global _LOG_FILE
    _LOG_FILE = os.path.join(sa_dir, 'output.txt')
    log('=' * 80)
    log(f"{title or 'Sensitivity analysis'}: job {job_name}.inp, {len(cases)} cases, options: "
        f"dry_run={dry_run}, workers={workers}, cpus_per_job={cpus_per_job}, force={force}")

    # A completed job is skipped (and its folder left untouched) only if its inputs are identical
    # to the ones that would be written now; otherwise the folder is refreshed and the job rerun.
    to_run = []
    for case in cases:
        set_case_inputs(case, sa_dir)
        completed = job_completed(case['run_dir'], job_name)
        if completed and not force and inputs_unchanged(case, base_dir, files):
            case['status'] = 'completed (skipped)'
            continue
        if completed and not force:
            log(f"{case['folder']}: inputs changed since the completed run -> will be rerun")
        prepare_case(case, base_dir, files)
        clean_previous_attempt(case['run_dir'], job_name)
        case['status'] = 'not run' if dry_run else 'pending'
        to_run.append(case)
    skipped = [c['folder'] for c in cases if c not in to_run]
    log(f"{len(cases)} cases in {sa_dir}: {len(skipped)} completed and unchanged (skipped), "
        f"{len(to_run)} to run")
    if skipped:
        log(f"Skipped: {', '.join(skipped)}")
    if to_run:
        log(f"{'Prepared (dry run)' if dry_run else 'To run'}: {', '.join(c['folder'] for c in to_run)}")

    if not dry_run and to_run:
        log(f"Running {len(to_run)} jobs, {workers} at a time...")
        t0 = time.time()
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(run_case, c, sa_dir, job_name, cpus_per_job): c for c in to_run}
            for n, future in enumerate(concurrent.futures.as_completed(futures), 1):
                case = futures[future]
                code, completed = future.result()
                case['status'] = 'completed' if completed else f'FAILED (return code {code})'
                log(f"Progress {n}/{len(to_run)} finished ({(time.time() - t0) / 60:.1f} min elapsed)")

    with open(os.path.join(sa_dir, 'sa_cases.csv'), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['folder', 'label', 'overrides', 'INITMU', 'status'])
        for c in cases:
            overrides = '; '.join(f'{k}={v:g}' for k, v in c['overrides'].items())
            writer.writerow([c['folder'], c['label'], overrides, f"{c['initmu']:.6f}", c['status']])

    failed = [c['folder'] for c in cases if c['status'].startswith('FAILED')]
    log(f"Summary written to {os.path.join(sa_dir, 'sa_cases.csv')}")
    if failed:
        log(f"{len(failed)} jobs failed: {', '.join(failed)}")
    return cases


# ------------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--dry-run', action='store_true', help='only create the run folders')
    parser.add_argument('--only', nargs='+', help='variables to include (default: all)')
    parser.add_argument('--workers', type=int, default=MAX_WORKERS, help='simultaneous Abaqus jobs')
    parser.add_argument('--force', action='store_true', help='rerun completed jobs')
    args = parser.parse_args()

    base_props = base_properties(BASE_OVERRIDES)
    cases = build_cases(base_props, STUDY_PROPS_INFO, args.only)
    run_study(cases, JOB_NAME, SA_DIR_NAME, FILES_TO_COPY, workers=args.workers,
              cpus_per_job=CPUS_PER_JOB, dry_run=args.dry_run, force=args.force,
              title=f"One-at-a-time sensitivity analysis (only={args.only or 'all'})")


if __name__ == '__main__':
    main()
