"""One-at-a-time sensitivity analysis of the UEL model in Abaqus.

For every variable in STUDY_PROPS_INFO and every value swept for it, a run folder
SA/<VAR>/<VAR>_<value>/ is created with:
  - properties.inp : base properties with only <VAR> changed, and INITMU recomputed for
                     that set of properties (initial chemical equilibrium, as get_initmu.py)
  - the files needed to run the job: JOB_NAME.inp, sec_uel_cube.inp, prefdir.inp (read by
    uexternaldb at start-up), uel.f90, abaqus_v6.env, aba_param.inc
All jobs are then run in parallel (MAX_WORKERS at a time), each in its own folder.

Base properties are read from ./properties.inp (single source of truth), with optional
overrides in BASE_OVERRIDES. A summary of all cases is written to SA/sa_cases.csv.

Usage (from test_in_abaqus/):
  python3 run_sa.py                  # create all folders and run all jobs
  python3 run_sa.py --dry-run        # only create the folders (no Abaqus)
  python3 run_sa.py --only DX KEQ    # restrict to some variables
  python3 run_sa.py --workers 8      # number of simultaneous Abaqus jobs
  python3 run_sa.py --force          # rerun jobs that already completed
Completed jobs (.sta containing "COMPLETED") are skipped, so the script can be relaunched
to finish or retry a study.
"""
import argparse
import ast
import concurrent.futures
import csv
import math
import operator
import os
import shutil
import subprocess
import time

# ------------------------------------------------------------------------------------------
# Study definition
# ------------------------------------------------------------------------------------------
JOB_NAME = 'cube_indent_uel_relax_auto'
SA_DIR_NAME = 'SA'
MAX_WORKERS = 4        # simultaneous Abaqus jobs (limited by CPUs and license tokens)
CPUS_PER_JOB = 1
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

# Files copied into every run folder (properties.inp is written, not copied)
FILES_TO_COPY = [f'{JOB_NAME}.inp', 'sec_uel_cube.inp', 'prefdir.inp', 'uel.f90',
                 'abaqus_v6.env', 'aba_param.inc']

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
    jc = 1.0 + p['VMOL'] * cr
    return p['MU0'] + rgas * p['THETA'] * (
        math.log(thetaf0 / (1.0 - thetaf0)) + p['CHI'] * (1.0 - 2.0 * thetaf0)
        - (p['K'] * p['VMOL'] / (rgas * p['THETA'])) * (math.log(1.0 / jc) / jc))


def write_properties(path, props, initmu, header):
    with open(path, 'w') as f:
        f.write('*parameter\n')
        f.write(f'** {header}\n')
        for key, value in props.items():
            f.write(f'{key} = {value:.12g}\n')
        f.write('*' * 79 + '\n')
        f.write(f'INITMU = {initmu:.6f}\n')


# ------------------------------------------------------------------------------------------
# Cases and folders
# ------------------------------------------------------------------------------------------
def build_cases(base_props, study, only=None):
    cases = []
    if RUN_BASE_CASE and not only:
        cases.append(dict(var='base', value='', folder=os.path.join('base'), props=dict(base_props)))
    for var, values in study.items():
        if only and var not in only:
            continue
        if var not in base_props:
            raise KeyError(f"{var} is not a property in properties.inp")
        if var == 'DXC' and base_props.get('KCATCH0', 0.0) == 0.0:
            print("WARNING: sweeping DXC with KCATCH0 = 0 has no effect on the results "
                  "(no catch pathway). Set KCATCH0 in BASE_OVERRIDES.")
        for value in values:
            props = dict(base_props)
            props[var] = float(value)
            cases.append(dict(var=var, value=float(value),
                              folder=os.path.join(var, f'{var}_{float(value):g}'), props=props))
    return cases


def prepare_case(case, base_dir, sa_dir):
    run_dir = os.path.join(sa_dir, case['folder'])
    os.makedirs(run_dir, exist_ok=True)
    case['run_dir'] = run_dir
    case['initmu'] = compute_initmu(case['props'])
    header = 'base properties' if case['var'] == 'base' else f"{case['var']} = {case['value']:g}"
    write_properties(os.path.join(run_dir, 'properties.inp'), case['props'], case['initmu'], header)
    for name in FILES_TO_COPY:
        shutil.copy(os.path.join(base_dir, name), os.path.join(run_dir, name))


def job_completed(run_dir):
    sta = os.path.join(run_dir, f'{JOB_NAME}.sta')
    if not os.path.exists(sta):
        return False
    with open(sta, errors='ignore') as f:
        return 'COMPLETED SUCCESSFULLY' in f.read()


def clean_previous_attempt(run_dir):
    """Remove the files of an unfinished attempt (the .lck file blocks a new run)."""
    for ext in ('.lck', '.odb', '.sta', '.msg', '.dat', '.com', '.prt', '.sim', '.log', '.exception'):
        path = os.path.join(run_dir, f'{JOB_NAME}{ext}')
        if os.path.exists(path):
            os.remove(path)


def run_abaqus(run_dir, poll_seconds=10):
    """Run one job in the foreground ('interactive'), with its output written to a file.

    The worker only returns once the job has finished, so that at most MAX_WORKERS jobs run
    at the same time: 'interactive' keeps the launcher in the foreground, and as a safeguard
    (in case the launcher ever returns early, e.g. queued or background submission) the
    worker also waits while the job lock file <job>.lck exists (Abaqus removes it at the end).
    stdin is closed so that no Abaqus prompt can block a worker.
    """
    cmd = ['abaqus', f'job={JOB_NAME}', 'user=uel.f90', f'cpus={CPUS_PER_JOB}',
           'ask_delete=off', 'interactive']
    with open(os.path.join(run_dir, f'{JOB_NAME}_output.txt'), 'w') as out:
        result = subprocess.run(cmd, cwd=run_dir, stdin=subprocess.DEVNULL,
                                stdout=out, stderr=subprocess.STDOUT)
    lock = os.path.join(run_dir, f'{JOB_NAME}.lck')
    while os.path.exists(lock):
        time.sleep(poll_seconds)
    return result.returncode, job_completed(run_dir)


# ------------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--dry-run', action='store_true', help='only create the run folders')
    parser.add_argument('--only', nargs='+', help='variables to include (default: all)')
    parser.add_argument('--workers', type=int, default=MAX_WORKERS, help='simultaneous Abaqus jobs')
    parser.add_argument('--force', action='store_true', help='rerun completed jobs')
    args = parser.parse_args()

    base_dir = os.path.dirname(os.path.abspath(__file__))
    sa_dir = os.path.join(base_dir, SA_DIR_NAME)
    os.makedirs(sa_dir, exist_ok=True)

    base_props = read_properties(os.path.join(base_dir, 'properties.inp'))
    base_props.update(BASE_OVERRIDES)
    cases = build_cases(base_props, STUDY_PROPS_INFO, args.only)

    for case in cases:
        prepare_case(case, base_dir, sa_dir)
    print(f"Prepared {len(cases)} run folders in {sa_dir}")

    to_run = []
    for case in cases:
        if job_completed(case['run_dir']) and not args.force:
            case['status'] = 'completed (skipped)'
        else:
            clean_previous_attempt(case['run_dir'])
            case['status'] = 'not run' if args.dry_run else 'pending'
            to_run.append(case)

    if not args.dry_run and to_run:
        print(f"Running {len(to_run)} jobs, {args.workers} at a time...")
        t0 = time.time()
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(run_abaqus, c['run_dir']): c for c in to_run}
            for n, future in enumerate(concurrent.futures.as_completed(futures), 1):
                case = futures[future]
                code, completed = future.result()
                case['status'] = 'completed' if completed else f'FAILED (return code {code})'
                print(f"[{n}/{len(to_run)}] {case['folder']}: {case['status']} "
                      f"({(time.time() - t0) / 60:.1f} min elapsed)")

    with open(os.path.join(sa_dir, 'sa_cases.csv'), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['variable', 'value', 'folder', 'INITMU', 'status'])
        for c in cases:
            writer.writerow([c['var'], c['value'], c['folder'], f"{c['initmu']:.6f}", c['status']])

    failed = [c['folder'] for c in cases if c['status'].startswith('FAILED')]
    print(f"Summary written to {os.path.join(sa_dir, 'sa_cases.csv')}")
    if failed:
        print(f"{len(failed)} jobs failed: {failed}")


if __name__ == '__main__':
    main()
