"""Stage 1 study: single-element simple shear (step + hold), catch-slip vs slip kinetics.

Input: stage1_shear_hold.inp (generated with 1_build_mesh.py + 3_uel_cube.py, stage-1 settings):
one element, all displacements prescribed (homogeneous simple shear, top u_x = <GAMMA>, ramped
over 0.5 s, then held for 200 s) and mu = INITMU at all nodes (open system: free crosslinker
concentration held at the bath value), so only the kinetics change the response.

Base properties: ./properties.inp with the stage-1 overrides in STAGE1_BASE (pure network, no
matrix; ETAC = 2/3; physical reactive distances of 1 nm). Every case sets the kinetic law
explicitly (KOFF0, KCATCH0) and the shear amplitude GAMMA. Laws with the same unloaded off-rate
koff(0) = KOFF0 + KCATCH0 (hence same KEQ and unloaded kinetics) are compared, so that only the
force dependence differs.

Force window: with dx = dxc = 1 nm and a catch fraction of 0.9, the catch optimum is
f* = kT/(dx + dxc) ln 9 ~ 4.5 pN. In this network, filament forces of 1-10 pN only appear in the
most stretched directions, in the nonlinear (stiffening) regime: GAMMA ~ 0.45-0.55 for
ETAC = 2/3 (~0.3-0.4 for ETAC = 1, ~0.8-1.0 for ETAC = 1/3). The amplitudes below span that
window, plus linear-regime controls (GAMMA = 0.1, 0.3) where both laws should coincide.

Sweeps (folders SA_stage1/<sweep>/<case>/):
  A_core          law {slip, catch} x GAMMA {0.1, 0.3, 0.4, 0.45, 0.5, 0.55}          (12)
  B_catch_frac    catch fraction KCATCH0/koff(0) {0, 0.5, 0.8, 0.9, 0.95} x GAMMA     (20)
  C_dx            DX = DXC {0.5, 1, 2} nm x law x GAMMA                              (24)
  C2_dxc_ratio    DXC/DX {0.5, 1, 2} (DX = 1 nm, catch law) x GAMMA                   (12)
  D_etac          ETAC {1/3, 2/3, 1} x law x GAMMA window of each ETAC               (18)
  E_rate          koff(0) {0.01, 0.05, 0.25} x law, GAMMA = 0.5                       (6)
  (GAMMA for B, C, C2: {0.3, 0.45, 0.5, 0.55})
Some cases coincide across sweeps (e.g. A catch GAMMA=0.5 and C_dx catch DX=1 nm GAMMA=0.5);
they are run in each sweep folder (single-element runs are cheap).

Usage (from test_in_abaqus/):
  python3 run_stage1.py                       # all sweeps
  python3 run_stage1.py --only A_core D_etac  # some sweeps
  python3 run_stage1.py --dry-run             # only create the run folders
  python3 run_stage1.py --workers 8 --force
"""
import argparse
import os

import run_sa

JOB_NAME = 'stage1_shear_hold'
SA_DIR_NAME = 'SA_stage1'
MAX_WORKERS = 8
CPUS_PER_JOB = 1       # one element: parallel threads only add overhead

# Stage-1 base: changes to ./properties.inp applied to every case. The values that define the
# stage-1 design are fixed here, so the study does not depend on the current state of
# properties.inp (the remaining properties are taken from it)
DX_BASE = 0.001        # 1 nm (in um): physical Bell reactive distance
STAGE1_BASE = {
    'PHINET': 1.0, 'C10': 0.0, 'C01': 0.0,   # pure network (no neo-Hookean matrix)
    'ETAC': 2.0 / 3.0,                       # GAMMA windows below chosen for this ETAC
    'DX': DX_BASE, 'DXC': DX_BASE,           # catch optimum f* ~ 4.5 pN with catch fraction 0.9
    'KEQ': 0.25,
}

# Kinetic laws, same unloaded off-rate koff(0) = KOFF0 + KCATCH0 = 0.05 1/s
KOFF_TOTAL = 0.05
LAWS = {
    'slip': {'KOFF0': KOFF_TOTAL, 'KCATCH0': 0.0},
    'catch': {'KOFF0': 0.1 * KOFF_TOTAL, 'KCATCH0': 0.9 * KOFF_TOTAL},   # f* = 4.5 pN, lifetime x1.67
}
GAMMAS_CORE = [0.1, 0.3, 0.4, 0.45, 0.5, 0.55]   # linear controls + force window (ETAC = 2/3)
GAMMAS = [0.3, 0.45, 0.5, 0.55]                   # one linear control + force window (ETAC = 2/3)
GAMMA_WINDOW_ETAC = {1.0 / 3.0: [0.8, 0.9, 1.0],  # amplitudes reaching 1-10 pN for each ETAC
                     2.0 / 3.0: [0.45, 0.5, 0.55],
                     1.0: [0.3, 0.35, 0.4]}
GAMMA_RATE = 0.5


def g(x):
    """Compact number format for folder names and labels."""
    return f'{float(x):g}'


def nm(x_um):
    """Reactive distance in nm (from um) for folder names and labels."""
    return f'{1000.0 * x_um:g}'


def sweep_cases():
    """{sweep name: [(case folder, label, overrides)]} for all stage-1 sweeps."""
    sweeps = {}

    sweeps['A_core'] = [
        (f'{law}_gamma{g(gm)}', f'A: {law}, GAMMA = {g(gm)}', {**LAWS[law], 'GAMMA': gm})
        for law in LAWS for gm in GAMMAS_CORE]

    sweeps['B_catch_frac'] = [
        (f'cf{g(cf)}_gamma{g(gm)}', f'B: catch fraction = {g(cf)}, GAMMA = {g(gm)}',
         {'KOFF0': (1.0 - cf) * KOFF_TOTAL, 'KCATCH0': cf * KOFF_TOTAL, 'GAMMA': gm})
        for cf in [0.0, 0.5, 0.8, 0.9, 0.95] for gm in GAMMAS]

    sweeps['C_dx'] = [
        (f'{law}_dx{nm(dx)}nm_gamma{g(gm)}', f'C: {law}, DX = DXC = {nm(dx)} nm, GAMMA = {g(gm)}',
         {**LAWS[law], 'DX': dx, 'DXC': dx, 'GAMMA': gm})
        for law in LAWS for dx in [0.0005, 0.001, 0.002] for gm in GAMMAS]

    sweeps['C2_dxc_ratio'] = [
        (f'ratio{g(r)}_gamma{g(gm)}', f'C2: catch, DXC/DX = {g(r)} (DX = {nm(DX_BASE)} nm), GAMMA = {g(gm)}',
         {**LAWS['catch'], 'DX': DX_BASE, 'DXC': DX_BASE * r, 'GAMMA': gm})
        for r in [0.5, 1.0, 2.0] for gm in GAMMAS]

    sweeps['D_etac'] = [
        (f'{law}_etac{eta:.4g}_gamma{g(gm)}', f'D: {law}, ETAC = {eta:.4g}, GAMMA = {g(gm)}',
         {**LAWS[law], 'ETAC': eta, 'GAMMA': gm})
        for law in LAWS for eta, gammas in GAMMA_WINDOW_ETAC.items() for gm in gammas]

    sweeps['E_rate'] = [
        (f'{law}_koff{g(k)}', f'E: {law}, koff(0) = {g(k)}, GAMMA = {g(GAMMA_RATE)}',
         {'KOFF0': LAWS[law]['KOFF0'] / KOFF_TOTAL * k, 'KCATCH0': LAWS[law]['KCATCH0'] / KOFF_TOTAL * k,
          'GAMMA': GAMMA_RATE})
        for law in LAWS for k in [0.01, 0.05, 0.25]]

    return sweeps


def main():
    sweeps = sweep_cases()
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--only', nargs='+', choices=list(sweeps), help='sweeps to run (default: all)')
    parser.add_argument('--dry-run', action='store_true', help='only create the run folders')
    parser.add_argument('--workers', type=int, default=MAX_WORKERS, help='simultaneous Abaqus jobs')
    parser.add_argument('--force', action='store_true', help='rerun completed jobs')
    args = parser.parse_args()

    base_props = run_sa.base_properties(STAGE1_BASE)
    cases = [run_sa.make_case(os.path.join(sweep, folder), label, base_props, overrides)
             for sweep, items in sweeps.items() if not args.only or sweep in args.only
             for folder, label, overrides in items]

    run_sa.run_study(cases, JOB_NAME, SA_DIR_NAME, workers=args.workers, cpus_per_job=CPUS_PER_JOB,
                     dry_run=args.dry_run, force=args.force,
                     title=f"Stage 1 (single-element shear-hold), sweeps: {', '.join(args.only or sweeps)}")


if __name__ == '__main__':
    main()
