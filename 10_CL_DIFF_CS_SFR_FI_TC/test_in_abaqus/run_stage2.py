"""Stage 2 study: column along x with a bath on one end face, simple shear (step + hold).

Input: stage2_column_bath.inp (generated with 1_build_mesh.py + 3_uel_cube.py, stage-2 presets):
column 50 x 1 x 1 um along x, 15 elements graded towards the bath face x = 50 (1.8 to 5.5 um),
homogeneous simple shear prescribed on all nodes (u_x = GAMMA*y over the height 1, ramped over
0.5 s, then held for 1000 s), mu = INITMU on the bath face x = 50 only (all other faces no-flux).
The mechanics is the same as in stage 1; the crosslinkers now have to diffuse to/from the bath.

Base properties: ./properties.inp with the stage-1 overrides (run_stage1.STAGE1_BASE: pure
network, ETAC = 2/3, dx = dxc = 1 nm, KEQ = 0.25) and a realistic diffusivity D = 3 um^2/s
(reaction-diffusion length sqrt(D/koff) ~ 8 um; diffusion distance in the hold ~ sqrt(D_eff t) ~ 50 um).
Kinetic laws as in stage 1 (same unloaded off-rate koff(0) = 0.05 1/s).

Sweeps (folders SA_stage2/<sweep>/<case>/):
  A_core     law {slip, catch} x GAMMA {0.1, 0.45, 0.5, 0.55}                   (8)
  B_D        D {1, 3, 10, 30} um^2/s x law, GAMMA = 0.5                         (8)
  C_KEQ      KEQ {0.0625, 0.25, 1, 5} x law, GAMMA = 0.5                         (8)
  D_RBMAX    RBMAX {0.25, 0.5, 1} x law, GAMMA = 0.5                             (6)
Reference limits (no extra runs): the closed column is identical to the closed single element
(stage 1 without bath); the open single element (stage 1) is the limit of fast diffusion.
Mesh convergence (N = 10, 15, 25) is checked separately, with one input file per mesh.

Usage (from test_in_abaqus/):
  python3 run_stage2.py                       # all sweeps
  python3 run_stage2.py --only A_core B_D     # some sweeps
  python3 run_stage2.py --dry-run             # only create the run folders
  python3 run_stage2.py --workers 6 --force
"""
import argparse
import os

import run_sa
from run_stage1 import STAGE1_BASE, LAWS, g

JOB_NAME = 'stage2_column_bath'
SA_DIR_NAME = 'SA_stage2'
MAX_WORKERS = 15
CPUS_PER_JOB = 1       # 15 elements: one CPU per job, several jobs in parallel

D_BASE = 3.0           # um^2/s, free crosslinker diffusivity in the network
STAGE2_BASE = {**STAGE1_BASE, 'D': D_BASE}

GAMMAS_CORE = [0.1, 0.45, 0.5, 0.55]   # linear control + force window (ETAC = 2/3)
GAMMA_SWEEPS = 0.5                     # amplitude for the B, C, D sweeps


def sweep_cases():
    """{sweep name: [(case folder, label, overrides)]} for all stage-2 sweeps."""
    sweeps = {}

    sweeps['A_core'] = [
        (f'{law}_gamma{g(gm)}', f'A: {law}, GAMMA = {g(gm)}', {**LAWS[law], 'GAMMA': gm})
        for law in LAWS for gm in GAMMAS_CORE]

    sweeps['B_D'] = [
        (f'{law}_D{g(d)}', f'B: {law}, D = {g(d)} um^2/s, GAMMA = {g(GAMMA_SWEEPS)}',
         {**LAWS[law], 'D': d, 'GAMMA': GAMMA_SWEEPS})
        for law in LAWS for d in [1.0, 3.0, 10.0, 30.0]]

    sweeps['C_KEQ'] = [
        (f'{law}_keq{g(k)}', f'C: {law}, KEQ = {g(k)}, GAMMA = {g(GAMMA_SWEEPS)}',
         {**LAWS[law], 'KEQ': k, 'GAMMA': GAMMA_SWEEPS})
        for law in LAWS for k in [0.0625, 0.25, 1.0, 5.0]]

    # sweeps['D_RBMAX'] = [
    #     (f'{law}_rbmax{g(rb)}', f'D: {law}, RBMAX = {g(rb)}, GAMMA = {g(GAMMA_SWEEPS)}',
    #      {**LAWS[law], 'RBMAX': rb, 'GAMMA': GAMMA_SWEEPS})
    #     for law in LAWS for rb in [0.25, 0.5, 1.0]]
    sweeps['D_RFMAX'] = [
        (f'{law}_rfmax{g(rb)}', f'D: {law}, RFMAX = {g(rb)}, GAMMA = {g(GAMMA_SWEEPS)}',
         {**LAWS[law], 'RFMAX': rb, 'GAMMA': GAMMA_SWEEPS})
        for law in LAWS for rb in [0.25, 0.5, 1.0, 2.]]
    

    return sweeps


def main():
    sweeps = sweep_cases()
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--only', nargs='+', choices=list(sweeps), help='sweeps to run (default: all)')
    parser.add_argument('--dry-run', action='store_true', help='only create the run folders')
    parser.add_argument('--workers', type=int, default=MAX_WORKERS, help='simultaneous Abaqus jobs')
    parser.add_argument('--force', action='store_true', help='rerun completed jobs')
    args = parser.parse_args()

    base_props = run_sa.base_properties(STAGE2_BASE)
    cases = [run_sa.make_case(os.path.join(sweep, folder), label, base_props, overrides)
             for sweep, items in sweeps.items() if not args.only or sweep in args.only
             for folder, label, overrides in items]

    run_sa.run_study(cases, JOB_NAME, SA_DIR_NAME, workers=args.workers, cpus_per_job=CPUS_PER_JOB,
                     dry_run=args.dry_run, force=args.force,
                     title=f"Stage 2 (column with bath, shear-hold), sweeps: {', '.join(args.only or sweeps)}")


if __name__ == '__main__':
    main()
