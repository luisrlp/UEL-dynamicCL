"""Stage 3 study: AFM indentation (spherical indenter, quarter model), catch-slip vs slip kinetics,
at two length scales, with filament prestress and deep indentation.

Inputs (one per scale and chemical boundary), generated with 1_build_mesh.py (stage-3 preset:
mesh 3 x 2.5 x 3 in units of the indenter radius R, graded towards the indented corner,
'base_mesh_afm') and 2_uel_afm.py (stage-3 preset), with length_scale = R in um:
  python3 2_uel_afm.py length_scale=20 output_file=stage3_afm_large_closed.inp
  python3 2_uel_afm.py length_scale=0.666666666667 output_file=stage3_afm_small_closed.inp
  (open faces, bath on the far side walls x = 0 and z = 0:
  python3 2_uel_afm.py length_scale=20 side_boundary_condition=open output_file=stage3_afm_large_open.inp
  python3 2_uel_afm.py length_scale=0.666666666667 side_boundary_condition=open output_file=stage3_afm_small_open.inp)
  large: R = 20 um, gel 60 x 50 x 60 um (quarter), reconstituted gel probed with a colloidal tip
  small: R = 0.67 um, sample 2 x 1.67 x 2 um (quarter), the same problem 30 times smaller
numElem in global.f90 must be >= the number of elements (3584), and uel.f90 regenerated
(compile_elem.sh) before running.

Protocol (2_uel_afm.py stage-3 preset): Equilibrate 200 s (indenter fixed: chemical equilibrium
under the filament prestress when LAMBDA0 > 1; nothing happens when LAMBDA0 = 1), Indentation
2 s (linear ramp to the depth), Hold 600 s. Chemistry closed on all faces (the symmetry planes
x = max, z = max are always no-flux).

Length scales: the mechanics has no intrinsic length (affine network, frictionless contact), so
the only length scale is that of transport, through Da = k l^2 / D_eff (l ~ contact radius
a = sqrt(R delta), k ~ koff ~ 0.05 1/s, D = 3 um^2/s): large scale, a ~ 4.5-11 um, diffusion
time a^2/D ~ 7-40 s, comparable to the kinetics (coupled reaction-diffusion); small scale, a ~
0.15-0.37 um, a^2/D ~ 0.01-0.05 s (diffusion instantaneous, local kinetics only). The force
scales as R^2 (equal stresses) when the transport regime is the same: the small sample with
D * (R_small/R_large)^2 must reproduce the large gel (D_scaling check).

Prestress: filament prestretch LAMBDA0 (extensible filament with compliant crosslinkers,
ETAC = 2/3, dx = dxc = 1 nm). Rest force f0 per filament, ratio to the catch optimum
f* ~ 4.5 pN, and small-strain shear modulus relative to LAMBDA0 = 1 (standalone driver):
  LAMBDA0 = 1.1: f0 = 0.07 pN, x9.5;  1.2: f0 = 0.46 pN, x86;  1.24: f0 = 1.4 pN (0.3 f*), x410
(the filament contour limit is reached at LAMBDA0 ~ 1.28, so prestressed cases are indented
shallowly; LAMBDA0 = 1.24 reaches the contour limit at shear strains ~ 0.3).

Base properties: ./properties.inp with the stage-1 design (run_stage1.STAGE1_BASE: ETAC = 2/3,
dx = dxc = 1 nm, KEQ = 0.25), D = 3 um^2/s and a weak neo-Hookean matrix (PHINET = 0.5,
C10 = 0.008: matrix shear modulus 2 (1 - PHINET) C10 = 0.008 Pa, ~5% of the network modulus)
that keeps the compressed region under the indenter stable (compressed filaments carry no load);
E_nomatrix checks its influence. Kinetic laws as in stages 1-2 (same koff(0) = 0.05 1/s).

Sweeps (folders SA_stage3/<scale>_<chem>/<sweep>/<case>/), delta = indentation depth:
  A_core       law x delta/R {0.05, 0.15, 0.3}                    large, small    (12)
  B_prestress  law x LAMBDA0 {1.1, 1.2, 1.24}, delta/R = 0.05     large, small    (12)
  B2_prestress_deep  law x LAMBDA0 {1.1, 1.2}, delta/R = 0.15     large            (4)
  C_D          law x D {0.3, 30} um^2/s, delta/R = 0.15            large            (4)
  D_scaling    law, D = 3 (R_small/R_large)^2, delta/R = 0.15      small            (2)
  E_nomatrix   law, pure network (PHINET = 1, C10 = 0), delta/R = 0.15   large      (2)
  F_bath       law x delta/R {0.15}, open far walls (--chem open)  large, small    (4)
(A_core with LAMBDA0 = 1 is the reference of B; A_core large, delta/R = 0.15, D = 3 that of C, D, E.)
Mesh and time-step convergence are checked separately (one input file per mesh).

Per case, properties.inp gets the indentation as Abaqus parameters (um): DELTA (depth),
RIND (indenter radius) and UIND = -(initial gap + DELTA) (indenter displacement, <UIND> in the
input); the radius and the gap are read from the '** INDENTER:' line of the input file.

Usage (from test_in_abaqus/):
  python3 run_stage3.py                              # all closed-face sweeps, both scales
  python3 run_stage3.py --scale large --only A_core B_prestress
  python3 run_stage3.py --chem closed open           # also the open-face cases (F_bath)
  python3 run_stage3.py --dry-run                    # only create the run folders
  python3 run_stage3.py --workers 4 --cpus 4 --force
"""
import argparse
import os
import re

import run_sa
from run_stage1 import STAGE1_BASE, LAWS, g

SA_DIR_NAME = 'SA_stage3'
MAX_WORKERS = 8
CPUS_PER_JOB = 4

D_BASE = 3.0                                    # um^2/s, as in stage 2
MATRIX = {'PHINET': 0.5, 'C10': 0.008}          # weak matrix, ~5% of the network shear modulus
STAGE3_BASE = {**STAGE1_BASE, 'D': D_BASE, 'LAMBDA0': 1.0, **MATRIX}

SCALES = ('large',) #, 'small')
CHEMS = ('closed', 'open')
DEPTHS = [0.05, 0.15, 0.3]                      # delta / R
DEPTH_REF = 0.15
LAMBDA0S = [1.1, 1.2, 1.24]


def job_name(scale, chem):
    return f'stage3_afm_{scale}_{chem}'


def read_indenter(inp_path):
    """Indenter radius and initial gap (um) from the '** INDENTER:' line written by 2_uel_afm.py."""
    with open(inp_path) as f:
        for line in f:
            if line.startswith('** INDENTER:'):
                vals = dict(re.findall(r'(\w+) = ([^,\s]+)', line))
                return float(vals['radius']), float(vals['initial_gap'])
    raise ValueError(f"No '** INDENTER:' line in {inp_path}: regenerate it with 2_uel_afm.py")


def sweep_cases(radius):
    """{sweep: [(scales, chem, case folder, label, overrides without the indentation, delta/R)]}."""
    r_small, r_large = radius.get('small'), radius.get('large')
    sweeps = {}
    sweeps['A_core'] = [
        (SCALES, 'closed', f'{law}_d{g(d)}', f'A: {law}, delta/R = {g(d)}', {**LAWS[law]}, d)
        for law in LAWS for d in DEPTHS]
    sweeps['B_prestress'] = [
        (SCALES, 'closed', f'{law}_lambda0{g(l0)}', f'B: {law}, LAMBDA0 = {g(l0)}, delta/R = {g(DEPTHS[0])}',
         {**LAWS[law], 'LAMBDA0': l0}, DEPTHS[0])
        for law in LAWS for l0 in LAMBDA0S]
    sweeps['B2_prestress_deep'] = [
        (('large', ), 'closed', f'{law}_lambda0{g(l0)}', f'B2: {law}, LAMBDA0 = {g(l0)}, delta/R = {g(DEPTH_REF)}',
         {**LAWS[law], 'LAMBDA0': l0}, DEPTH_REF)
        for law in LAWS for l0 in LAMBDA0S[:2]]
    sweeps['C_D'] = [
        (('large', ), 'closed', f'{law}_D{g(dd)}', f'C: {law}, D = {g(dd)} um^2/s, delta/R = {g(DEPTH_REF)}',
         {**LAWS[law], 'D': dd}, DEPTH_REF)
        for law in LAWS for dd in [0.3, 30.0]]
    if r_small and r_large:
        d_scaled = D_BASE * (r_small / r_large) ** 2
        sweeps['D_scaling'] = [
            (('small', ), 'closed', f'{law}_D{g(d_scaled)}',
             f'D: {law}, D = {g(d_scaled)} um^2/s (= D_large (R_small/R_large)^2), delta/R = {g(DEPTH_REF)}',
             {**LAWS[law], 'D': d_scaled}, DEPTH_REF)
            for law in LAWS]
    sweeps['E_nomatrix'] = [
        (('large', ), 'closed', f'{law}_nomatrix', f'E: {law}, pure network, delta/R = {g(DEPTH_REF)}',
         {**LAWS[law], 'PHINET': 1.0, 'C10': 0.0}, DEPTH_REF)
        for law in LAWS]
    sweeps['F_bath'] = [
        (SCALES, 'open', f'{law}_d{g(DEPTH_REF)}', f'F: {law}, open far walls, delta/R = {g(DEPTH_REF)}',
         {**LAWS[law]}, DEPTH_REF)
        for law in LAWS]
    return sweeps


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--scale', nargs='+', choices=SCALES, default=list(SCALES), help='length scales')
    parser.add_argument('--chem', nargs='+', choices=CHEMS, default=['closed'], help='chemical boundaries')
    parser.add_argument('--only', nargs='+', help='sweeps to run (default: all)')
    parser.add_argument('--dry-run', action='store_true', help='only create the run folders')
    parser.add_argument('--workers', type=int, default=MAX_WORKERS, help='simultaneous Abaqus jobs')
    parser.add_argument('--cpus', type=int, default=CPUS_PER_JOB, help='CPUs per Abaqus job')
    parser.add_argument('--force', action='store_true', help='rerun completed jobs')
    args = parser.parse_args()

    base_dir = os.path.dirname(os.path.abspath(__file__))
    jobs = [(s, c) for s in args.scale for c in args.chem]
    indenter = {}
    for scale, chem in jobs:
        inp = os.path.join(base_dir, job_name(scale, chem) + '.inp')
        if not os.path.exists(inp):
            raise FileNotFoundError(f"{inp} not found: generate it with 2_uel_afm.py (see the usage above)")
        indenter[scale, chem] = read_indenter(inp)
    radius = {s: r for (s, c), (r, gap) in indenter.items()}
    for s in SCALES:     # the scaling sweep needs both radii: read the other scale too if present
        inp = os.path.join(base_dir, job_name(s, 'closed') + '.inp')
        if s not in radius and os.path.exists(inp):
            radius[s] = read_indenter(inp)[0]

    sweeps = sweep_cases(radius)
    if args.only:
        unknown = [s for s in args.only if s not in sweeps]
        if unknown:
            raise SystemExit(f"Unknown sweeps {unknown}; available: {list(sweeps)}")

    base_props = run_sa.base_properties(STAGE3_BASE)
    for scale, chem in jobs:
        r_ind, gap = indenter[scale, chem]
        cases = []
        for sweep, items in sweeps.items():
            if args.only and sweep not in args.only:
                continue
            for scales, case_chem, folder, label, overrides, d_r in items:
                if scale not in scales or case_chem != chem:
                    continue
                delta = d_r * r_ind
                indentation = {'DELTA': delta, 'RIND': r_ind, 'UIND': -(gap + delta)}
                cases.append(run_sa.make_case(os.path.join(sweep, folder), f'{label} [{scale}, R = {g(r_ind)} um]',
                                              base_props, {**overrides, **indentation}))
        if not cases:
            continue
        run_sa.run_study(cases, job_name(scale, chem), os.path.join(SA_DIR_NAME, f'{scale}_{chem}'),
                         workers=args.workers, cpus_per_job=args.cpus, dry_run=args.dry_run, force=args.force,
                         title=f"Stage 3 (AFM indentation, {scale} scale R = {g(r_ind)} um, {chem} faces), "
                               f"sweeps: {', '.join(args.only or sweeps)}")


if __name__ == '__main__':
    main()
