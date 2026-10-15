"""Verification checks for the catch-slip crosslinker kinetics (CATCH_BOND_PLAN.md, step 6).

Checks
  1  Regression : k_catch_0 = 0 reproduces the original slip-only Bell model
                  (unclipped Bell law, plain Newton kinetics, no sub-stepping).
  2  Derivatives: dKoff/df (_koff) and dH/dcb (_kinetics_residual) vs central finite differences.
  3  Single direction: koff(f) and lifetime, steady-state cb(lambda), dH/dcb at the root vs dt.
  4  Network response: step and cyclic shear, slip-only vs catch-slip, dt-convergence and sub-steps.
  5  Prestress trend: relaxation time of a small shear step superposed on a held prestrain
                  (catch-slip: relaxation should slow down with prestrain; slip-only: speed up).

Usage
  python catch_bond_checks.py            # all checks
  python catch_bond_checks.py 1 2        # selected checks
Figures are written to figures/catch_bond_checks/.

The catch-slip parameters below are illustrative (not calibrated): the catch pathway is a split of
the slip-only k_off_0, so both parameter sets have the same unloaded off-rate and KEQ.
"""
import io
import os
import sys
import time
import contextlib
import numpy as np
import matplotlib.pyplot as plt

from material import DynamicCrosslinking
from kinematics import Kinematics
from main import solve_shear, solve_uniaxial, profile_step, profile_cyclic

plt.rcParams['text.usetex'] = False          # plain labels for the diagnostic figures
FIGDIR = os.path.join('figures', 'catch_bond_checks')

THETA = 25 + 273.0
KT = 1.38065e-5 * THETA                      # pN.um

# Kinetic parameter sets: slip-only (main.py defaults) and catch-slip (split of the same k_off_0)
SLIP = dict(ks0=0.05, kc0=0.0, dxs=0.1, dxc=0.0)
CATCH = dict(ks0=0.01, kc0=0.04, dxs=0.1, dxc=0.1)


def build(kin, cls=DynamicCrosslinking, factor=3, **sim):
    """Material with the main.py parameters and the kinetic set `kin` (dx in um)."""
    diffusion_props = {
        'k_off_0': kin['ks0'], 'k_on_0': 0.0125, 'k_catch_0': kin['kc0'],
        'CHI': 0.1, 'THETA': THETA,
        'dx_kT': kin['dxs'] / KT, 'dx_c_kT': kin['dxc'] / KT,
        'D': 0.1, 'MU0': 0.0, 'VMOL': 0.15, 'RGAS': 8.31446261815324,
        'CACTIN': 0.0095, 'R': 0.1, 'RBMAX': 0.25, 'RFMAX': 1.,
    }
    mech_props = {
        'KBULK': 1.e3, 'LAMBDA0': 1.0, 'MU0STRETCH': 1.e8, 'Lp': 16.0, 'BETA': 0.5, 'R0C': 0.014,
        'ETA': 2. / 3., 'a_ratio': 1.2, 'b_disp': 0.000001, 'prefdir': np.array([1.0, 0.0, 0.0]),
    }
    sim_props = {'deformation_type': 'shear', 'factor': factor, 'max_iter': 50, 'tol': 1.0e-12}
    sim_props.update(sim)
    with contextlib.redirect_stdout(io.StringIO()):
        return cls(diffusion_props, mech_props, sim_props)


def quiet():
    return contextlib.redirect_stdout(io.StringIO())


def cb_total(mat):
    return float(np.sum(mat.cb_history * mat.weights))


def savefig(fig, name):
    os.makedirs(FIGDIR, exist_ok=True)
    path = os.path.join(FIGDIR, name)
    fig.tight_layout()
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    print(f"   figure: {path}")


# ---------------------------------------------------------------------------------------------
# Check 1: regression against the original slip-only implementation
# ---------------------------------------------------------------------------------------------
class OriginalBell(DynamicCrosslinking):
    """Original model: unclipped slip-only Bell law and plain (clipped) Newton kinetics."""

    def _koff(self, f):
        k = self.diff_props['k_off_0'] * np.exp(self.diff_props['dx_kT'] * f)
        return k, self.diff_props['dx_kT'] * k

    def _solve_local_kinetics(self, kinematics, dir, cb_old, Kon, B, dt, warn=True, lamn=1.0):
        li = np.linalg.norm(kinematics.F_iso @ dir)
        compression = li <= 1.0
        cbtau, converged = cb_old, False
        for _ in range(self.sim_props['max_iter']):
            r0f = 1.6 * (cbtau * 1.e3) ** (-2 / 5)
            if compression:
                f, DfDcb = 0.0, 0.0
            else:
                lf = self._lambdaf(r0f, li)
                f = self._f_inext(r0f, lf)
                DfDcb = self.dfdcb(f, cbtau, r0f, li, lf)
            Koff = self.diff_props['k_off_0'] * np.exp(self.diff_props['dx_kT'] * f)
            thetab = cbtau / self.CB_MAX
            Ri = B - Koff * self.CB_MAX * thetab / (1 - thetab)
            res = cbtau - cb_old - Ri * dt
            if abs(res) < self.sim_props['tol']:
                converged = True
                break
            DRiDcb = - Koff * (cbtau / (1 - thetab) * self.diff_props['dx_kT'] * DfDcb + 1 / (1 - thetab) ** 2)
            cbtau = float(np.clip(cbtau - res / (1 - DRiDcb * dt), 1e-10, self.CB_MAX - 1e-10))
        dH = self._kinetics_residual(cbtau, li, compression, cb_old, B, dt)[1]
        return cbtau, li, Kon, Koff, Ri, converged, dH, f


def run_profile(mat, mode, profile, dt, T):
    """Strain/stretch-controlled run; returns time, sigma (3x3 list), cb_total, nsub per step."""
    out = {'t': [], 'sigma': [], 'cb_tot': [], 'cb_max': [], 'cb_min': [], 'nsub': []}
    mu = mat.INITMU
    t_array = np.arange(0.0, T + 0.5 * dt, dt)
    with quiet():
        for i, t in enumerate(t_array):
            step_dt = 0.0 if i == 0 else dt
            val = profile(t)
            if mode == 'shear':
                sigma, *_ = solve_shear(mat, val, mu, step_dt)
            else:
                sigma, *_ = solve_uniaxial(mat, val, mu, step_dt)
            out['t'].append(t)
            out['sigma'].append(np.array(sigma))
            out['cb_tot'].append(cb_total(mat))
            out['cb_max'].append(float(np.max(mat.cb_history)))
            out['cb_min'].append(float(np.min(mat.cb_history)))
            out['nsub'].append(mat.nsub)
    return {k: np.array(v) for k, v in out.items()}


def check1():
    print("\n=== Check 1: regression, k_catch_0 = 0 vs original slip-only implementation")
    print("   (PASS/FAIL with sub-stepping off in both models; the effect of the default sub-stepping is reported separately)")
    ok_all = True
    for mode, ramp in [('shear', lambda t: 0.2 * min(t / 2.0, 1.0)),
                       ('uniaxial', lambda t: 1.0 + 0.2 * min(t / 2.0, 1.0))]:
        ref = run_profile(build(SLIP, cls=OriginalBell, max_substeps=1), mode, ramp, 0.05, 4.0)
        scale = np.max(np.abs(ref['sigma']))
        for label, sim in [('sub-stepping off', {'max_substeps': 1}), ('default sub-stepping', {})]:
            new = run_profile(build(SLIP, **sim), mode, ramp, 0.05, 4.0)
            dsig = np.max(np.abs(new['sigma'] - ref['sigma'])) / scale
            dcb = np.max(np.abs(new['cb_tot'] - ref['cb_tot']) / ref['cb_tot'])
            line = (f"   {mode:8s} {label:20s}: max rel diff stress {dsig:.1e}, cb_tot {dcb:.1e}, "
                    f"max sub-steps {new['nsub'].max()}")
            if label == 'sub-stepping off':
                ok = dsig < 1e-6 and dcb < 1e-6
                ok_all &= ok
                print(line + f"  -> {'PASS' if ok else 'FAIL'}")
            else:
                print(line + "  (differences = accuracy gained by sub-stepping)")
    return ok_all


# ---------------------------------------------------------------------------------------------
# Check 2: derivatives
# ---------------------------------------------------------------------------------------------
def check2():
    print("\n=== Check 2: derivatives vs central finite differences (catch-slip parameters)")
    mat = build(CATCH)
    worst = 0.0
    for f in [0.0, 0.01, 0.03, 0.1, 0.3]:
        h = 1e-6 * max(f, 1e-3)
        k, dk = mat._koff(f)
        fd = (mat._koff(f + h)[0] - mat._koff(f - h)[0]) / (2 * h)
        err = abs(dk - fd) / max(abs(fd), 1e-30)
        worst = max(worst, err)
        print(f"   dKoff/df  f={f:5.2f} pN: analytic {dk:+.6e}, FD {fd:+.6e}, rel.err {err:.1e}")

    thetaf = mat.cf0 / mat.CF_MAX
    B = mat._binding_rate(thetaf)[1]
    for li in [1.05, 1.10, 1.15]:
        for dt in [0.1, 10.0]:
            cb = 1.5 * mat.cb0
            h = 1e-6 * cb
            H, dH, *_ = mat._kinetics_residual(cb, li, False, mat.cb0, B, dt)
            Hp = mat._kinetics_residual(cb + h, li, False, mat.cb0, B, dt)[0]
            Hm = mat._kinetics_residual(cb - h, li, False, mat.cb0, B, dt)[0]
            fd = (Hp - Hm) / (2 * h)
            err = abs(dH - fd) / max(abs(fd), 1e-30)
            worst = max(worst, err)
            print(f"   dH/dcb    lambda={li:.2f} dt={dt:5.1f}: analytic {dH:+.6e}, FD {fd:+.6e}, rel.err {err:.1e}")
    ok = worst < 1e-5
    print(f"   worst relative error {worst:.1e} -> {'PASS' if ok else 'FAIL'}")
    return ok


# ---------------------------------------------------------------------------------------------
# Check 3: single direction
# ---------------------------------------------------------------------------------------------
def direction_state(mat, lam):
    """Kinematics of isochoric uniaxial stretch lam; the x direction then has stretch lam."""
    kin = Kinematics(np.diag([lam, lam ** -0.5, lam ** -0.5]))
    thetaf = mat._solve_thetaf(kin, mat.INITMU)
    return kin, thetaf


def check3():
    print("\n=== Check 3: single direction (x, aligned with uniaxial stretch)")
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    f = np.linspace(0.0, 0.2, 400)
    xdir = np.array([1.0, 0.0, 0.0])
    for name, kinp, col in [('slip only', SLIP, 'tab:blue'), ('catch-slip', CATCH, 'tab:red')]:
        mat = build(kinp)
        koff = np.array([mat._koff(fi)[0] for fi in f])
        ax[0].semilogy(f, 1.0 / koff, color=col, label=name)
        if kinp['kc0'] > 0:
            fstar = KT / (kinp['dxc'] + kinp['dxs']) * np.log(kinp['kc0'] * kinp['dxc'] / (kinp['ks0'] * kinp['dxs']))
            print(f"   catch-slip lifetime peak: numerical f* = {f[np.argmin(koff)]:.4f} pN, analytic f* = {fstar:.4f} pN")
            ax[0].axvline(fstar, color=col, ls=':', lw=1)

        lams = np.linspace(1.0, 1.17, 35)
        cb_ss, f_ss = [], []
        with quiet():
            for lam in lams:
                kin, thetaf = direction_state(mat, lam)
                out = mat._solve_local_kinetics(kin, xdir, mat.cb0, *mat._binding_rate(thetaf), 1.0e6, warn=False)
                cb_ss.append(out[0])
                f_ss.append(out[7])
        cb_ss = np.array(cb_ss)
        ax[1].plot(lams, cb_ss / mat.cb0, color=col, label=name)
        imax = int(np.argmax(cb_ss))
        print(f"   {name:10s}: steady cb/cb0 over lambda 1.00-1.17: min {cb_ss.min()/mat.cb0:.3f}, "
              f"max {cb_ss.max()/mat.cb0:.3f} at lambda {lams[imax]:.3f} (f = {f_ss[imax]:.4f} pN)")

        dts = np.logspace(-2, 3, 21)
        for lam, ls in [(1.05, '-'), (1.10, '--'), (1.15, ':')]:
            dH = []
            with quiet():
                kin, thetaf = direction_state(mat, lam)
                for dt in dts:
                    dH.append(mat._solve_local_kinetics(kin, xdir, mat.cb0, *mat._binding_rate(thetaf), dt, warn=False)[6])
            ax[2].semilogx(dts, dH, color=col, ls=ls, label=f"{name}, lambda={lam}")
            print(f"   {name:10s}: lambda={lam:.2f}: min dH/dcb over dt 0.01-1000 s = {min(dH):.3f}")
    ax[0].set(xlabel='filament force f (pN)', ylabel='bond lifetime 1/koff (s)', title='Lifetime vs force')
    ax[1].set(xlabel='stretch lambda', ylabel='steady-state cb / cb0', title='Steady bound concentration')
    ax[2].axhline(0.1, color='k', lw=0.8, ls='-.', label='dH_min')
    ax[2].set(xlabel='time step dt (s)', ylabel='dH/dcb at root', title='Uniqueness indicator (from cb0)')
    for a in ax:
        a.grid(alpha=0.3)
        a.legend(fontsize=7)
    savefig(fig, 'check3_single_direction.pdf')
    return True


# ---------------------------------------------------------------------------------------------
# Check 4: network response
# ---------------------------------------------------------------------------------------------
def check4():
    print("\n=== Check 4: network response in shear (step and cyclic), slip-only vs catch-slip")
    T = 10.0
    step = lambda t: profile_step(t, t_step=0.5, target_stretch=0.2, deformation_type='shear')
    cyc = lambda t: profile_cyclic(t, amplitude=0.2, frequency=0.2, deformation_type='shear')
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    for name, kinp, col in [('slip only', SLIP, 'tab:blue'), ('catch-slip', CATCH, 'tab:red')]:
        runs = {}
        for dt, ls in [(0.01, '-'), (0.05, '--'), (0.5, ':')]:
            t0 = time.time()
            r = run_profile(build(kinp), 'shear', step, dt, T)
            runs[dt] = r
            ax[0, 0].plot(r['t'], r['sigma'][:, 0, 1], color=col, ls=ls, label=f"{name}, dt={dt}")
            ax[0, 1].plot(r['t'], r['cb_tot'] / r['cb_tot'][0], color=col, ls=ls, label=f"{name}, dt={dt}")
            ax[0, 2].plot(r['t'], r['nsub'], color=col, ls=ls, label=f"{name}, dt={dt}")
            print(f"   step  {name:10s} dt={dt:4.2f}: sigma12(T) = {r['sigma'][-1, 0, 1]:.5f} Pa, "
                  f"cb_tot(T)/cb0 = {r['cb_tot'][-1]/r['cb_tot'][0]:.4f}, max cb_i/cb0 = {r['cb_max'][-1]/r['cb_tot'][0]:.4f}, "
                  f"max sub-steps {r['nsub'].max()}  ({time.time()-t0:.1f} s)")
        ref = runs[0.01]['sigma'][-1, 0, 1]
        for dt in (0.05, 0.5):
            print(f"         {name:10s} dt-convergence: |sigma12(T; dt={dt}) - sigma12(T; dt=0.01)| / |.| = "
                  f"{abs(runs[dt]['sigma'][-1, 0, 1] - ref) / abs(ref):.2e}")

        r = run_profile(build(kinp), 'shear', cyc, 0.05, T)
        gam = np.array([cyc(t) for t in r['t']])
        ax[1, 0].plot(gam, r['sigma'][:, 0, 1], color=col, label=name)
        ax[1, 1].plot(r['t'], r['cb_tot'] / r['cb_tot'][0], color=col, label=f"{name}, total")
        ax[1, 1].plot(r['t'], r['cb_max'] / r['cb_tot'][0], color=col, ls='--', label=f"{name}, max cb_i")
        ax[1, 1].plot(r['t'], r['cb_min'] / r['cb_tot'][0], color=col, ls=':', label=f"{name}, min cb_i")
        ax[1, 2].plot(r['t'], r['nsub'], color=col, label=name)
        peaks = [np.max(r['sigma'][(r['t'] >= k * 5.0) & (r['t'] < (k + 1) * 5.0), 0, 1]) for k in range(2)]
        print(f"   cyclic {name:10s}: peak sigma12 cycle 1 = {peaks[0]:.5f} Pa, cycle 2 = {peaks[1]:.5f} Pa, "
              f"max sub-steps {r['nsub'].max()}")
    ax[0, 0].set(xlabel='time (s)', ylabel='sigma12 (Pa)', title='Step shear (gamma = 0.2 at t = 0.5 s)')
    ax[0, 1].set(xlabel='time (s)', ylabel='cb_tot / cb_tot(0)', title='Total bound crosslinkers')
    ax[0, 2].set(xlabel='time (s)', ylabel='kinetics sub-steps', title='Sub-steps per time step')
    ax[1, 0].set(xlabel='gamma', ylabel='sigma12 (Pa)', title='Cyclic shear (0.2, 0.2 Hz), dt = 0.05 s')
    ax[1, 1].set(xlabel='time (s)', ylabel='cb / cb_tot(0)', title='Bound crosslinkers, cyclic')
    ax[1, 2].set(xlabel='time (s)', ylabel='kinetics sub-steps', title='Sub-steps per time step, cyclic')
    for a in ax.ravel():
        a.grid(alpha=0.3)
        a.legend(fontsize=7)
    savefig(fig, 'check4_network_response.pdf')
    return True


# ---------------------------------------------------------------------------------------------
# Check 5: prestress trend
# ---------------------------------------------------------------------------------------------
def prestrain_run(kinp, gamma0, dgamma, t_hold=60.0, t_obs=30.0, dt_hold=0.5, dt_obs=0.1):
    """Hold shear prestrain gamma0, then (optionally) add a small step dgamma; return t and sigma12."""
    mat = build(kinp)
    mu = mat.INITMU
    t, s = [], []
    with quiet():
        solve_shear(mat, gamma0, mu, 0.0)
        for _ in range(int(round(t_hold / dt_hold))):
            solve_shear(mat, gamma0, mu, dt_hold)
        for k in range(int(round(t_obs / dt_obs)) + 1):
            sigma, *_ = solve_shear(mat, gamma0 + dgamma, mu, 0.0 if k == 0 else dt_obs)
            t.append(k * dt_obs)
            s.append(sigma[0, 1])
    return np.array(t), np.array(s)


def check5():
    print("\n=== Check 5: relaxation of a small step (dgamma = 0.005) superposed on a held shear prestrain")
    dgamma = 0.005
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for name, kinp, col in [('slip only', SLIP, 'tab:blue'), ('catch-slip', CATCH, 'tab:red')]:
        taus = []
        for gamma0, ls in [(0.0, '-'), (0.1, '--'), (0.15, '-.'), (0.2, ':')]:
            t, s_step = prestrain_run(kinp, gamma0, dgamma)
            _, s_base = prestrain_run(kinp, gamma0, 0.0)
            ds = s_step - s_base
            frac = (ds[0] - ds[-1]) / ds[0]
            if abs(frac) < 1e-3:
                tau = np.nan
                txt = "negligible relaxation"
            else:
                target = ds[0] - 0.632 * (ds[0] - ds[-1])
                tau = float(t[np.argmax(ds <= target)]) if frac > 0 else float(t[np.argmax(ds >= target)])
                txt = f"tau63 = {tau:.2f} s"
            taus.append(tau)
            ax[0].plot(t, ds / ds[0], color=col, ls=ls, label=f"{name}, gamma0={gamma0}")
            print(f"   {name:10s} gamma0={gamma0:4.2f}: relaxed fraction of the step stress {frac:+.4f}, {txt}")
        ax[1].plot([0.0, 0.1, 0.15, 0.2], taus, 'o-', color=col, label=name)
        valid = [x for x in taus if np.isfinite(x)]
        if len(valid) >= 2:
            d = np.diff(valid)
            if np.all(d > 0):
                trend = 'slows down monotonically (tau increases)'
            elif np.all(d < 0):
                trend = 'speeds up monotonically (tau decreases)'
            else:
                trend = 'is non-monotonic'
            print(f"   {name:10s}: with increasing prestrain the relaxation {trend}: "
                  f"tau63 = {', '.join(f'{x:.2f}' for x in valid)} s")
    ax[0].set(xlabel='time after step (s)', ylabel='d sigma12 / d sigma12(0)', title='Relaxation of the superposed step')
    ax[1].set(xlabel='prestrain gamma0', ylabel='tau63 (s)', title='Relaxation time vs prestrain')
    for a in ax:
        a.grid(alpha=0.3)
        a.legend(fontsize=7)
    savefig(fig, 'check5_prestress_trend.pdf')
    return True


CHECKS = {1: check1, 2: check2, 3: check3, 4: check4, 5: check5}

if __name__ == '__main__':
    selected = [int(a) for a in sys.argv[1:]] or sorted(CHECKS)
    results = {}
    for k in selected:
        t0 = time.time()
        results[k] = CHECKS[k]()
        print(f"   (check {k} took {time.time() - t0:.1f} s)")
    print("\nSummary:", ", ".join(f"check {k}: {'PASS' if v else 'FAIL'}" if k in (1, 2) else f"check {k}: done"
                                   for k, v in results.items()))
