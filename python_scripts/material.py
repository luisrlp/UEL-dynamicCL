import numpy as np
import math                                                                         
from scipy.optimize import brentq  
from dfd import generate_icosahedron, apply_dispersion
from kinematics import Kinematics

class DynamicCrosslinking:
    def __init__(self, diff_props, mech_props, sim_props):
        self.diff_props = diff_props
        self.mech_props = mech_props
        self.sim_props = sim_props

        # Icosahedron directions and weights
        self.dirs, self.areas = generate_icosahedron(self.sim_props['factor'])
        self.weights = apply_dispersion(self.dirs, self.areas, mech_props['prefdir'], mech_props['b_disp'])
        print(f"Initialized {len(self.dirs)} integration directions.")

        # Equilibrium at t=0
        self.cb0, self.cf0, self.INITMU = self._init_equilibrium()

        # Statev
        self.cb_history = np.full(len(self.dirs), self.cb0)
        self.cf_history = self.cf0
        self.mu_history = self.INITMU
        self.kon_history = np.full(len(self.dirs), self._binding_rate(self.cf0 / self.CF_MAX)[0])
        self.koff_history = np.full(len(self.dirs), self._koff(0.0)[0])
        self.Rbind_history = np.zeros(len(self.dirs))  # Initialize Rbind_history to zeros
        self.f_history = np.zeros(len(self.dirs))
        # Stretch-free binding
        self.stretch_free_binding = self.mech_props.get('stretch_free_binding', False)
        self.xi_t = self.mech_props.get('xi_t', 1.0) if self.stretch_free_binding else 0.0
        if self.stretch_free_binding and self.mech_props.get('filament', 'inext') != 'ext':
            raise NotImplementedError("Stretch-free binding is only implemented for extensible filaments.")
        self.lamn_history = np.ones(len(self.dirs))
        # Local kinetics diagnostics of the last evaluation (convergence flag and dH/dcb per direction)
        self.kin_converged = np.ones(len(self.dirs), dtype=bool)
        self.dHdcb = np.ones(len(self.dirs))
        self.nsub = 1                     # kinetics sub-steps used in the last evaluation
        self.dcb_rel = np.zeros(len(self.dirs))   # max relative cb change per sub-step (last evaluation)

    def _init_equilibrium(self):
        """Chemical equilibrium of the undeformed network (t = 0, F = I; 
        with LAMBDA0 > 1 the extensible filament carries its rest force f0).

        KEQ is defined with the total unloaded off-rate koff(0) = k_off_0 + k_catch_0
        (slip + catch pathways). 
        To keep the same KEQ as a slip-only parameter set,
        split the old k_off_0 between the two pathways instead of adding k_catch_0 on top.
        """

        self.B0 = 1.38065e-5 * self.diff_props['THETA'] * self.mech_props['Lp']
        self.KACTIN = self.diff_props['CACTIN'] * 42.e-3 * 6.022e5 / 16.
        self.CR = self.diff_props['CACTIN'] * self.diff_props['R']
        self.CF_MAX = self.diff_props['CACTIN'] * self.diff_props['RFMAX']
        self.CB_MAX = self.diff_props['CACTIN'] * self.diff_props['RBMAX']
        self.KEQ = self.diff_props['k_on_0'] / self._koff(0.0)[0]
        ext = self.mech_props.get('filament', 'inext') == 'ext'

        def evalh(cb0):
            cf0 = self.CR - cb0
            thetaf0 = cf0 / self.CF_MAX
            # rest force of the prestretched filament (lambda_bar = 1); zero for LAMBDA0 = 1 or the inextensible filament
            f0 = self._f_ext(1.6 * (cb0*1.e3)**(-2/5), 1.0) if (ext and cb0 > 0.0) else 0.0
            # lhs - rhs = 0
            lhs = self.diff_props['k_on_0'] / self._koff(f0)[0] * cf0 * (1.0 - cb0 / self.CB_MAX)
            rhs = cb0 * (1.0 - thetaf0) * math.exp(-self.diff_props['CHI'] * (1.0 - 2.0 * thetaf0))            
            return lhs - rhs

        # Solve for cb0 using Brent's method (same as pullchem.f90)
        cb0_root = brentq(evalh, 0.0, min(self.CR, self.CB_MAX))

        # --- CALCULATE INITMU ---
        cf0 = self.CR - cb0_root
        thetaf0 = cf0 / self.CF_MAX
        JC = 1.0    # J^c = 1 + VMOL*(c - c0) is measured from the initial content c0 = CR (stress-free reference)
        JE = 1.0 / JC

        INITMU = self.diff_props['MU0'] + self.diff_props['RGAS'] * self.diff_props['THETA'] * (
            math.log(thetaf0 / (1.0 - thetaf0)) + 
            self.diff_props['CHI'] * (1.0 - 2.0 * thetaf0) - 
            (self.mech_props['KBULK'] * self.diff_props['VMOL'] / (self.diff_props['RGAS'] * self.diff_props['THETA'])) * (math.log(JE) / JC)
        )

        return cb0_root, cf0, INITMU

    def _binding_rate(self, thetaf):
        """Kon = Kon0 *  exp(chi*(1 - 2 * thetaf)) and 
        Binding rate B = kon(thetaf) * CF_MAX * thetaf/(1 - thetaf) (independent of cb and force)."""
        Kon = self.diff_props['k_on_0'] * np.exp(self.diff_props['CHI'] * (1 - 2*thetaf))
        return Kon, Kon * self.CF_MAX * thetaf / (1 - thetaf)


    def _koff(self, f):
        """Catch-slip off-rate and its force derivative.

        koff(f) = k_catch_0 * exp(-dx_c_kT * f) + k_off_0 * exp(dx_kT * f)
        k_off_0 / dx_kT define the slip pathway, k_catch_0 / dx_c_kT the catch pathway.
        k_catch_0 = 0 reduces to the slip-only Bell law.
        """
        k_s0 = self.diff_props['k_off_0']
        dx_s = self.diff_props['dx_kT']
        k_c0 = self.diff_props.get('k_catch_0', 0.0)
        dx_c = self.diff_props.get('dx_c_kT', 0.0)

        arg_s = dx_s * f
        arg_c = -dx_c * f
        k_s = k_s0 * math.exp(min(max(arg_s, -700.0), 700.0))
        k_c = k_c0 * math.exp(min(max(arg_c, -700.0), 700.0))

        dk_s = dx_s * k_s if abs(arg_s) < 700.0 else 0.0
        dk_c = -dx_c * k_c if abs(arg_c) < 700.0 else 0.0

        return k_s + k_c, dk_s + dk_c

    def _solve_thetaf(self, kinematics, mu):
        J = kinematics.J
        cb_tot = np.sum(self.cb_history * self.weights) # self.weights already accounts for dispersion, area and the consideration of only 1 hemisphere 

        def thetaf_eq(thetaf):
            RT = self.diff_props['RGAS'] * self.diff_props['THETA']
            Jc = 1 + self.diff_props['VMOL'] * (cb_tot + self.CF_MAX * thetaf - self.CR)   # relative to c0 = CR
            Je = J / Jc
            resH = (self.diff_props['MU0'] - mu) / RT + math.log(thetaf / (1. - thetaf)) + self.diff_props['CHI'] * (1. - 2 * thetaf) - ((self.mech_props['KBULK'] * self.diff_props['VMOL']) / RT) * (math.log(Je) / Jc)
            dH = (1. / thetaf) + (1. / (1. - thetaf)) - 2. * self.diff_props['CHI'] + ((self.mech_props['KBULK'] * self.diff_props['VMOL']) / RT) * (self.diff_props['VMOL'] * self.CF_MAX) * (1. + math.log(Je)) / (Jc * Jc)  
            return resH, dH

        iter = 0
        thetaf = self.cf_history / self.CF_MAX 

        rootmin = 1.0e-8
        rootmax = 1.0 - 1.0e-8

        # Check if root is safely bracketed
        resmin, _ = thetaf_eq(rootmin)
        resmax, _ = thetaf_eq(rootmax)
        if resmin * resmax > 0:
            print("Warning: Root is not bracketed. Using initial guess for thetaf.")
            print(f"resmin = {resmin:.6e}, resmax = {resmax:.6e}")

        # Orient the search so that f(xl) < 0 (makes it easier for bracket updating)
        if (resmin < 0.0):
            xl = rootmin
            xh = rootmax
        else:
            xh = rootmin
            xl = rootmax
            swap = resmin
            resmin = resmax
            resmax = swap

        thetaf = np.clip(thetaf, rootmin, rootmax)  # Ensure initial guess is within bounds

        dxold = abs(rootmax - rootmin)
        dx = dxold

        resH, dH = thetaf_eq(thetaf) 

        while iter < self.sim_props['max_iter']:
            # Conditions for activating bisection method
            # 1) Newton is out of bounds. Equivalent to (newroot - xh) * (newroot - xl) > 0.0
            cond1 = ((thetaf - xh) * dH - resH) * ((thetaf - xl) * dH - resH) > 0.0
            # 2) Newton not decreasing fast enough. Equivalent to new_step > 0.5 * old_step
            cond2 = abs(2.0 * resH) > abs(dxold * dH)
            if cond1 or cond2:
                dxold = dx
                dx = 0.5 * (xh - xl)
                thetaf = xl + dx
                if xl == thetaf:
                    break

            else:
                dxold = dx
                dx = resH / dH
                temp = thetaf
                thetaf -= dx
                if temp == thetaf:
                    break

            if abs(dx) < self.sim_props['tol']:
                break

            resH, dH = thetaf_eq(thetaf)

            if resH < 0.0:
                xl = thetaf
            else:
                xh = thetaf

            iter += 1

        # print(f"Converged thetaf = {thetaf:.6f} after {iter} iterations with residual = {resH:.6e} and initial thetaf0 = {self.cf_history/self.CF_MAX:.6f}")
        return thetaf

    def _lambdaf(self, r0f, lambdai):
        return self.mech_props['ETA'] * (r0f + self.mech_props['R0C'])/ r0f * (- 1 + lambdai) + 1

    def _f_inext(self, r0f,lambdaf):
        # If needed, explicitly ensure that lambdaf < a_ratio
        a = self.mech_props['a_ratio']
        aux = self.B0 * np.pi**2 / (r0f * a)**2
        return aux * (((a - 1)/(a - lambdaf))**(1/self.mech_props['BETA']) - 1)

    def _G(self, f, lf, l0f, r0f, L) -> float:
        """
        Computes the difference between the left-hand side (lhs) and right-hand side (rhs)
        of the force equation for an extensible filament.
        """
        mu0 = self.mech_props['MU0STRETCH']
        beta = self.mech_props['BETA']
        lhs = lf * l0f * r0f / L
        aux1 = f / mu0
        aux2 = L**2 / (np.pi**2 * self.B0)
        numerator = (1 + 2 * aux1) * (1 + aux1)**beta * (1 - r0f / L)
        denominator = (1 + f * aux2 + f * aux1 * aux2)**beta
        return lhs - (1 + aux1 - numerator / denominator)
    
    def _f_ext(self, r0f, lf) -> float:
        """
        Computes the force for an extensible filament (with compliant crosslinker), 
        given material properties.
        """
        L = self.mech_props['a_ratio'] * r0f
        l0f = self._lambdaf(r0f, self.mech_props['LAMBDA0'])
        if lf * l0f <= 1.0:
            return 0.0
        f_hi = self.mech_props['MU0STRETCH'] * (1.0 + lf * l0f * r0f / L)
        try:
            return brentq(self._G, 0.0, f_hi, args=(lf, l0f, r0f, L),
                        xtol=1e-14, rtol=1e-13, maxiter=200)
        except ValueError:
            print(f"Root not found for lambdaf = {lf}")
            return np.nan

    def dfdcb(self, f,cbtau,r0f,li,lf):
        a = self.mech_props['a_ratio']
        etac = self.mech_props['ETA']
        r0c = self.mech_props['R0C']
        aux0 = 2/5 * cbtau**-1
        aux1 = 2 * f
        aux2 = etac * r0c / r0f * (li - 1)
        aux3 = self.B0 * np.pi**2 / (r0f * a)**2
        aux4 = (a - lf) * self.mech_props['BETA']
        aux5 = (f + aux3) / aux4
        return aux0 * (aux1 + aux2 * aux5)

    def dfdcb_ext(self,f,cb,r0f,li,lf,l0f,l0,r0):
        ll = self.mech_props['a_ratio'] * r0f
        mu0 = self.mech_props['MU0STRETCH']
        etac = self.mech_props['ETA']
        beta = self.mech_props['BETA']
        a = self.mech_props['a_ratio']
        dfdcb = (2*(-(self.B0*etac*ll*l0f*(-1 + li)*mu0**2*(r0 - r0f)) - self.B0*etac*ll*(-1 +
        l0)*lf*mu0**2*(r0 - r0f) + self.B0*ll*mu0*(l0f*lf*mu0 - (((f +
        mu0)/mu0)**beta*(2*f + mu0)*np.pi**(2*beta))/((f*ll**2*(f +
        mu0))/(self.B0*mu0) + np.pi**2)**beta)*r0f - a*r0f*(2*beta*f*ll**2*(f +
        mu0)*((f + mu0)/mu0)**beta*(2*f + mu0)*np.pi**(2*beta)*((f*ll**2*(f +
        mu0))/(self.B0*mu0) + np.pi**2)**(-1 - beta)*(ll - r0f) + self.B0*l0f*lf*mu0**2*r0f
        - (self.B0*mu0*((f + mu0)/mu0)**beta*(2*f +
        mu0)*np.pi**(2*beta)*r0f)/((f*ll**2*(f + mu0))/(self.B0*mu0) +
        np.pi**2)**beta)))/(5.*self.B0*cb*ll**2*(-mu0 - (beta*ll*((f +
        mu0)/mu0)**beta*(2*f + mu0)**2*np.pi**(2*beta)*((f*ll**2*(f +
        mu0))/(self.B0*mu0) + np.pi**2)**(-1 - beta)*(ll - r0f))/self.B0 + (2*mu0*((f +
        mu0)/mu0)**beta*np.pi**(2*beta)*(ll - r0f))/(ll*((f*ll**2*(f +
        mu0))/(self.B0*mu0) + np.pi**2)**beta) + (beta*((f + mu0)/mu0)**(-1 +
        beta)*(2*f + mu0)*np.pi**(2*beta)*(ll - r0f))/(ll*((f*ll**2*(f +
        mu0))/(self.B0*mu0) + np.pi**2)**beta)))
        return dfdcb

    def _kinetics_residual(self, cb, li, compression, cbt, B, dt):
        """Backward-Euler residual of one direction, H(cb) = cb - cbt - dt*R(cb), and dH/dcb.

        Returns (H, dH, f, Koff, Ri). Beyond the locking stretch of the inextensible filament
        (lambda_f >= a) the force is unbounded, so unbinding is infinitely fast: H is returned
        as +inf, which the bracketed solver reads as 'the root lies at lower cb'.
        """
        ext = self.mech_props.get('filament', 'inext') == 'ext'
        r0f = 1.6 * (cb*1.e3)**(-2/5)
        if ext:
            # Slack filament: total filament stretch (prestretch included) <= 1, evaluated at the
            # trial cb as in kineticsFunc.f90. With LAMBDA0 = 1 this is li <= 1; with LAMBDA0 > 1
            # the filament is loaded at rest (the 'compression' flag is not used here).
            lf = self._lambdaf(r0f, li)
            l0 = self.mech_props['LAMBDA0']
            l0f = self._lambdaf(r0f, l0)
            if lf * l0f <= 1.0:
                f = 0.0
                DfDcb = 0.0
            else:
                f = self._f_ext(r0f, lf)
                r0 = r0f + self.mech_props['R0C']
                DfDcb = self.dfdcb_ext(f, cb, r0f, li, lf, l0f, l0, r0)
        elif compression:
            # Inextensible filament: no prestretch in the force law, slack for li <= 1
            f = 0.0
            DfDcb = 0.0
        else:
            lf = self._lambdaf(r0f, li)
            if lf >= self.mech_props['a_ratio']:
                print(f"Warning: lambdaf = {lf:.6f} exceeds a_ratio = {self.mech_props['a_ratio']:.6f}. Returning H = +inf.")
                return np.inf, np.inf, np.inf, np.nan, np.nan
            f = self._f_inext(r0f, lf)
            DfDcb = self.dfdcb(f, cb, r0f, li, lf)

        # Two-pathway (catch-slip) off-rate; overflow guard is inside _koff
        Koff, dKoff_df = self._koff(f)
        thetab = cb / self.CB_MAX

        Ri = B - Koff * self.CB_MAX * thetab / (1-thetab)
        H = cb - cbt - Ri * dt

        DRiDcb = - (dKoff_df * cb / (1 - thetab) * DfDcb + Koff / (1 - thetab)**2)
        dH = 1 - DRiDcb * dt
        return H, dH, f, Koff, Ri

    def _rtsafe(self, fun, x1, x2, x0, tol, maxiter, H1=None, F0=None):
        """Safeguarded Newton (Numerical Recipes rtsafe, as in solveKinetics.f90).

        fun(x) -> (H, dH). Requires a sign change of H on [x1, x2]. Newton steps from x0 are
        taken while they stay inside the current bracket and decrease fast enough; otherwise
        the bracket is bisected. Returns (x, converged).
        """
        f1 = fun(x1)[0]
        if f1 < 0.0:
            xl, xh = x1, x2
        else:
            xl, xh = x2, x1
        x = min(max(x0, min(x1, x2)), max(x1, x2))
        dxold = abs(x2 - x1)
        dx = dxold
        H, dH = fun(x) if F0 is None else F0
        for _ in range(maxiter):
            if abs(H) < tol:
                return x, True
            newton_ok = (np.isfinite(H) and np.isfinite(dH) and dH != 0.0
                         and ((x - xh)*dH - H) * ((x - xl)*dH - H) < 0.0
                         and abs(2.0*H) <= abs(dxold*dH))
            dxold = dx
            if newton_ok:
                dx = H / dH
                x = x - dx
            else:
                dx = 0.5 * (xh - xl)
                x = xl + dx
            if abs(dx) < tol:
                H, dH = fun(x)
                return x, bool(np.isfinite(H))
            H, dH = fun(x)
            if H < 0.0:
                xl = x
            else:
                xh = x
        return x, False

    def _solve_local_kinetics(self, kinematics, dir, cb_old, Kon, B, dt, warn=True, lamn=1.0):
        """Implicit (Backward-Euler) update of the bound concentration of one direction.

        Catch bonds can make H(cb) non-monotonic (dH/dcb -> 0 or < 0), with several roots.
        The solution must continue the previous state: starting at cb_old, H is marched in
        the direction of the reaction (upwards if R(cb_old) > 0, i.e. H(cb_old) < 0) with
        growing steps up to the first sign change, and the root is then found in that bracket
        with a safeguarded Newton. No intermediate root (other branch) can be skipped, up to
        the resolution of the marching steps. dH/dcb at the root is returned so that the
        caller can reduce the time step where uniqueness is being lost. warn=False silences
        the non-convergence warning (used for attempts that the caller will retry).

        Returns (cb, li_bar, Kon, Koff, Ri, converged, dHdcb, f).
        """
        tol = self.sim_props['tol']
        maxiter = self.sim_props['max_iter']

        li_bar = np.linalg.norm(kinematics.F_iso @ dir)
        li = li_bar / lamn
        ext = self.mech_props.get('filament', 'inext') == 'ext'
        if (li > 1.2) and (not ext):
            print(f"li = {li}")
            print(f"dir = {dir}")
            print(f"F_iso = {kinematics.F_iso}")
        compression = li <= 1.0   # slack test of the inextensible filament (the extensible one tests lf*l0f <= 1)

        last = {}
        def fun(cb):
            out = self._kinetics_residual(cb, li, compression, cb_old, B, dt)
            last['cb'], last['out'] = cb, out
            return out[:2]
        
        lo, hi = 1e-10, self.CB_MAX - 1e-10
        x0 = min(max(cb_old, lo), hi)
        H0, dH0 = fun(x0)
        converged = False
        cbtau = x0
        if abs(H0) < tol:
            converged = True
        elif np.isfinite(H0):
            # March from cb_old in the direction of the reaction up to the first sign change of H
            upwards = H0 < 0.0
            end = hi if upwards else lo
            a, Ha, dHa, s = x0, H0, dH0, 0.05
            while True:
                b = min(a * (1.0 + s), end) if upwards else max(a / (1.0 + s), end)
                Hb, dHb = fun(b)
                if np.isnan(Hb):
                    break
                if Ha * Hb <= 0.0:
                    H1 = Ha if a < b else Hb
                    cbtau, converged = self._rtsafe(fun, min(a, b), max(a, b), a, 
                                                    tol, maxiter, H1, F0=(Ha, dHa))
                    break
                if b == end:
                    break
                a, Ha, dHa, s = b, Hb, dHb, 1.5 * s
        if last.get('cb') == cbtau:
            H, dHdcb, f, Koff, Ri = last['out']
        else:
            H, dHdcb, f, Koff, Ri = self._kinetics_residual(cbtau, li, compression, cb_old, B, dt)
        if not converged and warn:
            print(f"!!! WARNING !!! Local kinetics not converged. cbtau = {cbtau:.10e}, residual = {H:.6e}, "
                  f"dH/dcb = {dHdcb:.3e}, li = {li:.6f}, dir = {dir}")

        return cbtau, li_bar, Kon, Koff, Ri, converged, dHdcb, f

    def _assemble_stress(self, kinematics, cb_array, li_array, thetaf, f_array=None, lamn_array = None):
        # Isochoric contribution
        # SFIC = J^-1 <n w' lf^-1 m otimes m>
        ext = self.mech_props.get('filament', 'inext') == 'ext'
        F_iso = kinematics.F_iso
        detF = kinematics.J

        if ext:
            # NON-VECTORIZED
            sfic_iso = np.zeros((3, 3))
            if f_array is None:
                f_array = np.zeros(len(self.dirs))
                compute_f = True
            else:
                compute_f = False
            lamn_array = np.ones(len(self.dirs)) if lamn_array is None else lamn_array
            xi = self.xi_t if self.stretch_free_binding else 0.0
            l0 = self.mech_props['LAMBDA0']
            for i, dir in enumerate(self.dirs):
                r0f = 1.6 * (cb_array[i]*1.e3)**(-2/5)
                r0 = r0f + self.mech_props['R0C']
                l0f = self._lambdaf(r0f, l0)
                m = F_iso @ dir
                n = self.KACTIN / (self.mech_props['a_ratio'] * r0f)
                dw = 0.0
                # Stretch-free rebinding fraction: relative stretch lambda_bar/lambda_n, dw/dlambda_bar = w'(lambda_rel)/lambda_n
                if xi > 0.0:
                    lf = self._lambdaf(r0f, li_array[i] / lamn_array[i])
                    if lf * l0f > 1.0:
                        if compute_f:
                            f_array[i] = self._f_ext(r0f, lf)
                        dw += xi * f_array[i] * l0 * r0 / lamn_array[i]
                        # dw = f_array[i] * l0 * (r0f + self.mech_props['R0C'])
                # Permanent fraction (without reference update): current stretch
                if xi < 1.0:
                    lf = self._lambdaf(r0f, li_array[i])
                    if lf * l0f > 1.0:
                        if self.stretch_free_binding:
                            f_perm = self._f_ext(r0f, lf)
                        else:
                            if compute_f:
                                f_array[i] = self._f_ext(r0f, lf)
                            f_perm = f_array[i]
                        dw += (1.0 - xi) * f_perm * l0 * r0  
                if dw != 0.0:
                    sfic_iso += (1.0 / detF) * dw * n / li_array[i] * np.outer(m, m) * self.weights[i]

        else:
            r0f_array = 1.6 * (cb_array*1.e3)**(-2/5)
            lf_array = self._lambdaf(r0f_array, li_array)
            f_array = self._f_inext(r0f_array, lf_array) * (li_array > 1.0)
            n_array = self.KACTIN / (self.mech_props['a_ratio'] * r0f_array)
            dw_array = f_array * self.mech_props['LAMBDA0'] * (r0f_array + self.mech_props['R0C'])
            m_array = self.dirs @ F_iso.T
            scalars = (1.0 / detF) * dw_array * n_array / li_array * self.weights
            sfic_iso = np.einsum('i,ij,ik->jk', scalars, m_array, m_array)
        
        # Fictitious to cauchy
        siso = np.einsum('ijkl,kl->ij', kinematics.PROJEUL, sfic_iso)

        # Volumetric contribution
        cb_tot = np.sum(cb_array * self.weights)
        Jc = 1 + self.diff_props['VMOL'] * (cb_tot + self.CF_MAX * thetaf - self.CR)   # relative to c0 = CR
        Je = detF / Jc
        pv = self.mech_props['KBULK'] * math.log(Je) / detF
        svol = pv * np.eye(3)

        return siso, svol, f_array
    
    def evaluate_state(self, F: np.ndarray, mu, dt: float):
        """Takes F and mu, updates internal history, returns Stress.

        mu = None: closed system (no exchange of crosslinkers with the surroundings, e.g. a
        homogeneous rheometer sample). The total crosslinker content is conserved over the step,
            CF_MAX*thetaf + cb_tot^tau(thetaf) = cf^t + cb_tot^t,
        and thetaf is solved from this condition instead of from mu (the corresponding mu is
        stored in self.mu_history). This is what a single element with zero flux enforces in the
        UEL (capacity + reaction = 0). Otherwise (open system), thetaf follows from the given mu.
        """
        kinematics = Kinematics(F)

        # Local kinetics with automatic sub-stepping.
        # The full step is attempted first. The kinetics are recomputed from the start-of-step
        # state with n sub-steps of dt/n (n = 2, 4, 8, ... <= max_substeps) if, in any sub-step,
        #  - any direction does not converge,
        #  - dH/dcb at a root falls below dH_min (uniqueness of the Backward-Euler root being lost),
        #  - any direction changes by more than dcb_max relative to its sub-step start value
        #    (accuracy control: large relative changes mean a large time-discretisation error).
        # F and thetaf are held at their end-of-step values, so every sub-step is still
        # Backward Euler. With n = 1 this is the original single step.
        dH_min = self.sim_props.get('dH_min', 0.1)
        dcb_max = self.sim_props.get('dcb_max', 0.1)
        max_substeps = self.sim_props.get('max_substeps', 64)
        ndir = len(self.dirs)
        li_bar_all = np.linalg.norm(kinematics.F_iso @ self.dirs.T, axis=0)

        def kinetics_sweep(cb_start, dt_sub, Kon, B, lamn, warn):
            """One Backward-Euler kinetics update of all directions over dt_sub."""
            cb = np.zeros_like(cb_start)
            li = np.zeros_like(cb_start)
            kon = np.zeros_like(cb_start)
            koff = np.zeros_like(cb_start)
            R = np.zeros_like(cb_start)
            conv = np.ones(ndir, dtype=bool)
            dH = np.ones(ndir)
            f = np.zeros_like(cb_start)
            for i in range(ndir):
                # print(f"Evaluating direction {i+1}/{len(self.dirs)}: dir = {self.dirs[i]}")
                (cb[i], li[i], kon[i], koff[i], R[i], conv[i], dH[i], f[i]) = self._solve_local_kinetics(
                    kinematics, self.dirs[i], cb_start[i], Kon, B, dt_sub, warn=warn, lamn=lamn[i])
            return cb, li, kon, koff, R, conv, dH, f

        def kinetics_update(thetaf):
            """Kinetics update of all directions over dt (with sub-stepping) for a given thetaf,
            starting from the history at time t. Does not modify the history."""
            Kon, B = self._binding_rate(thetaf)
            nsub = 1
            while True:
                last_try = nsub >= max_substeps
                dt_sub = dt / nsub
                cb_sub = self.cb_history.copy()
                lamn_sub = self.lamn_history.copy()
                conv_all = np.ones(ndir, dtype=bool)
                dH_all = np.full(ndir, np.inf)
                dcb_all = np.zeros(ndir)
                for k in range(nsub):
                    cb_prev = cb_sub
                    if self.stretch_free_binding:
                        # reference stretch: new bonds form at the current stretch
                        lamn_sub = (cb_prev * lamn_sub + dt_sub * B * li_bar_all) / (cb_prev + dt_sub * B)
                    (cb_sub, li_arr, kon_arr, koff_arr, R_arr,
                     conv, dH, f_arr) = kinetics_sweep(cb_sub, dt_sub, Kon, B, lamn_sub, warn=last_try)
                    dcb = np.abs(cb_sub - cb_prev) / cb_prev
                    conv_all &= conv
                    dH_all = np.minimum(dH_all, dH)
                    dcb_all = np.maximum(dcb_all, dcb)
                    if not last_try and (not conv.all() or dH.min() < dH_min or dcb.max() > dcb_max):
                        break                       # abandon this attempt and refine
                ok = conv_all.all() and dH_all.min() >= dH_min and dcb_all.max() <= dcb_max
                if ok or last_try:
                    break
                nsub *= 2
            return (cb_sub, li_arr, kon_arr, koff_arr, R_arr, f_arr, conv_all, dH_all, dcb_all, nsub, ok, lamn_sub)

        if mu is None:
            # Closed system: CF_MAX*thetaf + cb_tot^tau(thetaf) = c_tot (conserved). The left-hand
            # side increases with thetaf (slope CF_MAX + dcb_tot/dthetaf > 0) and is bracketed by
            # thetaf -> 0 (only unbinding: cb_tot^tau <= cb_tot^t < c_tot) and
            # thetaf = c_tot/CF_MAX (left-hand side >= c_tot).
            c_tot = self.cf_history + np.sum(self.cb_history * self.weights)
            memo = {}
            def mass_balance(thetaf):
                # Bound crosslinkers after thetaf update
                memo[thetaf] = kinetics_update(thetaf)
                return self.CF_MAX * thetaf + np.sum(memo[thetaf][0] * self.weights) - c_tot
            lo = 1.0e-12
            hi = min(1.0 - 1.0e-8, c_tot / self.CF_MAX)
            if mass_balance(hi) < 0.0:
                raise ValueError(f"Closed system: total content {c_tot:.6e} exceeds the capacity "
                                 f"(CF_MAX + bound content); increase RFMAX.")
            # Solveing for thetaf based on the mass balance equation
            thetaf = brentq(mass_balance, lo, hi, xtol=1.0e-15, rtol=1.0e-13, maxiter=200)
            update = memo[thetaf] if thetaf in memo else kinetics_update(thetaf)
            # Chemical potential consistent with the solved thetaf (thetaf equation, lagged cb_tot in J^c)
            RT = self.diff_props['RGAS'] * self.diff_props['THETA']
            Jc = 1 + self.diff_props['VMOL'] * (np.sum(self.cb_history * self.weights) + self.CF_MAX * thetaf - self.CR)
            mu_out = self.diff_props['MU0'] + RT * (math.log(thetaf / (1. - thetaf)) + self.diff_props['CHI'] * (1. - 2 * thetaf)) \
                     - self.mech_props['KBULK'] * self.diff_props['VMOL'] * math.log(kinematics.J / Jc) / Jc
        else:
            # Open system: thetaf from the prescribed chemical potential
            thetaf = self._solve_thetaf(kinematics, mu)
            update = kinetics_update(thetaf)
            mu_out = mu

        (cb_sub, new_li_array, new_kon_array, new_koff_array, new_Rbind_array, new_f_array,
         conv_all, dH_all, dcb_all, nsub, ok, new_lamn_array) = update

    
        if nsub > 1 and ok and self.sim_props.get('verbose_substep', False):
            print(f"Sub-stepping: {nsub} sub-steps (min dH/dcb = {dH_all.min():.3e}, "f"max dcb_rel = {dcb_all.max():.3e})")
        if not ok:
            print(f"!!! WARNING !!! Local kinetics not resolved with {nsub} sub-steps: "
                  f"{np.sum(~conv_all)} non-converged directions, min dH/dcb = {dH_all.min():.3e}, "
                  f"max relative cb change per sub-step = {dcb_all.max():.3e}")
        new_cb_array = cb_sub

        # Local kinetics diagnostics of this evaluation used for time-step control
        # Boolean, false for any direction that did not converge (in any sub-step)
        self.kin_converged = conv_all
        # Minimum dH/dcb over the sub-steps, per direction
        self.dHdcb = dH_all
        # Maximum relative cb change in a single sub-step, per direction
        self.dcb_rel = dcb_all
        # Number of kinetics sub-steps used in this evaluation
        self.nsub = nsub
        # print(f"min cb = {np.min(new_cb_array):.10f}, max cb = {np.max(new_cb_array):.6f}, max li = {np.max(new_li_array):.6f}")                                                                                                                                                                     
        # Assemble stress
        siso, svol, f_array = self._assemble_stress(kinematics, new_cb_array, new_li_array, thetaf, 
                                                    f_array=new_f_array, lamn_array = new_lamn_array)
        
        # Update history
        self.cb_history = new_cb_array
        self.cf_history = thetaf * self.CF_MAX
        self.mu_history = mu_out
        self.kon_history = new_kon_array
        self.koff_history = new_koff_array
        self.Rbind_history = new_Rbind_array
        self.f_history = f_array
        self.lamn_history = new_lamn_array
        return siso + svol