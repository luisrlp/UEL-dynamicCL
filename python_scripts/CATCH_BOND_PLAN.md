# Plan: catch–slip crosslinker kinetics in the Python (Gauss-point) implementation

## Context
The Python driver (`python_scripts/material.py` + `main.py`) models force-dependent unbinding with
Bell's slip law only: `Koff = k_off_0 * exp(dx_kT * f)`. The literature reviewed (Lieleg 2009 prestress,
Mulla 2022 α-actinin-4, Schiffhauer 2016, Hosseini 2020) indicates that crosslinkers stabilise under
moderate force (catch) before slipping. main.pdf Eq. 37 introduces a catch–slip extension. The goal is to
add a two-pathway catch–slip off-rate to the Python model first (fast to test, no Abaqus), keep the current
slip model as an exact special case, and make the per-direction implicit kinetics robust to the loss of
monotonicity that catch bonds introduce. This validates the formulation before porting to the Fortran UEL.

Deliverable of this request: this plan saved as `python_scripts/CATCH_BOND_PLAN.md`. No existing file is
edited until implementation is requested.

## Formulation (force form of main.pdf Eq. 37, two-pathway Bell–Evans)
k_off(f) = k_c0·exp(−f·dx_c/kBT) + k_s0·exp(+f·dx_s/kBT)
dk_off/df = −(dx_c/kBT)·k_c0·exp(−f·dx_c/kBT) + (dx_s/kBT)·k_s0·exp(+f·dx_s/kBT)
- k_off(0) = k_c0 + k_s0; lifetime peak at f* = kBT/(dx_c+dx_s)·ln(k_c0·dx_c / (k_s0·dx_s)).
- k_c0 = 0 ⇒ identical to the current Bell law (regression guarantee).
- Units as now: f in pN, dx in µm, kBT = 1.38065e-5·THETA pN·µm.

## Implementation steps (when requested)

### 1. Parameters (`main.py`, `diffusion_props`) [**DONE**]
- Keep `k_off_0` and `dx_kT` as the **slip** pathway (backward compatible), add:
  `k_catch_0` (default 0.0) and `dx_c_kT` (catch reactive distance / kBT).
- `material.py` reads them with `.get(..., 0.0)` so existing scripts run unchanged.

### 2. Off-rate law in one place (`material.py`) [**DONE**]
- New method `_koff(self, f)` → `(Koff, dKoff_df)` implementing the formula above, with exponent
  clipping (`np.clip(..., -100, 100)`, the commented guard already at material.py:192).
- Replace the two hard-coded Bell uses in `_solve_local_kinetics` (material.py:193 and :210):
  - `Koff, dKoff_df = self._koff(f)`
  - `DRiDcb = -(dKoff_df * cbtau/(1-thetab) * DfDcb + Koff/(1-thetab)**2)`
    (generalises the current `Koff*dx_kT*DfDcb` term; compression f=0 ⇒ Koff = k_c0+k_s0).

### 3. Equilibrium / on-rate consistency (`material.py::_init_equilibrium`, line 36) [**DONE**]
- `KEQ = k_on_0 / (k_off_0 + k_catch_0)` (unloaded total off-rate), so the t=0 equilibrium
  (brentq on `evalh`) and INITMU stay consistent. `Kon` (line 166) unchanged.
- Note in docstring: for the same KEQ as before, split the old k_off_0 between the two pathways.

### 4. Robust local kinetics solve (`_solve_local_kinetics`) [**DONE**]
Catch bonds make H(cb) = cb − cb_t − Δt·R non-monotonic (dH/dcb can → 0 or < 0), so plain Newton
(current loop, material.py:180–222) can oscillate or jump branches.
- Replace with safeguarded Newton (rtsafe-style, same logic as Fortran `solveKinetics.f90`):
  bracket [1e-10, CB_MAX − 1e-10]; start from `cb_old`; bisect when Newton leaves the bracket.
- If H has the same sign at both ends or multiple sign changes, pick the root in the sub-interval that
  contains `cb_old` (scan H on a coarse grid ~20 points once; only on failure).
- Return a convergence flag and the value of dH/dcb at the root.

### 5. Automatic sub-stepping (`material.py::evaluate_state`) [**DONE**]
- If any direction reports non-convergence or dH/dcb < ε (e.g. 0.1), split `dt` into n sub-steps
  (n = 2, 4, 8 …, max 64) with F held at the end-of-step value (still backward Euler), updating
  `cb` per sub-step. Keeps `main.py` drivers (secant solves in `solve_uniaxial/biaxial`) untouched.
- Record the number of sub-steps used (new attribute) for diagnostics.

### 6. Diagnostics / outputs [**DONE**]
- Keep `koff_history`, `f_history`; add `dHdcb_history` (min over directions) and `nsub_history`.
- New standalone script `python_scripts/catch_bond_checks.py` (new file) for the verification below;
  no changes to `main.py` plotting blocks.

## Critical files
- `python_scripts/material.py` — `_init_equilibrium` (L29–58), `_solve_local_kinetics` (L163–224),
  `evaluate_state` (L284–311); add `_koff`.
- `python_scripts/main.py` — `diffusion_props` dict (L222–237) only (new keys with defaults).
- Reuse: `dfdcb` (material.py:151, already corrected inextensible derivative), `_f_inext`, `_lambdaf`,
  `Kinematics`, `generate_icosahedron`/`apply_dispersion` (dfd.py). Fortran `solveKinetics.f90` as the
  template for the safeguarded Newton.

## Verification (`catch_bond_checks.py`)
1. **Regression:** `k_catch_0 = 0` ⇒ stress, cb, koff histories identical (to 1e-12) to current code for
   the existing shear/uniaxial/cyclic profiles in `main.py`.
2. **Derivative checks:** finite differences of `_koff` (dKoff/df) and of the full `DRiDcb` vs analytic.
3. **Single-direction kinetics:** constant stretch λ sweep → steady cb,i(λ) and k_off(f): expect
   non-monotonic lifetime (peak at f*) and cb rising then falling; map dH/dcb vs Δt to confirm the
   sub-stepping triggers where uniqueness is lost.
4. **Network response (GP):** step shear / cyclic shear with slip-only vs catch–slip: catch should show
   cb increase in stretched directions and less softening; check stability over dt ∈ {0.01, 0.05, 0.5}.
5. **Qualitative targets:** prestress + small oscillation (Lieleg 2009 Fig. 5, Mulla Ext. Fig. 4e):
   apparent relaxation slows with prestress under catch–slip, speeds up under slip-only.

## Out of scope (later)
- Energy-derivative form of Eq. 37 (∂ψ/∂cb) — needs ∂ψ_i/∂c_b,i at fixed stretch.
- Parameter calibration (Mulla Fig. 1f fit; force-scale reconciliation with extensible filaments).

---

# Part B — Port to the Fortran UEL (Abaqus)

## Goal
Reproduce the validated Python model (Part A) in **`6_CL_DIFF_CS`**. This is a copy of `5_CL_DIFF` made for the catch–slip work (commit 932ade1), and it already uses the **extensible filament** (`fil` + `pullforce`) in both the kinetics residual and the final-state/tangent block. All paths below are relative to `6_CL_DIFF_CS/`.

The port covers:
- the catch–slip off-rate;
- the marching and safeguarded-Newton kinetics solve;
- sub-stepping.

The UEL must stay consistent in two respects:
- **Tangent consistency:** Abaqus needs dc_b/dλ for the Kuu terms, and it must stay exact when sub-stepping is active.
- **Regression:** with `KCATCH0 = 0` and no sub-stepping triggered, results are unchanged with respect to `6_CL_DIFF_CS` **after** the prerequisite fixes F0.

## Current Fortran chain (for reference)

| Stage | Routine | Current behaviour |
|---|---|---|
| Property read | `_umat_.f90` | 24 props; `Koff0 = props(23)`, `Keq = props(24)`, `Kon0 = Koff0*Keq` |
| Initial equilibrium | `pullchem`, `get_initmu.py` | depend only on `Keq` |
| Network loop | `affclnetfic_discrete_imp.f90` | loops over directions; `kon` at line 208 |
| Kinetics solve | `solveKinetics(cbtau_i, args, 17, cbt_i)` | `rtsafe` on the fixed bracket `[1e-8, cbmax − 1e-10]`; calls `exit` when that bracket has no sign change |
| Residual | `kineticsFunc.f90` | calls `fil` (extensible, Brent solve in `pullforce`); `koff = koff0*exp(dx_kT*fi)`; Bell `dRiDcb` |
| Final state and tangent | `affclnetfic_discrete_imp.f90` lines 234–293 | recalls `fil` (line 249, writing `DfDcb` into `dummy_DfDcb`, see F0); Bell `koff_i` (251); `dHdlambda` (269), `dRiDcb` (271), `auxchem = dHdlambda/dHdcb` (273) |
| Filament | `fil.f90` | force from `pullforce` on `[0, 1e9]`; `ddw` = extensible tangent (commit 359335c); `DfDcb` = Mathematica expression (Eq. 80) |

Build facts that matter:
- `compile_elem.sh` and `compile_mat.sh` concatenate every `src/**/*.f90`, so a new file in `src/material_AFFCL/` is compiled automatically.
- The standalone debug build uses `-ffpe-trap=invalid,zero,overflow`: an `exp` overflow or an `Inf` stops the program. No `Inf` may be produced intentionally; the Python `+inf` convention must become a large finite number.

## Parameter semantics
- `KOFF0` (props 23) and `DX` (props 12, a length in µm) stay the **slip** pathway.
- New `KCATCH0` (props 25, 1/s) and `DXC` (props 26, µm) are the catch pathway. `DXC/(kb·θ)` is formed inside the code, as for `DX`.
- `KEQ` stays an **input**, defined with the total unloaded off-rate: `Kon0 = (KOFF0 + KCATCH0)·KEQ`. This is the same definition as Python's `KEQ = k_on_0 / koff(0)`.
- `pullchem`, `get_initmu.py` and `INITMU` do not change, because they depend only on `KEQ`.
- To keep the same kinetics as a slip-only run, split the old `KOFF0` between the two pathways; do not add `KCATCH0` on top of it.

## Steps

### F0. Prerequisite fixes in the extensible-filament path (before any catch–slip change) [**DONE**]
Both bugs are on the active path of `6_CL_DIFF_CS`. Fix them first, and use the fixed code as the regression reference.
1. **`src/material_AFFCL/fil.f90`, `DfDcb` (around line 86):** the Mathematica expression contains `a*r0f`, where `a` was meant as the aspect ratio. Inside `fil`, however, `a` is the local `pullforce` bracket bound, set to `zero`. **Replace `a*r0f` with `ll`** (L = a·r₀f).
   - Python check at c_b = 0.9 µM, λᵢ = 1.02: the correct value is +21.6, the current code gives −71.5.
   - With the fix, the expression matches finite differences of the extensible force, including past the inextensible locking stretch (`material.py::dfdcb_ext`).
   - The wrong `DfDcb` currently enters dH/dc_b (the Newton solve for c_b) and the chemical tangent.
2. **`src/material_AFFCL/affclnetfic_discrete_imp.f90:249`:** `CALL fil(..., cb(node_num), dummy_DfDcb)` writes the derivative into `dummy_DfDcb`. `DfDcb` therefore stays at the zero set on line 202, and it is then used in `ddwdcb` (line 264) and `dRiDcb` (line 271). **Pass `DfDcb`** instead, as the former `fil_inext` call did. `dummy_DfDcb` can then be removed.

### F1. Properties 25–26 (all property-carrying files) 
- **`src/material_AFFCL/_umat_.f90`:**
  - read `Kcatch0 = PROPS(25)` and `dxc = PROPS(26)`;
  - set `Kon0 = (Koff0 + Kcatch0) * Keq` (line 233);
  - pass `Kcatch0` and `dxc` to `affclnetfic_discrete` (new dummy arguments, or append them to `filprops`, whose size changes from 10 to 12). [CREATE NEW CHEM PROPS VECTOR]
- **`affclnetfic_discrete_imp.f90`:** receive them, and set line 208 to `kon = (Koff0 + Kcatch0) * Keq * exp(CHI*(1 − 2θf))`.
- **`test_in_abaqus/sec_uel_cube.inp`:** add a 4th data line `<KCATCH0>,<DXC>`. Abaqus reads 8 values per line, so the lines hold 8 + 8 + 8 + 2 values.
- **`test_in_abaqus/properties.inp` and `get_initmu.py`:** add `KCATCH0 = 0.0` and `DXC = ...`, and keep both files in sync.
- **Every input with `Properties=24`**: change it to `26`. These are the 12 files `uel_cube.inp`, `uel_cube_rotation.inp`, `cube_*_uel*.inp` and `cube_indent_new.inp`. Any input that inlines `*uel property` instead of including `sec_uel_cube.inp` must also get the two extra values.
- **Generators** `2_uel_afm.py`, `3_uel_cube.py`, `4_uel_cube_simple_afm.py`: set `num_properties = 26`.
- **Standalone drivers:**
  - `main.f90`: set `NPROPS=26` and add `PROPS(25) = 0.d0` and `PROPS(26)`;
  - `main_umat.f90`: `NPROPS=30` is already large enough; just set `PROPS(25:26)`.

### F2. Off-rate law in one place (new file `src/material_AFFCL/koffcs.f90`)
```fortran
SUBROUTINE koffcs(koff, dkoffdf, fi, koff0, kcatch0, dxs_kT, dxc_kT)
!  koff = kcatch0*exp(-dxc_kT*fi) + koff0*exp(dxs_kT*fi), and dkoff/dfi
!  Exponents clipped to [-700, 700] (no overflow: -ffpe-trap=overflow);
!  derivative of a clipped pathway set to zero (as material.py::_koff).
```
The clip is essential with the extensible filament: there is no locking stretch, so the force keeps growing with stretch.
- With `DX = 0.005` µm, `dx_kT` ≈ 1.2 /pN, and `exp` overflows above about 580 pN.
- With the Python value `DX = 0.1` µm, it overflows above about 29 pN.
- With `MU0STRETCH = 3.86e13` (the current `properties.inp`), forces past λ_f ≈ a become huge.
Every `koff` in the code must come from this routine, together with `dkoffdf`, which replaces every `dx/(kb*theta) * koff` product. The call sites are:
- `kineticsFunc.f90`: the `koff` line and `dRiDcb`;
- `affclnetfic_discrete_imp.f90`: line 251 (`koff_i`), line 269 (`dHdlambda`) and line 271 (`dRiDcb`).

### F3. Residual (`kineticsFunc.f90`)
- **Arguments:** `nargs` goes from 17 to 19, with `args(18) = kcatch0` and `args(19) = dxc_kT`. Update `nargs` in `affclnetfic_discrete_imp.f90` (line 39) and fill the two new slots next to lines 210–226.
- **Off-rate:** replace the Bell line with `CALL koffcs(koff, dkoffdf, fi, koff0, kcatch0, dx_kT, dxc_kT)`.
- **Derivative:** `dRiDcb = -(dkoffdf*cbtau/(1-thetab)*DfDcb + koff/(1-thetab)**2)`.
- **New output `dHdlambda` (∂H/∂λᵢ):** `dt·cbmax·θb/(1−θb)·dkoffdf·ddwi/(lambda0·r0)`, with `ddwi` from `fil`, the extensible tangent. It is needed for the sub-step tangent recursion (F5). It relies on `ddwi = λ0·r0·df/dλᵢ`, so verify `fil`'s `ddw` against a finite difference of `dw` (Verification 3).
- **No locking guard:** the extensible law has no locking stretch, so `lambdaif >= aratio` needs no special case. A guard is only needed if an inextensible option is re-introduced (see F8).
- **Diagnostics:** remove the `|f| > 1e50` error print, or keep it only behind a debug flag.
- **Cost:** every residual evaluation runs a Brent solve (`pullforce`) and the long `DfDcb` expression. Keep the number of evaluations per direction low (F4: no re-evaluation of known values).

### F4. Kinetics solve (`solvekinetics.f90`), same algorithm as `material.py`
- **New interface:** `solveKinetics(root, args, nargs, rootOld, converged, dHdcb, dHdlambda)`.
- **Algorithm:**
  1. Clip `rootOld` to `[1e-10, cbmax − 1e-10]`, the same bounds as Python, and evaluate H₀. If |H₀| < tol, the direction has converged.
  2. **March** from `rootOld` in the direction of the reaction: upwards if H₀ < 0, downwards otherwise. Use multiplicative steps `s = 0.05`, growing by ×1.5, up to the **first** sign change.
  3. Run `rtsafe` on that bracket, starting from `a`, with its H and dH already known (no re-evaluation, as in the Python fix 3).
  4. Return `converged` and the values of `dHdcb` and `dHdlambda` at the root, taken from the last evaluation.
- **Remove `call exit`:** if no sign change is found before the bound, return `converged = .false.` and `root = rootOld`. With catch bonds, H can have the same sign at both ends of the full interval, so the current global bracket check would abort valid runs.

### F5. Per-direction sub-stepping and a consistent tangent (`affclnetfic_discrete_imp.f90`)
Within one Abaqus increment each direction is independent: λᵢ, θ_f (end of increment) and cb_t are fixed. Sub-stepping can therefore be done **per direction**, which is cheaper than the global Python loop and gives the same Backward-Euler sub-steps. The variables `isub`, `n_sub`, `dtime_sub` and `cb_sub` are already declared.

**New parameters:** `DH_MIN = 0.1d0`, `DCB_MAX = 0.1d0`, `MAX_SUBSTEPS = 64` (implemented as local `PARAMETER`s in `affclnetfic_discrete`, plus `CB_SCALE = 1.0d-3`; nothing added to `global.f90`).

**Sub-step loop** (replaces the single `solveKinetics` call at line 228):
```
n_sub = 1
do                                             ! refinement loop
  cb_sub = cbt_i ; S = 0 ; ok = .true.         ! S = dcb/dlambda_i (sensitivity)
  do isub = 1, n_sub
    args(12) = dtime/n_sub ; args(16) = cb_sub
    CALL solveKinetics(cb_new, args, nargs, cb_sub, conv, dHdcb, dHdlambda)
    S = (S - dHdlambda) / dHdcb                     ! sensitivity recursion (below)
    if (.not.conv .or. dHdcb < DH_MIN .or. abs(cb_new-cb_sub)/cb_sub > DCB_MAX) ok = .false.
    cb_sub = cb_new
    if (.not.ok .and. n_sub < MAX_SUBSTEPS) exit
  end do
  if (ok .or. n_sub >= MAX_SUBSTEPS) exit
  n_sub = 2*n_sub
end do
cbtau_i = cb_sub
```
Keep `min dHdcb`, `max n_sub` and a non-convergence flag over all directions for F6 and F7.

**Consistent tangent: sensitivity recursion.** Sub-step k solves H_k(cb_k, cb_{k−1}, λ) = cb_k − cb_{k−1} − Δt_s·R(cb_k, λ) = 0. Differentiating gives

  S_k = dcb_k/dλ = (S_{k−1} − ∂H_k/∂λ) / (∂H_k/∂cb), with S_0 = 0 (cb_t does not depend on λ at t + Δt).

- For n_sub = 1 this reduces to S = −dHdlambda/dHdcb, which is today's `auxchem` with the opposite sign.
- So set `auxchem = −S` (line 273), keep the existing sign conventions of the `cfic` and `dcbdC` terms, and stop recomputing `dHdlambda` and `dHdcb` with the full `dtime` (lines 269–272).
- Without this, sub-stepped increments have a wrong Kuu and lose quadratic convergence.

**Final-state block** (lines 234–293): keep it, with two changes:
- `koff_i` and its derivative come from `koffcs`;
- `DfDcb` (returned by `fil` after F0.2) is still needed for `ddwdcb`.

Also remove the `lambdaif >= 1.18` debug print (line 252). With the extensible filament it fires in normal operation.

**Accuracy criterion: problem found and fix (implemented, Fortran only).**
- **Problem:** the relative-change check `|Δcb|/cb_start > DCB_MAX` misfires when the kinetics are much faster than the sub-step (Δt_sub·k_off ≫ 1). Example: λᵢ ≥ 1.3, or complete unbinding at the lower bound.
  - c_b then drops to its (near-zero) equilibrium within any sub-step, so the relative change is always ~100%, and every such direction hits the 64 cap at every Δt.
  - Backward Euler is actually accurate there, since it returns the quasi-equilibrium. But the cap would trigger endless `PNEWDT` cut-backs (F6).
  - Exempting roots at the bound (option a) or using a mixed tolerance `DCB_MAX·max(cb_start, CB_SCALE·cbmax)` (option b) does not fix it. A drop from a large to a near-zero value fails the first sub-step at any n_sub.
- **Fix:** estimate the Backward-Euler error of each sub-step from the stiffness the solver already provides, kh = Δt_sub·|∂R/∂cb| = |dH/dcb − 1|:
  - `err = |Δcb|·min(kh/2, 1/kh)`, which is exact in order for linear kinetics: ~|Δcb|·kh/2 for slow kinetics, ~|Δcb|/kh → 0 at quasi-equilibrium;
  - kh is evaluated at the **start and end** of the sub-step, and the larger factor is used (one extra `kineticsFunc` call). The end-of-step rate alone accepted steps with errors up to 400% in strongly nonlinear transients;
  - sub-step if `err > DCB_MAX·max(min(cb_start, cb_end), CB_SCALE·cbmax)`. The tolerance is relative to the **smaller** value, because normalising by the start value accepted a 13%-of-start error that was 30% of the result.
- **Verification:** 112 cases against a 4096-sub-step reference (slip and strong catch, Δt 0.05–1000 s, λᵢ 1.02–1.5):
  - every accepted result is within 9.8% of the reference;
  - the fast-equilibrium cases now need `n_sub = 1`;
  - moderate cases still need no sub-stepping at Δt ≤ 1 s;
  - the tangent sensitivity still matches finite differences (≤ 1.5e-4).
- **Remaining:** 11 of the 112 cases hit the cap even though their results are accurate (errors 1e-7 to 1.4e-2), because per-sub-step error control is conservative in the fast transient after a stretch jump. With `MAX_SUBSTEPS = 1024`, 5 remain. So F6 should cut back only on hard failures (no convergence, or dH/dcb < `DH_MIN`) and only warn when the cap is reached on the accuracy check.
- The Python (`material.py`) still uses the old relative-change criterion, so the two implementations differ whenever sub-stepping triggers.

### F6. Increment cut-back (PNEWDT)
- **Pass `PNEWDT`** from `_umat_.f90` (it is already a dummy argument of `material`) to `affclnetfic_discrete` as a new `INTENT(INOUT)` argument.
- **When a direction is unresolved** (not converged at `MAX_SUBSTEPS`): set `PNEWDT = min(PNEWDT, 0.5d0)` and write **one** warning per call, giving the element, the integration point and the number of unresolved directions. Abaqus then repeats the increment with a smaller Δt. This replaces the old `call exit`.
- **Do not** cut back just because sub-stepping was needed. `PNEWDT < 1` discards a converged increment, and sub-stepping exists to avoid that.

### F7. Diagnostics in STATEV (optional, recommended) [SKIPPED]
- **`global.f90`:** `NSDV = 16 + NDIR` becomes `18 + NDIR`. `sdvread` and `initialize` address c_b as `nsdv − ndir + i`, so the c_b slots move automatically.
- **`sdvwrite.f90`:** slot 17 = max `n_sub` over directions, slot 18 = min `dHdcb`. Update the layout comment.
- **Inputs:** since F1 already touches all of them, change these at the same time:
  - `Variables = 8·NSDV`: 848 becomes 864 for FACTOR = 3;
  - `*User output variables`: 106 becomes 108;
  - the `nlSdv, ngSdv` line in `sec_uel_cube.inp`: `106, 106` becomes `108, 108`;
  - `uel_variables` in the generators.
- `uvarm.f90` loops over `nsdv`, so it needs no change.
 
### F8. Extensible-filament follow-ups (remaining items) [SKIPPED]
`6_CL_DIFF_CS` already uses `fil`, and its `DfDcb` bug is fixed in F0. What remains:
- **Force bracket (`fil.f90`):** `pullforce` is called on `[0, 1e9]`.
  - When the root lies above 1e9, the bracket fails silently. This is possible with a large `MU0STRETCH` (the current `properties.inp` has 3.86e13) past λ_f ≈ a.
  - Python uses `[0, μ₀(1 + λf·λ0f·r0f/L)]`, on which G < 0 is guaranteed. Use the same upper bound, or at least check that G(b) < 0.
  - Also return f = 0 when `λf·λ0f ≤ 1` (G(0) ≤ 0, no bracket).
- **μ₀ value:** `properties.inp` has `MU0STRETCH = 3.86e13`, which is effectively inextensible, while `main.f90` has `PROPS(8) = 38600`. If 3.86e13 is intentional, to mimic the inextensible law, the extensible behaviour is effectively off in the Abaqus runs. Decide the value before the Python comparison (Verification 2).
- **Optional law switch:** props 27 `FILMODEL` (0 = inextensible, 1 = extensible), or a PARAMETER in `global.f90`, to choose `fil` or `fil_inext` in `kineticsFunc.f90` and in the final-state block. It is only needed for an inextensible regression inside the UEL. With `FILMODEL = 0`, add the locking guard:
  - if `lambdaif >= aratio`, return `f = 1.0d30`, `df = 0` and `dHdlambda = 0`, without calling `fil_inext`;
  - the large finite residual means "root at lower c_b" without producing `Inf`, and `df = 0` forces bisection.

## Critical files (Part B, all in `6_CL_DIFF_CS/`)
| File | Steps |
|---|---|
| `src/material_AFFCL/fil.f90` | F0.1, F8 |
| `src/material_AFFCL/_umat_.f90` | F1, F6 |
| `src/material_AFFCL/affclnetfic_discrete_imp.f90` | F0.2, F1, F2, F3 (args), F5, F6 |
| `src/material_AFFCL/kineticsFunc.f90` | F2, F3 (F8) |
| `src/material_AFFCL/solvekinetics.f90` | F4 |
| `src/material_AFFCL/koffcs.f90` (new) | F2 |
| `global.f90` | F5 parameters, F7 `NSDV` |
| `src/material_AFFCL/sdvwrite.f90` | F7 |
| `main.f90`, `main_umat.f90` | F1 |
| `test_in_abaqus/sec_uel_cube.inp`, `properties.inp`, `get_initmu.py` | F1 (F7) |
| All `test_in_abaqus/*.inp` with `Properties=24`, and generators `2_/3_/4_*.py` | F1 (F7) |

## Verification (Part B)
1. **Regression:** first apply F0 and re-run `uel_cube.inp` and `cube_shear_uel_auto.inp` to create the reference UVARM histories. Outputs from before F0 are not a valid reference, since F0 changes the tangent and dH/dc_b. After F1–F7, with `KCATCH0 = 0` and sub-stepping not triggered, the results must match that reference to round-off.
2. **Single element vs Python:** run a homogeneous cube in shear or uniaxial with identical parameters, using `properties.inp` values equal to `catch_bond_checks.build`, the same factor and the same Δt.
   - Python must use **`'filament': 'ext'`** with the same `MU0STRETCH` as the UEL.
   - c_b(t) per direction and σ must match `material.py` (slip and catch–slip) to solver tolerance.
   - Differences can come only from Abaqus's own iteration on F, from the reference swelling (J꜀) at t = 0, or from the `pullforce` tolerance (`t = 1e-14`) against `brentq`.
3. **Derivative and tangent checks:**
   - standalone finite-difference checks, in `main_umat.f90` or a small driver, of `fil`'s `DfDcb` (after F0.1) and `ddw` (against d(dw)/dλᵢ), and of `kineticsFunc`'s `df` and `dHdlambda`;
   - watch the Abaqus Newton convergence rate in the `.msg` file. It must stay quadratic, including in increments where `n_sub > 1` (slot 17 > 1).
4. **Robustness:**
   - with a strong catch parameter set and large increments, the run must complete, through sub-stepping or `PNEWDT` cut-backs, without calling `exit`;
   - run the standalone debug build (`-ffpe-trap`) at stretches above the old locking value (λ_f > a, force in the hundreds of pN) without triggering an FPE (`koffcs` clip, `pullforce` bracket).
5. **Objectivity:** `uel_cube_rotation.inp` with the catch–slip parameters must keep all quantities constant during the rotation.

# Part C — Consistent chemical tangent ($K_{cu}$, $K_{cc}$) in the UEL

## Finding
A central finite-difference check of `AMATRX` showed that one block is essentially missing. The check perturbed each of the 32 element DOFs, resetting `SVARS` before every call, in the standalone driver. The state was one element rotated 13.5–90° about $z$, with a small general stretch, a $\mu$ gradient, and the `properties.inp` parameters (catch law, matrix on). `AMATRX` $= -\partial \mathrm{RHS}/\partial U$, the Abaqus convention, block by block:

| Block | $\lVert \mathrm{FD} \rVert$ | Relative error of `AMATRX` |
|---|---|---|
| $K_{uu}$ | $1.5\times10^{3}$ | $5\times10^{-9}$ |
| $K_{uc}$ | $1.5\times10^{-5}$ | $3\times10^{-4}$ |
| $K_{cc}$ | $6.2\times10^{-7}$ | $2.5\times10^{-3}$ |
| $\mathbf{K_{cu}}$ | $\mathbf{2.7\times10^{-4}}$ | **1.00**: `AMATRX` returns $6\times10^{-9}$ |

- **Converged solutions are unaffected,** since the tangent only changes the Newton path.
- **But the coupled Newton iteration is not quadratic.** In zero-flux problems, where Abaqus accepts increments without a flux reference, the leftover error builds up. Example: the rigid-rotation test, where $\mu$ drifts by 0.1 J/mol and gives an isotropic stress of $3\times10^{-6}$ Pa, decreasing with the rotation increment.

## Residual and its dependencies (per integration point)
Chemical residual of node $a$, integrated over the current volume (`detMapJC`; `dshC` $= \partial N/\partial \mathbf{x}$ at the integration point):

$$
R_{c,a} = \int_v \left[ N_a\,\mathrm{ResFac} - M\,\nabla_x N_a \cdot \nabla_x \mu \right] dv,
\qquad
\mathrm{ResFac} = -\frac{c_{f,\max}\,\dot\theta_f + R_{\mathrm{macro}}}{\bar J}
$$

**F-bar.** The element is the fully integrated 8-node brick, so `u3d8.f90` (lines 344–355) replaces $\mathbf F$ at each integration point by $\bar{\mathbf F} = (\bar J/J)^{1/3}\,\mathbf F$, where $J = \det\mathbf F$ is local and $\bar J$ = `detFc_tau` is the determinant at the element centre. The material receives $\bar{\mathbf F}$ and recomputes its determinant from it (`fslip`, `_umat_.f90:269`), so every material quantity uses $\bar J$: $\theta_f$, the mobility, and `detF` in `ResFac`. But the residual is integrated over the local volume, $dv = J\,dV_0$. The isochoric part (`distgr`) is the same for $\bar{\mathbf F}$ and $\mathbf F$, so the stretches are local.

*Notation in Part C:* $\bar J$ and $\bar{\mathbf F}$ are F-bar quantities. $\mathbf F_{\mathrm{iso}} = J^{-1/3}\mathbf F$ is the isochoric part, and $\bar\lambda_i = \lvert\mathbf F_{\mathrm{iso}}\,\mathbf m_{0,i}\rvert$ keeps its usual overbar for the isochoric stretch.

- **The source term** depends on the displacements through $\theta_f(\mu,\bar J)$, $c_{b,\mathrm{tot}}(\mathbf F_{\mathrm{iso}},\theta_f)$ and the volume ratio $J/\bar J$:
  $$
  \int_v N_a\,\mathrm{ResFac}\,dv = -\int_{V_0} N_a\,\frac{J}{\bar J}\,\frac{c_{f,\max}\,\Delta\theta_f + \Delta c_{b,\mathrm{tot}}}{\Delta t}\,dV_0
  $$
  - $\theta_f$ is solved first, in `_umat_.f90`, with $c_{b,\mathrm{tot}}$ at time $t$ and `DET` $= \bar J$. Its only dependence is $\theta_f(\mu, \bar J)$, with
    $$
    \frac{\partial\theta_f}{\partial \bar J} = \frac{K\,V_{\mathrm{mol}}}{RT\,J^c\,\bar J\,d_f}
    \quad\text{(already in the DMDJ comment)},
    \qquad
    \frac{\partial\theta_f}{\partial\mu} = \texttt{DTHETAFDMU}.
    $$
  - $c_{b,\mathrm{tot},\tau} = \sum_i \rho\,a_i\,c_{b,i}$ depends on the isochoric stretches $\bar\lambda_i = \lvert \mathbf F_{\mathrm{iso}}\,\mathbf m_{0,i} \rvert$, with $\mathbf F_{\mathrm{iso}} = J^{-1/3}\mathbf F$ (`distgr`), and on $\theta_f$ (through $k_{on}$ and the free fraction).
- **Time integration (Backward Euler, with one lag).** The residual imposes $(c^{\tau}-c^{t})/\Delta t = -\mathrm{Div}\,\mathbf j_R^{\tau}$ for the total content $c = c_f + c_{b,\mathrm{tot}}$. Here $\theta_f^{\tau}$ is solved with $\mu^{\tau}$ and $\bar J^{\tau}$, $c_{b,\mathrm{tot}}^{\tau}$ comes from the Backward-Euler kinetics with $\theta_f^{\tau}$, and the flux uses $\mu^{\tau}$ and $m^{\tau}$. The update is conservative in time: the total change telescopes over increments. The exception is that $\theta_f^{\tau}$ is solved first, with $J^c = 1 + V_{\mathrm{mol}}(c_{b,\mathrm{tot}}^{\,t} + c_{f,\max}\theta_f^{\tau} - c_0)$, i.e. with the bound content at time $t$ (`ARGS(9) = CB_TOT`). $J^c$ is not recomputed after the kinetics update. So the $\mu$–$\theta_f$ relation and the volumetric pressure lag $c_b$ by one increment: a staggered, first-order-consistent update, with an error of order $K\,V_{\mathrm{mol}}\,\Delta c_{b}$ per increment ($\sim10^{-4}$ Pa here). The tangent terms below differentiate exactly this algorithm. In particular there is no $\partial J^c/\partial c_{b,\mathrm{tot}}^{\tau}$ term, because $J^c$ uses $c_{b,\mathrm{tot}}^{t}$. A fully implicit variant is described at the end of Part C.
- **The flux term** $-\int_v M\,\nabla_x N \cdot \nabla_x \mu\,dv$ is already in $K_{cu}$ (`SpCUModFac`, the geometric and mobility part) and in $K_{cc}$.

## Missing terms
With `AMATRX` $= -\partial R/\partial U$, and DOF $b$ in direction $k$.

**Relation to `literature/implementation_plan.md` (§6, "Residuals and Tangents").** That plan neglects most of these terms on purpose, as a "Modified-Newton approximation", to keep the tangent symmetric:

| Term in `implementation_plan.md` | Treatment there | Here |
|---|---|---|
| $\mathbf K_{\mu_f\mathbf u}$ (I) Capacity, $\partial\dot c_f/\partial\mathbf u$ | neglected | 1a, plus the capacity part of 1b (F-bar, not considered there) |
| $\mathbf K_{\mu_f\mathbf u}$ (II) Diffusion | partly: the dilation and shape-distortion parts are dropped, and the mobility shift is projected onto $\nabla N^A\cdot\nabla N^B$ | item 5 |
| $\mathbf K_{\mu_f\mathbf u}$ (III) Reaction, $\partial\mathcal R/\partial\mathbf u$ | neglected | 1c (stretches), plus the $G$ part of 1a ($\theta_f$ through $\bar J$) and the reaction part of 1b |
| $\mathbf K_{\mu_f\mathbf u}$, surface flux on a deforming face | not discussed ("No change to the tangent matrix" in the code) | item 6 |
| $\mathbf K_{\mu_f\mu_f}$ Term 2, Reaction, $\partial\mathcal R/\partial\mu_f \approx 0$ | neglected | item 3 |
| $\mathbf K_{\mathbf u\mu_f}$, kinetic–mechanical coupling, $\partial c_b/\partial\theta_f \approx 0$ | neglected | item 4 |

The symmetry argument does not apply: the element is declared `Unsymm` (`*User Element, …, Unsymm` in every input), so Abaqus already uses the unsymmetric solver.

**Capacity (I) and reaction (III) in `ResFac`.** The code combines terms (I) and (III) of `implementation_plan.md` into one integrand:
$$
\underbrace{-\int_v N_a\,\frac{\dot c_f}{\bar J}\,dv}_{(I)}\;\underbrace{-\int_v N_a\,\frac{\mathcal R}{\bar J}\,dv}_{(III)}
= \int_v N_a\,\mathrm{ResFac}\,dv,
\qquad
\mathrm{ResFac} = -\frac{\dot c_f + \mathcal R}{\bar J},
\quad \dot c_f = \texttt{CFMAX}\cdot\texttt{DTHETAFDT},
\quad \mathcal R = \texttt{RMACRO}.
$$
Both the residual and its derivative are linear in the two contributions, so the code keeps them combined (`ResFac`, and a single coefficient `DSRCDJ`). Their origin in each tangent term below:

| Term | From (I), capacity | From (III), reaction |
|---|---|---|
| 1a ($\theta_f$ through $\bar J$) | $\dfrac{c_{f,\max}}{\Delta t}\,\dfrac{\partial\theta_f}{\partial\bar J}$ | $\dfrac{G}{\Delta t}\,\dfrac{\partial\theta_f}{\partial\bar J}$ (binding responds to $\theta_f$) |
| 1b (F-bar ratio $J/\bar J$) | $\dfrac{\dot c_f}{\bar J}\left(\dfrac{\partial N_b}{\partial x_k} - \left.\dfrac{\partial N_b}{\partial x_k}\right\vert_c\right)$ | $\dfrac{\mathcal R}{\bar J}\left(\dfrac{\partial N_b}{\partial x_k} - \left.\dfrac{\partial N_b}{\partial x_k}\right\vert_c\right)$ |
| 1c (stretches) | none | all of it ($\mathbf Q$) |
| Item 3 ($K_{cc}$, `TanFac`) | $\dfrac{c_{f,\max}}{\Delta t\,\bar J}\,\dfrac{\partial\theta_f}{\partial\mu}$ (already in the code) | $\dfrac{G}{\Delta t\,\bar J}\,\dfrac{\partial\theta_f}{\partial\mu}$ (new) |

In 1b, $-\mathrm{ResFac} = (\dot c_f + \mathcal R)/\bar J$, so the two parts are its capacity and reaction contributions.

1. **$K_{cu}$, source term (`ResFac`: capacity + reaction).** One derivation (below) gives three parts:
   $$
   K^{(1a)}_{cu,\,a(bk)} = \int_v N_a\,\frac{c_{f,\max} + G}{\Delta t}\,\frac{\partial\theta_f}{\partial \bar J}\,\left.\frac{\partial N_b}{\partial x_k}\right|_{c}\,dv,
   \qquad
   K^{(1b)}_{cu,\,a(bk)} = -\int_v N_a\,\mathrm{ResFac}\left(\frac{\partial N_b}{\partial x_k} - \left.\frac{\partial N_b}{\partial x_k}\right|_{c}\right) dv,
   $$
   $$
   K^{(1c)}_{cu,\,a(bk)} = \int_v N_a\,\frac{1}{\bar J\,\Delta t}\,\sum_l Q_{kl}\,\frac{\partial N_b}{\partial x_l}\,dv,
   $$
   where $\partial N_b/\partial x_k$ is taken at the integration point (`dshC`) and $\left.\partial N_b/\partial x_k\right|_c$ at the element centre, in the current configuration (`dshC0`, already computed in `u3d8.f90`, line 215).
   - **1a, free fraction through $\bar J$** (the dominant part at small strain): a volume change changes $\theta_f$, hence the free content (capacity, $c_{f,\max}$) and the binding (reaction, $G$).
   - **1b, F-bar volume ratio:** the derivative of the factor $J/\bar J$. It vanishes without F-bar (where $J = \bar J$) and at chemical equilibrium ($\mathrm{ResFac} = 0$).
   - **1c, binding through the stretches** (important at large strain): $\mathbf Q = \partial c_{b,\mathrm{tot}}^{\tau}/\partial(\nabla_x\mathbf u)$ is the binding sensitivity tensor derived in step 4 (not a stress):
     $$
     \mathbf Q = \sum_i \rho\,a_i\,S_i\,\bar\lambda_i \left( \hat{\mathbf m}_i \otimes \hat{\mathbf m}_i - \tfrac13\,\mathbf I \right),
     $$
     with $S_i = \partial c_{b,i}/\partial\bar\lambda_i$ the existing sub-step sensitivity (`sens_i`). This is exactly `−dcbdC` as already accumulated in `affclnetfic_discrete`: $\lambda^{-1}\cdot$`auxchem`$\cdot$`pfdlambdadcfil`, with `auxchem = −sens_i` and `pfdlambdadc` $= \rho\,a\,(\mathbf m \otimes \mathbf m - \tfrac13\lambda^2\mathbf I)$ for the non-normalised $\mathbf m = \mathbf F_{\mathrm{iso}}\,\mathbf m_0$ (`mfi`, from `distgr`, which is local and unaffected by F-bar). Confirm the sign with the FD check.

   The three parts are listed separately because they are tested in different FD states (1c only at large strain, 1b only under non-homogeneous deformation) and use different coefficients and gradients in the code (C3).

   **Derivation.**
   1. *Source part of the residual.* With $\dot\theta_f$ and $R_{\mathrm{macro}}$ discretised by Backward Euler, and $dv = J\,dV_0$ ($J$ local):
      $$
      R^{\mathrm{src}}_{c,a} = \int_v N_a\,\mathrm{ResFac}\,dv
      = -\int_v N_a\,\frac{s}{\bar J}\,dv
      = -\int_{V_0} N_a\,\frac{J}{\bar J}\,s\,dV_0,
      \qquad
      s = \frac{c_{f,\max}\,(\theta_f^{\tau}-\theta_f^{t}) + (c_{b,\mathrm{tot}}^{\tau}-c_{b,\mathrm{tot}}^{t})}{\Delta t}.
      $$
      On the reference volume, $N_a$, $V_0$ and the values at $t$ do not depend on the current displacements. Only $J/\bar J$, $\theta_f^{\tau}$ and $c_{b,\mathrm{tot}}^{\tau}$ do.
   2. *Tangent.* With `AMATRX` $= -\partial R/\partial U$ and the nodal displacement $u_{bk}$ (node $b$, direction $k$):
      $$
      K^{\mathrm{src}}_{cu,\,a(bk)} = -\frac{\partial R^{\mathrm{src}}_{c,a}}{\partial u_{bk}}
      = \int_{V_0} N_a \left[ \frac{J}{\bar J}\,\frac{\partial s}{\partial u_{bk}} + s\,\frac{\partial}{\partial u_{bk}}\!\left(\frac{J}{\bar J}\right) \right] dV_0,
      \qquad
      \frac{\partial s}{\partial u_{bk}} = \frac{1}{\Delta t}\left[ c_{f,\max}\,\frac{\partial\theta_f^{\tau}}{\partial u_{bk}} + \frac{\partial c_{b,\mathrm{tot}}^{\tau}}{\partial u_{bk}} \right].
      $$
   3. *Free fraction.* In `_umat_.f90`, $\theta_f^{\tau}$ solves $H(\theta_f;\mu,\bar J) = 0$ (`thetafFunc`, with `DET` $= \bar J$ and $c_{b,\mathrm{tot}}$ held at time $t$):
      $$
      H = \frac{\mu_0-\mu}{RT} + \ln\frac{\theta_f}{1-\theta_f} + \chi\,(1-2\theta_f) - \frac{K\,V_{\mathrm{mol}}}{RT}\,\frac{\ln(\bar J/J^c)}{J^c},
      \qquad J^c = 1 + V_{\mathrm{mol}}\left(c_{b,\mathrm{tot}}^{t} + c_{f,\max}\,\theta_f - c_0\right).
      $$
      At fixed nodal $\mu$, $\theta_f$ depends on the displacements only through $\bar J$. By the implicit function theorem, with $J^c$ independent of $\bar J$ at fixed $\theta_f$:
      $$
      \frac{\partial H}{\partial \bar J} = -\frac{K\,V_{\mathrm{mol}}}{RT}\,\frac{1}{\bar J\,J^c},
      \qquad
      \frac{\partial H}{\partial\theta_f} = d_f \ (\texttt{df}),
      \qquad
      \frac{\partial\theta_f}{\partial \bar J} = -\frac{\partial H/\partial \bar J}{\partial H/\partial\theta_f} = \frac{K\,V_{\mathrm{mol}}}{RT\,\bar J\,J^c\,d_f}.
      $$
   4. *Bound content.* The total bound content at the end of the increment is the weighted sum over directions, and each direction's bound content depends on two quantities that change with the displacements: the free fraction $\theta_f$ (binding rate, through $k_{on}$ and $\theta_f/(1-\theta_f)$) and its isochoric stretch $\bar\lambda_i$ (force-dependent unbinding):
      $$
      c_{b,\mathrm{tot}}^{\tau} = \sum_i \rho\,a_i\,c_{b,i}^{\tau}\big(\theta_f,\ \bar\lambda_i\big).
      $$
      *Chain rule.*
      $$
      \delta c_{b,\mathrm{tot}}^{\tau}
      = \underbrace{\sum_i \rho\,a_i\,\underbrace{\frac{\partial c_{b,i}^{\tau}}{\partial\theta_f}}_{S^{\theta}_i}}_{G}\;\delta\theta_f
      \;+\; \sum_i \rho\,a_i\,\underbrace{\frac{\partial c_{b,i}^{\tau}}{\partial\bar\lambda_i}}_{S_i}\;\delta\bar\lambda_i ,
      $$
      where $S^{\theta}_i$ is the sensitivity of each direction to the free fraction (new, `sensthetaf_i`) and $G$ its weighted sum (`dcbdthetaf`); $S_i$ is the existing stretch sensitivity of each direction (`sens_i`). Both are computed by the same sub-step recursion (item 3).
      - *First term (free fraction).* $\delta\theta_f = (\partial\theta_f/\partial\bar J)\,\delta\bar J$ (step 3), with $\delta\bar J$ from step 5. This is the $G$ part of 1a.
      - *Second term (stretches).* The isochoric stretch of direction $i$ is $\bar\lambda_i^2 = \mathbf m_{0,i}\cdot\bar{\mathbf C}\,\mathbf m_{0,i}$, with $\bar{\mathbf C} = J^{-2/3}\mathbf C$ and the local $J$, since `distgr` is local. Its variation is
        $$
        \delta(\bar\lambda_i^2) = J^{-2/3}\,\mathbf m_{0,i}\cdot\delta\mathbf C\,\mathbf m_{0,i} - \tfrac23\,\bar\lambda_i^2\,\frac{\delta J}{J},
        \qquad
        \delta\mathbf C = \mathbf F^{T}\left(\nabla_x\delta\mathbf u + \nabla_x\delta\mathbf u^{T}\right)\mathbf F .
        $$
        With $\mathbf F\,\mathbf m_{0,i} = J^{1/3}\,\bar\lambda_i\,\hat{\mathbf m}_i$ ($\hat{\mathbf m}_i$ the deformed unit direction) and $\delta J/J = \mathrm{div}_x\,\delta\mathbf u$:
        $$
        \delta\bar\lambda_i = \bar\lambda_i\left(\hat{\mathbf m}_i\cdot\nabla_x\delta\mathbf u\,\hat{\mathbf m}_i - \tfrac13\,\mathrm{div}_x\,\delta\mathbf u\right)
        = \bar\lambda_i\left(\hat{\mathbf m}_i\otimes\hat{\mathbf m}_i - \tfrac13\,\mathbf I\right) : \nabla_x\delta\mathbf u .
        $$
        Summing over the directions:
        $$
        \sum_i \rho\,a_i\,S_i\,\delta\bar\lambda_i = \mathbf Q : \nabla_x\delta\mathbf u,
        \qquad
        \mathbf Q = \frac{\partial c_{b,\mathrm{tot}}^{\tau}}{\partial(\nabla_x\mathbf u)} = \sum_i \rho\,a_i\,S_i\,\bar\lambda_i\left(\hat{\mathbf m}_i\otimes\hat{\mathbf m}_i - \tfrac13\,\mathbf I\right).
        $$
        $\mathbf Q$ is the **binding sensitivity tensor**: the change of the bound content per unit displacement gradient, in content units. It is not a stress. It is symmetric and deviatoric (a pure volume change does not change the isochoric stretches). Each direction contributes in proportion to its stretch sensitivity $S_i$, projected onto the gradient components that stretch it. The code already accumulates it, with the opposite sign, as `dcbdC` (part 1c).
      - *Nodal form.* With $\delta\mathbf u = N_b\,\delta\mathbf u_b$, the component $(\nabla_x\delta\mathbf u)_{kl} = (\partial N_b/\partial x_l)\,\delta u_{bk}$, so
      $$
      \frac{\partial c_{b,\mathrm{tot}}^{\tau}}{\partial u_{bk}}
      = G\,\frac{\partial\theta_f}{\partial \bar J}\,\frac{\partial \bar J}{\partial u_{bk}} + \sum_l Q_{kl}\,\frac{\partial N_b}{\partial x_l} .
      $$
      The first term enters 1a (reaction through $\theta_f$); the second is 1c (reaction through the stretches).
   5. *Volume changes.* With $\delta\mathbf u = N_b\,\delta\mathbf u_b$ and $\delta J = J\,\mathbf F^{-T} : \delta\mathbf F = J\,\mathrm{tr}(\delta\mathbf F\,\mathbf F^{-1}) = J\,\mathrm{div}_x\,\delta\mathbf u$, evaluated at the integration point for $J$ and at the centre for $\bar J$:
      $$
      \frac{\partial J}{\partial u_{bk}} = J\,\frac{\partial N_b}{\partial x_k},
      \qquad
      \frac{\partial \bar J}{\partial u_{bk}} = \bar J\,\left.\frac{\partial N_b}{\partial x_k}\right|_{c},
      \qquad
      \frac{\partial}{\partial u_{bk}}\!\left(\frac{J}{\bar J}\right) = \frac{J}{\bar J}\left(\frac{\partial N_b}{\partial x_k} - \left.\frac{\partial N_b}{\partial x_k}\right|_{c}\right).
      $$
   6. *Result.* Substituting steps 3–5 into step 2, and going back to the current volume ($dV_0 = dv/J$):
      $$
      K^{\mathrm{src}}_{cu,\,a(bk)}
      = \int_v \frac{N_a}{\bar J}\left[ \frac{c_{f,\max} + G}{\Delta t}\,\frac{\partial\theta_f}{\partial \bar J}\,\bar J\left.\frac{\partial N_b}{\partial x_k}\right|_{c}
      + \frac{1}{\Delta t}\sum_l Q_{kl}\,\frac{\partial N_b}{\partial x_l}
      + s\left(\frac{\partial N_b}{\partial x_k} - \left.\frac{\partial N_b}{\partial x_k}\right|_{c}\right) \right] dv .
      $$
      With $s/\bar J = -\mathrm{ResFac}$, the three terms are $K^{(1a)}$, $K^{(1c)}$ and $K^{(1b)}$. In the code, the coefficient of $\left.\partial N_b/\partial x_k\right|_c$ in $K^{(1a)}$ is `DSRCDJ` $= (c_{f,\max} + G)\,(\partial\theta_f/\partial \bar J)/\Delta t$, integrated with `detMapJC`·$w$ (C2, C3).

*(Item 2 is no longer used: the stretch term is now part 1c of item 1. Items 3–6 keep their numbers.)*

3. **$K_{cc}$, binding through $\theta_f$** (term 2, Reaction, of $\mathbf K_{\mu_f\mu_f}$ in `implementation_plan.md`, which was set to zero there).
   - **Origin.** The reaction part of the chemical residual is
     $$
     R^{(III)}_a = -\int_{\Omega^e} N_a\,\frac{\mathcal R}{\bar J}\,dv,
     \qquad
     \mathcal R = \frac{c_{b,\mathrm{tot}}^{\tau} - c_{b,\mathrm{tot}}^{t}}{\Delta t}.
     $$
     At fixed displacements, $\mu$ reaches $c_{b,\mathrm{tot}}^{\tau}$ only through $\theta_f$, via $k_{on}(\theta_f)$ in the kinetics. $J^c$ uses the lagged $c_b^t$, so there is no feedback. With $\partial\mu/\partial\mu_b = N_b$:
     $$
     \frac{\partial \mathcal R}{\partial \mu_b} = \frac{G}{\Delta t}\,\frac{\partial\theta_f}{\partial\mu}\,N_b .
     $$
     $\bar J$ and $dv$ do not depend on $\mu$.
   - **Term to add** ($K = -\partial R/\partial\mu$):
     $$
     \Delta K_{cc,\,ab} = \int_{\Omega^e} N_a\,\frac{G}{\Delta t\,\bar J}\,\frac{\partial\theta_f}{\partial\mu}\,N_b\,dv .
     $$
     This has the same structure as the existing capacity term (term 1), $\int N_a\,\frac{c_{f,\max}}{\Delta t\,\bar J}\frac{\partial\theta_f}{\partial\mu}N_b\,dv$. The two combine into a single factor, so in the code only `TanFac` changes, from
     $$
     \frac{c_{f,\max}}{\Delta t\,\bar J}\,\frac{\partial\theta_f}{\partial\mu}
     \quad\text{to}\quad
     \frac{c_{f,\max} + G}{\Delta t\,\bar J}\,\frac{\partial\theta_f}{\partial\mu}
     = \frac{\texttt{DSRCDMU}}{\texttt{detF}},
     $$
     with $\bar J$ = `detF`, the same F-bar determinant as in `ResFac`. The integration (`detMapJC`·$w$, i.e. the local $dv$), the matrix $N_aN_b$ and the diffusion terms of $K_{cc}$ stay as they are.
   - F-bar does not add terms here: $\mu$ enters only through its value at the integration point, and $\bar J$ does not depend on $\mu$.
   - $G = \partial c_{b,\mathrm{tot}}/\partial\theta_f = \sum_i \rho\,a_i\,S^{\theta}_i$ is the same $G$ as in part 1a; it was implemented with item 1 (C1, `dcbdthetaf`), so item 3 adds nothing to the kinetics. $S^{\theta}_i$ comes from a second sensitivity recursion in the existing sub-step loop, at no extra kinetics solves. Its input $\partial H_k/\partial\theta_f$ is not returned by `kineticsFunc` (only $\partial H/\partial c_b$ and $\partial H/\partial\bar\lambda_i$ are), so it is evaluated inline next to the `sens_i` update:
     $$
     S^{\theta}_k = \frac{S^{\theta}_{k-1} - \partial H_k/\partial\theta_f}{\partial H_k/\partial c_b},
     \qquad S^{\theta}_0 = 0,
     $$
     $$
     \frac{\partial H_k}{\partial\theta_f} = -\Delta t_{\mathrm{sub}}\,c_{f,\max}\,k_{on}\left[\frac{1}{(1-\theta_f)^2} - \frac{2\chi\,\theta_f}{1-\theta_f}\right],
     \qquad k_{on} = k_0\,e^{\chi(1-2\theta_f)}
     $$
     ($k_{off}$ does not depend on $\theta_f$).
   - The FD error of $K_{cc}$ ($2.5\times10^{-3}$) matches
     $$
     \frac{G}{c_{f,\max}} \approx \frac{\Delta t\,k_{on}}{(1-\theta_f)^2\,\partial H/\partial c_b}
     \quad\text{for } \Delta t = 0.05\ \mathrm{s}.
     $$
4. **$K_{uc}$, network stress through $\theta_f$** (optional, $2.5\times10^{-4}$): `DSIGDMU` lacks
   $$
   \sum_i \frac{\partial\boldsymbol\sigma_{\mathrm{net}}}{\partial c_{b,i}}\,S^{\theta}_i\,\frac{\partial\theta_f}{\partial\mu},
   $$
   the change of the filament lengths with the bound content. Add it with an accumulator in the direction loop (the same structure as `dPK2ficdcb`, weighted by $S^{\theta}_i$) only if the FD check shows it matters at large strain.
   - **How the stress is built.** $\boldsymbol\sigma = \boldsymbol\sigma_{\mathrm{vol}} + \boldsymbol\sigma_{\mathrm{iso}}$, with
     $$
     \boldsymbol\sigma_{\mathrm{vol}} = \frac{K\ln(\bar J/J^c)}{\bar J}\,\mathbf I,
     \qquad
     \boldsymbol\sigma_{\mathrm{iso}} = \mathbb P : \bar{\boldsymbol\sigma},
     \qquad
     \bar{\boldsymbol\sigma}_{\mathrm{net}} = \sum_i \frac{n_i}{\bar J}\,\rho\,a_i\,\frac{w'_i}{\bar\lambda_i}\,\mathbf m_i\otimes\mathbf m_i
     $$
     (`svol` from `vol`/`sigvol`; `siso = proje : sfic` in `sigiso`; $\bar{\boldsymbol\sigma}_{\mathrm{net}}$ = `snetficaf`, accumulated as `aux*sfilfic` with `aux` $= n_i/\bar J$ and $\mathbf m_i$ = `mfi`). $\theta_f$ enters $\boldsymbol\sigma_{\mathrm{vol}}$ directly through $J^c$, and $\bar{\boldsymbol\sigma}_{\mathrm{net}}$ only through the bound contents $c_{b,i}^{\tau}$. Each direction depends on its own $c_{b,i}$ through $n_i \propto c_{b,i}^{2/5}$ and $w'_i = \lambda_0\,r_0\,f_i$, at fixed $\bar\lambda_i$.
   - **Total derivative** (at fixed $\mathbf u$, so $\bar J$, $\bar\lambda_i$, $\mathbf m_i$ and $\mathbb P$ are fixed):
     $$
     \frac{\partial\boldsymbol\sigma}{\partial\mu}
     = \Bigg[\;
     \underbrace{-\frac{K\,V_{mol}\,c_{f,\max}}{\bar J\,J^c}\,\mathbf I}_{\text{existing: }\partial\boldsymbol\sigma_{\mathrm{vol}}/\partial\theta_f}
     \;+\;
     \underbrace{\mathbb P : \sum_i \frac{\partial\bar{\boldsymbol\sigma}_{\mathrm{net},i}}{\partial c_{b,i}}\,S^{\theta}_i}_{\text{new: }\partial\boldsymbol\sigma_{\mathrm{iso}}/\partial\theta_f}
     \;\Bigg]\frac{\partial\theta_f}{\partial\mu},
     \qquad
     \frac{\partial\bar{\boldsymbol\sigma}_{\mathrm{net},i}}{\partial c_{b,i}}
     = \frac{n_i}{\bar J}\,\rho\,a_i\,\frac{1}{\bar\lambda_i}\left(\frac{2\,w'_i}{5\,c_{b,i}} + \frac{\partial w'_i}{\partial c_{b,i}}\right)\mathbf m_i\otimes\mathbf m_i .
     $$
     The bracket $\left(\tfrac{2w'_i}{5c_{b,i}} + \partial w'_i/\partial c_{b,i}\right)$ is the existing `auxdwdcb` (the factor $2/5$ comes from $n_i$), as used in the $K_{uu}$ term 3 (`cfilficchem`). Both parts are multiplied by the same $\partial\theta_f/\partial\mu$ (`DTHETAFDMU`), so the new term is simply added to `DSIGDMU` before the multiplication. It is symmetric and deviatoric, and it makes `DSIGDMU` anisotropic (non-zero shear components).
   - **Code.**
     1. `affclnetfic_discrete`: new output `dsficdthetaf(3,3)` ($= \sum_i \partial\bar{\boldsymbol\sigma}_{\mathrm{net},i}/\partial c_{b,i}\,S^{\theta}_i$), initialised to zero. In the loaded branch, after `auxdwdcb` is computed:
        ```fortran
        call sigfilfic(dsfilficdcb, rho, lambdai, auxdwdcb, mfi, ai, ndi)   ! deformed mfi (Cauchy), not mf0i (PK2)
        dsficdthetaf = dsficdthetaf + aux * sensthetaf_i * dsfilficdcb
        ```
        Slack directions carry no stress, so they add nothing (their $S^{\theta}_i$ still enters $G$).
     2. `_umat_.f90`: project and add to the existing isotropic term:
        ```fortran
        CALL contraction42(dsigisodthetaf, proje, dsficdthetaf, ndi)        ! P : (...), as in sigiso
        DSIGDMU(I1,J1) = ( -((K*VMOL*CFMAX)/(DET*Jc))*UNIT2(I1,J1) + dsigisodthetaf(I1,J1) ) * DTHETAFDMU
        ```
        (no `phinet` factor: `sfic = (1-phinet)*sisomatfic + snetfic` adds the network stress unweighted).
     3. `u3d8.f90`: no change. `AmatUC` already takes all six components of `DSIGDMU` in the `Smat`/`Bmat` order (11, 22, 33, 12, 23, 13).
   - **Implemented and verified.** FD error of $K_{uc}$ (large strain, non-homogeneous, `BB` $=10^{-9}$): $2.5\times10^{-4}$ before, $5\times10^{-6}$ after with $h_\mu = 10^{-3}$. The remaining error scales as $1/h_\mu$ (round-off, since $\|K_{uc}\| \approx 1.5\times10^{-5}$): $5\times10^{-8}$ with $h_\mu = 0.1$. Use $h_\mu \approx 0.1$ when checking $K_{uc}$.
   - FD error of $K_{uc}$ without it: $4\times10^{-5}$ (state a) and $2.5\times10^{-4}$ (large strain, non-homogeneous); it grows with the strain, as expected of a network term.
   - The "bound-swelling" term $-\frac{K V_{mol}}{\bar J J^c}\,G\,\frac{\partial\theta_f}{\partial\mu}\,\mathbf I$ of the exact expression in `implementation_plan.md` ($\mathbf K_{\mathbf u\mu_f}$) must **not** be added. $J^c$ uses the lagged $c_{b,\mathrm{tot}}^{t}$, so $\boldsymbol\sigma_{\mathrm{vol}}$ does not depend on $c_b^{\tau}$ (same reason as the disabled `cvolchem`). FD check: adding it raises the $K_{uc}$ error to $2.6\times10^{-3}$ ($= G/c_{f,\max}$).
5. **$K_{cu}$, exact diffusion term (II).** The residual contains $-\int_v m\,\nabla_x N_a \cdot \nabla_x\mu\,dv$. At fixed nodal $\mu$, the current-configuration quantities vary as
   $$
   \delta(dv) = \mathrm{div}\,\delta\mathbf u\,dv,\qquad
   \delta\!\left(\frac{\partial N_a}{\partial x_i}\right) = -\frac{\partial N_a}{\partial x_j}\frac{\partial\delta u_j}{\partial x_i},\qquad
   \delta\!\left(\frac{\partial \mu}{\partial x_i}\right) = -\frac{\partial \mu}{\partial x_j}\frac{\partial\delta u_j}{\partial x_i},\qquad
   \delta m = \bar J\,\frac{\partial m}{\partial \bar J}\,\left.\mathrm{div}\,\delta\mathbf u\right|_{c},
   $$
   with $\partial m/\partial \bar J$ = `DMDJ` (the mobility is evaluated with $\theta_f(\bar J)$, so it varies with the centre volume, while $dv$ and the gradients are local). With $\delta u_j = N_b\,\delta u_{bj}$ this gives
   $$
   K^{(5)}_{cu,\,a(bk)} = \int_v \left[
   \underbrace{m\left(\nabla_x N_a\cdot\nabla_x\mu\right)\frac{\partial N_b}{\partial x_k}}_{\text{dilation}}
   \;+\; \underbrace{\bar J\,\frac{\partial m}{\partial \bar J}\left(\nabla_x N_a\cdot\nabla_x\mu\right)\left.\frac{\partial N_b}{\partial x_k}\right|_{c}}_{\text{mobility shift}}
   \;-\; \underbrace{m\,\frac{\partial N_a}{\partial x_k}\left(\nabla_x\mu\cdot\nabla_x N_b\right)}_{\text{distortion of }\nabla N_a}
   \;-\; \underbrace{m\,\frac{\partial \mu}{\partial x_k}\left(\nabla_x N_a\cdot\nabla_x N_b\right)}_{\text{distortion of }\nabla\mu}
   \right] dv .
   $$
   Without F-bar, the dilation and mobility-shift terms combine into $(m + J\,\partial m/\partial J)(\nabla_x N_a\cdot\nabla_x\mu)\,\partial N_b/\partial x_k$.
   - The code (`SpCUMod` built from `dMUdX` and `SpCUModFac` = $(m + \bar J\,\partial m/\partial \bar J)\,\delta_{il}$) keeps only a term of the last form, with $m + \bar J\,\partial m/\partial \bar J$ in place of $m$; this is the projection described in `implementation_plan.md`. The dilation, mobility-shift and $\nabla N_a$-distortion terms are missing.
   - Replace the `SpCUMod`/`AmatCU` assembly with $K^{(5)}$ written directly with `dshC`, `dshC0`, `dMUdX`, `MFLUID` and `DMDJ`. Check its sign and the column ordering against the FD driver; the existing assembly may already encode part of it with a different sign convention.
   - All four terms are proportional to $\nabla_x\mu$. In the FD state of the "Finding" (a gradient of a few J/mol per µm) they were smaller than part 1a. They matter at a diffusion front (stage 2) and under the indenter (stage 3).
6. **$K_{cu}$, prescribed surface flux on a deforming face.** A `*Dflux` with label `U1`–`U6` on face $n$ adds to the residual (`u3d8.f90`, lines 725–793)
   $$
   R^{\mathrm{flux}}_{c,a} = -\int_{a_n} N_a\,q\,da,
   $$
   where $q$ = `ADLMAG` is the outward normal flux per **current** area: positive means crosslinkers leaving, the opposite of the usual Abaqus `*Dflux` and `*Cflux` convention. The code computes $da$ from `coordsC` and states "no change to the tangent matrix", but the current area depends on the displacements. With the surface divergence
   $$
   \delta(da) = \left(\mathrm{div}_x\,\delta\mathbf u - \mathbf n\cdot\nabla_x\delta\mathbf u\,\mathbf n\right) da
   $$
   ($\mathbf n$ the outward unit normal in the current configuration) and $\delta\mathbf u = N_b\,\delta\mathbf u_b$:
   $$
   K^{(6)}_{cu,\,a(bk)} = -\frac{\partial R^{\mathrm{flux}}_{c,a}}{\partial u_{bk}}
   = \int_{a_n} N_a\,q\left(\frac{\partial N_b}{\partial x_k} - n_k \sum_l n_l\,\frac{\partial N_b}{\partial x_l}\right) da .
   $$
   - *Face-parametric form, convenient in `computeSurf3D`.* With $da = \lvert \mathbf x_{,\alpha}\times\mathbf x_{,\beta}\rvert\,d\alpha\,d\beta$, where $(\alpha,\beta)$ are the two local coordinates of the face and $\mathbf x_{,\alpha} = \sum_b N_{b,\alpha}\,\mathbf x_b$:
     $$
     \frac{\partial(da)}{\partial u_{bk}} = \mathbf n\cdot\left( N_{b,\alpha}\,\mathbf e_k\times\mathbf x_{,\beta} + \mathbf x_{,\alpha}\times N_{b,\beta}\,\mathbf e_k \right) d\alpha\,d\beta ,
     $$
     so that $K^{(6)}_{cu,\,a(bk)} = \sum_{\text{face points}} w\,N_a\,q\,\partial(da)/\partial u_{bk}$. This form only needs the in-face derivatives already computed in `computeSurf3D` and avoids the 3D gradient at the face.
   - It needs the outward unit normal. `computeSurf3D` computes one but does not use it, and its $y$-component has the wrong sign: the cross product gives $(\mathbf x_{,\alpha}\times\mathbf x_{,\beta})_2 = z_{,\alpha}\,x_{,\beta} - x_{,\alpha}\,z_{,\beta}$, but the code computes `dXdXi*dZdEta - dXdEta*dZdXi`, which is the negative. The same applies to the other two faces. Fix it, normalise it, and check the outward orientation of each face.
   - If $q$ were defined per **reference** area instead (`coords` instead of `coordsC`), $R^{\mathrm{flux}}$ would not depend on $\mathbf u$ and item 6 would vanish. Choose and document the convention.
   - Only relevant when `*Dflux` is used. None of the current inputs do; the bath is a prescribed $\mu$.
   - **If item 6 is implemented, also update the equations in `paper/model/res_tang.tex`:** the surface-flux term of the chemical residual (sign: the code uses $-\int N^A\,\check j\,\mathrm ds$ with $\check j = \mathbf j\cdot\mathbf n$ outward, while the document currently writes $+\int N^A\,\check j\,\mathrm ds$), its quadrature form, and the new $K_{cu}$ surface term.
   - FD (driver with `NDLOAD = 1`, faces 2 and 4): without item 6 the $K_{cu}$ error is $0.23$ for $q = 10^{-4}$ and $0.82$ for $q = 10^{-3}$, so it is required as soon as `*Dflux` is used.
7. **$K_{uu}$, network stress through $\theta_f(\bar J)$** (term 4 of $\mathbf K_{\mathbf{uu}}$ in `implementation_plan.md`, neglected there; found by the validation pass). At fixed $\mu$, a change of $\bar J$ changes $\theta_f$, hence the binding $c_{b,i}^{\tau}$ and the network stress:
   $$
   \delta\boldsymbol\sigma_{\mathrm{iso}} = \frac{\partial\boldsymbol\sigma_{\mathrm{iso}}}{\partial\theta_f}\,\frac{\partial\theta_f}{\partial\bar J}\,\delta\bar J,
   \qquad
   \frac{\partial\boldsymbol\sigma_{\mathrm{iso}}}{\partial\theta_f} = \mathbb P : \sum_i \frac{\partial\bar{\boldsymbol\sigma}_{\mathrm{net},i}}{\partial c_{b,i}}\,S^{\theta}_i
   $$
   (the same tensor `dsigisodthetaf` as item 4). With $\delta\bar J = \bar J\,\nabla_{x}\!\cdot\delta\mathbf u|_c$, the spatial tangent gains
   $$
   \Delta\mathbb c_{ijkl} = \bar J\,\frac{\partial\theta_f}{\partial\bar J}\left(\frac{\partial\boldsymbol\sigma_{\mathrm{iso}}}{\partial\theta_f}\right)_{ij}\delta_{kl}
   \qquad(\texttt{ddsigdde(i,j,k,k) += DET*DTHETAFDJ*dsigisodthetaf(i,j)}).
   $$
   This is the push-forward of the Lagrangian (PK2) form of term 4 in `implementation_plan.md`, restricted to the isochoric part and written per direction:
   $$
   \mathbb C^{(7)} = 2\,\frac{\partial\mathbf S_{\mathrm{iso}}}{\partial\theta_f}\otimes\frac{\partial\theta_f}{\partial\mathbf C},
   \qquad
   \frac{\partial\theta_f}{\partial\mathbf C} = \frac{\bar J}{2}\,\frac{\partial\theta_f}{\partial\bar J}\,\mathbf C^{-1},
   $$
   and $\frac{1}{\bar J}\bar{\mathbf F}\bar{\mathbf F}\,\mathbb C^{(7)}\,\bar{\mathbf F}^T\bar{\mathbf F}^T$ with $\bar{\mathbf F}\,\frac{\partial\mathbf S_{\mathrm{iso}}}{\partial\theta_f}\,\bar{\mathbf F}^T = \bar J\,\frac{\partial\boldsymbol\sigma_{\mathrm{iso}}}{\partial\theta_f}$ and $\bar{\mathbf F}\mathbf C^{-1}\bar{\mathbf F}^T = \mathbf I$ gives $\Delta\mathbb c$ above. The code builds it directly in spatial form, since `dsigisodthetaf` is already the spatial derivative.
   Added to `ddsigdde`, it passes through the F-bar correction exactly: `Gmat`$^T$`Amat Gmat` gives the local divergence, and the `Qmat` term ($\tfrac13 A_{ijmm}$) replaces it by the centroid divergence, i.e. $\delta\bar J/\bar J$. The bound-swelling counterpart ($\partial\boldsymbol\sigma_{\mathrm{vol}}/\partial c_b$) is zero (lagged $J^c$), as in `cvolchem`.
   - FD: the missing term scales with $K\,V_{mol}$ and the kinetics (it vanishes with $V_{mol}\to0$ or kinetics off). Force-sensitive catch state with $\Delta t = 0.5$: $K_{uu}$ error $1.1\times10^{-7}$ relative at $K = 1000$ ($1.6\times10^{-4}$ absolute) and $1.05\times10^{-7}$ at $K = 1$; with the term, $5.5\times10^{-9}$ in all states.
   - **Implemented** in `_umat_.f90`, after `DSIGDMU` (`ddsigdde(I1,J1,K1,K1) += DET*DTHETAFDJ*dsigisodthetaf(I1,J1)`); no change in `u3d8.f90`.

## Code changes
- **C1. `affclnetfic_discrete_imp.f90`:**
  - new `INTENT(OUT)` argument `dcbdthetaf` ($G$);
  - in the sub-step loop, alongside `sens_i`, accumulate `sensthetaf_i = (sensthetaf_i - dHdthetaf)/dHdcb_i`, with
    `dHdthetaf = -args(12)*cfmax*kon*(1/(1-thetaf)**2 - 2*CHI*thetaf/(1-thetaf))`;
  - at the end, `dcbdthetaf = dcbdthetaf + rho*ai*sensthetaf_i` (the weights used for `cbtau_tot`);
  - `dcbdC` is unchanged ($\mathbf Q$ = `−dcbdC`).
- **C2. `_umat_.f90` (`MATERIAL`):** new outputs, set to zero when `DTIME ≤ 1e-12`:
  - `DSRCDMU = (CFMAX + G)*DTHETAFDMU/DTIME` (replaces the $c_{f,\max}$-only rate in `TanFac`);
  - `DSRCDJ = (CFMAX + G)*DTHETAFDJ/DTIME`, with `DTHETAFDJ = K*VMOL/(RGAS*THETA*Jc*DET*df)` (factor out of the existing `DMDJ` expression; inside `MATERIAL`, `DET` $= \bar J$ after `fslip`);
  - `DSRCDGU(3,3) = -dcbdc/DTIME`.

  No new globals (local variables only); `MATERIAL` keeps its other arguments.
- **C3. `u3d8.f90`:**
  - pass the three new outputs from `call material(...)`;
  - `TanFac = DSRCDMU/detF`;
  - after the existing $K_{cu}$ update, add the source part
    $$
    K_{cu}\big(a,\ n_{\mathrm{dim}}(b-1)+k\big) \mathrel{+}= \texttt{detMapJC}\cdot w\cdot N_a \left[ \texttt{DSRCDJ}\cdot\texttt{dshC0}(b,k) + \frac{1}{\texttt{detF}} \sum_l \texttt{DSRCDGU}_{kl}\cdot\texttt{dshC}(b,l) - \texttt{ResFac}\,\big(\texttt{dshC}(b,k) - \texttt{dshC0}(b,k)\big) \right],
    $$
    with `detF` $= \bar J$ (`u3d8.f90`, line 355). The three terms are parts 1a, 1c and 1b. Check the column ordering against the existing `Gmat`/`Kuc` convention before coding;
  - replace the `SpCUMod`/`AmatCU`/`Gmat` diffusion part of $K_{cu}$ with the exact item 5, assembled directly from `dshC`, `dshC0`, `dMUdX`, `MFLUID` and `DMDJ` (`MATERIAL` already returns `MFLUID`, `DMDJ` and `DMUDX`; `SPCUMODFAC` is then no longer needed). The mobility-shift term uses `dshC0`, all the others `dshC`;
  - in the surface-flux loop (lines 725–793), add $K^{(6)}$ to `Kcu` at each face integration point (item 6). `computeSurf3D` must also return the in-face derivatives and the corrected outward unit normal. The surface-integration order `nInttS` is hard-coded to 1 in `uel.f90` (line 137); use 4 (`xintSurf3D4pt`) for distorted faces.
- **C4. Same changes in `5_CL_DIFF`** (same element and UMAT structure) after verification in `6_CL_DIFF_CS`.
- **Clean-up (independent):** `dpk2voldcb` (`_umat_.f90:108`) is never assigned and `dpk2dcb` is never used. Remove both, or initialise them, to avoid NaN with `-init=snan`.

## Verification
1. **FD tangent driver:** keep the scratch driver as a repository tool, e.g. `6_CL_DIFF_CS/test_tangent_fd.f90`. It wraps `main.f90` and does central differences of RHS with `SVARS` reset, with $h = 10^{-7}$ for $\mathbf u$ and $10^{-3}$ for $\mu$, and reports block errors. Targets: $K_{cu}$ and $K_{cc}$ within $10^{-5}$–$10^{-4}$ relative, $K_{uu}$ unchanged. States:
   - (a) the rotated, slightly deformed state above;
   - (b) a large-strain catch state ($\gamma \approx 0.5$, forces of several pN, $n_{\mathrm{sub}} > 1$ in some directions), which tests $\mathbf Q$ and both sensitivity recursions;
   - (c) prestress `LAMBDA0` $= 1.2$;
   - (d) a large chemical-potential gradient (e.g. 100–1000 J/mol across the element), which tests item 5;
   - (e) a non-homogeneous deformation (e.g. one node displaced, so that $J$ differs between integration points and from $\bar J$), with the chemistry out of equilibrium ($\mathrm{ResFac} \neq 0$). This tests the F-bar term 1b and the use of `dshC0` in part 1a and item 5. States (a)–(d) are homogeneous, so $J = \bar J$ everywhere and these terms are not tested there;
   - (f) a prescribed surface flux on one face (`NDLOAD = 1`, `JDLTYP = n`, `ADLMAG = q` passed to the UEL by the driver), on a deformed element, which tests item 6.

   Check the parts one at a time: the $K_{cu}$ error should fall after 1a (and 1b in state (e)), again after 1c in state (b), and after item 5 in state (d).

   **Capacity-only step (kinetics off).** Before the full checks, run states (a) and (e) with `KOFF0 = KCATCH0 = 0`. Then $k_{on} = (k_{off,0}+k_{catch,0})\,K_{eq} = 0$, the bound content is frozen, and $\mathcal R = 0$, $G = 0$, $\mathbf Q = 0$. Only the capacity parts (I) of 1a and 1b, and the existing capacity part of `TanFac`, remain, so the check tests them in isolation. Any remaining $K_{cu}$ or $K_{cc}$ error then lies in the $\theta_f(\bar J)$ derivative, the centre gradients or the F-bar ratio. With the kinetics back on, any additional error belongs to the reaction parts (III): $G$ in 1a and item 3, part 1c, and the $\mathcal R$ part of 1b. The code needs no flag for this; the properties alone switch the reaction off.
2. **Regression:** the residual is unchanged, so the standalone driver histories (stress, $c_b$, $\theta_f$) must match the current code to round-off.
3. **Abaqus:**
   - rotation test with the default increments: the $\mu$ drift and the isotropic stress should drop by orders of magnitude;
   - a stage-2 and a stage-3 case: fewer iterations per increment and a quadratic residual decrease in the `.msg`. Results should be unchanged to solver tolerance (spot-check one case per stage).
4. **Mass conservation:** in the closed stage-1 and stage-3 runs, the drift of the total content (currently $10^{-5}$ to $10^{-4}$) should decrease.
5. **Surface-flux check (Abaqus, one element):** a closed cube, with no other chemical boundary conditions, and a uniform `*Dflux` of magnitude $q$ (label `U`$n$) on one face of area $A$:
   - the total content must change at the rate $-q\,A$ (it decreases for $q > 0$, the outward-flux convention). This confirms the sign, the face numbering and the magnitude;
   - repeat for two faces (e.g. `U1` and `U4`), with `nInttS` = 1 and 4, and with the face stretched beforehand. The rate must follow the **current** area if $q$ is per current area;
   - compare with an equivalent `*Cflux` distribution (each of the face's 4 nodes getting $-q\,A/4$, positive into the body) to confirm the opposite sign convention of the two keywords.

### FD validation results (items 1, 3, 4, 5, 7 implemented; 2026-10-08)
One element, central differences of the residual with `SVARS` reset; `BB` $= 10^{-9}$ unless stated (to remove the neglected $\partial\rho/\partial\mathbf F$ term).
- **States covered:** rotation + general deformation (strain scale 1, 5, 10, 20), node perturbation ($J \neq \bar J$), $\mu$ gradient (100 J/mol), $D = 0$, kinetics off, $V_{mol} \to 0$, prestress `LAMBDA0` $= 1.1$–$1.2$, first increment (KINC = 1), $\Delta t = 0.05$ and $0.5$ s, bulk modulus $K = 1$, 10, 1000 (to resolve the chemical parts of $K_{uu}$).
- **Force-sensitive kinetics:** the driver ramps the deformation as $s = t/T$, so at the checked increment (KINC = 3) only 15% of the nominal strain is applied, and the forces there were small (< 0.12 pN). At full strain with the study properties the forces are 0.06 pN (max isochoric stretch 1.18), 0.4 pN (1.38) and 1.9 pN (1.38 with `LAMBDA0` $= 1.1$); with $\Delta x/k_BT = 1.2$ pN$^{-1}$ the kinetics are moderately force-sensitive. These full-strain states were added (`tshear` = 0.15, also with $\Delta t = 0.5$ and $K = 1$): all blocks at noise level. The kinetics were also tested with `DX` $= 0.5$, `KCATCH0` $= 0.5$, `DXC` $= 1.0$ (catch peak at $f^\ast \approx 0.008$ pN); there $n_{sub}$ reached 2 and 4.
- **Results (item 7 included):** $K_{uu}$ $5.5\times10^{-9}$, $K_{cu}$ $\le 4\times10^{-8}$, $K_{cc}$ $\le 2\times10^{-9}$, $K_{uc}$ $5\times10^{-9}$ (with $h_\mu = 0.1$). All remaining differences scale as $1/h$ (round-off).
- **Solver tolerance:** in the force-sensitive states $K_{cu}$ and $K_{cc}$ showed $h$-independent differences of $\sim10^{-6}$ and $2\times10^{-8}$. They disappear when the kinetics root is polished by two extra Newton steps: the residual uses a root accurate to $|H| < 10^{-12}$ (`solveKinetics`), not the exact root. This is a property of the FD test, not of the tangent, and is too small to affect convergence.
- **Discontinuities:** with prestress, FD steps $h_u \ge 10^{-6}$ cross the slack threshold of some direction (error $\sim 1$); use $h_u \le 10^{-7}$.
- **Building blocks** (`fil`, `kineticsFunc`) checked in isolation: $\partial w'/\partial\lambda$, $\partial f/\partial c_b$, $\partial H/\partial\lambda$, $\partial H/\partial c_b$ match FD to $\le 3\times10^{-8}$ for forces up to 230 pN.
- **Not covered:** item 6 (not implemented; see its FD numbers), multi-element assembly and Abaqus runs.

## Mobility convention (reference vs. spatial concentration) — RESOLVED: spatial (2026-10-08)
**Implemented** in `_umat_.f90`: `MFLUID` and `DMDMU` divided by `det` ($\bar J$), `DMDJ` divided by `det` with the extra $-m/\bar J$ term; documented in `res_tang.tex` ($K_{\mu_f\mu_f}$ mobility definition, $K_{\mu_f u}$ Step 1, listings). FD check (polished kinetics root, `BB` $=10^{-9}$): all blocks at the noise floor ($K_{uu}$ 4–6e-9, $K_{cu}$ ≤ 4e-8, $K_{cc}$ ≤ 1.4e-9, $K_{uc}$ round-off of $h_\mu$) in the reference, non-homogeneous, full-strain/prestress and force-sensitive states. Without the $-m/\bar J$ term the $K_{cu}$ error is 5e-3 (1.05e-2 with the sign flipped), so the term is exercised. In the reference state the residual derivative $K_{cc}$ changes by 7e-5 relative ($\mathcal O(\bar J - 1)$). The discussion that led to the decision is kept below.
- **Issue.** The residual uses $J^{-1}(\dot c_f + \mathcal R)\,dv = (\dot c_f + \mathcal R)\,dV$, so $c_f$ and $c_b$ are contents **per unit reference volume**. The flux is **spatial** (per current area, current gradient): $\mathbf j = -m\,\nabla_x\mu$ (`dshC`, `dMUdX`). But the mobility uses the reference content:
  $$
  m = \frac{D}{RT}\, c_f\,(1-\theta_f) \qquad (\texttt{MFLUID} = D/(RT)\cdot c_f\cdot(1-\theta_f)).
  $$
  A spatial Fickian flux needs the current concentration $c_f/\bar J$, so the current $m$ is too large by a factor $\bar J$.
- **Literature.** Chester, Di Leo & Anand (2015), *A finite element implementation of a coupled diffusion-deformation theory for elastomeric gels*, IJSS 52, 1–18 (`literature/`):
  - Eq. (2.11) and the Remark after it: the flux is taken as spatial, $\mathbf j = -m\,\mathrm{grad}\,\mu$. Their earlier papers (Chester & Anand 2010, 2011; Chester 2012) used a referential flux $\mathbf j_R = -m\nabla\mu$; following Hong et al. (2008) and Duda et al. (2010) they consider the spatial form more appropriate for isotropic gels.
  - Eq. (4.1): $m = Dc/(R\vartheta)$ with $c = c_R/J$ the concentration **per unit deformed volume**, i.e. the mobility is spatial.
  - Their UEL had the $1/J$: the commented-out line in `_umat_.f90`, `MFLUID = (D*(ONE/PHI_TAU - ONE))/(DET*VMOL*RGAS*THETA)`, is $Dc_R/(J\Omega R\vartheta) = Dc/(R\vartheta)$ (with $\Omega c_R = 1/\phi - 1$). The division by `DET` was lost when the mobility was rewritten for the crosslinker model.

  | | Chester et al. 2015 | current code |
  |---|---|---|
  | flux | spatial, $-m\,\mathrm{grad}\,\mu$ | spatial, $-m\nabla_x\mu$ ✓ |
  | content variable | $c_R$ per reference volume | $c_f$ per reference volume |
  | mobility | $m = Dc/(R\vartheta)$, $c = c_R/J$ (spatial) | $m = D/(RT)\,c_f(1-\theta_f)$ (reference content) ✗ |
- **Consistent version** (only free crosslinkers diffuse; $\bar J$ for F-bar consistency):
  $$
  m = \frac{D}{RT}\, \frac{c_f}{\bar J}\,(1-\theta_f) .
  $$
- **Code changes if adopted** (`_umat_.f90`):
  - `MFLUID` $\to$ `MFLUID/det`;
  - `DMDMU` $\to$ `DMDMU/det`;
  - `DMDJ` $\to$ (current expression)`/det` $-$ `MFLUID/det` (new $-m/\bar J$ term from $\partial(1/\bar J)/\partial\bar J$).

  No structural change in `u3d8.f90`: items 5 and 1 already use `MFLUID`, `DmDmu` and `detF*DmDJ`. Rerun the FD check afterwards, and update the mobility definition (and $\partial m/\partial\mu_f$, $\partial m/\partial\bar J$) in `paper/model/res_tang.tex` (K_μμ and K_μu sections). State the chosen convention explicitly in the paper either way.
- **Magnitude.** $\mathcal O(\bar J - 1)$: small for the nearly isochoric network ($\bar J \approx 1$; $v_{mol}\Delta c$ small), unlike Chester's swelling gels ($J \sim 2$–$10$).

## Neglected
The angular density $\rho_i$ depends on the deformed direction through `bangle` (fixed `prefdir`). With `BB` $= 10^{-6}$ its derivative is negligible; revisit only for anisotropic networks.
FD check (after item 1): through $c_{b,\mathrm{tot}} = \sum_i c_{b,i}\,\rho_i\,a_i$ this is the only remaining $K_{cu}$ error, and it scales linearly with `BB`. With `BB` $= 10^{-6}$ it is $3\times10^{-6}$ relative ($8.6\times10^{-10}$ absolute, also with kinetics off, $V_{mol}\to 0$ and $D = 0$). With `BB` $= 10^{-9}$ it falls to $4\times10^{-8}$, the FD noise floor. With `BB` $= 10^{-3}$ it rises to $3\times10^{-3}$.

## Optional: fully implicit $\theta_f$–$c_b$ update (not planned)
To remove the one-increment lag of $c_b$ in $J^c$, $\theta_f^{\tau}$ and the bound contents $c_{b,i}^{\tau}$ would have to satisfy, simultaneously,
$$
H\big(\theta_f;\,\mu^{\tau}, \bar J^{\tau}, c_{b,\mathrm{tot}}^{\tau}\big) = 0,
\qquad
H_i\big(c_{b,i};\,\theta_f, \bar\lambda_i^{\tau}\big) = 0 \quad (i = 1,\dots,N_{\mathrm{dir}}),
\qquad
c_{b,\mathrm{tot}}^{\tau} = \sum_i \rho\,a_i\,c_{b,i}^{\tau} .
$$
- **Coupling strength.** The two parts are coupled only through the scalar $c_{b,\mathrm{tot}}$ in $J^c$. The contraction factor of a fixed-point iteration is
  $$
  a\,G,
  \qquad
  a = \frac{\partial\theta_f}{\partial c_{b,\mathrm{tot}}} = -\frac{\partial H/\partial c_{b,\mathrm{tot}}}{d_f} = -\frac{\texttt{dHdcb}}{\texttt{df}},
  $$
  with `dHdcb` already computed in `_umat_.f90` (line 318). For the current parameters $a \sim 2\times10^{-3}$ mM⁻¹ and $G \sim 6\times10^{-6}$ mM, so $aG \sim 10^{-8}$. One extra pass (solve $\theta_f$ again with $c_{b,\mathrm{tot}}^{\tau}$, then repeat the kinetics) converges to round-off.
- **Algorithm.** Iterate {solve $\theta_f$ with the latest $c_{b,\mathrm{tot}}$ → kinetics (sub-stepped Backward Euler) → update $c_{b,\mathrm{tot}}$} until $|\Delta c_{b,\mathrm{tot}}|$ is below a tolerance, then evaluate $J^c$, the pressure and the stresses once with the converged state.
- **Code impact.**
  - `affclnetfic_discrete` currently does the kinetics and the stress assembly in one direction loop, so it must be split into a kinetics-only pass and a final stress and tangent pass. Otherwise every fixed-point iteration also pays for the stress assembly.
  - Expected cost: about 1.5–2 times the current material cost (the kinetics are the expensive part).
- **Tangent impact.**
  - $\theta_f$ then also depends on $c_{b,\mathrm{tot}}^{\tau}$. The derivatives of item 1 (1a–1c) and item 3 gain the scalar factor $1/(1-aG) \approx 1$, plus the stretch-driven shift $a\,\mathbf Q$ of $\theta_f$.
  - The volumetric stress depends on $c_{b,\mathrm{tot}}^{\tau}$ through $J^c$, so the term $\partial\boldsymbol\sigma_{\mathrm{vol}}/\partial c_{b,\mathrm{tot}}\cdot\partial c_{b,\mathrm{tot}}/\partial\mathbf C$ (`cvolchem`, currently commented out in `_umat_.f90`) becomes necessary in $K_{uu}$.
- **Assessment.** It removes an error of about $10^{-4}$ Pa per increment, far below the solver tolerances and the effects studied, at roughly twice the material cost. Not worth it for the current parameters. Reconsider only if $V_{\mathrm{mol}}$ or the bound content were orders of magnitude larger, i.e. if $aG$ were no longer small.

# Part D — Transient network (stress-free rebinding)

## D.1 Motivation and target
- **Finding (2026-10-09).** In the LAOS loop of Unterberger et al. (2013, Fig. 4a; digitised in `validation_data/unterberger2013_fig4a_laos.csv`) the stress at zero strain is $\sigma_{xy}(\gamma=0) \approx \pm 0.5$ Pa, about 12% of the peak: the stress leads the strain ($\sigma(\gamma=0) \approx \pm G''\gamma_A$). The current model cannot produce this. Its stress depends only on the current state $(\mathbf F, \theta_f, \{c_{b,i}\})$, and at $\mathbf F = \mathbf I$ every filament has $\lambda_f = 1$ and carries no force, whatever its bound content. Hence $\sigma_{xy}(\gamma = 0) = 0$ exactly, and the loop is pinched at the origin.
- **Python test** (12% LAOS at 0.05 Hz, 3rd cycle, study kinetics): $|\sigma_{xy}(\gamma=0)|/\sigma_{\max}$ = 5e-16 for $\lambda_0 = 1$; 4e-6 for $\lambda_0 = 1.012$ (Unterberger's prestretch); 4e-5 for $\lambda_0 = 1.05$; 4e-3 for $\lambda_0 = 1.05$ with 10× force sensitivity. A prestretch converts the anisotropy of $c_{b,i}$ into a residual shear stress of the right sign, but 3–5 orders of magnitude too small (experiment: about 0.12).
- **Missing physics: stress-free rebinding.** Crosslinks that rebind while the network is deformed take the current configuration as their reference. When the strain returns to zero those bonds are deformed, so the network carries stress; at fixed strain the stress relaxes on the turnover time scale.
- **Targets:** Fig. 4a (loop opening) → Fig. 4b (amplitude sweep, prediction) → Fig. 5 (frequency sweep). Unterberger's single Maxwell element ($\tau_1 = 2$ s, $\theta_1 = 0.835$) is the phenomenological baseline that the kinetics-based relaxation should replace.
- **References:** Tanaka & Edwards (1992), *Macromolecules* (transient networks); Vernerey, Long & Brighenti (2017), *J. Mech. Phys. Solids* (statistically based continuum theory for transient networks).

## D.2 Model: mean-field natural stretch per direction
- **New internal variable** per direction: the natural (reference) stretch $\lambda_{n,i}$ ($\lambda_{n,i} = 1$ at $t = 0$), the population-averaged stretch at which the currently bound crosslinks of direction $i$ were formed.
- **Split of the reaction** $\mathcal R_i = B - U_i$:
  $$
  B = k_{on}(\theta_f)\, c_f^{max}\, \frac{\theta_f}{1-\theta_f} \quad\text{(binding; independent of } c_{b,i} \text{ and of force)}, \qquad
  U_i = k_{off}(f_i)\, c_b^{max}\, \frac{\theta_{b,i}}{1-\theta_{b,i}} \quad\text{(unbinding)} .
  $$
- **Balance of the natural stretch.** New bonds form at the current stretch; leaving bonds carry the population average:
  $$
  \frac{\mathrm d}{\mathrm dt}\big(c_{b,i}\,\lambda_{n,i}\big) = B\, \bar\lambda_i - U_i\, \lambda_{n,i}
  \quad\Longrightarrow\quad
  c_{b,i}\, \frac{\mathrm d\lambda_{n,i}}{\mathrm dt} = B\, \big(\bar\lambda_i - \lambda_{n,i}\big) .
  $$
  *Where it comes from.* The quantity balanced is not $\lambda_{n,i}$ itself but the product $\Lambda_i = c_{b,i}\,\lambda_{n,i}$: the **sum of the reference stretches of all bound crosslinkers** of direction $i$ (per unit volume). It plays the same role for the reference stretch as $c_{b,i}$ plays for the number of bonds, and $\lambda_{n,i} = \Lambda_i/c_{b,i}$ is its average per bond.
  1. *Populations.* Let $c_{b,i}^{s}(t)$ be the concentration of bound crosslinkers of direction $i$ **per unit reference stretch** $s$: $c_{b,i}^{s}\,\mathrm ds$ is the content of bonds of direction $i$ whose reference stretch lies in $(s, s+\mathrm ds)$, i.e. bonds formed when the stretch of that direction was $s$ (a histogram over $s$ whose area is $c_{b,i}$). Then
     $$
     c_{b,i} = \int c_{b,i}^{s}\, \mathrm ds , \qquad \Lambda_i = \int s\, c_{b,i}^{s}\, \mathrm ds , \qquad \lambda_{n,i} = \frac{\Lambda_i}{c_{b,i}} .
     $$
     Equivalently, with discrete populations (cohorts) $k$ of content $c_{b,i}^{k}$ and reference stretch $\lambda_{n,i}^{k}$: $c_{b,i} = \sum_k c_{b,i}^{k}$, $\Lambda_i = \sum_k c_{b,i}^{k}\,\lambda_{n,i}^{k}$, so $\lambda_{n,i}$ is the content-weighted average of the cohorts' reference stretches. The bonds present at $t = 0$ form the cohort $k = 0$ with $\lambda_{n,i}^{0} = 1$.
  2. *Kinetics of the populations.* Bonds form at rate $B$, all with the current stretch as reference ($s = \bar\lambda_i(t)$, stress-free rebinding); every bound crosslinker of direction $i$ unbinds with the same rate per bond, $U_i/c_{b,i}$ (the mean-field assumption: all bonds of the direction see the same force):
     $$
     \frac{\partial c_{b,i}^{s}}{\partial t} = B\, \delta\big(s - \bar\lambda_i(t)\big) - \frac{U_i}{c_{b,i}}\, c_{b,i}^{s}(t) .
     $$
     The Dirac delta places all new bonds at the single reference stretch $s = \bar\lambda_i(t)$: it vanishes for $s \neq \bar\lambda_i$ and has unit area, so $\int B\,\delta\,\mathrm ds = B$ (rate of new bonds) and $\int s\,B\,\delta\,\mathrm ds = B\,\bar\lambda_i$ (each brings reference $\bar\lambda_i$). In the cohort picture it means "each time step adds a cohort of content $\Delta t\,B$ with $\lambda_{n,i}^{k} = \bar\lambda_i$", while the unbinding term reduces every cohort by the same factor $(1 - \Delta t\,U_i/c_{b,i})$. Because $\bar\lambda_i(t)$ changes with the loading, the spike moves and the cohorts spread over the stretches visited, older ones fading exponentially. (Tanaka & Edwards and Vernerey et al. insert new chains with the stress-free Gaussian distribution $P_0(\mathbf r)$; the delta is its per-direction counterpart for a single affine stretch. Only the mean of the new bonds' reference enters $\Lambda_i$, so a spread around $\bar\lambda_i$ would not change the result.)
  3. *Zeroth moment* ($\int \cdot\, \mathrm ds$): $\dot c_{b,i} = B - U_i = \mathcal R_i$, the existing kinetics (unchanged).
  4. *First moment* ($\int s\, \cdot\, \mathrm ds$): the new bonds add $B\,\bar\lambda_i$ (each brings its reference $\bar\lambda_i$); the leaving bonds remove $(U_i/c_{b,i})\,\Lambda_i = U_i\,\lambda_{n,i}$ (they carry, on average, the population's reference):
     $$
     \dot\Lambda_i = \frac{\mathrm d}{\mathrm dt}\big(c_{b,i}\,\lambda_{n,i}\big) = B\,\bar\lambda_i - U_i\,\lambda_{n,i} .
     $$
  5. *Equation for the average.* With the product rule, $\dot\Lambda_i = c_{b,i}\,\dot\lambda_{n,i} + \lambda_{n,i}\,\dot c_{b,i} = c_{b,i}\,\dot\lambda_{n,i} + \lambda_{n,i}\,(B - U_i)$, so the unbinding terms cancel:
     $$
     c_{b,i}\,\dot\lambda_{n,i} = B\,\bar\lambda_i - U_i\,\lambda_{n,i} - \lambda_{n,i}\,(B - U_i) = B\,\big(\bar\lambda_i - \lambda_{n,i}\big) .
     $$
  Interpretation: unbinding removes bonds but does not shift the average reference (bonds leave in proportion to the population); binding adds bonds with reference $\bar\lambda_i$ and pulls the average towards it, with a weight $B/c_{b,i}$ (each new bond counts more when fewer bonds remain). It is a mixing rule, like the temperature of a tank where water at temperature $\bar\lambda$ flows in at rate $B$ and the mixture flows out at rate $U$. The full distribution $c_{b,i}^{s}(t)$ (or the individual cohorts $c_{b,i}^{k}$, $\lambda_{n,i}^{k}$) never needs to be stored: only its first two moments ($c_{b,i}$, $\Lambda_i$) enter, which is the closure used by Tanaka & Edwards (distribution function) and Vernerey et al. (2017, distribution tensor $\boldsymbol\mu$, their Eq. 37). The closure is exact for the reference-stretch average; the stress, being nonlinear in $\lambda_{rel}$, is evaluated with this average (mean-field, see D.7).

  The relaxation rate is $B/c_{b,i}$, equal to $k_{off}/(1-\theta_b)$ at equilibrium: $\tau \approx (1-\theta_b)/k_{off}$, with **no new time parameter**. For HMM ($k_{off} \approx 0.3$ s⁻¹) this is a few seconds, consistent with Unterberger's $\tau_1 = 2$ s.
- **Filament stretch measure.** The relative stretch $\lambda_{rel,i} = \bar\lambda_i/\lambda_{n,i}$ replaces $\bar\lambda_i$ in the filament model: $\lambda_f$, the force, $w'$, $w''$, the slack test $\lambda_f(\lambda_{rel})\,\lambda_{0f} \le 1$, and the force in the kinetics. The orientation stays affine ($\bar{\mathbf m} = \bar{\mathbf F}\mathbf m_0$); only the stretch reference relaxes.
- **Stress.** With $w(\lambda_{rel})$, $\partial w/\partial\bar\lambda = w'(\lambda_{rel})/\lambda_n$ and $\partial^2 w/\partial\bar\lambda^2 = w''(\lambda_{rel})/\lambda_n^2$ enter $\bar{\boldsymbol\sigma}^{net}$ and $\bar{\mathbb c}^{net}$ in place of $w'$ and $w''$.
- **Partial relaxation.** Unterberger's data relax only partly ($\Delta G/G_\infty \approx 0.45$; G′ nearly constant over the frequency sweep). A transient fraction $\xi_t \in [0,1]$ of the network stress uses $\lambda_{rel}$; the permanent part $(1-\xi_t)$ uses $\bar\lambda$, with the same $c_{b,i}$. $\xi_t = 0$ recovers the current model exactly (regression). $\xi_t$ plays the role of Unterberger's $\theta_1$. (Alternative, not recommended first: let the isotropic matrix carry the permanent part.)
- **Expected behaviour:** at fixed strain the transient stress relaxes to zero at rate $B/c_b$; after shear history $\sigma_{xy}(\gamma = 0) \neq 0$; in the small-strain limit the response is close to one Maxwell element with modulus $\xi_t G$ and time $c_b/B$.

## D.3 Time integration (Backward Euler, explicit update)
- Backward Euler of the balance, combined with the kinetics update $c_{b,i}^\tau = c_{b,i}^t + \Delta t\,(B - U_i^\tau)$, gives
  $$
  \lambda_{n,i}^\tau = \frac{c_{b,i}^t\, \lambda_{n,i}^t + \Delta t\, B^\tau\, \bar\lambda_i^\tau}{c_{b,i}^t + \Delta t\, B^\tau} .
  $$
  It is a convex combination of $\lambda_{n,i}^t$ and $\bar\lambda_i^\tau$ (bounded), and it does not depend on $c_{b,i}^\tau$ because $B$ does not depend on $c_{b,i}$: no extra solve is needed.
- **Order per increment:** $\theta_f^\tau$ (as now) → $\lambda_{n,i}^\tau$ (explicit) → kinetics with $\lambda_{rel} = \bar\lambda/\lambda_n^\tau$ (fully implicit, existing solver) → stress.
- **Sub-stepping:** the same formula per sub-step, with $c_{k-1}$ and $\Delta t_{sub}$ ($\bar\lambda$ and $\theta_f$ frozen at $\tau$, as now).

## D.4 Python implementation (first; Gauss point)
- `material.py`:
  - new history array `lamn_history` (initialised to ones), committed with the other histories;
  - $\lambda_n^\tau$ computed in the kinetics sub-step loop of `evaluate_state` (`kinetics_update`);
  - $\lambda_{rel}$ passed to `_solve_local_kinetics` / `_kinetics_residual` (force, `DfDcb`, slack test);
  - `_assemble_stress` with $w'/\lambda_n$ and the $\xi_t$ split.
- New parameter `xi_t` (default 0.0) in `mech_props`, read with `.get`.
- `main.py`: the drivers that save/restore the history around trial evaluations (`solve_uniaxial`) must also save/restore `lamn_history`.

## D.5 Fortran port (after the Python validation)
- `affclnetfic_discrete_imp.f90`:
  - $\lambda_{n,i}$ per direction in the state variables: the `NSDV` value in `global.f90` grows by `NDIR` (106 → 196; changing the value only, no new global variable);
  - the same explicit update; $\lambda_{rel}$ in the `kineticsFunc` arguments and in `fil`; $w'/\lambda_n$ and $w''/\lambda_n^2$ in `sigfilfic`/`csfilfic`; the $\xi_t$ split.
- Input: `uel_variables` 848 → 1568 and `*User output variables` 106 → 196 in the `test_in_abaqus` builders; a new property (e.g. `PROPS(27)` = `XI_T`, `NPROPS` 26 → 27).
- **Tangents** (Part C machinery):
  - $\partial\lambda_n^\tau/\partial\bar\lambda = \Delta t\,B/(c_b^t + \Delta t\,B)$, so $\mathrm d\lambda_{rel}/\mathrm d\bar\lambda = \big(1 - \lambda_{rel}\,\partial\lambda_n/\partial\bar\lambda\big)/\lambda_n$, which enters Term 1 ($\bar{\mathbb c}^{net}$) and $S_i$ (Term 3, and $\mathbf Q$ in $K_{\mu u}$ 1c);
  - $\partial\lambda_n/\partial\theta_f$ through $B(\theta_f)$ adds contributions to the G-type terms (Term 4, item 4 $K_{u\mu}$, item 1a, item 3).
  - Verify with the FD kit and the sign-flip test, as for Part C.

## D.6 Verification
1. **Regression:** $\xi_t = 0$ is bit-identical to the current code (Python and Fortran).
2. **Relaxation test** (step shear, small strain): $\sigma(t)/\sigma(0^+) = 1 - \xi_t\big(1 - e^{-t/\tau}\big)$ with $\tau = c_b/B$; compare with the analytic single-Maxwell response. Check positive dissipation (loop area > 0).
3. **LAOS** (Python, closed system, $\gamma_A \approx 0.126$, 0.05 Hz, 3 cycles): $\sigma_{xy}(\gamma = 0) \neq 0$ with the experimental sign; compare with `validation_data/unterberger2013_fig4a_laos.csv`.
4. **Fit sequence:**
   - elastic recalibration first: the model's peak stress is currently about 150–300× too low; fit the linear plateau of Fig. 4b with $c_b$ frozen;
   - fit Fig. 4a (geometry, $\lambda_0$, $k_{off}$/$K_{eq}$, $\Delta x$, $\xi_t$);
   - predict Fig. 4b;
   - Fig. 5: with a single time scale, expect to need a relaxation spectrum later (e.g. CGD-type).
5. **Fortran:** FD tangent kit and sign-flip test; Python–Fortran parity on the LAOS case.

## D.7 Open choices (to confirm before coding)
- **Mean-field vs cohorts.** One natural stretch per direction (mean-field) first; several bound populations (cohorts) per direction only if the force nonlinearity makes the averaging inaccurate.
- **Permanent part:** the $\xi_t$ split is recommended.
- **Stretch-selective unbinding** (bonds with larger $\lambda_{rel}$ leave first) is not represented in the mean-field model. Note it as a limitation.

## D.8 Related work: Qian et al. (2026)
Qian, Tian, Liu, Lin & Fang (2026), *A visco-elasto-plastic constitutive theory of cytoskeletal deformation governed by dissipative crosslinkers and inter-filament slip*, Int. J. Eng. Sci. 228, 104640 (`literature/1-s2.0-S0020722526001783-main.pdf`).
- **Same physical assumption.** "Newly associated crosslinkers are considered stress-free": a crosslinker that rebinds at time $\tau$ takes the network configuration at $\tau$ as its reference and later sees only $\tilde{\mathbf F}_{t,\tau} = \mathbf F_t\mathbf F_\tau^{-1}$. This is the stress-free rebinding of D.2.
- **Same treatment of the history.** They sort crosslinkers into generations (number of rebinding events), each surviving with $e^{-k_{off}(t-\tau)}$, and show (their Eq. 13, Appendix B) that the sum over generations collapses into one convolution:
  $$
  \Delta\Psi_{clk}(t) = n_{sb}\, e^{-k_{off} t}\, \varphi_{sc}(\mathbf F_t) + n_{sb}\, k_{off} \int_0^t e^{-k_{off}(t-\tau)}\, \varphi_{sc}\big(\tilde{\mathbf F}_{t,\tau}\big)\, \mathrm d\tau .
  $$
  Our balance $\mathrm d(c_b\lambda_n)/\mathrm dt = B\bar\lambda - U\lambda_n$ is the differential form of the same exponential-kernel convolution, and the explicit update of D.3 is its recursive evaluation. In both cases the exponential kernel removes the need to store the history (answering the concern of a growing number of cohorts per direction and Gauss point).
- **Differences:**

  | | Qian et al. 2026 | Part D |
  |---|---|---|
  | carrier of the memory | crosslinker as a 3D neo-Hookean spring (tensor $\tilde{\mathbf F}_{t,\tau}$) | filament segment stretch per microsphere direction (scalar $\lambda_{n,i}$) |
  | kinetics | constant $k_{on}$, $k_{off}$; bound density at steady state ($n_{sb}$) | force-dependent catch–slip $k_{off}(f_i)$; $c_{b,i}$ varies in time and between directions |
  | exactness | exact: the neo-Hookean stress is linear in $\mathbf F_t\mathbf C_\tau^{-1}\mathbf F_t^T$, so the convolution reduces to one tensor-valued history variable (6 components) | mean-field: exact evolution of the average; the filament force is nonlinear, so stress of the average ≈ average of the stresses (exact at small strain) |
  | relaxation spectrum | several crosslinker types with different $k_{off}$ (Prony series / generalised Maxwell) | one time scale per direction (slightly spread by the force dependence) |
  | other mechanisms | Gaussian filament chains, inter-filament slip (Eyring-type plasticity), myosin activity | affine semiflexible filaments, diffusion of free crosslinkers, swelling, catch–slip |
- **What to take from it:**
  1. Independent support for the stress-free-rebinding assumption; without slip their model reduces to a finite-strain Maxwell-type viscoelastic response, as expected for Part D at small strain.
  2. Where the response is linear in the history, no mean-field approximation is needed; for us this holds only at small strain (nonlinear filament force near locking).
  3. A physically motivated route to a relaxation spectrum (Fig. 5, where Unterberger needed three Maxwell elements): several bound populations per direction with different $k_{off}$ (e.g. fast and slow), as an alternative to a CGD-type spectrum (D.6).
  4. Their irreversible deformation comes from inter-filament slip, absent in our model; not needed for Unterberger's in vitro network, but relevant for cell-level unloading tests.
