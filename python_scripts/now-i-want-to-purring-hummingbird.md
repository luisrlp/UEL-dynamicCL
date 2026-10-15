# Plan: staged optimisation of the LAOS fit (Unterberger 2013, Fig. 4a)

## Context
`python_scripts/optimise_laos.py` fits the stretch-free-binding model (`material.py`, `DynamicCrosslinking`) to the
3rd LAOS cycle (`validation_data/unterberger2013_fig4a_laos_dataset0.csv`) with `scipy.optimize.least_squares`
(TRF, normalised bounded parameters, parallel forward-difference Jacobian, cache with settings signature,
`setup.json`/`result.json`/`evals.csv`). Best known set: raw RMS 0.21 Pa.

Two issues were found that the current objective hides:
- the first cycle from rest is exaggerated (best set: cycle 1 σ_min = −5.3 Pa vs −3.7 Pa steady), because the
  reference λ_n is displaced during the first quarter cycle;
- the model's 3rd cycle is not yet steady (cycle 2→3 change 0.16 Pa RMS, cycle 4→5 0.007 Pa), while the paper
  states the measured 3rd cycle was steady.

Goal: a staged plan — Step 1 baseline fit of the 3rd cycle with a parameter add/remove procedure; Step 2a the same
fit with a penalty on the first-cycle overshoot; Step 2b a fit of the steady-state cycle. 2a and 2b are separate
fits started from the Step 1 optimum and compared with Step 1 on common metrics. All changes are in
`python_scripts/optimise_laos.py`; `material.py` is not modified. Each fit writes to its own `--out` folder
(`python_scripts/optimisation_runs/<step>_<label>/`).

## Common code changes (before Step 1)
1. **Full-history simulation.** `simulate(p, cfg)` returns γ(t), σ(t) of all cycles (not only the last); a new helper
   `split_cycles(t, g, s)` returns per-cycle arrays. `residuals_from_loop` is reused unchanged on the selected cycle.
2. **Per-cycle diagnostics in `evals.csv`.** For every run log, per cycle k: σ_max,k, σ_min,k, and the RMS change of
   cycle k vs k−1 (at the data strains, via `residuals_from_loop` with unit weights). Pure logging; the objective is
   unchanged. Since the column set changes, use new `--out` folders.
3. **Candidate parameters with default-fixed status.** Add `keq` (bounds [1, 20], x0 5), `beta` ([0.30, 0.70], x0 0.438),
   `mu0stretch` ([2.0e4, 6.0e4] pN, x0 38600) to `PARAMS`, and a `DEFAULT_FIXED = {'keq', 'beta', 'mu0stretch'}` set
   merged with `--fix` (new option `--free name` releases a default-fixed one). Map them in `material_props`:
   `k_on_0 = p['keq'] * koff0`, `BETA = p['beta']`, `MU0STRETCH = p['mu0stretch']`; remove `BETA`, `MU0STRETCH`
   from `FIXED_MECH` and the module constant `KEQ` (replace its uses in `main()` result output and `setup.json`).
   With these default-fixed at their old values, residuals are bit-identical to the current script.
4. **Start from the best known set**: run Step 1 with `--x0 lambda0=1.0467 --x0 koff0=0.08 --x0 dxc_nm=2 --x0 dxs_nm=0.5`
   (current `PARAMS` x0 values differ) or update `PARAMS` x0 accordingly.

## Step 1 — baseline: fit the 3rd cycle, add/remove parameters
1.0 **Resolution check** (`--evaluate-only` at the best set): (dt 0.1, factor 3) vs (dt 0.05, factor 3) vs
    (dt 0.1, factor 4). Keep the cheapest setting whose raw RMS differs < 0.01 Pa from the finest.
1.1 **Fit the core set**: λ₀, s_taut, koff0, φ, Δx_c, Δx_s (`--max-nfev 20`, 7 workers).
1.2 **Pruning rule** (applied after every fit, from `result.json`): fix a parameter at its fitted value when
    - it sits at a bound (`at_bound`) and the bound is physically acceptable; or
    - its standard error > 25 % of its bound range (not identified); or
    - |correlation| > 0.95 with another free parameter (fix the less physically constrained one); or
    - it changed < 1 % of its bound range between two consecutive fits (converged to a constant).
    Refit with the reduced set; record each decision in a short `notes.md` in the step folder.
1.3 **Add candidates one at a time**, starting from the pruned optimum: `keq`, then `beta`, later `mu0stretch`
    (`--free keq`, …). Keep a new parameter only if the raw RMS drops by > 5 % and it passes the pruning rule;
    otherwise fix it back to its default. Note: `keq` changes c_b0 and r_0f (elastic stiffness) and, with λ₀ > 1,
    the prestress-aware equilibrium in `_init_equilibrium`, so expect correlation with λ₀ / s_taut.
1.4 **Step 1 result** = final parameter set, `result.json`, diagnostics of cycles 1–5 at the optimum.

## Step 2a — penalise the first-cycle overshoot (start from rest, fit cycle 3)
- New option `--overshoot-weight w` (default 0 → identical to Step 1). Append 4 residual rows per early cycle
  k = 1, 2: `w_eff · max(0, |σ_max,k| − |σ_max,3|)` and the same for σ_min, so only overshoot beyond the cycle-3
  magnitude is penalised (the softer virgin loading of cycle 1 is physical and not penalised).
  `w_eff = w · sqrt(N_data / N_pen)` so that w = 1 gives the penalty the same RMS weight as the data rows.
- Include `overshoot_weight` in `cfg` so the cache signature and `setup.json` record it.
- Fits: w = 0.5, 1, 2 from the Step 1 optimum with the Step 1 free set.

## Step 2b — fit a steady-state cycle
- Use a fixed number of cycles N (not a stopping criterion, which would make the residual jump between Jacobian
  runs): choose N from the cycle diagnostics at the Step 1 optimum as the first cycle with change < 0.01 Pa RMS
  vs the previous (≈ 5 from the current best set; take N + 1 for margin). Run with `--cycles N`; the fit uses
  the last cycle (existing behaviour of `simulate`).
- Add a warning (and an `evals.csv` column) when the last-cycle change exceeds 0.01 Pa, i.e. the fitted cycle is
  not steady for that parameter set.
- Fit from the Step 1 optimum with the Step 1 free set. Cost ≈ N/3 × Step 1 per run.

## Comparison (Steps 1, 2a, 2b)
Table with: raw and weighted RMS vs data; σ_max, σ_min, σ(γ=0) of the fitted cycle; cycle-1 peak overshoot;
cycle 2→3 and last-cycle change; fitted parameters ± standard error; strongest correlations. Note the data's
asymmetric part (RMS 0.07 Pa) as an irreducible floor for a symmetric model.

## Verification
- With all new options at defaults (`--overshoot-weight 0`, candidates fixed at old values, `--cycles 3`) the
  residual vector at the best set is bit-identical to the current script (`--evaluate-only`, compare `cache.json`
  values / printed RMS 0.2111 raw).
- Smoke tests with cheap settings (`--dt 1 --factor 2 --cycles 1–3 --max-nfev 1`) for each option: penalty rows
  zero when w = 0 and positive when cycle 1 overshoots; `--free keq` changes `k_on_0` in `setup.json`/`result.json`.
- Step 2b: `--cycles 3` reproduces Step 1; `--cycles 6` at the best set gives raw RMS ≈ 0.19 Pa (cycles 4–5 value).
- Every fit folder contains `setup.json`, `setup_<time>.json`, `evals.csv` with per-cycle columns, `result.json`.
