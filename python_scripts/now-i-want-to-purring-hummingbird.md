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


--------------------------------------------------------------------------------

Proof [CATCH-SLIP, NO BELL]

From the previous step: if k_on does not depend on the mechanical state, consistency requires the ratio

$$\frac{k_{on}}{k_{off,i}} = K_0,e^{\chi(1-2\theta_f)},e^{-x_i}, \qquad x_i \equiv \frac{\mu^{mech}_{b,i}}{RT}$$

With k_on = k_on⁰e^{χ(1−2θ_f)}, that means

$$k_{off,i} = \frac{k_{on}^0}{K_0},e^{x_i} \qquad \text{for every attainable } x_i$$

Now take the seminar's catch–slip form (Eq. 5.36):

$$k_{off,i} = k_{catch}^0,e^{-\alpha_c x_i} + k_{slip}^0,e^{\alpha_s x_i}, \qquad \alpha_c, \alpha_s > 0$$

Consistency would need k_off,i e^{−x_i} to be constant:

$$k_{catch}^0,e^{-(\alpha_c+1)x_i} + k_{slip}^0,e^{(\alpha_s-1)x_i} = \text{const} \quad \forall x_i$$

Since α_c + 1 > 0, the first term varies with x_i unless k_catch⁰ = 0. The second term is constant only if α_s = 1. So the only consistent case is pure slip with α_s = 1, which is exactly Eq. (5.35). Any catch pathway (k_catch⁰ > 0) violates the condition. ∎

The same argument applies to the force-based law in the code, k_c0 e^{−Δx_c f/k_BT} + k_s0 e^{Δx_s f/k_BT}. It cannot equal const·e^{μ^mech/RT}, because it is non-monotonic in load while e^{μ^mech/RT} increases monotonically with stretch in our Ψ.

Consequence for your seminar: Eq. (5.36), the catch–slip extension with the stated constant k_on, is not thermodynamically consistent as written. Only Eq. (5.35) is.

How to make catch–slip consistent

Put the catch–slip shape into both rates, and keep their ratio set by the free energy. For any chosen unbinding law g(x_i), for example the two-pathway form in μ^mech or in f:

$$k_{off,i} = g(x_i), \qquad k_{on,i} = K_0,e^{\chi(1-2\theta_f)},e^{-x_i},g(x_i)$$

Then κ_i ≡ 1 by construction, so 𝒜_i w_i ≥ 0. Physically:
- The catch–slip shape sets how fast bonds turn over: slower in the catch range, faster in the slip range. That affects remodelling and dissipation.
- The free energy alone sets the occupancy of each direction, through e^{−x_i}.

This is the paired form with Γ_i = g(x_i)e^{−φx_i} and φ = 1, written out.

[SLIP ONLY, BELL]

The condition

With k_off = k_off⁰ e^{fΔx/k_BT} and a mechanics-independent k_on, consistency requires (from before)

$$k_{off,i} \propto e^{\mu^{mech}{b,i}/RT} \quad\Longleftrightarrow\quad \frac{\mu^{mech}{b,i}}{RT} = \frac{f_i,\Delta x}{k_BT} + \text{const}$$

So a slip-only Bell law with constant k_on is consistent if and only if the free energy's mechanical chemical potential per bond equals the Bell work, N_A f_i Δx per mole. With our current Ψ_mech, μ^mech_b = ∂Ψ_mech/∂c_b comes from the r_0f and n scaling. That is not of the form f·Δx with a free Δx, so Bell in force with a constant k_on is not consistent as the model stands.

Three consistent ways to use slip-only Bell

1. Bell in the mechanical chemical potential (the seminar's Eq. 5.35).

$$k_{off,i} = k_{off}^0,e^{\mu^{mech}{b,i}/RT}, \qquad k{on} = k_{on}^0,e^{\chi(1-2\theta_f)}$$

This is consistent with a constant k_on. It has Bell's exponential form, but the "force × distance" is replaced by μ^mech_b from Ψ. There is no free Δx: the exponent must be exactly μ^mech/RT. A tunable sensitivity is available only through the split:

$$k_{off,i} = k_{off}^0,e^{\phi x_i}, \qquad k_{on,i} = k_{on}^0,e^{\chi(1-2\theta_f)},e^{-(1-\phi)x_i}, \qquad 0 \le \phi \le 1$$

φ plays the role of Bell's x_b/(x_b + x_u). The sum of the two exponents is fixed by Ψ.

2. Bell in force, with Ψ extended to contain the Bell work. Add a free-energy contribution whose c_b-derivative gives N_A f_i Δx. Physically this is the work tied to the extension released when a crosslinker unbinds under filament tension: the compliant crosslinker's own energy, as part of the series compound. Then:
- k_off = k_off⁰ e^{fΔx/k_BT} with a constant k_on is consistent (φ = 1);
- or, with the split, k_off ∝ e^{f x_b/k_BT} and k_on ∝ e^{−f x_u/k_BT} with x_b + x_u = Δx.

The catch: the new term also enters the stress through ∂Ψ/∂C, so it changes the mechanics too. It has to be derived carefully, not just added.

3. Bell in force, Ψ unchanged, binding compensated. Keep k_off = k_off⁰ e^{fΔx/k_BT}, and set

$$k_{on,i} = k_{on}^0,e^{\chi(1-2\theta_f)},e^{,f_i\Delta x/k_BT, -, \mu^{mech}_{b,i}/RT}$$

This is consistent, but the binding rate then depends on the difference between the Bell term and the free energy's mechanical potential. That has no clear physical meaning. It is the same construction as for catch–slip, with g = Bell.

Recommendation

- The clean choice is Option 1: Bell written in μ^mech_b, optionally split by φ. It is already derived in your seminar framework. It needs no new energy term and leaves no free Δx. The only physical input is Ψ.
- Use Option 2 only if you want the classic force-based Bell with a measured Δx. That requires extending Ψ with the crosslinker work term, which changes the stress as well.


[CATCH-SLIP, BELL]

Detailed balance per pathway

Write the unbinding law as two Bell pathways with signed distances (x_p < 0 for the catch pathway):

$$k_{off} = \sum_{p=c,s} k_{p0},e^{f x_p/k_BT}, \qquad x_c = -\Delta x_c < 0,\quad x_s = \Delta x_s > 0$$

Microscopic reversibility means each pathway must also be travelled in the binding direction, with its own detailed balance. In Bell's tilted-landscape picture, the transition state of pathway p lies at x_p from the bound well, and the unbound well lies at Δ, the extension released when the bond breaks. Seen from the unbound side, the barrier of pathway p tilts by f(x_p − Δ). So:

$$k_{on} = K_0 \sum_{p} k_{p0},e^{f(x_p-\Delta)/k_BT} = K_0\Big[k_{c0},e^{-f(\Delta x_c+\Delta)/k_BT} + k_{s0},e^{f(\Delta x_s-\Delta)/k_BT}\Big]$$

Each pathway then satisfies k_on,p/k_off,p = K₀e^{−fΔ/k_BT}, and so does the total:

$$\frac{k_{on}}{k_{off}} = K_0,e^{-f\Delta/k_BT}$$

That is consistent provided the free energy contains the matching work per bond, μ^mech_b = N_A fΔ (Option 2 from before, extended to two pathways).

Two consistent versions for our model

A. Fully Bell. Extend Ψ with the crosslinker work term, so that μ^mech_b = N_A fΔ, and use both rates above. The parameters are k_c0, k_s0, Δx_c, Δx_s (as now) plus one more physical length, Δ. As noted before, the extra Ψ term also enters the stress and has to be derived carefully.

B. Bell catch–slip shape for turnover, Ψ unchanged for occupancy. Keep the current Ψ and use the catch–slip Bell law as the shared factor:
$$k_{off,i} = g(f_i), \qquad k_{on,i} = K_0,e^{\chi(1-2\theta_f)},e^{-\mu^{mech}_{b,i}/RT},g(f_i)$$

with g the current two-pathway law in force. This is consistent without changing Ψ. It is the shared-prefactor construction again, with Bell's catch–slip law as g.

What consistency implies physically

- Binding becomes force-dependent automatically. In version A, the catch pathway's binding rate drops as e^{−f(Δx_c + Δ)/k_BT}. The slip pathway's binding rate drops whenever Δx_s < Δ. So a consistent catch–slip Bell model contains force-suppressed binding. This is not an extra assumption but a consequence of microscopic reversibility.
- Occupancy falls under load. In version A, the equilibrium bound fraction drops as e^{−fΔ/k_BT}: the Wolff-type bond loss. That is the same kind of effect that worsened the fit in our probe (redistribution to the directions carrying G_M). Its strength now depends on Δ, which is a physical length rather than a free Δx_u.
- In version B, the occupancy follows μ^mech_b from our Ψ. How strongly that depends on load I haven't checked yet. It may be much weaker than e^{−fΔ}.

Summary

Version: Current code
Unbinding: Two-pathway Bell in force
Binding: Constant
Ψ: Unchanged
Consistent?: ✗
────────────────────────────────────────
Version: A
Unbinding: Two-pathway Bell in force
Binding: Same pathways, shifted by Δ
Ψ: + crosslinker work term
Consistent?: ✓
────────────────────────────────────────
Version: B
Unbinding: Two-pathway Bell in force (g)
Binding: K₀e^{χ}e^{−μ^mech/RT}·g
Ψ: Unchanged
Consistent?: ✓
## SFR dissipation with prestretch

### 1. Ingredients
For each microsphere direction $i$ (weight $\rho_i$), stretch-free rebinding (SFR) introduces a reference stretch $\lambda_{n,i}$, and the filament energy depends on the relative stretch
$$\lambda_{e,i} = \frac{\bar\lambda_i}{\lambda_{n,i}}, \qquad \bar\lambda_i = \sqrt{\mathbf M_i\cdot\bar{\mathbf C}\,\mathbf M_i}.$$
The mechanical free energy is (Ch. 4, with SFR)
$$\Psi_{mech} = \sum_i \rho_i\, n_i(c_{b,i})\, w(\lambda_{e,i}; c_{b,i}),$$
where $w$ is the compound-filament energy with prestretch $\lambda_0$ (Eqs. 4.11, 4.16). Its derivative with respect to the relative stretch is (Eq. 4.17)
$$w'(\lambda_e) \equiv \frac{\partial w}{\partial\lambda_e} = f\,\lambda_0\, r_0 \;\ge 0,$$
with the tension–compression switch: $f = 0$ (and $w' = 0$) when the filament is slack, $\lambda_f(\lambda_e)\,\lambda_{0,f} \le 1$. Define the slack threshold $\lambda_e^*$ by $\lambda_f(\lambda_e^*)\,\lambda_{0,f} = 1$:
- without prestretch ($\lambda_0 = 1$, so $\lambda_{0,f} = 1$): $\lambda_e^* = 1$;
- with prestretch ($\lambda_0 > 1$): $\lambda_e^* < 1$, and the reference state is taut, $w'(1) = f_0\lambda_0 r_0 > 0$.

The SFR update used in `material.py` is
$$\lambda_n^{t+\Delta t} = \frac{c_b\,\lambda_n^t + \Delta t\,B\,\bar\lambda}{c_b + \Delta t\,B}
\;\Rightarrow\;
\frac{\lambda_n^{t+\Delta t}-\lambda_n^t}{\Delta t} = \frac{B\,(\bar\lambda - \lambda_n^t)}{c_b + \Delta t B}
\;\xrightarrow{\;\Delta t\to 0\;}\;
\dot\lambda_{n,i} = \Gamma_i\,(\bar\lambda_i - \lambda_{n,i}), \qquad \Gamma_i = \frac{B}{c_{b,i}} > 0 .$$

### 2. The dissipation term
With $\Psi = \Psi(\mathbf C, \lambda_{n,i}, c_{b,i}, c_f)$, the isothermal dissipation inequality (Eq. 5.10) gives the state relation for the stress, $\mathbf S = 2\,\partial\Psi/\partial\mathbf C$, so the stress power cancels. The chemical terms ($c_{b,i}$, $c_f$) give the reaction dissipation $\sum_i\rho_i\mathcal A_i w_i$ treated before. The new internal variable leaves the term
$$\mathcal D_\lambda = -\sum_i \frac{\partial\Psi}{\partial\lambda_{n,i}}\,\dot\lambda_{n,i} \;\ge 0 \quad\text{(required)}.$$
With $\partial\lambda_{e,i}/\partial\lambda_{n,i} = -\bar\lambda_i/\lambda_{n,i}^2$:
$$\frac{\partial\Psi}{\partial\lambda_{n,i}} = \rho_i n_i\, w'(\lambda_{e,i})\,\frac{\partial\lambda_{e,i}}{\partial\lambda_{n,i}} = -\rho_i n_i\, w'(\lambda_{e,i})\,\frac{\bar\lambda_i}{\lambda_{n,i}^2}.$$
Inserting the SFR evolution law:
$$\mathcal D_{\lambda,i} = \rho_i n_i\, w'(\lambda_{e,i})\,\frac{\bar\lambda_i}{\lambda_{n,i}^2}\,\Gamma_i(\bar\lambda_i - \lambda_{n,i}).$$
Since $\dfrac{\bar\lambda}{\lambda_n^2}(\bar\lambda - \lambda_n) = \dfrac{\bar\lambda}{\lambda_n}\Big(\dfrac{\bar\lambda}{\lambda_n} - 1\Big) = \lambda_e(\lambda_e - 1)$:
$$\boxed{\;\mathcal D_{\lambda,i} = \rho_i\, n_i\,\Gamma_i\; w'(\lambda_{e,i})\;\lambda_{e,i}\,(\lambda_{e,i} - 1)\;}$$

### 3. Sign analysis
All prefactors are non-negative ($\rho_i, n_i, \Gamma_i, \lambda_{e,i} > 0$, $w' \ge 0$), so $\operatorname{sign}\mathcal D_{\lambda,i} = \operatorname{sign}\big[w'(\lambda_{e,i})\,(\lambda_{e,i}-1)\big]$:

| State of direction $i$ | $w'$ | $\lambda_e - 1$ | $\mathcal D_{\lambda,i}$ |
|---|---|---|---|
| stretched beyond its reference, $\lambda_e > 1$ | $> 0$ | $> 0$ | $> 0$ ✓ |
| at its reference, $\lambda_e = 1$ | $\ge 0$ | $0$ | $0$ ✓ |
| **taut but below its reference, $\lambda_e^* < \lambda_e < 1$** | $> 0$ | $< 0$ | **$< 0$ ✗** |
| slack, $\lambda_e \le \lambda_e^*$ | $0$ | $< 0$ | $0$ ✓ |

- **Without prestretch** ($\lambda_0 = 1$): $\lambda_e^* = 1$, the third row is empty, and $\mathcal D_\lambda \ge 0$ in every state. SFR is then a purely dissipative relaxation toward the force-free state. ✓
- **With prestretch** ($\lambda_0 > 1$): in the band $\lambda_e^* < \lambda_e < 1$ the update moves $\lambda_n$ down toward $\bar\lambda$. This raises $\lambda_e$ toward 1 and re-tensions the filament toward $f_0$, so the stored energy increases with no stress power supplied: $\mathcal D_\lambda < 0$. ✗
- **Size of the band** at the reference parameters (`plate_guo_dxs`: $\lambda_0 = 1.027$, $r_{0,f} = 1.694\,\mu$m, $f_0 = 0.185$ pN): $\lambda_e^* = 0.9737$, i.e. a band of about 2.6 % in relative stretch. Directions in mild compression during LAOS, and at the $\gamma = 0$ crossings, fall in it.
- **Most negative value:** linearising $w' \approx \lambda_0 r_0\,[f_0 + k(\lambda_e - 1)]$ with $k = df/d\lambda_e$ at the reference, and $u = 1 - \lambda_e \in (0, f_0/k)$, gives $\mathcal D_{\lambda,i} \approx -\rho_i n_i\Gamma_i\lambda_0 r_0\,(f_0 - k u)\,u$. Its minimum is at $u = f_0/(2k)$:
$$\mathcal D_{\lambda,i}^{min} \approx -\rho_i n_i\,\Gamma_i\,\lambda_0 r_0\,\frac{f_0^2}{4k}.$$
It vanishes as $f_0 \to 0$.

### 4. Why it cannot be fixed by the kinetics
The update is driven by the gross binding flux $B$, not the net reaction rate. Consider one exchange in a direction at steady turnover: an existing bond unbinds (at the mean state $\lambda_e$) and a new bond binds at the reference ($\lambda_e = 1$, force $f_0$).
- **Chemically, the exchange is neutral:** the binding free energy gained on binding equals the amount lost on unbinding, so its affinity is zero, whatever the rate laws.
- **Mechanically,** it changes the stored energy by $w(1) - w(\lambda_e)$ per bond. In the band, $w$ increases between $\lambda_e$ and 1, so this is positive.

A chemically neutral cycle that raises the stored free energy of a passive system violates the second law. No choice of $k_{on}$ and $k_{off}$ (paired or not) can fix it: it is the reference target $\lambda_e \to 1$ (force $f_0$) that is inconsistent, because with $f_0 > 0$ that target is not a minimum of $w$ ($w'(1) > 0$; the minimum is the slack state). Charging the cost to the binding affinity only works for **net** binding, e.g. transient network formation, as in Ma & Berro. It does not work for steady turnover.

### 5. Consistent alternatives
1. **No prestretch** ($\lambda_0 = 1$). Consistent, but rejected by the fits: $G'_{lin}$ collapses (3.6 vs 13.5 Pa) because the rest state then lies on the slack switch.
2. **An energy source.** Prestress regeneration needs active work (myosin, ATP). That is absent in rigor-HMM networks during the rheology; in Unterberger's protocol, motors act only before the rigor transition.
3. **Force-free rebinding with slowly relaxing prestress.** Let new bonds form at the force-free state, $\lambda_e^{target} = \lambda_e^*$, instead of at $f_0$. The evolution becomes $\dot\lambda_{n,i} = \Gamma_i\,(\bar\lambda_i/\lambda_e^* - \lambda_{n,i})$, and the same steps give
$$\mathcal D_{\lambda,i} = \rho_i n_i\Gamma_i\, w'(\lambda_{e,i})\,\frac{\lambda_{e,i}}{\lambda_e^*}\,(\lambda_{e,i} - \lambda_e^*) \;\ge 0 \quad\text{for every state and any } \lambda_0,$$
because $w' = 0$ wherever $\lambda_e \le \lambda_e^*$. This is consistent with prestress: the prestress built during network formation becomes stored energy that relaxes, instead of a state that is regenerated.
   - **With the current single-rate SFR,** the prestress would then decay within $1/\Gamma \approx 10$ s, and $G'_{lin}$ would collapse as in alternative 1.
   - **With crosslink-governed dynamics (CGD, Broedersz 2010),** a uniform tension along a filament is the longest-wavelength mode, so it relaxes only at the terminal time $\tau_r \approx \tau_{off}(L/\ell_c)^2$, i.e. tens of minutes. The prestress then persists over the rheology experiments (consistent with the slow internal-stress relaxation in Lieleg 2011), while local modes still relax at $\tau_{off}$ and give the dissipation.

   So the thermodynamically consistent form of SFR (force-free rebinding) is compatible with the measured linear modulus only if the relaxation is crosslink-governed. This is an independent argument for the CGD extension (Step 1).

### 6. Summary
- $\mathcal D_{\lambda,i} = \rho_i n_i\Gamma_i\,w'(\lambda_{e,i})\,\lambda_{e,i}(\lambda_{e,i}-1)$: non-negative for $\lambda_0 = 1$, negative in the taut-below-reference band ($\approx 2.6$ % wide at the reference set) for $\lambda_0 > 1$.
- The violation comes from the SFR target (regenerating $f_0$), not from the kinetics, so paired $k_{on}/k_{off}$ cannot remove it.
- Consistent options: no prestretch (rejected by the fits), active regeneration (not in rigor-HMM), or force-free rebinding with crosslink-governed, slowly relaxing prestress (the CGD route).

## CGD extension in our constitutive framework (proposal)

### 1. Where the split comes from
The split of the direction stretch is not a phenomenological choice of elements. It is the mode decomposition of Broedersz et al. (2010), Eq. (4), applied to the filament of each microsphere direction.

- **Broedersz Eq. (4):** between crosslink unbinding events, the transverse shape $\mathbf r_\perp(x)$ of a filament (arclength $x$) relaxes by diffusion along the filament,
$$\tau_{off}\,\partial_t\mathbf r_\perp = \frac{\ell_c^2}{2}\,\partial_x^2\mathbf r_\perp + \boldsymbol\xi_\perp .$$
- **Fourier modes:** writing $\mathbf r_\perp = \sum_q \mathbf u_q e^{iqx}$ (with $q$ in units of $1/\ell_c$) decouples it into modes, each relaxing at
$$\Gamma(q) = \frac{q^2}{2\,\tau_{off}}, \qquad q \in \Big[\frac{\pi\ell_c}{L},\ \pi\Big]$$
  The upper cut-off is the crosslink spacing $\ell_c$; the lower cut-off is the filament length $L$.
- **Length stored in each mode:** the end-to-end length of the filament, which is the stretch the network sees, changes when undulations of wavelength $2\pi\ell_c/q$ are pulled out or let in. Each mode therefore contributes a relative length change $\varepsilon_q$, and these add along the filament (series, compliances add; Broedersz: "contributions to the fluctuation spectra add").
- **Weights:** the mean-field compliance of the modes, Broedersz Eq. (5), is
$$\chi_{cgd}(\omega) \propto \frac{k_BT\ell_c^3}{\kappa^2}\int\frac{dq}{q^2 + 2i\omega\tau_{off}},$$
  i.e. each mode carries a compliance proportional to $\Delta q/q^2$.
- **Elastic part:** wiggles shorter than $\ell_c$ equilibrate in milliseconds. They are the entropic spring of a segment between crosslinks, which is our β-model $w(\lambda_e)$ (Eqs. 4.11, 4.16).
- **Terminal part:** relaxing the whole filament takes $\tau_r \approx \tau_{off}(L/\ell_c)^2$ (Broedersz, text).

Discretising the mode integral into $K$ bands gives the split below. $K$ is a numerical resolution parameter (convergence), not a material parameter.

### 2. Kinematics (per direction $i$)
The current SFR reference $\lambda_{n,i}$ is generalised into a product of modal and terminal stretches:
$$\bar\lambda_i = \lambda_{e,i}\;\lambda_{r,i}\prod_{k=1}^{K}\lambda_{k,i}, \qquad \varepsilon_{k,i} = \ln\lambda_{k,i},\quad \varepsilon_{r,i} = \ln\lambda_{r,i},$$
$$\lambda_{n,i} \equiv \lambda_{r,i}\prod_k\lambda_{k,i}, \qquad \lambda_{e,i} = \bar\lambda_i/\lambda_{n,i} \quad\text{(as now)} .$$

### 3. Free energy
$$\Psi_{mech} = \sum_i\rho_i\,n_i(c_{b,i})\Big[w(\lambda_{e,i}; c_{b,i}) + \sum_k\frac{\varepsilon_{k,i}^2}{2\,j_{k,i}}\Big]$$
- $w$: the β-model compound-filament energy with prestretch $\lambda_0$, unchanged.
- $j_{k,i}$: modal compliance (units 1/energy, per filament), from Broedersz Eq. (5):
$$j_{k,i} = R\,a_k\,j_{\infty,i}, \qquad a_k = \frac{\Delta q_k/q_k^2}{\sum_l \Delta q_l/q_l^2}, \qquad j_{\infty,i} = \Big[\frac{d(w'\lambda_e)}{d\ln\lambda_e}\Big]^{-1}_{\lambda_e = 1},$$
  where $j_\infty$ is the plateau compliance of the elastic part and $R$ is the ratio of total modal compliance to plateau compliance. A first estimate from Broedersz Eq. (5) with the segment modulus $\mu = 90\kappa^2/(k_BT\ell_c^3)$ is $R \approx 0.0036\cdot 90\,\big(L/(\pi\ell_c) - 1/\pi\big)$, about 1.1 for $L = 20\,\mu$m. The prefactors are uncertain, which Step 0 checks.
- The terminal mode stores no energy; it only relaxes.

### 4. Stress
At fixed internal variables, $\partial\lambda_{e,i}/\partial\bar\lambda_i = 1/\lambda_{n,i}$, so the stress keeps its current form (Eqs. 4.24–4.27 with SFR):
$$\boldsymbol\sigma^{affine} = J^{-1}\sum_i\rho_i\,n_i\,\frac{w'(\lambda_{e,i})}{\lambda_{n,i}}\,\frac{1}{\bar\lambda_i}\,\mathbf m_i\otimes\mathbf m_i ,$$
with $\lambda_{n,i} = \lambda_{r,i}\prod_k\lambda_{k,i}$. The code's `dw = f*l0*r0/lamn` is unchanged; only the update of `lamn` changes. The tangent's $\lambda_n$-path term (Term 1b) becomes a sum over the modes.

### 5. Thermodynamic forces and evolution (gradient flow)
With $\partial\lambda_e/\partial\varepsilon_k = \partial\lambda_e/\partial\varepsilon_r = -\lambda_e$, the forces conjugate to the internal variables are
$$X_{k,i} = -\frac{\partial\Psi}{\partial\varepsilon_{k,i}} = \rho_i n_i\Big[w'(\lambda_{e,i})\lambda_{e,i} - \frac{\varepsilon_{k,i}}{j_{k,i}}\Big], \qquad X_{r,i} = -\frac{\partial\Psi}{\partial\varepsilon_{r,i}} = \rho_i n_i\,w'(\lambda_{e,i})\lambda_{e,i} \ge 0 .$$
Evolution with non-negative mobilities set by the crosslink kinetics:
$$\dot\varepsilon_{k,i} = \Gamma_{k,i}\big[j_{k,i}\,w'(\lambda_{e,i})\lambda_{e,i} - \varepsilon_{k,i}\big], \qquad \Gamma_{k,i} = \Gamma_i\,\frac{q_k^2}{2},$$
$$\dot\varepsilon_{r,i} = \Gamma_{r,i}\,j_{r,i}\,w'(\lambda_{e,i})\lambda_{e,i}, \qquad \Gamma_{r,i} = \Gamma_i\,\frac{q_{min}^2}{2},\quad j_{r,i} = (1+R)\,j_{\infty,i},$$
with $\Gamma_i = B/c_{b,i}$, the current per-bond turnover rate (the SFR rate), which carries the catch–slip $k_{off}(f_i)$ and the bound content.

- **Modes $k$** relax toward their equilibrium share of the filament force, $\varepsilon_k \to j_k\,w'\lambda_e$. On their own they add compliance (creep) but never relax the force to zero.
- **The terminal mode** flows while the filament carries force and stops when it is slack ($w' = 0$). This is force-free rebinding of the whole filament, at the slow rate $\Gamma_r \approx 1/\tau_r$.
- **Microstructural mapping:** $\ell_c = r_{0,f}(c_{b,i})$ (Eq. 5.1, already in the model), $\tau_{off} = 1/\Gamma_i$, $L$ = filament length (the one new constant, from the literature). With the dimensionless $q$, $\Gamma(q) = q^2/(2\tau_{off})$, which is Broedersz's rate.
- **Initial state:** with prestretch, the modes must start in equilibrium with the rest force, $\varepsilon_{k,i}(0) = j_{k,i}\,w'(\lambda_{e,i})\lambda_{e,i}$ at $\bar\lambda = 1$, solved together with $\lambda_e$. Otherwise the first seconds show a spurious creep.

### 6. Dissipation
Using $\dot\varepsilon_k = (\Gamma_k j_k/(\rho_i n_i))\,X_k$ and $\dot\varepsilon_r = (\Gamma_r j_r/(\rho_i n_i))\,X_r$:
$$\mathcal D_{k,i} = X_{k,i}\,\dot\varepsilon_{k,i} = \frac{\Gamma_{k,i}\,j_{k,i}}{\rho_i n_i}\,X_{k,i}^2 \ge 0, \qquad \mathcal D_{r,i} = \frac{\Gamma_{r,i}\,j_{r,i}}{\rho_i n_i}\,X_{r,i}^2 \ge 0$$
for every state and any prestretch, and for any force dependence of $\Gamma_i$ (catch–slip enters only the positive mobilities). Together with paired reaction kinetics (ratio $k_{on}/k_{off}$ from $\Psi$, with $\partial\Psi_{mech}/\partial c_b$ now including the modal energies through $j_k(c_b)$ and $n(c_b)$):
$$\mathcal D = \sum_i\rho_i\,\mathcal A_i w_i + \sum_i\Big[\sum_k\mathcal D_{k,i} + \mathcal D_{r,i}\Big] \ge 0 .$$

### 7. Limits and checks
- **Linear limit:** linearising around rest gives the compliance used in Step 0,
$$J(\omega) = J_\infty\Big[1 + R\sum_k\frac{a_k}{1 + i\omega\tau_k}\Big] + \frac{(1+R)\,J_\infty}{i\omega\tau_r}, \qquad \tau_k = \frac{2\tau_{off}}{q_k^2},\quad \tau_r = \frac{2\tau_{off}}{q_{min}^2},$$
  with $G^* = 1/J$, i.e. Broedersz's mean-field result discretised.
- **Current model as a special case:** $K = 0$ and a terminal mode with rate $\Gamma_i$ and target $\lambda_e = 1$ (force $f_0$) recovers the present SFR rule. The changes are the mode spectrum and the force-free terminal target (thermodynamics).
- **Optional fast solvent modes** (high-frequency $G''$): the same structure with mobilities $\propto 1/\eta$ and times $\sim\eta\ell_c^4/\kappa$, added in series.
- **Optional stress-dependent spectrum** (Mulla 2019): to be built into $\Psi$, e.g. modal stiffnesses from $\kappa q^4 + f q^2$, not into the rates.

### 8. What is microstructural and what is not

| Ingredient | Source |
|---|---|
| Modal decomposition, rates $\Gamma(q) = q^2/(2\tau_{off})$ | Broedersz Eq. (4) |
| Modal weights $\propto\Delta q/q^2$ | Broedersz Eq. (5) |
| Cut-offs $q_{max} = \pi$, $q_{min} = \pi\ell_c/L$ | crosslink spacing $\ell_c = r_{0,f}(c_b)$; filament length $L$ |
| Turnover time $\tau_{off} = 1/\Gamma = c_b/B$ | our kinetics (catch–slip, bound content) |
| Elastic part $w(\lambda_e)$ | β-model (Holzapfel & Ogden; Unterberger) |
| Series combination, affine averaging | Broedersz (spectra add; $G^* = \rho/15\chi$); microsphere |
| Force-free terminal target | thermodynamic consistency (see SFR dissipation section) |
| Amplitude ratio $R$ | Broedersz Eq. (5) prefactor and the WLC segment modulus: uncertain, to be checked in Step 0 |
| Number of modes $K$ | numerical resolution only |
