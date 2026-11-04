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

## Consistent SFR with an isotensional prestress reservoir

This section supersedes the conclusions of "SFR dissipation with prestretch", §5–6. The two-population alternative (permanent prestressed fraction plus force-free relaxing fraction) was tested and fails (§5 below).

### 1. Physical picture
When a crosslink rebinds, the freed segment re-equilibrates with its neighbours (Broedersz 2010, local re-equilibration). The neighbours are still prestressed, so the segment returns to the prevailing filament tension, i.e. $f_0$, not zero. The energy to re-tension it is supplied by the rest of the filament: the long-wavelength, filament-scale prestress. On the time scale of the experiments this prestress is large and slow. Rigor-HMM networks keep their internal stress over hours (Lieleg 2011), and Fig. 5 shows no flow down to 0.01 Hz.

The negative dissipation of the plain SFR derivation appears because the free energy contained no such reservoir. It is added as an isotensional (constant-tension) reservoir, the mechanical analogue of a chemostat at fixed chemical potential.

### 2. Free energy
With $F(\lambda_e) \equiv w'(\lambda_e)\,\lambda_e = f\lambda_0 r_0\lambda_e \ge 0$ (energy units per filament) and $F_0 = F(1) = f_0\lambda_0 r_0$:
$$\Psi = \Psi_{chem}(c_b, c_f) + \sum_i\rho_i\,n_i(c_{b,i})\Big[w(\lambda_{e,i}; c_{b,i}) + F_0(c_{b,i})\,\ln\lambda_{n,i}\Big], \qquad \lambda_{e,i} = \bar\lambda_i/\lambda_{n,i}.$$
The reservoir term depends only on $\lambda_{n,i}$ (and $c_{b,i}$), not on the deformation.

### 3. Stress
At fixed internal variables, $\partial\Psi_{res}/\partial\mathbf C = \mathbf 0$, so the stress (Eqs. 4.24–4.27 with SFR), the tangent, and the code's `dw = f*l0*r0/lamn` are unchanged.

### 4. Dissipation
**Thermodynamic force on $\ln\lambda_{n,i}$.** With $\partial\lambda_{e,i}/\partial\ln\lambda_{n,i} = -\lambda_{e,i}$:
$$X_i = -\frac{\partial\Psi}{\partial\ln\lambda_{n,i}} = \rho_i n_i\Big[w'(\lambda_{e,i})\,\lambda_{e,i} - F_0\Big] = \rho_i n_i\big[F(\lambda_{e,i}) - F(1)\big].$$
**SFR evolution** (unchanged), in logarithmic form:
$$\frac{d\ln\lambda_{n,i}}{dt} = \frac{\dot\lambda_{n,i}}{\lambda_{n,i}} = \Gamma_i\Big(\frac{\bar\lambda_i}{\lambda_{n,i}} - 1\Big) = \Gamma_i\,(\lambda_{e,i} - 1), \qquad \Gamma_i = \frac{B}{c_{b,i}} > 0 .$$
**Dissipation:**
$$\boxed{\;\mathcal D_i = X_i\,\frac{d\ln\lambda_{n,i}}{dt} = \rho_i n_i\,\Gamma_i\,\big[F(\lambda_{e,i}) - F(1)\big]\big(\lambda_{e,i} - 1\big) \;\ge 0\;}$$
**Proof of the sign.** $F$ is non-decreasing in $\lambda_e$: $F = 0$ when slack ($\lambda_e \le \lambda_e^*$), and for $\lambda_e > \lambda_e^*$ it is the product of the non-negative, increasing filament force and $\lambda_e > 0$. Hence $F(\lambda_e) - F(1)$ has the sign of $\lambda_e - 1$ (or vanishes), so $\mathcal D_i \ge 0$ in every state, including the band $\lambda_e^* < \lambda_e < 1$ that was negative without the reservoir. In that band, the reservoir supplies the energy that re-tensions the segment: $F_0\,|d\ln\lambda_n| > F(\lambda_e)\,|d\ln\lambda_n|$.

**Comparison with the plain derivation:** $\mathcal D_i^{plain} = \rho_i n_i\Gamma_i F(\lambda_{e,i})(\lambda_{e,i}-1)$, and the reservoir contributes $-\rho_i n_i\Gamma_i F_0(\lambda_{e,i}-1)$, i.e. the work exchanged with the tension bath.

### 5. Why the two-population alternative was dropped
Quick check (single point at the `plate_guo_dxs` parameters, not refitted; scratchpad `twopop/`):

| Variant | $G'_{lin}$ | $G''_{lin}$ | Loop area (data 0.30) | $G'(0.158)$/data |
|---|---|---|---|---|
| current SFR | 10.66 | 2.90 | 0.434 | 1.37 |
| permanent $(1-\xi)$ + force-free $\xi = 0.53$ | 6.38 | **0.31** | **0.031** | **7.48** |
| force-free for all bonds | **1.84** | 0.57 | 0.066 | 0.23 |

- The force-free relaxing population drains to zero force and dissipates almost nothing.
- The permanent population never relaxes or ruptures, so the large-amplitude stiffness explodes.

The linear dissipation of the model comes from re-referencing around the prestressed state, which is exactly what the reservoir formulation keeps.

### 6. Consequences and assumptions
- **No change to the model equations, the fit, or the Fortran implementation.** Only the free-energy bookkeeping (and its statement in the paper) changes. No new parameter.
- **Assumption:** the filament-scale prestress is an infinitely large, slow reservoir on the experimental time scale. A finite reservoir would drain slowly (terminal relaxation over minutes to hours), which is outside the data.
- **Chemistry:** $\Psi_{res}$ depends on $c_{b,i}$ through $n_i$ and $F_0$, so it contributes to
$$\mu^{mech}_{b,i} = \frac{1}{\rho_i}\frac{\partial}{\partial c_{b,i}}\Big[\rho_i n_i\big(w + F_0\ln\lambda_{n,i}\big)\Big],$$
  which is the quantity the paired reaction kinetics need.
- **Remaining inconsistency:** the reaction kinetics. The catch–slip $k_{off}$ with a constant $k_{on}$ violates detailed balance. The consistent form has the catch–slip shape $g$ in both rates, $k_{off,i} = g$, $k_{on,i} = K_0\,e^{\chi(1-2\theta_f)}\,e^{-\mu^{mech}_{b,i}/RT}\,g$, with $\mu^{mech}_{b,i}$ including the reservoir term above.

## Finite (slowly draining) prestress reservoir

### 1. Should the prestress decrease over long times?
For a passive network, yes in principle. Without an energy source, stored elastic energy can only be dissipated, so a prestress must eventually relax. The infinite reservoir of the previous section is the limit in which that relaxation is much slower than any simulated process.

What the literature says for long times:
- **Slow relaxation does occur:**
  - Lieleg 2011: fascin-bundled networks release internal stress over hours, with power-law-like dynamics.
  - Burla 2019: stresses frozen in by kinetic trapping during network formation relax over hours.
- **Mechanism (CGD, Broedersz 2010):** uniform tension along a filament relaxes only through many local unbinding events. The mean-field terminal time is $\tau_r \approx (2/\pi^2)\,\tau_{off}(L/\ell_c)^2$, about 2 min for $\tau_{off} \approx 4$ s, $L = 20\,\mu$m and $\ell_c = 1.69\,\mu$m. That is faster than observed, so real networks relax more slowly than this estimate.
- **Rigor-HMM networks stiffen with age** (Lieleg 2011), the opposite trend, probably because of ongoing consolidation (late rigor-bond formation, kinetic trapping) that acts as a source. For our system the net long-time trend is therefore not clearly a decrease.
- **Our data:** Fig. 5 shows no terminal flow down to 0.01 Hz, so any reservoir relaxation time is $\gg 100$ s.

### 2. Formulation
Replace the constant reservoir tension $F_0$ by a reservoir tension $F_{res,i}$ per direction that relaxes slowly. Its initial value is $F_{res,i}(0) = F_0$.

- **Local rebinding (unchanged rate $\Gamma_i = B/c_{b,i}$):** $\lambda_{n,i}$ relaxes toward the state with $F(\lambda_{e,i}) = F_{res,i}$ instead of $F(\lambda_{e,i}) = F_0$. In the free energy, $+F_0\ln\lambda_{n,i}$ becomes $+F_{res,i}\ln\lambda_{n,i}$, so the thermodynamic force on $\ln\lambda_{n,i}$ is
$$X_i = \rho_i n_i\big[F(\lambda_{e,i}) - F_{res,i}\big],$$
  and the evolution must have the same sign as $X_i$. For example, in logarithmic form,
$$\frac{d\ln\lambda_{n,i}}{dt} = \Gamma_i\,\big(\lambda_{e,i} - \lambda_{e,i}^{res}\big), \qquad F(\lambda_{e,i}^{res}) = F_{res,i},$$
  which reduces to the current SFR update when $F_{res,i} = F_0$ ($\lambda_e^{res} = 1$). Because $F$ is non-decreasing, $X_i$ and $\lambda_{e,i} - \lambda_{e,i}^{res}$ share a sign, so
$$\mathcal D_{\lambda,i} = \rho_i n_i\,\Gamma_i\,\big[F(\lambda_{e,i}) - F_{res,i}\big]\big(\lambda_{e,i} - \lambda_{e,i}^{res}\big) \ge 0 .$$
- **The reservoir itself:** two options.
  1. **Simple decay:** $\dot F_{res,i} = -F_{res,i}/\tau_{res}$. This is easy, but the reservoir energy is not explicit, and consistency is not guaranteed: with $F_{res,i}$ inside $\Psi$, the dissipation gains the term $-\rho_i n_i\,\ln\lambda_{n,i}\,\dot F_{res,i}$, which has no definite sign.
  2. **Explicit reservoir energy (preferred for consistency):** give the filament-scale reservoir a finite stored energy, e.g. a filament-scale spring of stretch $\lambda_{s,i}$ with energy $\psi_{res}(\lambda_{s,i})$ and tension $F_{res,i} = \psi'_{res}\lambda_{s,i}$. Let it relax by a gradient flow with a positive mobility $1/\tau_{res}$. Its dissipation is then $\ge 0$ by construction, and the energy it supplies to local re-tensioning is accounted for explicitly.
- **Limit:** $\tau_{res}\to\infty$ recovers the infinite reservoir, $F_{res} = F_0$ (previous section).

### 3. Microstructural meaning and parameter
- $\tau_{res}$ is the terminal relaxation time of filament-scale tension, i.e. the CGD terminal mode. It is set by the turnover time, the filament length and the crosslink spacing, but the mean-field estimate (about 2 min) appears too short compared with observations.
- **Constraint from data:** only a lower bound, $\tau_{res} \gg 100$ s (Fig. 5), possibly hours (Lieleg 2011, Burla 2019).
- **One new parameter,** weakly constrained.

### 4. Assessment
- **LAOS and amplitude sweep (100 s per test):** no effect. The fits are unchanged.
- **Long UEL simulations** (holds, creep or stress relaxation over many minutes, diffusion-dominated problems over $\sim 10^3$ s): the finite reservoir would add a slow loss of stiffness and prestress.
- **For rigor HMM the long-time trend is uncertain** (stiffening with age), so a decreasing reservoir is not clearly right for this system.

**Recommendation:** keep the infinite reservoir ($\tau_{res}\to\infty$) in the current model. State in the paper that it is the limit of a slowly draining filament-scale prestress, justified by Fig. 5 and the persistent prestress of rigor-HMM networks. Implement the finite version (option 2) only if a UEL application runs long enough for it to matter and data are available to constrain $\tau_{res}$.

## Thermodynamically consistent reaction kinetics: all admissible $k_{on}$/$k_{off}$ configurations

**Notation:** $r_i$ = net binding rate of direction $i$ (per unit volume; written $w_i$ or $R_i$ elsewhere, renamed here to avoid clashing with the filament energy $w$ and the gas constant $R$). $RT$ = gas constant times temperature.

### 1. Ingredients
**Free energy** (seminar Eqs. 5.25–5.28, with separate reference energies for free and bound crosslinkers, and the SFR prestress reservoir):
$$\Psi = \mu^0_f c_f + \sum_i\rho_i\mu^0_b c_{b,i} + \Psi_{mix}(\theta_f, \theta_{b,i}) + \sum_i\rho_i\,n_i(c_{b,i})\Big[w(\lambda_{e,i}; c_{b,i}) + F_0(c_{b,i})\ln\lambda_{n,i}\Big]$$

**Chemical potentials:**
$$\mu_f = \mu^0_f + RT\Big[\ln\frac{\theta_f}{1-\theta_f} + \chi(1-2\theta_f)\Big], \qquad
\mu_{b,i} = \mu^0_b + RT\ln\frac{\theta_{b,i}}{1-\theta_{b,i}} + \mu^{mech}_{b,i},$$
$$\mu^{mech}_{b,i} = \frac{1}{\rho_i}\frac{\partial}{\partial c_{b,i}}\Big[\rho_i n_i\big(w + F_0\ln\lambda_{n,i}\big)\Big] \quad\text{(at fixed } \bar\lambda_i, \lambda_{n,i}\text{)}, \qquad x_i \equiv \frac{\mu^{mech}_{b,i}}{RT}.$$

**Affinity and reaction dissipation:** $\mathcal A_i = \mu_f - \mu_{b,i}$, with $\mathcal D_{reac} = \sum_i\rho_i\,\mathcal A_i\,r_i \ge 0$ required.

**Mass-action rate law** (seminar Eq. 5.33):
$$r_i = k_{on,i}\,\frac{\theta_f}{1-\theta_f} - k_{off,i}\,\frac{\theta_{b,i}}{1-\theta_{b,i}} .$$

### 2. The consistency theorem
Shown before: $r_i = P_i\,(\kappa_i\,e^{\mathcal A_i/RT} - 1)$ with $P_i = k_{off,i}\,\theta_{b,i}/(1-\theta_{b,i}) > 0$, and $\mathcal A_i r_i \ge 0$ in all admissible states **if and only if** $\kappa_i \equiv 1$, i.e.
$$\boxed{\;\frac{k_{on,i}}{k_{off,i}} = K_0\;e^{\chi(1-2\theta_f)}\;e^{-x_i}, \qquad K_0 = e^{(\mu^0_f-\mu^0_b)/RT}\;}$$
In the dissociation convention, $K_d = 1/K_0 = k^0_{off}/k^0_{on}$. In the code's normalisation, $K_0 = k^0_{on}C_{F,max}/(k^0_{off}C_{B,max})$.

**Only the ratio is constrained.** The complete family of consistent mass-action kinetics is
$$k_{off,i} = h_i, \qquad k_{on,i} = K_0\,e^{\chi(1-2\theta_f)}\,e^{-x_i}\,h_i, \qquad h_i > 0 \text{ any function of the state}$$
(force, stretch, $\theta$'s, temperature, ...). Equivalently, with a split factor $\phi_i$ (any real; physically $0 \le \phi_i \le 1$):
$$k_{off,i} = \Gamma_i\,e^{\phi_i x_i}, \qquad k_{on,i} = \Gamma_i\,K_0\,e^{\chi(1-2\theta_f)}\,e^{-(1-\phi_i)x_i}, \qquad \Gamma_i > 0 .$$
- **The prefactor $\Gamma_i$ sets the turnover speed,** i.e. how fast equilibrium is approached.
- **The ratio sets the equilibrium occupancy:**
$$\frac{\theta_{b,i}}{1-\theta_{b,i}} = K_0\,e^{\chi(1-2\theta_f)}\,\frac{\theta_f}{1-\theta_f}\,e^{-x_i}.$$
  So mechanics shifts occupancy only through $x_i$, i.e. through the free energy.
- **The Flory–Huggins factor $e^{\chi(1-2\theta_f)}$** may likewise be put on $k_{on}$ (as in Eq. 5.34), on $k_{off}$, or split. Only its presence in the ratio matters.

### 3. Admissible configurations

Recall $x_i = \frac{\mu^{mech}_{b,i}}{RT}, \qquad \mu^{mech}_{b,i} = \frac{1}{\rho_i}\frac{\partial}{\partial c_{b,i}}\Big[\rho_i\,n_i\big(w(\lambda_{e,i}) + F_0\ln\lambda_{n,i}\big)\Big]$

| # | $k_{off,i}$ | $k_{on,i}$ | Consistent? | Remarks |
|---|---|---|---|---|
| C0 | catch–slip in filament force $g(f_i)$ | $k^0_{on}e^{\chi(1-2\theta_f)}$ (no mechanics) | **✗** | current implementation: ratio $\propto 1/g(f_i)$, not $e^{-x_i}$ |
| C1 | $k^0_{off}\,e^{x_i}$ | $k^0_{on}e^{\chi(1-2\theta_f)}$ | ✓ | seminar Eq. (5.35); $\phi = 1$; pure slip-type, set by $\Psi$ |
| C2 | $k^0_{off}\,e^{\phi x_i}$ | $k^0_{on}e^{\chi(1-2\theta_f)}e^{-(1-\phi)x_i}$ | ✓ | Bell-type split of the $\Psi$ bias; $\phi = 0$: all bias in binding (force-suppressed binding only) |
| C3 | $g_i\,e^{\phi x_i}$ | $g_i\,K_0e^{\chi(1-2\theta_f)}e^{-(1-\phi)x_i}$ | ✓ | **catch–slip (or any law) $g_i$ in the shared prefactor**: catch–slip controls turnover, $\Psi$ controls occupancy |
| C4 | $g_i$ (any, e.g. catch–slip or Bell in $f$) | $K_0e^{\chi(1-2\theta_f)}e^{-x_i}\,g_i$ | ✓ | C3 with $\phi = 0$ (all of the $\Psi$ bias in binding, $k_{off}$ unchanged); earlier "version B" |
| C5 | $k^0_{off}\,e^{f_i x_b/k_BT}$ | $k^0_{on}e^{\chi(1-2\theta_f)}\,e^{-f_i x_u/k_BT}$ | ✓ **only if** $\Psi$ gives $x_i = f_i(x_b+x_u)/k_BT$ | classic Bell in force; requires a crosslinker-work term in $\Psi$ (which also enters the stress); $x_u = 0$ gives Bell on $k_{off}$ only |
| C6 | $\sum_p k_{p0}\,e^{f_i x_p/k_BT}$ ($x_c<0$ catch, $x_s>0$ slip) | $K_0\sum_p k_{p0}\,e^{f_i(x_p-\Delta)/k_BT}$ | ✓ **only if** $\Psi$ gives $x_i = f_i\Delta/k_BT$ | two-pathway Bell with microscopic reversibility per pathway (earlier "version A"); force-suppressed binding is then automatic |
| C7 | any $g_i$ | $K_0e^{\chi(1-2\theta_f)}\,g_i$ (no $e^{-x}$) | ✓ **only if** $x_i \equiv 0$ | mechanics-independent equilibrium; consistent only if $\Psi_{mech}$ does not depend on $c_b$ (not our model) |
| ✗ | Bell or catch–slip in $f$ on $k_{off}$ only, free $\Delta x$ | constant | ✗ | unless $\Psi$ is changed as in C5 |
| ✗ | constant | Bell on $k_{on}$ only (free $\Delta x_u$; Wolff-type probe) | ✗ | unless $\Psi$ is changed as in C5 |
| ✗ | seminar Eq. (5.36) (catch–slip in $\mu^{mech}$) | constant | ✗ | only the pure-slip limit $k_{catch}^0 = 0$, $\alpha_s = 1$ (= C1) is consistent |

**Non-mass-action forms** that are consistent by construction:
- **Linear (Onsager):** $r_i = L_i\,\mathcal A_i$ with $L_i \ge 0$ any state function (Salvadori et al.). Exact near equilibrium, a simplification far from it.
- **General dissipation potential:** $r_i = \partial\varphi_i/\partial\mathcal A_i$ with $\varphi_i$ convex in $\mathcal A_i$ and $\varphi_i(0)$ minimal, or any $r_i = \eta_i(\mathcal A_i)$ with $\mathcal A_i\,\eta_i(\mathcal A_i) \ge 0$.
- **Butler–Volmer form:** $r_i = \Gamma'_i\,[e^{(1-\phi)\mathcal A_i/RT} - e^{-\phi\mathcal A_i/RT}]$. This is the mass-action family of §2 rewritten in the affinity; same content.

### 4. Size of $x_i$ in our model
Computed at the `plate_guo_dxs` parameters (scratchpad `kin/mu.py`), per direction, with the energy $w$ obtained by integrating $w' = f\lambda_0 r_0$ from the slack state and $\mu^{mech}_b = \partial(n w)/\partial c_b$ at fixed $\lambda_e$ (units: $n w$ in J/m³, $c_b$ in mol/m³, so $\mu$ in J/mol; $RT = 2444$ J/mol):

| $\lambda_e$ | 0.98 | 1.00 (rest) | 1.02 | 1.04 | 1.06 | 1.07 |
|---|---|---|---|---|---|---|
| $f$ (pN) | 0.03 | 0.185 | 0.57 | 2.0 | 15.8 | 97.7 |
| $x = \mu^{mech}_b/RT$ | 0.0005 | 0.011 | 0.048 | 0.17 | 0.81 | 2.64 |
| effective Bell length $x\,k_BT/f$ (nm) | | 0.24 | 0.34 | 0.34 | 0.21 | 0.11 |

- **Reservoir contribution** ($\lambda_n \ne 1$): $|x_{res}| \approx 0.01$ for $|\ln\lambda_n| = 0.01$, about 0.03 for $\lambda_n = 1.03$. Small.
- **Comparison:** Bell with $\Delta x = 1$ nm gives $f\Delta x/k_BT = 0.25$ at 1 pN and 2.5 at 10 pN. So the $\Psi$-based bias acts like a weak Bell factor with an effective length of about 0.1–0.35 nm.
- **Consequences:**
  - Near rest the consistent kinetics barely differ from the current ones ($e^{-x} \approx 0.99$).
  - In the most loaded directions at $\gamma \approx 0.126$–0.158 (peak forces about 7–18 pN, depending on the fit and on single-point vs plate-edge strain), $x \approx 0.5$–0.9, so $e^{-x} \approx 0.4$–0.6. The equilibrium occupancy there drops by a factor of up to about 2: a moderate, thermodynamically required bond loss under load. The force-dependent-binding probe showed that such redistribution tends to raise $G_M/G'_{lin}$, so its effect on the fit has to be checked.

### 5. Implementation notes
- **Needed quantity:** $x_i(\lambda_{e,i}, c_{b,i}, \lambda_{n,i})$ per direction and step. The code computes $w' = f\lambda_0 r_0$ and $w''$, but not $w$ itself. $w$ requires integrating the (implicit) β-model force from the slack state, $w(\lambda_e) = \lambda_0 r_0\int_{\lambda_e^*}^{\lambda_e} f\,d\lambda$, either numerically per evaluation or from a pre-tabulated function of $(\lambda_e, r_{0f})$. Its $c_b$-derivative includes the $n(c_b)$, $r_{0f}(c_b)$ (and $L = a\,r_{0f}$) dependences, which the code already differentiates for $f$ (`dfdcb_ext`).
- **Kinetics residual:** $k_{on}$ becomes direction-dependent and enters the backward-Euler residual $H(c_b)$ through $r_i$. Its derivative $\partial k_{on,i}/\partial c_{b,i}$ (via $x_i$) must be added to $dH/dc_b$, which affects convergence and uniqueness.
- **SFR update:** the binding flux $B$ that drives $\lambda_n$ becomes $B_i$ (direction-dependent), as in the force-dependent-binding probe.
- **Fortran (folders 7 and 9):** the same changes in the local kinetics (`kineticsFunc.f90`, `solveKinetics`), and new tangent contributions through $k_{on,i}$ in the chemo-mechanical coupling terms (Terms 3 and 5–6 of the tangent split). This adds no `global.f90` variables; per-direction quantities are state variables or are recomputed.
- **Recommended configuration:** **C4** (or C3 with $\phi$ chosen), with the current catch–slip law as $g_i$. It keeps all fitted kinetics as the turnover law, adds no parameter (with C4), and makes the occupancy follow $\Psi$. Before any Fortran work, test it in Python at the reference parameters, then refit $k^0_{off}$ and $\Delta x_s$.

### 6. Computation of $\mu^{mech}_b$ (chain rule)
Implemented in `material_TC._x_mech` (Python) and `xmech.f90` (Fortran, folder 10). Per direction (index $i$ dropped):
$$\mu^{mech}_b = \frac{\partial}{\partial c_b}\Big[n\,\big(w + F_0\ln\lambda_n\big)\Big]_{\bar\lambda,\,\lambda_n}, \qquad x = \frac{\mu^{mech}_b}{RT}.$$
Fixed $\bar\lambda$ and $\lambda_n$ (the other state variables) is the same as fixed $\lambda_e = \bar\lambda/\lambda_n$ and $\lambda_n$. By the product rule:
$$\mu^{mech}_b = \underbrace{n'\,w + n\,\frac{\partial w}{\partial c_b}\Big|_{\lambda_e}}_{\text{compound energy}} \;+\; \underbrace{\ln\lambda_n\,\big(n'\,F_0 + n\,F_0'\big)}_{\text{prestress reservoir}}, \qquad (\cdot)' = \frac{d}{dc_b}.$$

**Step 1: where $c_b$ enters.** Only through the crosslink spacing $r_{0f} = 1.6\,(10^3 c_b)^{-2/5}$ (Tharmann et al.). Everything else is a chain of geometric quantities:

| quantity | definition | derivative with respect to $c_b$ |
|---|---|---|
| $r_{0f}$ | $1.6\,(10^3 c_b)^{-2/5}$ | $r_{0f}' = -\tfrac25\,r_{0f}/c_b$ |
| $r_0$ | $r_{0f} + r_{0c}$ | $r_0' = r_{0f}'$ |
| $L$ | $a\,r_{0f}$ | $L' = -\tfrac25\,L/c_b$ |
| $n$ | $K_{actin}/L$ | $n' = +\tfrac25\,n/c_b$ |
| $\kappa$ | $\eta\,r_0/r_{0f}$ (filament share of the compound stretch) | $\kappa' = -\eta\,r_{0c}\,r_{0f}'/r_{0f}^2$ |
| $\lambda_{0f}$ | $\kappa(\lambda_0 - 1) + 1$ | $\lambda_{0f}' = \kappa'(\lambda_0 - 1)$ |
| $\lambda_f$ | $\kappa(\lambda_e - 1) + 1$ | $\lambda_f' = \kappa'(\lambda_e - 1)$ (at fixed $\lambda_e$) |
| $c$ | $L^2/(\pi^2 B_0)$ | $c' = -\tfrac45\,c/c_b$ |
| $A$ | $\lambda_0 r_0/\kappa$ | $A' = \lambda_0\,(r_{0f}'\kappa - r_0\kappa')/\kappa^2$ |

**Step 2: the compound energy $w$ in force space.** The code only has the force law, $w' = dw/d\lambda_e = f\lambda_0 r_0$, with $f$ given implicitly by $G(f, \lambda_f) = 0$ (`evalg.f90`). Integrating from the slack state ($f = 0$ at $\lambda_f\lambda_{0f} = 1$), with $d\lambda_e = d\lambda_f/\kappa$ and integration by parts:
$$w = \frac{\lambda_0 r_0}{\kappa}\int f\,d\lambda_f = A\,\Big[f\lambda_f - I\Big], \qquad I = \int_0^{f}\lambda_f(f')\,df'.$$
The inverse force law is explicit, which is why the force-space form is used:
$$\lambda_f(f') = \frac{a}{\lambda_{0f}}\,\mathrm{rhs}(f'), \qquad \mathrm{rhs} = u - K\,u^\beta\,v\,D^{-\beta},$$
$$u = 1 + f'/\mu_0, \quad v = 1 + 2f'/\mu_0, \quad D = 1 + c\,f'u, \quad K = 1 - 1/a .$$
Check: $\lambda_f(0) = 1/\lambda_{0f}$, the slack threshold.

**Step 3: differentiate $w$ at fixed $\lambda_e$.** The upper limit $f$ of $I$ depends on $c_b$. By Leibniz's rule:
$$I' = \lambda_f(f)\,f' + J, \qquad J = \int_0^{f}\frac{\partial\lambda_f(f')}{\partial c_b}\Big|_{f'}\,df'.$$
The term $\lambda_f(f)\,f'$ cancels the $\lambda_f f'$ from $(f\lambda_f)'$. So **no $\partial f/\partial c_b$ is needed for $w$**:
$$\frac{\partial w}{\partial c_b}\Big|_{\lambda_e} = A'\,\big[f\lambda_f - I\big] + A\,\big[f\,\kappa'(\lambda_e - 1) - J\big].$$
The integrand of $J$ follows from the explicit inverse law ($c_b$ enters through $c$ and $\lambda_{0f}$):
$$\frac{\partial\lambda_f}{\partial c_b}\Big|_{f'} = \frac{a}{\lambda_{0f}}\,\beta K\,u^{\beta+1}\,v\,f'\,D^{-\beta-1}\,c' \;-\; \lambda_f\,\frac{\lambda_{0f}'}{\lambda_{0f}} .$$

**Step 4: quadrature.** Both $I$ and $J$ use 16-point Gauss–Legendre in $s = \ln(1 + c f')$, $s \in [0, \ln(1 + c f)]$:
$$f' = \frac{e^s - 1}{c}, \qquad df' = \frac{e^s}{c}\,ds .$$
The logarithmic variable puts the nodes where $\lambda_f(f')$ varies fastest, in the entropic regime at small force. The relative error is about $10^{-10}$ up to more than 200 pN. One force solve at $\lambda_e$ is needed (`fil`), plus 16 explicit evaluations.

**Step 5: the reservoir term.** $F_0 = f_0\,\lambda_0 r_0$, with $f_0$ the force at $\lambda_e = 1$ ($\lambda_f = 1$). Here the force itself is differentiated, so the implicit-function derivative is needed:
$$F_0' = \lambda_0\big(r_{0f}'\,f_0 + r_0\,f_0'\big), \qquad f_0' = -\frac{\partial G/\partial c_b}{\partial G/\partial f}\Big|_{\lambda_e = 1},$$
with $f_0'$ from `fil` (`DfDcb`) or `dfdcb_ext`. The term is skipped when $\lambda_n = 1$.

**Step 6: assemble.**
$$\mu^{mech}_b = n'\,w + n\,\frac{\partial w}{\partial c_b}\Big|_{\lambda_e} + \ln\lambda_n\,\big(n'F_0 + n\,F_0'\big), \qquad x = \mu^{mech}_b/RT .$$
- **Units:** $n$ is the number of compounds per µm³, so $n w$ is in pN/µm² = J/m³. With $c_b$ in mM = mol/m³, $\mu^{mech}_b$ is in J/mol.
- **Normalisation:** C4 uses $x - x_{rest}$, with $x_{rest} = x(c_{b0}, \lambda_e = 1, \lambda_n = 1)$.

**Tangent (Fortran only).** The stretch sensitivity of the bias does not need the quadrature. Since $\partial w/\partial\lambda_e = w' = \lambda_0 r_0 f$ and the reservoir does not depend on $\lambda_e$:
$$\frac{\partial x}{\partial\lambda_e} = \frac{1}{RT}\,\frac{\partial(n\,w')}{\partial c_b}\Big|_{\lambda_e} = \frac{n}{RT}\Big[\tfrac25\,\frac{w'}{c_b} + \lambda_0\big(r_{0f}'\,f + r_0\,\partial f/\partial c_b\big)\Big],$$
the same combination as `auxdwdcb` in `affclnetfic_discrete_imp.f90`.

**Verification (2026-11-03, scratchpad):**
- `xmech.f90` against `_x_mech` at η = 1, 0.8, 0.65 ($c_b/c_{b0}$ = 0.8–1.2, $\lambda_e$ = 0.96–1.06, $\lambda_n$ = 0.98–1.02): maximum relative difference 5e-13.
- $\partial x/\partial\lambda_e$ and $\partial f/\partial\lambda_e$ against finite differences: about 1e-9.

## Future extensions (outside the current scope)

Parked items, with the reason each is not pursued now.

### 1. Stress-dependent actin (de)polymerisation per direction
- **Idea:** filament content per direction $n_i$ exchanged with a G-actin pool,
  $\dot n_i = k_{pol}c_G - k_{dep}(f_i)\,n_i$, with total actin conserved
  (constrained-mixture growth and remodelling; same structure as $c_b$/$c_f$).
- **Effect:** if new filaments reset the reference (deposited at the current stretch), the rates just add to $\Gamma$ (no new relaxation time). If turnover changes $n_i$, it is a separate channel with time $1/k_{actin}$, visible in the linear regime through the prestress $F_0$; its weight is set by the force sensitivity of $k_{dep}$.
- **Why not now:** in reconstituted rigor-HMM gels, turnover is far too slow (treadmilling about 0.1–0.3 subunits/s, of order 10 h per 21 µm filament; barbed ends gelsolin-capped; ATP depleted). Supplying the missing 30–60 s modes would need a fitted rate with no microstructural basis.
- **Where it fits:** cells and the cortex (fast, ATP-driven turnover). It needs the ATP free energy as an explicit source term for thermodynamic consistency (treadmilling is active).

### 2. High-frequency $G''$ rise (Fig. 5 above about 1 Hz)
- **Mechanism:** filament–solvent friction of the segment fluctuations, $G^* \propto (i\omega)^{3/4}$ (Gittes–MacKintosh, Morse; the $f^{3/4}$ term of Lieleg 2008), with $\tau \sim \eta L^4/\kappa$.
- **Form:** one fast relaxation channel (free Maxwell fits put it at about 7–9 ms; Unterberger's 7 ms element), or a Prony approximation of the $(i\omega)^{3/4}$ law with weights from $\eta$, $\kappa$ and $L$.
- **Why not now:** no effect at 0.05 Hz (LAOS, amplitude sweep). Needed only if the UEL is used for fast loading (high-rate indentation, impact).

### 3. Finite (slowly draining) prestress reservoir
- See the section "Finite (slowly draining) prestress reservoir" (option 2, explicit reservoir energy).
- **Why not now:** Fig. 5 shows no terminal flow down to 0.01 Hz ($\tau_{res} \gg 100$ s), and rigor-HMM networks stiffen with age (Lieleg 2011). Only relevant for UEL simulations over many minutes, with data to constrain $\tau_{res}$.

### 4. Chain relaxation: refinements and open points (`chain_relaxation_rules.md`)
- **Per-crosslink rates:** $k_{j,j+1}(f)$ from the force at each crosslink instead of one $k = B/(2c_b)$ per direction.
- **Populations following $c_b$:** tabulate $\kappa_p$, $\omega_p$ in $\langle L_{fil}\rangle/L(c_b)$ instead of freezing them at $c_{b0}$.
- **Stress-dependent spectrum (Mulla 2019):** tension-limited bending length $\sqrt{\kappa/f}$ in the mode energies.
- **Open:** the crosslink rate that Fig. 5 needs with $\langle L_{fil}\rangle = 21$ µm ($k \approx 1.5$ s⁻¹) is about 5 times the literature values (0.26–0.3 s⁻¹), and it makes the closed-system kinetics about 13 times more expensive per step.

### 5. Yield onset at large amplitude
- **Limitation:** the model yields at $\gamma_{taut} \approx s_{taut} - 1/s_{taut} \approx 0.17$ (directions reach the taut limit), the data at about 0.2–0.22. Same origin as the excess $G''$ at $\gamma = 0.158$ (about 1.7 times the data).
- **Parked test:** a plate-averaged single-point fit with $\gamma = 0.2$ included, $\beta$ fixed at 0.438 and $\Delta x_s$ free (plate averaging shifts the apparent yield; `c5_w15_plate_guo` reached $G'(0.2)/G'(0.158) = 1.31$ against 1.46 in the data).
- Plate averaging itself is parked for the current fits.

### 6. $G_M/G'_{lin}$ and first-cycle overshoot
- **$G_M/G'_{lin}$ diagnostics:** curvature of the force law near $f_0$ (scan $\beta$, $s_{taut}$ beyond their bounds); slack filaments at $\gamma = 0$ after large strain (Wolff 2012, Fig. 1B).
- **Overshoot:** slip plasticity (Qian 2026) or a slow shakedown variable (Wolff 2012). Unterberger's data show only the steady third cycle, so they do not constrain it; Wolff's staircase data would.

### 7. Fortran ports (after the Python decisions)
- **C4 consistent kinetics:** force-space quadrature for $w$, closed-form tangent terms reusing `dw`, `DfDcb`, `dkoffdf`, `ddw` (corrected reservoir sign).
- **Lumped populations (chain):** $P$ reference stretches per direction as state variables; the $\lambda_n$ term of the tangent summed over populations. No new `global.f90` variables.

## Constitutive assumptions with internal variables (Coleman–Noll, step 2)

Same structure as Sec. 2.3 of the paper draft (`literature/main.pdf`, Eqs. 15–22), starting from the free-energy imbalance (Eq. 14 there),
$$\dot\Psi_R - \boldsymbol\sigma^e_{act} : \dot{\mathbf F}^e - \mu_{act}\,\dot c_R + \mathbf j_R\cdot\nabla_X\mu \le 0 ,$$
but with the state extended by the bound content and the reference stretch of each direction $i = 1 \dots n_{dir}$ (weight $\rho_i$, unit vector $\mathbf m_i$).

### 1. State variables
| variable | meaning | evolves by |
|---|---|---|
| $\mathbf C^e = \mathbf F^{eT}\mathbf F^e$ | elastic right Cauchy–Green tensor | balance of momentum |
| $c_{fR}$ | free (mobile) crosslinker content | balance of mass (transport and reactions) |
| $c_{b,i}$ | bound crosslinker content of direction $i$ | local reaction law, $\dot c_{b,i} = R_i$ |
| $\lambda_{n,i}$ | reference stretch of direction $i$ | local remodelling law, $\dot\lambda_{n,i} = \Lambda_i$ |

The total content and the swelling remain those of the paper: $c_R = c_{fR} + \sum_i \rho_i\,c_{b,i}$, $\dot c_R = -\mathrm{Div}\,\mathbf j_R$, $J^c = 1 + v\,c_R$. Binding converts free into bound crosslinkers at constant $c_R$, so it does not change $J^c$.

### 2. Constitutive equations
At each material point $\mathbf X$ and time $t$, with $\mathbf c_b = \{c_{b,i}\}$ and $\boldsymbol\lambda_n = \{\lambda_{n,i}\}$:
$$\boldsymbol\sigma^e_{act} = \hat{\boldsymbol\sigma}^e_{act}(\mathbf F^e, c_{fR}, \mathbf c_b, \boldsymbol\lambda_n), \qquad
\Psi_R = \hat\Psi_R(\mathbf F^e, c_{fR}, \mathbf c_b, \boldsymbol\lambda_n), \qquad
\mu_{act} = \hat\mu_{act}(\mathbf F^e, c_{fR}, \mathbf c_b, \boldsymbol\lambda_n),$$
$$\mathbf j_R = \hat{\mathbf j}_R(\mathbf F^e, c_{fR}, \mathbf c_b, \boldsymbol\lambda_n, \nabla_X\mu), \qquad
R_i = \hat R_i(\mathbf F^e, c_{fR}, \mathbf c_b, \boldsymbol\lambda_n), \qquad
\Lambda_i = \hat\Lambda_i(\mathbf F^e, c_{fR}, \mathbf c_b, \boldsymbol\lambda_n),$$
subject to $\det\mathbf F^e = J^e = 1$. The new entries are the internal variables $\mathbf c_b$, $\boldsymbol\lambda_n$ in every response function and their two evolution laws.

**Frame indifference.** The dependence on $\mathbf F^e$ is replaced by $\mathbf C^e$. In the affine microsphere model it enters only through the directional stretches
$$\bar\lambda_i = \sqrt{\mathbf m_i\cdot\mathbf C^e\,\mathbf m_i}\,,$$
which are invariant under rotations of the current configuration:
$$\boldsymbol\sigma^e_{act} = \mathbf F^e\,\hat{\mathbf S}(\mathbf C^e, c_{fR}, \mathbf c_b, \boldsymbol\lambda_n), \qquad \Psi_R = \hat\Psi_R(\mathbf C^e, c_{fR}, \mathbf c_b, \boldsymbol\lambda_n), \ \dots$$
with $\hat{\mathbf S}$ symmetric. As in the paper (Eq. 16), $\boldsymbol\sigma^e_{act}$ is the two-point (first Piola–Kirchhoff-type) stress conjugate to $\dot{\mathbf F}^e$, hence a single $\mathbf F^e$; the Cauchy stress follows only at the end, by push-forward (Eq. 22).

**Free energy, summarised.** Three parts: reference chemical potentials, mechanics, and mixing.
$$\hat\Psi_R = \Psi_{ref}(c_{fR}, \mathbf c_b) + \Psi_{mech}(\mathbf C^e, \mathbf c_b, \boldsymbol\lambda_n) + \Psi_{mix}(c_{fR}, \mathbf c_b),$$
$$\Psi_{ref} = \mu_0\,c_{fR} + \sum_i\rho_i\,\mu_b^0\,c_{b,i}, \qquad
\Psi_{mech} = \sum_i\rho_i\,n(c_{b,i})\Big[w(\lambda_{e,i}) + F_0\ln\lambda_{n,i}\Big], \qquad
\Psi_{mix} = \psi_{mix}^{f}(\theta_f) + \sum_i\rho_i\,\psi_{mix}^{b}(\theta_{b,i}).$$
- $\Psi_{ref}$: energy of putting a crosslinker in the free or in the bound state; only $\mu_b^0 - \mu_0$ (the binding free energy) matters.
- $\Psi_{mech}$: elastic energy of the filament–crosslinker compounds ($\lambda_0$ inside $w$) and the prestress reservoir; it couples to the chemistry through $n(c_b)$ and $r_{0f}(c_b)$.
- $\Psi_{mix}$: configurational entropy (and, for the free species, the Flory–Huggins interaction $\chi$) of each species on its own lattice.

The expanded form below writes out $\psi_{mix}^f$ and $\psi_{mix}^b$.

**Free energy of the model, expanded** (the structure to which the restrictions below are applied), with $\lambda_{e,i} = \bar\lambda_i/\lambda_{n,i}$, $\theta_f = c_{fR}/C_{F,max}$ and $\theta_{b,i} = c_{b,i}/C_{B,max}$:
$$\hat\Psi_R = \sum_i \rho_i\Big[\,\underbrace{n(c_{b,i})\,w(\lambda_{e,i}; c_{b,i})}_{\text{compounds (}\lambda_0\text{ inside }w\text{)}}
+ \underbrace{n(c_{b,i})\,F_0(c_{b,i})\ln\lambda_{n,i}}_{\text{prestress reservoir}}
+ \underbrace{\mu_b^0\,c_{b,i} + RT\,C_{B,max}\big[\theta_{b,i}\ln\theta_{b,i} + (1-\theta_{b,i})\ln(1-\theta_{b,i})\big]}_{\text{bound species: binding energy and site mixing}}\Big]$$
$$\qquad + \underbrace{\mu_0\,c_{fR} + RT\,C_{F,max}\big[\theta_f\ln\theta_f + (1-\theta_f)\ln(1-\theta_f) + \chi\,\theta_f(1-\theta_f)\big]}_{\text{free species: reference potential and Flory–Huggins mixing}} ,$$
with $n = K_{actin}/(a\,r_{0f})$, $r_{0f} = 1.6\,(10^3 c_b)^{-2/5}$, $F = w'(\lambda_e)\lambda_e$, $F_0 = F(1)$, $C_{F,max} = c_A R_{F,max}$, $C_{B,max} = c_A R_{B,max}$.
- Free and bound crosslinkers occupy two separate lattices (free sites, capacity $C_{F,max}$; binding sites of each direction, capacity $C_{B,max}$), so the free mixing term depends only on $\theta_f$ and each direction has its own bound mixing term in $\theta_{b,i}$. The bound term is the one implied by the kinetics: the unbinding term $k_{off}\,C_{B,max}\,\theta_b/(1-\theta_b)$ of $R_i$ (`_kinetics_residual`) gives, at $R_i = 0$, the equilibrium of two Langmuir lattices (Oriani-type trapping).
- Only the difference $\mu_b^0 - \mu_0$ (the binding free energy, contained in $K_0$) is physical; $\mu_0 = 0$ (`MU0`) in the parameter set.
- The volumetric part depends on the total content $c_R$ through $J^c$: the incompressibility constraint $J^e = 1$ (as in the paper), or the `KBULK` penalty of the code, in which case $\Psi_{vol}(J^e)$ is added to $\hat\Psi_R$.

### 3. Chain rule
$$\dot\Psi_R = \frac{\partial\hat\Psi_R}{\partial\mathbf C^e} : \dot{\mathbf C}^e
+ \frac{\partial\hat\Psi_R}{\partial c_{fR}}\,\dot c_{fR}
+ \sum_i \frac{\partial\hat\Psi_R}{\partial c_{b,i}}\,\dot c_{b,i}
+ \sum_i \frac{\partial\hat\Psi_R}{\partial\lambda_{n,i}}\,\dot\lambda_{n,i},$$
each derivative at fixed values of the other state variables. With $\dot{\mathbf C}^e = \dot{\mathbf F}^{eT}\mathbf F^e + \mathbf F^{eT}\dot{\mathbf F}^e$, the first term is $2\mathbf F^e\,\partial\hat\Psi_R/\partial\mathbf C^e : \dot{\mathbf F}^e$. With $\dot c_R = \dot c_{fR} + \sum_i\rho_i\,\dot c_{b,i}$, the dissipation inequality becomes
$$\Big(2\mathbf F^e\frac{\partial\hat\Psi_R}{\partial\mathbf C^e} - \mathbf F^e\hat{\mathbf S}\Big) : \dot{\mathbf F}^e
+ \Big(\frac{\partial\hat\Psi_R}{\partial c_{fR}} - \hat\mu_{act}\Big)\dot c_{fR}
+ \sum_i\Big(\frac{\partial\hat\Psi_R}{\partial c_{b,i}} - \rho_i\,\hat\mu_{act}\Big)R_i
+ \sum_i\frac{\partial\hat\Psi_R}{\partial\lambda_{n,i}}\,\Lambda_i
+ \hat{\mathbf j}_R\cdot\nabla_X\mu \;\le\; 0 .$$

### 4. State relations
$\dot{\mathbf F}^e$ and $\dot c_{fR}$ can be prescribed arbitrarily (through the boundary conditions and the flux divergence), while the response functions do not depend on them. Their coefficients must therefore vanish:
$$\boldsymbol\sigma^e_{act} = 2\mathbf F^e\,\frac{\partial\hat\Psi_R}{\partial\mathbf C^e}, \qquad
\mu_{act} = \frac{\partial\hat\Psi_R}{\partial c_{fR}},$$
as in Eq. (19) of the paper, now evaluated at fixed $\mathbf c_b$ and $\boldsymbol\lambda_n$. With the constraints (Eqs. 12–13 of the paper) this gives Eqs. (21)–(22) unchanged in form:
$$\mathbf P = 2\mathbf F^e\frac{\partial\hat\Psi_R}{\partial\mathbf C^e}\mathbf F^{e-T} - P\,\mathbf F^{-T}, \qquad
\mu = \frac{\partial\hat\Psi_R}{\partial c_{fR}} + \bar p\,v .$$
In the microsphere form the stress uses $\partial\hat\Psi_R/\partial\bar\lambda_i = \rho_i\,n\,w'(\lambda_{e,i})/\lambda_{n,i}$, i.e. $w' = f\,r_0\lambda_0/\lambda_n$ (the $1/\lambda_n$ is the chain rule through $\lambda_e$).

### 5. Residual dissipation inequality
The rates $R_i$, $\Lambda_i$ and the flux $\mathbf j_R$ are given by response functions, so what remains is the restriction on them:
$$\underbrace{-\,\hat{\mathbf j}_R\cdot\nabla_X\mu}_{\text{transport}}
\;+\; \underbrace{\sum_i \rho_i\,\big(\mu_f - \mu_{b,i}\big)\,R_i}_{\text{reaction}}
\;+\; \underbrace{\sum_i \Big(-\frac{\partial\hat\Psi_R}{\partial\lambda_{n,i}}\Big)\Lambda_i}_{\text{remodelling}} \;\ge\; 0,$$
with the chemical potentials of the free and bound crosslinkers
$$\mu_f = \frac{\partial\hat\Psi_R}{\partial c_{fR}} = \mu_0 + RT\Big[\ln\frac{\theta_f}{1-\theta_f} + \chi(1-2\theta_f)\Big],$$
$$\mu_{b,i} = \frac{1}{\rho_i}\frac{\partial\hat\Psi_R}{\partial c_{b,i}}\Big|_{\lambda_{e,i},\,\lambda_{n,i}}
= \mu_b^0 + RT\ln\frac{\theta_{b,i}}{1-\theta_{b,i}} + \underbrace{\frac{\partial\big[n\,(w + F_0\ln\lambda_n)\big]}{\partial c_b}\Big|_{\lambda_{e,i},\,\lambda_{n,i}}}_{\mu_b^{mech}\ (\texttt{\_x\_mech})} .$$
The affinity $\mu_f - \mu_{b,i}$ therefore contains exactly the terms of the kinetics: the free and bound mixing logarithms, $\chi$, the binding free energy ($K_0$) and the mechanical part.
The pressure term $\bar p\,v$ cancels in the affinity $\mu_f - \mu_{b,i}$, since binding does not change $c_R$. As in the paper (Eq. 20), each process is required to dissipate on its own (sufficient condition):
- **transport:** $\hat{\mathbf j}_R\cdot\nabla_X\mu \le 0$;
- **reaction:** $\rho_i\,(\mu_f - \mu_{b,i})\,R_i \ge 0$ (the net binding rate follows the affinity);
- **remodelling:** $-\,\partial\hat\Psi_R/\partial\lambda_{n,i}\;\Lambda_i \ge 0$, with $-\partial\hat\Psi_R/\partial\lambda_{n,i} = \rho_i\,n\,\big(F(\lambda_{e,i}) - F_0\big)/\lambda_{n,i}$ (the reservoir term supplies the $-F_0$).

**Choices of the model that satisfy them** (step 5 of the procedure, derived in the sections above):
- transport: $\mathbf j_R = -\mathbf M\nabla_X\mu$, $\mathbf M = D\,c_f(1-\theta_f)/(J\,RT)\,\mathbf I \ge 0$ (Chester–Anand);
- reaction: $k_{on}/k_{off} = K_0\,e^{\chi(1-2\theta_f)}\,e^{-\mu_b^{mech}/RT}$ (with the mixing factors $\theta_f/(1-\theta_f)$ in $B_i$ and $\theta_{b,i}/(1-\theta_{b,i})$ in the unbinding term of $R_i$) and the catch–slip law in a shared positive prefactor, so that $\mathrm{sign}(R_i) = \mathrm{sign}(\mu_f - \mu_{b,i})$ (mass-action kinetics with the rate ratio set by the chemical potentials, as in Arricca et al. 2023, JMPS 181:105425); the legacy kinetics ($k_{on}$ independent of $\mu_b^{mech}$) does not satisfy it;
- remodelling: $\Lambda_i = \Gamma_i\,\lambda_{n,i}(\lambda_{e,i} - 1)$, $\Gamma_i = B_i/c_{b,i}$; since $F$ is non-decreasing in $\lambda_e$, $(F - F_0)(\lambda_e - 1) \ge 0$ for every state (nonlinear, sign-preserving; not of Onsager form).

**To state explicitly:** the reaction and the remodelling share the binding flux $B_i$ (in $R_i$ and in $\Gamma_i$); requiring each term to be non-negative separately is sufficient, and it relies on the split of $\mu_{b,i}$ at fixed $\lambda_{e,i}$, $\lambda_{n,i}$ (new bonds form at the current stretch, which is accounted for in $\Lambda_i$).


---------------------------------------------

$$\rho_i,(\mu_f - \mu_{b,i}),R_i ;\ge; 0 ,$$
so the net binding rate R_i must have the same sign as the affinity A_i = μ_f − μ_b,i:
- free crosslinkers "higher" (A_i > 0): net binding, R_i > 0;
- bound crosslinkers "higher" (A_i < 0): net unbinding, R_i < 0;
- A_i = 0 (equilibrium): R_i = 0.

Two ways to guarantee it

1. Linear (Onsager): R_i = L A_i with L ≥ 0, so R_i A_i = L A_i² ≥ 0. This is the same trick as the flux, but it is only accurate close to equilibrium.
2. Mass action with detailed balance (Arricca et al.; your C4): write the rate as (unbinding flux) × (ratio − 1),
   $$R_i = \underbrace{k_{off},C_{B,max}\frac{\theta_{b,i}}{1-\theta_{b,i}}}_{>0};\Big(e^{A_i/RT} - 1\Big).$$
   Since e^x − 1 has the sign of x, R_i·A_i ≥ 0 for every state, however far from equilibrium.

Why C4 has this form and the legacy kinetics doesn't

Your rate is R_i = B − k_off C_B,max θ_b/(1 − θ_b), with B = k_on C_F,max θ_f/(1 − θ_f)·e^{χ(1−2θ_f)}. The ratio of binding to unbinding flux is
$$\frac{B}{k_{off}C_{B,max},\theta_b/(1-\theta_b)} = \frac{k_{on}}{k_{off}},\frac{C_{F,max}}{C_{B,max}},\frac{\theta_f}{1-\theta_f},\frac{1-\theta_b}{\theta_b},e^{\chi(1-2\theta_f)} .$$
- C4: k_on/k_off = K₀e^{χ(1−2θ_f)}e^{−μ_b^mech/RT}, with K₀ set by μ₀ − μ_b⁰. The ratio is exactly e^{A_i/RT}, so the sign condition holds always.
- Legacy: k_on/k_off has no e^{−μ_b^mech/RT}. The ratio is e^{(A_i + μ_b^mech)/RT}, so R_i follows A_i + μ_b^mech instead of A_i. Whenever the affinity is small and opposite in sign to μ_b^mech (|A_i| < |μ_b^mech|), the reaction runs against its affinity and R_i A_i < 0.
The full condition

Write k_on = k_on0 e^G and k_off = k_off0 e^F. With your rate,
$$\frac{B}{k_{off}C_{B,max}\theta_b/(1-\theta_b)} = \frac{k_{on0}}{k_{off0}},e^{G-F},\frac{C_{F,max}}{C_{B,max}},\frac{\theta_f}{1-\theta_f},\frac{1-\theta_b}{\theta_b},$$
and the affinity is
$$\frac{A_i}{RT} = \frac{\mu_0 - \mu_b^0}{RT} + \ln\frac{\theta_f}{1-\theta_f} + \chi(1-2\theta_f) - \ln\frac{\theta_b}{1-\theta_b} - \frac{\mu_b^{mech}}{RT}.$$
Setting the ratio equal to e^{A_i/RT}, the mixing logarithms cancel, leaving
$$G - F = \chi(1-2\theta_f) - \frac{\mu_b^{mech}}{RT} ;+; \underbrace{\frac{\mu_0 - \mu_b^0}{RT} - \ln\frac{C_{F,max}}{C_{B,max}} - \ln\frac{k_{on0}}{k_{off0}}}_{\text{constant}} .$$

Your expression is this one with the constant reduced to −ln(k_on0/k_off0). That's what you get if μ₀ = μ_b⁰ and C_F,max = C_B,max. Neither holds in general: in your parameter set R_F,max = 1 and R_B,max = 0.25, so C_F,max/C_B,max = 4.


----


![alt text](image.png)