# Chain relaxation: nonlinear rules

Per microsphere direction (the direction index $i$ is dropped). Linear checks: scratchpad `chain/chain_linear.py`, `chain/chain_reduced.py`.

## 1. Idea

```
current model:   [seg]   [seg]   [seg]   [seg]        independent segments, each unbinding erases its stress
chain:         ~~[seg]x[seg]x[seg]x[seg]~~            one filament, crosslinks x pin it to the network
                 free end                  free end
```

- **Current model:** a crosslink event resets one segment toward the rest state.
- **Chain:** an interior event only lets two neighbouring segments share their stored length. Stress leaves the filament only at its ends.
- So interior stress has to *diffuse* to the ends, and long filaments relax slowly.

## 2. Variables

**Glossary.** Physical objects (they exist in the network):

| term | what it is | how many | symbol |
|---|---|---|---|
| **Filament** | one actin filament, a continuous polymer of length $L$ | many per direction; lengths random (exponential, mean $\langle L\rangle$, §6) | $L$ |
| **Junction** | a point where a bound crosslinker attaches a filament to the network (node) | $K$ per filament, random (Poisson, mean $L/\ell_c$, §6) | $K$ |
| **Segment** | the piece of a filament between two neighbouring junctions; one Holzapfel compound, carrying force $f_j$ | $N = K - 1$ load-bearing segments per filament (the two end pieces dangle) | $j = 1 \dots N$, $\lambda_{n,j}$ |

Groupings (bookkeeping over the physical objects):

| term | what it is | why it exists |
|---|---|---|
| **Class** | all filaments with the same $N$ | filaments with equal $N$ have identical dynamics (same matrix $A_N$), so they are treated together; weight $\omega_N$ (§6) |
| **Mode** | an eigenvector of $A_N$: a lag pattern along a filament of class $N$ that decays at one rate $k\mu_{N,q}$ | turns the coupled segments of a class into independent variables (§5); one Maxwell element each |
| **Population** | a group of modes, from all classes, with similar rates, lumped into one stretch-free reference $\lambda_{n,p}$ with weight $\omega_p$ and rate $\kappa_p$ | numerical compression of the thousands of modes into about 6 variables for the UMAT (§8) |

```
direction i
 └─ filaments (random lengths L)                         physical
     └─ segments j = 1..N  (between junctions)           physical, carry force

grouped as:
 class N       = all filaments with N segments           -> weight omega_N          (§6)
   └─ modes q  = eigenvectors of A_N, rate k mu_{N,q}    -> weight omega_N * w_{N,q} (§5)
 population p  = modes of ALL classes with similar rate  -> omega_p, kappa_p, lambda_n,p (§8)
```
- A class groups filaments; a population groups modes and cuts across classes (the slowest population holds the $q = 1$ modes of the long-filament classes; the fastest holds the $N = 1$–2 classes and the high-$q$ modes of the rest).
- Segments and modes are two descriptions of the same filament: $N$ coupled segment references or, equivalently, $N$ independent mode amplitudes (related by the eigenvectors).

In `material_chain.py`:

| `relaxation` | stored per direction | meaning |
|---|---|---|
| `'sfr'` | 1 reference | one segment per filament, one class, one population: the current model |
| `'chain'` | segment references of selected classes (`chain_N`, half-chains: 154 in total) | physical segments of representative classes; exact mean-field |
| `'pops'` | `n_pops` = 6 references | populations; approximate outside the linear regime |

**Variables of one filament** of $N$ segments, $j = 1 \dots N$:

| Symbol | Meaning |
|---|---|
| $\bar\lambda$ | directional stretch (affine, **the same for all segments**, since every crosslink is pinned to the network) |
| $\lambda_{n,j}$ | reference stretch of segment $j$; $\lambda_{n,j}\,r_0$ is its zero-force end-to-end length. Sliding moves contour length; with $a = L/r_{0f}$ fixed this is represented by moving reference end-to-end length (see §4, reference stretch vs contour length) |
| $\lambda_{e,j} = \bar\lambda/\lambda_{n,j}$ | elastic stretch of segment $j$ |
| $f_j = f(\lambda_{e,j})$ | compound force (filament + linker in series, $f_f = f_c$), the existing force law |
| $F_j = w'(\lambda_{e,j})\,\lambda_{e,j}$ | the same force in energy units, non-decreasing in $\lambda_{e}$; $F_0 = F(1)$ is the prestress |

Rest state: $\lambda_{n,j} = 1$ for all $j$, so every segment carries $f_0$.

## 3. Stress

The same expression as now, with the segment term averaged along the filament:
$$\frac{w'(\lambda_e)}{\lambda_n} \;\longrightarrow\; \Big\langle \frac{w'(\lambda_{e,j})}{\lambda_{n,j}} \Big\rangle_j = \frac{1}{N}\sum_{j=1}^{N}\frac{w'(\lambda_{e,j})}{\lambda_{n,j}} .$$
Every segment counts equally: each has the same $\bar\lambda$ and the same $n$.

## 4. Event rules (what one crosslink does)

**Topology and force balance at a junction.** Each segment is a Holzapfel compound: a filament piece in series with a linker, so inside segment $j$ the filament and the linker carry the same force, $f_f = f_c = f_j$. The linkers attach the segments to a network junction (node) that moves affinely.

```
 bound junction (one crosslinker)                         unbound junction

   segment j                     segment j+1                segment j          segment j+1
 ~~~~~~~~~~~~~~[ c | c ]~~~~~~~~~~~~~~                    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      f_j          |          f_j+1                            f_j    -->|<--    f_j+1
                   * node (affine, network)                    slides until equal
                   | F_node                                 (crosslinker released, F_node = 0)

 ~~~ filament piece   [ c | c ] one crosslinker: each half is in series with its own piece
                      (half of its compliance belongs to each compound)   * network junction
```

- **One crosslinker per junction:** a filament of $N$ segments has $N$ linker lengths ($r_0 = r_{0f} + r_{0c}$ per compound) and about $N$ junctions, so each compound owns half of the crosslinker at each of its ends. Both halves belong to one crosslinker, so one unbinding event releases the junction (the rate $k$ is per junction).
- **Bound:** node equilibrium $f_j - f_{j+1} + F_{node} = 0$. Each crosslinker half carries its own segment's full force; the difference $f_j - f_{j+1}$ is carried by the node, i.e. by the affine network. This is the same imbalance the affine microsphere model accepts at every junction, where compounds of different directions meet. So neighbouring segments can carry different forces (different $\lambda_{n,j}$) between events.
- **Unbound:** the crosslinker lets go (both halves at once), the two filament pieces are joined directly, and the node no longer acts on them ($F_{node} = 0$). The free-body diagram then gives $f_j = f_{j+1}$: the filament slides through the junction until the forces are equal, with the total stored length conserved.
- **Rebinding:** the crosslinker reattaches at the node, which moves with the network, so both halves again have the stretch $\bar\lambda$ and each compound again satisfies $f_f = f_c$, now with the equalised force.
- **Assumption (hybrid topology):** the elasticity follows Holzapfel (the filament effectively ends at each junction, and the linker carries the full force), while the relaxation follows CGD (the filament is continuous through the junction and slides when it unbinds). A real crosslink is a side attachment on a continuous filament; combining the two pictures is a modelling assumption.
- **Slack segments:** if both halves are slack ($f = 0$), force equality does not fix the split of stored length; the average below is then a choice, with no mechanical effect.
- **Reference stretch vs contour length:** sliding physically moves contour length between segments ($L_j$ changes). Exchanging $\lambda_n$ instead captures the change of extension ratio ($\delta\ln\lambda_e = -\delta L/L$) but keeps the force scale $B_0\pi^2/L^2$ fixed. With $d\ln f/d\ln\lambda_e \approx 56$–160 for the compound law (0.1–16 pN, reference parameters), the missing term is about $2/60 \approx 3\,\%$ of the force change near rest and less at higher force. Exchanging $L_j$ per segment would be exact but needs per-segment geometry in the force solve and the tangent.

**Interior crosslink between $j$ and $j+1$ unbinds.** The filament slides through it until the two forces are equal, and the stored length is conserved:
$$\lambda_{n,j} + \lambda_{n,j+1} = \text{const}, \qquad f(\bar\lambda/\lambda_{n,j}) = f(\bar\lambda/\lambda_{n,j+1}).$$
Both segments have the same $\bar\lambda$ and the same force law, so equal forces mean equal references:
$$\boxed{\lambda_{n,j},\ \lambda_{n,j+1} \;\to\; \tfrac12(\lambda_{n,j} + \lambda_{n,j+1})}$$
*Intuition:* sliding makes the two neighbours identical. This holds exactly in the nonlinear regime; no linearisation is needed.

**End crosslink unbinds.** The end segment slides into the free end and re-attaches to the prestressed network, i.e. the current stretch-free rebinding rule:
$$\boxed{\lambda_{n,1} \to \bar\lambda \quad(\lambda_{e,1} = 1,\ f = f_0)}$$
*Intuition:* the network outside the filament is the prestress reservoir (Section 7).

**A new crosslink binds inside segment $j$.** This splits the segment into two parts that inherit the same $\lambda_{n,j}$, so no force changes. The bound content $c_b$ therefore enters only through the spacing $\ell_c(c_b) = a\,r_{0f}(c_b)$, i.e. through the force law (as now) and the number of segments $N = L/\ell_c$.

## 5. Averaged (mean-field) equations

Each crosslink unbinds at rate $k$. Averaging the event rules gives
$$\text{interior:}\quad \frac{d\lambda_{n,j}}{dt} = \frac{k}{2}\big(\lambda_{n,j-1} - 2\lambda_{n,j} + \lambda_{n,j+1}\big),$$
$$\text{end:}\quad \frac{d\lambda_{n,1}}{dt} = k\,\big(\bar\lambda - \lambda_{n,1}\big) + \frac{k}{2}\big(\lambda_{n,2} - \lambda_{n,1}\big) \qquad(\text{same at } j = N).$$
- *Interior:* a discrete diffusion equation for the stored length. Neighbours pull each other toward their mean.
- *End:* the end segment is also pulled toward $\bar\lambda$, as in the current model.
- The equations are linear in $\lambda_{n,j}$ even though the forces are nonlinear, because all segments share $\bar\lambda$. The nonlinearity enters through the force law in the stress and through the rates.

**Matrix form.** Stacking the reference stretches of one filament, $\boldsymbol\lambda_n = (\lambda_{n,1},\dots,\lambda_{n,N})^T$, the rows above are linear in $\boldsymbol\lambda_n$ and $\bar\lambda$:
- interior row $j$: $\dot\lambda_{n,j} = k\big[\tfrac12\lambda_{n,j-1} - \lambda_{n,j} + \tfrac12\lambda_{n,j+1}\big]$;
- end row $1$: $\dot\lambda_{n,1} = k\big[(-\tfrac12 - 1)\lambda_{n,1} + \tfrac12\lambda_{n,2} + 1\cdot\bar\lambda\big]$ ($-\tfrac12$: exchange with the neighbour; $-1$: release at the end crosslinker; $+1\cdot\bar\lambda$: end target). Row $N$ is the mirror image.

Collecting the coefficients:
$$\frac{d\boldsymbol\lambda_n}{dt} = k\,\big(A\,\boldsymbol\lambda_n + \mathbf b\,\bar\lambda\big),
\qquad\text{e.g. } N = 4:\quad
A = \begin{pmatrix} -\tfrac32 & \tfrac12 & 0 & 0\\ \tfrac12 & -1 & \tfrac12 & 0\\ 0 & \tfrac12 & -1 & \tfrac12\\ 0 & 0 & \tfrac12 & -\tfrac32 \end{pmatrix},
\quad \mathbf b = \begin{pmatrix}1\\0\\0\\1\end{pmatrix}.$$
- $A$ is tridiagonal (each segment exchanges only with its neighbours) and symmetric; $\mathbf b$ is non-zero only at the ends, where the target $\bar\lambda$ enters.
- Consistency: every row of $A$ plus $\mathbf b$ sums to zero, $A\mathbf 1 + \mathbf b = \mathbf 0$, so the rest state $\boldsymbol\lambda_n = \bar\lambda\,\mathbf 1$ is stationary.
- Time stepping (prototype, `_chain_matrix`, `_relax_segments`): backward Euler, one small linear solve per direction and step,
$$\big(I - k\Delta t\,A\big)\,\boldsymbol\lambda_n^{t+\Delta t} = \boldsymbol\lambda_n^{t} + k\Delta t\,\mathbf b\,\bar\lambda^{t+\Delta t}.$$
- Spectrum: the eigenvalues of $-A$ (times $k$) are the relaxation rates of the $N$ modes, and the projection of the uniform initial state on the eigenvectors gives their stress weights. Doing this for every $N$ of the length distribution (§6) and binning the modes by rate gives the populations $\kappa_p$, $\omega_p$ (§8). This exact discrete spectrum, not the continuum $r_q$ below, is what the model uses.

**The current model is the case $N = 1$.** One segment between two end crosslinks gives $\dot\lambda_n = 2k(\bar\lambda - \lambda_n)$, which is the current rule $\dot\lambda_n = \Gamma(\bar\lambda - \lambda_n)$ with
$$k = \frac{\Gamma}{2} = \frac{B}{2\,c_b}.$$

**Load dependence:** in the current model $\Gamma = B/c_b$ is the same for all bonds of a direction; force enters through $c_b$ (the catch–slip $k_{off}(f)$ in the kinetics). The prototype keeps this: one $k = B/(2c_b)$ per direction, and the kinetics of $c_b$ use the segment-averaged off-rate $\langle k_{off}(f_j)\rangle$. A per-crosslink rate $k_{j,j+1} = \tfrac12\Gamma\big(\tfrac12(f_j + f_{j+1})\big)$ is a possible refinement.

**Continuum picture.** For long filaments ($N \gg 1$) the discrete chain becomes a diffusion equation for the stored length along the filament. It is used only to understand the scaling ($D$, $\tau_{term}$); the prototype uses the discrete chains or the populations.

*From the discrete chain to a diffusion equation.* Let $s$ be the arclength along the filament and $\ell_c$ the contour length between crosslinks ($\ell_c = L = a\,r_{0f}$ of the β-model). Segment $j$ occupies $s \in [(j-1)\ell_c,\ j\ell_c]$, and $\lambda_{n,j}$ is the value of a smooth profile $\lambda_n(s)$ at its centre. A Taylor expansion of the neighbours gives
$$\lambda_{n,j\pm1} = \lambda_n \pm \ell_c\,\partial_s\lambda_n + \tfrac12\,\ell_c^2\,\partial_s^2\lambda_n \pm \dots
\quad\Rightarrow\quad
\lambda_{n,j-1} - 2\lambda_{n,j} + \lambda_{n,j+1} \approx \ell_c^2\,\partial_s^2\lambda_n .$$
Inserted into the interior equation:
$$\frac{\partial\lambda_n}{\partial t} = D\,\frac{\partial^2\lambda_n}{\partial s^2}, \qquad D = \frac{k\,\ell_c^2}{2}\quad[\mu\text{m}^2/\text{s}].$$
Stored length spreads along the filament by about one segment length per unbinding time $1/k$: a random walk.

*Discrete model, continuous approximation.* The physical model is discrete: each segment has one $\lambda_{n,j}$ and a uniform tension, so along the filament $\lambda_n$ is a staircase of $N$ steps, and the equations above for that staircase are the model. The profile $\lambda_n(s)$ is a smooth curve through the step values, valid when $\lambda_n$ changes slowly from segment to segment. Placing each value at the segment centre $s_j = (j - \tfrac12)\ell_c$ is a labelling convention: it makes the neighbour differences symmetric (error $O(\ell_c^2)$) and puts the boundary exactly at the end crosslinker ($s = 0$). Unlike a finite-difference scheme, the grid here is physical: the continuum approximates the discrete chain, not the other way round. It holds for modes whose wavelength spans many segments ($q \ll N$); the rates
$$r_q^{discrete} \approx 2k\,\sin^2\!\Big(\frac{q\pi}{2N}\Big) \qquad\text{vs}\qquad r_q^{continuum} = \frac{k}{2}\Big(\frac{q\pi}{N}\Big)^2$$
agree for small $q/N$ (slowest rate for $N = 10$: $0.0489\,k$ against $0.0494\,k$), but the discrete rate saturates at about $2k$ for $q \to N$ while the continuum keeps growing. Hence the prototype uses the discrete chains or the populations, and the continuum only for scaling ($D$, $\tau_{term}$).

*Boundary condition.* The end equation is the interior equation with a ghost segment $0$ outside the filament, set to $\lambda_{n,0} = 2\bar\lambda - \lambda_{n,1}$:
$$\frac{k}{2}\big(\lambda_{n,0} - 2\lambda_{n,1} + \lambda_{n,2}\big) = k\,(\bar\lambda - \lambda_{n,1}) + \frac{k}{2}\,(\lambda_{n,2} - \lambda_{n,1}).$$
The ghost value makes $\lambda_n = \bar\lambda$ midway between segments $0$ and $1$, i.e. at the end crosslinker. In the continuum this is a fixed-value (Dirichlet) condition at both ends:
$$\lambda_n(0, t) = \lambda_n(L_{fil}, t) = \bar\lambda(t).$$

*Relaxation modes.* After a step strain, the deviation from $\bar\lambda$ decays as a sum of sine modes along the filament:
$$\lambda_n(s,t) - \bar\lambda = \sum_{q\ \text{odd}} a_q\,\sin\!\Big(\frac{q\pi s}{L_{fil}}\Big)\,e^{-r_q t},
\qquad r_q = D\Big(\frac{q\pi}{L_{fil}}\Big)^2 = \frac{k}{2}\Big(\frac{q\pi}{N}\Big)^2,
\qquad N = \frac{L_{fil}}{\ell_c}.$$
- The slowest mode ($q = 1$) gives the terminal time $\tau_{term} = 1/r_1 \approx \dfrac{2N^2}{\pi^2 k}$, against $\dfrac{1}{2k}$ for $N = 1$.
- For a uniform initial deviation, mode $q$ carries the weight $\dfrac{8}{\pi^2 q^2}$ (about 81 % in $q = 1$).
- Both agree with the discrete chain: slowest rate $0.0489\,k$ (discrete) against $0.0494\,k$ (continuum) for $N = 10$.

*Why it does not appear in the current model.* There each direction is a single segment ($N = 1$) between two crosslinkers: no interior, no neighbours, no Laplacian. The chain collapses to
$$\frac{d\lambda_n}{dt} = 2k\,(\bar\lambda - \lambda_n) = \Gamma\,(\bar\lambda - \lambda_n),$$
the stretch-free rebinding rule.

**Symmetry:** uniform affine loading keeps the chain symmetric about its midpoint. Only half the chain is needed, with no exchange across the centre ($\lambda_{n,\text{mid}+1} = \lambda_{n,\text{mid}}$ in the stencil).

**From the matrix to the relaxation spectrum** (linear regime, one filament of $N$ segments).

1. *Deviation form.* Because $A\mathbf 1 + \mathbf b = \mathbf 0$, the right-hand side depends only on the lag of the references behind the current stretch, $\mathbf e = \bar\lambda\,\mathbf 1 - \boldsymbol\lambda_n$:
$$A\boldsymbol\lambda_n + \mathbf b\,\bar\lambda = -A\,\mathbf e
\qquad\Rightarrow\qquad
\frac{d\mathbf e}{dt} = k\,A\,\mathbf e + \dot{\bar\lambda}\,\mathbf 1 .$$
A change of stretch creates a uniform lag in every segment; $kA$ relaxes it.

2. *Modes.* $A$ is symmetric, so it has $N$ orthonormal eigenvectors, $A\mathbf v_q = -\mu_q\mathbf v_q$ ($\mu_q > 0$, $\mathbf v_r^T\mathbf v_q = \delta_{rq}$). They form a basis (like $\hat{\mathbf i},\hat{\mathbf j},\hat{\mathbf k}$, but chosen so that $A$ is diagonal in it), so both $\mathbf e$ and the forcing vector $\mathbf 1$ can be expanded in it:
$$\mathbf e(t) = \sum_q a_q(t)\,\mathbf v_q, \quad a_q = \mathbf v_q^T\mathbf e,
\qquad\qquad
\mathbf 1 = \sum_q (\mathbf v_q^T\mathbf 1)\,\mathbf v_q .$$
The $a_q$ are time-dependent amplitudes (components), not eigenvalues. Substituting into $\dot{\mathbf e} = kA\mathbf e + \dot{\bar\lambda}\,\mathbf 1$ and using $A\mathbf v_q = -\mu_q\mathbf v_q$:
$$\sum_q \dot a_q\,\mathbf v_q = \sum_q \Big[-k\mu_q\,a_q + \dot{\bar\lambda}\,(\mathbf v_q^T\mathbf 1)\Big]\mathbf v_q .$$
The $\mathbf v_q$ are independent, so the coefficients of each $\mathbf v_q$ must match (equivalently, multiply by $\mathbf v_r^T$):
$$\dot a_q = -k\mu_q\,a_q + \dot{\bar\lambda}\,(\mathbf v_q^T\mathbf 1).$$
The coupled $N\times N$ system becomes $N$ independent scalar equations. Mode $q$ is a pattern of stored length along the filament ($q = 1$: one half-wave over the whole filament); it relaxes on its own at the rate $k\mu_q$, and the stretch rate feeds it in proportion to $\mathbf v_q^T\mathbf 1$.

   *Example, $N = 2$* (both segments are end segments):
$$A = \begin{pmatrix} -\tfrac32 & \tfrac12\\[2pt] \tfrac12 & -\tfrac32\end{pmatrix}, \qquad
\mathbf v_1 = \tfrac{1}{\sqrt2}\begin{pmatrix}1\\ 1\end{pmatrix}\ (\mu_1 = 1), \qquad
\mathbf v_2 = \tfrac{1}{\sqrt2}\begin{pmatrix}1\\ -1\end{pmatrix}\ (\mu_2 = 2).$$
   - $a_1 = (e_1 + e_2)/\sqrt2$ is the common lag of both segments. It can only relax through the ends, at rate $k$.
   - $a_2 = (e_1 - e_2)/\sqrt2$ is their difference. It also relaxes through the shared junction, so it is faster, at rate $2k$.
   - After a step, $\mathbf e = (\Delta\lambda, \Delta\lambda)$, so $a_1 = \sqrt2\,\Delta\lambda$ and $a_2 = 0$. Only the symmetric mode is excited, with weight $w_1 = (\mathbf 1^T\mathbf v_1)^2/2 = 1$. A uniformly stretched 2-segment filament therefore relaxes at $k$, half the rate of a single segment ($2k$): the middle junction can only redistribute the lag, not release it.

3. *Stress weights.* Question: after a sudden stretch, how does the stress of one filament decay in time? For $N = 1$ (the current model) it is one exponential, $e^{-2kt}$. For several segments it is a sum of exponentials, and the *weights* say what fraction of the stress decays at each rate.

**Lag and stress.** The lag $e_j = \bar\lambda - \lambda_{n,j}$ measures how stretched segment $j$ still is, i.e. how much force it carries. The filament's stress is proportional to the average lag, $\bar e = \tfrac1N\sum_j e_j$. The force of segment $j$ depends on the ratio $\bar\lambda/\lambda_{n,j} = 1 + e_j/\lambda_{n,j} \approx 1 + e_j$ (with $\lambda_{n,j} \approx 1$), so $f_j \approx f_0 + f'_0\,e_j$, with $f'_0$ absorbed into $G_0$. This is an approximation valid only for small amplitudes (errors of order $e_j^2$). For larger deformations it must be corrected: the stress has to use the exact ratio and the full nonlinear force, $\sum_j f(\bar\lambda/\lambda_{n,j})/\lambda_{n,j}$, as the prototype does, and the weights below are then no longer exact.

**Example, $N = 3$.**
```
 free end ── seg 1 ──[J1]── seg 2 ──[J2]── seg 3 ── free end
```
*At $t = 0$* the stretch jumps by $\Delta\lambda$. The references cannot move instantly, so every segment lags by the same amount, $\mathbf e = (1, 1, 1)\,\Delta\lambda$, and the stress jumps to $G_0\Delta\lambda$.

*Afterwards* the rules of §4 decide who can lose lag:
- the end segments (1 and 3) lose lag through the free end at rate $k$ (filament slides in from the unstressed reservoir), and exchange lag with their neighbour;
- the middle segment (2) has no free end: it can only exchange lag with its neighbours, and only when they lag less than it does.

Each event gives one term in the equations of the lags:
- end crosslink unbinds (rate $k$): the end segment's natural stretch jumps to $\bar\lambda$, so its lag is lost: $-k\,e_1$;
- interior crosslink unbinds (rate $k$): the two neighbours average their natural stretches, so both lags become $(e_1 + e_2)/2$, and segment 1 changes by $\tfrac12(e_2 - e_1)$.

$$\frac{de_1}{dt} = -k\,e_1 + \tfrac{k}{2}(e_2 - e_1), \qquad
\frac{de_2}{dt} = \tfrac{k}{2}(e_1 - e_2) + \tfrac{k}{2}(e_3 - e_2), \qquad e_3 = e_1 .$$
These are the rows of $A$ ($-\tfrac32,\ \tfrac12$ for an end segment; $\tfrac12,\ -1,\ \tfrac12$ for the middle one). At $t = 0$, with $\mathbf e = (1, 1, 1)$: $de_1/dt = -k$ and $de_2/dt = 0$. Nothing is exchanged yet: only the end segments relax, the middle one is stuck.

**Where the two rates come from.** A *mode* is a lag pattern that only shrinks and never changes shape: $\mathbf e(t) = c(t)\,\mathbf v$. Inserting this into $\dot{\mathbf e} = kA\mathbf e$ requires $A\mathbf v = -\mu\,\mathbf v$, so **the modes are the eigenvectors of $A$** (the $\mathbf v_q$ of step 2) and **their rates are $k$ times minus the eigenvalues** ($k\mu_q$). For $N = 3$:

| eigenvector of $A$ (mode) | eigenvalue | rate |
|---|---|---|
| $(1, 2, 1)$ | $-\tfrac12$ | $\tfrac12 k$ |
| $(1, 0, -1)$ | $-\tfrac32$ | $1.5\,k$ |
| $(1, -1, 1)$ | $-2$ | $2k$ |

By hand: the two ends are identical, so a mode is symmetric $(a, b, a)$ or antisymmetric $(a, 0, -a)$. For $(a, b, a)$, segment 2 gives $a - b = -\mu b$, i.e. $a = (1-\mu)\,b$, and segment 1 gives $-\tfrac32 a + \tfrac12 b = -\mu a$; together $\mu^2 - \tfrac52\mu + 1 = 0$, so $\mu = \tfrac12$ ($a = \tfrac12 b$: $(1, 2, 1)$) or $\mu = 2$ ($a = -b$: $(1, -1, 1)$). For $(a, 0, -a)$, segment 1 gives $\mu = \tfrac32$. Numerically, `np.linalg.eigh(A)` gives the same pairs (as `_setup_segments` in `material_chain.py` does for every $N$). The scaling of an eigenvector is arbitrary: the integer forms are used here, the normalised ones in step 2; both give the same weights. Check of the two symmetric modes:
- $(1, 2, 1)$ is a mode. Segment 1 loses 1 to the end at rate $k$ and gains $\tfrac12(2 - 1)$ from the middle: net $-\tfrac12 k$. Segment 2 exchanges with both neighbours, $\tfrac12(1-2) + \tfrac12(1-2) = -1$, i.e. $-k = -\tfrac12 k \times 2$. Every segment shrinks by the same fraction, so the pattern decays at rate $\tfrac12 k$.
- $(1, -1, 1)$ is also a mode, at rate $2k$: neighbours differ strongly, so exchange is fast.

**Splitting the step into modes.** The eigenvectors of the symmetric $A$ are orthogonal and form a basis, so any lag vector, at any time, is a combination of them; the coefficients are the amplitudes $a_q$ of step 2. Each coefficient is a projection, $c_q = (\mathbf v_q\cdot\mathbf e)/(\mathbf v_q\cdot\mathbf v_q)$ (with normalised eigenvectors, $a_q = \mathbf v_q^T\mathbf e$). For the uniform lag of the step, $\mathbf e = (1, 1, 1)$:

| mode | $\mathbf v\cdot\mathbf e$ | $\mathbf v\cdot\mathbf v$ | coefficient |
|---|---|---|---|
| $(1, 2, 1)$ | 4 | 6 | $\tfrac23$ |
| $(1, 0, -1)$ | 0 | 2 | $0$ |
| $(1, -1, 1)$ | 1 | 3 | $\tfrac13$ |

$$(1, 1, 1) = \tfrac23\,(1, 2, 1) + 0\cdot(1, 0, -1) + \tfrac13\,(1, -1, 1).$$
The antisymmetric mode is not discarded: its coefficient is zero, because a symmetric vector has no projection on an antisymmetric one. It also stays zero: its forcing $\mathbf v^T\mathbf 1$ in step 2 is zero, so nothing ever feeds it (both ends behave identically, so a uniformly stretched filament stays symmetric), and its average is zero, so the stress would not see it anyway. It would appear only for a loading that is not uniform along the filament (e.g. different rates or initial lags of the two end segments); in the model all segments of a direction share $\bar\lambda$ and $k$, which is why the prototype stores half-chains only.

**From the split to the exponentials** (valid only for a step stretch that is then held, $\dot{\bar\lambda} = 0$ for $t > 0$: the relaxation test that defines $G(t)$; for any other history the forcing term of step 2 stays, and the code integrates the rate equation incrementally with backward Euler instead of using these exponentials).
1. *Modes at every time.* The split holds at every time, with coefficients that change in time; the split of $(1, 1, 1)$ only gives their starting values, $c_1(0) = \tfrac23$, $c_2(0) = 0$, $c_3(0) = \tfrac13$:
$$\mathbf e(t) = c_1(t)\,(1, 2, 1) + c_2(t)\,(1, 0, -1) + c_3(t)\,(1, -1, 1).$$
2. *Substitute.* For $t > 0$ the stretch is held after the step ($\dot{\bar\lambda} = 0$), so $\dot{\mathbf e} = kA\mathbf e$. Each $A\mathbf v$ is a multiple of $\mathbf v$ (eigenvectors): $A(1,2,1) = -\tfrac12(1,2,1)$, $A(1,0,-1) = -\tfrac32(1,0,-1)$, $A(1,-1,1) = -2\,(1,-1,1)$. Hence
$$\dot c_1\,(1, 2, 1) + \dot c_2\,(1, 0, -1) + \dot c_3\,(1, -1, 1) = -\tfrac{k}{2}c_1\,(1, 2, 1) - \tfrac{3k}{2}c_2\,(1, 0, -1) - 2k\,c_3\,(1, -1, 1).$$
3. *Match coefficients.* Both sides use the same three independent vectors, so each coefficient must match: one scalar equation per mode (step 2's $\dot a_q = -k\mu_q a_q$, without forcing since $\dot{\bar\lambda} = 0$):
$$\dot c_1 = -\tfrac{k}{2}\,c_1, \qquad \dot c_2 = -\tfrac{3k}{2}\,c_2, \qquad \dot c_3 = -2k\,c_3 .$$
4. *Solve.* Each amplitude loses a fixed fraction of itself per unit time (as a Maxwell element at constant strain, $\dot\sigma = -\sigma/\tau$):
$$\frac{dc}{dt} = -r\,c \;\Rightarrow\; \frac{dc}{c} = -r\,dt \;\Rightarrow\; \ln c = -rt + \text{const} \;\Rightarrow\; c(t) = c(0)\,e^{-rt},$$
so $c_1 = \tfrac23\,e^{-kt/2}$, $c_2 = 0$, $c_3 = \tfrac13\,e^{-2kt}$, and
$$\mathbf e(t) = \Delta\lambda\Big[\tfrac23\,(1, 2, 1)\,e^{-kt/2} + \tfrac13\,(1, -1, 1)\,e^{-2kt}\Big].$$
Check: the time derivative of $\tfrac23(1,2,1)\,e^{-kt/2}$ is $-\tfrac{k}{2}$ times itself, which is $kA$ times itself because $A(1,2,1) = -\tfrac12(1,2,1)$; each term satisfies $\dot{\mathbf e} = kA\mathbf e$ on its own.

**Lag profile in time.** Written per segment (units of $\Delta\lambda$):
$$e_1 = e_3 = \tfrac23\,e^{-kt/2} + \tfrac13\,e^{-2kt}, \qquad
e_2 = \tfrac43\,e^{-kt/2} - \tfrac13\,e^{-2kt}.$$
Check: at $t = 0$ they give $(1, 1, 1)$, with slopes $-k$ (ends) and $0$ (middle), as above. Example, $kt = 1$ ($e^{-0.5} = 0.607$, $e^{-2} = 0.135$): $e_1 = 0.404 + 0.045 = 0.45$, $e_2 = 0.809 - 0.045 = 0.76$, $\bar e = (0.45 + 0.76 + 0.45)/3 = 0.55$; a single segment ($N = 1$) would be at $e^{-2} = 0.14$. The other rows use the same formulas (a direct numerical integration of the lag equations gives the same values):

| $kt$ | $e_1 = e_3$ (ends) | $e_2$ (middle) | $e_2/e_1$ | $\bar e$ ($\propto$ stress) | $N = 1$ |
|---|---|---|---|---|---|
| 0 | 1.00 | 1.00 | 1.0 | 1.00 | 1.00 |
| 0.5 | 0.64 | 0.92 | 1.4 | 0.73 | 0.37 |
| 1 | 0.45 | 0.76 | 1.7 | 0.55 | 0.14 |
| 2 | 0.25 | 0.48 | 1.9 | 0.33 | 0.02 |
| 4 | 0.09 | 0.18 | 2.0 | 0.12 | 0.00 |

Two phases:
1. *Fast phase* (the $(1, -1, 1)$ part, rate $2k$): the end segments drop quickly while the middle barely moves. The profile changes shape, from flat $(1, 1, 1)$ to peaked $(1, 2, 1)$.
2. *Slow phase* (the $(1, 2, 1)$ part, rate $\tfrac12 k$): once the fast part has gone, the profile stops changing shape and only shrinks ($e_2/e_1 \to 2$). The middle passes lag to the ends, the ends pass it to the reservoir. This bottleneck is slow.

**From lags to stress.** The stress only needs the average lag. The average of $(1, 2, 1)$ is $\tfrac43$, the average of $(1, -1, 1)$ is $\tfrac13$:
$$\bar e(t) = \Delta\lambda\Big[\tfrac23\cdot\tfrac43\,e^{-kt/2} + \tfrac13\cdot\tfrac13\,e^{-2kt}\Big] = \Delta\lambda\Big[\tfrac89\,e^{-kt/2} + \tfrac19\,e^{-2kt}\Big].$$
- 89 % of the stress relaxes slowly (rate $\tfrac12 k$): lag stored in the middle that must drain through the ends.
- 11 % relaxes quickly (rate $2k$): the early drop of the end segments.

The weights sum to 1 because at $t = 0$ the average lag is exactly $\Delta\lambda$, i.e. all of the stress. Mechanically, a 3-segment filament behaves as two Maxwell elements in parallel, $(\tfrac89 G_0,\ 2/k)$ and $(\tfrac19 G_0,\ 1/2k)$, against one element $(G_0,\ 1/2k)$ for a single segment: most of its stress lives 4 times longer. The antisymmetric pattern $(1, 0, -1)$ (rate $1.5k$) is also a mode, but it never appears: a uniform step does not create it, and its average is zero.

**Any $N$.** The same bookkeeping gives
$$G(t) = G_0\sum_q w_q\,e^{-k\mu_q t}, \qquad G^*(\omega) = G_0\sum_q w_q\,\frac{i\omega}{k\mu_q + i\omega},$$
with $w_q$ = (share of the uniform lag in mode $q$) $\times$ (average of mode $q$). With the normalised eigenvectors of step 2 both factors are the same projection, so $w_q = (\mathbf 1^T\mathbf v_q)^2/N$. Longer filaments repeat the same story with more modes: the slowest one (one bump, peaked in the middle) slows down roughly as $1/N^2$ and carries most of the weight ($\to 8/\pi^2 \approx 81\,\%$). A network mixes many filament lengths, so it relaxes with many times. Check: $N = 1$ gives one mode, rate $2k$, $w = 1$: $G(t) = G_0e^{-2kt} = G_0e^{-\Gamma t}$, the current model.

4. *All filament lengths* (§6): each class $N$ contributes its modes with its stress weight $\omega_N$,
$$G(t) = G_0\sum_N \omega_N \sum_q w_{N,q}\,e^{-k\mu_{N,q}t}.$$
This is the spectrum compared with Fig. 5 (`chain_linear.py`). The populations of §8 are bins of these $(N, q)$ modes.

## 6. Filament lengths

Two assumptions:
1. Filament lengths are exponential: $P(L) = e^{-L/\langle L\rangle}/\langle L\rangle$.
2. Crosslinks sit at random along each filament (Poisson process): a filament of length $L$ carries $K$ crosslinks with probability $e^{-x}x^K/K!$, $x = L/\ell_c$.

Averaging 2 over 1, with $m = \langle L\rangle/\ell_c$:
$$P(K) = \int_0^\infty \frac{e^{-x}x^K}{K!}\,\frac{e^{-x/m}}{m}\,dx = \frac1m\Big(1+\frac1m\Big)^{-(K+1)} = \theta_K\,(1-\theta_K)^K,
\qquad \theta_K = \frac{1}{1 + \langle L\rangle/\ell_c},$$
i.e. the number of crosslinks per filament is geometric. The parts beyond the outer crosslinks dangle, so a filament has $N = K - 1$ load-bearing segments.
Each filament contributes in proportion to its number of segments, so the stress weight of a class with $N$ segments is
$$\omega_N = \frac{P(N+1)\,N}{\sum_{N'} P(N'+1)\,N'} .$$
The direction's stress uses $\sum_N \omega_N \langle w'(\lambda_{e,j})/\lambda_{n,j}\rangle_{j \in N}$.

*Only new constant:* $\langle L\rangle$. Linear check: $\langle L\rangle = 10$–20 µm and $k = 0.4$–1.5 s⁻¹ reproduce Fig. 5 (shape error 0.018–0.033, log10 RMS).

## 7. Thermodynamics

**Free energy per direction** (one filament class shown):
$$\Psi = \rho\,n\,\frac{1}{N}\sum_j w(\lambda_{e,j}) \;+\; \Psi_{res}.$$

**Interior exchange: no reservoir is involved.** Use $Y_j = -\partial\Psi/\partial\lambda_{n,j} \propto F_j/\lambda_{n,j}$ as the driving force on the stored length. An exchange moves stored length from $j+1$ to $j$ at rate $J = \frac{k}{2}(\lambda_{n,j+1} - \lambda_{n,j})$, with dissipation
$$\mathcal D_{j,j+1} \propto J\,\Big(\frac{F_j}{\lambda_{n,j}} - \frac{F_{j+1}}{\lambda_{n,j+1}}\Big) \ge 0 .$$
*Why $\ge 0$:* $F(\bar\lambda/\lambda_n)/\lambda_n$ decreases with $\lambda_n$, so the segment with the smaller reference has the larger $F/\lambda_n$, and stored length flows into it. Equivalently, stored length flows from the slack side to the taut side and never the other way. This holds for any state, including slack segments ($F = 0$ on both sides gives $\mathcal D = 0$).

**End exchange: the prestress reservoir.** The network outside the filament holds the tension $F_0$. When the end segment rebinds, it exchanges stored length with this isotensional reservoir. Its energy changes by $d\Psi_{res} = +\rho n\,\frac{1}{N}\,F_0\,d\ln\lambda_{n,1}$ during end events. Then
$$\mathcal D_{end} \propto k\,\big[F(\lambda_{e,1}) - F_0\big]\big(\lambda_{e,1} - 1\big) \ge 0,$$
the same result as the single-segment reservoir. $F$ is non-decreasing, so both factors have the same sign.

*Intuition:* a compressed-but-taut end segment ($\lambda_e^* < \lambda_e < 1$) is re-tensioned to $f_0$; the reservoir pays for that, and the energy balance stays positive.

**Sign note.** The reservoir enters $\Psi$ as $+F_0\ln\lambda_n$, which gives $X = -\partial\Psi/\partial\ln\lambda_n = F - F_0$. The earlier note (`now-i-want-to-purring-hummingbird.md`, reservoir sections) writes $-F_0\ln\lambda_n$ in $\Psi$ and in $\mu^{mech}_b$, which is inconsistent with its own $X = F - F_0$. Only the $+$ sign gives $\mathcal D \ge 0$.

## 8. Cheap version for the UMAT (lumped populations)

The exact chains need about $N/2$ variables per class and many classes, which is fine in Python (the reference) but too many for a UMAT. The linear spectrum of the full mixture is reproduced by $P \approx 6$ parallel populations per direction (max error about 15 % against the full mixture). Each population is the current rule with its own rate and weight:
$$\frac{d\lambda_{n,p}}{dt} = \kappa_p\,k\,\big(\bar\lambda - \lambda_{n,p}\big), \qquad \frac{w'(\lambda_e)}{\lambda_n} \;\to\; \sum_{p=1}^{P}\omega_p\,\frac{w'(\lambda_{e,p})}{\lambda_{n,p}} .$$
- $\kappa_p$, $\omega_p$ are the binned mode rates and weights of the chain mixture. They depend only on $\langle L\rangle/\ell_c$ and are not fitted.
- **What a population is:** a bundle of $(N, q)$ modes of §5 (steps 2–4) whose rates $k\mu_{N,q}$ fall in one bin: $\omega_p = \sum_{\text{bin}}\omega_N w_{N,q}$, $\kappa_p$ = weighted log-mean of $\mu_{N,q}$. Its reference $\lambda_{n,p}$ is the modal coordinate of that bundle, not the reference of a real segment; each bundle follows step 2 with a single rate, which is the current rule with rate $\kappa_p k$.
- **Thermodynamics:** each population is a single-segment reservoir model, so $\mathcal D_p \propto \omega_p\,\kappa_p\,k\,[F(\lambda_{e,p}) - F_0](\lambda_{e,p} - 1) \ge 0$.
- **Caveat:** a population is a lumped relaxation channel, not a piece of filament. It is exact in the linear regime; in the nonlinear regime, using the full force law per population is an assumption. The LAOS check compares it with the exact chains.

## 9. What is fixed and what is open

| Item | Status |
|---|---|
| Interior rule (average of neighbouring references) | exact consequence of force equalisation at conserved stored length; $\lambda_n$ exchange approximates contour exchange (about 3 %, §4) |
| End rule ($\lambda_n \to \bar\lambda$) and $k = \Gamma/2$ | chosen so that $N = 1$ is the current model |
| Rates $k$ | one $k = B/(2c_b)$ per direction, as the current $\Gamma$ (prototype); a per-crosslink $k(f)$ is a refinement |
| $c_b$ | enters only through $\ell_c(c_b)$: force law, $N = L/\ell_c$ (exact chains), $\kappa_p,\ \omega_p$ (populations) |
| Filament lengths | exponential, one constant $\langle L\rangle$ (literature value needed) |
| Open: populations with $c_b$ | freeze $\kappa_p,\ \omega_p$ at $c_{b0}$, or tabulate them in $\langle L\rangle/\ell_c(c_b)$ |
| Open: kinetic consistency | the paired $k_{on}/k_{off}$ (C4) needs $\mu^{mech}_b$ with the reservoir term of the correct sign |
| Open: nonlinear reduction | 6 populations are accurate in the linear regime but not at large amplitude: against the exact chain at $\gamma_A$ = 0.126, $G'$ +8 %, $G''$ +36 % (§11). Next: 8–12 populations or fixed-weight populations against the exact chain; otherwise exact half-chains of a few length classes |

## 10. Python prototype (`material_chain.py`)

Clean copy of `material.py` with `mech_props['relaxation']`:

| Mode | References per direction | Use |
|---|---|---|
| `'sfr'` (default) | 1 | current model, bit-identical to `material.py` |
| `'chain'` | half-chains of the classes $N \in$ `chain_N` (default 1, 2, 3, 4, 6, 8, 11, 16, 22, 32, 45, 64, 90; 154 references at $\langle L\rangle = 21$ µm) | nonlinear reference, expensive |
| `'pops'` | `n_pops` lumped populations (default 6) | UMAT candidate |

- **Inputs:** `L_mean` ($\langle L\rangle$, µm). Classes, populations and weights are fixed at the rest spacing $\ell_c = a\,r_{0f}(c_{b0})$.
- **Update per (sub)step:** backward Euler $(I - k\Delta t\,A)\lambda_n = \lambda_n^{old} + k\Delta t\,b\,\bar\lambda$ for all directions at once, with $k = B/(2c_b^{prev})$ (same lagging as the current rule).
- **Kinetics:** $c_b$ with $\langle k_{off}(f_j)\rangle$ (stress weights) and $\langle k_{off}'\,\partial f_j/\partial c_b\rangle$ in the Newton derivative.
- **Stress:** $w'(\lambda_e)/\lambda_n \to \sum_j \omega_j\,w'(\lambda_{e,j})/\lambda_{n,j}$.

**Checks done:**
- `'sfr'` reproduces `material.py` exactly (one LAOS cycle at $\gamma = 0.126$: difference 0 Pa); `'chain'` with `chain_N = [1]` to $5\cdot10^{-13}$ Pa.
- Vectorised and scalar force solvers agree with `_f_ext` (brentq) to $2\cdot10^{-13}$.
- Linear spectra against the full mixture (max log10 error of $G'$, $G''$ in the Fig. 5 band, $\omega/k = 0.03$–4): chain classes 0.006, 6 populations 0.019.
- Multi-reference code path: 2 and 3 identical populations (rate $2k$, any weights) reproduce `material.py` to $5\cdot10^{-13}$ Pa.
- Full model at $\gamma = 0.5$ % (6 populations, $k = 1.47$ s⁻¹) against a linear twin with the same stepping: $G'$ within 2 % (0.01–0.3 Hz), $G''$ 6–17 % lower. The twin lacks two linear channels of the full model: the prestress geometric stiffness (non-relaxing, about 0.5–0.65 Pa, seen in the single-rate control) and the bound-content kinetics (fast at this rate scaling).

**Inputs from the literature and their consequence:**
- $\langle L\rangle = 21$ µm: the Bausch-lab rigor-HMM preparations at the same actin concentration set the mean filament length to 21 µm with gelsolin (Tharmann 2007; Lieleg 2007, 2009). Unterberger 2013 do not state it.
- With 21 µm, Fig. 5 needs $k \approx 1.47$ s⁻¹ per crosslink (shape error 0.033–0.036). That is about 31 times the current model's rest turnover ($k = \Gamma/2 = 0.047$ s⁻¹) and about 5 times the literature crosslink rates (Lieleg 2008: about 0.3 s⁻¹; Guo: 0.26 s⁻¹). With $k = 0.3$ s⁻¹ the linear dissipation at 0.05 Hz is too low ($\tan\delta$ 0.09 against 0.21).


## 11. LAOS check of the chain model (`laos_chain_check.py`, VM, 2026-10-30)

**Question:** does the strain stiffening of Fig. 4 survive when the reference stretch relaxes through the chain mechanism at the rate Fig. 5 needs?

**Setup:** Fig. 4 parameter set (`c5_w15_plate_guo_dxs`), **not refitted**, point stress, dt 0.1 s, factor 3, 5 cycles from rest. The rate is raised by scaling $k_{off,0}$, $k_{catch,0}$, $k_{on,0}$ together ($k_{eq}$ fixed), so the rest state is unchanged. Cases: `ref` (current SFR), `sfr_k035` (control: SFR with the rate ×7.4), `pops7` (6 populations, $\langle L\rangle$ = 7 µm, ×7.4, $k \approx 0.35$ s⁻¹), `pops15` (15 µm, ×18.5, $k \approx 0.87$ s⁻¹), `chain7` (exact half-chains, 7 µm, ×7.4; at 0.0126 and 0.126 only). Point stress over-stiffens all cases at 0.126, so the cases are compared with each other rather than with the data.

**1. The decisive test passes.** A faster single rate makes the network almost fluid; the chain at the same speed keeps stiffness and stiffening:

| | ref | single rate ×7.4 | chain 7 µm | chain 15 µm | data |
|---|---|---|---|---|---|
| $G'_{lin}$ (Pa) | 10.66 | **2.47** | 9.78 | 9.90 | 13.55 |
| $G'(0.126)/G'_{lin}$ | 4.51 | **1.41** | 4.94 | 4.98 | 2.19 |

Unlike the Broedersz-type modes in series (which killed the stiffening), parallel references only shift the zero-force point of the unchanged compound law.

**2. Improvements at large amplitude (no refit):**

| | ref | chain 7 µm | chain 15 µm | data |
|---|---|---|---|---|
| $G''(0.126)$ (Pa) | 15.1 | 9.9 | 8.3 | 6.6 |
| $G''(0.158)$ (Pa) | 41.2 | **15.4** | **14.3** | 11.7 |
| $G'(0.158)$ (Pa) | 65.4 | 51.4 | 48.8 | 47.8 |
| $G_M$ at 0.126 (Pa) | 12.6 | 11.6 | 11.4 | 11.1 |

The over-dissipation at 0.158 (3.5× the data with the current model) drops to about 1.3×.

**3. Open problems:**
1. *Linear dissipation lower:* $G''_{lin}$ 2.29 Pa (7 µm), 1.79 Pa (15 µm), against 2.90 (ref) and 4.44 (data). The rate came from the shape of Fig. 5; the Fig. 4 sample is more dissipative at 0.05 Hz ($\tan\delta \approx 0.33$ against about 0.21 in Fig. 5). A refit on Fig. 4 will move $k$ (the Fig. 4 / Fig. 5 sample difference again).
2. *Stiffening too strong:* $G'(0.1)/G'_{lin}$ 2.9 against 2.3 (ref) and 1.6 (data); refit target ($s_{taut}$, $\lambda_0$), also reduced by plate averaging.
3. *Six populations are not accurate enough in the nonlinear regime* (`pops7` against the exact `chain7`): linear $G'$ +0.3 %, $G''$ +4 % (fine); at 0.126 $G'$ +8 %, **$G''$ +36 %**, $\sigma_{max}$ +9 %. Expected from lumping (a nonlinear force of one lumped reference is not the sum over the segments it replaces), but larger than acceptable. The exact chain's loop is also the closest to the data (loop RMS 0.97 Pa against 1.14 for ref and 1.20 for `pops7`).
4. *Not steady at 0.158 after 5 cycles:* cycle-to-cycle change 0.13–0.18 Pa (ref: 0.002). The slow populations (minutes) have not settled; a fit needs more cycles or a fixed-cycle convention.
5. *Large swings of the bound content* (last cycle): $c_b/c_{b0}$ = 0.55–1.30 at 0.126 and 0.25–1.38 at 0.158, well beyond the ±10 % where frozen populations are accurate. Fixed-weight populations (rates following $c_b$) were accurate up to about −23 % / +27 %: to be checked at 0.126, out of range at 0.158.
6. *Cost per run:* ref 0.07 h; `pops7` 1.2–1.7 h (17–25×); `pops15` 2.2–2.9 h; `chain7` 3.1–5.6 h. A refit with populations at these settings takes days on the VM.

**Verdict:** valid as a mechanism (it fixes stiffness and stiffening under fast relaxation, and the large-amplitude over-dissipation); not yet fit-ready.

**Next steps, in order:**
1. *Nonlinear reduction:* 0.126 with 8, 10, 12 populations and with fixed-weight populations, against `chain7` (about 2 h each on the VM). If none gets $G''$ within about 10 %, implement exact half-chains of a few length classes instead (more expensive, exact).
2. *Speed-up:* profile the populations step; most of the cost is probably the closed-system mass balance re-running the kinetics, which scales with the faster rate.
3. *Refit on Fig. 4* with the chain ($s_{taut}$, $\lambda_0$, rate scale at $\langle L\rangle$ = 7 µm), plate-averaged, with enough cycles; this is where the Fig. 4 / Fig. 5 trade-off shows.

---------------------------------

Read the formula term by term and ask which factors depend on the segment. For one direction:
$$\sigma_{\text{dir}} = J^{-1},n,w',\lambda^{-1},\mathbf m\otimes\mathbf m .$$

Step 1: only w′ knows about the segments

- J, λ and m are kinematics. All segments in a direction share them, because the junctions are pinned to the affine network.
- n and r₀ depend on c_b of that direction, which is also shared.
- The λn,j (and so the lags) appear only inside w′.

So at a given deformation, σ_dir is a fixed prefactor times w′. Everything that relaxes is in w′.

Step 2: w′ is an average over segments

$$w' = r_0\lambda_0\sum_j w_j,g(e_j), \qquad g(e_j) = \frac{f(\lambda/\lambda_{n,j})}{\lambda_{n,j}}, \quad \lambda_{n,j} = \lambda - e_j, \quad \sum_j w_j = 1 .$$
w′ is r₀λ₀ times the average of one function g, evaluated at each segment's lag.

Step 3: linearise g, then average

For small lags, expand g about the relaxed state e_j = 0 (that is, λn,j = λ):
$$g(e_j) \approx g(0) + g'(0),e_j ,\qquad g(0) = \frac{f_0}{\lambda}, \qquad g'(0) = \frac{f'_0 + f_0}{\lambda^2}.$$
The average of a linear function is the same linear function of the average, so
$$w' \approx r_0\lambda_0\Big[g(0) + g'(0),\bar e\Big], \qquad \bar e = \sum_j w_j e_j .$$

Result

$$\sigma_{\text{dir}} \approx \underbrace{J^{-1} n, r_0\lambda_0, g(0),\lambda^{-1}\mathbf m\otimes\mathbf m}{\text{relaxed state: no lag}} ;+; \underbrace{J^{-1} n, r_0\lambda_0, g'(0),\lambda^{-1}\mathbf m\otimes\mathbf m}{\text{fixed prefactor}};\bar e .$$
- The first term is the stress when every segment has caught up, λn,j = λ. It carries the prestress force f₀, does not relax, and is the source of G_geo.
- The second term is the part that relaxes. It is a prefactor (kinematics, n, f′₀ + f₀) times ē. This is the meaning of "proportional to the average lag": its prefactor, linearised at rest, is G₀.