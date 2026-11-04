# Abaqus UEL Implementation Breakdown: `5_CL_DIFF/src/element/`

This directory contains the Fortran source files responsible for the **macroscopic element formulation** (the "User Element" or UEL). While the `material/` directory handles the local physics at a single mathematical point, the `element/` directory handles the actual 3D geometry, shape functions, integration over the element volume, and communication with the global Abaqus solver.

Here is a breakdown of the purpose of each file and their exact sequence of execution.

---

## 1. File Breakdown

### `uel.f90`
*   **Purpose:** This is the main entry point. Abaqus always looks for a subroutine named `UEL` to execute user-defined elements. 
*   **Action:** It performs global safety checks (verifying the procedure is a coupled analysis, ensuring large deformations `nlgeom` is on, etc.) and routes the execution to the correct dimensionality-specific subroutine (e.g., calling the 3D element or the 2D element).

**Available Subroutines:**
1. `UEL(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, ...)`
    *   **Functionality:** This is the standard interface required by the Abaqus solver. It acts strictly as a logical gateway and error handler before executing any finite element math.
    *   **Procedure Verification:** Validates that the active Abaqus step is compatible with a fully coupled, large-deformation analysis by checking the internal `LFLAGS` array:
        *   $\text{LFLAGS}(1) \in \{64, 65, 72, 73\}$: Ensures the step type supports coupled chemo-mechanics (e.g., mass diffusion or pore-fluid flow coupled with stress).
        *   $\text{LFLAGS}(2) = 1$: Ensures non-linear geometric kinematics (`nlgeom=yes`) are active.
        *   $\text{LFLAGS}(4) = 0$: Ensures the current increment is a general step, not a linear perturbation.
    *   **Element Routing:** Reads the user-defined element integer (`JTYPE`). If `JTYPE = 3` (representing a 3D block element), it hands over all raw Abaqus arrays to the `U3D8` subroutine for actual element integration.

### `u3d8.f90`
*   **Purpose:** The core engine of the 8-node, 3D hexahedral element (a coupled displacement/chemical potential element).
*   **Action:** This massive subroutine drives the actual finite element math. It loops over the nodes to get nodal displacements and chemical potentials. It then loops over the 8 Gauss integration points, calculating the Deformation Gradient ($\mathbf{F}$) and chemical gradients. It calls the `MATERIAL` subroutine to get the local stresses, and then numerically integrates those stresses over the volume to construct the element residual vectors (`Ru`, `Rc`) and the partitioned Jacobian stiffness matrices (`Kuu`, `Kcc`, `Kuc`, `Kcu`).

### `elmroutines.f90`
*   **Purpose:** A library of geometric and mathematical helper functions.
*   **Action:** Contains the exact definitions for the Gauss integration points and weights (e.g., `xint3D8pt`), the isoparametric shape functions (e.g., `calcShape3DLinear`), their spatial derivatives, and the mappings from the local "master" element coordinates to the global real-world coordinates.

**Available Subroutines:**
1. **Gauss Quadrature Routines** (e.g., `xint3D8pt`, `xint3D1pt`, `xint2D4pt`)
    *   **Functionality:** Defines the standard optimal coordinates and weights for Gaussian numerical volume integration.
    *   **Mathematical Action:** For an 8-point fully integrated 3D element, it sets the local master element coordinates $(\xi, \eta, \zeta)$ to $\pm 1/\sqrt{3}$ and assigns a uniform weight of $w_i = 1.0$ for all 8 points.
2. **Surface Quadrature Routines** (e.g., `xintSurf3D4pt`, `xintSurf3D1pt`)
    *   **Functionality:** Similar to the volume integration routines, but defines the 2D coordinates and weights for numerical integration specifically over the 6 exterior faces of the element (used for applying boundary fluid fluxes).
    *   **Mathematical Action:** When Abaqus signals a flux on a specific face (e.g., Face 3), this routine logically locks the local coordinate perpendicular to that face (e.g., freezing $\eta = -1$) and distributes 4 optimal Gauss points across the remaining 2D plane (e.g., varying $\xi$ and $\zeta$ at $\pm 1/\sqrt{3}$). It assigns a surface integration weight of $w_i = 1.0$ to each of these 4 points.
3. **Shape Function Routines** (`calcShape3DLinear`, `calcShape2DLinear`)
    *   **Functionality:** Evaluates the isoparametric shape functions and their local derivatives at a given Gauss integration point.
    *   **Action for 3D Elements:** For an 8-node brick element, it mathematically evaluates the trilinear shape functions $N^k$ for node $k$:
        $$ N^k(\xi, \eta, \zeta) = \frac{1}{8} (1 \pm \xi)(1 \pm \eta)(1 \pm \zeta) $$
    *   It also evaluates the first-order local derivatives (e.g., $\frac{\partial N^k}{\partial \xi}$, $\frac{\partial N^k}{\partial \eta}$), as well as all 9 cross-terms of the second-order derivatives $\frac{\partial^2 N^k}{\partial \xi_i \partial \xi_j}$.
4. **Surface Mapping Routines** (`computeSurf3D`, `computeSurf`)
    *   **Functionality:** Computes the shape functions and the mapping Jacobian strictly for the 2D exterior faces of the element.
    *   **Mathematical Action:** Used to calculate the true physical differential surface area $dA$ corresponding to the master element face. For example, if a boundary flux is applied on a face where $\zeta$ is constant, the physical area Jacobian is:
        $$ dA = \sqrt{ \left( \frac{\partial y}{\partial \xi} \frac{\partial z}{\partial \eta} - \frac{\partial y}{\partial \eta} \frac{\partial z}{\partial \xi} \right)^2 + \left( \frac{\partial x}{\partial \xi} \frac{\partial z}{\partial \eta} - \frac{\partial x}{\partial \eta} \frac{\partial z}{\partial \xi} \right)^2 + \left( \frac{\partial x}{\partial \xi} \frac{\partial y}{\partial \eta} - \frac{\partial x}{\partial \eta} \frac{\partial y}{\partial \xi} \right)^2 } $$
        This mapping ensures that any user-defined chemical or fluid fluxes (like `*Dsload` commands) applied to the element "skin" in Abaqus are correctly integrated into the global equilibrium equations.

### `assembleelement.f90`
*   **Purpose:** Matrix stitching.
*   **Action:** Abaqus requires a single, monolithic stiffness matrix (`amatrx`) and a single residual vector (`rhs`). However, `u3d8.f90` calculates the mechanical components (`Kuu`, `Ru`) and chemical components (`Kcc`, `Rc`) separately. This subroutine systematically stitches those sub-matrices together into the final monolithic arrays, ensuring the mechanical DOFs (1, 2, 3) and the chemical DOF (4) are in the exact order Abaqus expects.

**Available Subroutines:**
1. `AssembleElement(nDim, nNode, ndofel, Ru, Rc, Kuu, Kuc, Kcu, Kcc, rhs, amatrx)`
    *   **Functionality:** Maps the partitioned local element residuals and stiffness matrices into the global monolithic arrays required by Abaqus, based on the dimensionality (`nDim`).
    *   **Residual Assembly:** For a 3D element ($nDim = 3$), the global residual vector $\{ \text{rhs} \}$ for a given node $i$ is assembled by stacking the displacement residual $\mathbf{R}_u$ and the chemical potential residual $R_c$:
        $$ \{ \text{rhs} \}_{I} = \begin{Bmatrix} R_u^{x}(i) \\ R_u^{y}(i) \\ R_u^{z}(i) \\ R_c(i) \end{Bmatrix} $$
        where the global row index $I$ spans from $4(i-1)+1$ to $4(i-1)+4$.
    *   **Tangent Matrix Assembly:** The global coupled tangent stiffness matrix $[ \text{amatrx} ]$ between node $i$ and node $j$ is assembled by inserting the 4 decoupled sub-blocks into their correct row/column positions:
        $$ [ \text{amatrx} ]_{I,J} = \begin{bmatrix} \mathbf{K}_{uu}(i,j) & \mathbf{K}_{uc}(i,j) \\ \mathbf{K}_{cu}(i,j) & K_{cc}(i,j) \end{bmatrix} $$
        where $\mathbf{K}_{uu}$ is a $3 \times 3$ mechanical matrix, $\mathbf{K}_{uc}$ is a $3 \times 1$ chemomechanical coupling vector, $\mathbf{K}_{cu}$ is a $1 \times 3$ mechanochemical coupling vector, and $K_{cc}$ is the scalar chemical stiffness.

### `uvarm.f90`
*   **Purpose:** Visualization and Post-Processing.
*   **Action:** Abaqus Viewer cannot natively plot internal variables from a `UEL`. To fix this, your model uses a "dummy mesh" of standard Abaqus elements overlaid on top of the UELs. Abaqus calls `UVARM` for the dummy mesh, and this subroutine simply reads the saved state variables (like $c_b$, $c_f$, $J$) out of the `globalSdv` array and hands them to the viewer.

**Available Subroutines:**
1. `UVARM(UVAR, DIRECT, T, TIME, DTIME, ...)`
    *   **Functionality:** Solves the "black box" problem of Abaqus UELs by bridging the data computed inside the hidden UEL directly to the visible dummy element mesh for post-processing.
    *   **Data Mapping:** It loops through every single state variable index $i \in \{1 \dots \text{nsdv}\}$ at the current integration point (`NPT`) and populates the visualization array `UVAR` by fetching the data from the shared Fortran `globalSdv` module array:
        $$ \text{UVAR}_i = \text{globalSdv}(\text{NOEL} - \text{ElemOffset}, \text{NPT}, i) $$
        The parameter `ElemOffset` acts as a numerical shift to perfectly map the active dummy element ID (`NOEL`) back to the ID of the mathematical User Element that actually computed the mechanics.

---

## 2. Sequence of Execution (The Call Order)

During a single global Newton-Raphson iteration, the execution flows exactly like this:

1.  **Abaqus Core Solver** calls $\rightarrow$ `UEL` (`uel.f90`)
2.  `UEL` verifies the setup and calls $\rightarrow$ `U3D8` (`u3d8.f90`)
3.  `U3D8` calls $\rightarrow$ `elmroutines.f90` to get the Gauss points, shape functions, and spatial gradients.
4.  `U3D8` loops over the 8 integration points. Inside the loop, it calls $\rightarrow$ `MATERIAL` (located in the `material_AFFCL/_umat_.f90` file) to compute the physical stress and chemical flux at that specific point.
5.  After the loop finishes, `U3D8` has the full element stiffness sub-matrices. It calls $\rightarrow$ `AssembleElement` (`assembleelement.f90`).
6.  `AssembleElement` stitches the matrices and passes `amatrx` and `rhs` back to `U3D8`, which passes them back to `UEL`, which returns them to the **Abaqus Core Solver**.
7.  *...Abaqus solves the global system...*
8.  At the end of the increment, Abaqus calls $\rightarrow$ `UVARM` (`uvarm.f90`) to pull the data for the visualization dummy mesh.
