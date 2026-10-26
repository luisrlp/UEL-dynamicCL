# 1) Adapt notebook kinetics.ipynb to allow for incremental applications of the
# deformation gradient with time being an important parameter as well
# 2) Note that an open system is being simulated: we prescribe the chemical potential!
# 3) Ensure correct enforcement of uniaxial, biaxial and shear 
# boundary conditions (might be done just by summing the stress contributions from 
# directions where we know the stress must be null)

# 4) Add possibility for cyclic stretching and stretch-relaxation

# OUTSIDE THE CLASS
# def solve_uniaxial_step(material, lxambda_1, mu, dt):
#     """Uses Scipy to find transverse stretch and returns final stress"""
#     def objective(lambda_trans):
#         F_trial = np.diag([lambda_1, lambda_trans[0], lambda_trans[0]])
#         # Notice how we don't have to pass history variables!
#         stress = material.evaluate_state(F_trial, mu, dt) 
#         return stress[1, 1]  # Enforce sigma_22 = 0
        
#     # Run scipy root finder here...
#     return converged_stress

from matplotlib.lines import Line2D
import numpy as np
import math                                                                         
from scipy.optimize import brentq, root_scalar
from dfd import generate_icosahedron, apply_dispersion
from material import DynamicCrosslinking
from plotting import plot_microsphere
import matplotlib.pyplot as plt
import contextlib, os
from concurrent.futures import ProcessPoolExecutor
import sys


tex_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 12,
    "font.size": 12,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10
}
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams.update(tex_fonts)


# Deformation profiles

def profile_linear(t, rate, deformation_type: str):                                           
    if deformation_type in ['uniaxial', 'biaxial']:
        return 1.0 + rate * t
    elif deformation_type == 'shear':
        return rate * t
                                                                        
def profile_step(t, t_step, target_stretch, deformation_type: str):
    if deformation_type in ['uniaxial', 'biaxial']:
        return target_stretch if t >= t_step else 1.0                      
    elif deformation_type == 'shear':
        return target_stretch if t >= t_step else 0.0
                                                                        
def profile_cyclic(t, amplitude, frequency, deformation_type: str):
    if deformation_type in ['uniaxial', 'biaxial']:                          
        return 1.0 + amplitude * np.sin(2 * np.pi * frequency * t)
    elif deformation_type == 'shear':
        return amplitude * np.sin(2 * np.pi * frequency * t)

# Boundary condition solvers

def solve_uniaxial(material, lambda_1, mu, dt):                        
    """                                                                
    Enforces sigma_22 = sigma_33 = 0 by finding the correct transverse stretch.                                                                 
    """                                                                
        
    # Save history at the start of the time step           
    saved_cb_history = material.cb_history.copy()
    saved_cf_history = material.cf_history
    saved_kon_history = material.kon_history.copy()
    saved_koff_history = material.koff_history.copy()
    saved_lamn_history = material.lamn_history.copy()
    
    def objective(lambda_trans):
        # Trial deformation gradient
        F_trial = np.diag([lambda_1, lambda_trans, lambda_trans])      
        
        # Reset history to time `t` before the trial evaluation        
        material.cb_history = saved_cb_history.copy()
        material.cf_history = saved_cf_history
        material.kon_history = saved_kon_history.copy()
        material.koff_history = saved_koff_history.copy()
        material.lamn_history = saved_lamn_history.copy()
        
        # Evaluate stress
        stress = material.evaluate_state(F_trial, mu, dt)
        
        # We want the transverse stress sigma_22 to be exactly 0       
        return stress[1, 1]
        
    # Initial guess for transverse    
    guess = 1.0 / np.sqrt(lambda_1)
    try:
        sol = root_scalar(objective, x0=guess, x1=guess*1.01, method='secant', xtol=1e-9)
        lambda_trans_converged = sol.root
    except ValueError:
        print(f"Warning: Uniaxial solver failed to find a root within bounds at lambda_1 = {lambda_1}!")
        lambda_trans_converged = 1.0 / np.sqrt(lambda_1)
        
    F_converged = np.diag([lambda_1, lambda_trans_converged, lambda_trans_converged])
        
    # Reset history 
    material.cb_history = saved_cb_history.copy()
    material.cf_history = saved_cf_history
    material.kon_history = saved_kon_history.copy()
    material.koff_history = saved_koff_history.copy()
    material.lamn_history = saved_lamn_history.copy()
    final_stress = material.evaluate_state(F_converged, mu, dt)

    # print(f"Uniaxial: lambda_1 = {lambda_1:.4f}, lambda_trans = {lambda_trans_converged:.4f}, sigma_11 = {final_stress[0, 0]:.4f}, sigma_22 = {final_stress[1, 1]:.4f}, sigma_33 = {final_stress[2, 2]:.4f}")
    return final_stress, F_converged, material.cb_history, material.cf_history, material.kon_history, material.koff_history, material.Rbind_history, material.f_history
                                                                        

def solve_biaxial(material, lambda_1, lambda_2, mu, dt):
    """
    Enforces sigma_33 = 0 by finding the correct out-of-plane stretch (lambda_3).
    """
    # Save history at the start of the time step           
    saved_cb_history = material.cb_history.copy()
    saved_cf_history = material.cf_history
    saved_kon_history = material.kon_history.copy()
    saved_koff_history = material.koff_history.copy()
    saved_lamn_history = material.lamn_history.copy()
    
    def objective(lambda_3):
        # Trial deformation gradient
        F_trial = np.diag([lambda_1, lambda_2, lambda_3])      
        
        # Reset history to time `t` before the trial evaluation        
        material.cb_history = saved_cb_history.copy()
        material.cf_history = saved_cf_history
        material.kon_history = saved_kon_history.copy()
        material.koff_history = saved_koff_history.copy()
        material.lamn_history = saved_lamn_history.copy()
        
        # Evaluate stress
        stress = material.evaluate_state(F_trial, mu, dt)
        
        # We want the out-of-plane stress sigma_33 to be exactly 0       
        return stress[2, 2]

    guess = 1.0 / (lambda_1 * lambda_2) 
    try:
        sol = root_scalar(objective, x0=guess, x1=guess*1.01, method='secant', xtol=1e-6)
        lambda_3_converged = sol.root
    except ValueError:
        print(f"Warning: Biaxial solver failed to find a root within bounds at lambda_1 = {lambda_1}!")
        lambda_3_converged = 1.0 / (lambda_1 * lambda_2)
            
    F_converged = np.diag([lambda_1, lambda_2, lambda_3_converged])
        
    # Reset history 
    material.cb_history = saved_cb_history.copy()
    material.cf_history = saved_cf_history
    material.kon_history = saved_kon_history.copy()
    material.koff_history = saved_koff_history.copy()
    material.lamn_history = saved_lamn_history.copy()
    final_stress = material.evaluate_state(F_converged, mu, dt)      
    
    return final_stress, F_converged, material.cb_history, material.cf_history, material.kon_history, material.koff_history, material.Rbind_history, material.f_history
                                                                        
def solve_shear(material, gamma, mu, dt):
    """
    Evaluates Simple Shear. 
    Because simple shear is fully constrained (constant volume, J=1), 
    we do not need a Scipy solver to find unknown stretches.
    """
    # Construct the simple shear deformation gradient
    # print(f"Evaluating simple shear with gamma = {gamma:.4f}")
    F_converged = np.array([
        [1.0, gamma, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0]
    ])
    
    # Evaluate the material (history updates automatically)
    final_stress = material.evaluate_state(F_converged, mu, dt)
    
    return final_stress, F_converged, material.cb_history, material.cf_history, material.kon_history, material.koff_history, material.Rbind_history, material.f_history

# Main function
def run_simulation(material, time_array, profile_func, mode, mu):      
    history = {'time': [], 'def_val': [], 'stress': [], 'F': [], 'cb': [], 'cf': [], 'kon': [], 'koff': [], 'Rbind': [], 'f': [],'lamn':[],
               # Local kinetics diagnostics of the accepted evaluation of each time step
               'nsub': [], 'min_dHdcb': [], 'max_dcb_rel': [], 'kin_converged': []}

    for i in range(len(time_array)):                                     
        if i == 0:                                                       
            dt = 0.0                                                     
        else:                                                            
            dt = time_array[i] - time_array[i-1]                         
                                                                             
        t = time_array[i]  
        
        # Prescribed stretch/shear for time = t     
        def_val = profile_func(t)
        print(f"----------Time: {t:.2f}, Prescribed deformation: {def_val:.4f}----------")
        # Enforce BCs and evaluate material
        if mode == 'uniaxial':
            stress, F, cb_tot, cf, kon, koff, Rbind, f = solve_uniaxial(material, def_val, mu, dt)  
        elif mode == 'biaxial':
            stress, F, cb_tot, cf, kon, koff, Rbind, f = solve_biaxial(material, def_val, def_val, mu, dt) # Equibiaxial
        elif mode == 'shear':
            stress, F, cb_tot, cf, kon, koff, Rbind, f = solve_shear(material, def_val, mu, dt)     
            
        # 3. Log results
        history['time'].append(t)
        history['def_val'].append(def_val)
        history['stress'].append(stress)
        history['F'].append(F)
        history['cb'].append(cb_tot)
        history['cf'].append(cf)
        history['kon'].append(kon)
        history['koff'].append(koff)
        history['Rbind'].append(Rbind)
        history['f'].append(f)
        history['lamn'].append(material.lamn_history.copy())
        history['nsub'].append(material.nsub)                        # kinetics sub-steps used
        history['min_dHdcb'].append(float(np.min(material.dHdcb)))   # min dH/dcb over directions
        history['max_dcb_rel'].append(float(np.max(material.dcb_rel)))  # max relative cb change per sub-step
        history['kin_converged'].append(bool(np.all(material.kin_converged)))  # all directions converged
    return history

def moduli(t, gamma, sigma, frequency, both_halves=False):
    """ 
    G' and G'' from Unterberger et al. 2013:
    G' = (sigma_A/gamma_A) * cos(delta)
    G'' = (sigma_A/gamma_A) * sin(delta)
    Positive peaks are used, with the option for averaging with the negative half cycle.
    """
    def peak(i, y):
        '''Returns the time and value of the peak using parabolic interpolation'''
        y0, y1, y2 = y[i - 1], y[i], y[i + 1]
        # Parabola: y(s) = a s² + b s + y1
        # b = (y2 − y0)/2, a = (y0 − 2 y1 + y2)/2 = den/2
        den = y0 - 2.0 * y1 + y2
        s = 0.5 * (y0 - y2) / den if den != 0.0 else 0.0
        return t[i] + s * (t[i + 1] - t[i]), y1 - 0.25 * (y0 - y2) * s
    out = []
    for sgn in ((1.0, -1.0) if both_halves else (1.0,)):
        # Gamma peaks
        tg, gA = peak(np.argmax(sgn * gamma[1:-1]) + 1, sgn * gamma)
        # Sigma peaks
        ts, sA = peak(np.argmax(sgn * sigma[1:-1]) + 1, sgn * sigma)
        out.append((gA, sA, 2.0 * np.pi * frequency * (tg - ts)))
    gA, sA, delta = np.mean(out, axis=0)
    return sA / gA * np.cos(delta), sA / gA * np.sin(delta), delta

def sweep_worker(args):
    """One amplitude of the sweep from rest (fresh material): gamma_A, G', G'', delta of the last cycle."""
    gA, diffusion_props, mech_props, sim_props, frequency, n_cycles, dt, mu = args
    t = np.arange(0.0, n_cycles / frequency + 0.5 * dt, dt)
    with open(os.devnull, 'w') as devnull, contextlib.redirect_stdout(devnull):     # silence the per-step prints
        mat = DynamicCrosslinking(diffusion_props, mech_props, sim_props)
        prof = lambda tt: profile_cyclic(tt, amplitude=gA, frequency=frequency,
                                         deformation_type=sim_props['deformation_type'])
        res = run_simulation(mat, t, prof, mode=sim_props['deformation_type'], mu=mu)
    th = np.array(res['time'])
    gh = np.array(res['def_val'])
    sh = np.array([s[0, 1] for s in res['stress']])
    m = th >= th[-1] - 1.0 / frequency - 0.5 * dt                                    # last full cycle
    G1, G2, delta = moduli(th[m], gh[m], sh[m], frequency)
    return gA, G1, G2, delta
    

if __name__ == "__main__":

    RUN = sys.argv[1] if len(sys.argv) > 1 else None

    # --- MATERIAL PARAMETERS ---
    THETA = 294. # 21 + 273.0  # Temperature (K)
    diffusion_props = {
        'k_off_0' : 0.17686159107698207*(1-0.998999999440339),    # 1/s (baseline off-rate) [0.05]
        'k_catch_0' : 0.17686159107698207*0.998999999440339,  # 1/s (baseline catch-rate) [0.0]
        'k_on_0' : 30.9044970*0.17686159107698207,   # 1/s (baseline on-rate) [0.0125]
        'CHI' : 0.1,     # Flory-Huggins parameter
        'THETA' : THETA,  # Temperature (K)
        'dx_kT' : 0.00050198 / (1.38065e-5 * (THETA)),       # Reactive distance for slip bonds over thermal energy: Delta x / (k_B T)
        'dx_c_kT': 0.009333 / (1.38065e-5 * (THETA)),  # Reactive distance for catch bonds over thermal energy: Delta x_c / (k_B T)
        'D' : 0.1,
        'MU0' : 0.0,
        'VMOL' : 0.15,
        'RGAS' : 8.31446261815324,  # J/(mol*K)
        # Concentration constraints
        'CACTIN' : 0.0095,     # Reference actin concentration (mM)
        'R' : 0.1,      # Initial total crosslinker ratio (cb + cf)
        'RBMAX' : 0.25,     # Max bound ratio
        'RFMAX' : 1.,       # Max free ratio (dummy large)
    }

    mech_props = {
        'filament': 'ext',  # 'inext' for inextensible, 'ext' for extensible
        'KBULK': 1.e3,
        'LAMBDA0': 1.0113904611496767,
        'MU0STRETCH': 38600.0, # 38600
        'Lp': 16.0,
        'BETA': 0.438, # 0.5,
        'R0C': 0.014,
        'ETA': 1., # 2./3.,
        'a_ratio': 1.0113904611496767*1.080591339670614, # 1.2,  # L / r0f ratio
        'b_disp': 0.000001, # Dispersion parameter
        'prefdir': np.array([1.0, 0.0, 0.0]),  # Preferred direction for anisotropy
        'stretch_free_binding': True, # True: reference state updated each increment
        'xi_t': 1.0, # fraction of the network stress that binds un-stretched (used only if stretch_free_binding=True)
        }

    ### Unterberger props
    ## Direct
    # T = 294 K
    # Lp = 16 microm
    # beta = 0.438
    # mu0 = 38600
    # lambda0 = 1.012
    ## Indirect
    # p = 13.12
    # r0 (= r0f) = 1 microm
    # L = 1.048 microm
    # n = 14.31
    ## Visco
    # tau1 = 2 s
    # theta1 = 0.835
    # tau2 = 0.007 s
    # theta2 = 1.6
    # tau3 = 15 s
    # theta3 = 1.4


    # --- SIMULATION PARAMETERS ---
    sim_props = {
        'deformation_type': 'shear',  # Choose from ['uniaxial', 'biaxial', 'shear']
        'def_init': 0.0,  # Initial deformation magnitude
        'def_max': 0.126,  # Maximum deformation magnitude
        't_init': 0.0,  # Initial time
        't_max': 100.0,  # Maximum time
        'delta_t': 0.05,  # Time step
        'factor': 3,  # Sphere discretization factor
        'max_iter': 50,  # Maximum iterations for solvers
        'tol': 1.0e-12  # Tolerance for solvers
    }

    # ----------- LAOS (Fig. 4a) -----------
    run_laos = False
    mu_fixed = None # material.INITMU  # Keeping mu constant for now          
    frequency = 0.05
    tag = '0'                                                   # prefix of the saved figures and CSV
    if RUN is not None:
        from amplitude_sweep_runs import load_run
        diffusion_props, mech_props, cfg = load_run(RUN)        # optimum + fixed constants of that run
        sim_props.update(def_max=cfg['gamma_A'], delta_t=cfg['dt'], factor=cfg['factor'],
                         t_max=cfg['cycles'] / frequency)       # same protocol as the fit (cycle 3)
        tag = RUN
    if run_laos:
        factor = sim_props['factor']
        material = DynamicCrosslinking(diffusion_props, 
                                    mech_props,
                                    sim_props)
        # Time array
        t_array = np.arange(sim_props['t_init'], sim_props['t_max'] + sim_props['delta_t'], sim_props['delta_t'])
        # Define stretch rate and profile     
        # stretch_rate = (sim_props['def_max'] - sim_props['def_init']) / sim_props['t_max']
        my_profile = lambda t: profile_cyclic(t, 
                                            amplitude=sim_props['def_max'], 
                                            frequency=frequency, 
                                            deformation_type=sim_props['deformation_type'])      
        # Get results
        print(f"\n--- Running {sim_props['deformation_type']} Simulation ---")
        results = run_simulation(material, 
                                t_array, 
                                my_profile,
                                mode=sim_props['deformation_type'], mu=mu_fixed)
        # Extract data and plot
        stretches = results['def_val']
        # Extract the P11 component of stress (sigma_11)
        stress_11 = [s[0, 0] for s in results['stress']]
        stress_12 = [s[0, 1] for s in results['stress']]
        plt.figure(figsize=(8, 6))
        plt.plot(stretches, stress_12, linewidth=2)
        d = np.loadtxt('validation_data/unterberger2013_fig4a_laos.csv', delimiter=',')
        gamma, sigma, branch = d[:, 0], d[:, 1], d[:, 2]
        plt.plot(gamma, sigma, 'kx', ms=5, label='Unterberger et al. 2013 (LAOS)')
        plt.xlabel('Shear strain ($\\gamma$)')
        # plt.xlabel('Shear Strain ($\\gamma$)')
        plt.ylabel('Shear stress ($\\sigma_{12}$)')
        plt.title(f'Dynamic Crosslinking: {sim_props["deformation_type"].capitalize()}')
        plt.legend()
        plt.grid(True)
        # plt.savefig(f'figures/{sim_props["deformation_type"]}/uni_results.pdf', bbox_inches='tight')
        plt.savefig(f'figures/validation/{tag}.pdf', bbox_inches='tight')
        # plt.show()

        # LAST CYCLE PLOT
        t_hist = np.array(results['time'])
        period = 1.0 / frequency
        last = t_hist >= t_hist[-1] - period - 0.5 * sim_props['delta_t']   # last full cycle (both end points included)
        gamma_last = np.array(stretches)[last]
        sigma_last = np.array(stress_12)[last]

        d0 = np.loadtxt('validation_data/unterberger2013_fig4a_laos_dataset0.csv', delimiter=',', comments='#')
        g0, s0, b0 = d0[:, 0], d0[:, 1], d0[:, 2]

        dg_last = np.gradient(gamma_last)
        s_model = np.empty_like(s0)
        for b in (1, -1):
            msk = dg_last > 0 if b == 1 else dg_last < 0
            o = np.argsort(gamma_last[msk])
            s_model[b0 == b] = np.interp(g0[b0 == b], gamma_last[msk][o], sigma_last[msk][o])
        rmse = np.sqrt(np.mean((s_model - s0)**2))
        r2 = 1.0 - np.sum((s_model - s0)**2) / np.sum((s0 - np.mean(s0))**2)

        plt.figure(figsize=(8, 6))
        plt.plot(gamma_last, sigma_last, linewidth=2,
                label=f'Model (cycle {int(round(t_hist[-1] / period))}): RMSE = {rmse:.3f} Pa, $R^2$ = {r2:.4f}')    
        plt.plot(g0, s0, 'kx', ms=5, label='Unterberger et al. 2013')
        # plt.plot(g0[b0 == 1], s0[b0 == 1], 'kx', ms=5, label='Unterberger et al. 2013 (loading)')
        # plt.plot(g0[b0 == -1], s0[b0 == -1], 'k+', ms=6, label='Unterberger et al. 2013 (unloading)')
        plt.xlabel('Shear strain ($\\gamma$)')
        plt.ylabel('Shear stress ($\\sigma_{12}$)')
        plt.title(f'Dynamic Crosslinking: {sim_props["deformation_type"].capitalize()} (last cycle)')
        plt.legend()
        plt.grid(True)
        Gp_A, Gpp_A, delta_A = moduli(t_hist[last], gamma_last, sigma_last, frequency)
        print(f"gamma_A = {sim_props['def_max']}: G' = {Gp_A:.2f} Pa, G'' = {Gpp_A:.2f} Pa, "
            f"delta = {np.degrees(delta_A):.1f} deg  (Fig. 4b at 0.126: G' = 29.6, G'' = 6.6)")
        plt.savefig(f'figures/validation/{tag}_last_cycle.pdf', bbox_inches='tight')

    # ----------- AMPLITUDE SWEEP (Fig. 4b) -----------
    run_sweep = True
    if run_sweep:
        dG1 = np.loadtxt('validation_data/4b_G_.csv', delimiter=',')
        dG2 = np.loadtxt('validation_data/4b_G__.csv', delimiter=',')
        dG2 = dG2[np.argsort(dG2[:, 0])]
        max_gamma = 0.2
        gamma_As = dG1[dG1[:, 0] <= max_gamma, 0]
        n_cycles = cfg['cycles'] if RUN else 3
        dt_sweep = cfg['dt'] if RUN else 0.1
        jobs = [(gA, diffusion_props, mech_props, sim_props, frequency, n_cycles, dt_sweep, mu_fixed)
                for gA in gamma_As]
        n_workers = 16
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            out = list(ex.map(sweep_worker, jobs))                 # same order as gamma_As
        Gp = [o[1] for o in out]
        Gpp = [o[2] for o in out]

        np.savetxt(f'validation/{tag}_amplitude_sweep.csv', np.column_stack([gamma_As, Gp, Gpp]),
           delimiter=',', header="gamma_A, G' (Pa), G'' (Pa)")
        plt.figure(figsize=(8, 6))
        plt.loglog(gamma_As, Gp, 'b-', lw=2, label="Model $G'$")
        plt.loglog(gamma_As, Gpp, 'r-', lw=2, label="Model $G''$")
        plt.loglog(dG1[:, 0], dG1[:, 1], 'bx', ms=6, label="Unterberger et al. 2013 $G'$")
        plt.loglog(dG2[:, 0], dG2[:, 1], 'r+', ms=7, label="Unterberger et al. 2013 $G''$")
        plt.xlabel('Amplitude of shear strain $\\gamma_A$')
        plt.ylabel("Moduli $G'$, $G''$ (Pa)")
        plt.title(f'Amplitude sweep ({tag})'.replace('_', '\\_'))
        plt.legend()
        plt.grid(True, which='both')
        plt.savefig(f'figures/validation/{tag}_amplitude_sweep.pdf', bbox_inches='tight')


    # --------------------------------------------
    # ----------- SENSITIVITY ANALYSIS -----------
    # --------------------------------------------

    # The values you want to test                                                                                                                                                
    sensitivity_array = [] # np.array([0.01, 0.05, 0.1])  # Frequencies in Hz for cyclic loading                                                                                             
                                                                                                                                                                                    
    # fig, ax = plt.subplots(3,2,figsize=(14, 16))
    # fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    # colors = ['b', 'orange','g', 'r', 'c', 'm', 'y', 'k']  # Add more colors if needed                                                                                                                                          
    colors = ['midnightblue', 'mediumblue', 'steelblue', 'firebrick', 'darkcyan', 'darkmagenta', 'goldenrod', 'black']  # Add more colors if needed
    colors = ['saddlebrown', 'darkorange', 'gold']
    labels_list = ['$\\mathbf{1}$ \\textbf{s}', '$\\mathbf{10}$ \\textbf{s}', '$\\mathbf{60}$ \\textbf{s}']
    # y_pos = [-0.035, -0.01, -0.08]                                                                                                                                                                     
    for i, param in enumerate(sensitivity_array):
        # break                                                                                                                                                        
        print(f"Running simulation for t = {param}...")                                                                                                                       
                                                                                                                                                  
        # Update the parameter
        diffusion_props['k_catch_0'] = param
        # sim_props['t_max'] = param                                                                                                                                
                                                                                                                                                                                    
        # Instantiate material                                                                                                                                    
        material = DynamicCrosslinking(diffusion_props, mech_props, sim_props)
        mu_fixed = material.INITMU  # Keeping mu constant for now
        print(f"Material initialized with mu = {mu_fixed:.4f}")                                                                              
                                                                                                                                                                                    
        # Time array and deformation profile
        # ----- LINEAR -----                                                                                                                                  
        # t_array = np.linspace(sim_props['t_init'], sim_props['t_max'], 75)
        # t_array = np.arange(sim_props['t_init'], sim_props['t_max'] + sim_props['delta_t'], sim_props['delta_t']) 
        # stretch_rate = (sim_props['def_max'] - sim_props['def_init']) / sim_props['t_max']                                                                                       
        # my_profile = lambda t: profile_linear(t, rate=stretch_rate, deformation_type=sim_props['deformation_type'])
        # ----- STEP -----
        # t_step = 1.0                                                                                                                                
        # t_array = np.linspace(sim_props['t_init'], sim_props['t_max'], 100)
        # t_array = np.arange(sim_props['t_init'], sim_props['t_max'] + sim_props['delta_t'], sim_props['delta_t'])                                                                                                       
        # my_profile = lambda t: profile_step(t, t_step=t_step, 
                                            # target_stretch=sim_props['def_max'], 
                                            # deformation_type=sim_props['deformation_type'])                                                                                                    
        # ----- CYCLIC -----
        frequency = 1.0 # param  # Frequency of cyclic loading (Hz)
        t_array = np.arange(sim_props['t_init'], sim_props['t_max'] + sim_props['delta_t'], sim_props['delta_t'])                                                                                                       
        my_profile = lambda t: profile_cyclic(t, amplitude=sim_props['def_max'], 
                                              frequency=frequency, 
                                              deformation_type=sim_props['deformation_type'])

        # Run and Plot                                                                                                                                                        
        results = run_simulation(material, t_array, my_profile, mode=sim_props['deformation_type'], mu=mu_fixed)                                                                    

        stretches = results['def_val']
        cf = results['cf']                                                                                                                                       
        cb = results['cb']
        kon_mean = [np.mean(kon_i) for kon_i in results['kon']]
        koff_mean = [np.mean(koff_i) for koff_i in results['koff']]                                                               
        stress_11 = [s[0, 0] for s in results['stress']]
        stress_12 = [s[0, 1] for s in results['stress']]
        thetaf = [c / 0.0095 for c in results['cf']]
        cb_total = [sum(cb_i *1000 * material.weights) for cb_i in results['cb']]
        max_cb = [max(cb_i*1000) for cb_i in results['cb']]
        min_cb = [min(cb_i*1000) for cb_i in results['cb']]

        # Deformation in x-axis          
        ax[0,0].plot(stretches, stress_11, linewidth=2, label=f"$k_{{catch,0}} = {param}$", color=colors[i])
        # ax.text(stretches[-1] - 0.03, stress_11[-1] + y_pos[i], labels_list[i], color=colors[i], fontsize=14, fontweight='bold', va='center')
        ax[0,1].plot(stretches, thetaf, linewidth=2, color=colors[i])
        ax[1,0].plot(stretches, cb_total, linewidth=2, color=colors[i])
        ax[1,0].plot(stretches, max_cb, linewidth=2, linestyle='--', label=f"Max $c_{{b,i}}$", color=colors[i])
        ax[1,0].plot(stretches, min_cb, linewidth=2, linestyle='-.', label=f"Min $c_{{b,i}}$", color=colors[i])
        ax[1,1].plot(stretches, np.array(cb_total) + np.array(cf), linewidth=2, color=colors[i])
        ax[2,1].plot(stretches, kon_mean, linewidth=2, color=colors[i])
        ax[2,0].plot(stretches, koff_mean, linewidth=2, color=colors[i])

        # Time in x-axis
        # ax[0,0].plot(t_array, stress_12, linewidth=2, label=f"$R_{{b}}^{{max}} = {param}$", color=colors[i])
        # ax[0,1].plot(t_array, thetaf, linewidth=2, color=colors[i])
        # ax.plot(t_array, cb_total, linewidth=2, color=colors[i])
        # ax.plot(t_array, max_cb, linewidth=2, linestyle='--', color=colors[i])
        # ax.plot(t_array, min_cb, linewidth=2, linestyle='dotted', color=colors[i])
        # ax[1,1].plot(t_array, np.array(cb_total) + np.array(cf), linewidth=2, color=colors[i])
        # ax[2,1].plot(t_array, kon_mean, linewidth=2, color=colors[i])
        # ax[2,0].plot(t_array, koff_mean, linewidth=2, color=colors[i])

        # Cyclic-specific plots

        # 1) Histeresis loops: stress vs deformation
        # ax[i,0].plot(stretches, stress_12, linewidth=2, label=f"$\\omega = {param}$", color=colors[i])
        # ax[i,0].set_xlabel('Shear ($\\gamma$)')
        # ax[i,0].set_ylabel('Cauchy Stress ($\\sigma_{12}$)')
        # ax[i,0].legend()
        # ax[i,0].grid(True, alpha=0.3)
        # ax[i,1].plot(t_array, cb_total, linewidth=2, color=colors[i])
        # ax[i,1].set_xlabel('Time ($t$)')
        # ax[i,1].set_ylabel('Total Bound Crosslinkers ($c_b$)')
        # ax[i,1].grid(True, alpha=0.3)
        # 2) Peak stress vs frequency


    # ----- PLOT FORMATTING -----

    # # (0,0) - Stress vs deformation/time
    # # ax[0,0].set_xlabel('Stretch ($\\lambda_1$)')
    # ax.set_xlabel('Shear Strain $\\gamma$', fontsize=15)
    # ax.set_xlabel('Stretch $\\lambda$', fontsize=15)
    # ax.set_xticks([0, 0.1, 0.2, 0.3])
    # ax.set_xticks([1.0, 1.1, 1.2])
    # # ax[0,0].set_xlabel('Time ($t$)')
    # ax.set_ylabel('Normal Stress $\\sigma_{11}$ (Pa)', fontsize=15)
    # ax.set_yticks(3*np.array([0, 0.05, 0.1, 0.15]))
    # ax.tick_params(labelsize=12)
    # # ax[0,0].set_title('Sensitivity to Unbinding Rate')
    # ax[0,0].legend()
    # ax[0,0].grid(True, alpha=0.3)

    # # (0,1) - Thetaf vs deformation/time
    # # ax[0,1].set_xlabel('Stretch ($\\lambda_1$)')
    # ax[0,1].set_xlabel('Shear ($\\gamma$)')
    # # ax[0,1].set_xlabel('Time ($t$)')
    # ax[0,1].set_ylabel('Crosslinker Fraction ($\\theta_f$)')
    # ax[0,1].set_ylim([0.05, 0.1])
    # # ax[0,1].legend()
    # ax[0,1].grid(True, alpha=0.3)

    # # (1,0) - Total Bound Crosslinkers vs deformation/time
    # # ax[1,0].set_xlabel('Stretch ($\\lambda_1$)')
    # ax[1,0].set_xlabel('Shear ($\\gamma$)')
    # ax.set_xlabel('Time ($t$)')
    # ax.set_ylabel('Bound Crosslinker Concentration $c_b$ ($\\mu$M)')
    # # ax[1,0].set_ylim([0.0, 0.25])
    # ax[1,0].legend()
    # ax[1,0].grid(True, alpha=0.3)

    # custom_lines = [
    #     Line2D([0], [0], color='black', lw=2, linestyle='-'),
    #     Line2D([0], [0], color='black', lw=2, linestyle='--'),
    #     Line2D([0], [0], color='black', lw=2, linestyle='dotted')
    #                 ]
  
    # Apply the custom legend to the axis
    # ax.legend(custom_lines, ['Homogenized $c_b$', 'Max $c_{b,i}$', 'Min $c_{b,i}$'], 
    #           loc='best', fontsize=12)

    # ax.grid(True, alpha=0.2)

    # # (1,1) - Total Crosslinkers vs deformation/time
    # # ax[1,1].set_xlabel('Stretch ($\\lambda_1$)')
    # ax[1,1].set_xlabel('Shear ($\\gamma$)')
    # # ax[1,1].set_xlabel('Time ($t$)')
    # ax[1,1].set_ylabel('Total Crosslinkers ($c_b + c_f$)')
    # # ax[1,1].set_ylim([0.0, 0.25])
    # # ax[1,1].legend()
    # ax[1,1].grid(True, alpha=0.3)

    # # (2,0) - Mean koff vs deformation/time
    # # ax[2,0].set_xlabel('Stretch ($\\lambda_1$)')
    # ax[2,0].set_xlabel('Shear ($\\gamma$)')
    # # ax[2,0].set_xlabel('Time ($t$)')
    # ax[2,0].set_ylabel('Mean $k_{{off}}$')
    # # ax[2,0].legend()
    # ax[2,0].grid(True, alpha=0.3)

    # # (2,1) - Mean kon vs deformation/time
    # # ax[2,1].set_xlabel('Stretch ($\\lambda_1$)')
    # ax[2,1].set_xlabel('Shear ($\\gamma$)')
    # # ax[2,1].set_xlabel('Time ($t$)')
    # ax[2,1].set_ylabel('Mean $k_{{on}}$')
    # # ax[2,1].legend()
    # ax[2,1].grid(True, alpha=0.3)

    # ----- SAVE PLOTS -----
    # fig.savefig(f'figures/catch-slip/{sim_props["deformation_type"]}_kcatch_20s.pdf', bbox_inches='tight')
    # plt.show()

    # ---------------------------------------------------------
    # ----------- SWEEP & CYCLIC SOFTENING SETTINGS -----------
    # ---------------------------------------------------------
                                                                                                                 
    # frequencies = np.logspace(-3, 1, 15)
    # frequencies = np.array([1.])  # For testing purposes                                                                                                                                        
                                                                                                                                                                                    
    # We will pick the 1st, 8th, and 15th frequencies to plot side-by-side                                                                                                       
    # indices_to_plot = [0] #[0, 7, 14]                                                                                                                                                 
    # colors = ['saddlebrown', 'darkorange', 'gold','midnightblue', 'mediumblue', 'steelblue']# ['r', 'g', 'b']                                                                                                                                                     
    # cycle_color_idx = 0                                                                                                                                                          
                                                                                                                                                                                    
    # num_cycles = 10  # 10 cycles gives the material time to plateau to steady-state                                                                                              
    # pts_per_cycle = 100                                                                                                                                                          
    # num_points = num_cycles * pts_per_cycle                                                                                                                                      
                                                                                                                                                                                    
    # Create a 1x2 side-by-side figure                                                                                                                                           
    # fig, ax = plt.subplots(1, 1, figsize=(6, 4))                                                                                                                                

    # for ik, koff0 in enumerate([0.01, 0.1, 1.]):  # Test different baseline off-rates
    #     peak_stresses_steady = []                                                                                                                                                    
    #     diffusion_props['k_off_0'] = koff0                                                                                                                                                              
    #     for freq_idx, freq in enumerate(frequencies):                                                                                                                                
    #         print(f"Running frequency: {freq:.3f} Hz...")                                                                                                                            
                                                                                                                                                                                        
    #         # Time Array                                                                                                                                                    
    #         period = 1.0 / freq                                                                                                                                                      
    #         t_max = num_cycles * period                                                                                                                                              
    #         t_array = np.linspace(0, t_max, num_points)                                                                                                                              
                                                                                                                                                                                        
    #         # Instantiate material                                                                                                                                             
    #         material = DynamicCrosslinking(diffusion_props, mech_props, sim_props)                                                                                                   
    #         mu_fixed = material.INITMU                                                                                                                                               
                                                                                                                                                                                        
    #         # Deformation Profile                                                                                                                                                       
    #         amplitude = sim_props['def_max']                                                                                                                                         
    #         my_profile = lambda t: profile_cyclic(t, amplitude, freq, sim_props['deformation_type'])                                                                                 
                                                                                                                                                                                        
    #         # Run Simulation                                                                                                                                                      
    #         results = run_simulation(material, t_array, my_profile, mode=sim_props['deformation_type'], mu=mu_fixed)

    #         # Microsphere anisotropy plotting
    #         # if freq_idx == 0 or freq_idx == len(frequencies) - 1:
    #         #     # 1. Start of the last cycle 
    #         #     idx_start = -pts_per_cycle
    #         #     # 2. End of the last cycle 
    #         #     idx_end = -1
    #         #     # 3. Peak of the last cycle
    #         #     last_cycle_defs = results['def_val'][-pts_per_cycle:]
    #         #     idx_peak = -pts_per_cycle + np.argmax(last_cycle_defs)
    #         #     # Now call the plot function 3 times to see the evolution!
    #         #     title_str = ""# f"Anisotropy at Peak Deformation ({freq:.3f} Hz)"
    #         #     path1 = f'figures/{sim_props["deformation_type"]}/microsphere_{sim_props["deformation_type"]}_{freq_idx}Hz_t0.pdf'
    #         #     path2 = f'figures/{sim_props["deformation_type"]}/microsphere_{sim_props["deformation_type"]}_{freq_idx}Hz_tdefmax.pdf'
    #         #     path3 = f'figures/{sim_props["deformation_type"]}/microsphere_{sim_props["deformation_type"]}_{freq_idx}Hz_tend.pdf'
    #         #     plot_microsphere(material.dirs, results['cb'][idx_start], material.cb0, 
    #         #                         title=f"1. Start of Cycle ({freq:.3f} Hz)",
    #         #                         path=path1)
    #         #     plot_microsphere(material.dirs, results['cb'][idx_peak], material.cb0, 
    #         #                         title=f"2. Peak Stretch ({freq:.3f} Hz)",
    #         #                         path=path2)
    #         #     plot_microsphere(material.dirs, results['cb'][idx_end], material.cb0, 
    #         #                         title=f"3. End of Cycle ({freq:.3f} Hz)",
    #         #                         path=path3)

            
    #         # Reshape the 1D stress array into a 2D matrix: (10 cycles, 100 points)
    #         stress_12 = [s[0, 1] for s in results['stress']]
    #         stress_matrix = np.array(stress_12).reshape((num_cycles, pts_per_cycle))
            
    #         # # Get the maximum positive stress for every single cycle
    #         peaks_per_cycle = np.max(stress_matrix, axis=1)
            
    #         # # The steady-state peak is the peak of the very last cycle
    #         steady_peak = peaks_per_cycle[-1]
    #         peak_stresses_steady.append(steady_peak)
            
    #         # # If this frequency is one of our chosen 3, plot its cyclic softening!
    #         if freq_idx in indices_to_plot:
    #             cycle_numbers = np.arange(1, num_cycles + 1)
    #             ax.plot(cycle_numbers, peaks_per_cycle, marker='o', linewidth=2, 
    #                         color=colors[ik], )# label=f'{freq:.3f} Hz')
                # cycle_color_idx += 1
        # ax.plot(frequencies, peak_stresses_steady, 'k-o', linewidth=2, label=koff0, color=colors[ik])   

    # ----------- PLOT 1: FREQUENCY SWEEP -----------
    # ax.set_xscale('log')
    # ax.set_xlabel('Frequency $f$ (Hz)',fontsize=15)
    # ax.set_ylabel('Peak Shear Stress $\\sigma_{12}$ (Pa)',fontsize=15)
    # ax.set_yticks([0.1, 0.2, 0.3])
    # ax.tick_params(labelsize=12)
    # # # ax.set_title('Dynamic Stiffness vs. Frequency')
    # ax.grid(True, which="both", alpha=0.2)
    # ax.legend(title='$k_\\mathrm{{off}}^0$ (s$^{-1}$)', title_fontsize=14, fontsize=13)

    # ----------- PLOT 2: CYCLIC SOFTENING -----------
    # ax.set_xlabel('Cycle Number')
    # ax.set_ylabel('Peak Stress per Cycle ($\\sigma_{12}$)')
    # ax[1].set_title('Cyclic Softening (Mullins-like effect)')
    # ax.set_xticks(np.arange(1, num_cycles + 1))  # Force integer ticks for cycles
    # ax.grid(True, ls="--")
    # ax.grid(True, which="both", alpha=0.2)
    # ax.legend()
    
    # plt.tight_layout()
    # fig.savefig(f'figures/{sim_props["deformation_type"]}/{sim_props["deformation_type"]}_imgforppt22_new.svg', bbox_inches='tight')
    # fig.savefig(f'figures/{sim_props["deformation_type"]}/{sim_props["deformation_type"]}_imgforppt22_new.pdf', bbox_inches='tight')

    # plt.show()
