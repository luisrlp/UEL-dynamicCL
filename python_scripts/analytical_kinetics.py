import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import os
import sys

# Ensure the script can import filament.py
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from filament import Filament

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
plt.rcParams.update(tex_fonts)

# ==========================================
# 1. PARAMETERS (Refined against main.pdf)
# ==========================================
# Kinetic parameters
k_off_0 = 0.05      # 1/s (baseline off-rate) [0.05]
K_eq = 0.1        # Equilibrium constant [0.25]
chi_FH = 0.1      # Flory-Huggins parameter
dx_kT = 0.1 / (1.38065e-5 * 298)       # Reactive distance over thermal energy: Delta x / (k_B T)

# Concentration constraints
c_actin = 0.0095     # Reference actin concentration
R_total = 0.1      # Initial total crosslinker ratio (cb + cf)
Rb_max = 0.25      # Max bound ratio
Rf_max = 1.0       # Max free ratio (dummy large)

c_total = c_actin * R_total
cb_max = c_actin * Rb_max
cf_max = c_actin * Rf_max

# Filament properties
mat_props_base = {
    'LAMBDA0': 1.0,
    'MU0': 38600.0,
    'Lp': 16.0,
    'BETA': 0.5,
    'R0C': 0.014,
    'ETA': 0.5,
    'a_ratio': 1.2  # L / r0f ratio
}

# Instantiate Filament class
def_info = {'def_initial': 1.0, 'def_max': 1.2, 'increments': 100}
sim_dir = os.path.dirname(os.path.abspath(__file__))
fil = Filament(def_info, sim_dir=sim_dir)

# ==========================================
# 2. CORE FUNCTIONS
# ==========================================
def get_geometry_and_force(cb, stretch):
    """Calculates dynamically updated geometry and filament force."""
    props = mat_props_base.copy()
    props['STRETCH'] = stretch
    # print(f"Calculating geometry and force for cb={cb:.6f}, stretch={stretch:.3f}")
    
    # Eq 62: r0,f = 1.6 * (10^3 * cb)^(-2/5)
    # Assuming cb is in mM, so 10^3 cb is in microM. 
    # Let's assume our cb is already properly scaled for this empirical law.
    cb_in_mM = max(cb * 1000.0, 1e-6) # prevent div by zero
    r0f = 1.6 * (cb_in_mM)**(-0.4)
    
    # Eq 63: L = a * r0,f
    L = props['a_ratio'] * r0f
    
    props['R0F'] = r0f
    props['L'] = L
    
    try:
        res = fil.fil_force(props, cl=True)
        force = res['force'][0]
        if np.isnan(force):
            force = 0.0
    except Exception as e:
        force = 0.0
        
    # Eq 65: n \propto L^{-1}. 
    # The true macroscopic stress will scale with n * f \propto f / L
    # We will return the scaling factor 1/L for the network stress.
    n_scaling = 1.0 / L
    
    return force, n_scaling

def reaction_rate(t, cb_array, stretch_func):
    """ODE for bound crosslinkers."""
    cb = cb_array[0]
    
    # Constrain concentrations
    cb = max(1e-6, min(cb, c_total - 1e-6))
    cf = c_total - cb
    
    # Eq 66: theta_b and theta_f
    theta_b = cb / cb_max
    theta_f = cf / cf_max
    
    # Get current stretch and force
    lam = stretch_func(t)
    f, _ = get_geometry_and_force(cb, lam)
    
    # Eq 67: Forward binding rate
    k_on = k_off_0 * K_eq * np.exp(chi_FH * (1.0 - 2.0 * theta_f))
    
    # Eq 68: Slip-bond unbinding rate
    k_off_i = k_off_0 * np.exp(min(f * dx_kT, 200.0))
    if f * dx_kT > 200.0:
        print(f"Warning: Force-induced unbinding rate is very high at t={t:.2f}s, f={f:.2f}pN.")
    
    # Eq 69: Reaction rate equation (scaled by max concentrations)
    term1 = k_on * cf_max * (theta_f / (1.0 - theta_f))
    term2 = k_off_i * cb_max * (theta_b / (1.0 - theta_b))
    
    dcb_dt = term1 - term2
    return [dcb_dt]

from scipy.optimize import brentq

# ==========================================
# 3. INITIALIZATION (Eq 59)
# ==========================================
def get_equilibrium_cb0():
    """Finds the thermodynamic equilibrium cb_0 at t=0 (undeformed state)."""
    def H(cb0):
        cf0 = c_total - cb0
        theta_f0 = cf0 / cf_max
        
        # Eq 59: H(cb0) = 0
        lhs = K_eq * cf0 * (1.0 - cb0 / cb_max)
        rhs = cb0 * (1.0 - theta_f0) * np.exp(-chi_FH * (1.0 - 2.0 * theta_f0))
        return lhs - rhs

    # cb0 must be strictly between 0 and the maximum physically possible bound
    upper_bound = min(c_total, cb_max) - 1e-8
    cb0_eq = brentq(H, 1e-8, upper_bound)
    return cb0_eq

# ==========================================
# 4. ANALYTICAL EXPERIMENTS
# ==========================================

def experiment_constant_stretch():
    """Applies a sudden step stretch and plots relaxation."""
    print("Running Constant Stretch Experiment...")
    t_span = (0, 10)
    t_eval = np.linspace(t_span[0], t_span[1], 200)
    
    lam_step = 1.15
    stretch_func = lambda t: lam_step if t > 0.5 else 1.0
    
    # Calculate exact equilibrium initial condition
    cb_eq = get_equilibrium_cb0()
    print(f"Calculated equilibrium cb_0: {cb_eq:.4f}")
    cb_0 = [cb_eq]
    
    sol = solve_ivp(reaction_rate, t_span, cb_0, args=(stretch_func,), t_eval=t_eval, method='BDF')
    
    # Post-process forces and stress
    forces = []
    macroscopic_stress = []
    k_offs = []
    
    for cb, t in zip(sol.y[0], sol.t):
        f, n_scaling = get_geometry_and_force(cb, stretch_func(t))
        forces.append(f)
        macroscopic_stress.append(n_scaling * f)
        k_offs.append(k_off_0 * np.exp(min(f * dx_kT, 100.0)))
    
    # Plotting
    fig, axs = plt.subplots(3, 1, figsize=(8, 10))
    axs[0].plot(sol.t, sol.y[0], label='Bound ($c_b$)')
    axs[0].plot(sol.t, c_total - sol.y[0], label='Free ($c_f$)')
    axs[0].set_ylabel('Concentration')
    axs[0].legend()
    axs[0].set_title('Crosslinker Evolution (Step Stretch)')
    
    axs[1].plot(sol.t, macroscopic_stress, color='red')
    axs[1].set_ylabel('Macroscopic Stress ($\\propto f_i / L$)')
    axs[1].set_title('Stress Relaxation')
    
    axs[2].plot(sol.t, k_offs, color='green')
    axs[2].set_ylabel('Unbinding Rate ($k_\\mathrm{off}$)')
    axs[2].set_xlabel('Time (s)')
    
    plt.tight_layout()
    plt.savefig(os.path.join(sim_dir, 'exp1_relaxation.pdf'))
    plt.close()

def experiment_cyclic_loading():
    """Applies a sinusoidal stretch and plots hysteresis."""
    print("Running Cyclic Loading Experiment...")
    t_span = (0, 600)
    t_eval = np.linspace(t_span[0], t_span[1], 10000)
    
    omegas = [0.01, 2.] # , 1.0]
    
    # Calculate exact equilibrium initial condition
    cb_eq = get_equilibrium_cb0()
    cb_0 = [cb_eq]
    
    fig, axs = plt.subplots(3, 2, figsize=(12, 14))
    
    for omega in omegas:
        stretch_func = lambda t: 1.05 - 0.05 * np.cos(omega * t)
        sol = solve_ivp(reaction_rate, 
                        t_span, 
                        cb_0, 
                        args=(stretch_func,), 
                        t_eval=t_eval, 
                        method='BDF',
                        rtol=1e-8, atol=1e-10)
        
        stretches = []
        macroscopic_stress = []
        net_rates = []
        
        for cb, t in zip(sol.y[0], sol.t):
            lam = stretch_func(t)
            f, n_scaling = get_geometry_and_force(cb, lam)
            stretches.append(lam)
            macroscopic_stress.append(n_scaling * f)
            rate = reaction_rate(t, [cb], stretch_func)[0]               
            net_rates.append(rate)  
        
        axs[0,0].plot(sol.t, macroscopic_stress, label=f'$\\omega$ = {omega}')
        
        idx_steady = len(sol.t) // 2
        axs[0,1].plot(stretches[idx_steady:], macroscopic_stress[idx_steady:], label=f'$\\omega$ = {omega}')
        axs[1,0].plot(sol.t, sol.y[0], label=f'Bound ($c_b$), $\\omega$ = {omega}')

        # 1. Find the indices of the peaks in the stress array           
        peak_indices, _ = find_peaks(macroscopic_stress)                 
                                                                            
        # 2. Extract the time and stress values at those exact peaks     
        peak_times = [sol.t[i] for i in peak_indices]
        peak_stresses = [macroscopic_stress[i] for i in peak_indices]    
        
        # 3. Plot Peak Stress vs Time to clearly see the transient decay 
        axs[1,1].plot(peak_times, peak_stresses, marker='o', linestyle='--',
                    label=f'$\\omega$ = {omega}')
        
        axs[2,0].plot(sol.t, net_rates, label=f'$\\omega$ = {omega}')

        axs[2,1].plot(sol.t, sol.y[0]/cb_max, label=f'Bound ($\\theta_b$), $\\omega$ = {omega}')


    axs[0,0].set_xlabel('Time (s)')
    axs[0,0].set_ylabel('Macroscopic Stress')
    axs[0,0].legend()
    axs[0,0].set_title('Stress vs. Time')
    
    axs[0,1].set_xlabel('Stretch ($\\lambda$)')
    axs[0,1].set_ylabel('Macroscopic Stress')
    axs[0,1].legend()
    axs[0,1].set_title('Hysteresis Loops')

    # axs[1,0].plot(sol.t, c_total - sol.y[0], label='Free ($c_f$)')
    axs[1,0].set_xlabel('Time (s)')
    axs[1,0].set_ylabel('Concentration')
    axs[1,0].legend()
    axs[1,0].set_title('Crosslinker Evolution (Cyclic Loading)')

    axs[1,1].set_xlabel('Time (s)')
    axs[1,1].set_ylabel('Peak Macroscopic Stress')
    axs[1,1].legend()
    axs[1,1].set_title('Peak Stress Evolution')

    axs[2,0].axhline(0, color='black', linestyle='--', linewidth=0.8) 
    axs[2,0].set_xlabel('Time (s)')
    axs[2,0].set_ylabel('Net Reaction Rate $\\mathcal{R}_i$')     
    axs[2,0].set_title('Reaction Rate over Time')
    axs[2,0].legend()

    axs[2,1].set_xlabel('Time (s)')
    axs[2,1].set_ylabel('Bound Fraction ($\\theta_b$)')
    axs[2,1].set_title('Bound Fraction over Time')
    axs[2,1].legend()


    plt.tight_layout()
    plt.savefig(os.path.join(sim_dir, 'exp2_hysteresis.pdf'))
    plt.close()

if __name__ == '__main__':
    # experiment_constant_stretch()
    experiment_cyclic_loading()
    print("Plots saved in:", sim_dir)
