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
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams.update(tex_fonts)

# ==========================================
# 1. PARAMETERS
# ==========================================
# Kinetic parameters
k_off_0 = 0.05    # 1/s (baseline off-rate) [0.05]
k_on_0 = 0.0125   # 1/s (baseline on-rate) [0.0125]
chi_FH = 1.     # Flory-Huggins parameter
dx_kT = 0.1 / (1.38065e-5 * 298)       # Reactive distance over thermal energy: Delta x / (k_B T)

# Concentration constraints
c_actin = 0.0095     # Reference actin concentration
R_total = 0.1      # Initial total crosslinker ratio (cb + cf)
Rb_max = 0.025     # Max bound ratio
Rf_max = 0.1       # Max free ratio (dummy large)

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
        res = fil.fil_force(props, dw=True, cl=True)
        force = res['force'][0]
        dw = res['dw'][0]
        if np.isnan(force):
            force = 0.0
        if np.isnan(dw):
            dw = 0.0
    except Exception as e:
        print(f"Error in fil_force calculation: {e}")
        force = 0.0
        dw = 0.0
        
    # Eq 65: n \propto L^{-1}. 
    # The true macroscopic stress will scale with n * f \propto f / L
    # We will return the scaling factor 1/L for the network stress.
    n_scaling = 1.0 / L
    
    return force, n_scaling, dw

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
    f, _, _ = get_geometry_and_force(cb, lam)
    
    # Eq 67: Forward binding rate
    k_on = k_on_0 * np.exp(chi_FH * (1.0 - 2.0 * theta_f))
    
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
        lhs = (k_on_0 / k_off_0) * cf0 * (1.0 - cb0 / cb_max)
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
    
    lam_step = 1.1
    stretch_func = lambda t: lam_step if t > 0.5 else 1.0
    
    # Calculate exact equilibrium initial condition
    cb_eq = get_equilibrium_cb0()
    print(f"Calculated equilibrium cb_0: {cb_eq:.4f}")
    cb_0 = [cb_eq]
    
    sol = solve_ivp(reaction_rate, t_span, cb_0, args=(stretch_func,), t_eval=t_eval, method='BDF')
    
    # Post-process forces and stress
    forces = []
    dws = []
    macroscopic_stress = []
    k_offs = []
    k_ons = []
    
    for cb, t in zip(sol.y[0], sol.t):
        f, n_scaling, dw = get_geometry_and_force(cb, stretch_func(t))
        forces.append(f)
        dws.append(dw)
        # macroscopic_stress.append(n_scaling * f)
        macroscopic_stress.append(n_scaling * dw)
        k_offs.append(k_off_0 * np.exp(min(f * dx_kT, 500.0)))
        print((c_total - cb) / cf_max)
        k_ons.append(k_on_0 * np.exp(chi_FH * (1.0 - 2.0 * (c_total - cb) / cf_max)))
    
    # Plotting
    fig, axs = plt.subplots(2, 2, figsize=(8, 6))
    axs[0,0].plot(sol.t, sol.y[0])
    axs[0,0].set_ylabel('$c_b$')
    # axs[0,0].legend()
    # axs[0,0].set_title('Crosslinker Evolution (Step Stretch)')
    axs[0,0].grid(True,alpha=0.3)
    axs[0,1].plot(sol.t, c_total - sol.y[0])
    axs[0,1].set_ylabel('$c_f$')
    # axs[0,1].legend()
    # axs[0,1].set_title('Crosslinker Evolution (Step Stretch)')
    axs[0,1].grid(True,alpha=0.3)
    
    axs[1,0].plot(sol.t, k_offs, color='green')
    # axs[1,0].set_ylabel('Unbinding Rate ($k_\\mathrm{off}$)')
    axs[1,0].set_ylabel(r'$k_\mathrm{off} (1/\mathrm{s})$')    
    axs[1,0].set_xlabel('Time (s)')
    axs[1,0].grid(True,alpha=0.3)

    axs[1,1].plot(sol.t, k_ons, color='green')
    # axs[1,1].set_ylabel('Unbinding Rate ($k_\\mathrm{off}$)')
    axs[1,1].set_ylabel(r'$k_\mathrm{on} (1/\mathrm{s})$')    
    axs[1,1].set_xlabel('Time (s)')
    axs[1,1].grid(True,alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(sim_dir, 'exp1_relaxation2.pdf'))
    plt.close()


    plt.figure(figsize=(8, 4))
    plt.plot(sol.t, macroscopic_stress, color='red')
    # plt.set_ylabel('Macroscopic Stress ($\\propto f_i / L$)')
    # plt.set_ylabel('$w_i^{\\prime} n_i$')
    plt.ylabel('$\\sigma_i^*$')
    plt.title('Stress Relaxation')
    plt.grid(True,alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(sim_dir, 'exp1_relaxation3.pdf'))
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
            f, n_scaling, _ = get_geometry_and_force(cb, lam)
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

def experiment_sensitivity(param_name='k_off_0', param_values=[0.005, 0.05, 0.5]):
    """Runs a sensitivity analysis for a given parameter, outputting multi-panel plots."""
    print(f"Running Sensitivity Analysis for {param_name}...")
    global k_off_0, k_on_0, dx_kT, Rb_max, R_total, cb_max, c_total
    
    # Store originals
    orig_k_off_0 = k_off_0
    orig_k_on_0 = k_on_0
    orig_dx_kT = dx_kT
    orig_Rb_max = Rb_max
    orig_R_total = R_total
    
    t_span = (0, 10)
    t_eval = np.linspace(t_span[0], t_span[1], 250)
    
    lam_step = 1.2
    stretch_func = lambda t: lam_step if t > 0.5 else 1.0
    
    # Create the figures
    fig1, axs1 = plt.subplots(2, 2, figsize=(8, 6))
    fig2, ax2 = plt.subplots(figsize=(8, 4))
    
    for val in param_values:
        # Set parameter
        if param_name == 'k_off_0': k_off_0 = val
        elif param_name == 'k_on_0': k_on_0 = val
        elif param_name == 'dx_kT': dx_kT = val * 1.38065e-5 * 298
        elif param_name == 'dx': dx_kT = val / (1.38065e-5 * 298)
        elif param_name == 'Rb_max': 
            Rb_max = val
            cb_max = c_actin * Rb_max
        elif param_name == 'R_total':
            R_total = val
            c_total = c_actin * R_total
            
        cb_eq = get_equilibrium_cb0()
        cb_0 = [cb_eq]
        
        sol = solve_ivp(reaction_rate, t_span, cb_0, args=(stretch_func,), t_eval=t_eval, method='BDF')
        
        forces = []
        dws = []
        macroscopic_stress = []
        k_offs = []
        k_ons = []
        
        for cb, t in zip(sol.y[0], sol.t):
            f, n_scaling, dw = get_geometry_and_force(cb, stretch_func(t))
            forces.append(f)
            dws.append(dw)
            macroscopic_stress.append(n_scaling * dw)
            k_offs.append(k_off_0 * np.exp(min(f * dx_kT, 500.0)))
            print((c_total - cb) / cf_max)
            k_ons.append(k_on_0 * np.exp(chi_FH * (1.0 - 2.0 * (c_total - cb) / cf_max)))
            
        if param_name == 'k_off_0':
            label = f'$k_{{\\mathrm{{off}}}}^0 = {val}$'
        elif param_name == 'k_on_0':
            label = f'$k_{{\\mathrm{{on}}}}^0 = {val}$'
        elif param_name == 'dx':
            label = f'$\\Delta x = {val}$'
        elif param_name == 'Rb_max':
            label = f'$R_{{b,\\mathrm{{max}}}} = {val}$'
        elif param_name == 'R_total':
            label = r'$R = ' + str(val) + '$'
        else:
            label = f'${param_name} = {val}$'
        
        # Fig 1 plots
        axs1[0,0].plot(sol.t, sol.y[0] * 1000.0, label=label)
        axs1[0,1].plot(sol.t, (c_total - sol.y[0]) * 1000.0, label=label)
        axs1[1,0].plot(sol.t, k_offs, label=label)
        axs1[1,1].plot(sol.t, k_ons, label=label)
        
        # Fig 2 plot
        ax2.plot(sol.t, macroscopic_stress, label=label)
        
    # Format Fig 1
    import matplotlib.ticker as ticker
    
    class CustomScalarFormatter(ticker.ScalarFormatter):
        def __init__(self, decimals=2, **kwargs):
            super().__init__(**kwargs)
            self.decimals = decimals
            
        def __call__(self, x, pos=None):
            s = super().__call__(x, pos)
            try:
                # Handle matplotlib's unicode minus
                val = float(s.replace('\u2212', '-'))
                res = f"{val:.{self.decimals}f}"
                return res
            except ValueError:
                return s
                
        def get_offset(self):
            s = super().get_offset()
            if not s: return s
            import re
            m = re.search(r'10\^\{.*?(\d+)\}', s)
            if m:
                exponent = m.group(1)
                if '-' in s or '\u2212' in s:
                    return f'$\\times 10^{{-{exponent}}}$'
                else:
                    return f'$\\times 10^{{{exponent}}}$'
            return s
                
    def apply_formatting(ax, scientific=True, decimals=2):
        formatter = CustomScalarFormatter(decimals=decimals, useMathText=True)
        if scientific:
            formatter.set_scientific(True)
            formatter.set_powerlimits((-2, 3))
        else:
            formatter.set_scientific(False)
        ax.yaxis.set_major_formatter(formatter)
    
    axs1[0,0].set_ylabel(r'$c_b\ (\mu\mathrm{M})$')
    # axs1[0,0].set_title('Crosslinker Evolution (Step Stretch)')
    axs1[0,0].grid(True, alpha=0.3)
    apply_formatting(axs1[0,0])
    
    axs1[0,1].set_ylabel(r'$c_f\ (\mu\mathrm{M})$')
    # axs1[0,1].set_title('Crosslinker Evolution (Step Stretch)')
    axs1[0,1].grid(True, alpha=0.3)
    apply_formatting(axs1[0,1])
    
    axs1[1,0].set_ylabel(r'$k_\mathrm{off} (1/\mathrm{s})$')    
    axs1[1,0].set_xlabel('Time (s)')
    axs1[1,0].grid(True, alpha=0.3)
    apply_formatting(axs1[1,0], scientific=False)
    
    axs1[1,1].set_ylabel(r'$k_\mathrm{on} (1/\mathrm{s})$')    
    axs1[1,1].set_xlabel('Time (s)')
    axs1[1,1].grid(True, alpha=0.3)
    axs1[1,1].set_yticks([0.005, 0.006, 0.007])
    axs1[1,1].set_yticklabels(['0.005', '0.006', '0.007'])
    
    # Legend outside for Fig 1
    handles, labels = axs1[0,0].get_legend_handles_labels()
    fig1.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=len(param_values), fontsize=10)
    
    fig1.tight_layout()
    fig1.savefig(os.path.join(sim_dir, f'sensitivity_{param_name}_kinetics.pdf'), bbox_inches='tight')
    plt.close(fig1)
    
    # Format Fig 2
    ax2.set_ylabel('$\\sigma_i^*$',fontsize=16)
    ax2.set_xlabel('Time (s)')
    # ax2.set_title(f'Stress Relaxation Sensitivity: {param_name}')
    ax2.grid(True, alpha=0.3)
    apply_formatting(ax2)
    
    # Legend outside for Fig 2
    # ax2.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), fontsize=10)
    
    fig2.tight_layout()
    fig2.savefig(os.path.join(sim_dir, f'sensitivity_{param_name}_stress.pdf'), bbox_inches='tight')
    plt.close(fig2)
    
    # Restore originals
    k_off_0 = orig_k_off_0
    k_on_0 = orig_k_on_0
    dx_kT = orig_dx_kT
    Rb_max = orig_Rb_max
    cb_max = c_actin * Rb_max
    R_total = orig_R_total
    c_total = c_actin * R_total

if __name__ == '__main__':
    # experiment_constant_stretch()
    # experiment_cyclic_loading()
    
    # Easily change the variable name here to run for other parameters
    # Options: 'k_off_0', 'k_on_0', 'dx_kT', 'Rb_max', 'R_total'
    experiment_sensitivity(param_name='R_total', param_values=[0.02, 0.05, 0.1])
    print("Plots saved in:", sim_dir)
