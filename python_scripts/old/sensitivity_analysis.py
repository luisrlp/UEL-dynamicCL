import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add the directory to sys.path so we can import analytical_kinetics
scripts_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(scripts_dir)

import analytical_kinetics as ak

# Disable TeX to avoid missing latex errors on user's machine
ak.plt.rcParams['text.usetex'] = False

def run_relaxation(k_off_0, K_eq, dx_kT, Rb_max):
    # Set the globals in the module
    ak.k_off_0 = k_off_0
    ak.K_eq = K_eq
    ak.dx_kT = dx_kT
    ak.Rb_max = Rb_max
    
    # Update dependent globals
    ak.cb_max = ak.c_actin * ak.Rb_max
    
    # Run the equilibrium to get cb0
    cb0 = ak.get_equilibrium_cb0()
    
    # Setup time array and stretch function (Step to 1.5 over 1s, hold for 60s)
    t_end = 60.0
    t_eval = np.linspace(0, t_end, 500)
    
    def stretch_func(t):
        t_ramp = 1.0
        stretch_max = 1.5
        if t < t_ramp:
            return 1.0 + (stretch_max - 1.0) * (t / t_ramp)
        return stretch_max
        
    # Solve ODE
    sol = ak.solve_ivp(ak.reaction_rate, [0, t_end], [cb0], 
                       args=(stretch_func,), max_step=0.1, t_eval=t_eval)
                       
    # Calculate macroscopic stress
    stresses = []
    for i in range(len(sol.t)):
        cb = sol.y[0][i]
        stretch = stretch_func(sol.t[i])
        force, n_scaling, _ = ak.get_geometry_and_force(cb, stretch)
        
        # Nominal stress proportional to network density * force
        # n_scaling acts as 1/L for the network chain density
        stress = cb * force * n_scaling * 1e12 # scaled for visibility
        stresses.append(stress)
        
    return sol.t, np.array(stresses)

# Baseline parameters
baseline = {
    'k_off_0': 0.05,
    'K_eq': 0.25,
    'dx_kT': 24.3, # Equivalent to chi ~ 5.0 in some scaling, baseline is usually 0.1 / kT = 24.3
    'Rb_max': 0.25
}

# Define variations (Low, Baseline, High)
variations = {
    'k_off_0': [0.005, 0.05, 0.5],
    'K_eq': [0.05, 0.25, 1.25],
    'dx_kT': [10.0, 24.3, 50.0],
    'Rb_max': [0.1, 0.25, 0.5]
}

titles = {
    'k_off_0': 'Baseline Unbinding Rate ($k_{off,0}$)',
    'K_eq': 'Equilibrium Constant ($K_{eq} \propto k_{on}$)',
    'dx_kT': 'Stretch Sensitivity Parameter ($\Delta x / k_B T$)',
    'Rb_max': 'Bound Saturation Limit ($R_{b,max}$)'
}

fig, axs = plt.subplots(2, 2, figsize=(14, 10))
axs = axs.flatten()

for i, param in enumerate(variations.keys()):
    ax = axs[i]
    ax.set_title(titles[param], fontsize=14, fontweight='bold')
    ax.set_xlabel('Time (s)', fontsize=12)
    ax.set_ylabel('Macroscopic Stress (A.U.)', fontsize=12)
    
    for val in variations[param]:
        # Reset to baseline
        params = baseline.copy()
        params[param] = val
        
        # Run test
        t, stress = run_relaxation(params['k_off_0'], params['K_eq'], params['dx_kT'], params['Rb_max'])
        
        ax.plot(t, stress, linewidth=2, label=f'{param} = {val:.3g}')
        
    ax.legend(fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
out_path = os.path.join(scripts_dir, 'sensitivity_results.png')
plt.savefig(out_path, dpi=300)
print(f'Sensitivity analysis complete. Plot saved to {out_path}')
