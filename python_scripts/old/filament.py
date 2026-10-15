import numpy as np
import pandas as pd
import os
from scipy.optimize import brentq

class Filament:
    def __init__(self, def_info, sim_dir):
        self.sim_dir = sim_dir
        self.def_initial = def_info['def_initial']
        self.def_max = def_info['def_max']
        self.increments = def_info['increments']

    def inextensible(self, mat_props):
        """Returns the force for an inextensible filament"""
        a = mat_props['L'] / mat_props['R0F']
        mat_props['B0'] = 294.0 * 1.38065e-5 * mat_props['Lp']
        aux = mat_props['L']**2 / (np.pi**2 * mat_props['B0'])
        force_array = np.zeros((self.increments))
        def_array = np.linspace(self.def_initial, self.def_max, self.increments)
        force_array = 1/aux * ( ((a-1)/(a - def_array)) ** (1/mat_props['BETA']) - 1 )
        # force_array = 2 * def_array + 5 * np.exp(a) + mat_props['B0']**mat_props['BETA'] * mat_props['MU0'] / mat_props['L']**2 
        f_star = force_array * mat_props['L']**2 / (np.pi**2*mat_props['B0'])
        return {'force': np.array(force_array),
            'deformation': np.array(def_array),
            'f_star': np.array(f_star)
            }
    
    @staticmethod
    def G(force: float, stretch: float, mat_props: dict) -> float:
        """
        Computes the difference between the left-hand side (lhs) and right-hand side (rhs)
        of the force equation for an extensible filament.
        
        Parameters:
        - force (float): Force to evaluate.
        - stretch (float): Stretch value.
        - mat_props (dict): Material properties containing keys:
            LAMBDA0, R0F, L, MU0, Lp, BETA, #B0
        
        Returns:
        - float: The result of the equation (lhs - rhs).
        """
        lhs = (stretch * mat_props['LAMBDA0'] * mat_props['R0F']) / mat_props['L']
        aux1 = force / mat_props['MU0']
        aux2 = mat_props['L']**2 / (np.pi**2 * mat_props['B0'])
        numerator = (1 + 2 * aux1) * (1 + aux1)**mat_props['BETA'] * (1 - mat_props['R0F'] / mat_props['L'])
        denominator = (1 + force * aux2 + force * aux1 * aux2)**mat_props['BETA']
        rhs = 1 + aux1 - numerator / denominator
        return lhs - rhs
    
    def fil_force(self, mat_props: dict, dw: bool = False, ddw: bool = False, cl: bool = True) -> dict:
        """
        Computes the force for an extensible filament (with compliant crosslinker), 
        given material properties.

        Parameters:
        - mat_props (dict): A dictionary containing material properties. Expected keys:
            LAMBDA0, R0F, L, MU0, B0, BETA, R0C, ETA
        - dw (bool): Flag indicating whether to compute the first derivative of the strain energy.
        - ddw (bool): Flag indicating whether to compute the second derivative.
        - cl (bool): Flag indicating whether to include the compliant crosslinker effect.

        Returns:
        - dict: A dictionary containing:
            - 'force': Array of computed forces.
            - 'deformation': Array of corresponding stretches.
            - 'lambda_f': Array of filament stretch values.
        """
        # Validate input
        required_keys = ['LAMBDA0', 'R0F', 'L', 'MU0', 'Lp', 'BETA', 'R0C', 'ETA'] #,'B0']
        for key in required_keys:
            if key not in mat_props:
                raise ValueError(f"Missing required material property: {key}")

        # Generate deformation and force arrays
        if 'STRETCH' in mat_props.keys():
            # Dummy deformation array for consistency with pce pipeline (requires a deformation array)
            def_array = np.array([mat_props['STRETCH']])
            force_array = np.array([0.0])
            lambdaf_array = np.array([0.0])
        else:
            def_array = np.linspace(self.def_initial, self.def_max, self.increments)
            force_array = np.zeros((self.increments))
            lambdaf_array = np.zeros((self.increments))
        mat_props['B0'] = 298.0 * 1.38065e-5 * mat_props['Lp']
        if cl:
            R0 = mat_props['R0F'] + mat_props['R0C']
            LAMBDA0F = mat_props['ETA'] * (R0/mat_props['R0F']) * (mat_props['LAMBDA0']-1) + 1
            # Compute forces
            for i, stretch in enumerate(def_array):
                if 0 < mat_props['ETA'] <= 1:
                    stretch_pre = stretch
                    stretch = mat_props['ETA'] * (R0/mat_props['R0F']) * (stretch-1) + 1
                try:
                    lambdaf_array[i] = stretch
                    force_array[i] = brentq(
                        self.G, a=0.0, b=1e20, 
                        args=(stretch, mat_props), xtol=2e-12, maxiter=100
                    )
                except ValueError:
                    r0f = mat_props['R0F']
                    l = mat_props['L']
                    print(f"Root not found for r0f {r0f}, l {l}, stretch {stretch}, stretch_pre {stretch_pre}")
                    force_array[i] = np.nan
        else:
            R0 = mat_props['R0F']
            for i, stretch in enumerate(def_array):
                try:
                    lambdaf_array[i] = stretch
                    force_array[i] = brentq(
                        self.G, a=0.0, b=1e20, 
                        args=(stretch, mat_props), xtol=2e-12, maxiter=100
                    )
                except ValueError:
                    print(f"Root not found for stretch {stretch}")
                    force_array[i] = np.nan
        # Get results dictionary
        f_star = force_array * mat_props['L']**2 / (np.pi**2*mat_props['B0'])
        result = {
            'force': np.array(force_array),
            'deformation': np.array(def_array),
            'f_star': np.array(f_star),
            # 'lambda_f': np.array(lambdaf_array),
        }

        # If dw or ddw are True, compute derivatives
        if dw:
            ###OLD
            # dw_array = mat_props['LAMBDA0'] * mat_props['R0F'] * force_array
            dw_array = mat_props['LAMBDA0'] * R0 * force_array
            result['dw'] = np.array(dw_array)
        if ddw:
            alpha = np.pi**2 * mat_props['B0'] / (mat_props['L']**2 * mat_props['MU0'])
            f_star = force_array * mat_props['L']**2 / (np.pi**2*mat_props['B0'])
            ### OLD
            # num = mat_props['LAMBDA0']**2 * mat_props['R0F']**2 * mat_props['MU0'] / mat_props['L']
            num = mat_props['ETA'] * LAMBDA0F * mat_props['LAMBDA0'] * R0**2 * mat_props['MU0'] / mat_props['L']
            Y = mat_props['BETA']/alpha * (1 + 2*f_star*alpha)**2 / (1 + f_star + alpha * f_star**2) - mat_props['BETA'] * (1+2*alpha*f_star)/(1+alpha*f_star) - 2
            den = 1 + Y * ( (1+alpha*f_star) / (1+f_star+alpha*f_star**2) )**mat_props['BETA'] * (1-mat_props['R0F']/mat_props['L'])
            ddw_array = num / den
            result['ddw'] = np.array(ddw_array)
        
        return result
    

    def save_results(self, results: pd.DataFrame):
        """Save results to csv file"""
        # Save dataframe to pkl
        results.to_pickle(os.path.join(self.sim_dir, 'results.pkl'))
