import numpy as np

class Kinematics:
    """
    A class to compute and store all kinematic tensors and invariants
    from a given deformation gradient F.
    """
    def __init__(self, F):
        # Deformation gradient
        self.F = np.array(F, dtype=float)
        
        # Identity tensors
        self.I2 = np.eye(3)
        self.I4 = np.einsum('ik,jl->ijkl', self.I2, self.I2) # I_ijkl = delta_ik * delta_jl
        self.I4S = 0.5 * (self.I4 + np.einsum('il,jk->ijkl', self.I2, self.I2))
        
        # Determinant
        self.J = np.linalg.det(self.F)
        if self.J <= 0:
            raise ValueError(f"Non-positive Jacobian: J = {self.J}")
            
        # Isochoric part of the deformation gradient
        self.F_iso = (self.J ** (-1.0 / 3.0)) * self.F
        
        # Inverses
        self.F_inv = np.linalg.inv(self.F)
        self.F_iso_inv = np.linalg.inv(self.F_iso)
        
        # Right Cauchy-Green
        self.C = self.F.T @ self.F
        self.C_inv = np.linalg.inv(self.C)
        
        # Isochoric right Cauchy-Green
        self.C_iso = self.F_iso.T @ self.F_iso
        self.C_iso_inv = np.linalg.inv(self.C_iso)
        
        # Left Cauchy-Green
        self.b_cg = self.F @ self.F.T
        
        # Isochoric left Cauchy-Green
        self.b_cg_iso = self.F_iso @ self.F_iso.T
        
        # Invariants of the isochoric right Cauchy-Green tensor
        self.I1_bar = np.trace(self.C_iso)
        self.I2_bar = 0.5 * (self.I1_bar ** 2 - np.trace(self.C_iso @ self.C_iso))
        
        # Stretch tensors via eigenvalue decomposition (symmetric matrices)
        eigvl_Cbar, eigvc_Cbar = np.linalg.eigh(self.C_iso)
        eigvl_bbar, eigvc_bbar = np.linalg.eigh(self.b_cg_iso)
        
        # Right stretch tensor (U)
        # We use np.abs() inside sqrt to protect against numerical noise causing tiny negative eigenvalues
        self.U_bar = eigvc_Cbar @ np.diag(np.sqrt(np.abs(eigvl_Cbar))) @ eigvc_Cbar.T
        
        # Left stretch tensor (V)
        # Note: In your notebook this was called 'b_bar', but mathematically it is V_bar.
        # I added the np.sqrt() here which was missing in the notebook's reconstruction!
        self.V_bar = eigvc_bbar @ np.diag(np.sqrt(np.abs(eigvl_bbar))) @ eigvc_bbar.T
        
        # Rotation tensor
        self.R = self.F_iso @ np.linalg.inv(self.U_bar)
        
        # Deviatoric projection tensors
        self.PROJEUL = self.I4S - (1.0 / 3.0) * np.einsum('ij,kl->ijkl', self.I2, self.I2)
        self.PROJLAG = self.I4 - (1.0 / 3.0) * np.einsum('ij,kl->ijkl', self.C_inv, self.C)

if __name__ == "__main__":
    # --- QUICK TEST ---
    F_test = np.array([[1.0, 0.2, 0.06],
                       [0.15, 1.1, 0.2],
                       [0.1, 0.05, 1.0]])
    
    kin = Kinematics(F_test)
    
    print(f"Jacobian (J) = {kin.J:.6f}")
    print("\nIsochoric Right Cauchy-Green (C_iso):")
    print(np.round(kin.C_iso, 4))
    print(f"\nFirst Invariant (I1_bar) = {kin.I1_bar:.6f}")
    print(f"\nShape of Eulerian Projection Tensor (PROJEUL) = {kin.PROJEUL.shape}")
