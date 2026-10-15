import numpy as np
from scipy.special import erfi

### Check consistency with the dfd notebook

def generate_icosahedron(factor=1):
    """
    Generates integration directions and weights for a subdivided icosahedron.
    Only one hemisphere (10 faces) is integrated due to filament symmetry.
    
    Args:
        factor (int): The subdivision factor (NDIR = 10 * factor^2)
        
    Returns:
        directions: (N, 3) numpy array of unit vectors (centroids of subtriangles)
        weights: (N,) numpy array of integration weights (spherical area / 4pi)
    """
    # 1. Base Icosahedron Vertices
    phi = 0.5 * (np.sqrt(5.0) + 1.0)
    a = phi / np.sqrt(1.0 + phi**2)
    b = 1.0 / np.sqrt(1.0 + phi**2)
    z = 0.0
    
    verts = np.array([
        [ a,  b,  z], [ a, -b,  z], [ b,  z,  a], [ b,  z, -a],
        [ z,  a,  b], [ z,  a, -b], [ z, -a,  b], [ z, -a, -b],
        [-b,  z,  a], [-b,  z, -a], [-a,  b,  z], [-a, -b,  z]
    ])
    
    # 2. Base Icosahedron Faces (10 faces for one hemisphere)
    # Matches the Fortran face_point array (0-indexed for Python)
    faces = np.array([
        [0, 1, 3], [0, 2, 1], [0, 3, 5], [0, 4, 2], [0, 5, 4],
        [1, 2, 6], [1, 6, 7], [1, 7, 3], [2, 4, 8], [2, 8, 6]
    ])
    
    directions = []
    weights = []
    
    # Helper: Project flat planar point to unit sphere
    def project(pt):
        return pt / np.linalg.norm(pt)
        
    # Helper: Spherical triangle area (Girard's formula)
    def spherical_area(v1, v2, v3):
        # Cosines of the sides (clipped for numerical stability)
        ca = np.clip(np.dot(v2, v3), -1.0, 1.0)
        cb = np.clip(np.dot(v3, v1), -1.0, 1.0)
        cc = np.clip(np.dot(v1, v2), -1.0, 1.0)
        
        sa, sb, sc = np.arccos(ca), np.arccos(cb), np.arccos(cc)
        ss = (sa + sb + sc) / 2.0
        
        # Spherical angles using half-angle formulas
        # The max(0, ...) protects against floating point negative zeros
        tan_A2 = np.sqrt(max(0.0, np.sin(ss - sb)*np.sin(ss - sc) / (np.sin(ss)*np.sin(ss - sa))))
        tan_B2 = np.sqrt(max(0.0, np.sin(ss - sa)*np.sin(ss - sc) / (np.sin(ss)*np.sin(ss - sb))))
        tan_C2 = np.sqrt(max(0.0, np.sin(ss - sa)*np.sin(ss - sb) / (np.sin(ss)*np.sin(ss - sc))))
        
        A = 2.0 * np.arctan(tan_A2)
        B = 2.0 * np.arctan(tan_B2)
        C = 2.0 * np.arctan(tan_C2)
        
        return A + B + C - np.pi
    
    # 3. Subdivide each face using Barycentric coordinates
    for face in faces:
        vA, vB, vC = verts[face[0]], verts[face[1]], verts[face[2]]
        
        for i in range(factor):
            for j in range(factor - i):
                k = factor - i - j
                
                # Subtriangle Type 1 (pointing "up")
                p1, p2, p3 = (i+1, j, k-1), (i, j+1, k-1), (i, j, k)
                
                c1 = project((p1[0]*vA + p1[1]*vB + p1[2]*vC) / factor)
                c2 = project((p2[0]*vA + p2[1]*vB + p2[2]*vC) / factor)
                c3 = project((p3[0]*vA + p3[1]*vB + p3[2]*vC) / factor)
                
                centroid = project(c1 + c2 + c3)
                area = spherical_area(c1, c2, c3)
                
                directions.append(centroid)
                weights.append(area / (4.0 * np.pi))
                
                # Subtriangle Type 2 (pointing "down")
                if k - 2 >= 0:
                    p1, p2, p3 = (i+1, j, k-1), (i+1, j+1, k-2), (i, j+1, k-1)
                    
                    c1 = project((p1[0]*vA + p1[1]*vB + p1[2]*vC) / factor)
                    c2 = project((p2[0]*vA + p2[1]*vB + p2[2]*vC) / factor)
                    c3 = project((p3[0]*vA + p3[1]*vB + p3[2]*vC) / factor)
                    
                    centroid = project(c1 + c2 + c3)
                    area = spherical_area(c1, c2, c3)
                    
                    directions.append(centroid)
                    weights.append(area / (4.0 * np.pi))
                    
    return np.array(directions), np.array(weights)

def mises_density(b, cos_theta):
    """
    Computes the Von Mises-Fisher density for a given dispersion parameter b and cos(theta).
    
    Args:
        b: dispersion parameter (b > 0). If 0, perfectly isotropic.
        cos_theta: cosine of the angle between the direction and preferred direction.

    Returns:
        rho: the Von Mises-Fisher density.
    """
    if b <= 0:
        raise ValueError("Dispersion parameter 'b' must be positive.")

    aux1 = np.sqrt(b / (2.0 * np.pi))
    aux2 = np.exp(2.0 * b * cos_theta**2)
    efi = erfi(np.sqrt(2.0 * b))

    rho = 4.0 * aux1 * aux2 / efi
    return rho

def apply_dispersion(dirs, areas, pref_dir, b_disp):
    """
    Applies the Discrete Fibre Dispersion (DFD) model to the integration sphere.
    
    Args:
        dirs: (N, 3) numpy array of integration directions.
        areas: (N,) numpy array of isotropic weights (A_i / 4pi).
        pref_dir: (3,) array or list representing the preferred direction.
        b_disp: dispersion parameter (b > 0). If 0, perfectly isotropic.
        
    Returns:
        new_weights: (N,) numpy array of anisotropic weights. 
                     These sum to exactly 1.0 over the hemisphere!
    """
    pref_dir = np.array(pref_dir, dtype=float)
    pref_dir = pref_dir / np.linalg.norm(pref_dir)
    
    if b_disp <= 1e-10:
        # Isotropic case: Multiply by 2 so the hemisphere sum is 1.0
        return areas * 2.0
        
    # Calculate cos(theta) for each direction
    cos_theta = np.dot(dirs, pref_dir)
    
    # Calculate the Von Mises-Fisher density (matching density.f90)
    rho = mises_density(b_disp, cos_theta)
    
    # Normalization (Divide by 2*pi according to Li et al. 2018)
    rho = rho / (2.0 * np.pi)
    
    # The 'areas' array from icos.py contains (A_i / 4pi)
    # We extract the true geometric surface area A_i:
    actual_triangle_areas = areas * 4.0 * np.pi
    
    # The final integration weight for each direction is rho_i * A_i
    new_weights = rho * actual_triangle_areas
    
    return new_weights

if __name__ == "__main__":
    # Quick test
    dirs, areas = generate_icosahedron(factor=3)

    print(f"Number of (unique) directions: {len(dirs)}")
    print(f"Isotropic sum (areas * 2): {np.sum(areas * 2.0):.6f}")
    
    pref_dir = [1.0, 0.0, 0.0]
    b_disp = 0.1
    
    new_weights = apply_dispersion(dirs, areas, pref_dir, b_disp)
    print(f"Number of weights: {len(new_weights)}")
    print(f"Anisotropic sum (b={b_disp}): {np.sum(new_weights):.6f}")