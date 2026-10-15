import numpy as np                                                       
import matplotlib.pyplot as plt                                          
from scipy.interpolate import griddata                                   
import matplotlib.cm as cm

# PLOTTING THE MICROSPHERE

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
                                                                            
def subdivide(faces):                                                    
    new_faces = []                                                       
    for (v1, v2, v3) in faces:                                           
        m12 = v1 + v2                                                    
        m12 /= np.linalg.norm(m12)                                       
        m23 = v2 + v3                                                    
        m23 /= np.linalg.norm(m23)                                       
        m31 = v3 + v1                                                    
        m31 /= np.linalg.norm(m31)                                       
                                                                            
        new_faces.append([v1, m12, m31])                                 
        new_faces.append([v2, m23, m12])                                 
        new_faces.append([v3, m31, m23])                                 
        new_faces.append([m12, m23, m31])                                
    return new_faces                                                     
                                                                            
def create_icosahedron():                                                
    phi = (1.0 + np.sqrt(5.0)) / 2.0                                     
    verts = np.array([                                                   
        [-1,  phi, 0], [ 1,  phi, 0], [-1, -phi, 0], [ 1, -phi, 0],      
        [ 0, -1,  phi], [ 0,  1,  phi], [ 0, -1, -phi], [ 0,  1, -phi],  
        [ phi, 0, -1], [ phi, 0,  1], [-phi, 0, -1], [-phi, 0,  1]       
    ])                                                                   
    verts /= np.linalg.norm(verts, axis=1)[:, np.newaxis]                
                                                                            
    faces = [                                                            
        [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],       
        [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],       
        [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],           
        [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]          
    ]                                                                    
    return verts[faces]                                                  
                                                                            
def plot_microsphere_mathematica_style(dirs, cb_array, cb0, title="Bound Crosslinkers ($c_b$)"):                                                    
    fig = plt.figure(figsize=(8, 8))                                     
    ax = fig.add_subplot(111, projection='3d')                           
                                                                         
    # 1. Prepare Data                                                    
    full_dirs = np.vstack([dirs, -dirs])                                 
    full_cb = np.concatenate([cb_array, cb_array])                       
    cb_normalized = full_cb / cb0                                        
                                                                         
    # 2. Build a dense Icosahedral Triangular Mesh
    faces = create_icosahedron()                                         
    for _ in range(4):                                                   
        faces = subdivide(faces)                                         
                                                                         
    verts, triangles, v_dict = [], [], {}                                
    for face in faces:                                                   
        tri = []                                                         
        for pt in face:                                                  
            k = tuple(np.round(pt, 5))                                   
            if k not in v_dict:                                          
                v_dict[k] = len(verts)                                   
                verts.append(pt)                                         
            tri.append(v_dict[k])                                        
        triangles.append(tri)                                            
                                                                         
    verts = np.array(verts)                                              
    triangles = np.array(triangles)                                      
    X, Y, Z = verts[:, 0], verts[:, 1], verts[:, 2]                      
                                                                         
    # 3. Interpolate colors directly onto the triangle vertices          
    cb_verts = griddata(full_dirs, cb_normalized, verts, method='nearest')                                                          
    cb_faces = np.mean(cb_verts[triangles], axis=1)
                                                                         
    norm = plt.Normalize(vmin=0.0, vmax=1.0)                             
                                                                         
    # 4. Plot using Triangles! 
    # --> Notice the coordinate swap here: (Z, X, Y)
    # This places the original Y-data on Matplotlib's vertical Z-axis
    surf = ax.plot_trisurf(Z, X, Y, triangles=triangles, cmap='coolwarm',
                           alpha=0.9, antialiased=True, vmin=0.0, vmax=1.0,
                           shade=False)                                                                         
    surf.set_array(cb_faces) 
                                                                         
    # 5. Add the geometric wireframe on top
    wire_faces = create_icosahedron()                                    
    for _ in range(3):                                                   
        wire_faces = subdivide(wire_faces)                               
                                                                         
    edges = set()                                                        
    for face in wire_faces:                                              
        for i in range(3):                                               
            v1, v2 = tuple(face[i]), tuple(face[(i+1)%3])                
            if v1 > v2: v1, v2 = v2, v1                                  
            edges.add((v1, v2))                                          
                                                                         
    for (v1, v2) in edges:
        v1, v2 = np.array(v1), np.array(v2)
        l_x, l_y, l_z = zip(v1 * 1.01, v2 * 1.01)
        # --> Notice the coordinate swap here too: (l_z, l_x, l_y)
        ax.plot(l_z, l_x, l_y, color='black', linewidth=0.5, alpha=0.8)  
    
    # 6. Formatting
    m = cm.ScalarMappable(cmap='coolwarm', norm=norm)
    m.set_array([])
    cbar = plt.colorbar(m, ax=ax, shrink=0.5, pad=0.1)
    cbar.set_label('$c_b / c_{b0}$')
    
    # Ensure proportional scaling so the sphere doesn't look stretched
    ax.set_box_aspect([1, 1, 1]) 
    
    # ax.set_title(title)
    ax.set_axis_off() 
    plt.tight_layout()
    # plt.savefig("figures/microsphere_mathematica_style10s.svg", bbox_inches='tight')
    # plt.show()