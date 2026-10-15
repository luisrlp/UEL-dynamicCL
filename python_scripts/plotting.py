import numpy as np                                                       
import matplotlib.pyplot as plt                                          
from scipy.interpolate import griddata                                   
import matplotlib.cm as cm

# PLOTTING THE MICROSPHERE
                                                                             
def plot_microsphere(dirs, cb_array, cb0, title="Bound Crosslinkers ($c_b$)", path=None):                                                      
    """                                                                  
    Plots a 3D microsphere completely filled and colored by the          
    normalized bound crosslinker concentration.                          
    """                                                                  
    fig = plt.figure(figsize=(6, 6))                                     
    ax = fig.add_subplot(111, projection='3d')                           
                                                                            
    # 1. Duplicate directions to the opposite hemisphere                 
    full_dirs = np.vstack([dirs, -dirs])                                 
    full_cb = np.concatenate([cb_array, cb_array])                       
    cb_normalized = full_cb / cb0                                        
                                                                            
    # 2. Create a very fine, continuous spherical grid                   
    u = np.linspace(0, 2 * np.pi, 200)                                   
    v = np.linspace(0, np.pi, 100)                                       
    X = np.outer(np.cos(u), np.sin(v))                                   
    Y = np.outer(np.sin(u), np.sin(v))                                   
    Z = np.outer(np.ones(np.size(u)), np.cos(v))                         
                                                                            
    # 3. Interpolate the sparse 180 points onto the dense grid           
    # 'nearest' creates mathematically perfect equal-area patches!       
    grid_points = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T         
    cb_grid = griddata(full_dirs, cb_normalized, grid_points, method='linear')                                                          
    cb_grid = cb_grid.reshape(X.shape)                                   
                                                                            
    # 4. Map the concentrations to colors                                
    norm = plt.Normalize(vmin=0.0, vmax=1.0)                             
    colors = cm.coolwarm(norm(cb_grid))                                  
                                                                            
    # 5. Plot the fully filled continuous surface                        
    surf = ax.plot_surface(X, Y, Z, facecolors=colors, shade=False,      
antialiased=True)                                                          
                                                                            
    # 6. Formatting                                                      
    m = cm.ScalarMappable(cmap='coolwarm', norm=norm)                    
    m.set_array([])                                                      
    cbar = plt.colorbar(m, ax=ax, shrink=0.5, pad=0.1)                   
    cbar.set_label('$c_b / c_{b0}$')                                     
                                                                            
    ax.set_title(title)                                                  
    
    # Optional: Turn off the grid/axes box so only the sphere floats in  the image
    ax.set_axis_off() 
    
    if path:
        plt.savefig(path, bbox_inches='tight')
    else:
        plt.show()
