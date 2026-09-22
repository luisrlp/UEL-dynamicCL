import os

# ==============================================================================
# PARAMETERS
# ==============================================================================
input_file = "base_mesh.inp"
output_file = "cube_indent_uel_relax_auto.inp"

# Indenter Geometry (must match coordinates of top face)
cube_size = 2.0
indenter_radius = 0.5
indenter_center_x = 1.0  # Assumes indentation happens at the center of the top face
indenter_center_y = 1.0

# Step and Solver Parameters
t_load = 1.0
t_relax = 10.0
depth = -0.5
initial_dtime = 1e-4

# UEL Parameters
n_state_vars = 56 # 7 SDVs * 8 Integration points
num_properties = 14

# ==============================================================================
# SCRIPT
# ==============================================================================
print(f"Reading {input_file}...")
with open(input_file, 'r') as f:
    lines = f.readlines()

new_lines = []
node_lines = []
element_lines = []
in_nodes = False
in_elements = False
top_nodes = {} # node_id : (x, y, z)

# 1. Parse Nodes and Elements
for line in lines:
    if line.upper().startswith("*NODE"):
        in_nodes = True
        new_lines.append(line)
        continue
    elif line.upper().startswith("*ELEMENT"):
        in_nodes = False
        in_elements = True
        continue
    elif line.upper().startswith("*NSET") or line.upper().startswith("*SURFACE"):
        in_elements = False

    if in_nodes:
        # Extract coordinates to build top_nodes dict
        parts = line.split(',')
        if len(parts) >= 4:
            n_id = int(parts[0].strip())
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            # If node is on the top face
            if abs(z - cube_size) < 1e-6:
                top_nodes[n_id] = (x, y, z)
        new_lines.append(line)
        node_lines.append(line)
    elif in_elements:
        element_lines.append(line)
    else:
        # Append everything else (sets, surfaces, etc)
        new_lines.append(line)

# 2. Compute free_flux chemical boundary nodes
free_flux_nodes = []
for n_id, (x, y, z) in top_nodes.items():
    dist_to_center = ((x - indenter_center_x)**2 + (y - indenter_center_y)**2)**0.5
    if dist_to_center >= indenter_radius:
        free_flux_nodes.append(n_id)

# 3. Inject the UEL definition
uel_block = f"""
** ==============================================================================
** 1. UEL DEFINITION (MAIN MESH)
** ==============================================================================
*User Element, type=U3, Nodes=8, Coordinates=3, Properties={num_properties}, Iproperties=2, Variables={n_state_vars}, Unsymm
1, 2, 3, 11
*Element, type=U3, elset=main_mesh
"""
# Add the main mesh connectivity
for el in element_lines:
    uel_block += el

uel_block += f"""
** ==============================================================================
** 2. C3D8 DEFINITION (DUMMY MESH FOR VISUALIZATION)
** ==============================================================================
** Note: Dummy elements are offset by +10000 to avoid ID collisions
*Element, type=C3D8, elset=dummy_mesh
"""
# Add the dummy mesh connectivity (ID + 10000)
for el in element_lines:
    parts = el.split(',')
    if len(parts) > 1:
        parts[0] = str(int(parts[0]) + 10000)
        uel_block += ",".join(parts)

uel_block += f"""
** ==============================================================================
** 3. CHEMICAL BOUNDARY SET (FREE FLUX)
** ==============================================================================
*Nset, nset=free_flux
"""
# Chunk free_flux_nodes into 16 per line (Abaqus limit)
for i in range(0, len(free_flux_nodes), 16):
    uel_block += ", ".join(map(str, free_flux_nodes[i:i+16])) + "\n"

uel_block += f"""
** ==============================================================================
** 4. RIGID INDENTER DEFINITION
** ==============================================================================
*Node, nset=IndenterRef
99999, {indenter_center_x}, {indenter_center_y}, {cube_size}
*Surface, type=REVOLUTION, name=SURFACE-Probe
********
{indenter_center_x}, {indenter_center_y+0.75}, {cube_size},   {indenter_center_x}, {indenter_center_y+1.75}, {cube_size}
START,          {indenter_radius},           0.
 CIRCL,           0.,         -{indenter_radius},           0.,           0.
*Rigid Body, ref node=99999, analytical surface=SURFACE-Probe
*Surface Behavior, augmented Lagrange
*Contact Pair, interaction=INTPROP-Frictionless, type=SURFACE TO SURFACE
SURFACE-Top, SURFACE-Probe

** ==============================================================================
** 5. PROPERTIES AND SIMULATION STEPS
** ==============================================================================
*include, input=properties.inp

*Step, name=Indentation, nlgeom=YES, inc=10000
*Coupled Temperature-displacement, creep=none, deltmx=10.0
{initial_dtime}, {t_load}, 1e-15, 0.1
*Boundary
bottom_nodes, 1, 3, 0.0
IndenterRef, 1, 2, 0.0
IndenterRef, 3, 3, {depth}
free_flux, 11, 11, 0.0
*Output, field, frequency=10
*Node Output
U, NT
*Element Output
SDV,
*End Step

*Step, name=Relaxation, nlgeom=YES, inc=10000
*Coupled Temperature-displacement, creep=none, deltmx=10.0
{initial_dtime}, {t_relax}, 1e-15, 0.5
*Boundary
bottom_nodes, 1, 3, 0.0
IndenterRef, 1, 2, 0.0
IndenterRef, 3, 3, {depth}
free_flux, 11, 11, 0.0
*Output, field, frequency=10
*Node Output
U, NT
*Element Output
SDV,
*End Step
"""

# Find where to inject the UEL block (right after nodes, before materials/steps)
# We will just write a new file completely replacing everything below *Nset definitions
final_lines = []
for line in new_lines:
    if line.upper().startswith("*MATERIAL") or line.upper().startswith("*STEP"):
        break
    final_lines.append(line)

# Write to file
print(f"Writing {output_file}...")
with open(output_file, 'w') as f:
    f.writelines(final_lines)
    f.write(uel_block)

print("Done! You can now run the job with 'abaqus job=cube_indent_uel_relax_auto.inp'")
