import os
import math

# ==============================================================================
# PARAMETERS
# ==============================================================================
input_file = "base_mesh.inp"
output_file = "cube_indent_uel_relax_auto.inp"

# Geometry Parameters
cube_size = 2.0  # Size of the cube mesh (from 0 to 2.0)
top_face_y = 2.0 # The face where the indenter makes contact

# Indenter Location and Geometry
indenter_position = "vertex"  # Options: "vertex" or "center"
indenter_radius = 0.5
initial_gap = 0.25           # Gap between indenter tip and mesh surface before step 1

if indenter_position == "vertex":
    indenter_x = cube_size
    indenter_z = cube_size
elif indenter_position == "center":
    indenter_x = cube_size / 2.0
    indenter_z = cube_size / 2.0
else:
    raise ValueError("indenter_position must be 'vertex' or 'center'")

ref_point_y = top_face_y + initial_gap + indenter_radius

# UEL Parameters
uel_variables = 848
dummy_variables = 106
num_properties = 24

# ID Management (Increase these if you use a very dense mesh)
dummy_element_offset = 100000
ref_node_id = 100000

# ==============================================================================
# CHEMICAL BOUNDARY CONFIGURATION
# ==============================================================================
# "closed" -> Impermeable boundary. Fluid/proteins cannot cross this surface.
# "open"   -> Permeable boundary. Fluid/proteins can escape into the surrounding bath.
top_boundary_condition = "closed"   # Top surface (outside the indenter)
side_boundary_condition = "closed"  # Outer physical side walls

bath_chemical_potential = "<INITMU>" # The potential of the surrounding bath

# ==============================================================================
# STEP CONFIGURATIONS (Set time > 0.0 to enable step, 0.0 to disable)
# ==============================================================================
t_indent = 1.0
depth = -0.5
dtime_indent = 1e-3

t_hold = 10.0
dtime_hold = 1e-4

t_withdraw = 1.0
dtime_withdraw = 1e-4

t_relax = 10.0
dtime_relax = 1e-4

# ==============================================================================
# SCRIPT
# ==============================================================================
print(f"Reading {input_file}...")
try:
    with open(input_file, 'r') as f:
        lines = f.readlines()
except FileNotFoundError:
    print(f"Error: {input_file} not found. Please create it with raw *Node and *Element blocks.")
    exit()

new_lines = []
node_lines = []
element_lines = []
in_nodes = False
in_elements = False
top_nodes = {}
bottom_nodes = []
xmax_nodes = []
zmax_nodes = []
x0_nodes = []
z0_nodes = []
all_node_ids = []

skip_block = False
for line in lines:
    upper_line = line.upper()
    
    # Determine if we should skip this entire keyword block
    if upper_line.startswith("*"):
        if upper_line.startswith("*SURFACE") or upper_line.startswith("*CONTACT") or "_SURFACE" in upper_line:
            skip_block = True
        else:
            skip_block = False
            
    if skip_block:
        continue

    if upper_line.startswith("*NODE"):
        in_nodes = True
        new_lines.append(line)
        continue
    elif upper_line.startswith("*ELEMENT"):
        in_nodes = False
        in_elements = True
        continue
    elif upper_line.startswith("*NSET"):
        in_elements = False

    if in_nodes:
        parts = line.split(',')
        if len(parts) >= 4:
            n_id = int(parts[0].strip())
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            
            all_node_ids.append(n_id)
            if abs(y - top_face_y) < 1e-6:
                top_nodes[n_id] = (x, y, z)
            if abs(y - 0.0) < 1e-6:
                bottom_nodes.append(n_id)
            if abs(x - cube_size) < 1e-6:
                xmax_nodes.append(n_id)
            if abs(z - cube_size) < 1e-6:
                zmax_nodes.append(n_id)
            if abs(x - 0.0) < 1e-6:
                x0_nodes.append(n_id)
            if abs(z - 0.0) < 1e-6:
                z0_nodes.append(n_id)
                
        new_lines.append(line)
        node_lines.append(line)
    elif in_elements:
        element_lines.append(line)
    else:
        new_lines.append(line)

free_flux_nodes = []
for n_id, (x, y, z) in top_nodes.items():
    dist_to_center = math.sqrt((x - indenter_x)**2 + (z - indenter_z)**2)
    if dist_to_center > indenter_radius + 1e-4:
        free_flux_nodes.append(n_id)

uel_block = f"""
*Node, nset=extra_element
99992, 0.0,0.0,0.0
99993, 0.0001,0.0,0.0
99994, 0.0001,0.0001,0.0
99995, 0.0,0.0001,0.0
99996, 0.0,0.0,0.0001
99997, 0.0001,0.0,0.0001
99998, 0.0001,0.0001,0.0001
99999, 0.0,0.0001,0.0001
** ==============================================================================
** 1. UEL DEFINITION (MAIN MESH)
** ==============================================================================
*User Element, type=U3, Nodes=8, Coordinates=3, Properties={num_properties}, Iproperties=2, Variables={uel_variables}, Unsymm
1, 2, 3, 11
*Element, type=U3, elset=main_mesh
"""
for el in element_lines:
    uel_block += el

uel_block += f"""
** EXTRA ELEMENT
*Element,type=C3D8T,elset=extra_element
99999, 99992,99993,99994,99995,99996,99997,99998,99999
** ==============================================================================
** 2. C3D8 DEFINITION (DUMMY MESH FOR VISUALIZATION)
** ==============================================================================
*Element, type=C3D8, elset=dummy_mesh
"""
for el in element_lines:
    parts = el.split(',')
    if len(parts) > 1:
        # Offset dummy elements
        parts[0] = str(int(parts[0]) + dummy_element_offset)
        uel_block += ",".join(parts)

def write_nset(name, nodes):
    # Only write the set if it's not empty
    if not nodes:
        return ""
    block = f"*Nset, nset={name}\n"
    for i in range(0, len(nodes), 16):
        block += ", ".join(map(str, nodes[i:i+16])) + "\n"
    return block

top_dummy_surface_lines = []
for el in element_lines:
    parts = el.split(',')
    if len(parts) == 9:
        el_id = int(parts[0])
        dummy_el_id = el_id + dummy_element_offset
        nodes = [int(n) for n in parts[1:]]
        faces = {
            'S1': [nodes[0], nodes[1], nodes[2], nodes[3]],
            'S2': [nodes[4], nodes[7], nodes[6], nodes[5]],
            'S3': [nodes[0], nodes[4], nodes[5], nodes[1]],
            'S4': [nodes[1], nodes[5], nodes[6], nodes[2]],
            'S5': [nodes[2], nodes[6], nodes[7], nodes[3]],
            'S6': [nodes[3], nodes[7], nodes[4], nodes[0]]
        }
        for face_name, face_nodes in faces.items():
            if all(n in top_nodes for n in face_nodes):
                top_dummy_surface_lines.append(f"{dummy_el_id}, {face_name}\n")

uel_block += f"""
** ==============================================================================
** 3. NODE SETS & CHEMICAL BOUNDARY
** ==============================================================================
"""
uel_block += write_nset("all_nodes", all_node_ids)
if top_boundary_condition == "open":
    uel_block += write_nset("free_flux", free_flux_nodes)
uel_block += write_nset("bottom_nodes", bottom_nodes)
uel_block += write_nset("xmax_nodes", xmax_nodes)
uel_block += write_nset("zmax_nodes", zmax_nodes)
uel_block += write_nset("x0_nodes", x0_nodes)
uel_block += write_nset("z0_nodes", z0_nodes)

uel_block += f"""
** ==============================================================================
** 4. RIGID INDENTER AND CONTACT
** ==============================================================================
*Node, nset=IndenterRef
{ref_node_id}, {indenter_x}, {ref_point_y}, {indenter_z}
*Surface, type=REVOLUTION, name=SURFACE-Probe
********
{indenter_x}, {ref_point_y}, {indenter_z},   {indenter_x}, {ref_point_y + 1.0}, {indenter_z}
START,          {indenter_radius},           0.
 CIRCL,           0.,         -{indenter_radius},           0.,           0.
*Rigid Body, ref node={ref_node_id}, analytical surface=SURFACE-Probe

*Surface, type=ELEMENT, name=SURFACE-Top
"""
for line in top_dummy_surface_lines:
    uel_block += line

uel_block += f"""
*Surface Interaction, name=INTPROP-Frictionless
1.,
*Friction
0.,
*Surface Behavior, augmented Lagrange
*Contact Pair, interaction=INTPROP-Frictionless, type=SURFACE TO SURFACE
SURFACE-Top, SURFACE-Probe

** ==============================================================================
** 5. PROPERTIES AND SIMULATION STEPS
** ==============================================================================
*include, input=properties.inp
*INCLUDE, file=sec_uel_cube.inp
*Solid section, elset=extra_element, material=extra_material
*Solid section, elset=dummy_mesh, material=dummy_material
*Hourglass stiffness
250.0
*Material, name=extra_material
*Elastic
1.e-20
*Conductivity
1.0
*Density
1.0
*Specific heat
1.0
*Material, name=dummy_material
*elastic
1.e-20
*User output variables
{dummy_variables}
** ==============================================================================
** 6. INITIAL CONDITIONS
** ==============================================================================
*Initial conditions, type=temperature
extra_element, 0.0
*Initial conditions, type=temperature
all_nodes, <INITMU>
"""

output_block = """*Output, field, time marks=no
*Node Output, nset=all_nodes
U, NT, RF
*Node Output, nset=IndenterRef
U, RF
*Element Output, elset=dummy_mesh
UVARM, LE
*node output, nset=extra_element
u
*End Step
"""

uel_block += f"""
*Step, name=Indentation, nlgeom=YES, inc=10000
*Coupled Temperature-displacement, creep=none, deltmx=10.0
{dtime_indent}, {t_indent}, 1e-15, 0.02
*Boundary
** Fix the bottom of the gel
bottom_nodes, YSYMM
"""

if indenter_position == "vertex":
    uel_block += """** Apply Quarter Symmetry boundary conditions for vertex indentation
xmax_nodes, XSYMM
zmax_nodes, ZSYMM
** Extra Element Fix
*Boundary
extra_element, encastre
extra_element, 11, 11, 0.0
"""

uel_block += f"""** Indenter Constraints
IndenterRef, 1, 1
IndenterRef, 3, 3
IndenterRef, 4, 4
IndenterRef, 5, 5
IndenterRef, 6, 6
IndenterRef, 2, 2, {depth}
"""

if top_boundary_condition == "open":
    uel_block += f"""** Chemical Boundary (Open Top)
free_flux, 11, 11, {bath_chemical_potential}
"""

if side_boundary_condition == "open":
    uel_block += f"""** Chemical Boundary (Open Sides)
x0_nodes, 11, 11, {bath_chemical_potential}
z0_nodes, 11, 11, {bath_chemical_potential}
"""
    if indenter_position == "center":
        uel_block += f"""xmax_nodes, 11, 11, {bath_chemical_potential}
zmax_nodes, 11, 11, {bath_chemical_potential}
"""

uel_block += output_block

if t_hold > 0.0:
    uel_block += f"""
*Step, name=Hold, nlgeom=YES, inc=10000
*Coupled Temperature-displacement, creep=none, deltmx=10.0
{dtime_hold}, {t_hold}, 1e-15, 0.5
*Boundary
IndenterRef, 2, 2, {depth}
""" + output_block

if t_withdraw > 0.0:
    uel_block += f"""
*Step, name=Withdraw, nlgeom=YES, inc=10000
*Coupled Temperature-displacement, creep=none, deltmx=10.0
{dtime_withdraw}, {t_withdraw}, 1e-15, 0.02
*Boundary
IndenterRef, 2, 2, 0.0
""" + output_block

if t_relax > 0.0:
    uel_block += f"""
*Step, name=Relax, nlgeom=YES, inc=10000
*Coupled Temperature-displacement, creep=none, deltmx=10.0
{dtime_relax}, {t_relax}, 1e-15, 0.5
*Boundary
IndenterRef, 2, 2, 0.0
""" + output_block

final_lines = []
for line in new_lines:
    if line.upper().startswith("*MATERIAL") or line.upper().startswith("*STEP"):
        break
    if line.strip() == "":
        final_lines.append("**\n")
    else:
        final_lines.append(line)

# Clean up blank lines (Abaqus pre-processor crashes on blank lines)
cleaned_uel_block = ""
for line in uel_block.split('\n'):
    if line.strip() == "":
        cleaned_uel_block += "**\n"
    else:
        cleaned_uel_block += line + "\n"

print(f"Writing {output_file}...")
with open(output_file, 'w') as f:
    f.writelines(final_lines)
    f.write(cleaned_uel_block)

print("Done! You can now run the job with 'abaqus job=cube_indent_uel_relax_auto.inp'")
