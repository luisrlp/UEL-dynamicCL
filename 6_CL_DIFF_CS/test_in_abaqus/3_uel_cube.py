import os
import math

# ==============================================================================
# PARAMETERS
# ==============================================================================
input_file = "base_mesh.inp"

# Deformation Setup
deformation_type = "biaxial"  # Options: "uniaxial", "biaxial", "shear"
stretch_displacement = 0.3     # Displacement magnitude

output_file = f"cube_{deformation_type}_uel_auto.inp"

# Geometry Parameters
cube_size = 2.0  
top_face_y = 2.0 
dummy_element_offset = 100000

# Step Times
t_indent = 1.0;     dtime_indent = 0.001;    max_inc_indent = 0.05
t_hold = 10.0;      dtime_hold = 0.001;      max_inc_hold = 0.5
t_withdraw = 1.0;   dtime_withdraw = 0.001;  max_inc_withdraw = 0.05
t_relax = 10.0;     dtime_relax = 0.001;     max_inc_relax = 0.5

# UEL Parameters
uel_variables = 848
dummy_variables = 108
num_properties = 24

# ==============================================================================
# CHEMICAL BOUNDARY CONFIGURATION
# ==============================================================================
# "closed" -> Impermeable boundary. Fluid/proteins cannot cross this surface.
# "open"   -> Permeable boundary. Fluid/proteins can escape into the surrounding bath.
top_boundary_condition = "closed"   # Top surface (outside the indenter)
side_boundary_condition = "closed"  # Outer physical side walls

bath_chemical_potential = "<INITMU>" # The potential of the surrounding bath


# ==============================================================================
# SCRIPT LOGIC
# ==============================================================================
try:
    with open(input_file, 'r') as f:
        lines = f.readlines()
except FileNotFoundError:
    print(f"Error: {input_file} not found.")
    exit()

new_lines, element_lines = [], []
in_nodes, in_elements = False, False
top_nodes, bottom_nodes = {}, []
xmax_nodes, zmax_nodes, x0_nodes, z0_nodes = [], [], [], []
all_node_ids = []

for line in lines:
    upper_line = line.upper()
    
    if upper_line.startswith("*"):
        in_nodes = in_elements = False
        if upper_line.startswith("*NODE"):
            in_nodes = True
            new_lines.append(line)
            continue
        elif upper_line.startswith("*ELEMENT"):
            in_elements = True
            continue
        elif any(upper_line.startswith(x) for x in ["*SURFACE", "*CONTACT", "*BOUNDARY", "*NSET", "*ELSET", "*AMPLITUDE"]) or "_SURFACE" in upper_line:
            continue
        else:
            new_lines.append(line)
            continue

    if in_nodes:
        parts = line.split(',')
        if len(parts) >= 4:
            n_id = int(parts[0].strip())
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            
            all_node_ids.append(n_id)
            if abs(y - top_face_y) < 1e-6: top_nodes[n_id] = (x, y, z)
            if abs(y - 0.0) < 1e-6: bottom_nodes.append(n_id)
            if abs(x - cube_size) < 1e-6: xmax_nodes.append(n_id)
            if abs(z - cube_size) < 1e-6: zmax_nodes.append(n_id)
            if abs(x - 0.0) < 1e-6: x0_nodes.append(n_id)
            if abs(z - 0.0) < 1e-6: z0_nodes.append(n_id)
                
        new_lines.append(line)
    elif in_elements:
        element_lines.append(line)

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
        parts[0] = str(int(parts[0]) + dummy_element_offset)
        uel_block += ",".join(parts)

def write_nset(name, nodes):
    if not nodes: return ""
    block = f"*Nset, nset={name}\n"
    for i in range(0, len(nodes), 16):
        block += ", ".join(map(str, nodes[i:i+16])) + "\n"
    return block

uel_block += write_nset("all_nodes", all_node_ids)
uel_block += write_nset("top_nodes", list(top_nodes.keys()))
uel_block += write_nset("bottom_nodes", bottom_nodes)
uel_block += write_nset("xmax_nodes", xmax_nodes)
uel_block += write_nset("zmax_nodes", zmax_nodes)
uel_block += write_nset("x0_nodes", x0_nodes)
uel_block += write_nset("z0_nodes", z0_nodes)

open_chem_nodes = set()
if top_boundary_condition == "open": open_chem_nodes.update(list(top_nodes.keys()))
if side_boundary_condition == "open":
    open_chem_nodes.update(x0_nodes + z0_nodes + xmax_nodes + zmax_nodes)
uel_block += write_nset("open_chem_nodes", list(open_chem_nodes))

uel_block += f"""
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

# Dynamic Boundary Condition Generator
def get_fixed_bcs():
    if deformation_type == "uniaxial":
        return "bottom_nodes, 2, 2, 0.0\nx0_nodes, 1, 1, 0.0\nz0_nodes, 3, 3, 0.0\n*Boundary\nextra_element, encastre\nextra_element, 11, 11, 0.0\n"
    elif deformation_type == "biaxial":
        return "bottom_nodes, 2, 2, 0.0\nx0_nodes, 1, 1, 0.0\nz0_nodes, 3, 3, 0.0\n*Boundary\nextra_element, encastre\nextra_element, 11, 11, 0.0\n"
    elif deformation_type == "shear":
        # Isochoric Simple Shear (Plane Strain in Z)
        return "all_nodes, 3, 3, 0.0\nbottom_nodes, 1, 2, 0.0\ntop_nodes, 2, 2, 0.0\n*Boundary\nextra_element, encastre\nextra_element, 11, 11, 0.0\n"
    return ""

def get_active_bcs(disp):
    if deformation_type == "uniaxial":
        return f"top_nodes, 2, 2, {disp}\n"
    elif deformation_type == "biaxial":
        return f"top_nodes, 2, 2, {disp}\nxmax_nodes, 1, 1, {disp}\n"
    elif deformation_type == "shear":
        return f"top_nodes, 1, 1, {disp}\n"
    return ""

def write_step(name, dtime, t, max_inc, disp, is_first=False):
    if t <= 0.0: return ""
    
    bc_str = "*Boundary\n"
    if is_first:
        bc_str += get_fixed_bcs()
    bc_str += get_active_bcs(disp)
    
    step = f"""*Step, name={name}, nlgeom=YES, inc=10000
*Coupled Temperature-displacement, creep=none, deltmx=10.0
{dtime}, {t}, 1e-15, {max_inc}
{bc_str}"""
    if name == "Loading" and len(open_chem_nodes) > 0:
        step += f"open_chem_nodes, 11, 11, {bath_chemical_potential}\n"
    step += """*Output, field, time marks=no
*Node Output, nset=all_nodes
U, NT, RF
*Element Output, elset=dummy_mesh
UVARM, LE
*node output, nset=extra_element
u
*End Step
"""
    return step

uel_block += write_step("Loading", dtime_indent, t_indent, max_inc_indent, stretch_displacement, is_first=True)
uel_block += write_step("Hold", dtime_hold, t_hold, max_inc_hold, stretch_displacement)
uel_block += write_step("Unloading", dtime_withdraw, t_withdraw, max_inc_withdraw, 0.0)
uel_block += write_step("Relax", dtime_relax, t_relax, max_inc_relax, 0.0)

final_lines = []
for line in new_lines:
    if line.upper().startswith("*MATERIAL") or line.upper().startswith("*STEP"): break
    final_lines.append(line if line.strip() else "**\n")

cleaned_uel_block = "\n".join([line if line.strip() else "**" for line in uel_block.split('\n')])

with open(output_file, 'w') as f:
    f.writelines(final_lines)
    f.write(cleaned_uel_block)

print(f"Done! Generated {output_file} for {deformation_type} deformation.")
