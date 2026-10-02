import os
import ast
import sys
import math

# ==============================================================================
# PARAMETERS (defaults; the presets below may override them)
# ==============================================================================
input_file = "base_mesh.inp"
output_file = "cube_indent_uel_relax_auto.inp"

# Geometry Parameters (mesh units: every length below is multiplied by length_scale)
cube_size = 2.0  # Size of the cube mesh (from 0 to 2.0)
top_face_y = 2.0 # The face where the indenter makes contact
size_x = None    # box size in x (None -> cube_size)
size_z = None    # box size in z (None -> cube_size)
length_scale = 1.0  # multiplies the node coordinates, indenter radius, gap and a numeric depth
                    # (e.g. a mesh built in units of R, with length_scale = R in um)

# Indenter Location and Geometry
indenter_position = "vertex"  # Options: "vertex" or "center"
indenter_radius = 0.5
initial_gap = 0.25           # Gap between indenter tip and mesh surface before step 1

# UEL Parameters
uel_variables = 848
dummy_variables = 106
num_properties = 26

# ID Management (Increase these if you use a very dense mesh)
dummy_element_offset = 100000
ref_node_id = 100000

# ==============================================================================
# CHEMICAL BOUNDARY CONFIGURATION
# ==============================================================================
# "closed" -> Impermeable boundary. Fluid/proteins cannot cross this surface.
# "open"   -> Permeable boundary. Fluid/proteins can escape into the surrounding bath.
top_boundary_condition = "closed"     # Top surface (outside the indenter)
side_boundary_condition = "closed"    # Outer physical side walls (x = 0, z = 0; also x = max, z = max for "center")
bottom_boundary_condition = "closed"  # Bottom surface (y = 0)
# Single side faces (opened in addition to the above)
x0_boundary_condition = "closed"      # x = 0
z0_boundary_condition = "closed"      # z = 0
xmax_boundary_condition = "closed"    # x = size_x ("center" only: symmetry plane for "vertex")
zmax_boundary_condition = "closed"    # z = size_z ("center" only: symmetry plane for "vertex")

bath_chemical_potential = "<INITMU>" # The potential of the surrounding bath

# ==============================================================================
# STEP CONFIGURATIONS (Set time > 0.0 to enable step, 0.0 to disable)
# ==============================================================================
t_equil = 0.0     # Equilibrate step before the indentation, indenter fixed (e.g. chemical
dtime_equil = 1e-2  # equilibrium under a filament prestress, LAMBDA0 > 1)
max_inc_equil = 10.0

t_indent = 1.0
depth = -0.5      # Indenter displacement (< 0: into the gel, gap included): a number (mesh units,
                  # scaled by length_scale) or an Abaqus parameter, e.g. "<UIND>" (um, set per case)
dtime_indent = 1e-3
max_inc_indent = 0.02
ramp_indenter = False  # True: explicit linear ramp of the indenter displacement over the
                       # Indentation (0 -> depth) and Withdraw (depth -> 0) steps

t_hold = 10.0
dtime_hold = 1e-4
max_inc_hold = 0.5

t_withdraw = 1.0
dtime_withdraw = 1e-4
max_inc_withdraw = 0.02

t_relax = 10.0
dtime_relax = 1e-4
max_inc_relax = 0.5

# ==============================================================================
# OUTPUT CONFIGURATION
# ==============================================================================
node_output_vars = "U, NT, RF"     # add RFL: reaction flux of the chemical DOF (exchange with the bath)
element_output_vars = "UVARM, LE"  # dummy mesh (all elements)
probe_output_vars = None           # e.g. "UVARM": extra output on the elements along the indenter
                                   # axis (elset probe_mesh), e.g. all per-direction cb
field_number_interval = None       # None: field output at every increment; or {step name: number
                                   # of intervals}, e.g. {"Hold": 60} (steps not listed: every increment)

# ==============================================================================
# PRESETS (after the defaults above, so that they override them)
# ==============================================================================
# ------------------------------------------------------------------------------
# STAGE 3 (AFM indentation, quarter model; run_stage3.py). Mesh built in units of the indenter
# radius R (1_build_mesh.py stage-3 preset: 3 x 2.5 x 3, graded towards the indented corner,
# job_name = 'base_mesh_afm'); length_scale = R in um. Generate one input per scale with:
#   python3 2_uel_afm.py length_scale=20 output_file=stage3_afm_large_closed.inp
#   python3 2_uel_afm.py length_scale=0.666666666667 output_file=stage3_afm_small_closed.inp
# Open faces, e.g. a bath on the far side walls (x = 0, z = 0):
#   python3 2_uel_afm.py length_scale=20 side_boundary_condition=open output_file=stage3_afm_large_open.inp
input_file = "base_mesh_afm.inp"
output_file = "stage3_afm_large_closed.inp"
cube_size = 3.0;  top_face_y = 2.5;  size_x = 3.0;  size_z = 3.0
length_scale = 20.0                 # R = 20 um (60 x 50 x 60 um gel); 2/3 -> 2 x 1.67 x 2 um sample
indenter_position = "vertex";  indenter_radius = 1.0;  initial_gap = 0.005
depth = "<UIND>"                    # -(gap + delta) in um, set per case by run_stage3.py
ramp_indenter = True
t_equil = 200.0;  dtime_equil = 0.1;    max_inc_equil = 20.0
t_indent = 2.0;   dtime_indent = 0.002; max_inc_indent = 0.05
t_hold = 600.0;   dtime_hold = 0.01;    max_inc_hold = 10.0
t_withdraw = 0.0; t_relax = 0.0     # zero-length steps are skipped
node_output_vars = "U, NT, RF, RFL"
element_output_vars = "UVARM1, UVARM2, UVARM3, UVARM4, UVARM5, UVARM6, UVARM7, UVARM8, UVARM9, UVARM10, LE"
probe_output_vars = "UVARM"
field_number_interval = {"Equilibrate": 4, "Indentation": 20, "Hold": 120}
# ------------------------------------------------------------------------------

# Command-line overrides (after the presets): python3 2_uel_afm.py name=value ...
def _parse_value(text):
    try:
        return ast.literal_eval(text)
    except (ValueError, SyntaxError):
        pass
    try:
        return float(eval(text, {'__builtins__': {}}, {}))   # e.g. length_scale=2/3
    except Exception:
        return text                                           # plain string, e.g. ...=open

for arg in sys.argv[1:]:
    key, _, value = arg.partition('=')
    if key.startswith('_') or key not in globals() or not value:
        sys.exit(f"Unknown argument {arg!r}: use name=value with a parameter of this script")
    globals()[key] = _parse_value(value)
    print(f"Override: {key} = {globals()[key]!r}")

sx = size_x if size_x is not None else cube_size   # box sizes (mesh units)
sz = size_z if size_z is not None else cube_size

if indenter_position == "vertex":
    indenter_x = sx
    indenter_z = sz
elif indenter_position == "center":
    indenter_x = sx / 2.0
    indenter_z = sz / 2.0
else:
    raise ValueError("indenter_position must be 'vertex' or 'center'")

if indenter_position == "vertex" and "open" in (xmax_boundary_condition, zmax_boundary_condition):
    raise ValueError("x = max and z = max are symmetry planes for vertex indentation: keep them closed")

# Scaled geometry (um) written to the input file
ls = float(length_scale)
indenter_x_s, indenter_z_s = indenter_x * ls, indenter_z * ls
indenter_radius_s = indenter_radius * ls
initial_gap_s = initial_gap * ls
ref_point_y = (top_face_y + initial_gap + indenter_radius) * ls
depth_value = depth if isinstance(depth, str) else depth * ls

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
node_xyz = {}

skip_block = False
for line in lines:
    upper_line = line.upper()
    
    # Determine if we should skip this entire keyword block. The CAE node/element sets are
    # skipped too: their names (e.g. bottom_nodes = z = 0 in 1_build_mesh.py) are redefined below
    # from the coordinates, and Abaqus adds the nodes of a repeated *Nset to the same set
    if upper_line.startswith("*"):
        if not upper_line.startswith("**"):
            in_nodes = in_elements = False  # a keyword ends the node/element data
        if upper_line.startswith("*SURFACE") or upper_line.startswith("*CONTACT") or "_SURFACE" in upper_line \
                or upper_line.startswith("*NSET") or upper_line.startswith("*ELSET"):
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
            if abs(x - sx) < 1e-6:
                xmax_nodes.append(n_id)
            if abs(z - sz) < 1e-6:
                zmax_nodes.append(n_id)
            if abs(x - 0.0) < 1e-6:
                x0_nodes.append(n_id)
            if abs(z - 0.0) < 1e-6:
                z0_nodes.append(n_id)
            node_xyz[n_id] = (x, y, z)
            # Node sets are found in mesh units (above); the coordinates are written scaled
            if ls != 1.0:
                line = f"{n_id}, {x * ls:.10g}, {y * ls:.10g}, {z * ls:.10g}\n"

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
** INDENTER: radius = {indenter_radius_s:.10g}, initial_gap = {initial_gap_s:.10g}, length_scale = {ls:.10g}, position = {indenter_position}
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
uel_block += write_nset("top_nodes", list(top_nodes.keys()))
uel_block += write_nset("bottom_nodes", bottom_nodes)
uel_block += write_nset("xmax_nodes", xmax_nodes)
uel_block += write_nset("zmax_nodes", zmax_nodes)
uel_block += write_nset("x0_nodes", x0_nodes)
uel_block += write_nset("z0_nodes", z0_nodes)

# Chemical bath: all nodes of the open faces
open_chem_nodes = set()
if top_boundary_condition == "open":
    open_chem_nodes.update(free_flux_nodes)
if side_boundary_condition == "open":
    open_chem_nodes.update(x0_nodes + z0_nodes)
    if indenter_position == "center":
        open_chem_nodes.update(xmax_nodes + zmax_nodes)
if bottom_boundary_condition == "open":
    open_chem_nodes.update(bottom_nodes)
if x0_boundary_condition == "open":
    open_chem_nodes.update(x0_nodes)
if z0_boundary_condition == "open":
    open_chem_nodes.update(z0_nodes)
if xmax_boundary_condition == "open":
    open_chem_nodes.update(xmax_nodes)
if zmax_boundary_condition == "open":
    open_chem_nodes.update(zmax_nodes)
print(f"Chemical bath (mu = {bath_chemical_potential}) on {len(open_chem_nodes)} of {len(all_node_ids)} nodes")
uel_block += write_nset("open_chem_nodes", sorted(open_chem_nodes))

# Probe elements: the elements along the indenter axis (dummy mesh), for extra output
probe_elements = []
for el in element_lines:
    parts = el.split(',')
    if len(parts) == 9:
        nodes = [int(n) for n in parts[1:]]
        if any(abs(node_xyz[n][0] - indenter_x) < 1e-6 and abs(node_xyz[n][2] - indenter_z) < 1e-6
               for n in nodes):
            probe_elements.append(int(parts[0]) + dummy_element_offset)
if probe_output_vars and probe_elements:
    uel_block += f"*Elset, elset=probe_mesh\n"
    for i in range(0, len(probe_elements), 16):
        uel_block += ", ".join(map(str, probe_elements[i:i+16])) + "\n"
    print(f"Probe output ({probe_output_vars}) on {len(probe_elements)} elements along the indenter axis")

uel_block += f"""
** ==============================================================================
** 4. RIGID INDENTER AND CONTACT
** ==============================================================================
*Node, nset=IndenterRef
{ref_node_id}, {indenter_x_s:.10g}, {ref_point_y:.10g}, {indenter_z_s:.10g}
*Surface, type=REVOLUTION, name=SURFACE-Probe
********
{indenter_x_s:.10g}, {ref_point_y:.10g}, {indenter_z_s:.10g},   {indenter_x_s:.10g}, {ref_point_y + indenter_radius_s:.10g}, {indenter_z_s:.10g}
START,          {indenter_radius_s:.10g},           0.
 CIRCL,           0.,         -{indenter_radius_s:.10g},           0.,           0.
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

# Explicit ramps of the indenter displacement (model data, before the steps)
if ramp_indenter:
    uel_block += f"""*Amplitude, name=ramp_load, time=STEP TIME
0.0, 0.0, {t_indent}, 1.0
"""
    if t_withdraw > 0.0:
        uel_block += f"""*Amplitude, name=ramp_unload, time=STEP TIME
0.0, 1.0, {t_withdraw}, 0.0
"""

def output_block(step_name):
    """Field output (every increment, or field_number_interval[step_name] intervals), and the
    indenter history at every increment."""
    n_int = (field_number_interval or {}).get(step_name)
    freq = f"number interval={n_int}" if n_int else "frequency=1"
    block = f"""*Output, field, {freq}, time marks=no
*Node Output, nset=all_nodes
{node_output_vars}
*Node Output, nset=IndenterRef
U, RF
*Element Output, elset=dummy_mesh
{element_output_vars}
"""
    if probe_output_vars and probe_elements:
        block += f"""*Element Output, elset=probe_mesh
{probe_output_vars}
"""
    block += """*node output, nset=extra_element
u
*Output, history, frequency=1
*Node Output, nset=IndenterRef
U2, RF2
*End Step
"""
    return block

def fixed_bcs():
    """Constant boundary conditions, defined in the first step (they carry over)."""
    bcs = """*Boundary
** Fix the bottom of the gel
bottom_nodes, YSYMM
"""
    if indenter_position == "vertex":
        bcs += """** Apply Quarter Symmetry boundary conditions for vertex indentation
xmax_nodes, XSYMM
zmax_nodes, ZSYMM
"""
    bcs += """** Extra Element Fix
extra_element, encastre
extra_element, 11, 11, 0.0
** Indenter Constraints
IndenterRef, 1, 1
IndenterRef, 3, 3
IndenterRef, 4, 4
IndenterRef, 5, 5
IndenterRef, 6, 6
"""
    if open_chem_nodes:
        bcs += f"""** Chemical Boundary (bath on the open faces)
open_chem_nodes, 11, 11, {bath_chemical_potential}
"""
    return bcs

first_step = [True]

def write_step(name, dtime, t, max_inc, u_indenter, amplitude=None):
    """One step; the indenter displacement is in its own *Boundary block (ramped if amplitude)."""
    if t <= 0.0:
        return ""
    step = f"""
*Step, name={name}, nlgeom=YES, inc=10000
*Coupled Temperature-displacement, creep=none, deltmx=10.0
{dtime}, {t}, 1e-15, {max_inc}
"""
    if first_step[0]:
        step += fixed_bcs()
        first_step[0] = False
    amp = f", amplitude={amplitude}" if amplitude else ""
    step += f"""*Boundary{amp}
IndenterRef, 2, 2, {u_indenter}
"""
    return step + output_block(name)

uel_block += write_step("Equilibrate", dtime_equil, t_equil, max_inc_equil, 0.0)
uel_block += write_step("Indentation", dtime_indent, t_indent, max_inc_indent, depth_value,
                        "ramp_load" if ramp_indenter else None)
uel_block += write_step("Hold", dtime_hold, t_hold, max_inc_hold, depth_value)
if ramp_indenter:   # ramp_unload: depth -> 0
    uel_block += write_step("Withdraw", dtime_withdraw, t_withdraw, max_inc_withdraw, depth_value, "ramp_unload")
else:
    uel_block += write_step("Withdraw", dtime_withdraw, t_withdraw, max_inc_withdraw, 0.0)
uel_block += write_step("Relax", dtime_relax, t_relax, max_inc_relax, 0.0)

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

print(f"Done! You can now run the job with 'abaqus job={os.path.splitext(output_file)[0]}'")
