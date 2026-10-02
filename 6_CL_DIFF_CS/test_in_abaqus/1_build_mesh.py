from abaqus import *
from abaqusConstants import *
import regionToolset
import mesh

# ==============================================================================
# PARAMETERS
# ==============================================================================
# Geometry: box size_x x size_y x size_z (a cube when all equal); y is the "vertical" axis
# used by 3_uel_cube.py (bottom face y = 0, top face y = size_y)
size_x = 20.0
size_y = 1.0
size_z = 1.0

# Mesh: elements per direction (n_x, n_y, n_z); None -> from the global size mesh_size
mesh_size = 1.0
n_x = 20
n_y = 1
n_z = 1

# Grading: bias_axis is one axis ('x', 'y' or 'z') or several, e.g. ('x', 'y', 'z') to grade
# towards a corner; ratio largest/smallest element bias_ratio (1.0 = uniform); smallest elements
# at bias_end: 'max' (e.g. the bath face of the stage-2 column at x = size_x) or 'min'
# (coordinate 0). bias_ratio and bias_end apply to every graded axis
bias_axis = 'x'
bias_ratio = 1.0
bias_end = 'max'

job_name = 'base_mesh_column'

# ------------------------------------------------------------------------------
# Presets (replace the values above):
#   Stage 1 (single element):  size_x = size_y = size_z = 1.0; mesh_size = 1.0;
#                              job_name = 'base_mesh_1el'
#   Stage 2 (column along x, bath at x = size_x): size_x = 50.0; size_y = size_z = 1.0;
#                              n_x = 15; n_y = n_z = 1; bias_axis = 'x'; bias_ratio = 3.0;
#                              bias_end = 'max'; job_name = 'base_mesh_column'
#   Indentation cube (previous): size_x = size_y = size_z = 2.0; mesh_size = 0.2;
#                              job_name = 'base_mesh'
#   Stage 3 (AFM quarter model, in units of the indenter radius R; 2_uel_afm.py scales it):
#                              size_x = size_z = 3.0; size_y = 2.5; n_x = n_z = 16; n_y = 14;
#                              bias_axis = ('x', 'y', 'z'); bias_ratio = 10.0; bias_end = 'max';
#                              job_name = 'base_mesh_afm'
#                              (3584 elements, ~0.045 R at the indented corner, ~0.45 R far away;
#                              numElem in global.f90 must be >= the number of elements)
# ------------------------------------------------------------------------------

# ==============================================================================
# 1. SETUP MODEL
# ==============================================================================
mdb.models.changeKey(fromName='Model-1', toName='UEL_Model')
myModel = mdb.models['UEL_Model']

# CRITICAL for UELs: Export a flat input file (no *Part or *Instance blocks)
myModel.setValues(noPartsInputFile=ON)

# ==============================================================================
# 2. CREATE GEOMETRY
# ==============================================================================
myPart = myModel.Part(name='Cube', dimensionality=THREE_D, type=DEFORMABLE_BODY)
mySketch = myModel.ConstrainedSketch(name='sketch', sheetSize=2.0 * max(size_x, size_y, size_z, 5.0))
mySketch.rectangle(point1=(0.0, 0.0), point2=(size_x, size_y))
myPart.BaseSolidExtrude(sketch=mySketch, depth=size_z)

# ==============================================================================
# 3. MESH THE BOX
# ==============================================================================
def edges_at(points):
    """Edges of the part through the given points (edge midpoints)."""
    return myPart.edges.findAt(*[(p, ) for p in points])

def edge_start(edge):
    """Coordinates of the start (parameter 0) of an edge."""
    try:
        return edge.getCurvature(parameter=0.0)['evaluationPoint']
    except Exception:
        return myPart.vertices[edge.getVertices()[0]].pointOn[0]

SIZES = (size_x, size_y, size_z)
AXIS = {'x': 0, 'y': 1, 'z': 2}

def axis_edges(a):
    """The 4 edges parallel to axis a (0, 1, 2), found at their midpoints."""
    b, c = [k for k in range(3) if k != a]
    points = []
    for vb in (0.0, SIZES[b]):
        for vc in (0.0, SIZES[c]):
            p = [0.0, 0.0, 0.0]
            p[a], p[b], p[c] = 0.5 * SIZES[a], vb, vc
            points.append(tuple(p))
    return edges_at(points)

def axis_spacings(a):
    """Element spacing along axis a at coordinate 0 and at the max face, on its 4 edges."""
    b, c = [k for k in range(3) if k != a]
    out = []
    for vb in (0.0, SIZES[b]):
        for vc in (0.0, SIZES[c]):
            vals = sorted(n.coordinates[a] for n in myPart.nodes
                          if abs(n.coordinates[b] - vb) < 1e-6 and abs(n.coordinates[c] - vc) < 1e-6)
            out.append((vals[1] - vals[0], vals[-1] - vals[-2]))
    return out

bias_axes = [AXIS[b] for b in ((bias_axis, ) if isinstance(bias_axis, str) else bias_axis)]
uniform_global = n_x is None and n_y is None and n_z is None and bias_ratio == 1.0
if uniform_global:
    myPart.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
else:
    counts = [n or max(1, int(round(size / mesh_size))) for n, size in zip((n_x, n_y, n_z), SIZES)]
    for a in range(3):
        if a not in bias_axes or bias_ratio == 1.0:
            myPart.seedEdgeByNumber(edges=axis_edges(a), number=counts[a], constraint=FIXED)
    if bias_ratio != 1.0:
        for ab in bias_axes:
            # Smallest elements at the parameter start (end1) or end (end2) of each edge, chosen so
            # that they lie at the requested face (coordinate = size for 'max', 0 for 'min')
            fine = SIZES[ab] if bias_end == 'max' else 0.0
            b_edges = axis_edges(ab)
            end1 = [e for e in b_edges if abs(edge_start(e)[ab] - fine) < 1e-6]
            end2 = [e for e in b_edges if abs(edge_start(e)[ab] - fine) >= 1e-6]
            kwargs = dict(biasMethod=SINGLE, ratio=bias_ratio, number=counts[ab], constraint=FIXED)
            if end1:
                kwargs['end1Edges'] = end1
            if end2:
                kwargs['end2Edges'] = end2
            myPart.seedEdgeByBias(**kwargs)

elemType = mesh.ElemType(elemCode=C3D8, elemLibrary=STANDARD)
myPart.setElementType(regions=(myPart.cells,), elemTypes=(elemType,))
myPart.generateMesh()

# Check the grading: finest spacing at the requested face, on all 4 edges along each graded axis
if not uniform_global and bias_ratio != 1.0:
    for ab in bias_axes:
        name = 'xyz'[ab]
        spacings = axis_spacings(ab)
        want_max = bias_end == 'max'
        ok = all((d_max < d_min) if want_max else (d_min < d_max) for d_min, d_max in spacings)
        print("Spacing along %s at %s = 0 / %s = %g on its 4 edges: %s" % (
            name, name, name, SIZES[ab], ', '.join('%.4g / %.4g' % s for s in spacings)))
        if not ok:
            raise RuntimeError("Grading along %s is not finest at bias_end = '%s' on all edges: "
                               "check the edge directions (end1Edges/end2Edges)" % (name, bias_end))

print("Mesh: %d nodes, %d elements" % (len(myPart.nodes), len(myPart.elements)))

# ==============================================================================
# 4. CREATE NODE SETS AND SURFACES
# ==============================================================================
# We define sets on the ASSEMBLY so they appear globally in the flat .inp
# (3_uel_cube.py rebuilds its own node sets from the coordinates)
myAssembly = myModel.rootAssembly
myInstance = myAssembly.Instance(name='Cube-1', part=myPart, dependent=ON)

# All Nodes Set
all_nodes = myInstance.nodes
myAssembly.Set(nodes=all_nodes, name='all_nodes')

# Bottom Nodes Set (Z = 0)
bottom_nodes = myInstance.nodes.getByBoundingBox(
    xMin=-0.1, xMax=size_x+0.1,
    yMin=-0.1, yMax=size_y+0.1,
    zMin=-0.1, zMax=0.1
)
myAssembly.Set(nodes=bottom_nodes, name='bottom_nodes')

# Top Nodes Set (Z = size_z)
top_nodes = myInstance.nodes.getByBoundingBox(
    xMin=-0.1, xMax=size_x+0.1,
    yMin=-0.1, yMax=size_y+0.1,
    zMin=size_z-0.1, zMax=size_z+0.1
)
myAssembly.Set(nodes=top_nodes, name='top_nodes')

# Top Surface for Contact
top_faces = myInstance.faces.getByBoundingBox(
    xMin=-0.1, xMax=size_x+0.1,
    yMin=-0.1, yMax=size_y+0.1,
    zMin=size_z-0.1, zMax=size_z+0.1
)
myAssembly.Surface(side1Faces=top_faces, name='SURFACE-Top')

# ==============================================================================
# 5. EXPORT RAW INP FILE
# ==============================================================================
mdb.Job(name=job_name, model='UEL_Model', type=ANALYSIS)
mdb.jobs[job_name].writeInput(consistencyChecking=OFF)
print("Successfully generated %s.inp!" % job_name)
