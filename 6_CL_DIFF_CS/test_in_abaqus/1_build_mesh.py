from abaqus import *
from abaqusConstants import *
import regionToolset
import mesh

# ==============================================================================
# PARAMETERS
# ==============================================================================
cube_size = 2.0
mesh_size = 0.2
job_name = 'base_mesh'

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
mySketch = myModel.ConstrainedSketch(name='sketch', sheetSize=10.0)
mySketch.rectangle(point1=(0.0, 0.0), point2=(cube_size, cube_size))
myPart.BaseSolidExtrude(sketch=mySketch, depth=cube_size)

# ==============================================================================
# 3. MESH THE CUBE
# ==============================================================================
myPart.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
elemType = mesh.ElemType(elemCode=C3D8, elemLibrary=STANDARD)
myPart.setElementType(regions=(myPart.cells,), elemTypes=(elemType,))
myPart.generateMesh()

# ==============================================================================
# 4. CREATE NODE SETS AND SURFACES
# ==============================================================================
# We define sets on the ASSEMBLY so they appear globally in the flat .inp
myAssembly = myModel.rootAssembly
myInstance = myAssembly.Instance(name='Cube-1', part=myPart, dependent=ON)

# All Nodes Set
all_nodes = myInstance.nodes
myAssembly.Set(nodes=all_nodes, name='all_nodes')

# Bottom Nodes Set (Z = 0)
bottom_nodes = myInstance.nodes.getByBoundingBox(
    xMin=-0.1, xMax=cube_size+0.1, 
    yMin=-0.1, yMax=cube_size+0.1, 
    zMin=-0.1, zMax=0.1
)
myAssembly.Set(nodes=bottom_nodes, name='bottom_nodes')

# Top Nodes Set (Z = cube_size)
top_nodes = myInstance.nodes.getByBoundingBox(
    xMin=-0.1, xMax=cube_size+0.1, 
    yMin=-0.1, yMax=cube_size+0.1, 
    zMin=cube_size-0.1, zMax=cube_size+0.1
)
myAssembly.Set(nodes=top_nodes, name='top_nodes')

# Top Surface for Contact
top_faces = myInstance.faces.getByBoundingBox(
    xMin=-0.1, xMax=cube_size+0.1, 
    yMin=-0.1, yMax=cube_size+0.1, 
    zMin=cube_size-0.1, zMax=cube_size+0.1
)
myAssembly.Surface(side1Faces=top_faces, name='SURFACE-Top')

# ==============================================================================
# 5. EXPORT RAW INP FILE
# ==============================================================================
mdb.Job(name=job_name, model='UEL_Model', type=ANALYSIS)
mdb.jobs[job_name].writeInput(consistencyChecking=OFF)
print("Successfully generated base_mesh.inp!")
