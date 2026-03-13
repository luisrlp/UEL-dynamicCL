# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2024 replay file
# Internal Version: 2023_09_21-13.55.25 RELr426 190762
# Run by lpacheco on Fri Mar 13 11:13:00 2026
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=423.389434814453, 
    height=283.897888183594)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=4.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=(0.0, 0.0), point2=(2.0, 2.0))
p = mdb.models['Model-1'].Part(name='cube', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['cube']
p.BaseSolidExtrude(sketch=s, depth=2.0)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['cube']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=4.0)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=STANDALONE)
s1.ConstructionLine(point1=(0.0, -2.0), point2=(0.0, 2.0))
s1.FixedConstraint(entity=g[2])
s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(1.0, 0.0))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.60926, 
    farPlane=3.93321, width=2.33251, height=1.26549, cameraPosition=(0.0132254, 
    -0.00141003, 3.77124), cameraTarget=(0.0132254, -0.00141003, 0))
s1.Line(point1=(0.0, 0.0), point2=(1.0, 0.0))
s1.HorizontalConstraint(entity=g[4], addUndoState=False)
s1.Line(point1=(0.0, 0.0), point2=(0.0, -1.0))
s1.VerticalConstraint(entity=g[5], addUndoState=False)
s1.autoTrimCurve(curve1=g[3], point1=(-0.961063504219055, 0.270868510007858))
s1.autoTrimCurve(curve1=g[5], point1=(-0.00166182965040207, 
    -0.139202788472176))
s1.autoTrimCurve(curve1=g[4], point1=(0.203451573848724, -0.00361468270421028))
s1.undo()
s1.autoTrimCurve(curve1=g[6], point1=(0.988617122173309, 0.15842966735363))
s1.autoTrimCurve(curve1=g[4], point1=(0.482450067996979, -0.00251235440373421))
session.viewports['Viewport: 1'].view.fitView()
p = mdb.models['Model-1'].Part(name='probe', dimensionality=THREE_D, 
    type=ANALYTIC_RIGID_SURFACE)
p = mdb.models['Model-1'].parts['probe']
p.AnalyticRigidSurfRevolve(sketch=s1)
s1.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['probe']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
mdb.models['Model-1'].Material(name='Material-1')
mdb.models['Model-1'].materials['Material-1'].Elastic(table=((1.0, 0.45), ))
mdb.models['Model-1'].HomogeneousSolidSection(name='SECTION_elastic', 
    material='Material-1', thickness=None)
p = mdb.models['Model-1'].parts['cube']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['cube']
c = p.cells
cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
region = p.Set(cells=cells, name='Set-Cube')
p = mdb.models['Model-1'].parts['cube']
p.SectionAssignment(region=region, sectionName='SECTION_elastic', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['cube']
a.Instance(name='cube-1', part=p, dependent=OFF)
p = mdb.models['Model-1'].parts['probe']
a.Instance(name='probe-1', part=p, dependent=OFF)
p1 = a.instances['probe-1']
p1.translate(vector=(3.2, 0.0, 0.0))
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF, 
    engineeringFeatures=OFF)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p = mdb.models['Model-1'].parts['cube']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['probe']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.07876, 
    farPlane=8.38344, width=3.99877, height=2.11659, cameraPosition=(2.75561, 
    -5.51049, 3.55102), cameraUpVector=(-0.177375, 0.779525, 0.600732), 
    cameraTarget=(0.0539496, -0.607898, 0.0539494))
p = mdb.models['Model-1'].parts['probe']
v1, e, d1, n = p.vertices, p.edges, p.datums, p.nodes
p.ReferencePoint(point=v1[0])
a = mdb.models['Model-1'].rootAssembly
a.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
a = mdb.models['Model-1'].rootAssembly
a.translate(instanceList=('probe-1', ), vector=(-1.2, 3.0, 2.0))
#: The instance probe-1 was translated by -1.2, 3., 2. with respect to the assembly coordinate system
a = mdb.models['Model-1'].rootAssembly
a.translate(instanceList=('probe-1', ), vector=(0.0, 1.0, 0.0))
#: The instance probe-1 was translated by 0., 1., 0. with respect to the assembly coordinate system
session.viewports['Viewport: 1'].view.fitView()
p = mdb.models['Model-1'].parts['probe']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['probe']
s = p.features['3D Analytic rigid shell-1'].sketch
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=s)
s2 = mdb.models['Model-1'].sketches['__edit__']
g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
s2.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s2, 
    upToFeature=p.features['3D Analytic rigid shell-1'], filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1.71409, 
    farPlane=3.94277, width=1.64651, height=0.893301, cameraPosition=(0.476719, 
    -0.546417, 2.82843), cameraTarget=(0.476719, -0.546417, 0))
s2.offset(distance=0.75, objectList=(g[7], ), side=LEFT)
s2.autoTrimCurve(curve1=g[7], point1=(0.537047803401947, -0.845222055912018))
s2.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['probe']
p.features['3D Analytic rigid shell-1'].setValues(sketch=s2)
del mdb.models['Model-1'].sketches['__edit__']
p = mdb.models['Model-1'].parts['probe']
p.regenerate()
#* FeatureError: Regeneration failed
p = mdb.models['Model-1'].parts['probe']
p.backup()
p = mdb.models['Model-1'].parts['probe']
del p.features['RP']
session.viewports['Viewport: 1'].view.fitView()
p = mdb.models['Model-1'].parts['probe']
v2, e1, d2, n1 = p.vertices, p.edges, p.datums, p.nodes
p.ReferencePoint(point=v2[0])
a = mdb.models['Model-1'].rootAssembly
a.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
a = mdb.models['Model-1'].rootAssembly
a.translate(instanceList=('probe-1', ), vector=(0.0, -1.75, 0.0))
#: The instance probe-1 was translated by 0., -1.75, 0. with respect to the assembly coordinate system
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
a = mdb.models['Model-1'].rootAssembly
a.translate(instanceList=('probe-1', ), vector=(0.0, 0.25, 0.0))
#: The instance probe-1 was translated by 0., 250.E-03, 0. with respect to the assembly coordinate system
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)
mdb.models['Model-1'].StaticStep(name='indentation', previous='Initial', 
    nlgeom=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='indentation')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
    constraints=ON, connectors=ON, engineeringFeatures=ON, 
    adaptiveMeshConstraints=OFF)
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.80385, 
    farPlane=10.7141, width=4.30473, height=2.27854, cameraPosition=(3.35694, 
    6.95394, 6.66593), cameraUpVector=(-0.215648, 0.714875, -0.665169))
mdb.models['Model-1'].ContactProperty('IntProp-1')
mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
    formulation=FRICTIONLESS)
#: The interaction property "IntProp-1" has been created.
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['probe-1'].faces
side1Faces1 = s1.getSequenceFromMask(mask=('[#1 ]', ), )
region1=a.Surface(side1Faces=side1Faces1, name='Surface-Probe')
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['cube-1'].faces
side1Faces1 = s1.getSequenceFromMask(mask=('[#2 ]', ), )
region2=a.Surface(side1Faces=side1Faces1, name='SURFACE-Top')
mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='INTERACTION-Contact', 
    createStepName='indentation', main=region1, secondary=region2, 
    sliding=FINITE, thickness=ON, interactionProperty='IntProp-1', 
    adjustMethod=NONE, initialClearance=OMIT, datumAxis=None, 
    clearanceRegion=None)
#: The interaction "INTERACTION-Contact" has been created.
mdb.models['Model-1'].interactionProperties.changeKey(fromName='IntProp-1', 
    toName='INTPROP-Frictionless')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, interactions=OFF, constraints=OFF, 
    engineeringFeatures=OFF)
a = mdb.models['Model-1'].rootAssembly
r1 = a.instances['probe-1'].referencePoints
refPoints1=(r1[3], )
region = a.Set(referencePoints=refPoints1, name='Set-RP')
mdb.models['Model-1'].DisplacementBC(name='BC-ProbeDisplacement', 
    createStepName='indentation', region=region, u1=0.0, u2=-0.5, u3=0.0, 
    ur1=0.0, ur2=0.0, ur3=0.0, amplitude=UNSET, fixed=OFF, 
    distributionType=UNIFORM, fieldName='', localCsys=None)
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.04732, 
    farPlane=10.6467, width=4.48531, height=2.37413, cameraPosition=(2.20875, 
    -5.48423, 5.93631), cameraUpVector=(-0.325662, 0.513825, 0.79368), 
    cameraTarget=(1.13714, 1.25647, 1.13271))
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['cube-1'].faces
faces1 = f1.getSequenceFromMask(mask=('[#8 ]', ), )
region = a.Set(faces=faces1, name='Set-Bottom')
mdb.models['Model-1'].EncastreBC(name='BC-Fix', createStepName='indentation', 
    region=region, localCsys=None)
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.82537, 
    farPlane=11.6773, width=5.94817, height=3.14843, cameraPosition=(2.29502, 
    -6.30558, 6.33577), cameraUpVector=(0.700347, 0.484735, 0.523972), 
    cameraTarget=(1.1072, 1.16611, 1.01125))
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.68931, 
    farPlane=11.9141, width=6.01726, height=3.185, cameraPosition=(6.70335, 
    -4.89163, 5.33282), cameraUpVector=(0.46488, 0.755185, 0.46215), 
    cameraTarget=(1.15864, 1.16998, 1.0052))
session.viewports['Viewport: 1'].view.fitView()
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['cube-1'].faces
faces1 = f1.getSequenceFromMask(mask=('[#4 ]', ), )
region = a.Set(faces=faces1, name='Set-XSYMM')
mdb.models['Model-1'].XsymmBC(name='BC-XSYMM', createStepName='indentation', 
    region=region, localCsys=None)
#: Warning: Cannot continue yet--complete the step or cancel the procedure.
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['cube-1'].faces
faces1 = f1.getSequenceFromMask(mask=('[#10 ]', ), )
region = a.Set(faces=faces1, name='Set-ZSYMM')
mdb.models['Model-1'].ZsymmBC(name='BC-ZSYMM', createStepName='indentation', 
    region=region, localCsys=None)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON, loads=OFF, 
    bcs=OFF, predefinedFields=OFF, connectors=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
p = mdb.models['Model-1'].parts['probe']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
p = mdb.models['Model-1'].parts['cube']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
a1 = mdb.models['Model-1'].rootAssembly
a1.makeDependent(instances=(a1.instances['cube-1'], ))
p = mdb.models['Model-1'].parts['cube']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['cube']
p.seedPart(size=1.0, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['cube']
p.generateMesh()
p = mdb.models['Model-1'].parts['cube']
p.deleteMesh()
p = mdb.models['Model-1'].parts['cube']
p.seedPart(size=0.5, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['cube']
p.generateMesh()
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
mdb.Job(name='cube_indent', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
    multiprocessingMode=DEFAULT, numCpus=2, numDomains=2, numGPUs=0)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
    constraints=ON, connectors=ON, engineeringFeatures=ON)
mdb.models['Model-1'].interactions['INTERACTION-Contact'].setValues(
    initialClearance=OMIT, adjustMethod=NONE, sliding=FINITE, 
    enforcement=SURFACE_TO_SURFACE, thickness=ON, contactTracking=TWO_CONFIG, 
    interactionProperty='INTPROP-Frictionless', bondingSet=None)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=OFF, 
    constraints=OFF, connectors=OFF, engineeringFeatures=OFF)
mdb.jobs['cube_indent'].submit(consistencyChecking=OFF)
#: The job input file "cube_indent.inp" has been submitted for analysis.
#: Job cube_indent: Analysis Input File Processor completed successfully.
#: Job cube_indent: Abaqus/Standard completed successfully.
#: Job cube_indent completed successfully. 
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)
mdb.models['Model-1'].steps['indentation'].setValues(initialInc=0.1, 
    maxInc=0.2)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=OFF)
mdb.jobs['cube_indent'].submit(consistencyChecking=OFF)
#: The job input file "cube_indent.inp" has been submitted for analysis.
#: Job cube_indent: Analysis Input File Processor completed successfully.
#: Job cube_indent: Abaqus/Standard completed successfully.
#: Job cube_indent completed successfully. 
o3 = session.openOdb(
    name='/home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb')
#: Model: /home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     2
#: Number of Meshes:             2
#: Number of Element Sets:       5
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
p = mdb.models['Model-1'].parts['cube']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
elemType1 = mesh.ElemType(elemCode=C3D8, elemLibrary=STANDARD, 
    secondOrderAccuracy=OFF, distortionControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
p = mdb.models['Model-1'].parts['cube']
c = p.cells
cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
pickedRegions =(cells, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
    elemType3))
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
mdb.jobs['cube_indent'].submit(consistencyChecking=OFF)
#: The job input file "cube_indent.inp" has been submitted for analysis.
#: Job cube_indent: Analysis Input File Processor completed successfully.
#: Job cube_indent: Abaqus/Standard completed successfully.
#: Job cube_indent completed successfully. 
o3 = session.openOdb(
    name='/home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb')
#: Model: /home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     2
#: Number of Meshes:             2
#: Number of Element Sets:       5
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON, mesh=OFF)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=OFF)
p = mdb.models['Model-1'].parts['cube']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
mdb.models['Model-1'].Material(name='Material-Neo-Hooke')
mdb.models['Model-1'].materials['Material-Neo-Hooke'].Hyperelastic(
    materialType=ISOTROPIC, testData=OFF, type=NEO_HOOKE, 
    volumetricResponse=VOLUMETRIC_DATA, table=((1.0, 0.002), ))
mdb.models['Model-1'].materials['Material-Neo-Hooke'].Density(table=((1.0, ), 
    ))
mdb.models['Model-1'].sections['SECTION_elastic'].setValues(
    material='Material-Neo-Hooke', thickness=None)
mdb.models['Model-1'].sections.changeKey(fromName='SECTION_elastic', 
    toName='SECTION_hyperelastic')
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
mdb.jobs['cube_indent'].submit(consistencyChecking=OFF)
#: The job input file "cube_indent.inp" has been submitted for analysis.
#: Error in job cube_indent: 64 elements have missing property definitions. The elements have been identified in element set ErrElemMissingSection.
p = mdb.models['Model-1'].parts['cube']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
#: Job cube_indent: Analysis Input File Processor aborted due to errors.
#: Error in job cube_indent: Analysis Input File Processor exited with an error - Please see the  cube_indent.dat file for possible error messages if the file exists.
#: Job cube_indent aborted due to errors.
mdb.models['Model-1'].parts['cube'].sectionAssignments[0].setValues(
    sectionName='SECTION_hyperelastic')
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
mdb.jobs['cube_indent'].submit(consistencyChecking=OFF)
#: The job input file "cube_indent.inp" has been submitted for analysis.
#: Job cube_indent: Analysis Input File Processor aborted due to errors.
#: Error in job cube_indent: Analysis Input File Processor exited with an error - Please see the  cube_indent.dat file for possible error messages if the file exists.
#: Job cube_indent aborted due to errors.
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
p = mdb.models['Model-1'].parts['cube']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF, 
    engineeringFeatures=OFF, mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
elemType1 = mesh.ElemType(elemCode=C3D8H, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
p = mdb.models['Model-1'].parts['cube']
c = p.cells
cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
pickedRegions =(cells, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
    elemType3))
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
mdb.jobs['cube_indent'].submit(consistencyChecking=OFF)
#: The job input file "cube_indent.inp" has been submitted for analysis.
#: Job cube_indent: Analysis Input File Processor completed successfully.
#: Job cube_indent: Abaqus/Standard completed successfully.
#: Job cube_indent completed successfully. 
o3 = session.openOdb(
    name='/home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb')
#: Model: /home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     2
#: Number of Meshes:             2
#: Number of Element Sets:       5
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
    constraints=ON, connectors=ON, engineeringFeatures=ON)
mdb.models['Model-1'].interactions['INTERACTION-Contact'].setValues(
    initialClearance=OMIT, surfaceSmoothing=NONE, adjustMethod=NONE, 
    sliding=FINITE, enforcement=NODE_TO_SURFACE, thickness=OFF, 
    supplementaryContact=SELECTIVE, smooth=0.2, bondingSet=None)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=OFF, 
    constraints=OFF, connectors=OFF, engineeringFeatures=OFF)
mdb.jobs['cube_indent'].submit(consistencyChecking=OFF)
#: The job input file "cube_indent.inp" has been submitted for analysis.
#: Job cube_indent: Analysis Input File Processor completed successfully.
#: Job cube_indent: Abaqus/Standard completed successfully.
#: Job cube_indent completed successfully. 
o3 = session.openOdb(
    name='/home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb')
#: Model: /home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     2
#: Number of Meshes:             2
#: Number of Element Sets:       5
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.99472, 
    farPlane=10.026, width=1.95241, height=1.00575, viewOffsetX=0.00628274, 
    viewOffsetY=0.133262)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
    constraints=ON, connectors=ON, engineeringFeatures=ON)
mdb.models['Model-1'].interactions['INTERACTION-Contact'].setValues(
    initialClearance=OMIT, surfaceSmoothing=NONE, adjustMethod=NONE, 
    sliding=FINITE, enforcement=SURFACE_TO_SURFACE, thickness=ON, 
    contactTracking=TWO_CONFIG, smooth=0.2, bondingSet=None)
mdb.models['Model-1'].interactionProperties['INTPROP-Frictionless'].tangentialBehavior.setValues(
    formulation=FRICTIONLESS)
mdb.models['Model-1'].interactionProperties['INTPROP-Frictionless'].NormalBehavior(
    pressureOverclosure=HARD, allowSeparation=ON, 
    constraintEnforcementMethod=DEFAULT)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=OFF, 
    constraints=OFF, connectors=OFF, engineeringFeatures=OFF)
mdb.jobs['cube_indent'].submit(consistencyChecking=OFF)
#: The job input file "cube_indent.inp" has been submitted for analysis.
#: Job cube_indent: Analysis Input File Processor completed successfully.
#: Job cube_indent: Abaqus/Standard completed successfully.
#: Job cube_indent completed successfully. 
o3 = session.openOdb(
    name='/home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb')
#: Model: /home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     2
#: Number of Meshes:             2
#: Number of Element Sets:       5
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
    constraints=ON, connectors=ON, engineeringFeatures=ON)
mdb.models['Model-1'].interactionProperties['INTPROP-Frictionless'].tangentialBehavior.setValues(
    formulation=FRICTIONLESS)
mdb.models['Model-1'].interactionProperties['INTPROP-Frictionless'].normalBehavior.setValues(
    pressureOverclosure=HARD, allowSeparation=ON, contactStiffness=DEFAULT, 
    contactStiffnessScaleFactor=1.0, clearanceAtZeroContactPressure=0.0, 
    constraintEnforcementMethod=AUGMENTED_LAGRANGE)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=OFF, 
    constraints=OFF, connectors=OFF, engineeringFeatures=OFF)
mdb.jobs['cube_indent'].submit(consistencyChecking=OFF)
#: The job input file "cube_indent.inp" has been submitted for analysis.
#: Job cube_indent: Analysis Input File Processor completed successfully.
#: Job cube_indent: Abaqus/Standard completed successfully.
#: Job cube_indent completed successfully. 
o3 = session.openOdb(
    name='/home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb')
#: Model: /home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     2
#: Number of Meshes:             2
#: Number of Element Sets:       5
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
    constraints=ON, connectors=ON, engineeringFeatures=ON)
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['probe-1'].faces
side2Faces1 = s1.getSequenceFromMask(mask=('[#1 ]', ), )
region1=a.Surface(side2Faces=side2Faces1, name='m_Surf-3')
mdb.models['Model-1'].interactions['INTERACTION-Contact'].setValues(
    main=region1, initialClearance=OMIT, adjustMethod=NONE, sliding=FINITE, 
    enforcement=SURFACE_TO_SURFACE, thickness=ON, contactTracking=TWO_CONFIG, 
    bondingSet=None)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=OFF, 
    constraints=OFF, connectors=OFF, engineeringFeatures=OFF)
del mdb.models['Model-1'].rootAssembly.surfaces['Surface-Probe']
mdb.jobs['cube_indent'].submit(consistencyChecking=OFF)
#: The job input file "cube_indent.inp" has been submitted for analysis.
#: Job cube_indent: Analysis Input File Processor completed successfully.
#: Job cube_indent: Abaqus/Standard completed successfully.
#: Job cube_indent completed successfully. 
o3 = session.openOdb(
    name='/home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb')
#: Model: /home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     2
#: Number of Meshes:             2
#: Number of Element Sets:       5
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
    constraints=ON, connectors=ON, engineeringFeatures=ON)
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['probe-1'].faces
side1Faces1 = s1.getSequenceFromMask(mask=('[#1 ]', ), )
region1=a.Surface(side1Faces=side1Faces1, name='m_Surf-4')
mdb.models['Model-1'].interactions['INTERACTION-Contact'].setValues(
    main=region1, initialClearance=OMIT, adjustMethod=NONE, sliding=FINITE, 
    enforcement=SURFACE_TO_SURFACE, thickness=ON, contactTracking=TWO_CONFIG, 
    bondingSet=None)
del mdb.models['Model-1'].rootAssembly.surfaces['m_Surf-3']
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=OFF, 
    constraints=OFF, connectors=OFF, engineeringFeatures=OFF)
mdb.jobs['cube_indent'].submit(consistencyChecking=OFF)
#: The job input file "cube_indent.inp" has been submitted for analysis.
#: Job cube_indent: Analysis Input File Processor completed successfully.
#: Job cube_indent: Abaqus/Standard completed successfully.
#: Job cube_indent completed successfully. 
o3 = session.openOdb(
    name='/home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb')
#: Model: /home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     2
#: Number of Meshes:             2
#: Number of Element Sets:       5
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON, mesh=OFF)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=OFF)
p = mdb.models['Model-1'].parts['cube']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
mdb.models['Model-1'].materials['Material-Neo-Hooke'].hyperelastic.setValues(
    table=((1.0, 0.005), ))
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
mdb.jobs['cube_indent'].submit(consistencyChecking=OFF)
#: The job input file "cube_indent.inp" has been submitted for analysis.
#: Job cube_indent: Analysis Input File Processor completed successfully.
#: Job cube_indent: Abaqus/Standard completed successfully.
#: Job cube_indent completed successfully. 
o3 = session.openOdb(
    name='/home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb')
#: Model: /home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     2
#: Number of Meshes:             2
#: Number of Element Sets:       5
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF, 
    engineeringFeatures=OFF, mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
p = mdb.models['Model-1'].parts['cube']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['cube']
p.deleteMesh()
p = mdb.models['Model-1'].parts['cube']
p.seedPart(size=0.2, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['cube']
p.generateMesh()
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
mdb.jobs['cube_indent'].submit(consistencyChecking=OFF)
#: The job input file "cube_indent.inp" has been submitted for analysis.
#: Job cube_indent: Analysis Input File Processor completed successfully.
#: Job cube_indent: Abaqus/Standard completed successfully.
#: Job cube_indent completed successfully. 
o3 = session.openOdb(
    name='/home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb')
#: Model: /home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     2
#: Number of Meshes:             2
#: Number of Element Sets:       5
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='indentation', frame=1)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='indentation', frame=2)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='indentation', frame=3)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='indentation', frame=4)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='indentation', frame=3)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='indentation', frame=4)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='indentation', frame=5)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='indentation', frame=6)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='indentation', frame=7)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='indentation', frame=8)
session.viewports['Viewport: 1'].view.setValues(session.views['Left'])
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.setValues(session.views['Right'])
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=OFF)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p = mdb.models['Model-1'].parts['cube']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p1 = mdb.models['Model-1'].parts['probe']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['Model-1'].parts['probe']
s = p.features['3D Analytic rigid shell-1'].sketch
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=s)
s1 = mdb.models['Model-1'].sketches['__edit__']
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s1, 
    upToFeature=p.features['3D Analytic rigid shell-1'], filter=COPLANAR_EDGES)
s1.offset(distance=0.05, objectList=(g[8], ), side=RIGHT)
s1.autoTrimCurve(curve1=g[8], point1=(0.0949200093746185, -0.231783509254456))
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.429094, 
    farPlane=0.98512, width=0.490031, height=0.265863, cameraPosition=(
    0.124383, -0.124002, 0.707107), cameraTarget=(0.124383, -0.124002, 0))
s1.offset(distance=0.2, objectList=(g[9], ), side=RIGHT)
s1.autoTrimCurve(curve1=g[9], point1=(0.202805683016777, -0.221501022577286))
s1.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['probe']
p.features['3D Analytic rigid shell-1'].setValues(sketch=s1)
del mdb.models['Model-1'].sketches['__edit__']
p = mdb.models['Model-1'].parts['probe']
p.regenerate()
#* FeatureError: Regeneration failed
p = mdb.models['Model-1'].parts['probe']
p.backup()
p = mdb.models['Model-1'].parts['probe']
del p.features['RP']
p = mdb.models['Model-1'].parts['probe']
v1, e, d1, n = p.vertices, p.edges, p.datums, p.nodes
p.ReferencePoint(point=v1[0])
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
a1 = mdb.models['Model-1'].rootAssembly
a1.translate(instanceList=('probe-1', ), vector=(0.0, 0.25, 0.0))
#: The instance probe-1 was translated by 0., 250.E-03, 0. with respect to the assembly coordinate system
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
    constraints=ON, connectors=ON, engineeringFeatures=ON)
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['probe-1'].faces
side1Faces1 = s1.getSequenceFromMask(mask=('[#1 ]', ), )
region1=a.Surface(side1Faces=side1Faces1, name='SURFACE-Probe')
mdb.models['Model-1'].interactions['INTERACTION-Contact'].setValues(
    main=region1, initialClearance=OMIT, adjustMethod=NONE, sliding=FINITE, 
    enforcement=SURFACE_TO_SURFACE, thickness=ON, contactTracking=TWO_CONFIG, 
    bondingSet=None)
del mdb.models['Model-1'].rootAssembly.surfaces['m_Surf-4']
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, interactions=OFF, constraints=OFF, 
    engineeringFeatures=OFF)
a = mdb.models['Model-1'].rootAssembly
region = a.sets['Set-RP']
mdb.models['Model-1'].boundaryConditions['BC-ProbeDisplacement'].setValues(
    region=region)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON, loads=OFF, 
    bcs=OFF, predefinedFields=OFF, connectors=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
p = mdb.models['Model-1'].parts['probe']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
p = mdb.models['Model-1'].parts['cube']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['cube']
p.deleteMesh()
p = mdb.models['Model-1'].parts['cube']
p.seedPart(size=0.5, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['cube']
p.generateMesh()
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
del mdb.models['Model-1'].rootAssembly.sets['Set-RP']
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON)
a = mdb.models['Model-1'].rootAssembly
r1 = a.instances['probe-1'].referencePoints
refPoints1=(r1[4], )
region = a.Set(referencePoints=refPoints1, name='Set-RP')
mdb.models['Model-1'].boundaryConditions['BC-ProbeDisplacement'].setValues(
    region=region)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
mdb.jobs['cube_indent'].submit(consistencyChecking=OFF)
#: The job input file "cube_indent.inp" has been submitted for analysis.
#: Job cube_indent: Analysis Input File Processor completed successfully.
#: Job cube_indent: Abaqus/Standard completed successfully.
#: Job cube_indent completed successfully. 
o3 = session.openOdb(
    name='/home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb')
#: Model: /home/lpacheco/UEL-ABAQUS/4_AFFCL_COUPLED/test_in_abaqus/cube_indent.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     2
#: Number of Meshes:             2
#: Number of Element Sets:       5
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.setValues(session.views['Right'])
session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
