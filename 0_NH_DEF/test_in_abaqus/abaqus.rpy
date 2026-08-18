# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 2024 replay file
# Internal Version: 2023_09_21-13.55.25 RELr426 190762
# Run by lpacheco on Mon Aug 17 22:14:44 2026
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=435.818725585938, 
    height=569.240173339844)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
o1 = session.openOdb(
    name='/home/lpacheco/UEL-dynamicCL/0_NH_DEF/test_in_abaqus/uel_cube.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: /home/lpacheco/UEL-dynamicCL/0_NH_DEF/test_in_abaqus/uel_cube.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       4
#: Number of Node Sets:          12
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
#: Warning: Results for the current deformed variable are not available for one or more nodes contained in the model. Deformations at such nodes are assumed to be zero.
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='RF', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
#: Warning: Results for the current deformed variable are not available for one or more nodes contained in the model. Deformations at such nodes are assumed to be zero.
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='static1', frame=9)
#: Warning: Results for the current deformed variable are not available for one or more nodes contained in the model. Deformations at such nodes are assumed to be zero.
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='static1', frame=8)
