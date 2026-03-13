# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 2024 replay file
# Internal Version: 2023_09_21-13.55.25 RELr426 190762
# Run by lpacheco on Fri Feb 27 11:41:57 2026
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=296.451934814453, 
    height=198.172912597656)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
o1 = session.openOdb(
    name='/home/lpacheco/UEL-ABAQUS/3_NH_COUPLED/test_in_abaqus/uel_cube.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: /home/lpacheco/UEL-ABAQUS/3_NH_COUPLED/test_in_abaqus/uel_cube.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       5
#: Number of Node Sets:          12
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Swell', frame=56)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Swell', frame=17)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Swell', frame=18)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Swell', frame=4)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Swell', frame=18)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='UVARM3', outputPosition=INTEGRATION_POINT, )
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Swell', frame=37)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Swell', frame=44)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Swell', frame=45)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Swell', frame=46)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Swell', frame=47)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Swell', frame=48)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Swell', frame=49)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Swell', frame=50)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Swell', frame=49)
session.viewports['Viewport: 1'].view.fitView()
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Swell', frame=50)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Swell', frame=51)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Swell', frame=52)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Swell', frame=53)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Swell', frame=54)
