# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 2025 replay file
# Internal Version: 2024_09_20-14.00.46 RELr427 198590
# Run by lpacheco on Mon May 18 16:32:55 2026
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=278.469116210938, 
    height=198.031066894531)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
o1 = session.openOdb(
    name='/home/lpacheco/Desktop/UEL-dynamicCL/4_AFFCL_COUPLED/test_in_abaqus/cube_indent_uel.odb', 
    readOnly=False)
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: /home/lpacheco/Desktop/UEL-dynamicCL/4_AFFCL_COUPLED/test_in_abaqus/cube_indent_uel.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     2
#: Number of Meshes:             2
#: Number of Element Sets:       8
#: Number of Node Sets:          8
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='MUGRAD', outputPosition=INTEGRATION_POINT, refinement=(
    INVARIANT, 'Magnitude'), )
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    SYMBOLS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(
    scratchCoordSystemDisplay=OFF, pointElements=OFF)
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    renderStyle=WIREFRAME)
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(triad=OFF, 
    title=OFF, state=OFF, annotations=OFF, compass=OFF)
session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(
    spectrum='Blue to red', maxValue=4.36038E-07, minValue=9.48758E-09)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    SYMBOLS_ON_DEF, ))
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    renderStyle=WIREFRAME, )
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    renderStyle=WIREFRAME, )
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    SYMBOLS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.symbolOptions.setValues(
    vectorColorSpectrum='Blue to red', vectorMaxValue=4.65989E-07, 
    vectorMinValue=1.16574E-08)
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.fitView()
session.printOptions.setValues(vpDecorations=OFF)
session.printToFile(fileName='mugrad_maxind', format=EPS, canvasObjects=(
    session.viewports['Viewport: 1'], ))
session.printToFile(fileName='mugrad_maxind', format=PNG, canvasObjects=(
    session.viewports['Viewport: 1'], ))
session.printToFile(fileName='mugrad_maxind', format=SVG, canvasObjects=(
    session.viewports['Viewport: 1'], ))
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].odbDisplay.symbolOptions.setValues(
    vectorMaxValueAutoCompute=OFF, vectorMinValueAutoCompute=OFF, 
    vectorMinValue=0)
session.printToFile(fileName='mugrad_maxind', format=SVG, canvasObjects=(
    session.viewports['Viewport: 1'], ))
session.printToFile(fileName='mugrad_maxind', format=PNG, canvasObjects=(
    session.viewports['Viewport: 1'], ))
session.printToFile(fileName='mugrad_maxind', format=EPS, canvasObjects=(
    session.viewports['Viewport: 1'], ))
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='indentation', frame=18)
session.viewports['Viewport: 1'].view.fitView()
session.printToFile(fileName='mugrad_half_ind', format=EPS, canvasObjects=(
    session.viewports['Viewport: 1'], ))
session.printToFile(fileName='mugrad_half_ind', format=SVG, canvasObjects=(
    session.viewports['Viewport: 1'], ))
session.printToFile(fileName='mugrad_half_ind', format=PNG, canvasObjects=(
    session.viewports['Viewport: 1'], ))
