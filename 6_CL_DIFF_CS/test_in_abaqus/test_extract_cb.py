# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 2024 replay file
# Internal Version: 2023_09_21-13.55.25 RELr426 190762
# Run by lpacheco on Tue Sep  1 01:16:45 2026
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=408.580047607422, 
    height=311.98486328125)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
import numpy as np
executeOnCaeStartup()
o1 = session.openOdb(
    name='/home/lpacheco/UEL-dynamicCL/5_CL_DIFF/test_in_abaqus/cube_indent_uel_relax_eta067.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: /home/lpacheco/UEL-dynamicCL/5_CL_DIFF/test_in_abaqus/cube_indent_uel_relax_eta067.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     2
#: Number of Meshes:             2
#: Number of Element Sets:       8
#: Number of Node Sets:          8
#: Number of Steps:              2
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.fitView()
odb = session.odbs['/home/lpacheco/UEL-dynamicCL/5_CL_DIFF/test_in_abaqus/cube_indent_uel_relax_eta067.odb']
session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=((
    'UVARM100', INTEGRATION_POINT), ), nodePick=(('PART-1-1', 1, ('[#1 ]', )), 
    ), )
# xyp = session.XYPlot('XYPlot-1')
xy1 = session.xyDataObjects['UVARM100 (Avg: 100%) PI: PART-1-1 N: 1']
c1 = session.Curve(xyData=xy1)

cb_array = []
time_array = [point[0] for point in xy1]
for i in range(17, 107):
    uvarm_name = 'UVARM'+str(i)
    session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=((
    uvarm_name, INTEGRATION_POINT), ), nodePick=(('PART-1-1', 1, ('[#1 ]', )), 
    ), )
    aux_str = uvarm_name + ' (Avg: 100%) PI: PART-1-1 N: 1'
    xy1 = session.xyDataObjects[aux_str]
    cb_aux_array = [point[1] for point in xy1]
    cb_array.append(cb_aux_array)

cb_array = np.array(cb_array)

cb_array1s = cb_array[:, 29]
cb_array10s = cb_array[:, 47]
np.save('cb_array1s.npy', cb_array1s)
np.save('cb_array10s.npy', cb_array10s)
np.save('time_array.npy', time_array)
print(cb_array.shape)
print(time_array)


