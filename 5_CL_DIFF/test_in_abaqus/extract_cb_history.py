# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 2024 replay file
# Internal Version: 2023_09_21-13.55.25 RELr426 190762
# Run by lpacheco on Tue Sep  1 2026
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *                                                     
from abaqusConstants import *                                            
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()                                              
                                                                            
# 1. Open the specific ODB
o1 = session.openOdb(
    name='/home/lpacheco/UEL-dynamicCL/5_CL_DIFF/test_in_abaqus/cube_indent_uel_relax_eta067.odb')                                              
odb = session.odbs['/home/lpacheco/UEL-dynamicCL/5_CL_DIFF/test_in_abaqus/cube_indent_uel_relax_eta067.odb']

# odb_path = '/home/lpacheco/UEL-dynamicCL/5_CL_DIFF/test_in_abaqus/cube_indent_uel_relax_eta067.odb'       
# odb = session.openOdb(name=odb_path)                                     
                                                                            
time_array = []                                                          
uvarm_data = {}                                                          
                                                                            
print("Extracting UVARM17 to UVARM106 history for Node 1...")            
                                                                            
# 2. Loop through directions 17 to 106                                   
for i in range(17, 107):                                                 
    uvarm_name = 'UVARM'+str(i)                                           
    xy = None                                                            
                                                                            
    try:                                                                 
        session.xyDataListFromField(odb=odb, outputPosition=NODAL, 
        variable=(('UVARM100', INTEGRATION_POINT), ), nodePick=(('PART-1-1', 1, ('[#1 ]', )), ), )
    except Exception as e:
        print("  Warning: Could not extract " + uvarm_name)
    
    if xy:
        data_tuples = xy[0].data
        
        # Grab the time column only once on the first successful pass    
        if not time_array:
            time_array = [pt[0] for pt in data_tuples]
            
        # Store the values
        uvarm_data[uvarm_name] = [pt[1] for pt in data_tuples]       
                                                                            
# 3. Save directly to CSV                                                
if time_array and uvarm_data:                                            
    out_file = odb_path.replace('.odb', '_cb_node1.csv')                 
                                                                            
    # Sort keys to guarantee exact order from 17 to 106                  
    uvarm_keys = ['UVARM%d' % i for i in range(17, 107) if 'UVARM%d' % i in uvarm_data]                                                             
                                                                            
    with open(out_file, 'w') as f:                                       
        # Write header                                                   
        header = ['Time'] + uvarm_keys                                   
        f.write(','.join(header) + '\n')                                 
                                                                            
        # Write data rows                                                
        for row_idx in range(len(time_array)):                           
            row = [str(time_array[row_idx])]                             
            for name in uvarm_keys:
                row.append(str(uvarm_data[name][row_idx]))
            f.write(','.join(row) + '\n')
            
    print("Successfully saved to: " + out_file)
else:
    print("Error: No data extracted.")

odb.close()
