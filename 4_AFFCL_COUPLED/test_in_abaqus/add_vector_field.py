# import sys, getopt, os, string
# import math
# from odbAccess import *
# from abaqusConstants import *
# import numpy as np
# from numpy import linalg as LA
# from collections import defaultdict
import sys, os
from odbAccess import openOdb, isUpgradeRequiredForOdb, upgradeOdb
from abaqusConstants import VECTOR, MAGNITUDE, INTEGRATION_POINT

def upgrade_if_necessary(job_id):
    odbPath = job_id + '.odb'
    print(f"Checking ODB: {odbPath}")
    if isUpgradeRequiredForOdb(upgradeRequiredOdbPath=odbPath):
        print("Upgrade required. Upgrading...")
        path, file = os.path.split(odbPath)
        new_odbPath = os.path.join(path, 'upgraded_' + file)
        upgradeOdb(existingOdbPath=odbPath, upgradedOdbPath=new_odbPath)
        return new_odbPath
    return odbPath

def add_vector_fields(job_id):
    odbPath = upgrade_if_necessary(job_id)
    odb = openOdb(path=odbPath, readOnly=False)
    root_assembly = odb.rootAssembly

    # 1. Print Instance Info (Optional, kept for your debugging)
    if root_assembly.instances:
        for inst_name, instance in root_assembly.instances.items():
            print(f"Instance: '{inst_name}', Total Elements: {len(instance.elements)}")
    else:
        print("No instances found in the root assembly.")

    # 2. Define the new fields to be created
    new_fields = [
        {
            'name': 'MUGRAD',
            'desc': 'Gradient of the Chemical Potential',
            'uvarms': ('UVARM10', 'UVARM11', 'UVARM12')
        },
        {
            'name': 'JFLUX',
            'desc': 'Fluid Flux',
            'uvarms': ('UVARM13', 'UVARM14', 'UVARM15')
        }
    ]
    for step_name, step in odb.steps.items():
        print(f"\n==================================================")
        print(f"Processing Step: '{step_name}'")
        print(f"==================================================")
    # step = odb.steps['indentation']
    
        # 3. Process each frame
        for frame in step.frames:
            fo = frame.fieldOutputs
            
            # Fast check: ensure all required UVARMs exist in this frame
            all_req_uvarms = [uvarm for field in new_fields for uvarm in field['uvarms']]
            if not all(name in fo for name in all_req_uvarms):
                continue

            print(f"\nProcessing Frame Id = {frame.frameId}, Time = {frame.frameValue:.4f}")

            # Process each vector field defined in our configuration
            for field in new_fields:
                f_name = field['name']
                
                if f_name in fo:
                    print(f"  Skipping {f_name}: already exists.")
                    continue

                # Extract the three component arrays
                c1, c2, c3 = [fo[uv].values for uv in field['uvarms']]

                if not (len(c1) == len(c2) == len(c3)):
                    raise RuntimeError(f"Mismatched UVARM counts for {f_name} in frame {frame.frameId}.")

                if not c1:
                    continue

                # Dictionaries to group data by instance
                labels_by_inst = {}
                data_by_inst = {}

                # Map the data
                for v1, v2, v3 in zip(c1, c2, c3):
                    inst_name = v1.instance.name
                    
                    if inst_name not in labels_by_inst:
                        labels_by_inst[inst_name] = []
                        data_by_inst[inst_name] = []

                    # Append element label only once per element
                    if not labels_by_inst[inst_name] or labels_by_inst[inst_name][-1] != v1.elementLabel:
                        labels_by_inst[inst_name].append(v1.elementLabel)
                    
                    # Append the vector tuple for this integration point
                    data_by_inst[inst_name].append((float(v1.data), float(v2.data), float(v3.data)))

                # Create the new field in the ODB
                new_vector_field = frame.FieldOutput(
                    name=f_name,
                    description=field['desc'],
                    type=VECTOR,
                    validInvariants=[MAGNITUDE]
                )

                # Write data back grouped by instance
                for inst_name, labels in labels_by_inst.items():
                    data = data_by_inst[inst_name]
                    print(f"  Adding {f_name} to '{inst_name}': {len(labels)} elements, {len(data)} vectors.")
                    
                    new_vector_field.addData(
                        position=INTEGRATION_POINT,
                        instance=root_assembly.instances[inst_name],
                        labels=labels,
                        data=data
                    )

    odb.save()
    odb.close()
    print("\nScript completed successfully.")

if __name__ == '__main__':
    usage = "usage: abaqus python add_vector_field.py <job name>"
    
    if len(sys.argv) < 2:
        print(usage)
        sys.exit(0)
        
    # Safely handle the job name whether the user includes ".odb" or not
    JobID = sys.argv[1].replace('.odb', '')
    odbPath = JobID + '.odb'
    
    if not os.path.exists(odbPath):
        print(f"Error: {odbPath} does not exist!")
        sys.exit(0)
        
    add_vector_fields(JobID)