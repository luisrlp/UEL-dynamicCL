import sys, getopt, os, string
import math
from odbAccess import *
from abaqusConstants import *
import numpy as np
from numpy import linalg as LA
from collections import defaultdict

#####
## Create vectors for dMudx
#####

def get_p1_vector(s):
    '''
    11, 22, 33, 12, 13, 23
    '''
    s = np.array([[s[0],s[3],s[4]],[s[3],s[1],s[5]],[s[4],s[5],s[2]]])
    w, v = LA.eigh(s)
    return v[:,np.argsort(w)[-1]]
def get_p0_scalar(s):
    '''
    11, 22, 33, 12, 13, 23
    '''
    #s = np.ScalarType([s[0]])
    return s
def norm(vector):
    """ Returns the norm (length) of the vector."""
    # note: this is a very hot function, hence the odd optimization
    # Unoptimized it is: return np.sqrt(np.sum(np.square(vector)))
    return np.sqrt(np.dot(vector, vector))

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / norm(vector)
def upgrade_if_necessary(job_id):
    odbPath = job_id + '.odb'
    new_odbPath = None
    print odbPath
    if isUpgradeRequiredForOdb(upgradeRequiredOdbPath=odbPath):
        print "Upgrade required"
        path,file = os.path.split(odbPath)
        file = 'upgraded_'+file
        new_odbPath = os.path.join(path,file)
        upgradeOdb(existingOdbPath=odbPath, upgradedOdbPath=new_odbPath)
        odbPath = new_odbPath
    else:
        print "Upgrade not required"
    return odbPath
def NewField(job_id):
        # retrieve steps from the odb
    odbPath = upgrade_if_necessary(job_id)
    odb = openOdb(path=odbPath)
    part_instance = odb.rootAssembly.instances['PART-1-1']#.elementSets['DUMMY_MESH']
    keys = odb.steps.keys()
        #if v1.sectionCategory.name['DUMMY_MATERIAL']:
        #    print v1.label
    #for instanceName in odb.rootAssembly.instances['PART-1-1'].elementSets.keys():
    #    print instanceName
    #    For each field output value in the last frame,
    #print the name, description, and type members.
    # lastFrame = odb.steps['static'].frames[-1]
    # for f in lastFrame.fieldOutputs.values():
    #     print f.name, ':', f.description
    #     print 'Type: ', f.type
    #     # For each location value, print the position.
    #     for loc in f.locations:
    #         print 'Position:',loc.position
    #     print
    # for v1 in part_instance.elements:
    #     if v1.sectionCategory.name == 'solid < DUMMY_MATERIAL >':
    #         print v1.label
    for key in keys:
        step = odb.steps[key]
        frameRepository = step.frames
        if len(frameRepository):
            for frame in frameRepository:
                print 'Id = %d, Time = %f\n'%(frame.frameId,frame.frameValue)
                fo = frame.fieldOutputs
                # Get fields from output database.
                #f0 = fo['SDV36']
                f1 = fo['UVARM33']
                f2 = fo['UVARM34']
                f3 = fo['UVARM35']
                #field0 = f0.getSubset(region=part_instance, position=INTEGRATION_POINT)
                #field1 = f1.getSubset(region=part_instance, position=INTEGRATION_POINT)
                #field2 = f2.getSubset(region=part_instance, position=INTEGRATION_POINT)
                #field3 = f3.getSubset(region=part_instance, position=INTEGRATION_POINT)
                # Compute new fields fields.
                FieldData = defaultdict(list)
                for v1 in f1.values:
                    FieldData.setdefault(v1.elementLabel,[]).append(v1.data)
                #print FieldData
                fData = FieldData.items()
                elementLabels = []
                elementData = []
                filDir = frame.FieldOutput(name='FDIR3', description='Filament direction', type=VECTOR,validInvariants=[MAGNITUDE,])
                elData = []
                for v1,v2,v3 in zip(f1.values,f2.values,f3.values):
                    elData.append([v1.data,v2.data,v3.data])
                #elLabels = [ v.label for v in part_instance.elements ]
                elLabels = [ v.label for v in part_instance.elements if v.sectionCategory.name == 'solid < DUMMY_MATERIAL >']
                #print elLabels
                #print elData
                filDir.addData(position=INTEGRATION_POINT, instance=part_instance, labels=elLabels, data=elData)
    odb.save()
    odb.close()
if __name__ == '__main__':
        # Get command line arguments.
        usage = "usage: abaqus python NewField.py <job name>"
        optlist, args = getopt.getopt(sys.argv[1:],'')
        JobID = args[0]
        odbPath = JobID + '.odb'
        print JobID
        print odbPath
        if not odbPath:
                print usage
                sys.exit(0)
        if not os.path.exists(odbPath):
                print "odb %s does not exist!" % odbPath
                sys.exit(0)
        NewField(JobID)
