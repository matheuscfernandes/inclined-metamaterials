from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import numpy as np
import random
random.seed(1111)

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

def DeleteAbaqusFiles(Job):
	try:
		os.remove(Job+'.odb')
	except: pass
	
	try:
		os.remove(Job+'.dat')
	except: pass
	
	try:
		os.remove(Job+'.com')
	except: pass
	
	try:
		os.remove(Job+'.ipm')
	except: pass
	
	try:
		os.remove(Job+'.log')
	except: pass
	
	try:
		os.remove(Job+'.prt')
	except: pass
	
	try:
		os.remove(Job+'.sim')
	except: pass
	
	try:
		os.remove(Job+'.sta')
	except: pass
		
	try:
		os.remove(Job+'.msg')
	except: pass
	
	try:
		os.remove(Job+'.lck')
	except: pass

	try:
		os.remove(Job+'.inp')
	except: pass

def DeleteAbaqusFilesButODB(Job):
	try:
		os.remove(Job+'.dat')
	except: pass
	
	try:
		os.remove(Job+'.com')
	except: pass
	
	try:
		os.remove(Job+'.ipm')
	except: pass
	
	try:
		os.remove(Job+'.log')
	except: pass
	
	try:
		os.remove(Job+'.prt')
	except: pass
	
	try:
		os.remove(Job+'.sim')
	except: pass
	
	try:
		os.remove(Job+'.sta')
	except: pass
		
	try:
		os.remove(Job+'.msg')
	except: pass
	
	try:
		os.remove(Job+'.lck')
	except: pass
	
	try:
		os.remove(Job+'.inp')
	except: pass

		
		
def ExtractValues(setname,JobName,OutputFile):  
	NameStep='Step-1'
	instancename='PART-1-1'
	text_file = open(OutputFile, 'w')
	odb=openOdb(path=JobName+'.odb')
	FrameLength=len(odb.steps[NameStep].frames)
	Frame=odb.steps[NameStep].frames[-1]
	Time=odb.steps[NameStep].frames[-1].frameValue
	set=odb.rootAssembly.instances[instancename].nodeSets[setname] 
	displacement=Frame.fieldOutputs['U'].getSubset(region=set).values
	
	for nod in range(0,len(set.nodes)):
				xDisp=displacement[nod].data[0]
				yDisp=displacement[nod].data[1]
				zDisp=displacement[nod].data[2]
				
				xCoor=set.nodes[nod].coordinates[0]
				yCoor=set.nodes[nod].coordinates[1]
				
				
				text_file.write('%e %e %e %e %e %e\r\n' % (Time, xCoor, yCoor, xDisp, yDisp, zDisp))
	text_file.close()
			
	
def ReturnRF(instancename,setname,JobName):  
	NameStep='Step-1'
	odb=openOdb(path=JobName+'.odb')
	FrameLength=len(odb.steps[NameStep].frames)
	Frame=odb.steps[NameStep].frames[FrameLength-1]
	Time=odb.steps[NameStep].frames[FrameLength-1].frameValue
	set=odb.rootAssembly.nodeSets[setname] 
	displacement=Frame.fieldOutputs['RF'].getSubset(region=set).values
	xRF=[]
	yRF=[]
	zRF=[]

	xCoor=[]
	yCoor=[]
	for nod in range(0,len(set.nodes)):
		xRF.append(displacement[nod].data[0])
		yRF.append(displacement[nod].data[1])
		zRF.append(displacement[nod].data[2])

	odb.close()
	return np.sum(xRF), np.sum(yRF), np.sum(zRF)

def ExtractRF(setname,JobName,OutputFile):  
	NameStep='Step-1'
	instancename='PART-1-1'
	text_file = open(OutputFile, 'w')
	odb=openOdb(path=JobName+'.odb')
	FrameLength=len(odb.steps[NameStep].frames)
	Frame=odb.steps[NameStep].frames[FrameLength-1]
	Time=odb.steps[NameStep].frames[FrameLength-1].frameValue
	set=odb.rootAssembly.instances[instancename].nodeSets[setname] 
	displacement=Frame.fieldOutputs['RF'].getSubset(region=set).values
	
	for nod in range(0,len(set.nodes)):
				xRF=displacement[nod].data[0]
				yRF=displacement[nod].data[1]
				zRF=displacement[nod].data[2]
				
				xCoor=set.nodes[nod].coordinates[0]
				yCoor=set.nodes[nod].coordinates[1]
				
				
				text_file.write('%e %e %e %e %e %e\r\n' % (Time, xCoor, yCoor, xRF, yRF, zRF))	
	text_file.close()
	odb.close()

def ExtractVirtualPointRF(instanceVP,nodeSetVP,stepName,JobName,OutputFile):
	text_file = open(OutputFile, 'w')
	odb=openOdb(path=JobName+'.odb')
	FrameLength=len(odb.steps[stepName].frames)
	set = odb.rootAssembly.instances[instanceVP].nodeSets[nodeSetVP]
	for i in range(FrameLength):
		Frame = odb.steps[stepName].frames[i].frameValue
		RF = odb.steps[stepName].frames[i].fieldOutputs['RF'].getSubset(region = set).values[0].data[1]
		text_file.write('%e %e\r\n' % (Frame, RF))
	text_file.close()
	odb.close()

def ApplyBuckling(mdbName,odbName,ImpFrames,ImpWeights,StepName,PartName,ModelName,InstanceName,AllNodes):
	########## DESCRIPTION OF VARIABLES FOR THE DEFINED FUNCTION
	'''
	mdbName -- 	THIS IS THE NAME FOR THE MDB PATH FOR THE CAE
	odbName -- THIS IS THE PATH FOR THE ODB FILE
	ImpFrames -- THIS IS THE FRAME NUMBERS OF THE FREQUENCY NUMBER OF THE ANALYSIS
	ImpWeigts -- THIS IS THE WEIGHT FOR THE FREQUENCY DISPLACEMENT B.C.
	StepName -- THIS IS THE NAME OF THE STEP
	PartName -- THIS IS THE NAME OF THE PART WHERE THE BUCKLING WILL BE APPLIED
	ModelName -- THIS IS THE MODEL NAME IN THE CAE
	InstanceName -- THIS IS THE INSTANCE NAME IN THE CAE
	AllNodes -- THIS IS THE SET NAME THAT CONTAINS ALL OF THE NODES
	'''
	##########
	# BEGIN CODE BY IMPORTING THE MDB MODEL OF FREQUENCY ANALYSIS
	mdb=openMdb(pathName=mdbName)
	odb=openOdb(path=odbName)
	# THIS OPENS THE RESULTS ODB WITH THE NODESET OF ALL OF THE NODES
	pbpPartNodes=odb.rootAssembly.instances[InstanceName.upper()].nodeSets[AllNodes.upper()]

	#CREATE A MATRIX TO SAVE THE NEW COORDINATES OF ALL NODES
	NewCoord=np.zeros((len(mdb.models[ModelName].parts[PartName].nodes), 3))	
	
	for CImp in range(len(ImpFrames)):
		cframe = ImpFrames[CImp]
		firstFrame = odb.steps[StepName].frames[cframe]
		displacement = firstFrame.fieldOutputs['U']
		pbpDispField = displacement.getSubset(region=pbpPartNodes)
		pbpDisp = pbpDispField.values


		# Imperfection Using Buckling Analysis Results
		#---------------------------------------------------------------
		ind=0;
		IMP = ImpWeights[CImp]
		#CALCULATE THE MODIFIED COORDINATES
		for i in mdb.models[ModelName].parts[PartName].nodes:
			NewCoord[ind][0]=i.coordinates[0]+IMP*pbpDisp[ind].data[0]
			NewCoord[ind][1]=i.coordinates[1]+IMP*pbpDisp[ind].data[1]
			#NewCoord[ind][2]=i.coordinates[2]+IMP*pbpDisp[ind].data[2]
			ind=ind+1

		#SET THE NEW COORDINATES
		mdb.models[ModelName].parts[PartName].editNode(
			nodes=mdb.models[ModelName].parts[PartName].nodes,
			coordinates=NewCoord)

	mdb.models[ModelName].rootAssembly.regenerate()

	print 'Original Coordinates Modified Successfully!'
	
def ExtractEigenMode(JobName, NumberOfModes):
    odb = openOdb(path=JobName + '.odb')
    Freq = []
    for i in xrange(1, NumberOfModes + 1):
        Desc = odb.getFrame(i).description
        Desc = Desc.split("=")
        Freq.append(float(Desc[1]))
    odb.close()
    return Freq

def ApplyRandomImperfection(mdb,IMP,PartName,ModelName):
	########## DESCRIPTION OF VARIABLES FOR THE DEFINED FUNCTION
	'''
	IMP -- THIS IS THE WEIGHT FOR THE FREQUENCY DISPLACEMENT B.C.
	StepName -- THIS IS THE NAME OF THE STEP
	PartName -- THIS IS THE NAME OF THE PART WHERE THE BUCKLING WILL BE APPLIED
	ModelName -- THIS IS THE MODEL NAME IN THE CAE
	InstanceName -- THIS IS THE INSTANCE NAME IN THE CAE
	AllNodes -- THIS IS THE SET NAME THAT CONTAINS ALL OF THE NODES
	'''
	##########
	# BEGIN CODE BY IMPORTING THE MDB MODEL OF FREQUENCY ANALYSIS

	#CREATE A MATRIX TO SAVE THE NEW COORDINATES OF ALL NODES
	NewCoord=np.zeros((len(mdb.models[ModelName].parts[PartName].nodes), 3))	
	
	# Imperfection Using Buckling Analysis Results
	#---------------------------------------------------------------
	ind=0
	#CALCULATE THE MODIFIED COORDINATES
	for i in mdb.models[ModelName].parts[PartName].nodes:
		NewCoord[ind][0]=i.coordinates[0]+IMP*(0.5-random.random())
		NewCoord[ind][1]=i.coordinates[1]+IMP*(0.5-random.random())
		NewCoord[ind][2]=i.coordinates[2]#+IMP*(0.5-random.random())
		ind=ind+1

	#SET THE NEW COORDINATES
	mdb.models[ModelName].parts[PartName].editNode(
		nodes=mdb.models[ModelName].parts[PartName].nodes,
		coordinates=NewCoord)

	mdb.models[ModelName].rootAssembly.regenerate()

	print 'Original Coordinates Modified Successfully!'

def ExtractSetRFAndU3(instance,nodeSet,stepName,JobName,OutputFile):
	text_file = open(OutputFile, 'w')
	odb=openOdb(path=JobName+'.odb')
	FrameLength=len(odb.steps[stepName].frames)
	set = odb.rootAssembly.instances[instance].nodeSets[nodeSet]
	for i in range(FrameLength):
		U3 = odb.steps[stepName].frames[i].fieldOutputs['U'].getSubset(region = set).values[0].data[2]
		RF3 = odb.steps[stepName].frames[i].fieldOutputs['RF'].getSubset(region = set).values[0].data[2]
		text_file.write('%e %e\r\n' % (U3, RF3))
	text_file.close()
	odb.close()

def ExtractSetMaxDistance(instance,stepName,JobName,nodeSet1,nodeSet2,OutputFile,Direction):
	text_file = open(OutputFile+'_Output.txt', 'w')
	odb=openOdb(path=JobName+'.odb')

	set1 = odb.rootAssembly.nodeSets[nodeSet1]
	set2 = odb.rootAssembly.nodeSets[nodeSet2]
	
	i=0
	for s1 in odb.steps[stepName].frames[-1].fieldOutputs['U'].getSubset(region = set1).values:
		vals=s1.data[Direction]
		coor=odb.rootAssembly.nodeSets[nodeSet1].nodes[0][i].coordinates
		text_file.write('%s %e %e %e %e\r\n' % (nodeSet1, vals,coor[0],coor[1],coor[2]))
		i+=1
	i2=0
	for s2 in odb.steps[stepName].frames[-1].fieldOutputs['U'].getSubset(region = set2).values:
		vals2=s2.data[Direction]
		coor2=odb.rootAssembly.nodeSets[nodeSet2].nodes[0][i2].coordinates
		text_file.write('%s %e %e %e %e\r\n' % (nodeSet2, vals2,coor2[0],coor2[1],coor2[2]))
		i2+=1


	text_file.close()
	
	odb.close()

def PeriodicBound2D(mdb, NameModel, NameSet, LatticeVec):
    from part import TWO_D_PLANAR, DEFORMABLE_BODY
    # Create reference parts and assemble
    NameRef1 = 'RefPoint-0'
    NameRef2 = 'RefPoint-1'
    mdb.models[NameModel].Part(
        dimensionality=TWO_D_PLANAR, name=NameRef1, type=DEFORMABLE_BODY)
    mdb.models[NameModel].parts[NameRef1].ReferencePoint(point=(0.0, 0.0, 0.0))
    mdb.models[NameModel].Part(
        dimensionality=TWO_D_PLANAR, name=NameRef2, type=DEFORMABLE_BODY)
    mdb.models[NameModel].parts[NameRef2].ReferencePoint(point=(0.0, 0.0, 0.0))
    mdb.models[NameModel].rootAssembly.Instance(dependent=ON, name=NameRef1,
                                                part=mdb.models[NameModel].parts[NameRef1])
    mdb.models[NameModel].rootAssembly.Instance(dependent=ON, name=NameRef2,
                                                part=mdb.models[NameModel].parts[NameRef2])

    # Create set of reference points
    mdb.models[NameModel].rootAssembly.Set(name=NameRef1, referencePoints=(
        mdb.models[NameModel].rootAssembly.instances[NameRef1].referencePoints[1],))
    mdb.models[NameModel].rootAssembly.Set(name=NameRef2, referencePoints=(
        mdb.models[NameModel].rootAssembly.instances[NameRef2].referencePoints[1],))

    # Get all nodes
    nodesAll = mdb.models[NameModel].rootAssembly.sets[NameSet].nodes
    nodesAllCoor = []
    for nod in mdb.models[NameModel].rootAssembly.sets[NameSet].nodes:
        nodesAllCoor.append(nod.coordinates)
    repConst = 0
    # Find periodically located nodes and apply equation constraints
    # Index array of nodes not used in equations constraint
    ranNodes = range(0, len(nodesAll))
    for repnod1 in xrange(0, len(nodesAll)):
        stop = False  # Stop will become true when equation constraint is made between nodes
        # Select Node 1 for possible equation constraint
        nod1 = nodesAll[repnod1]
        Coor1 = nodesAllCoor[repnod1]  # Coordinates of Node 1
        for repnod2 in ranNodes:  # Loop over all available nodes
            # Select Node 2 for possible equation constraint
            nod2 = nodesAll[repnod2]
            Coor2 = nodesAllCoor[repnod2]  # Coordinates of Node 2
            dx = Coor2[0] - Coor1[0]
            dy = Coor2[1] - Coor1[1]  # X and Y Distance between nodes
            # Check if nodes are located exactly the vector lattice apart
            for comb in xrange(0, len(LatticeVec)):
                if int(round(1000.0 * (LatticeVec[comb][0] - dx))) == 0:
                    if stop:
                        break
                    elif int(round(1000.0 * (LatticeVec[comb][1] - dy))) == 0:
                        # Correct combination found
                        # Create sets for use in equations constraints
                        mdb.models[NameModel].rootAssembly.Set(name='Node-1-' + str(
                            repConst), nodes=mdb.models[NameModel].rootAssembly.sets[NameSet].nodes[repnod1:repnod1 + 1])
                        mdb.models[NameModel].rootAssembly.Set(name='Node-2-' + str(
                            repConst), nodes=mdb.models[NameModel].rootAssembly.sets[NameSet].nodes[repnod2:repnod2 + 1])
                        # Create equations constraints for each dof
                        for Dim1 in [1, 2]:
                            mdb.models[NameModel].Equation(name='PerConst' + str(Dim1) + '-' + str(repConst),
                                                           terms=((1.0, 'Node-1-' + str(repConst), Dim1), (-1.0, 'Node-2-' + str(repConst), Dim1),
                                                                  (dx, 'RefPoint-' + str(Dim1 - 1), 1), (dy, 'RefPoint-' + str(Dim1 - 1), 2)))
                        mdb.models[NameModel].Equation(name='Rotation' + '-' + str(repConst), terms=(
                            (1.0, 'Node-1-' + str(repConst), 6),
                            (-1.0, 'Node-2-' + str(repConst), 6)))
                        repConst = repConst + 1  # Increase integer for naming equation constraint
                        # Remove used node from available list
                        ranNodes.remove(repnod1)
                        # Don't look further, go to following node.
                        stop = True
                        break

    # Return coordinates of free node so that it can be fixed
    return (NameRef1, NameRef2, repConst)

def ApplyAverageConstraint(NameModel,Dim,NameSet):
	nodesAll = mdb.models[NameModel].rootAssembly.sets[NameSet].nodes
	ct=0
	for i in xrange(0,len(nodesAll)):
		mdb.models[NameModel].rootAssembly.Set(name='AVGNode'+NameSet+str(i), 
			nodes=mdb.models[NameModel].rootAssembly.sets[NameSet].nodes[i:i + 1])
	listt=[]
	for i in xrange(0,len(nodesAll)):
		listt.append((1.0, 'AVGNode'+NameSet+str(i), Dim))

	mdb.models[NameModel].Equation(name='AVGAll'+NameSet, terms=tuple(listt))


# ApplyAverageConstraint('Model-1',3,'BOTTOM')	
# ExtractSetMaxDistance('PART-3-1','Step-1','Test','BOTTOM','TOP','Test',2)