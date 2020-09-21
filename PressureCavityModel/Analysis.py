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
import math

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

execfile('Functions.py')

def RunSimulation(ANGLECUT):
    Mdb()

    TOL=10E-3
    
    #PARAMETERS
    RHOLES=11.66/2.
    
    XNHOLES=8
    YNHOLES=8

    XSEPHOLES=3.*11.66/20.
    YSEPHOLES=3.*11.66/20.

    CHANNELDIAMETER=6.

    CHANNELHEIGHT=4.
    CHANNELWIDTH=1.5
    FILLETRADIUS=0.5

    FILMTHICKNESS=.5

    # ANGLECUT=15
    
    THICKNESS=11.66*2.#/np.sin(np.deg2rad(ANGLECUT))

    #MATERIAL PARAMETERS
    DENSITY=1.07E-9
    ADAMPING=15.
    BDAMPING=0.
    C10MAT=0.1875
    D1MAT=0.002061

    #COMPUTED PARAMETERS - Doumentation on these paramters are found in Blue Notebook
    XDOMAIN=XNHOLES*(RHOLES*2+XSEPHOLES)-2*RHOLES-XSEPHOLES
    YDOMAIN=YNHOLES*(RHOLES*2+YSEPHOLES)-2*RHOLES-YSEPHOLES
    # ANGLECUT=90-ANGLECUT

    LX=XDOMAIN/np.tan(np.deg2rad(ANGLECUT))
    LT=THICKNESS/np.sin(np.deg2rad(ANGLECUT))
    L=LT+LX
    LS=np.sqrt(LT**2-THICKNESS**2)
    XNEW=XDOMAIN/np.sin(np.deg2rad(ANGLECUT))


    BY=XDOMAIN*np.cos(np.deg2rad(ANGLECUT))
    BX=sqrt(XDOMAIN**2-BY**2)+sqrt(LX**2-BY**2)

    # DEFINING MATERIAL PROPERTIES
    mdb.models['Model-1'].Material(name='Material-1')
    mdb.models['Model-1'].materials['Material-1'].Hyperelastic(materialType=
        ISOTROPIC, table=((C10MAT, D1MAT), ), testData=OFF, type=NEO_HOOKE, 
        volumetricResponse=VOLUMETRIC_DATA)
    mdb.models['Model-1'].materials['Material-1'].Density(table=((DENSITY, ), ))
    mdb.models['Model-1'].materials['Material-1'].Damping(alpha=ADAMPING, beta=BDAMPING)

    #CREATING THE GEOMETRY
    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)

    mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(0, 0), 
        point2=(XDOMAIN, YDOMAIN))
    mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-Square', type=
        DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['Part-Square'].BaseSolidExtrude(depth=L, sketch=
        mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']

    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)

    circ=mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
        RHOLES+XSEPHOLES, RHOLES+YSEPHOLES), point1=(XSEPHOLES, RHOLES+YSEPHOLES))

    mdb.models['Model-1'].sketches['__profile__'].linearPattern(angle1=0.0, angle2=
        90.0, geomList=(circ, 
        ), number1=XNHOLES, number2=YNHOLES, spacing1=XSEPHOLES+2*RHOLES, spacing2=YSEPHOLES+2*RHOLES, vertexList=())
    
    mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-Holes', type=
        DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['Part-Holes'].BaseSolidExtrude(depth=L, sketch=
        mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']

    mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models['Model-1'].rootAssembly.Instance(dependent=OFF, name='Part-Holes-1', 
    part=mdb.models['Model-1'].parts['Part-Holes'])
    mdb.models['Model-1'].rootAssembly.Instance(dependent=OFF, name='Part-Square-1', 
    part=mdb.models['Model-1'].parts['Part-Square'])
    mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-Holes-1', ), 
        vector=(-RHOLES-XSEPHOLES, -RHOLES-YSEPHOLES, 0.0))  
    mdb.models['Model-1'].rootAssembly.InstanceFromBooleanCut(cuttingInstances=(
        mdb.models['Model-1'].rootAssembly.instances['Part-Holes-1'], ), 
            instanceToBeCut=
        mdb.models['Model-1'].rootAssembly.instances['Part-Square-1'], name=
            'Part-1', originalInstances=SUPPRESS)
        

    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)

    factor=BY

    mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(0, 0), 
        point2=(BX, YDOMAIN))
    
    mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-2', type=
        DEFORMABLE_BODY)

    mdb.models['Model-1'].parts['Part-2'].BaseSolidExtrude(depth=factor, sketch=
        mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']


    mdb.models['Model-1'].rootAssembly.Instance(dependent=OFF, name='Part-2-1', 
        part=mdb.models['Model-1'].parts['Part-2'])

    mdb.models['Model-1'].rootAssembly.Instance(dependent=OFF, name='Part-2-2', 
        part=mdb.models['Model-1'].parts['Part-2'])

    mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-2-1', ), 
        vector=((XDOMAIN-BX), 0.0, -BY)) 
   
    mdb.models['Model-1'].rootAssembly.rotate(angle=90-ANGLECUT, axisDirection=(0.0, 1.0, 
        0.0), axisPoint=(XDOMAIN, 0.0, 0.0), instanceList=('Part-2-1',))

    mdb.models['Model-1'].rootAssembly.rotate(angle=90-ANGLECUT, axisDirection=(0.0, 1.0, 
        0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('Part-2-2',))

    mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-2-2', ), 
        vector=(0.0, 0.0, L))   
    
    mdb.models['Model-1'].rootAssembly.InstanceFromBooleanCut(cuttingInstances=(
        mdb.models['Model-1'].rootAssembly.instances['Part-2-1'], 
        mdb.models['Model-1'].rootAssembly.instances['Part-2-2']), instanceToBeCut=
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'], name='Part-3', 
        originalInstances=DELETE)

    mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-3-1', ), 
        vector=(0.0, 0.0, -LX))   

    mdb.models['Model-1'].rootAssembly.rotate(angle=90-ANGLECUT, axisDirection=(0.0, -1.0, 
        0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('Part-3-1',))

    mdb.models['Model-1'].rootAssembly.features.changeKey(fromName='Part-3-1', toName='Part-3-2')
    

    #CREATING TOP AND BOTTOM FILMS
    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
    mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(0.0, 0.0), 
        point2=(XNEW, YDOMAIN))
    mdb.models['Model-1'].Part(dimensionality=THREE_D, name='FILM', type=
        DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['FILM'].BaseShell(sketch=
        mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']
    
    mdb.models['Model-1'].HomogeneousShellSection(idealization=NO_IDEALIZATION, 
        integrationRule=SIMPSON, material='Material-1', name='Section-film', 
        nodalThicknessField='', numIntPts=5, poissonDefinition=DEFAULT, 
        preIntegrate=OFF, temperature=GRADIENT, thickness=FILMTHICKNESS, thicknessField='', 
        thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
    mdb.models['Model-1'].parts['FILM'].Set(faces=
        (mdb.models['Model-1'].parts['FILM'].faces, ), name='Set-1')
    mdb.models['Model-1'].parts['FILM'].SectionAssignment(offset=0.0, offsetField=
        '', offsetType=MIDDLE_SURFACE, region=
        mdb.models['Model-1'].parts['FILM'].sets['Set-1'], sectionName='Section-film'
        , thicknessAssignment=FROM_SECTION)

    # mdb.models['Model-1'].MembraneSection(material='Material-1', name=
    #     'Section-membrane', poissonDefinition=DEFAULT, thickness=FILMTHICKNESS, 
    #     thicknessField='', thicknessType=UNIFORM)
    # mdb.models['Model-1'].parts['FILM'].SectionAssignment(offset=0.0, offsetField=
    #     '', offsetType=MIDDLE_SURFACE, region=
    #     mdb.models['Model-1'].parts['FILM'].sets['Set-1'], sectionName='Section-membrane'
    #     , thicknessAssignment=FROM_SECTION)

    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='FILM-1', part=
        mdb.models['Model-1'].parts['FILM'])
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='FILM-2', part=
        mdb.models['Model-1'].parts['FILM'])
    mdb.models['Model-1'].rootAssembly.translate(instanceList=('FILM-2', ), vector=
        (-LS, 0.0, THICKNESS))
    
    # mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
    #     instances=(mdb.models['Model-1'].rootAssembly.instances['FILM-2'], 
    #     mdb.models['Model-1'].rootAssembly.instances['Part-3-2'], 
    #     mdb.models['Model-1'].rootAssembly.instances['FILM-1']), name='Part-4', 
    #     originalInstances=DELETE)
    mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
        instances=(mdb.models['Model-1'].rootAssembly.instances['Part-3-2'], 
        mdb.models['Model-1'].rootAssembly.instances['FILM-1'], 
        mdb.models['Model-1'].rootAssembly.instances['FILM-2']), name='Part-4', 
        originalInstances=SUPPRESS)
    
    mdb.models['Model-1'].parts.changeKey(fromName='Part-3', toName='Part-3-old')
    mdb.models['Model-1'].parts.changeKey(fromName='Part-4', toName='Part-3')
    
    # CREATING SETS
    xSHIFT=THICKNESS/np.tan(np.deg2rad(ANGLECUT))

    mdb.models['Model-1'].rootAssembly.Set(name='Faces-Top', faces=
        mdb.models['Model-1'].rootAssembly.instances['Part-4-1'].faces.getByBoundingBox(
         -xSHIFT+TOL,  + TOL, THICKNESS-TOL, XNEW-xSHIFT - TOL, YDOMAIN- TOL, THICKNESS+TOL))

    mdb.models['Model-1'].rootAssembly.rotate(angle=90-ANGLECUT, axisDirection=(0.0, 1.0, 
        0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('Part-4-1', ))

    mdb.models['Model-1'].rootAssembly.Set(name='Holes-All', faces=
        mdb.models['Model-1'].rootAssembly.instances['Part-4-1'].faces.getByBoundingBox(
         +TOL,  + TOL, -100000, XDOMAIN-TOL, YDOMAIN-TOL, 100000))
    
    mdb.models['Model-1'].rootAssembly.SetByBoolean(name='Holes-OK', operation=
        DIFFERENCE, sets=(mdb.models['Model-1'].rootAssembly.sets['Holes-All'], 
        mdb.models['Model-1'].rootAssembly.sets['Faces-Top']))

    mdb.models['Model-1'].rootAssembly.Surface(name='Surf-1', side1Faces=
        mdb.models['Model-1'].rootAssembly.sets['Holes-OK'].faces,
         side2Faces=mdb.models['Model-1'].rootAssembly.sets['Faces-Top'].faces)
    
    mdb.models['Model-1'].rootAssembly.rotate(angle=-90+ANGLECUT, axisDirection=(0.0, 1.0, 
        0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('Part-4-1', ))

    # MESHING
    # mdb.models['Model-1'].parts['Part-3'].seedPart(minSizeFactor=0.01, size=2.*(1.+(85.-ANGLECUT)/(85.-10.)))
    mdb.models['Model-1'].parts['Part-3'].seedPart(minSizeFactor=0.005, size=1.6)
    mdb.models['Model-1'].parts['Part-3'].setMeshControls(elemShape=TET, regions=
        (mdb.models['Model-1'].parts['Part-3'].cells[0],), technique=FREE)

    # mdb.models['Model-1'].parts['Part-3'].setMeshControls(algorithm=ADVANCING_FRONT
    #     , elemShape=HEX_DOMINATED, regions=
    #     (mdb.models['Model-1'].parts['Part-3'].cells[0],), technique=SYSTEM_ASSIGN)

    # mdb.models['Model-1'].parts['Part-3'].setElementType(elemTypes=(ElemType(
    #     elemCode=C3D8R, elemLibrary=EXPLICIT), ElemType(elemCode=C3D6, 
    #     elemLibrary=EXPLICIT), ElemType(elemCode=C3D4, elemLibrary=EXPLICIT, 
    #     secondOrderAccuracy=OFF, distortionControl=DEFAULT)), regions=(
    #     mdb.models['Model-1'].parts['Part-3'].cells[0], ))
    # mdb.models['Model-1'].parts['Part-3'].setElementType(elemTypes=(ElemType(
    #     elemCode=C3D20R, elemLibrary=STANDARD), ElemType(elemCode=C3D15, 
    #     elemLibrary=STANDARD), ElemType(elemCode=C3D10, elemLibrary=STANDARD)), 
    #     regions=(mdb.models['Model-1'].parts['Part-3'].cells[0], ))
    mdb.models['Model-1'].parts['Part-3'].setElementType(elemTypes=(ElemType(elemCode=C3D4H, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, distortionControl=DEFAULT),), regions=(
        mdb.models['Model-1'].parts['Part-3'].cells[0],))

    mdb.models['Model-1'].parts['Part-3'].generateMesh()
    
    #CREATING STEP
    # mdb.models['Model-1'].ExplicitDynamicsStep(name='Step-1', previous='Initial', 
    #     timePeriod=1.0)
    # mdb.models['Model-1'].steps['Step-1'].setValues(improvedDtMethod=ON, 
    #     massScaling=((SEMI_AUTOMATIC, MODEL, AT_BEGINNING, 1000.0, 0.0, None, 0, 0, 
    #     0.0, 0.0, 0, None), ))
    mdb.models['Model-1'].ImplicitDynamicsStep(initialInc=0.01, maxInc=0.01, maxNumInc=1000, name='Step-1', nlgeom=ON, previous='Initial')


    ## CREATING AND WORKING WITH THE VIRTUAL NODES
    #CREATING SETS
    mdb.models['Model-1'].rootAssembly.regenerate()
    mdb.models['Model-1'].parts['Part-3'].Set(name='ALL',cells=mdb.models['Model-1'].parts['Part-3'].cells[:])

    mdb.models['Model-1'].rootAssembly.Set(name='BOTTOM', nodes=mdb.models['Model-1'].rootAssembly.instances['Part-4-1'].nodes.getByBoundingBox(
        0 - TOL, 0 - TOL, 0 -TOL, XNEW + TOL, YDOMAIN + TOL, 0 + TOL))

    mdb.models['Model-1'].rootAssembly.Set(name='TOP', nodes=mdb.models['Model-1'].rootAssembly.instances['Part-4-1'].nodes.getByBoundingBox(
        -LS - TOL, 0 - TOL, THICKNESS -TOL, XNEW -LS + TOL, YDOMAIN + TOL, THICKNESS + TOL))
    mdb.models['Model-1'].rootAssembly.Set(name='RIGHT', nodes=mdb.models['Model-1'].rootAssembly.instances['Part-4-1'].nodes.getByBoundingBox(
        -0 + TOL, 0 - TOL, -TOL, XNEW- TOL, 0 + TOL, 0+ TOL))
    mdb.models['Model-1'].rootAssembly.Set(name='LEFT', nodes=mdb.models['Model-1'].rootAssembly.instances['Part-4-1'].nodes.getByBoundingBox(
        -0 +TOL, YDOMAIN - TOL, -TOL, XNEW- TOL, YDOMAIN + TOL, 0+ TOL))


    #CREATE EQUATION CONSTRAING FOR ALL NODES AT BOTTOM FACE
    mdb.models['Model-1'].rootAssembly.Set(name='BOTTOM-AVG', nodes=mdb.models['Model-1'].rootAssembly.instances['Part-4-1'].nodes.getByBoundingBox(
        0 + TOL, 0 + TOL, 0 -TOL, XNEW - TOL, YDOMAIN - TOL, 0 + TOL))
    mdb.models['Model-1'].rootAssembly.Set(name='LEFT-AVG', nodes=mdb.models['Model-1'].rootAssembly.instances['Part-4-1'].nodes.getByBoundingBox(
        -LS -TOL, YDOMAIN - TOL, -TOL, XNEW+ TOL, YDOMAIN + TOL, THICKNESS+ TOL))
    
    ApplyRandomImperfection(mdb,0.2,'Part-3','Model-1')

    mdb.models['Model-1'].rootAssembly.regenerate()

    
    #APPLY BC
    # mdb.models['Model-1'].SmoothStepAmplitude(data=((0.0, 0.0), (1.0, 1.0)), name=
    #     'Amp-1', timeSpan=STEP)
    mdb.models['Model-1'].TabularAmplitude(data=((0.0, 0.0), (1.0, 1.0)), name=
        'Amp-1', smooth=SOLVER_DEFAULT, timeSpan=STEP)

    mdb.models['Model-1'].rootAssembly.ReferencePoint(point=(0.0, 0.0, 0.0))
    mdb.models['Model-1'].rootAssembly.Set(name='RP', referencePoints=(
        mdb.models['Model-1'].rootAssembly.referencePoints[mdb.models['Model-1'].rootAssembly.referencePoints.keys()[0]], ))

    mdb.models['Model-1'].Temperature(createStepName='Initial', 
        crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=
        UNIFORM, magnitudes=(0.0, ), name='Predefined Field-1', region=
        mdb.models['Model-1'].rootAssembly.sets['RP'])

    mdb.models['Model-1'].Temperature(amplitude='Amp-1', createStepName='Step-1', 
        crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=
        UNIFORM, magnitudes=(-0.9, ), name='Predefined Field-2', region=
        mdb.models['Model-1'].rootAssembly.sets['RP'])
    
    # DEFINING SECTION PROPERTIES
    mdb.models['Model-1'].HomogeneousSolidSection(material='Material-1', name=
        'Section-1', thickness=None)
    mdb.models['Model-1'].parts['Part-3'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=
        mdb.models['Model-1'].parts['Part-3'].sets['ALL'], sectionName=
        'Section-1', thicknessAssignment=FROM_SECTION)

    # HARD CONTACT DEFINITION
    mdb.models['Model-1'].ContactProperty('IntProp-1')
    mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
        allowSeparation=ON, constraintEnforcementMethod=DEFAULT, 
        pressureOverclosure=HARD)
    # mdb.models['Model-1'].ContactExp(createStepName='Step-1', name='Int-1')
    # mdb.models['Model-1'].interactions['Int-1'].includedPairs.setValuesInStep(
    #     stepName='Step-1', useAllstar=ON)
    # mdb.models['Model-1'].interactions['Int-1'].contactPropertyAssignments.appendInStep(
    #     assignments=((GLOBAL, SELF, 'IntProp-1'), ), stepName='Step-1')
    mdb.models['Model-1'].ContactStd(createStepName='Initial', name='Int-1')
    mdb.models['Model-1'].interactions['Int-1'].includedPairs.setValuesInStep(
        stepName='Initial', useAllstar=ON)
    mdb.models['Model-1'].interactions['Int-1'].contactPropertyAssignments.appendInStep(
        assignments=((GLOBAL, SELF, 'IntProp-1'), ), stepName='Initial')

    # FLUID CAVITY DEFINITION
    mdb.models['Model-1'].HistoryOutputRequest(createStepName='Step-1', name=
        'H-Output-2', frequency=1, rebar=EXCLUDE, region=
        mdb.models['Model-1'].rootAssembly.sets['RP'], sectionPoints=DEFAULT, 
        variables=('PCAV', 'CVOL'))
    mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(timeInterval=
        0.01)
    mdb.models['Model-1'].FluidCavityProperty(bulkModulusTable=((1000.0, ), ), 
        expansionTable=((0.33333333, ), ), fluidDensity=1.0, name='IntProp-2', 
        useBulkModulus=True, useExpansion=True)

    mdb.models['Model-1'].FluidCavity(cavityPoint=
        mdb.models['Model-1'].rootAssembly.sets['RP'], cavitySurface=
        mdb.models['Model-1'].rootAssembly.surfaces['Surf-1'], createStepName=
        'Initial', interactionProperty='IntProp-2', name='Int-2')

    # JOB INFORMATION
    JobName='Analysis_{}'.format(ANGLECUT)
    DeleteAbaqusFiles(JobName)
    
    mdb.models['Model-1'].rootAssembly.regenerate()
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
        multiprocessingMode=THREADS, name=JobName, nodalOutputPrecision=SINGLE, 
        numCpus=20, numDomains=20, numGPUs=0, queue=None, scratch='', type=ANALYSIS, 
        userSubroutine='', waitHours=0, waitMinutes=0)
    # mdb.jobs[JobName].setValues(explicitPrecision=DOUBLE_PLUS_PACK)
    mdb.saveAs(JobName+'.cae')

    asdf
    mdb.jobs[JobName].submit(consistencyChecking=OFF)
    mdb.jobs[JobName].waitForCompletion()

    ExtractSetMaxDistance('Part-3-1','Step-1',JobName,'BOTTOM','TOP',JobName,2)

    DeleteAbaqusFilesButODB(JobName)



RunSimulation(ANGLECUT)

# JobName='Analysis_{}'.format(ANGLECUT)
# ExtractSetMaxDistance('Part-3-1','Step-1',JobName,'BOTTOM','TOP',JobName,2)
