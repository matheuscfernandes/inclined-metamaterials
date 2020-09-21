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
import os
import glob

def list_files(directory, extension):
    return list((f for f in os.listdir(directory) if f.endswith('.' + extension)))

def ExtractPV(odbName):
    odb=openOdb(path=odbName)

    VOL = session.XYDataFromHistory(name='VOL',
        odb=odb,
        outputVariableName='Fluid cavity volume: CVOL PI: rootAssembly Node 1 in NSET RP' )
    PRESS = session.XYDataFromHistory(name='PRESS',
        odb=odb,
        outputVariableName='Fluid cavity pressure: PCAV PI: rootAssembly Node 1 in NSET RP' )

    Volume_aux = VOL.data
    Pressure_aux = PRESS.data
    Volume = np.array([Volume_aux[j][1] for j in range(len(Volume_aux))])
    Pressure = np.array([Pressure_aux[j][1] for j in range(len(Pressure_aux))])
    PV = np.column_stack((Volume, Pressure))
    np.savetxt(file.split('.')[0]+'.pv',PV, delimiter=',')
    odb.close()

path='.'
Files=list_files(path,'odb')
for file in Files:
    ExtractPV(file)