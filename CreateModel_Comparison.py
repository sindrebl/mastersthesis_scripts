'''
This script is an example for how createModel_BasicScript.py can be used. 
It was used to do FE simulations on multiple designs of actual ICA, and
the data was used in section "Comparing FE model with ICA".
There are two ways to use it:

	1. Open commandwindow in the same folder and write:
		"abaqus cae noGUI=CreateModel_BasicScript.py"
		this will run the file without opening the GUI of ABAQUS.
	
	2. Open commandwindow in the same folder and write:
		"abaqus cae script=CreateModel_BasicScript.py"
		this will run the file inside the ABAQUS GUI.
'''


#======================================================================================
# 		INITIALIZE 
#======================================================================================
import os 

dir_path = 'F:/masteroppgave/script/scripts-tobe-published/ABAQUS/'
os.chdir(dir_path)
import myUtilities as myUt

dir_path = 'F:/masteroppgave/script/'
workDirectory = dir_path + 'workdir/';
if not os.path.exists(workDirectory):
    os.makedirs(workDirectory)	
os.chdir(workDirectory)

import numpy as np

'''===================================================================================
Remember to change textOutDirectory!!!
======================================================================================'''
textOutDirectory = dir_path + 'textOutDirectory/testPublish/';
if not os.path.exists(textOutDirectory):
    os.makedirs(textOutDirectory)	
	

#=======================================================================================
# 		GEOMETRIC PARAMTERER SPECIFICATION
#=======================================================================================
radiusVec = np.array([5.0, 5.0, 5.0, 5.0, 10.0, 15.0, 20.0]); 
shellThicknessVec = np.array([0.07, 0.14, 0.14, 0.2, 0.14, 0.14, 0.14]);
partConcVec = np.array([0.48, 0.48, 0.50, 0.48, 0.48, 0.48, 0.48]);


areRandomSpheres = False;
contactWidth = 0.5;#0.5;
maxLengthChain = 5; #nrOfSpheres
randomChains = False; #True or False;

meshSize = 0.25;
runInput = True;

#====================================================================================
#		MATERIAL
#====================================================================================

silver = myUt.Material('silver',10.49e-12,76.0e+6,0.37);
PMMA = myUt.Material('PMMA', 1.16e-12,5.37e+06,0.343); 
epoxy = myUt.Material('epoxy',1.12e-12,4.35e+06,0.368);
polystyrene = myUt.Material('polystyrene',1.05e-12,3.45e+06,0.358);


#coreMaterial = PMMA;
corePS = np.array([1,1,1,1,1,0,0])
shellMaterial = silver;
matrixMaterial = epoxy;

#=====================================================================================
# 		RUN FUNCTION
#=====================================================================================
for ss in range(len(radiusVec)):
	radius = radiusVec[ss];
	shellThickness = shellThicknessVec[ss];
	partConc = partConcVec[ss];
	
	if corePS[ss] == 1:
		coreMaterial = polystyrene;
		nrOfSpheres = int(np.ceil(215.0/(2*(2.0/(3*partConc))**(1.0/3)*(radius+shellThickness))));
	else:
		coreMaterial = PMMA;
		nrOfSpheres = int(np.ceil(257.0/(2*(2.0/(3*partConc))**(1.0/3)*(radius+shellThickness))));		
	
	sphereChainVector = (nrOfSpheres,); #has to start with value >1, 
										#(nrOfSpheres,)& (areRandomSpheres = True): creates a continuous connected model  
	
	myUt.commandPrint('---------------------CORE =' + coreMaterial.name + '----------------------------------');
	myUt.commandPrint('-----------------------N = ' + str(nrOfSpheres) + '-----------------------------------');		
	myUt.commandPrint('------------R = ' + str(radius) + '-----ST = ' + str(shellThickness) + '-----PC = ' + str(partConc) + '-----------');
	
	odbName, cell_rho, shellFraction, nrOfElem, randomSphereChainVector = myUt.createModel(coreMaterial,
		shellMaterial,matrixMaterial,radius,shellThickness,partConc,nrOfSpheres,
		contactWidth,maxLengthChain,randomChains,sphereChainVector, runInput,meshSize,
		areRandomSpheres,workDirectory);
	
	if runInput:
		tStrsMtrx, bStrsMtrx, distVect, timeVect, widthVect = myUt.odbWork(workDirectory,odbName,cell_rho);
		myUt.writeOutputTofile(textOutDirectory, odbName, shellFraction, tStrsMtrx, bStrsMtrx, distVect, timeVect,nrOfElem,widthVect, randomSphereChainVector);
	
print 'finished'
myUt.commandPrint('-----======================-----finished------======================------');
