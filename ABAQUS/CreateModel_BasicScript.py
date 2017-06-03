'''
This script is the basic script for all the ABAQUS models.It uses the funtions in
myUtilities.py to create material, generate an ABAQUS model, operate on the ODB-file
and to write the extracted data to file. There are two ways to use it:

	1. Open commandwindow in the same folder and write:
		"abaqus cae noGUI=CreateModel_BasicScript.py"
		this will run the file without opening the GUI of ABAQUS.
	
	2. Open commandwindow in the same folder and write:
		"abaqus cae script=CreateModel_BasicScript.py"
		this will run the file inside the ABAQUS GUI.
		
This will run an ABAQUS analysis and, if specified, write out the output in .txt files 
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
Remember to change textOutDirectory!!! This is where the output txt-file are written
======================================================================================'''
textOutDirectory = dir_path + 'textOutDirectory/testPublish/';
if not os.path.exists(textOutDirectory):
    os.makedirs(textOutDirectory)	
	

#=======================================================================================
# 		GEOMETRIC PARAMTERER SPECIFICATION
#=======================================================================================
radius = 15.0;
shellThickness = 0.14;
partConc = 0.50; #[%vol]

areRandomSpheres = False;
contactWidth = 0.5;
maxLengthChain = 5; #nrOfSpheres
randomChains = False; #True or False;


meshSize = 0.25;
runInput = True;

#====================================================================================
#		MATERIAL
#====================================================================================

#myUt.Material('name',density,Young's Modulus, Poisson's ratio)
silver = myUt.Material('silver',10.49e-12,76.0e+6,0.37);
PMMA = myUt.Material('PMMA', 1.16e-12,5.37e+06,0.343); 
epoxy = myUt.Material('epoxy',1.12e-12,4.35e+06,0.368);
polystyrene = myUt.Material('polystyrene',1.05e-12,3.45e+06,0.358);


coreMaterial = PMMA;
shellMaterial = silver;
matrixMaterial = epoxy;



if (coreMaterial.name == 'polystyrene'):
	nrOfSpheres = int(np.ceil(215.0/(2*(2.0/(3*partConc))**(1.0/3)*(radius+shellThickness))));
else:
	nrOfSpheres = int(np.ceil(257.0/(2*(2.0/(3*partConc))**(1.0/3)*(radius+shellThickness))));

sphereChainVector = (nrOfSpheres,); #(3,2,3,1,3,5,4,); #If randomChain = False. 
													   #Decides the ordering of chains
													   #has to start with value >1

#=====================================================================================
# 		RUN FUNCTION
#=====================================================================================

odbName, cell_rho, shellFraction, nrOfElem, randomSphereChainVector = myUt.createModel(coreMaterial,
	shellMaterial,matrixMaterial,radius,shellThickness,partConc,nrOfSpheres,
	contactWidth,maxLengthChain,randomChains,sphereChainVector, runInput,meshSize,
	areRandomSpheres,workDirectory);
	
if runInput:
	tStrsMtrx, bStrsMtrx, distVect, timeVect, widthVect = myUt.odbWork(workDirectory,odbName,cell_rho);
	myUt.writeOutputTofile(textOutDirectory, odbName, shellFraction, tStrsMtrx, bStrsMtrx, distVect, timeVect,nrOfElem,widthVect, randomSphereChainVector);
	
print 'finished' #Writes to ABAQUS GUI
myUt.commandPrint('----------finished------'); #Writes to command window
