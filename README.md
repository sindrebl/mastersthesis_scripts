# mastersthesis_scripts
The python-scripts used to create FE models are found in the folder "ABAQUS", and the 
scripts used for post processing are found in the folder "MATLAB". 

 ## ABAQUS 

If ABAQUS is installed on the computer, the scripts can be run by opening a command window 
in the folder (shift + right mouse button -> "Open command window here") and using the
commands:

	abaqus cae nogui=filename.py 
		This command will run the file without opening the ABAQUS GUI.	
	abaqus cae script=filename.py
		This command will run the file inside the ABAQUS GUI.

Make sure to edit the directories used in the script before running them. The author of the 
scripts has used ABAQUS 6.14-1, and cannot ensure that the files will work for other editions.

### The scripts:

	- CreateModel_BasicScript.py
		  The most basic script used to create a model and run simulations. 	
		  It can be used as a starting point for every simulation.
	- CreateModel_Comparison.py
		  An example of how the basic script can be used to create multiple models. 		
		  It is the actual script used in the comparison study in the master's 		
		  thesis.	
	- myUtilities.py
		  The functions used by "CreateModel_BasicScript.py". An explanation for every 
		  function can be found inside the script. 

## MATLAB 
Make sure to edit directories.

### The script:

	- ABAQUS_study_NoPlot_Basic.m
		  The basic script used for post processing of ABAQUS data. The user chooses 
		  the directory that has the ABAQUS output files. It then reads the infoFile.txt
		  and calculates the acoustic impedance.

### The supporting functions:
An explanation for every function are found inside the functions.

	- density_eff.m
		  Calculates the effective density of a two-phase composite.
	- homogenization_threePhase.m
		  Approximates the effective homogeneous three-phase composite using 
		  "threePhaseModel_bulk.m ", "density_eff.m" and "ThreePhaseModel_shear.m".
	- threePhaseModel_bulk.m 
		  Calculates the effective bulk modulus of a two-phase composite using the 
		  TPM as described in section 2.4.3. 	  
	- threePhaseModel_shear.m
		  Calculates the effective shear modulus of a two-phase composite using the
		  TPM as described in section 2.5.4.
	- xcorr_shift.m
		  Calculates the lag used in the Time Delay Estimation described in section 3.6.
		  It returns what is inside the brackets of eq. (3.17). 		
