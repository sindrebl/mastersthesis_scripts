""" myUtilities - a module of functions to use in ABAQUS.
	FEM variables can be altered in "createModel - FEM variables".
"""

#======================================================================================
# 		initialization
#======================================================================================
from sketch import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from mesh import *
from load import *
from job import *
from visualization import *
from odbAccess import openOdb
import math
import numpy as np;
import time
import random
import sys

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

#==================================================================================
# 		functions and classes
#==================================================================================	

class Material:
	#Creates materials to be used as core, shell or matrix material in "createModel"
	def __init__(self,name,rho,youngs,poisson):
		self.name = name;
		self.rho = rho;
		self.youngs = youngs;
		self.poisson = poisson;
		self.bulk = youngs/(3*(1-2*poisson))
		self.shear = youngs/(2*(1+poisson));
		self.cL = math.sqrt((self.bulk+(4.0/3)*self.shear)/self.rho);
		

def commandPrint(toPrint):
	# Writes output to command window
	print >> sys.__stdout__, toPrint;
	
def createModel(coreMaterial,shellMaterial,matrixMaterial,radius,shellThickness,partConc,nrOfSpheres,contactWidth,maxLengthChain,randomChains,sphereChainVector, runInput,meshSize,areRandomSpheres,workDirectory):
	"""
	INPUT:
		coreMaterial,shellMaterial,matrixMaterial are Material objects
		radius: Radius of the core										[micro metre]
		shellThickness: Thickness of the shell 							[micro metre]
		partConc: volume percent of particles							[vol%]	
		nrOfSpheres: Particles included in the FE model 				[int]
		contactWidth: Radius of the contact zone 						[micro metre]
		maxLengthChain: Maximum chain length IF randomchain = True 		[int]
		randomChains: Particles are ordered in random chains			[True/False]
		sphereChainVector: Particle chains IF randomchain = False		(int,int,int,)
		runInput: Create and run job/analysis in ABAQUS					[True/False]	
		meshSize: Global mesh size										[micro metre]
		areRandomSpheres: Random (true) or ordered model(false)			[True/False]			
		workDirectory: 	Directory for ABAQUS files						[path]	
	OUTPUT:
		odbName: File name, used in "odbWork", includes geometry		[R_10_ST_01.....]
		cell_rho: density of composite, used in development				[g/cm^3]*1e-12 
		shellFraction: Shell fraction, used in development				[vol%]
		nrOfElements: Nr of elements in the FEmodel						[int]
		randomSphereChainVector: Vector of particle chains				(2,3,4,1,)
	"""
	#==================================================================================
	# 		FEM variables
	#==================================================================================	
	quadraticElements = 0; 				# set 1 or 0 for quad or linear triangular elements
	linBulkViscosity = 0.0; 			# used to reduce oscillations of stress
	quadrBulkViscosity = 0.0;			
	stressOutput = ('S33',);			# normal stress in z-direction
	infLayerHeight = 5.0;				
	totalTimePeriod = 8e-7;				
	partitionWedgeThreshold = 2.0; 		# Threshold for when to mesh wedge seperatly
	
	#==================================================================================
	# 		Geometric calculations
	#==================================================================================	
	volumeSphere = (4.0/3)*pi*(radius+shellThickness)**3;
	volumeCore = (4.0/3)*pi*radius**3;
	volumeShell = volumeSphere-volumeCore;
	width =  (volumeSphere/(2*partConc*pi))**(1 / 3.0); 
	cellHeight = width;
	sphere_rho = (volumeCore/volumeSphere)*coreMaterial.rho + (1-(volumeCore/volumeSphere))*shellMaterial.rho;
	cell_rho = partConc*sphere_rho + (1-partConc)*matrixMaterial.rho;
	shellFraction =(1-(volumeCore/volumeSphere))*partConc;
	
	#==================================================================================
	# 		Naming
	#==================================================================================	
	if quadraticElements == 1:
		type = 'AXQUAD'
	else:
		type = 'AXLIN'
		
	if areRandomSpheres:
		type += '_random';
		odbName = 'R'+str(radius)+'ST'+str(shellThickness)+'PC'+str(partConc)+'N'+str(nrOfSpheres) +'MS'+ str(meshSize)+'CW'+str(contactWidth)+'mCL'+str(maxLengthChain) + type;
	else:
		type += '_ordered';
		odbName = 'R'+str(radius)+'ST'+str(shellThickness)+'PC'+str(partConc)+'N'+str(nrOfSpheres) +'MS' +str(meshSize) + type;
	
	modelName = 'Model ' + type; 
	sketchName = 'Sketch ' + type;
	partName = 'Part ' + type;
	instanceName = partName + '-1';
	stepName = 'Loading';
	jobName = 'Job_' + type;
	#odbName = 'R'+str(radius)+'ST'+str(shellThickness)+'PC'+str(partConc)+'N'+str(nrOfSpheres)+'CW'+str(contactWidth)+'mCL'+str(maxLengthChain) + type;
	print odbName
	odbName = odbName.replace('.','_');
	print odbName

	#==================================================================================
	# 		Vector initialization
	#==================================================================================	
	timeVector=[];
	chainSystem=[];
	sphereCenterYVector = [];
	wedgeCoordYVector = [];

	#==================================================================================
	# 		Create model, materials and section
	#==================================================================================	
	
	# Model
	myModel = mdb.Model(name=modelName)

	# Materials 
	myCore = myModel.Material(name=coreMaterial.name)
	myCore.Elastic(table=((coreMaterial.youngs, coreMaterial.poisson), ))
	myCore.Density(table=((coreMaterial.rho, ), ))

	myShell = myModel.Material(name=shellMaterial.name)
	myShell.Elastic(table=((shellMaterial.youngs,shellMaterial.poisson), ))
	myShell.Density(table=((shellMaterial.rho, ), ))

	myMatrix = myModel.Material(name=matrixMaterial.name)
	myMatrix.Elastic(table=((matrixMaterial.youngs, matrixMaterial.poisson), ))
	myMatrix.Density(table=((matrixMaterial.rho, ), ))
	# Sections
	myModel.HomogeneousSolidSection(material=coreMaterial.name, 
									name='Core', thickness=None)
	myModel.HomogeneousSolidSection(material=matrixMaterial.name, 
									name='Matrix', thickness=None)
	myModel.HomogeneousSolidSection(material=shellMaterial.name, 
									name='Shell', thickness=None)
									
	#==================================================================================
	# 		Create parts
	#==================================================================================	
	
	# Initialization 
	mySketch = myModel.ConstrainedSketch(name='Sketch '+type, sheetSize=50.0)
	mySketch.sketchOptions.setValues(viewStyle=AXISYM)
	mySketch.ConstructionLine(point1=(0.0,-25.0), 
							  point2=(0.0, 25.0)
							  )					
	mySketch.rectangle(point1=(0.0,-(cellHeight+ infLayerHeight)), 
					   point2=(width,2*cellHeight*(nrOfSpheres-1)+cellHeight)
					   )					
	myPart = myModel.Part(name=partName, dimensionality=AXISYMMETRIC, type=DEFORMABLE_BODY)
	myPart.BaseShell(sketch=mySketch)
	
	
	#========================================================================================================
	# 		Partition spheres with shell
	#========================================================================================================
	sphereSketch = myModel.ConstrainedSketch(name='sphere', sheetSize=50.0)
	sphereSketch.sketchOptions.setValues(viewStyle=AXISYM)
	sphereSketch.ConstructionLine(point1=(0.0,-25.0), point2=(0.0, 25.0))


	if areRandomSpheres:
		#========================================================================================================
		# 		Partition random spheres with shell
		#========================================================================================================

		if contactWidth + meshSize/2 > partitionWedgeThreshold: # Wedge partition that gives elements of ~ normal mesh size 
			partitionWedgeThreshold = contactWidth
		elif partitionWedgeThreshold > contactWidth + meshSize:
			partitionWedgeThreshold = contactWidth + meshSize;
		partitionWedgeWidth = partitionWedgeThreshold + meshSize;

		matrixSpacingPerSphere = 2*(cellHeight-(radius+shellThickness));
		remainingSpheres = nrOfSpheres;
		bottomOfSystem = -cellHeight;
		nSphereChain = 0;
		randomSphereChainVector = [];
		while remainingSpheres > 0:
			if randomChains == True:
				if remainingSpheres < maxLengthChain: #makes sure generated spheres are not more than specified
					maxLengthChain = remainingSpheres;
				sphereChainLength = random.randint(1,maxLengthChain);
				if remainingSpheres - sphereChainLength == 1:
					sphereChainLength += 1; # last sphere is not a single sphere, to have spacing above
				while sphereChainLength == 1 and remainingSpheres == nrOfSpheres:
					# prevent first sphere to be a single sphere. This is because of the matrix layering 
					sphereChainLength = random.randint(1,maxLengthChain);
				randomSphereChainVector.append(sphereChainLength);
				
			else:
				randomSphereChainVector = sphereChainVector;
				sphereChainLength = sphereChainVector[nSphereChain];
				nSphereChain += 1;
			chainSystem.append(sphereChainLength);
			for n in np.arange(1,sphereChainLength+1,1):
				nrOfMatrixSpacings = sphereChainLength;
				
				if remainingSpheres == nrOfSpheres: # The first sphere. 
					# This has to have matrix spacing underneath,
					# to prevent reflection between sphere and inf layer
					nrOfMatrixSpacings -= 1;
					bottomOfSystem += matrixSpacingPerSphere;
					
				sphereCenterY = bottomOfSystem + (radius + shellThickness);
				sphereCenterYVector.append(sphereCenterY);
				
				#------------partition the core--------------------------------
				sphereSketch.ArcByCenterEnds(center=(0.0 ,sphereCenterY), 
											 direction=COUNTERCLOCKWISE, 
											 point1=(0.0, sphereCenterY-radius), 
											 point2=(0.0, sphereCenterY+radius))
				
				#------------------partitioning the shell--------------------------
				if sphereChainLength == 1:
					print 'single sphere'
					# a single sphere, not a chain
					sphereSketch.ArcByCenterEnds(center=(0.0, sphereCenterY),				
													direction=COUNTERCLOCKWISE, 
													point1=(0.0, sphereCenterY-(radius+shellThickness)), 
													point2=(0.0, sphereCenterY+(radius+shellThickness)))
					sphereSketch.Line(point1=(0.0, sphereCenterY-(radius+shellThickness)), 
													point2=(0.0, sphereCenterY+(radius+shellThickness)))
					# updating the start of next chain system
					bottomOfSystem += 2*(radius+shellThickness) + matrixSpacingPerSphere;
					remainingSpheres -= 1; #keeping track of total sphere count
				elif n == 1 and sphereChainLength != 1:
					print 'bottom sphere'
					# bottom sphere of chain system
					sphereSketch.ArcByCenterEnds(center=(0.0, sphereCenterY), 
												 direction=COUNTERCLOCKWISE, 
												 point1=(0.0, sphereCenterY-(radius+shellThickness)), 
												 point2=(contactWidth, sphereCenterY + sqrt((radius+shellThickness)**2 - contactWidth**2)))
					# upper connection zone
					sphereSketch.Line(point1=(contactWidth, sphereCenterY + sqrt((radius+shellThickness)**2 - contactWidth**2)), 
									  point2=(contactWidth, sphereCenterY +(radius+shellThickness)))
					
					#partition for wedge mesh
					if contactWidth < partitionWedgeThreshold:
						sphereSketch.Line (point1=(partitionWedgeWidth, sphereCenterY +sqrt((radius+shellThickness)**2 - partitionWedgeWidth**2)),
										point2=(partitionWedgeWidth, sphereCenterY +(radius+shellThickness)))
						wedgeCoordYVector.append(sphereCenterY + (radius + shellThickness));
					
					# axissymm symmetry line
					sphereSketch.Line(point1=(0.0, sphereCenterY-(radius+shellThickness)), 
									  point2=(0.0, sphereCenterY+(radius+shellThickness)))
					# updating the start of next sphere
					bottomOfSystem += 2*(radius+shellThickness);
					remainingSpheres -= 1; #keeping track of total sphere count
				elif n == sphereChainLength:
					print 'top sphere'
					#top sphere of chain system
					sphereSketch.ArcByCenterEnds(center=(0.0, sphereCenterY), 
												 direction=CLOCKWISE, 
												 point1=(0.0, sphereCenterY+(radius+shellThickness)), 
												 point2=(contactWidth, sphereCenterY - sqrt((radius+shellThickness)**2 - contactWidth**2)))
					#lower connection zone
					sphereSketch.Line(point1=(contactWidth, sphereCenterY - sqrt((radius+shellThickness)**2 - contactWidth**2)), 
									  point2=(contactWidth, sphereCenterY - (radius+shellThickness)))
					#partition for mesh
					if contactWidth < partitionWedgeThreshold:
						sphereSketch.Line (point1=(partitionWedgeWidth, sphereCenterY  - sqrt((radius+shellThickness)**2 - partitionWedgeWidth**2)),
										point2=(partitionWedgeWidth, sphereCenterY  - (radius+shellThickness)))
					#axissymm symmetry line
					sphereSketch.Line(point1=(0.0, sphereCenterY-(radius+shellThickness)), 
									  point2=(0.0, sphereCenterY+(radius+shellThickness)))
					# updating the start of next sphere chain, adding on matrix spacing
					bottomOfSystem += 2*(radius+shellThickness) + nrOfMatrixSpacings*matrixSpacingPerSphere;
					remainingSpheres -= 1; #keeping track of total sphere count
				else:
					print 'middel sphere'
					#midle sphere of chain system
					sphereSketch.ArcByCenterEnds(center=(0.0, sphereCenterY), 
												 direction=CLOCKWISE, 
												 point1=(contactWidth, sphereCenterY + sqrt((radius+shellThickness)**2 - contactWidth**2)), 
												 point2=(contactWidth, sphereCenterY - sqrt((radius+shellThickness)**2 - contactWidth**2)))
					#upper connection zone
					sphereSketch.Line(point1=(contactWidth, sphereCenterY + sqrt((radius+shellThickness)**2 - contactWidth**2)), 
									  point2=(contactWidth, sphereCenterY +(radius+shellThickness)))
					#lower connection zone
					sphereSketch.Line(point1=(contactWidth, sphereCenterY - sqrt((radius+shellThickness)**2 - contactWidth**2)), 
									  point2=(contactWidth, sphereCenterY - (radius+shellThickness)))
					#partition for mesh
					if contactWidth < partitionWedgeThreshold:
						sphereSketch.Line (point1=(partitionWedgeWidth, sphereCenterY  - sqrt((radius+shellThickness)**2 - partitionWedgeWidth**2)),
										point2=(partitionWedgeWidth, sphereCenterY  - (radius+shellThickness)))
						sphereSketch.Line (point1=(partitionWedgeWidth, sphereCenterY  + sqrt((radius+shellThickness)**2 - partitionWedgeWidth**2)),
										point2=(partitionWedgeWidth, sphereCenterY  + (radius+shellThickness)))
						wedgeCoordYVector.append(sphereCenterY + (radius + shellThickness));
					#axissymm symmetry line
					sphereSketch.Line(point1=(0.0, sphereCenterY-(radius+shellThickness)), 
									  point2=(0.0, sphereCenterY+(radius+shellThickness)))
					# updating the start of next sphere
					bottomOfSystem += 2*(radius+shellThickness);
					remainingSpheres -= 1; #keeping track of total sphere count
		print chainSystem;
	else:
		#========================================================================================================
		# 		Partition ordered spheres with shell
		#========================================================================================================
		randomSphereChainVector = np.ones(nrOfSpheres);
		for n in range(nrOfSpheres): 
			
			sphereCenterY = 0.0 + 2*cellHeight*n;
			sphereCenterYVector.append(sphereCenterY);
			
			#partition sphere
			sphereSketch.ArcByCenterEnds(center=(0.0 ,sphereCenterY), 
										 direction=COUNTERCLOCKWISE, 
										 point1=(0.0, sphereCenterY-radius), 
										 point2=(0.0, sphereCenterY+radius))

			#partition shell		
			sphereSketch.ArcByCenterEnds(center=(0.0, sphereCenterY), 
										 direction=COUNTERCLOCKWISE, 
										 point1=(0.0, sphereCenterY-(radius+shellThickness)), 
										 point2=(0.0, sphereCenterY+(radius+shellThickness)))
										 
			sphereSketch.Line(point1=(0.0, sphereCenterY-(radius+shellThickness)), 
							  point2=(0.0, sphereCenterY+(radius+shellThickness)))
			
	myPart.PartitionFaceBySketch(faces=myPart.faces, sketch=sphereSketch)
	#========================================================================================================
	# 		Partition infinite layer
	#========================================================================================================	
	infSketch = myModel.ConstrainedSketch(name='inf', sheetSize=50.0)
	infSketch.sketchOptions.setValues(viewStyle=AXISYM)
	infSketch.ConstructionLine(point1=(0.0,-25.0), point2=(0.0, 25.0))
	infSketch.Line(point1=(0.0, -cellHeight), point2=(width, -cellHeight))
	myPart.PartitionFaceBySketch(faces=myPart.faces, sketch=infSketch)
	
	#========================================================================================================
	# 		Sets and faces
	#========================================================================================================	 
	
	# Faces
	for nr in range(nrOfSpheres):
		if nr == 0:
			coreFaceMesh = myPart.faces.findAt(((radius/2,radius/2+sphereCenterYVector[nr],0.0),),);
			shellFaceMesh = myPart.faces.findAt(((0.0,radius + (shellThickness/2)+sphereCenterYVector[nr],0.0),),);
		else:
			coreFaceMesh += myPart.faces.findAt(((radius/2,radius/2+sphereCenterYVector[nr],0.0),),);
			shellFaceMesh += myPart.faces.findAt(((0.0,radius + (shellThickness/2)+sphereCenterYVector[nr],0.0),),);
	
	if areRandomSpheres and contactWidth < partitionWedgeThreshold:
		for nr in range(len(wedgeCoordYVector)):
			if nr == 0:
				wedgeFaceMesh = myPart.faces.findAt(((partitionWedgeThreshold+(partitionWedgeWidth-partitionWedgeThreshold)/2,wedgeCoordYVector[nr],0.0),),);
			else:
				wedgeFaceMesh += myPart.faces.findAt(((partitionWedgeThreshold+(partitionWedgeWidth-partitionWedgeThreshold)/2,wedgeCoordYVector[nr],0.0),),);

				
	matrixFace = myPart.faces.findAt((((width-radius)/2 + radius, 0.0, 0.0), ))	
	infFace = myPart.faces.findAt(((width/2, - (cellHeight + infLayerHeight/2), 0.0), ))
	
	
	# Sets
	
	coreSet = myPart.Set(faces=coreFaceMesh, name='Set-Sphere')	
	shellSet = myPart.Set(faces=shellFaceMesh, name='Set-Shell')
	infSet = myPart.Set(faces=infFace, name='Set-Inf')
	if areRandomSpheres and contactWidth < partitionWedgeThreshold:
		matrixSet = myPart.Set(faces=(matrixFace,wedgeFaceMesh,), name='Set-Matrix');
	else:
		matrixSet = myPart.Set(faces=matrixFace, name='Set-Matrix');
	
	#========================================================================================================
	# 		Assign sections
	#========================================================================================================	 	
	# Core
	myPart.SectionAssignment(offset=
		0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=coreSet, sectionName='Core', #'Core' 
		thicknessAssignment=FROM_SECTION)

	# Matrix
	myPart.SectionAssignment(offset=
		0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=matrixSet, sectionName='Matrix', 
		thicknessAssignment=FROM_SECTION)

	# Shell
	myPart.SectionAssignment(offset=
		0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=shellSet, sectionName='Shell', #'Shell'
		thicknessAssignment=FROM_SECTION)
		
	# InfElement
	myPart.SectionAssignment(offset=
		0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=infSet, sectionName='Matrix', 
		thicknessAssignment=FROM_SECTION)	


	#========================================================================================================
	# 		Assembly
	#========================================================================================================	 	
	myModel.rootAssembly.DatumCsysByThreePoints(coordSysType=CYLINDRICAL, 
												origin=(0.0, 0.0, 0.0), 
												point1=(1.0, 0.0, 0.0), 
												point2=(0.0, 0.0, -1.0))
	myInstance = myModel.rootAssembly.Instance(dependent=ON, name=instanceName, part=myPart)
	myModel.rootAssembly.regenerate()



	#========================================================================================================
	# 		Element sets for output data
	#========================================================================================================	 	
	bottomEdgeSet = myModel.rootAssembly.Set(edges=
		myInstance.edges.findAt(((width/3.0,-cellHeight,0.0),),), 
								name='bottomEdge');
		
	unitCellEdgeSet = myModel.rootAssembly.Set(edges=
		myInstance.edges.findAt(((width,-cellHeight/3.0,0.0),),), 
								name='unitCellEdge');
		
	loadEdgeSet = myModel.rootAssembly.Set(edges=
		myInstance.edges.findAt(((width/3.0,2*cellHeight*(nrOfSpheres-1)+cellHeight,0.0),),), 
								name='loadEdge');
		
	horisontalEdgesSet = myModel.rootAssembly.Set(edges=
		myInstance.edges.findAt(((width/3.0,-cellHeight,0.0),),
								((width/3.0,2*cellHeight*(nrOfSpheres-1)+cellHeight,0.0),),), 
								name='horisontalEdges')

	#========================================================================================================
	# 		Mesh
	#========================================================================================================
	# Shell
	
	myPart.setMeshControls(elemShape=QUAD, regions=shellFaceMesh, technique=FREE, algorithm=MEDIAL_AXIS)
	myPart.setElementType(
		elemTypes=(ElemType(elemCode=CAX4R, elemLibrary=EXPLICIT, 
		secondOrderAccuracy=OFF, hourglassControl=ENHANCED, 
		distortionControl=DEFAULT), ElemType(elemCode=CAX3, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, distortionControl=DEFAULT)), 
		regions=shellSet)

		
	# Infinite section 
	
	myPart.setMeshControls(regions=infFace, elemShape=QUAD, technique=SWEEP) 
	myPart.setElementType(
		elemTypes=(ElemType(elemCode=ACAX4R, elemLibrary=EXPLICIT), ElemType(elemCode=ACAX3, 
		elemLibrary=EXPLICIT)), 
		regions=infSet) # will be replaced by CINAX4 elements in inp file

	
	# Core and matrix
	
	myPart.setMeshControls(regions=coreFaceMesh, elemShape=TRI, technique=FREE)
	myPart.setMeshControls(regions=matrixFace, elemShape=TRI, technique=FREE)
	
	if areRandomSpheres and contactWidth < partitionWedgeThreshold:
		myPart.setMeshControls(regions=wedgeFaceMesh, elemShape=QUAD_DOMINATED, technique=FREE, algorithm=MEDIAL_AXIS)

	if quadraticElements == 1:
		print '-----------quadraticElements-----------'
		myPart.setElementType(
			elemTypes=(ElemType(elemCode=UNKNOWN_QUAD, elemLibrary=EXPLICIT), ElemType(
			elemCode=CAX6M, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, 
			distortionControl=DEFAULT)),
			regions=coreSet
			)
		myPart.setElementType(
			elemTypes=(ElemType(elemCode=UNKNOWN_QUAD, elemLibrary=EXPLICIT), ElemType(
			elemCode=CAX6M, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, 
			distortionControl=DEFAULT)), 
			regions=matrixSet
			) #CAX6M
	else:
		print '------------Linear Elements------------'
		myPart.setElementType(
			elemTypes=(ElemType(elemCode=CAX4R, elemLibrary=EXPLICIT), ElemType(elemCode=CAX3, 
			elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, distortionControl=DEFAULT)), 
			regions=coreSet
			)
		myPart.setElementType(
			elemTypes=(ElemType(elemCode=CAX4R, elemLibrary=EXPLICIT), ElemType(elemCode=CAX3, 
			elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, distortionControl=DEFAULT)), 
			regions=matrixSet
			) #CAX6M

	
	# seed part and generate mesh

	myPart.seedEdgeByNumber(constraint=FINER, edges=
		myPart.edges.findAt(((0.0, 
		-(cellHeight+infLayerHeight/2.0), 0.0), ), ((width, -(cellHeight+infLayerHeight/2.0), 0.0), ), ), number=1)
	myPart.seedPart(size=meshSize, deviationFactor=0.1, minSizeFactor=0.1) #constraint=FINER
	myPart.generateMesh()
	nrOfElements = len(myPart.elements);


	#========================================================================================================
	# 		Step
	#========================================================================================================
	myModel.ExplicitDynamicsStep(description=
		'Loading at .... frequency', linearBulkViscosity=linBulkViscosity, 
		maxIncrement=None, name=stepName, nlgeom=OFF, previous='Initial',  
		timePeriod=totalTimePeriod, quadBulkViscosity=quadrBulkViscosity, scaleFactor=1.0, 
		timeIncrementationMethod=AUTOMATIC_EBE) #EBE=element-by-element

		
	#--------------------------------------BC-----------------------------	
	myModel.DisplacementBC(amplitude=UNSET, createStepName='Initial', 
		distributionType=UNIFORM, fieldName='', localCsys=None, 
		name='UnitCellBC', region=unitCellEdgeSet, u1=SET, u2=UNSET, ur3=UNSET)
		
	#========================================================================================================
	# 		Pressure pulse
	#========================================================================================================
	fc = 5e6; 						# center frequency in [MHz]
	n = 4;							# number of half cycles included
	t0 = n*np.pi/(4*np.pi*fc);		# time shift of envelope curve
	t = np.arange(0,2*t0,1e-9); 	# time vector
	u = np.sin(2*np.pi*fc*t)*(np.cos(2*np.pi*fc*(t-t0)/n))**2; # amplitude
	amp = zip(t,u);
	amplitudeName = '5MHz'
	myModel.TabularAmplitude(data=amp,name=amplitudeName) #,timeSpan=TOTAL							
	
	myModel.rootAssembly.Surface(name='loadSurf', side1Edges=myInstance.edges.findAt(((width/3.0,2*cellHeight*(nrOfSpheres-1)+cellHeight,0.0), )));
	myModel.Pressure(amplitude=amplitudeName, createStepName='Loading', 
			distributionType=UNIFORM, field='', magnitude=1000.0, 
			name='pressure', 
			region=myModel.rootAssembly.surfaces['loadSurf'])
			
	#========================================================================================================
	# 		Output request
	#========================================================================================================
	myModel.HistoryOutputRequest(createStepName='Loading', name=
		'H-Output-2', rebar=EXCLUDE, region=horisontalEdgesSet, sectionPoints=
		DEFAULT, variables=stressOutput) #timeInterval=4e-09
		
		
		
	#========================================================================================================
	# 		Create and edit input file
	#========================================================================================================
	mdb.Job(activateLoadBalancing=False, atTime=None, contactPrint=OFF, 
		description='', echoPrint=OFF, explicitPrecision=SINGLE, 
		getMemoryFromAnalysis=True, historyPrint=OFF, memory=90, memoryUnits=
		PERCENTAGE, model=modelName, modelPrint=OFF, multiprocessingMode=DEFAULT, 
		name=jobName, nodalOutputPrecision=SINGLE, numCpus=1, numDomains=1, 
		parallelizationMethodExplicit=DOMAIN, queue=None, resultsFormat=ODB, 
		scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

	mdb.jobs[jobName].writeInput()
	writeFile = None
	filename = workDirectory + jobName + '.inp';
	readFile = open(filename, 'r')
	# Read in the file
	writeFile = readFile.read()
	# Replace the string
	writeFile = writeFile.replace('ACAX4R', 'CINAX4')
	# Write the file out again
	outName = 'infinite'+jobName;
	outFile = open(workDirectory + outName + '.inp', 'w')
	outFile.write(writeFile)
	readFile.close();
	outFile.close();
	print 'input file closed'
	del mdb.jobs[jobName]

	#========================================================================================================
	# 		Run job
	#========================================================================================================
	if runInput:
		myJob = mdb.JobFromInputFile(inputFileName= workDirectory + outName + '.inp', 
			name=odbName, type=ANALYSIS, numCpus=2, parallelizationMethodExplicit=DOMAIN,
			numDomains=2, multiprocessingMode=THREADS, explicitPrecision=SINGLE, nodalOutputPrecision=FULL,) #SINGLE/DOUBLE
		myJob.submit()
		myJob.waitForCompletion()
		
	return odbName, cell_rho, shellFraction, nrOfElements, randomSphereChainVector;
	
def odbWork(workDirectory,odbName,cell_rho):
	"""Operates on ODB data created by the function CREATEMODEL. 
		INPUT:
			workDirectory: Directory for ABAQUS files
			odbName: File name, used in "odbWork", includes geometry		
			cell_rho: Composite density. Used in development.
		OUTPUT: 
			topStressMatrix: Stress matrix of TOP elements ordered from inner to outer element 
			bottomStressMatrix: Stress matrix of BOTTOM elements ordered from inner to outer element
			distanceVector: Distance between integration points ordered from inner to outer
			timeVector: Time vector of elements ordered from inner to outer for the whole time periode.
			widthVector: Vector containing the width of every element form inner to outer.
	"""
	#==================================================================================
	# 		Operate on odb data
	#==================================================================================	
	myodbpath = workDirectory + odbName + '.odb';
	odb = openOdb(myodbpath);
	stressOutput = ('S33',);#('S33', 'MISES',); 	#'S33', 'MISES' etc.

	elementSetNameList = odb.rootAssembly.elementSets.keys()

	stepKey = odb.steps.keys() # vector of keys
	step = odb.steps[stepKey[0]] # using first key, as we only have one step
	regionString = step.historyRegions.keys()[1]; #using 1, as the [0] is for energy plots
	
	#==================================================================================
	# 		Sort bottom and top elements on their coordinates
	#==================================================================================	
	
	# Bottom Elements
	bottomElementList = {};
	bottomNodeList = {};
	for bottomElement in odb.rootAssembly.elementSets['BOTTOMEDGE'].elements[0]:
		if bottomElement.type != 'CINAX4':
			centerx = 0;
			centery = 0;
			centerz = 0;
			elementNodes = {};
			for node in bottomElement.connectivity:
				coord = odb.steps['Loading'].frames[0].fieldOutputs['S'].values[0].instance.getNodeFromLabel(node).coordinates
				elementNodes[node]=coord;
				centerx += coord[0]/3;
				centery += coord[1]/3;
				centerz += coord[2]/3;	
			center = (centerx,centery,centerz);
			bottomElementList[bottomElement.label] = center;
			#--------------------------------------------------------------------------
			#		Create a width vector of elements along bottom
			#--------------------------------------------------------------------------
			for node in elementNodes.keys():
				#coords = elementNodes[node];
				coords = (elementNodes[node][0],elementNodes[node][1], elementNodes[node][2]);
				if coords[0]<center[0] and coords[1]<center[1]:
					#found the first node 
					bottomNodeList[node] = coords 
				elif coords[0]>center[0] and coords[1]<center[1]:
					#found the last node
					bottomNodeList[node] = coords
	
	sortedBottomNodes = sorted(bottomNodeList, key=bottomNodeList.__getitem__)
	widthVector = np.zeros(len(sortedBottomNodes)-1);
	for i in range(len(widthVector)):
		widthVector[i] = bottomNodeList[sortedBottomNodes[i+1]][0]-bottomNodeList[sortedBottomNodes[i]][0];
			#---------------------------------------------------------------------------
	
	sortedBottomElements = sorted(bottomElementList, key=bottomElementList.__getitem__); # sorted by x-coordinates

	# Top Elements
	topElementList = {};
	for topElement in odb.rootAssembly.elementSets['LOADEDGE'].elements[0]:
		centerx = 0;
		centery = 0;
		centerz = 0;
		for node in topElement.connectivity:
			coord = odb.steps['Loading'].frames[0].fieldOutputs['S'].values[0].instance.getNodeFromLabel(node).coordinates
			centerx += coord[0]/3;
			centery += coord[1]/3;
			centerz += coord[2]/3;	
		center = (centerx,centery,centerz);
		topElementList[topElement.label] = center;				
	sortedTopElements = sorted(topElementList, key=topElementList.__getitem__) # sorted by x-coordinates
	
	
	#==================================================================================
	# 		Create a sorted stress matrix for top and bottom elemets
	#==================================================================================	
	
	# Top
	topStressMatrix = np.zeros((len(sortedTopElements),201));
	index = 0;
	for element in sortedTopElements:
		regionSubString = regionString[(regionString.find('-1.')+3):(regionString.find('Int Point')-1)];
		regionKey = regionString.replace(regionSubString,str(element))
		region = step.historyRegions[regionKey]
		for rk in region.historyOutputs.keys(): # taking out the data rk = 'MISES' or 'S33'
			if rk != stressOutput[0]: # and rk != stressOutput[1]: # Skip all plots of energy
				continue
			data = np.matrix(region.historyOutputs[rk].data); 
			topStressMatrix[index,:] = np.transpose(data[:,1]);
			index += 1;
	
	# Bottom
	bottomStressMatrix = np.zeros((len(sortedTopElements),201));
	index = 0;
	for element in sortedBottomElements:
		regionSubString = regionString[(regionString.find('-1.')+3):(regionString.find('Int Point')-1)];
		regionKey = regionString.replace(regionSubString,str(element))
		region = step.historyRegions[regionKey]
		for rk in region.historyOutputs.keys(): # taking out the data rk = 'MISES' or 'S33'
			if rk != stressOutput[0]: # and rk != stressOutput[1]: # Skip all plots of energy
				continue
			data = np.matrix(region.historyOutputs[rk].data); 
			bottomStressMatrix[index,:] = np.transpose(data[:,1]);
			index += 1;
	#==================================================================================
	# 		Calculate time shift, distance and velocity
	#==================================================================================	
	timeVector = data[:,0];	
	distanceVector = [];
	
	for elementsIndex in range(len(sortedBottomElements)):
		distance = (topElementList[sortedTopElements[elementsIndex]][1]
			- bottomElementList[sortedBottomElements[elementsIndex]][1]); #[1] gives y-coordinate
		distanceVector.append(distance);
	odb.close();
	
	#return np.mean(velocityVector)*10**-6, np.mean(velocityVector)*cell_rho*10**3, topStressMatrix, bottomStressMatrix, distanceVector, timeVector, widthVector;
	return topStressMatrix, bottomStressMatrix, distanceVector, timeVector, widthVector;

def writeOutputTofile(textOutDirectory, odbName, shellFraction, topStressMatrix, bottomStressMatrix, distanceVector, timeVector, nrOfElements, widthVector,randomSphereChainVector):		
	"""
	Writes output from "odbWork" to txt-file. stress, distance, time, width and randomChainVector is written out to seperate files.
	One infoFile is made to keep track of all the odbNames, shellfraction and number of elements used in the simulation.
	These txt-files are used for post processing in MATLAB.
	INPUT:
		textOutDirectory: Directory for the txt-files.
		odbName: File name, used in "odbWork", includes geometry.
		shellFraction: volume fraction of shell in the composite. 
		topStressMatrix: Stress matrix of TOP elements ordered from inner to outer element.
		bottomStressMatrix: Stress matrix of BOTTOM elements ordered from inner to outer element.
		distanceVector: Distance between integration points ordered from inner to outer.
		timeVector: Time vector of elements ordered from inner to outer for the whole time periode.
		nrOfElements: Total nr of elements used in the FE model. 
		widthVector: Vector containing the width of every element form inner to outer.
		randomSphereChainVector: Vector of particle chains	
	"""
	np.savetxt(textOutDirectory+'tStr__'+odbName+'.txt',np.transpose(topStressMatrix) , fmt='%.18e', delimiter=' ', newline='\n')
	np.savetxt(textOutDirectory+'bStr__'+odbName+'.txt',np.transpose(bottomStressMatrix) , fmt='%.18e', delimiter=' ', newline='\n')
	np.savetxt(textOutDirectory+'dist__'+odbName+'.txt', distanceVector, fmt='%.18e', delimiter=' ', newline='\n')
	np.savetxt(textOutDirectory+'time__'+odbName+'.txt',timeVector,fmt='%.18e', delimiter=' ', newline='\n')
	np.savetxt(textOutDirectory+'width_'+odbName+'.txt', widthVector,fmt='%f', delimiter=' ', newline='\n')
	np.savetxt(textOutDirectory+'randomChainVector_'+odbName+'.txt', randomSphereChainVector,fmt='%f', delimiter=' ', newline='\n')
	
	infoFileName = textOutDirectory + 'infoFile.txt';
	with open(infoFileName,'a') as infoFile:
			infoFile.write("%s,%f,%f\n" % (odbName, shellFraction, nrOfElements));