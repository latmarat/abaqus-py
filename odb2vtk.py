from abaqus import *
from abaqusConstants import *
from odbAccess import *
import __main__

import numpy as np
import sys

def mequit(f):
	f.close()
	sys.exit('Script exited with an error')

def writeVtkMesh(vtkFileName,myInstance,numNodes,numElements,myFrame):

	vtkFile = open(vtkFileName, 'w')
	vtkFile.write('# vtk DataFile Version 2.0\n')
	vtkFile.write('Reconstructed Lagrangian Field Data\n')
	vtkFile.write('ASCII\n\n')
	vtkFile.write('DATASET UNSTRUCTURED_GRID\n')

	# Get the initial nodal coordinates
	initialCoords = np.ndarray(shape=(3,numNodes))
	for nd in range(0, numNodes):
		coords = myInstance.nodes[nd].coordinates
		for iCoor, coord in enumerate(coords):
			initialCoords[iCoor,nd] = coord

	elementConnectivity = []
	for el in range(0, numElements):
		con = myInstance.elements[el].connectivity
		elementConnectivity.append(con)

	totalNumNodes = numNodes
	totalNumElements = numElements

	# Print out all the coordinates to the vtk file
	vtkFile.write('POINTS %i float\n' % (totalNumNodes))

	r = np.ndarray(shape=(3,numNodes))

	# Isolate the displacement field
	displacements = myFrame.fieldOutputs['U'].getSubset(region=myInstance).values

	# Add displacements to the initial coordinates
	for nd in range(0, numNodes):
		for iCoor in range(0,3):
			r[iCoor,nd] = initialCoords[iCoor,nd] + displacements[nd].data[iCoor]
			vtkFile.write('%f\t' % (r[iCoor,nd]))
		vtkFile.write('\n')

	# Print out all the elements to the vtk file
	vtkFile.write('\nCELLS %i %i\n' % (totalNumElements, 9*totalNumElements))
	for el in range(0,numElements):
		numNodesEl = len(elementConnectivity[el])
		n = [0]*numNodesEl
		vtkFile.write('%i\t' % numNodesEl)
		for iNode, node in enumerate(elementConnectivity[el]):
			n[iNode] = node - 1
			vtkFile.write('%i\t' % n[iNode])
		vtkFile.write('\n')

	# Print out all the cell types to the vtk file
	vtkFile.write('\nCELL_TYPES %i\n' % (totalNumElements))
	for el in range(0,totalNumElements):
		vtkFile.write('12\n')

	vtkFile.write('\nCELL_DATA %i' % (totalNumElements))

	vtkFile.close()

def getScalarData(varName, myFrame, varField, myInstance, numElements):

	# vtkFile.write('\nTENSORS %s float\n' % varName)
	aveScalar = np.zeros((numElements))

	numIntPts = len(varField.getSubset(region=myInstance.elements[0]).values)
	scalarData = np.zeros((numIntPts))
	for el in range(0,numElements):
		# Isolate current and previous element's stress field
		abaVar = varField.getSubset(region=myInstance.elements[el]).values

		for ip,ipValue in enumerate(abaVar):
			scalarData[ip] = ipValue.data

		# Average the variable values for the element
		aveScalar[el] = np.average(scalarData,axis=0)

	return aveScalar

def getSymTensorData(varName, myFrame, varField, myInstance, numElements):

	# vtkFile.write('\nTENSORS %s float\n' % varName)


	numIntPts = len(varField.getSubset(region=myInstance.elements[0]).values)
	numComps = len(varField.getSubset(region=myInstance.elements[0]).values[0].data)
	tensorData = np.zeros((numIntPts,numComps))
	aveTensor = np.zeros((numElements,numComps))

	for el in range(0,numElements):
		# Isolate current and previous element's stress field
		abaVar = varField.getSubset(region=myInstance.elements[el]).values

		for ip,ipValue in enumerate(abaVar):
			for icomp,component in enumerate(varAtPt.data)
				tensorData[ip,icomp] = ipValue.data[icomp]

		# Average the variable values for the element
		aveTensor[el,:] = np.average(tensorData,axis=0)
	return aveTensor


def writeTensorData(vtkFileName, tensorData, tensorName):

	with open(vtkFileName,'a') as vtkFile:
		vtkFile.write('\nTENSORS %s float\n' % tensorName)
		np.savetxt(vtkFile, tensorData, fmt='%.4f', delimiter=' ')

fileName = 'sdv-test'
# varNames = ['SDV14', 'SDV15', 'SDV16', 'SDV17', 'SDV18', 'SDV19', 'SDV20', 'SDV21', 'SDV22']
varNames = ['S','PEEQ']
instanceName = None

f = open('odb2vtk.out','w')

# Open the odb and request re-input if the odb is not found
try:
	odbPath = fileName + '.odb'
	myOdb = session.openOdb(name=odbPath)
except OdbError:
	f.write('ERROR! Correct path to odb was not provided\n')
	mequit(f)

msg = 'odb is successfully opened at %s\n' % odbPath
f.write(msg)

# The first instance is default if not provided
if instanceName == None:
	instanceName = myOdb.rootAssembly.instances.keys()[0]

# Name the Step object for convenience
myStepName = myOdb.steps.keys()[0]
myStep = myOdb.steps[myStepName]
msg = 'Working with step %s and instance %s\n' % (myStepName, instanceName)
f.write(msg)

# Check if the given instance exists and request re-input if it doesn't
instanceNames = myOdb.rootAssembly.instances.keys()
if instanceName not in instanceNames:
	instanceList = '\n'.join(iInstance for iInstance in instanceNames)
	msg = '%s was not found! \nRe-input the instance name. \nInstances in the odb:\n%s'  % (instanceName, instanceList)
	f.write(msg)
	mequit(f)

myInstance = myOdb.rootAssembly.instances[instanceName]
numElements = len(myInstance.elements)
numNodes = len(myInstance.nodes)
numFrames = len(myStep.frames)
time = [0]

f.write('Number of frames: %s\n' % str(numFrames))

mySet = myOdb.rootAssembly.instances[instanceName].elementSets['ALLELEMENTS']

velGrdData = np.zeros((numElements,9))

for iframe in range(1,numFrames):
	vtkFileName = '%s-%03d.vtk' % (fileName, iframe)

	myFrame = myStep.frames[iframe]
	time[0] = myFrame.frameValue

	f.write('%7.4f\n' % time[0])
	f.write(vtkFileName + '\n')

	writeVtkMesh(vtkFileName,myInstance,numNodes,numElements,myFrame)

	for i,ivar in enumerate(varNames):
		varField = myFrame.fieldOutputs[ivar]
		field = varField.getSubset(region=mySet, position=CENTROID)

		print str(field.type)

		if str(field.type) == 'SCALAR':
		 	scalar = getScalarData(ivar, myFrame, varField, myInstance, numElements)
			writeScalarData(vtkFileName, velGrdData, ivar)
		if str(field.type) == '3D_TENSOR':
			scalar = getScalarData(ivar, myFrame, varField, myInstance, numElements)
			writeTensorData(vtkFileName, velGrdData, ivar)

f.close()