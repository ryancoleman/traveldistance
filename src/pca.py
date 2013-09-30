#!/usr/bin/env python2.5

#use numeric utility to compute principal components of points

try: #to use numeric (old)
  from Matrix import Matrix
  from LinearAlgebra import eigenvectors
except ImportError: #use numpy (new)
  import numpy.oldnumeric as Numeric
  from numpy import matrix as Matrix
  from numpy.linalg import eig as eigenvectors
except ImportError:
  print "you do not have numpy or numeric installed under this python version"
  exit(1)
import geometry #for getAverage of points, and dot product
import operator #for sorting tricks

def pca2d(pointList):
  '''sets up the pca for a list of points in 2d. solves eig problem'''
  matrixList = [[0,0],[0,0]] #init to 0
  avgPt = geometry.getAverageArbitraryDimension(pointList, 2)
  diffs = [0,0]
  for point in pointList:
    for index in range(2):
      diffs[index] = point[index] - avgPt[index]
    #matrix is old plus a2 ab
    #                   ba b2 
    #          only compute upper diagonal now
    matrixList[0][0] += diffs[0]*diffs[0] #just hardcode it
    matrixList[0][1] += diffs[0]*diffs[1]
    matrixList[1][1] += diffs[1]*diffs[1]
  #make symmetric
  matrixList[1][0] = matrixList[0][1]
  actualMatrix = Matrix(matrixList)
  val,vec = eigenvectors(actualMatrix)
  return val, vec

def pca3d(pointList):
  '''sets up the pca for a list of points in 3d. solves eig problem'''
  matrixList = [[0,0,0],[0,0,0],[0,0,0]] #init to 0
  avgPt = geometry.getAverage(pointList)
  diffs = [0,0,0]
  for point in pointList:
    for index in range(3):
      diffs[index] = point[index] - avgPt[index]
    #matrix is old plus a2 ab ac
    #                   ba b2 bc
    #                   ca cb c2 only compute upper diagonal now
    matrixList[0][0] += diffs[0]*diffs[0] #just hardcode it
    matrixList[0][1] += diffs[0]*diffs[1]
    matrixList[0][2] += diffs[0]*diffs[2]
    matrixList[1][1] += diffs[1]*diffs[1]
    matrixList[1][2] += diffs[1]*diffs[2]
    matrixList[2][2] += diffs[2]*diffs[2]
  #make symmetric
  matrixList[1][0] = matrixList[0][1]
  matrixList[2][0] = matrixList[0][2]
  matrixList[2][1] = matrixList[1][2]
  actualMatrix = Matrix(matrixList)
  val,vec = eigenvectors(actualMatrix)
  return val, vec

def findLongestDirection(pointList):
  '''does pca, gets the eigenresult, find the longest direction, returns it'''
  eigenvalues, eigenvectors = pca3d(pointList)
  maxVal, maxIndex = 0,0
  try:
    for index in range(len(eigenvalues)):
      if maxVal < eigenvalues[index]:
        maxVal = eigenvalues[index]
        maxIndex = index
    maxVec = eigenvectors[maxIndex]
  except TypeError: #caused by complex numbers
    for index in range(len(eigenvalues)):
      if maxVal < eigenvalues[index].real: #real prevents imaginary problems
        maxVal = eigenvalues[index].real
        maxIndex = index
    maxVec = eigenvectors[maxIndex].real
  try:
    maxVecRet = maxVec.tolist() #stupid numpy new2010
  except AttributeError:
    maxVecRet = maxVec
  return maxVecRet
    
def sortDirections(pointList):
  '''finds the eigenvalues and eigenvectors, sorts based on eigenvalue and 
  returns list of direction vectorsin descending order'''
  eigenvalues, eigenvectors = pca3d(pointList)
  maxVal, maxIndex = 0,0
  #print eigenvalues, eigenvectors
  try:
    newEigList = []
    for index in xrange(len(eigenvalues)):
      newEigList.append((eigenvalues[index], eigenvectors[index]))
    newEigList.sort(key=operator.itemgetter(0))
    newEigList.reverse()
  except TypeError: #imaginary problem
    newEigList = []
    for index in xrange(len(eigenvalues)):
      newEigList.append((eigenvalues[index].real, eigenvectors[index].real))
    newEigList.sort(key=operator.itemgetter(0))
    newEigList.reverse()
  retEigVecList = []
  for eigenvalue, eigenvector in newEigList:
    try:
      newEig = eigenvector.tolist() #stupid numpy new 2010 [0] before?
      retEigVecList.append(newEig)
    except AttributeError: #numpy is actually not being run, not a problem
      retEigVecList.append(eigenvector)
  return retEigVecList

def findLongestDimension(pointList):
  '''calls direction, returns length in that direction'''
  direction = findLongestDirection(pointList)
  #direction is a unit vector, so can do scalar projection, i.e.
  #dot product of each point with direction gives the length in that direction
  #and is negative if in the opposite direction, so just find max-min and return
  firstPoint = pointList[0]
  try:
    if 1 == len(direction):
      direction = direction[0] #numpy/numeric return differently
  except TypeError:
    pass
  min = geometry.dot(firstPoint, direction)
  max = min #same for now
  for point in pointList[1:]: #already done first point
    newScalar = geometry.dot(point, direction)
    try:
      if newScalar < min:
        min = newScalar
      if newScalar > max:
        max = newScalar
    except TypeError:
      newScalar = newScalar.real
      min = min.real
      max = max.real
      if newScalar < min:
        min = newScalar
      if newScalar > max:
        max = newScalar
  return max-min   

def findDimensions(pointList):
  '''finds the dimension in the 3 principal directions, longest first'''
  dimensions = []
  directions = sortDirections(pointList)
  for direction in directions:
    firstPoint = pointList[0]
    try:
      if 1 == len(direction):
        direction = direction[0] #stupid numpy
    except TypeError:
      pass
    min = geometry.dot(firstPoint, direction)
    max = min #same for now
    for point in pointList[1:]: #already done first point
      newScalar = geometry.dot(point, direction)
      try:
        if newScalar < min:
          min = newScalar
        if newScalar > max:
          max = newScalar
      except TypeError:
        newScalar = newScalar.real
        min = min.real
        max = max.real
        if newScalar < min:
          min = newScalar
        if newScalar > max:
          max = newScalar
    dimensions.append(max - min)
  return dimensions 

#disable testing now
'''
import sys,string
if -1 != string.find(sys.argv[0], "pca.py"):
  print findLongestDimension([[1,1,1],[0,0,0],[-1,-1,-1]]) #stupid test
  print findLongestDimension([[10,10,10],[9,9,9],[8,8,8]]) #stupid test
  print findLongestDimension([[10,9,1],[9,9,9],[8,7,8]]) #stupid test
  print findDimensions([[10,9,1],[9,9,9],[8,7,8],[1,1,1]]) #stupid test
'''
