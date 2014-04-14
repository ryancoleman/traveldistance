#!/usr/bin/env python

#Ryan G. Coleman, Kim A. Sharp
#ryan.g.coleman@gmail.com http://crystal.med.upenn.edu

#this is a program to replace the excellent trigen (fortran) with python
#reads in a .tri file, outputs a .tst file

#2 kinds of input, first the pdb_record, then the triangles:
"""
PDB_RECORD
ATOM      1  N   HIS     4      10.859  -0.424   8.156  1.65 28.87
ATOM      2  CA  HIS     4      10.020   0.343   9.069  1.90 28.74
...
END PDB_RECORD
TRIANGLE_XYZ
   -30.627  -11.711   16.376
   -31.291  -11.711   17.039
   -30.627  -12.374   17.039
   -30.627  -12.374   17.039
   -31.291  -11.711   17.039
   -31.291  -11.711   18.020
   -30.627  -12.374   17.039
...
"""
#note that each 'triple' is an xyz coord and the first 3 lines are a triangle,
#as are the next 3, etc, and the points are repeated as necessary

import string, sys
import pdb, tstdata
import geometry

def readTri(triName):
  '''reads a tri file, returns the data inside'''
  pdbD = pdb.pdbData()
  triLineToNum = {}
  triPoints, oneTri = [], []
  triFile = open(triName, 'r')
  try:
    inPdb, inTri = False, False
    for line in triFile:
      if line.startswith('PDB_RECORD'):
        inPdb = True
      elif line.startswith('END PDB_RECORD'):
        inPdb = False
      elif line.startswith('TRIANGLE_XYZ'):
        inTri = True
      elif line.startswith('END TRIANGLE_XYZ'):
        inTri = False
      else: #actually do some processing of lines
        if inPdb:
          pdbD.processLine(line)
        elif inTri:
          strippedLine = string.strip(line)
          if strippedLine not in triLineToNum:
            triLineToNum[strippedLine] = len(triLineToNum) + 1
          triNum = triLineToNum[strippedLine]
          oneTri.append(triNum)
          if len(oneTri) == 3:
            triPoints.append(oneTri)
            oneTri = []
  except StopIteration: #EOF
    pass
  return pdbD, triLineToNum, triPoints

def makePointXyz(triLineToNum):
  '''takes a dictionary of pointString -> pointNum and makes a list of
  pointNum, pointX, pointY, pointZ'''
  #make the reverse dict
  reverseDict = {}
  for key, value in triLineToNum.items():
    reverseDict[value] = key
  pointNums = reverseDict.keys()
  pointNums.sort()
  outputList = []
  for pointNum in pointNums:
    pointStrings = string.split(reverseDict[pointNum])
    newInnerList = [pointNum]
    for pointStr in pointStrings:
      newInnerList.append(float(pointStr))
    outputList.append(newInnerList)
  return outputList

def makeTriPoints(triPoints):
  '''assigns numbers to the triangles, returns list'''
  outputList = []
  triNum = 0
  for triPointList in triPoints:
    triNum += 1
    triPointList.insert(0, triNum)
    outputList.append(triPointList)
  return outputList

def makeFakeTriNeighbors(trianglePoints):
  '''makes fake neighbors, like current trigen, TODO update to do it for real'''
  outputList = []
  for triPointList in trianglePoints:
    outputList.append([triPointList[0], 0,0,0])
  return outputList

def nearbyPoints(pointXyz, pdbD):
  '''finds the nearest atom to each surface point, makes list'''
  outputList = []
  for pointCoord in pointXyz:
    pointNum = pointCoord[0]
    xyz = pointCoord[1:4]
    closestIndex, closestDist = 0, 1000000000000.
    for index, atom in enumerate(pdbD.coords):
      thisDist = geometry.dist(xyz, atom, metric="L2SQUARED") #l2^2 monotonic
      if thisDist < closestDist:
        closestDist = thisDist
        closestIndex = index
    outputList.append([pointNum, closestIndex+1])
  return outputList

def nearbyTriangles(pointPdbRecord, trianglePoint):
  '''finds the consensus of nearest atoms to each triangle'''
  outputList = []
  for trianglePointRecord in trianglePoint:
    triNum = trianglePointRecord[0]
    nearbys = []
    for point in trianglePointRecord[1:]: #each point in this tri
      nearbys.append(pointPdbRecord[point-1][1])
    nearbys.sort()
    if nearbys[1] == nearbys[2]: #if these 2 match, use one of them
      outputList.append([triNum, nearbys[1]])
    else: #otherwise no matches or the first one matches, just use first one
      outputList.append([triNum, nearbys[0]])
  return outputList

def makeDicts(triPoints):
  '''finds each tri that each point is in'''
  pointTriDict = {}
  triPointDict = {}
  for triPointRec in triPoints:
    triNum = triPointRec[0]
    triPointDict[triNum] = triPointRec[1:] #kept clockwise
    for point in triPointRec[1:]:
      if point not in pointTriDict:
        pointTriDict[point] = [] #new empty list
      pointTriDict[point].append(triNum)
  return pointTriDict, triPointDict

def getNextPoint(curTri, point):
  '''helper that gets the next clockwise point from the triangle'''
  position = curTri.index(point)
  if position == 2:
    return curTri[0]
  else:
    return curTri[position+1]

def makePointTriAndNeighbor(pointTriDict, triPointDict):
  '''makes the list of point triangles and neighbors in clockwise direction'''
  outputPtTri, outputPtNeighbor = [],[]
  points = pointTriDict.keys()
  points.sort()
  for point in points:
    tris = pointTriDict[point][:] #must copy[:] since going to get destroyed
    outputPtTri.append([point, len(tris)]) #know num,but not order
    outputPtNeighbor.append([point, len(tris)])
    curTri = tris[0] #could pick any to start
    curPoint = getNextPoint(triPointDict[curTri], point)
    orderedTris, orderedPoints = [curTri], [curPoint] #start lists
    tris.remove(curTri)
    while len(tris) > 0:
      nextPt = getNextPoint(triPointDict[curTri], curPoint)
      possNextTris = pointTriDict[nextPt]
      nextTri = False
      for possNextTri in possNextTris:
        if possNextTri in tris:
          nextTri = possNextTri
          break
      tris.remove(nextTri) #destructively change this copy
      orderedTris.append(nextTri)
      orderedPoints.append(nextPt)
      curPoint, curTri = nextPt, nextTri
    outputPtTri[-1].extend(orderedTris)
    outputPtNeighbor[-1].extend(orderedPoints)
  return outputPtTri, outputPtNeighbor

def calcNormals(pointXyz, triPointDict, pointTriDict):
  '''calculate normals for each point by averaging the adjacent tri normals'''
  triNormals = {}
  for tri in triPointDict.keys():
    points = triPointDict[tri]
    coords = []
    for point in points:
      coords.append(pointXyz[point-1][1:])
    triNormals[tri] = geometry.getTriNormalList(coords)
  outputList = []
  for pointXyzRecord in pointXyz:
    pointNum = pointXyzRecord[0]
    normals = []
    for tri in pointTriDict[pointNum]:
      normals.append(triNormals[tri])
    thisNormal = geometry.getAverage(normals) #just average the normals
    outputList.append([pointNum, thisNormal[0], thisNormal[1], thisNormal[2]])
  return outputList

def fakeForNow(pointXyz):
  '''just list of xyz with 0s, fake it for now, TODO make real'''
  outputList = []
  for pointRec in pointXyz:
    outputList.append([pointRec[0], 0.])
  return outputList

def triToTst(triName, tstName):
  '''reads a tri file, writes the tstfile'''
  pdbD, triLineToNum, triPoints = readTri(triName)
  #print pdbD, triLineToNum, triPoints
  tstD = tstdata.tstDataWritable(False) #empty, now add stuff to it
  tstD.dict['PDB_RECORD'] = pdbD.rawData #this one is easy
  tstD.dict['POINT_XYZ'] = makePointXyz(triLineToNum)
  tstD.dict['TRIANGLE_POINT'] = makeTriPoints(triPoints)
  tstD.dict['TRIANGLE_NEIGHBOR'] = \
                        makeFakeTriNeighbors(tstD.dict['TRIANGLE_POINT'])
  tstD.dict['POINT_PDB_RECORD'] = nearbyPoints(tstD.dict['POINT_XYZ'], pdbD)
  tstD.dict['TRIANGLE_PDB_RECORD'] = \
                        nearbyTriangles(tstD.dict['POINT_PDB_RECORD'], \
                              tstD.dict['TRIANGLE_POINT'])
  pointTriDict, triPointDict = makeDicts(tstD.dict['TRIANGLE_POINT'])
  tstD.dict['POINT_TRIANGLE'],tstD.dict['POINT_NEIGHBOR'] = \
                             makePointTriAndNeighbor(pointTriDict, triPointDict)
  tstD.dict['NORM_XYZ'] = calcNormals(tstD.dict['POINT_XYZ'], \
                             triPointDict, pointTriDict)
  tstD.dict['CURVATURE_XYZ'] = fakeForNow(tstD.dict['POINT_XYZ'])
  tstD.dict['PROPERTY_XYZ']= fakeForNow(tstD.dict['POINT_XYZ'])
  #print tstD.dict.keys()
  tstD.write(tstName)

if __name__ == "__main__" and -1 != string.find(sys.argv[0], "trigenTst.py"):
  #otherwise this script has been imported, no reason to run
  if len(sys.argv) > 1:
    triName = sys.argv[1]
    if len(sys.argv) > 2:
      tstName = sys.argv[2]
    else: #assume how to make name
      tstName = triName[:-4] + ".tst"
    triToTst(triName, tstName)
  else:
    print "usage: trigenTst.py name.tri [name.tst]"
