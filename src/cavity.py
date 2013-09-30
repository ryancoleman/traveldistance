#!/usr/bin/env python2.5

#a collection of functions for dealing with (mostly removing) cavities
#Ryan G Coleman, Kim A Sharp, crystal.med.upennn.edu, ryan.g.coleman@gmail.com

from unionfind2 import unionFind
import phi
import grid
import tstdata
import geometry
import orstHelper
import string, sys #so can execute by itself

def findBiggestDisjointSets(pointList, triList, pointNeighborList):
  '''slightly improved code-- well 15 seconds faster on small stuff'''
  pointSetUF = unionFind()
  for nhbrPointsList in pointNeighborList:
    #first check to see if point is already in a list
    startPt = nhbrPointsList[0]
    for otherPt in nhbrPointsList[2:]:
      pointSetUF.union(startPt, otherPt)
  pointSets = pointSetUF.toLists()
  #remove points + tris not in the biggest disjoint set (cavities)
  largest, size = 0,0
  for index in xrange(len(pointSets)):
    if len(pointSets[index]) > size:
      largest, size = index,len(pointSets[index])
  allowedPoints = pointSets[largest]
  #figured it out, make sets
  allPoints,cavPoints = set(), set()
  for point in pointList:
    if int(point[0]) in allowedPoints:
      allPoints.update([int(point[0])])
    else:
      cavPoints.update([int(point[0])])
  allTris, cavTris = set(), set()
  for tri in triList:
    if int(tri[1]) in allPoints: #any triangle point is okay
      allTris.update([int(tri[0])])
    else:
      cavTris.update([int(tri[0])])
  #print len(allPoints), len(pointList)
  #print len(allTris), len(triList)
  return allPoints, allTris, cavPoints, cavTris

def findBiggestDisjointSetsBreakCavities(pointList, triList, pointNeighborList):
  '''breaks out each cavity separately. doesn't return tris, just points'''
  pointSetUF = unionFind()
  for nhbrPointsList in pointNeighborList:
    #first check to see if point is already in a list
    startPt = nhbrPointsList[0]
    for otherPt in nhbrPointsList[2:]:
      pointSetUF.union(startPt, otherPt)
  pointSets = pointSetUF.toLists()
  #remove points + tris not in the biggest disjoint set (cavities)
  largest, size = 0,0
  for index in xrange(len(pointSets)):
    if len(pointSets[index]) > size:
      largest, size = index,len(pointSets[index])
  allowedPoints = pointSets[largest]
  #figured it out, make sets
  allPoints,cavPoints = set(allowedPoints), set()
  pointSets.remove(allowedPoints)
  for cavPtSet in pointSets:
    cavPoints.update(cavPtSet)
  #print len(allPoints), len(pointList)
  return allPoints, cavPoints, pointSets

def assumeNoCavities(pointList, triList, pointNeighborList=False):
  '''just like findBiggestDisjoinUnions except assume only one set'''
  allPoints, allTris = [], []
  for point in pointList:
    allPoints.append(int(point[0]))
  for tri in triList:
    allTris.append(int(tri[0]))
  return allPoints,allTris,[],[]

def makeMapLists(fullList, shortList):
  '''helper function that makes a map from the full to the short lists... 1 indexing'''
  mapPoints = {}
  allPointsList = list(shortList)
  allPointsList.sort()
  for point in fullList:
    newIndex = -2
    try:
      newIndex = allPointsList.index(int(point[0]))
    except ValueError:
      pass
    mapPoints[int(point[0])] = newIndex+1 #1-indexing
  return mapPoints

def replaceEntry(data, primaryMap, otherMaps=False, extendToEnd=False):
  '''needs a map for each 'column' in data'''
  newData = []
  for row in data:
    newRow  = []  
    if False != primaryMap:
      try:
        newRow.append(primaryMap[row[0]])
      except KeyError:
        newRow.append(row[0])
    else:
      newRow.append(row[0])
    for index,eachMap in enumerate(otherMaps):
      if False != eachMap:
        try:
          newRow.append(eachMap[row[index+1]])
        except KeyError:
          newRow.append(row[index+1])
      else:
        newRow.append(row[index+1])
    if extendToEnd:
      #lastMap = otherMaps[len(otherMaps)-1] #not used??? who knows, code is old
      for index in xrange(len(otherMaps)+1,len(row)):
        if False != eachMap:
          try:
            newRow.append(eachMap[row[index]])
          except KeyError:
            newRow.append(row[index])
        else:
          newRow.append(row[index])
    if newRow[0] != -1: #this means removed from list
      newData.append(newRow)
  return newData

def tstCavityRemoval(tstFileName,tstFileNameOut, \
                     phiFileName=False,phiFileNameOut=False):
  oldTst = tstdata.tstDataWritable(tstFileName) #read the file into the struct
  allPoints, allTris, cavPoints, cavTris = findBiggestDisjointSets( \
      oldTst.dict['POINT_XYZ'], oldTst.dict['TRIANGLE_POINT'], \
      oldTst.dict['POINT_NEIGHBOR']) 
  #print len(allPoints), len(allTris), len(cavPoints), len(cavTris) #debug
  #now remove everything that isn't allPoints and allTris
  #first generate maps from old numbering to new numbering
  mapPoints = makeMapLists(oldTst.dict['POINT_XYZ'], allPoints)
  mapTris = makeMapLists(oldTst.dict['TRIANGLE_POINT'], allTris)
  #ugh. phi map.
  phiData = phi.phi(phiFileName)  #read in the phimap if possible
  gridD = grid.makeGridFromPhi(phiData,False) #the False just makes it copy values
  mins, maxs = phiData.getMinsMaxs()
  maxVal = phiData.getMaxValues()
  gridSize = 1.0/phiData.scale
  #this function marks cavities as interior
  cavTriTuples = geometry.cacheTriangle( \
            oldTst.dict['TRIANGLE_POINT'], oldTst.dict['POINT_XYZ'], cavTris)
  #have to do this now before oldTst gets modified
  orstHelper.decideInside(gridD, cavTriTuples, cavPoints, \
                    oldTst.dict['POINT_TRIANGLE'], oldTst.dict['POINT_XYZ'], \
                 oldTst.dict['TRIANGLE_POINT'], maxVal)  #change anything inside cavities to max..
  #now piece-by-piece fix each entry
  newTN = replaceEntry(oldTst.dict['TRIANGLE_NEIGHBOR'], mapTris, [mapTris,mapTris,mapTris])
  oldTst.dict['TRIANGLE_NEIGHBOR'] = newTN
  newPX = replaceEntry(oldTst.dict['POINT_XYZ'], mapPoints, [False, False, False])
  oldTst.dict['POINT_XYZ'] = newPX
  newTP = replaceEntry(oldTst.dict['TRIANGLE_POINT'], mapTris, [mapPoints, mapPoints, mapPoints])
  oldTst.dict['TRIANGLE_POINT'] = newTP
  newPT = replaceEntry(oldTst.dict['POINT_TRIANGLE'], mapPoints, [False, mapTris], extendToEnd=True)
  oldTst.dict['POINT_TRIANGLE'] = newPT
  newPN = replaceEntry(oldTst.dict['POINT_NEIGHBOR'], mapPoints, [False, mapPoints], extendToEnd=True)
  oldTst.dict['POINT_NEIGHBOR'] = newPN
  newNX = replaceEntry(oldTst.dict['NORM_XYZ'], mapPoints, [False,False,False])
  oldTst.dict['NORM_XYZ'] = newNX
  newPPR = replaceEntry(oldTst.dict['POINT_PDB_RECORD'], mapPoints, [False])
  oldTst.dict['POINT_PDB_RECORD'] = newPPR
  newTPR = replaceEntry(oldTst.dict['TRIANGLE_PDB_RECORD'], mapTris, [False])
  oldTst.dict['TRIANGLE_PDB_RECORD'] = newTPR
  #don't touch PDB_RECORD...
  newCX = replaceEntry(oldTst.dict['CURVATURE_XYZ'], mapPoints, [False])
  oldTst.dict['CURVATURE_XYZ'] = newCX
  newPX = replaceEntry(oldTst.dict['PROPERTY_XYZ'], mapPoints, [False])
  oldTst.dict['PROPERTY_XYZ'] = newPX
  #now write it back out...should probably tag it somehow....hmmm
  oldTst.write(tstFileNameOut)
  #write phi map
  phiDataOut = phi.phi()
  phiDataOut.createFromGrid(gridD,gridSize,toplabel=phiData.toplabel, \
                head=phiData.head,title=phiData.title,botlabel=phiData.botlabel)
  phiDataOut.write(phiFileNameOut)

#this is where main is... maybe add some other arguments like gridSize?
if -1 != string.find(sys.argv[0], "cavity.py"):
  if 4 < len(sys.argv):
    tstFile, tstOut, phiFile, phiOut = sys.argv[1:5]
    print tstFile, tstOut, phiFile, phiOut
    tstCavityRemoval(tstFile, tstOut, phiFile, phiOut)
  else:
    print "Usage: cavity.py tstFile tstFileOut phiFile phiFileOut"
