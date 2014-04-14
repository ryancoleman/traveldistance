#!/usr/bin/env python

import math
import string
import sys
import os  # for running commandline instead of through pymol
from os.path import isfile
#following custom module must be in same directory or in usual path
import tstdata
from sharp_phi import phi  # phi map reader/writer
import tstdebug  # debugging functions, moved out of codebase
import geometry   # geometric primitives
from pdb import pdbData
import grid   # moved primitive grid functions
import cavity
import orstHelper
import travelDistNoMesh

def tstTravelDepthNoPhi(tstFileName, gridSize=1.0):
  '''does old style (v1.0) travel depth with no phi-map, finds own inside and
  outside data, proceeds as normal'''
  tstD = tstdata.tstData(tstFileName)  # read the file into the data structure
  if 'CONVEX_HULL_TRI_POINT_LIST' not in tstD.dict.keys():
    print "Run tstConvexHull.py on this tst data file first."
    sys.exit(1)
  #these sets are useful to construct
  convexHullPoints = set()
  for record in tstD.dict['CONVEX_HULL_TRI_POINT_LIST']:
    convexHullPoints.update(record[1:])
  #do the biggest disjoint set of tris/points stuff
  allPoints, allTris, cavPoints, cavTris = cavity.assumeNoCavities(
      tstD.dict['POINT_XYZ'], tstD.dict['TRIANGLE_POINT'],
      tstD.dict['POINT_NEIGHBOR'])
  #find the minimum coordinates in each dimension
  mins = tstD.dict['POINT_XYZ'][0][1:]
  maxs = tstD.dict['POINT_XYZ'][0][1:]
  for point in convexHullPoints:
    xyz = tstD.dict['POINT_XYZ'][point-1][1:]
    for coord in range(3):
      mins[coord] = min(mins[coord], xyz[coord])
      maxs[coord] = max(maxs[coord], xyz[coord])
  #floor and ceiling everything to whole numbers
  minsSave = mins[:]
  maxsSave = maxs[:]
  #cache a bunch of computations on each triangle
  #triTuples does not contain the ones that are for cavities
  triTuples = geometry.cacheTriangle(
      tstD.dict['TRIANGLE_POINT'], tstD.dict['POINT_XYZ'], allTris)
  convexTriTuples = geometry.cacheTriangle(
      tstD.dict['CONVEX_HULL_TRI_POINT_LIST'], tstD.dict['POINT_XYZ'])
  #grid encoding -1 = outside ch, 0 = between ch, ms, -2 = inside ms
  mins, maxs = [], []  # so values computed are saved
  mins = [math.floor(x)-gridSize for x in minsSave]
  maxs = [math.ceil(x)+gridSize for x in maxsSave]
  gridD = grid.makeNewEmptyGrid(mins, maxs, gridSize, -1)  # set to outside ch
  #tstdebug.debugGridCountVals(gridD)
  #first step, check and set outside ch
  orstHelper.decideInside(
      gridD, convexTriTuples, convexHullPoints, 
      tstD.dict['CONVEX_HULL_POINT_TRI_LIST'], tstD.dict['POINT_XYZ'],
      tstD.dict['CONVEX_HULL_TRI_POINT_LIST'], 0, False, 2)
   #0 inside convex hull, False any value, 2=max tris
  #tstdebug.debugGridCountVals(gridD)
  #now find inside molecular surface
  orstHelper.decideInside(
      gridD, triTuples, allPoints, tstD.dict['POINT_TRIANGLE'],
      tstD.dict['POINT_XYZ'], tstD.dict['TRIANGLE_POINT'], -2)
  #-2 inside everything
  #tstdebug.debugGridCountVals(gridD)
  #following lines make files that allow visual debugging in pymol
  #tstdebug.debugGrid(gridD, "debug.grid." + str(gridSize) + ".py")
  #tstdebug.debugGridNoBlue(gridD, "debug.grid.nb." + str(gridSize) + ".py")
  #tstdebug.debugGridJustGreen(gridD, "debug.grid.jg." + str(gridSize) + ".py")
  #now...
  #here's the (relatively simple) surface travel distance calculation finally
  #  assign following encoding -1 = outside ch, 0 = on border,
  #   pos ints = dist from border, -2 = far inside ms,
  #   other neg ints = -(dist)-3
  #whole algorithm wrapped into big function...
  extraEdges, surfaceEdgeBoxes = grid.findLongSurfEdges(
     tstD.dict['POINT_XYZ'], tstD.dict['POINT_NEIGHBOR'], gridSize, mins, maxs)
  volumePoints = False
  if 'POINT_TRAVEL_DEPTH_CHECK' in tstD.dict.keys():
    volumePoints = tstD.dict['POINT_TRAVEL_DEPTH_CHECK']
  pointTravelDist,traceback,volumePointDepths = \
          travelDistNoMesh.calcTravelDist(gridD, \
            tstD.dict['POINT_XYZ'], gridSize, mins, maxs, allPoints, \
            extraEdges, surfaceEdgeBoxes, tstFileName,volumePoints)
  #tstdebug.debugGridCountVals(gridD)
  #reassign -1 to max+1 (cavities)
  maximumTD = max([xxx[1] for xxx in pointTravelDist])
  for onePoint in pointTravelDist:
    if onePoint[1] < 0.0:
      onePoint[1] = maximumTD + 1.0
  maxTD = grid.finalizeGridTravelDist(gridD, gridSize)
  phiDataOut = phi()
  phiDataOut.createFromGrid(gridD,gridSize,toplabel="travel depth")
  phiDataOut.write(tstFileName+".travel.oldmethod.phi")
  #more pymol debugging if desired
  #tstdebug.debugTravelGrid(gridD,"debug.travel.grid." + str(gridSize) + ".py",maximumTD)
  #tstdebug.debugTravelSurfGrid(gridD, "debug.travel.surf.grid." + str(gridSize) + ".py",extraEdges,mins,maxs,gridSize,maximumTD)
  #save data into tstData
  tstD.dict['DEPTH_TRAVEL_DIST'] = pointTravelDist
  #write data to file
  tstFile = open(tstFileName, 'a')
  tstFile.write("DEPTH_TRAVEL_DIST\n")
  for line in pointTravelDist:
    lineOut = "%8d" % line[0]
    for count in xrange(1,len(line)):
      lineOut += "%+9.4f " % line[count]
    noPlusLine = string.replace(lineOut, "+", " ")
    tstFile.write(noPlusLine)
    tstFile.write("\n")
  tstFile.write("END DEPTH_TRAVEL_DIST\n")
  tstFile.close()
  if volumePointDepths:
    tstD.dict['POINT_TRAVEL_DEPTH_REPORT'] = volumePointDepths
    tstFile = open(tstFileName, 'a')
    tstFile.write("POINT_TRAVEL_DEPTH_REPORT\n")
    for line in volumePointDepths:
      lineOut = "%8d" % line[0]
      for count in xrange(1,len(line)):
        lineOut += "%+9.4f " % line[count]
      noPlusLine = string.replace(lineOut, "+", " ")
      tstFile.write(noPlusLine)
      tstFile.write("\n")
    tstFile.write("END POINT_TRAVEL_DEPTH_REPORT\n")
    tstFile.close()
  #save tracebacks into tstData
  tracebackArray = [] #format is end, dist, list of starts
  tracebackKeys = traceback.keys()
  tracebackKeys.sort()
  for endkey in tracebackKeys:
    startList,end,dist = traceback[endkey]
    tbLine = [end[0], end[1], end[2], dist]
    for start in startList:
      tbLine.append(start[0])
      tbLine.append(start[1])
      tbLine.append(start[2])
    tracebackArray.append(tbLine)
  tstD.dict['TRACEBACK_LIST'] = tracebackArray
  #now write to file
  tstFile = open(tstFileName, 'a')
  tstFile.write("TRACEBACK_LIST\n")
  for line in tracebackArray:
    lineOut = " "
    for count in xrange(0,len(line)):
      lineOut += "%+9.4f " % line[count]
    noPlusLine = string.replace(lineOut, "+", " ")
    tstFile.write(noPlusLine)
    tstFile.write("\n")
  tstFile.write("END TRACEBACK_LIST\n")
  tstFile.close()

def tstTravelDepthRun(tstD, phiData, tstFileName="temp.tst"):
  if 'CONVEX_HULL_TRI_POINT_LIST' not in tstD.dict.keys():
    print "Run tstConvexHull.py on this tst data file first."
    sys.exit(1)
  #these sets are useful to construct
  convexHullPoints = set()
  for record in tstD.dict['CONVEX_HULL_TRI_POINT_LIST']:
    convexHullPoints.update(record[1:])
  gridD, mins, maxs =  grid.makeTrimmedGridFromPhi(
      phiData, tstD.dict['POINT_XYZ'], convexHullPoints, 0.6, -2.0, -1.0)
  gridSize = 1.0/phiData.scale
  del phiData  # no longer needed in this function, so delete this reference
  #tstdebug.debugGridCountVals(gridD)
  #do the biggest disjoint set of tris/points stuff
  allPoints, allTris, cavPoints, cavTris = cavity.assumeNoCavities(
      tstD.dict['POINT_XYZ'], tstD.dict['TRIANGLE_POINT'],
      tstD.dict['POINT_NEIGHBOR'])
  convexTriTuples = geometry.cacheTriangle( \
                tstD.dict['CONVEX_HULL_TRI_POINT_LIST'], tstD.dict['POINT_XYZ'])
  #grid encoding -1 = outside ch, 0 = between ch, ms, -2 = inside ms
  #this function marks cavities as interior
  #cavities now removed already
  #cavTriTuples = geometry.cacheTriangle( \
  #                tstD.dict['TRIANGLE_POINT'], tstD.dict['POINT_XYZ'], cavTris)
  #orstHelper.decideInside(gridD, cavTriTuples, cavPoints, \
  #               tstD.dict['POINT_TRIANGLE'], tstD.dict['POINT_XYZ'], \
  #               tstD.dict['TRIANGLE_POINT'], -2, -1) #-2 inside everything, -1 only valid to change
  #tstdebug.debugGridCountVals(gridD)
  #this marks the editable volume between the two surfaces
  orstHelper.decideInside(gridD, convexTriTuples, convexHullPoints, \
                 tstD.dict['CONVEX_HULL_POINT_TRI_LIST'], \
                 tstD.dict['POINT_XYZ'], \
                 tstD.dict['CONVEX_HULL_TRI_POINT_LIST'], \
                 0, -1, 2) #0 inside convex hull, -1 outside (only valid to change), 2 = max tris
  #tstdebug.debugGridCountVals(gridD)
  #following lines make files that allow visual debugging in pymol
  #tstdebug.debugGrid(gridD, "debug.grid." + str(gridSize) + ".py")
  #tstdebug.debugGridNoBlue(gridD, "debug.grid.nb." + str(gridSize) + ".py")
  #tstdebug.debugGridJustGreen(gridD, "debug.grid.jg." + str(gridSize) + ".py")
  #now...
  #here's the (relatively simple) surface travel distance calculation finally
  #  assign following encoding -1 = outside ch, 0 = on border,
  #   pos ints = dist from border, -2 = far inside ms,
  #   other neg ints = -(dist)-3
  #whole algorithm wrapped into big function...
  extraEdges, surfaceEdgeBoxes = grid.findLongSurfEdges(  \
     tstD.dict['POINT_XYZ'], tstD.dict['POINT_NEIGHBOR'], \
                    gridSize, mins, maxs)
  volumePoints = False
  if 'POINT_TRAVEL_DEPTH_CHECK' in tstD.dict.keys():
    volumePoints = tstD.dict['POINT_TRAVEL_DEPTH_CHECK']
  pointTravelDist,traceback,volumePointDepths = \
         travelDistNoMesh.calcTravelDist(gridD, \
           tstD.dict['POINT_XYZ'], gridSize, mins, maxs, allPoints, \
           extraEdges, surfaceEdgeBoxes, tstFileName, volumePoints)
  #tstdebug.debugGridCountVals(gridD)
  #transform grid to actual travel distance
  maxTD = grid.finalizeGridTravelDist(gridD, gridSize)
  #more pymol debugging if desired
  #tstdebug.debugTravelGrid(gridD, "debug.travel.grid." + str(gridSize) + ".py",maximumTD)
  #tstdebug.debugTravelSurfGrid(gridD, "debug.travel.surf.grid." + str(gridSize) + ".py",extraEdges,mins,maxs,gridSize,maximumTD)
  #save data into tstD
  tstD.dict['DEPTH_TRAVEL_DIST'] = pointTravelDist
  #save tracebacks into tstD
  tracebackArray = [] #format is end, dist, list of starts
  tracebackKeys = traceback.keys()
  tracebackKeys.sort()
  for endkey in tracebackKeys:
    startList,end,dist = traceback[endkey]
    tbLine = [end[0], end[1], end[2], dist]
    for start in startList:
      tbLine.append(start[0])
      tbLine.append(start[1])
      tbLine.append(start[2])
    tracebackArray.append(tbLine)
  tstD.dict['TRACEBACK_LIST'] = tracebackArray
  phiDataOut = phi()
  phiDataOut.createFromGrid(gridD,gridSize,toplabel="travel depth", \
                            defaultValue=maxTD)
  return gridD,phiDataOut,volumePointDepths

#took all the hole/path finding stuff out....
def tstTravelSurfOutside(tstFileName,phiFileName=False, \
                         tstDataIn=False,phiDataIn=False,borderSize=10):
  if tstDataIn: #don't want to read in if calling function already did it for us
    tstD = tstDataIn
  else:
    tstD = tstdata.tstData(tstFileName) #read the file into the data structure
  if phiDataIn:
    phiData = phiDataIn
  else:
    phiData = phi(phiFileName)  #read in the phimap if possible
  #depends on convex hull _only_ to decide how far out to put the border
  if 'CONVEX_HULL_TRI_POINT_LIST' not in tstD.dict.keys():
    print "Run tstConvexHull.py on this tst data file first."
    sys.exit(1)
  #these sets are useful to construct
  convexHullPoints = set()
  for record in tstD.dict['CONVEX_HULL_TRI_POINT_LIST']:
    convexHullPoints.update(record[1:])
  gridD,mins,maxs =  grid.makeTrimmedGridFromPhi(phiData,\
                tstD.dict['POINT_XYZ'], \
                convexHullPoints, 0.6, -2.0, 0, borderSize)
  gridSize = 1.0/phiData.scale
  del phiData #no longer needed in this function, so delete this reference
  #do the biggest disjoint set of tris/points stuff
  allPoints, allTris, cavPoints, cavTris = cavity.assumeNoCavities( \
      tstD.dict['POINT_XYZ'], tstD.dict['TRIANGLE_POINT'], \
      tstD.dict['POINT_NEIGHBOR'])
  #here's the (relatively simple) surface travel distance calculation finally
  #  assign following encoding -1 = outside ch, 0 = on border,
  #   pos ints = dist from border, -2 = far inside ms,
  #   other neg ints = -(dist)-3
  extraEdges, surfaceEdgeBoxes = grid.findLongSurfEdges(  \
     tstD.dict['POINT_XYZ'], tstD.dict['POINT_NEIGHBOR'], \
                    gridSize, mins, maxs)
  for surfaceEdgeBox in surfaceEdgeBoxes.keys():
    x,y,z=gridD[surfaceEdgeBox[0]][surfaceEdgeBox[1]][surfaceEdgeBox[2]][1:]
    gridD[surfaceEdgeBox[0]][surfaceEdgeBox[1]][surfaceEdgeBox[2]] = (-1.,x,y,z)
  pointTravelDist,traceback,volumePointDepths = \
        travelDistNoMesh.calcTravelDist(gridD, \
          tstD.dict['POINT_XYZ'], gridSize, mins, maxs, allPoints, \
          extraEdges, surfaceEdgeBoxes, tstFileName)
  #transform grid to actual travel distance
  maxTD = grid.finalizeGridTravelDist(gridD, gridSize)
  gridMaximaRanking = grid.getMaximaRanking(gridD)
  #output the grids and the gridmaxima....
  phiDataOut = phi()
  phiDataOut.createFromGrid(gridMaximaRanking,gridSize,toplabel="travel depth maxima rank")
  phiDataOut.write(tstFileName+".travel.out.max.rank.phi")
  phiDataOut = phi()
  phiDataOut.createFromGrid(gridD,gridSize,toplabel="travel depth surf-out",defaultValue=maxTD+1.0)
  phiDataOut.write(tstFileName+".travel.out.phi")
  return mins,maxs,gridSize,convexHullPoints,allPoints,allTris, \
      extraEdges,surfaceEdgeBoxes,pointTravelDist,traceback,maxTD,gridD, \
      gridMaximaRanking

#run tstConvexHull.py first...
#this runs all the other functions and sets up necessary stuff, writes to disk
def tstTravelDepth(tstFileName,phiFileName):
  tstD = tstdata.tstData(tstFileName) #read the file into the data structure
  phiData = phi(phiFileName)  #read in the phimap if possible
  phiTravelDepthGrid,phiTravelDepthData,volumePointDepths = \
           tstTravelDepthRun(tstD, phiData, tstFileName) #modifies tstD in place
  #transform grid to actual travel distance
  phiTravelDepthData.write(tstFileName+".travel.phi")
  #write data to file
  tstFile = open(tstFileName, 'a')
  tstFile.write("DEPTH_TRAVEL_DIST\n")
  for line in tstD.dict['DEPTH_TRAVEL_DIST']:
    lineOut = "%8d" % line[0]
    for count in xrange(1,len(line)):
      lineOut += "%+9.4f " % line[count]
    noPlusLine = string.replace(lineOut, "+", " ")
    tstFile.write(noPlusLine)
    tstFile.write("\n")
  tstFile.write("END DEPTH_TRAVEL_DIST\n")
  tstFile.close()
  if volumePointDepths:
    tstD.dict['POINT_TRAVEL_DEPTH_REPORT'] = volumePointDepths
    tstFile = open(tstFileName, 'a')
    tstFile.write("POINT_TRAVEL_DEPTH_REPORT\n")
    for line in volumePointDepths:
      lineOut = "%8d" % line[0]
      for count in xrange(1,len(line)):
        lineOut += "%+9.4f " % line[count]
      noPlusLine = string.replace(lineOut, "+", " ")
      tstFile.write(noPlusLine)
      tstFile.write("\n")
    tstFile.write("END POINT_TRAVEL_DEPTH_REPORT\n")
    tstFile.close()
  #now write to file
  tstFile = open(tstFileName, 'a')
  tstFile.write("TRACEBACK_LIST\n")
  for line in tstD.dict['TRACEBACK_LIST']:
    lineOut = " "
    for count in xrange(0,len(line)):
      lineOut += "%+9.4f " % line[count]
    noPlusLine = string.replace(lineOut, "+", " ")
    tstFile.write(noPlusLine)
    tstFile.write("\n")
  tstFile.write("END TRACEBACK_LIST\n")
  tstFile.close()

def tstTravelSurfInside(tstFileName,phiFileName=False):
  tstD = tstdata.tstData(tstFileName) #read the file into the data structure
  phiData = phi(phiFileName)  #read in the phimap if possible
  if 'CONVEX_HULL_TRI_POINT_LIST' not in tstD.dict.keys():
    print "Run tstConvexHull.py on this tst data file first."
    sys.exit(1)
  #these sets are useful to construct
  convexHullPoints = set()
  for record in tstD.dict['CONVEX_HULL_TRI_POINT_LIST']:
    convexHullPoints.update(record[1:])
  gridD,mins,maxs =  grid.makeTrimmedGridFromPhi(phiData,\
                tstD.dict['POINT_XYZ'], \
                convexHullPoints, 0.6, 0, -2.0, 2)
  gridSize = 1.0/phiData.scale
  del phiData #no longer needed in this function, so delete this reference
  #do the biggest disjoint set of tris/points stuff
  allPoints, allTris, cavPoints, cavTris = cavity.assumeNoCavities( \
      tstD.dict['POINT_XYZ'], tstD.dict['TRIANGLE_POINT'], \
      tstD.dict['POINT_NEIGHBOR'])
  #here's the (relatively simple) surface travel distance calculation finally
  #  assign following encoding -1 = outside ch, 0 = on border,
  #   pos ints = dist from border, -2 = far inside ms,
  #   other neg ints = -(dist)-3
  #whole algorithm wrapped into big function...
  extraEdges, surfaceEdgeBoxes = grid.findLongSurfEdges(  \
     tstD.dict['POINT_XYZ'], tstD.dict['POINT_NEIGHBOR'], \
                    gridSize, mins, maxs)
  for surfaceEdgeBox in surfaceEdgeBoxes.keys():
    x,y,z=gridD[surfaceEdgeBox[0]][surfaceEdgeBox[1]][surfaceEdgeBox[2]][1:]
    gridD[surfaceEdgeBox[0]][surfaceEdgeBox[1]][surfaceEdgeBox[2]] = (-1.,x,y,z)
  pointTravelDist,traceback,volumePointDepths = \
        travelDistNoMesh.calcTravelDist(gridD, \
           tstD.dict['POINT_XYZ'], gridSize, mins, maxs, allPoints, \
           extraEdges, surfaceEdgeBoxes, tstFileName)
  #transform grid to actual travel distance
  maxTD = grid.finalizeGridTravelDist(gridD, gridSize)
  phiDataOut = phi()
  phiDataOut.createFromGrid(gridD,gridSize,toplabel="travel depth surf-in")
  phiDataOut.write(tstFileName+".travel.in.phi")
  pdbD = pdbData()
  for line in tstD.dict['PDB_RECORD']:
    pdbD.processLine(line)
  atomTravelInDepths = grid.assignAtomDepths(gridD, gridSize, mins, maxs, pdbD)
  #make a pdb file with the bfactor replaced
  for index,atomTID in enumerate(atomTravelInDepths):
    pdbD.updateFactors(index, (pdbD.factors[index][0], atomTID))
  pdbD.write(tstFileName+".travelin.pdb")
  #also add record to tstdata
  atomTIDRecord = []
  for index,atomTID in enumerate(atomTravelInDepths):
    atomTIDRecord.append([index+1, atomTID])
  tstD.dict['ATOM_TRAVEL_IN'] = atomTIDRecord
  #write data into tst file
  tstFile = open(tstFileName, 'a')
  tstFile.write("ATOM_TRAVEL_IN\n")
  for line in tstD.dict['ATOM_TRAVEL_IN']:
    lineOut = "%8d" % line[0]
    for count in xrange(1,len(line)):
      lineOut += "%+9.4f " % line[count]
    noPlusLine = string.replace(lineOut, "+", " ")
    tstFile.write(noPlusLine)
    tstFile.write("\n")
  tstFile.write("END ATOM_TRAVEL_IN\n")
  tstFile.close()

def tstTravelSurfInsideOld(tstFileName,phiFileName=False):
  '''does the old algorithm of just computing the shortest distance to any
  surface point from any atom by going through both lists the hard way'''
  tstD = tstdata.tstData(tstFileName) #read the file into the data structure
  #do the biggest disjoint set of tris/points stuff
  allPoints, allTris, cavPoints, cavTris = cavity.assumeNoCavities( \
      tstD.dict['POINT_XYZ'], tstD.dict['TRIANGLE_POINT'], \
      tstD.dict['POINT_NEIGHBOR'])
  pointXyz = tstD.dict['POINT_XYZ']
  pdbD = pdbData()
  for line in tstD.dict['PDB_RECORD']:
    pdbD.processLine(line)
  atomTravelInDepths = []
  for coord in pdbD.coords:
    minDist = geometry.distL2(coord, pointXyz[allPoints[0]-1][1:]) #set to first
    for point in allPoints[1:]:
      thisDist = geometry.distL2(coord, pointXyz[point-1][1:])
      minDist = min(minDist, thisDist)
    atomTravelInDepths.append(minDist)
  #make a pdb file with the bfactor replaced
  for index,atomTID in enumerate(atomTravelInDepths):
    pdbD.updateFactors(index, (pdbD.factors[index][0], atomTID))
  pdbD.write(tstFileName+".old.atomdepth.pdb")
  #also add record to tstdata
  atomTIDRecord = []
  for index,atomTID in enumerate(atomTravelInDepths):
    atomTIDRecord.append([index+1, atomTID])
  tstD.dict['ATOM_DEPTH_OLD'] = atomTIDRecord
  #write data into tst file
  tstFile = open(tstFileName, 'a')
  tstFile.write("ATOM_TRAVEL_IN\n")
  for line in tstD.dict['ATOM_DEPTH_OLD']:
    lineOut = "%8d" % line[0]
    for count in xrange(1,len(line)):
      lineOut += "%+9.4f " % line[count]
    noPlusLine = string.replace(lineOut, "+", " ")
    tstFile.write(noPlusLine)
    tstFile.write("\n")
  tstFile.write("END ATOM_DEPTH_OLD\n")
  tstFile.close()

def copyPhiMap(phiFileIn, phiFileOut):
  '''simple, just copies the phi map and writes it, to check binary i/o'''
  phiData = phi(phiFileIn)  #read in the phimap if possible
  phiData.write(phiFileOut) #write out

def examinePhiMap(phiFileIn):
  phiData = phi(phiFileIn)  #read in the phimap if possible
  print phiData.getMinMaxValues()
  print phiData.histogramValues()

#this is where main is... maybe add some other arguments like gridSize?
if -1 != string.find(sys.argv[0], "oldTravelDist"):
  if 1 < len(sys.argv) and sys.argv[1] == "depth" and  len(sys.argv) > 3:
    tstFile, phiFile = sys.argv[2:4]
    print tstFile, phiFile
    tstTravelDepth(tstFile, phiFile)
  elif 1 < len(sys.argv) and sys.argv[1] == "depthold"  and len(sys.argv) > 2:
    tstFile = sys.argv[2]
    gridSize = 1.0
    if len(sys.argv) > 3:
      gridSize = sys.argv[3]
    print tstFile
    tstTravelDepthNoPhi(tstFile, gridSize)
  elif 1 < len(sys.argv) and sys.argv[1] == "surfout" and  len(sys.argv) > 3:
    tstFile, phiFile = sys.argv[2:4]
    print tstFile, phiFile
    tstTravelSurfOutside(tstFile, phiFile)
  elif 1 < len(sys.argv) and sys.argv[1] == "surfin" and  len(sys.argv) > 3:
    tstFile, phiFile = sys.argv[2:4]
    print tstFile, phiFile
    tstTravelSurfInside(tstFile, phiFile)
  elif 1 < len(sys.argv) and sys.argv[1] == "surfinold" and  len(sys.argv) > 3:
    tstFile, phiFile = sys.argv[2:4]
    print tstFile, phiFile
    tstTravelSurfInsideOld(tstFile, phiFile)
  elif 1 < len(sys.argv) and sys.argv[1] == "cavityremove" and  len(sys.argv) > 3:
    tstFile, tstOut, phiFile, phiOut = sys.argv[2:6]
    print tstFile, tstOut, phiFile, phiOut
    cavity.tstCavityRemoval(tstFile, tstOut, phiFile, phiOut)
  elif 1 < len(sys.argv) and sys.argv[1] == "copyphi" and  len(sys.argv) > 3:
    phiFileIn, phiFileOut = sys.argv[2:4]
    print phiFileIn, phiFileOut
    copyPhiMap(phiFileIn, phiFileOut)
  elif 1 < len(sys.argv) and sys.argv[1] == "examinephi" and  len(sys.argv) > 2:
    phiFileIn= sys.argv[2]
    print phiFileIn
    examinePhiMap(phiFileIn)
  else:
    print "Usage: tstTravelDist.py depth tstFile phiFile"
    print "Usage: tstTravelDist.py depthold tstFile"
    print "Usage: tstTravelDist.py surfout tstFile phiFile"
    print "Usage: tstTravelDist.py surfin tstFile phiFile"
    print "Usage: tstTravelDist.py surfinold tstFile phiFile"
    print "Usage: tstTravelDist.py cavityremove tstFile tstFileOut phiFile phiFileOut"
    print "Usage: tstTravelDist.py copyphi phiFileIn phiFileOut"
    print "Usage: tstTravelDist.py examinephi phiFileIn"
