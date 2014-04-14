#!/usr/bin/env python

import math
import string
import sys
import os #for running commandline instead of through pymol
from os.path import isfile
import cPickle
#following custom module must be in same directory or in usual path
import tstdata
from sharp_phi import phi #phi map reader/writer
import tstdebug #debugging functions, moved out of codebase
import tstTopology  # gets strip surrounding handles
from comparePaths import comparePathsManyMetrics  # does pRMSD etc
import mesh
import geometry  # geometric primitives
from pdb import pdbData
import grid  # moved primitive grid functions
import paths  # moved primitive path functions
import orstHelper
import cavity
import charge
import statistics
import priodict
import tm3
import tstCurvature

def meshConstruct(tstD, phiData, tstFileName="temp.tst", borderSize=2, \
                  threshold="auto", cavities=False):
  '''returns a phi and a grid, modifies tstD in place.'''
  if 'CONVEX_HULL_TRI_POINT_LIST' not in tstD.dict.keys():
    print "Run tstConvexHull.py on this tst data file first."
    sys.exit(1)
  #these sets are useful to construct
  convexHullPoints = set()
  for record in tstD.dict['CONVEX_HULL_TRI_POINT_LIST']:
    convexHullPoints.update(record[1:])
  maxPhi = phiData.getMaxValues()
  if threshold == "auto" and maxPhi == 1.0:
    threshold = 0.6
  if threshold == "auto" and maxPhi == 10.0:
    threshold = 6.0
  gridD,mins,maxs =  grid.makeTrimmedGridFromPhi(phiData,\
                tstD.dict['POINT_XYZ'], \
                convexHullPoints, threshold, -2.0, -1.0, borderSize)
  gridSize = 1.0/phiData.scale
  del phiData #no longer needed in this function, so delete this reference
  #needs border = 2 so all surface points have a grid point on either side...
  #tstdebug.debugGridCountVals(gridD)
  #tstdebug.debugGridNoBlue(gridD, "debug.phigrid.py")
  #do the biggest disjoint set of tris/points stuff
  if not cavities:
    allPoints, allTris, cavPoints, cavTris = cavity.assumeNoCavities( \
        tstD.dict['POINT_XYZ'], tstD.dict['TRIANGLE_POINT'], \
        tstD.dict['POINT_NEIGHBOR'])
    cavPointLists = False
  else:
    allPoints, allCavPoints, cavPointLists = \
        cavity.findBiggestDisjointSetsBreakCavities( \
          tstD.dict['POINT_XYZ'], tstD.dict['TRIANGLE_POINT'], \
          tstD.dict['POINT_NEIGHBOR'])
  convexTriTuples = geometry.cacheTriangle( \
           tstD.dict['CONVEX_HULL_TRI_POINT_LIST'], tstD.dict['POINT_XYZ'])
  #this marks the editable volume between the two surfaces
  orstHelper.decideInside(gridD, convexTriTuples, convexHullPoints, \
                 tstD.dict['CONVEX_HULL_POINT_TRI_LIST'], \
                 tstD.dict['POINT_XYZ'], \
                 tstD.dict['CONVEX_HULL_TRI_POINT_LIST'], \
                 0, -1, 2) #0 inside convex hull,
                           # -1 outside (only valid to change), 2 = max tris
  #grid encoding -1 = outside ch, 0 = between ch, ms, -2 = inside ms
  #tstdebug.debugGridNoBlue(gridD, "debug.phigrid.py")
  meshData = mesh.mesh(gridD, tstD.dict['POINT_XYZ'], \
                       tstD.dict['POINT_NEIGHBOR'], \
                       gridSize, -1, -2, 0, allPoints, cavPointLists)
  #now mesh is initialized
  return meshData, gridD

def tstTravelDepthMeshRun(tstD, phiData, tstFileName="temp.tst",borderSize=2, \
                          threshold="auto", cavities=False):
  '''returns a phi and a grid, modifies tstD in place.'''
  meshData, gridData = meshConstruct(tstD, phiData, tstFileName, \
                                   borderSize=borderSize, threshold=threshold, \
                                   cavities=cavities)
  if cavities:
    validList = [2,3,5]
  else:
    validList = [2,3]
  meshData.calculateTravelDistance("traveldepth", [0], validList)
  pointTravelDepth = meshData.getSurfaceTravelDistance("traveldepth")
  gridTravelDepth = meshData.getGridTravelDistance(gridData, "traveldepth")
  tstD.dict['DEPTH_TRAVEL_DIST'] = pointTravelDepth #save to tstD...
  phiDataOut = phi()
  gridSize = 1.0/phiData.scale
  phiDataOut.createFromGrid(gridTravelDepth,gridSize, \
                            toplabel="travel depth",defaultValue=-1)
  return gridTravelDepth,phiDataOut,meshData

def tstTravelDepthMesh(tstFileName,phiFileName, ligandFileName=None, \
                       cavities=False, threshold="auto"):
  '''sets up a normal travel depth mesh run, calls tstTravelDepthMeshRun.
  if cavities is set, it means we want the (new) travel depth of those too'''
  distanceName = 'traveldepth'
  tstD = tstdata.tstData(tstFileName) #read the file into the data structure
  phiData = phi(phiFileName)  #read in the phimap if possible
  phiGridSpacing = 1./phiData.scale
  phiTravelDepthGrid,phiTravelDepth,meshData = \
       tstTravelDepthMeshRun(tstD, phiData, tstFileName, cavities=cavities, \
                             threshold=threshold) #modifies tstD in place
  if ligandFileName is not None:
    ligand = pdbData(ligandFileName)
    ligandXYZR = ligand.getHeavyAtomXYZRadius()
    betweenList = meshData.getBetweenNodes()
    surfaceList = meshData.getSurfaceNodes()
    nodeWithinSet = meshData.getWithinNodesNoInside(ligandXYZR)
    #print len(nodeWithinSet), len(ligandXYZR)
    tracebackSet = meshData.getTracebackSet(nodeWithinSet, distanceName)
    #print len(tracebackSet)
    minW, maxW, meanW = meshData.getMinMaxMeanNodes(nodeWithinSet, distanceName)
    minT, maxT, meanT = meshData.getMinMaxMeanNodes(tracebackSet, distanceName)
    minB, maxB, meanB = meshData.getMinMaxMeanNodes(betweenList, distanceName)
    minS, maxS, meanS = meshData.getMinMaxMeanNodes(surfaceList, distanceName)
    #print all on one line, header printed earlier
    print minW, maxW, meanW,
    print minT, maxT, meanT,
    print minB, maxB, meanB,
    print minS, maxS, meanS,
    volumeWithin = len(nodeWithinSet)*phiGridSpacing**3.
    volumeTrace = len(tracebackSet)*phiGridSpacing**3.
    print volumeWithin, volumeTrace #newline wanted here so no comma
    #print phiGridSpacing
    listWithin, listTrace = [],[]
    for node in nodeWithinSet:
      listWithin.append(node.getXYZ())
    for node in tracebackSet:
      listTrace.append(node.getXYZ())
    tstdebug.pointDebug(listWithin, filename=ligandFileName+".within.py")
    tstdebug.pointDebug(listTrace, filename=ligandFileName+".trace.py")
  #tstdebug.debugGridCountVals(phiTravelDepthGrid)
  #transform grid to actual travel distance
  phiTravelDepth.write(tstFileName+".travel.phi")
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

def tstCountHoles(tstFileName):
  tstD = tstdata.tstData(tstFileName) #read the file into the data structure
  numberHandles = tstD.countHandles()
  return numberHandles

#need to finish making this do as little duplicate work as surfout as possible
def tstTravelFindHoles(tstFileName,phiFileName,debugOut=False,borderSize=2,\
                       nearbyDistance=4.):
  '''if debugout is set, additional files are created.
  bordersize can change the amount of extra space around the protein.
  nearbydistance changes the distance that nearby residues are gathered from.'''
  print "reading tst and phi files"
  tstD = tstdata.tstData(tstFileName, necessaryKeys= \
                    tstdata.tstData.necessaryKeysForHoles)
  phiData = phi(phiFileName)  #read in the phimap if possible
  gridSize = 1.0/phiData.scale #needed later
  numberHandles = tstD.countHandles()
  print "there are ", numberHandles, " holes in this structure"
  print "running travel depth now"
  phiTravelDepthGrid,phiTravelDepthData,meshData = tstTravelDepthMeshRun( \
        tstD, phiData, tstFileName,borderSize=borderSize,threshold="auto")
  del phiData,phiTravelDepthGrid,phiTravelDepthData #not needed, reclaim memory
  print "calculating travel out distance"
  meshData.calculateTravelDistance("surfout", [3], [0,2])
  print "finding holes"
  loopTrisSave, loopPointsSave, regLoopTris, regLoopPts, pointNeighbors, \
      pointNeighborsNodes, outsidePoints, outsidePointsNodes, possHoleStarts = \
                tstTopology.fillInHolesAndGrow(tstD.dict['TRIANGLE_POINT'], \
                     tstD.dict['POINT_TRIANGLE'], \
                     tstD.dict['POINT_XYZ'], tstD.dict['NORM_XYZ'], \
                     numberHandles, tstFileName, debugOut, meshData, "surfout")
  if debugOut:
    tstdebug.debugTriangleList(regLoopTris, tstD.dict['TRIANGLE_POINT'], tstD.dict['POINT_XYZ'], tstFileName + ".regular.loops.py")
  tstDPointXYZ = tstD.dict['POINT_XYZ'] #save for later
  tstDPdbRecord = tstD.dict['PDB_RECORD'] #save for later
  del tstD #rest of tstD isn't needed, so free memory
  #output the possible places where HOLE could start... centers of regular plugs
  print "writing output files"
  paths.outputNodesText(possHoleStarts, tstFileName+".HOLE.start.txt")
  pathsList = meshData.getPaths("surfout", pointNeighbors, outsidePoints, \
                            possHoleStarts)
  del meshData #deletes everything that has no other refs (from paths)
  #print len(pointNeighbors), len(outsidePoints), len(possHoleStarts), len(pathsList)
  allPoints,outsidePts = [],[]
  for point in pointNeighborsNodes.keys():
    allPoints.append(point.pathXyz)
  for point in outsidePointsNodes:
    outsidePts.append(point.pathXyz)
  tstdebug.pointDebug(allPoints, filename=tstFileName+".tree.py")
  tstdebug.pointDebug(outsidePts, filename=tstFileName+".outside.py", mainColor=(.9,.1,.1), radius=0.55)
  #tstdebug.debugSetGridSpheres(pointNeighbors.keys(),gridSize,tstFileName+".tree.radius.py",radius=True,mainColor=(0.01,0.9,0.05)) #testing new output of tree with radius
  foundPaths = 0
  pathFile = string.replace(tstFileName, ".nocav.tst", ".py")  #very very specific... probably bad but nice... can always rerun standalone later.
  knownPathExists = isfile(pathFile)
  logName = tstFileName + ".findholes.log"
  logFile = open(logName, 'w')
  logFile.write("number endNumOne endNumTwo plugs stepLength pathLength pathMinRadius pathMaxInsideRadius endMinOne endMinTwo minimaCount travelDepthMax windingMetric avgTheta ")
  if knownPathExists:
    logFile.write("pRMSD coverage span wrmsd less1 lessrad radiicomp")
  logFile.write("\n")
  if knownPathExists:
    bestStats = [ "err","err","err","err","err","err","err" ]
    bestStatsPaths = [0,0,0,0,0,0,0]
  sortedByMinRadiusPaths = []
  for pathIndex, (outsideOne,outsideTwo,plugs,nodePath) in enumerate(pathsList):
    pointPath,spheres = [],[]
    for node in nodePath:
      pointPath.append(list(node.pathXyz))
      spheres.append(list(node.pathXyz))
      pointRadius = node.distances["surfout"] #add 'radius' info
      if not pointRadius or pointRadius == 0.:
        pointRadius = .000000001 #very small number, on surface
      pointPath[-1].insert(0, pointRadius)
      spheres[-1].append(pointRadius)
    minRad = paths.pathMinRadius(pointPath)
    newTuple = (minRad, outsideOne, outsideTwo,plugs,nodePath,pointPath,spheres)
    #insertion sort into new list
    position = 0
    while position < len(sortedByMinRadiusPaths) and \
                       sortedByMinRadiusPaths[position][0] > minRad:
      position += 1 #only if list can handle it and already inserted are bigger
    sortedByMinRadiusPaths.insert(position, newTuple)
  print "output files for individual paths"
  for pathIndex, newTuple in enumerate(sortedByMinRadiusPaths):
    (minRad, outsideOne, outsideTwo, plugs, nodePath,pointPath,spheres)=newTuple
    throughTris,throughPts = paths.checkPath(pointPath, loopPointsSave, \
                                       tstDPointXYZ)
    #if True: #testing this idea... possibly something wrong
    if throughTris: #it worked... make some more debugging output
      foundPaths += 1
      outName = tstFileName + "." + str(foundPaths)
      if debugOut and throughTris:
        tstdebug.debugTrianglesNotOrig(throughTris,  \
                       tstDPointXYZ, outName+".through.loop.py", \
                       ptList=throughPts)
      #always do these 2
      tstdebug.debugSetGridSpheres(pointPath,gridSize,outName+".pore.py",radius=True,mainColor=(0.01,0.9,0.05))
      tstdebug.debugSetGridSpheres(pointPath,gridSize,outName+".path.py",mainColor=(.01,0.95,0.9))
      #mesh.meshFromSpheres(spheres, 0.5, outName+".points.py")
      paths.outputRadiiTxt(pointPath, outName+".radii.txt")
      paths.outputNearbyResidues(pointPath, outName, \
                           tstDPdbRecord, nearbyDistance)
      pathLen = paths.pathLength(pointPath)
      minimaCount = paths.pathMinimaCount(pointPath)
      maxRad, endMinOne, endMinTwo = paths.insideTwoMinimaRadiusMax(pointPath)
      travelDepthMax = paths.pathMaxDistance(nodePath, 'traveldepth')
      windingMetric = paths.computeWindingMetric(pointPath)
      thetas, avgTheta = paths.averageTheta(pointPath)
      #print endMinOne, endMinTwo, minimaCount, travelDepthMax, windingMetric,  avgTheta
      #attempt to do pRMSD if possible...
      prmsd, coverage, span, wrmsd, less1, lessrad, radiicomp = "err","err","err","err","err","err","err"
      if knownPathExists:
        try:
          sourcePath,sourceRadii = [],[]
          for pathPt in pointPath:
            sourcePath.append(pathPt[1:4])
            sourceRadii.append(pathPt[0])
          prmsd, coverage, span, wrmsd, less1, lessrad, radiicomp =  \
                 comparePathsManyMetrics(False,pathFile,sourcePath, sourceRadii)
          if bestStats[0] == 'err' or bestStats[0] > prmsd:
            bestStats[0] = prmsd
            bestStatsPaths[0] = foundPaths
          if bestStats[1] == 'err' or bestStats[1] < coverage:
            bestStats[1] = coverage
            bestStatsPaths[1] = foundPaths
          if bestStats[2] == 'err' or bestStats[2] < span:
            bestStats[2] = span
            bestStatsPaths[2] = foundPaths
          if bestStats[3] == 'err' or bestStats[3] > wrmsd:
            bestStats[3] = wrmsd
            bestStatsPaths[3] = foundPaths
          if bestStats[4] == 'err' or bestStats[4] < less1:
            bestStats[4] = less1
            bestStatsPaths[4] = foundPaths
          if bestStats[5] == 'err' or bestStats[5] < lessrad:
            bestStats[5] = lessrad
            bestStatsPaths[5] = foundPaths
          if bestStats[6] == 'err' or bestStats[6] < radiicomp:
            bestStats[6] = radiicomp
            bestStatsPaths[6] = foundPaths
        except (IOError, TypeError): #if there is no known path  file, this should be the error
          pass
      #now output data
      logFile.write(str(foundPaths) + " ")
      logFile.write(str(outsideOne) + " ")
      logFile.write(str(outsideTwo) + " ")
      logFile.write(str(plugs) + " ")
      logFile.write(str(len(pointPath)) + " ")
      logFile.write(str(pathLen) + " ")
      logFile.write(str(minRad) + " ")
      logFile.write(str(maxRad) + " ")
      logFile.write(str(endMinOne) + " ")
      logFile.write(str(endMinTwo) + " ")
      logFile.write(str(minimaCount) + " ")
      logFile.write(str(travelDepthMax) + " ")
      logFile.write(str(windingMetric) + " ")
      logFile.write(str(avgTheta) + " ")
      if knownPathExists:
        logFile.write(str(prmsd) + " ")
        logFile.write(str(coverage) + " ")
        logFile.write(str(span) + " ")
        logFile.write(str(wrmsd) + " ")
        logFile.write(str(less1) + " ")
        logFile.write(str(lessrad) + " ")
        logFile.write(str(radiicomp) + " ")
      logFile.write("\n") #that's all
  logFile.close()
  if knownPathExists: #output bestStats and bestStatsPaths
    bestName = tstFileName + ".known.best.txt"
    bestFile = open(bestName, 'w')
    bestFile.write("pRMSD coverage span wrmsd less1 lessrad radiicomp ")
    bestFile.write("pRMSD# coverage# span# wrmsd# less1# lessrad# radiicomp#\n")
    for stat in bestStats:
      bestFile.write(str(stat) + " ")
    for stat in bestStatsPaths:
      bestFile.write(str(stat) + " ")
    bestFile.write("\n")
    bestFile.close()
  print "done with chunnel"

#took all the hole/path finding stuff out....
def tstTravelSurfOutsideMesh(tstFileName,phiFileName=None, \
                         tstDataIn=False,phiDataIn=None,borderSize=10, \
                         threshold="auto"):
  if tstDataIn: #don't want to read in if calling function already did it for us
    tstD = tstDataIn
  else:
    tstD = tstdata.tstData(tstFileName, \
            necessaryKeys=tstdata.tstData.necessaryKeysForMesh)
  if phiDataIn is not None:
    phiData = phiDataIn
  else:
    phiData = phi(phiFileName)  #read in the phimap if possible
  if 'CONVEX_HULL_TRI_POINT_LIST' not in tstD.dict.keys():
    print "Run tstConvexHull.py on this tst data file first."
    sys.exit(1)
  #these sets are useful to construct
  convexHullPoints = set()
  for record in tstD.dict['CONVEX_HULL_TRI_POINT_LIST']:
    convexHullPoints.update(record[1:])
  maxPhi = phiData.getMaxValues()
  if threshold == "auto" and maxPhi == 1.0:
    threshold = 0.6
  if threshold == "auto" and maxPhi == 10.0:
    threshold = 6.0
  gridD,mins,maxs =  grid.makeTrimmedGridFromPhi(phiData,\
                tstD.dict['POINT_XYZ'], \
                convexHullPoints, threshold, -2, 0, borderSize)
  gridSize = 1.0/phiData.scale
  del phiData #no longer needed in this function, so delete this reference
  #do the biggest disjoint set of tris/points stuff
  allPoints, allTris, cavPoints, cavTris = cavity.assumeNoCavities( \
      tstD.dict['POINT_XYZ'], tstD.dict['TRIANGLE_POINT'], \
      tstD.dict['POINT_NEIGHBOR'])
  #here is where code is mesh-specific
  meshData = mesh.mesh(gridD, tstD.dict['POINT_XYZ'], \
                  tstD.dict['POINT_NEIGHBOR'], \
                  gridSize, 0, -2, "X") #no between
  meshData.calculateTravelDistance("surfout", [3], [0,2])
  #transform grid to actual travel distance
  maxTD = grid.finalizeGridTravelDist(gridD, gridSize)
  gridMaximaRanking = grid.getMaximaRanking(gridD)
  #output the grids and the gridmaxima....
  phiDataOut = phi()
  phiDataOut.createFromGrid(gridMaximaRanking,gridSize,toplabel="travel depth maxima rank")
  phiDataOut.write(tstFileName+".mesh.travel.out.max.rank.phi")
  phiDataOut = phi()
  phiDataOut.createFromGrid(gridD,gridSize,toplabel="travel depth surf-out",defaultValue=maxTD+1.0)
  phiDataOut.write(tstFileName+".mesh.travel.out.phi")

def tstTravelSurfInsideMesh(tstFileName,phiFileName,threshold="auto"):
  '''calculates the burial depth'''
  print "reading in tst and phi files"
  tstD = tstdata.tstData(tstFileName, \
            necessaryKeys=tstdata.tstData.necessaryKeysForMesh+['PDB_RECORD'])
  phiData = phi(phiFileName)  #read in the phimap if possible
  if 'CONVEX_HULL_TRI_POINT_LIST' not in tstD.dict.keys():
    print "Run tstConvexHull.py on this tst data file first."
    sys.exit(1)
  #these sets are useful to construct
  convexHullPoints = set()
  for record in tstD.dict['CONVEX_HULL_TRI_POINT_LIST']:
    convexHullPoints.update(record[1:])
  maxPhi = phiData.getMaxValues()
  if threshold == "auto" and maxPhi == 1.0:
    threshold = 0.6
  if threshold == "auto" and maxPhi == 10.0:
    threshold = 6.0
  gridD,mins,maxs =  grid.makeTrimmedGridFromPhi(phiData,\
                tstD.dict['POINT_XYZ'], \
                convexHullPoints, threshold, 0, -2, 2)
  gridSize = 1.0/phiData.scale
  del phiData #no longer needed in this function, so delete this reference
  #do the biggest disjoint set of tris/points stuff
  allPoints, allTris, cavPoints, cavTris = cavity.assumeNoCavities( \
      tstD.dict['POINT_XYZ'], tstD.dict['TRIANGLE_POINT'], \
      tstD.dict['POINT_NEIGHBOR'])
  #here is where code is mesh-specific
  print "setting up mesh data structure"
  meshData = mesh.mesh(gridD, tstD.dict['POINT_XYZ'], \
                  tstD.dict['POINT_NEIGHBOR'], \
                  gridSize, -2, 0, "X") #no between
  print "calculating burial depth"
  meshData.calculateTravelDistance("travelin", [3], [1])
  gridTravelInDepth = meshData.getGridTravelDistance(gridD, "travelin")
  #tstdebug.debugGridCountVals(gridTravelInDepth)
  print "writing phi file output"
  phiDataOut = phi()
  phiDataOut.createFromGrid(gridTravelInDepth,gridSize, \
                            toplabel="travel depth surf-in")
  phiDataOut.write(tstFileName+".mesh.travel.in.phi")
  print "writing pdb file output"
  pdbD = pdbData()
  for line in tstD.dict['PDB_RECORD']:
    pdbD.processLine(line)
  atomTravelInDepths = grid.assignAtomDepths(gridTravelInDepth, gridSize, \
                                        mins, maxs, pdbD)
  #make a pdb file with the bfactor replaced
  for index,atomTID in enumerate(atomTravelInDepths):
    pdbD.updateFactors(index, (pdbD.factors[index][0], atomTID))
  pdbD.write(tstFileName+".mesh.travelin.pdb")
  #also add record to tstD
  atomTIDRecord = []
  for index,atomTID in enumerate(atomTravelInDepths):
    atomTIDRecord.append([index+1, atomTID])
  print "updating tst file"
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
  print "burial depth done"

def repairPointPdbRecord(tstD=None, tstFileName=False):
  '''checks and repairs pointpdbrecord if it has no data in it'''
  same = True
  if tstD is None: #hasn't been read in already
    tstD = tstdata.tstData(tstFileName, \
              necessaryKeys=tstdata.tstData.necessaryKeysForPocket)
  lastPdbNum = tstD.dict['POINT_PDB_RECORD'][0][1]-1
  for pointPdb in tstD.dict['POINT_PDB_RECORD']:
    pdbNum = pointPdb[1]-1
    if pdbNum != lastPdbNum:
      same = False
  if same: #needs repaired, otherwise this is over
    pdbD = pdbData()
    for line in tstD.dict['PDB_RECORD']:
      pdbD.processLine(line)
    ptCoords = {}
    for pointXyz in tstD.dict['POINT_XYZ']:
      ptCoords[pointXyz[0]] = tuple(pointXyz[1:])
    coordToNearbyAtoms = pdbD.getNearbyAtoms(ptCoords.values())
    newPointPdbRec = []
    for pointPdb in tstD.dict['POINT_PDB_RECORD']:
      atomFound = coordToNearbyAtoms[ptCoords[pointPdb[0]]][0]
      newPointPdbRec.append([pointPdb[0], atomFound])
    #replace old record with new record
    tstD.dict['POINT_PDB_RECORD'] = newPointPdbRec
    if tstFileName:
      tstFile = open(tstFileName, 'a') #append into file
      tstdata.writeEntryIntegers(tstD.dict['POINT_PDB_RECORD'], \
        "POINT_PDB_RECORD", \
        "END POINT_PDB_RECORD", tstFile)
      tstFile.close()
  #done now, possibly replaced record, just return

def calculateCharges(tstD, chargeD):
  '''actually does charge assignment, returns 2 lists'''
  pdbD = pdbData()
  for line in tstD.dict['PDB_RECORD']:
    pdbD.processLine(line)
  pdbD.assignCharges(chargeD)
  #for index in xrange(len(pdbD.atoms)): #used to debug assignments
  #  print pdbD.atoms[index], pdbD.resNames[index], pdbD.charges[index]
  chargeXyz, hydroXyz = [],[] #new tst record list of [number, charge]
  for pointPdb in tstD.dict['POINT_PDB_RECORD']:
    pdbNum = pointPdb[1]-1
    tempCharge = pdbD.charges[pdbNum]
    tempHydroCharge = pdbD.hydroCharges[pdbNum]
    if tempCharge is None: #warn user that charges didn't get assigned
      print "warning: charge is not assigned for " + pdbD.atoms[pdbNum] + \
          " " + pdbD.resNames[pdbNum]
      print "charge is set to ZERO for now"
      tempCharge = 0.
      tempHydroCharge = 0.
    chargeXyz.append([pointPdb[0], tempCharge])
    hydroXyz.append([pointPdb[0], tempHydroCharge])
  return chargeXyz, hydroXyz

def tstAssignCharges(tstFileName, chargeD, appendTstFile=False):
  '''reads the tstfile, assigns charges based on the charge data, appends tst'''
  tstD = tstdata.tstData(tstFileName) #read the file into the data structure
  chargeXyz, hydroXyz = calculateCharges(tstD, chargeD)
  tstD.dict['CHARGE_XYZ'] = chargeXyz
  tstD.dict['HYDROPHOBIC_XYZ'] = hydroXyz
  if appendTstFile:
    tstFile = open(tstFileName, 'a')
    tstFile.write("CHARGE_XYZ\n")
    for line in tstD.dict['CHARGE_XYZ']:
      lineOut = "%8d" % line[0]
      for count in xrange(1,len(line)):
        lineOut += "%+9.4f " % line[count]
      noPlusLine = string.replace(lineOut, "+", " ")
      tstFile.write(noPlusLine)
      tstFile.write("\n")
    tstFile.write("END CHARGE_XYZ\n")
    tstFile.write("HYDROPHOBIC_XYZ\n")
    for line in tstD.dict['HYDROPHOBIC_XYZ']:
      lineOut = "%8d" % line[0]
      for count in xrange(1,len(line)):
        lineOut += "%+9.4f " % line[count]
      noPlusLine = string.replace(lineOut, "+", " ")
      tstFile.write(noPlusLine)
      tstFile.write("\n")
    tstFile.write("END HYDROPHOBIC_XYZ\n")
    tstFile.close()
  return tstD, chargeXyz, hydroXyz #in case caller just wanted that

def repairNearby(tstFileName):
  '''repairs the POINT_PDB_RECORD if it needs it'''
  if tstD is None: #hasn't been read in already
    tstD = tstdata.tstData(tstFileName, \
              necessaryKeys=tstdata.tstData.necessaryKeysForPocket)
  repairPointPdbRecord(tstD, tstFileName)
  #that's it

def tstPocketMap(tstFileName, phiFileName, tstD=None, \
                 ligandFileName=None, nearbyDistance=0., appendTst=True, \
                 doPCA=True):
  '''pocket mapping algorithm, finds all pockets on entire surface, puts in
  tree and graph data structure, various outputs'''
  print "read tst file"
  if tstD is None: #hasn't been read in already
    tstD = tstdata.tstData(tstFileName, \
              necessaryKeys=tstdata.tstData.necessaryKeysForPocket)
  print "repairing nearby points if necessary"
  repairPointPdbRecord(tstD, tstFileName)
  print "calculating charges"
  chargeXyz, hydroXyz = calculateCharges(tstD, charge.charge())
  print "calculating curvatures"
  edgeCurv, ptCurv, ptWeighCurv = tstCurvature.tstEdgeCurvature( \
                     tstD.dict['TRIANGLE_POINT'], tstD.dict['POINT_XYZ'], \
                     tstD.dict['POINT_TRIANGLE'], tstD.dict['POINT_NEIGHBOR'])
  tstD.dict['POINT_CURVATURE_EDGE'] = ptWeighCurv
  tstD.dict['CHARGE_XYZ'] = chargeXyz
  tstD.dict['HYDROPHOBIC_XYZ'] = hydroXyz
  print "read in phi file"
  phiData = phi(phiFileName)  #read in the phimap
  print "making mesh data structure"
  meshData, gridData = meshConstruct(tstD, phiData, tstFileName, \
                                     threshold="auto", cavities=True)
  meshData.setPtHydro(tstD.dict['HYDROPHOBIC_XYZ'])
  meshData.setPtCurvature(tstD.dict['POINT_CURVATURE_EDGE'])
  gridSize = 1.0/phiData.scale
  tstPdbRecord = tstD.dict['PDB_RECORD']
  meshData.setSurfaceArea(tstD.dict['TRIANGLE_POINT'])
  del tstD,phiData,gridData #not needed, reclaim memory
  pdbD = pdbData()
  pdbD.processLines(tstPdbRecord)
  pointAtomList = meshData.calculateNearbyAtoms(pdbD, nearbyDistance)
  meshData.setVolume(gridSize)
  print "calculating travel depth"
  meshData.calculateTravelDistance("traveldepth", [0], [2,3,5])
  pointTravelDepth = meshData.getSurfaceTravelDistance("traveldepth")
  if ligandFileName is not None: #if there is a ligand, read it
    ligand = pdbData(ligandFileName)
    ligandXYZR = ligand.getHeavyAtomXYZRadius()
    nodeWithinSet = meshData.getWithinNodesNoInside(ligandXYZR)
    bestIU = 0. #intersection / union, 1 is perfect
    #print nodeWithinSet, len(nodeWithinSet)
  print "pocket mapping starting"
  if ligandFileName is not None and len(nodeWithinSet) > 0:
    outFileName = ligandFileName
    #tstdebug.nodeDebug(nodeWithinSet, \
    #              filename = tstFileName+".within.ligand.py")
    localMaxima, borders, tm3tree, surfNodeToLeaf = meshData.pocketMapping( \
           'traveldepth',  \
           [2,3,5], pointAtomList, pdbD, \
           outName=outFileName+".", groupName='group', \
           ligandNodes=nodeWithinSet, doPCA=doPCA)
  else:
    outFileName = tstFileName
    localMaxima, borders, tm3tree, surfNodeToLeaf = meshData.pocketMapping( \
           'traveldepth',  \
           [2,3,5], pointAtomList, pdbD, \
           outName=outFileName+".", groupName='group', doPCA=doPCA)
  #print len(localMaxima), len(borders), tm3tree, len(surfNodeToLeaf)
  #tstdebug.nodeDebug(localMaxima, \
  #              filename=tstFileName+".localmaxima.pocketmap.py")
  #tstdebug.nodeDebug(borders, \
  #              filename=tstFileName+".borders.pocketmap.py")
  #tstdebug.nodeDebug(meshData.getSurfaceNodes(), \
  #              filename=tstFileName+".groups.pocketmap.py", name='group')
  tm3tree.write(outFileName+".tree.tm3")
  #tm3tree.writeTNV(tstFileName+".tree.tnv") #doesn't seem to import into treemap correctly
  if appendTst: #turn off sometimes since appends to tst file
    print "appending data to tst file"
    surfNodes = meshData.getSurfaceNodes()
    pointLeafList = []
    for aNode in surfNodes:
      if aNode not in surfNodeToLeaf:
        print aNode, aNode.distances
        leafNum = 0 #made up and wrong for testing
      else:
        leafNum = surfNodeToLeaf[aNode]
      pointLeafList.append([aNode, int(leafNum)])
    #print pointLeafList
    leafToGroup = tm3tree.getLeafToGroup()
    leafGroupList = []
    leafKeyMax = max(leafToGroup.keys())
    for leaf in xrange(leafKeyMax):
      tempList = [leaf + 1]
      try:
        tempList.extend(leafToGroup[leaf + 1])
      except KeyError:
        pass #means leaf doesn't exist
      leafGroupList.append(tempList)
    #print leafGroupList
    tstFile = open(tstFileName, 'a')
    tstdata.writeEntryIntegers(pointLeafList, "POINT_LEAF LIST", \
                            "END POINT_LEAF", tstFile)
    tstdata.writeEntryIntegers(leafGroupList, "LEAF_GROUP LIST", \
                            "END LEAF_GROUP", tstFile)
    tstdata.writeEntryIntegers(pointAtomList, "POINT_NEARBY_ATOM LIST", \
                            "END POINT_NEARBY_ATOM", tstFile)
    #also write curvature and charge data here
    tstdata.writeEntrySingleFloat(ptWeighCurv, "POINT_CURVATURE_EDGE LIST", \
                            "END POINT_CURVATURE_EDGE", tstFile)
    tstdata.writeEntrySingleFloat(chargeXyz, "CHARGE_XYZ", \
                            "END CHARGE_XYZ", tstFile)
    tstdata.writeEntrySingleFloat(hydroXyz, "HYDROPHOBIC_XYZ", \
                            "END HYDROPHOBIC_XYZ", tstFile)
    #write data to file
    tstFile.write("DEPTH_TRAVEL_DIST\n")
    for line in pointTravelDepth:
      lineOut = "%8d" % line[0]
      for count in xrange(1,len(line)):
        lineOut += "%+9.4f " % line[count]
      noPlusLine = string.replace(lineOut, "+", " ")
      tstFile.write(noPlusLine)
      tstFile.write("\n")
    tstFile.write("END DEPTH_TRAVEL_DIST\n")
    tstFile.close()
  print "pocket mapping complete"

def printHelpMessage():
  '''prints the usage requirements for this script'''
  print "Usage: tstTravelDist.py meshdepth tstFile phiFile [ligandPdbFile]"
  print "Usage: tstTravelDist.py meshdepthcav tstFile phiFile [ligandPdbFile]"
  print "Usage: tstTravelDist.py meshsurfout tstFile phiFile"
  print "Usage: tstTravelDist.py meshsurfin tstFile phiFile"
  print "Usage: tstTravelDist.py cavityremove tstFile tstFileOut phiFile phiFileOut"
  print "Usage: tstTravelDist.py countholes tstFile"
  print "Usage: tstTravelDist.py findholes tstFile phiFile"
  print "Usage: tstTravelDist.py findholesdebug tstFile phiFile"
  print "Usage: tstTravelDist.py assigncharge tstFile chargeFile"
  print "Usage: tstTravelDist.py repairnearby tstFile [more tstFiles]"
  print "Usage: tstTravelDist.py pocketmap tstFile phiFile [ligandFile]"

#this is where main is... maybe add some other arguments like gridSize?
if -1 != string.find(sys.argv[0], "tstTravelDist.py"):
  if 1 < len(sys.argv):
    if sys.argv[1] == "meshdepth" and  len(sys.argv) > 4:
      tstFile, phiFile, ligandPdbFile = sys.argv[2:5]
      print "tstname phiname ligandname ",
      print "within-min max mean traceback-min max mean ",
      print "between-min max mean surface-min max mean",
      print "within-volume trace-volume" #actually want newline
      print tstFile, phiFile, ligandPdbFile, #wait for more output
      tstTravelDepthMesh(tstFile, phiFile, ligandPdbFile)
    elif sys.argv[1] == "meshdepth" and  len(sys.argv) > 3:
      tstFile, phiFile = sys.argv[2:4]
      print tstFile, phiFile
      tstTravelDepthMesh(tstFile, phiFile)
    elif sys.argv[1] == "repairnearby" and  len(sys.argv) > 2:
      for tstFile in sys.argv[2:]:
        print tstFile
        tstTravelRepairNearby(tstFile)
    elif sys.argv[1] == "meshdepthcav" and len(sys.argv) > 4:
      tstFile, phiFile, ligandPdbFile = sys.argv[2:5]
      print "tstname phiname ligandname ",
      print "within-min max mean traceback-min max mean ",
      print "between-min max mean surface-min max mean",
      print "within-volume trace-volume" #actually want newline
      print tstFile, phiFile, ligandPdbFile, #wait for more output
      tstTravelDepthMesh(tstFile, phiFile, ligandPdbFile, \
                         cavities=True, threshold="auto")
    elif sys.argv[1] == "meshdepthcav" and  len(sys.argv) > 3:
      tstFile, phiFile = sys.argv[2:4]
      print tstFile, phiFile
      tstTravelDepthMesh(tstFile, phiFile, cavities=True, threshold="auto")
    elif sys.argv[1] == "meshsurfout" and  len(sys.argv) > 3:
      tstFile, phiFile = sys.argv[2:4]
      print tstFile, phiFile
      tstTravelSurfOutsideMesh(tstFile, phiFile)
    elif sys.argv[1] == "findholes" and  len(sys.argv) > 3:
      tstFile, phiFile = sys.argv[2:4]
      print tstFile, phiFile
      tstTravelFindHoles(tstFile, phiFile, debugOut=False)
    elif sys.argv[1] == "findholesdebug" and  len(sys.argv) > 3:
      tstFile, phiFile = sys.argv[2:4]
      print tstFile, phiFile
      tstTravelFindHoles(tstFile, phiFile, debugOut=True)
    elif sys.argv[1] == "meshsurfin" and  len(sys.argv) > 3:
      tstFile, phiFile = sys.argv[2:4]
      print tstFile, phiFile
      tstTravelSurfInsideMesh(tstFile, phiFile)
    elif sys.argv[1] == "cavityremove" and  len(sys.argv) > 3:
      tstFile, tstOut, phiFile, phiOut = sys.argv[2:6]
      print tstFile, tstOut, phiFile, phiOut
      cavity.tstCavityRemoval(tstFile, tstOut, phiFile, phiOut)
    elif sys.argv[1] == "countholes" and  len(sys.argv) > 2:
      tstFileIn = sys.argv[2]
      print tstFileIn,
      print tstCountHoles(tstFileIn)
    elif sys.argv[1] == "assigncharge" and len(sys.argv)> 2:
      tstFile = sys.argv[2]
      print tstFile,
      if len(sys.argv) > 3:
        chargeFile = sys.argv[3]
        print chargeFile
        chargeD = charge.charge(chargeFile)
      else:
        chargeD = charge.charge()
      tstAssignCharges(tstFile, chargeD)
    elif sys.argv[1] == "pocketmap" and len(sys.argv)> 3:
      tstFile, phiFile = sys.argv[2:4]
      if len(sys.argv) > 4:
        ligandPdbFile = sys.argv[4]
        print tstFile, phiFile, ligandPdbFile
        tstPocketMap(tstFile, phiFile, ligandFileName=ligandPdbFile)
      else:
        print tstFile, phiFile
        tstPocketMap(tstFile, phiFile)
    else:
      printHelpMessage()
  else:
    printHelpMessage()
