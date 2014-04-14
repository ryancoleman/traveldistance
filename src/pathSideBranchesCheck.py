#!/usr/bin/env python

#idea is to find paths that exit into the membrane barrier, and figure out
#what residues line them and what their radius is.
#this version only does the branch from a 'proper' path to the membrane interior
#and removes other points along the path and other residues.

import string, sys, comparePaths, os, pdb, geometry, tstdebug
import paths as pathsModule
import addFoundHoleStats

def countCrossingsZ(path, barrier):
  crossings = 0
  thisPtZ = path[0][2]
  for nextPt in path[1:]:
    nextPtZ = nextPt[2]
    if thisPtZ <= barrier and barrier <= nextPtZ:
      crossings += 1
    elif thisPtZ >= barrier and barrier >= nextPtZ:
      crossings += 1
    thisPtZ = nextPtZ #setup for next iteration
  return crossings

def checkPathBarriers(prefix):
  tstName = prefix + ".nocav.tst"
  findHolesName = tstName + ".findholes.log"
  findHolesFile = open(findHolesName, 'r')
  findHolesLines = findHolesFile.readlines()
  findHolesFile.close()
  HolesName = tstName + ".sideshole.log" #holds all the output
  goodHolesName = tstName + ".good.sideshole.log" #just the 1 1 0 1 1
  sideHolesName = tstName + ".side.sideshole.log" #just the * * 1 * *
  badHolesName = tstName + ".bad.sideshole.log" #all others
  pdbWithBarriersFileName = "planes_" + prefix + ".pdb"
  pdbBarriers = pdb.pdbData(pdbWithBarriersFileName)
  #get the barriers read in and defined
  barrierAtomList = [[],[]]
  for index,resName in enumerate(pdbBarriers.resNames):
    if resName == "DUM":
      if pdbBarriers.atoms[index][0] == "O":
        barrierAtomList[0].append(pdbBarriers.coords[index])
      elif pdbBarriers.atoms[index][0] == "N":
        barrierAtomList[1].append(pdbBarriers.coords[index])
  barrierZ = [barrierAtomList[0][0][2], barrierAtomList[1][0][2]]
  barrierZ.sort()
  barrierSep = geometry.distL2(barrierAtomList[0][0],  \
                             barrierAtomList[1][0])
  #barrier is just Z coordinate
  #setup for main loop over paths
  poreSuffix = ".pore.py"
  logFile = open(HolesName, 'w')
  goodLogFile = open(goodHolesName, 'w')
  sideLogFile = open(sideHolesName, 'w')
  badLogFile = open(badHolesName, 'w')
  #the following 5 things are calculated and written for each path, headers
  #the 6th, barrier separation, is really the same for each structure
  logFile.write("endsBeyond1count barrier1count endsBetweenCount ")
  logFile.write("barrier2count endsBeyond2count barrierSeparation\n")
  goodLogFile.write("prefix ")
  goodLogFile.write(string.strip(findHolesLines[0]) + " ")
  goodLogFile.write("endsBeyond1count barrier1count endsBetweenCount ")
  goodLogFile.write("barrier2count endsBeyond2count barrierSeparation\n")
  sideLogFile.write("prefix ")
  sideLogFile.write(string.strip(findHolesLines[0]) + " ")
  sideLogFile.write("endsBeyond1count barrier1count endsBetweenCount ")
  sideLogFile.write("barrier2count endsBeyond2count barrierSeparation\n")
  badLogFile.write("prefix ")
  badLogFile.write(string.strip(findHolesLines[0]) + " ")
  badLogFile.write("endsBeyond1count barrier1count endsBetweenCount ")
  badLogFile.write("barrier2count endsBeyond2count barrierSeparation\n")
  holeNumber = 1
  poreFile =  tstName + "." + str(holeNumber) + poreSuffix
  print poreFile
  paths = []
  sides,goods = [],[]
  endsToPaths = {}
  pathsToEnds = {}
  while os.path.exists(poreFile):
    path = comparePaths.readCGOPath(poreFile)
    pathRad = comparePaths.readCGOPathWithRadius(poreFile)
    paths.append(pathRad)
    pathNum = len(paths) - 1
    for end in string.split(findHolesLines[holeNumber])[1:3]:
      if pathNum not in pathsToEnds:
        pathsToEnds[pathNum] = []
      pathsToEnds[pathNum].append(end)
      if end not in endsToPaths:
        endsToPaths[end] = []
      endsToPaths[end].append(pathNum)
    intersections = [0,0]
    for index, barrier in enumerate(barrierZ):
      intersections[index] = countCrossingsZ(path, barrier)
    ends = [0,0,0]
    for endPoint in [path[0],path[-1]]:
      endPointZ = endPoint[2]
      if endPointZ < barrierZ[0] and endPointZ < barrierZ[1]:
        ends[0] += 1
      elif endPointZ >= barrierZ[0] and endPointZ <= barrierZ[1]:
        ends[1] += 1
      elif endPointZ > barrierZ[0] and endPointZ > barrierZ[1]:
        ends[2] += 1
    outputThisTime = str(ends[0]) + " " + str(intersections[0]) + " " + \
                     str(ends[1]) + " " + str(intersections[1]) + " " + \
                     str(ends[2]) + " " + str(barrierSep) + " "
    logFile.write(outputThisTime)
    logFile.write("\n")
    if ends[0] + ends[1] + ends[2] != 2:
      print "problems sorting out the ends"
    if ends[0] == 1 and ends[2] == 1 and intersections == [1,1]: #it is 'good'
      goods.append(pathNum)
      goodLogFile.write(prefix + " ")
      goodLogFile.write(string.strip(findHolesLines[holeNumber]) + " ")
      goodLogFile.write(outputThisTime + "\n")
    elif ends[1] == 1:
      sides.append(pathNum)
      sideLogFile.write(prefix + " ")
      sideLogFile.write(string.strip(findHolesLines[holeNumber]) + " ")
      sideLogFile.write(outputThisTime + "\n")
    else:
      badLogFile.write(prefix + " ")
      badLogFile.write(string.strip(findHolesLines[holeNumber]) + " ")
      badLogFile.write(outputThisTime + "\n")
    #and that is it for this path
    holeNumber += 1     #get set up for next pass
    poreFile =  tstName + "." + str(holeNumber) + poreSuffix
  logFile.close()
  goodLogFile.close()
  sideLogFile.close()
  badLogFile.close()
  #next lines are for debugging the new data structures
  '''
  print sides
  print goods
  print endsToPaths
  print pathsToEnds
  '''
  #now want to find side branches of good paths
  branches = 0
  branchSuffix = ".branch.py"
  branchFile =   tstName + "." + str(branches) + branchSuffix
  branchLog = open(tstName + ".branchholes.log", 'w')
  branchLog.write(string.strip(findHolesLines[0]) + "\n")
  for side in sides:
    foundGoods = []
    for sideEnd in pathsToEnds[side]:
      for good in goods:
        for goodEnd in pathsToEnds[good]:
          if goodEnd == sideEnd:
            foundGoods.append(good)
    if len(foundGoods) > 0:
      branchedPath = paths[side] #start with whole path
      for good in foundGoods: #remove physiological intersecting paths
        branchedPath = pathsModule.subtractPaths(branchedPath, paths[good])
      if len(branchedPath) > 0: #has to have some length remaining
        branches += 1
        branchFile = tstName + "." + str(branches) + branchSuffix
        print branches, side, foundGoods
        tstdebug.debugSetGridSpheres(branchedPath,0.5,branchFile,radius=True,mainColor=(0.01,0.9,0.35))
        branchLog.write(str(branches) + " ")
        branchLog.write(str(pathsToEnds[side][0]) + " ")
        branchLog.write(str(pathsToEnds[side][1]) + " ")
        branchLog.write("- ") #dummy, not real
        branchLog.write("0. 0. 0. 0. 0. 0. 0. 0. 0. 0. \n")
  branchLog.close()
  addFoundHoleStats.redoFindholes(prefix, nearbyDistance=4., \
          logExt=".branchholes.log", poreSuffix = ".branch.py", \
          nearbyName=".branch")

if -1 != string.find(sys.argv[0], "pathSideBranchesCheck"):
  if len(sys.argv) >= 2:
    prefix = sys.argv[1]
    prefix = prefix.replace(".nocav.tst.findholes.log", "")#just in case
    checkPathBarriers(prefix)
  else:
    print "pathSidesCheck.py prefix[.nocav.tst.findholes.log]"
