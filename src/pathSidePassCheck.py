#!/usr/bin/env python

#find paths that enter and exit in membrane barrier

import string
import sys
import comparePaths
import os
import pdb
import geometry

def countCrossingsZ(path, barrier):
  crossings = 0
  thisPtZ = path[0][2]
  for nextPt in path[1:]:
    nextPtZ = nextPt[2]
    if thisPtZ <= barrier and barrier <= nextPtZ:
      crossings += 1
    elif thisPtZ >= barrier and barrier >= nextPtZ:
      crossings += 1
    thisPtZ = nextPtZ  # setup for next iteration
  return crossings

def checkPathBarriers(prefix):
  tstName = prefix + ".nocav.tst"
  findHolesName = tstName + ".findholes.log"
  findHolesFile = open(findHolesName, 'r')
  findHolesLines = findHolesFile.readlines()
  findHolesFile.close()
  HolesName = tstName + ".membranehole.log"  # holds all the output
  goodHolesName = tstName + ".good.membranehole.log"  # just the 1 1 0 1 1
  sideHolesName = tstName + ".side.membranehole.log"  # just the * * 1 * *
  badHolesName = tstName + ".bad.membranehole.log"  # all others
  pdbWithBarriersFileName = "planes_" + prefix + ".pdb"
  pdbBarriers = pdb.pdbData(pdbWithBarriersFileName)
  #get the barriers read in and defined
  barrierAtomList = [[], []]
  for index,resName in enumerate(pdbBarriers.resNames):
    if resName == "DUM":
      if pdbBarriers.atoms[index][0] == "O":
        barrierAtomList[0].append(pdbBarriers.coords[index])
      elif pdbBarriers.atoms[index][0] == "N":
        barrierAtomList[1].append(pdbBarriers.coords[index])
  barrierZ = [barrierAtomList[0][0][2], barrierAtomList[1][0][2]]
  barrierZ.sort()
  barrierSep = geometry.distL2(barrierAtomList[0][0], barrierAtomList[1][0])
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
  sides = []
  while os.path.exists(poreFile):
    path = comparePaths.readCGOPath(poreFile)
    paths.append(path)
    intersections = [0, 0]
    for index, barrier in enumerate(barrierZ):
      intersections[index] = countCrossingsZ(path, barrier)
    ends = [0, 0, 0]
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
    if ends[0] == 1 and ends[2] == 1 and intersections == [1, 1]:  
      # it is 'good'
      goodLogFile.write(prefix + " ")
      goodLogFile.write(string.strip(findHolesLines[holeNumber]) + " ")
      goodLogFile.write(outputThisTime + "\n")
    elif ends[1] == 2:
      sides.append(len(paths) - 1)
      sideLogFile.write(prefix + " ")
      sideLogFile.write(string.strip(findHolesLines[holeNumber]) + " ")
      sideLogFile.write(outputThisTime + "\n")
    else:
      badLogFile.write(prefix + " ")
      badLogFile.write(string.strip(findHolesLines[holeNumber]) + " ")
      badLogFile.write(outputThisTime + "\n")
    #and that is it for this path
    holeNumber += 1     # get set up for next pass
    poreFile =  tstName + "." + str(holeNumber) + poreSuffix
  print sides
  logFile.close()
  goodLogFile.close()
  sideLogFile.close()
  badLogFile.close()

if -1 != string.find(sys.argv[0], "pathSidePassCheck"):
  if len(sys.argv) >= 2:
    prefix = sys.argv[1]
    prefix = prefix.replace(".nocav.tst.findholes.log", "")  # just in case
    checkPathBarriers(prefix)
  else:
    print "pathSidesCheck.py prefix[.nocav.tst.findholes.log]"
