#!/usr/bin/env python2.5
#ryan g. coleman ryangc@mail.med.upenn.edu
#does everything to do travel depth, starting from a pdb file and optionally a 
#gridsize and optionally a path to the trisrf/trigen execs

import sys,string,os

import tstCreate
import tstConvexHull
import oldTravelDist
import tstdata
import cavity
import tstTravelDist #for repair of nearby point->pdb mapping

def makeTstDoChecks(pdbFileName, gridSize, pathToExecs, rootName, \
                    whichSrf, probeSize, radScaleIn, perturb=.0041, maxTry=50,
                    check2trisPerEdge=True):
  '''checks to make sure tst building goes okay... and then tris per point <=12
   modifies gridSize slightly if failure... usually enough'''
  tstFileName = rootName + ".tst"
  everythingOkay = False
  gridSizeTry = float(gridSize)
  while not everythingOkay and (-1 == maxTry or maxTry > 0):
    if maxTry != -1:
      maxTry -= 1
    print "trying again", maxTry, pdbFileName,
    attempt = tstCreate.makeTst(pdbFileName, gridSizeTry, pathToExecs, \
                                whichSrf, probeSize, radScaleIn)
    print "attempt", gridSizeTry
    if not attempt:
      gridSizeTry -= perturb
      everythingOkay = False
    else:
      try:
        tstD = tstdata.tstData(tstFileName)
        everythingOkay = True
        #this checks for cases where 13 or more points are shared by a tri. tst code has problems outputting this
        for pointTriRec in tstD.dict['POINT_TRIANGLE']:
          if everythingOkay and pointTriRec[1] >= 13:
            everythingOkay = False
            print "more than 13 triangles at one vertex", pointTriRec[1]
        #makes sure each edge is in only 2 triangles!
        for pointTriRec in tstD.dict['POINT_NEIGHBOR']:
          if everythingOkay: #no need to keep checking otherwise
            onePt = pointTriRec[0]
            theseTris = tstD.dict['POINT_TRIANGLE'][onePt-1][2:]
            for otherPt in pointTriRec[2:]:
              if otherPt > onePt: #don't need to check each edge twice, just once
                triCount = 0
                for otherTri in tstD.dict['POINT_TRIANGLE'][otherPt-1][2:]:
                  try:
                    theseTris.index(otherTri) 
                    triCount += 1
                  except ValueError:
                    pass
                if check2trisPerEdge and triCount != 2:
                  everythingOkay = False #this is bad...
                  print "more than 2 triangles per edge...", triCount
        if not everythingOkay:
          gridSizeTry -= perturb
      except IOError: #file doesn't exist, usually a problem with fortran, try again
        everythingOkay = False
        print "ioerror, usually fortran problem, check ....tst.log"
        gridSizeTry -= perturb
  if 0 == maxTry:
    print pdbFileName + " exceeded number of attempts"
    return False
  return True #indicates success
  #done, files created, all that needs to happen, maybe add max iteration count

#default parameters, change for your system
qhullExec = "qhull" #standard code
khullExec = "khull" #kim's convex hull code
def runTravelDepthCompletely(pdbFileName, gridSize=1.0, \
                             whichSrf="mesh", probeSize=False, \
                             radScaleIn=False, \
                             pathToExecs="$TDHOME/bin/", \
                             deleteFiles=True, check2trisPerEdge=True, \
                             doTravelDepth=True):
  pathToExecs = os.path.expandvars(pathToExecs)
  rootNameTemp = string.replace(pdbFileName,".pdb","")
  rootName = string.replace(rootNameTemp,".PDB","")
  #now there are 4 more files
  print "running the tst creation code"
  okay = makeTstDoChecks(pdbFileName, gridSize, pathToExecs, rootName, \
                        whichSrf, probeSize, radScaleIn, \
                        check2trisPerEdge=check2trisPerEdge)
  if okay:
    tstFileName = rootName + ".tst"
    print "repairing nearby points if necessary"
    tstTravelDist.repairPointPdbRecord(tstD=None, tstFileName=tstFileName)
    phiFileName = rootName + ".phi"
    tstCavFileName = rootName + ".cav.tst"
    phiCavFileName = rootName + ".cav.phi"
    nocavTstFileName = rootName + ".nocav.tst"
    nocavPhiFileName = rootName + ".nocav.phi"
    print "removing cavities"
    cavity.tstCavityRemoval(tstFileName, nocavTstFileName, \
                                   phiFileName, nocavPhiFileName)
    #now no more cavities
    print "running convex hull"
    tstConvexHull.tstConvexHull(pathToExecs + khullExec, nocavTstFileName, \
                                pathToExecs + qhullExec)
    convexTempFile1 = nocavTstFileName + ".tempQHullFile.output"
    convexTempFile2 = nocavTstFileName + ".tempQHullFile"
    #now has convex hull data,next step is travel depth
    if doTravelDepth: #sometimes just skip for pockets, etc
      print "running travel depth"
      oldTravelDist.tstTravelDepth(nocavTstFileName, nocavPhiFileName)
    #now has travel depth data in it
    print "renaming files, deleting temporary files"
    os.rename(tstFileName, tstCavFileName) #keep the cavities
    os.rename(phiFileName, phiCavFileName)
    if deleteFiles: #delete unnecessary debugging files...
      os.unlink(convexTempFile1)
      os.unlink(convexTempFile2)
      os.unlink(rootName + ".tri")

#this is main
if -1 != string.find(sys.argv[0], "tstTravelDepth.py"):
  if 6 < len(sys.argv):
    pdbFileName = sys.argv[1]
    gridSize = sys.argv[2]
    whichSurface = sys.argv[3]
    probe = float(sys.argv[4])
    radscale = sys.argv[5]
    path = sys.argv[6]
    runTravelDepthCompletely(pdbFileName, gridSize, whichSurface, \
                             probe, radscale, path)
  elif 5 < len(sys.argv):
    pdbFileName = sys.argv[1]
    gridSize = sys.argv[2]
    whichSurface = sys.argv[3]
    probe = float(sys.argv[4])
    radscale = sys.argv[5]
    runTravelDepthCompletely(pdbFileName, gridSize, whichSurface, \
                             probe, radscale)
  elif 4 < len(sys.argv):
    pdbFileName = sys.argv[1]
    gridSize = sys.argv[2]
    whichSurface = sys.argv[3]
    probe = float(sys.argv[4])
    runTravelDepthCompletely(pdbFileName, gridSize, whichSurface, probe)
  elif 3 < len(sys.argv):
    pdbFileName = sys.argv[1]
    gridSize = sys.argv[2]
    whichSurface = sys.argv[3]
    runTravelDepthCompletely(pdbFileName, gridSize, whichSurface)
  elif 2 < len(sys.argv):
    pdbFileName = sys.argv[1]
    gridSize = sys.argv[2]
    runTravelDepthCompletely(pdbFileName, gridSize)
  elif 1 < len(sys.argv):
    pdbFileName = sys.argv[1]
    runTravelDepthCompletely(pdbFileName)
  else:
    print "Usage: tstTravelDepth.py file.pdb [gridsize] [tri|mesh] [probe] [radius scale] [pathToExecs]"
