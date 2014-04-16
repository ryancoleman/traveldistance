#!/usr/bin/env python

#ryan g coleman ryangc@mail.med.upenn.edu
#takes a pdb file and makes a tst file, possibly add later improvements
#to do common tasks automatically.

import sys
import os
import string  # necessary for file and string stuff
import pdb
#import trigenTst #alternate python version of fortran, no sizelimits but slower

#names of executables
#trisrf = {65: "trisrf65", 129: "trisrf129", 193: "trisrf193"}
trisrf = {193: "trisrflarge193"}
#meshsrf = {65: "meshsrfA_65", 129: "meshsrfA_129", 193: "meshsrfA_193"}
meshsrf = {193: "meshsrflarge193"}
#gen="trigen"
gen = "trigenlarge"
#probe radius, set in fortran code
meshprobe = 1.2  # set extra small because of ions in tunnels
triprobe = 1.8
contour = 6.0
radscale = 1.0

def makeTst(
    pdbFileName, gridSpacing, pathTo="$TDHOME/bin/", whichSrf="tri",
    probeSize=False, radScaleIn=False):
  '''
  figures out how to create proper grid spacing, calls fortran programs
  '''
  pathTo = os.path.expandvars(pathTo)
  if "tri" == whichSrf:  # pick which method to use, default
    srf = trisrf
    probe = triprobe
  elif "mesh" == whichSrf:  # alternate, better for tunnels
    srf = meshsrf
    probe = meshprobe
  if probeSize:  # means use non-default value...
    probe = probeSize
  if radScaleIn:  # same, use non-default value
    radScaleUse = float(radScaleIn)
  else:
    radScaleUse = radscale
  pdbEntry = pdb.pdbData(pdbFileName)
  rootNameTemp = string.replace(pdbFileName, ".pdb", "")
  rootName = string.replace(rootNameTemp, ".PDB", "")
  tstFileName = rootName + ".tst"
  triFileName = rootName + ".tri"
  phiFileName = rootName + ".phi"
  mins, maxs = [xVal for xVal in pdbEntry.coords[0]], [
      xVal for xVal in pdbEntry.coords[0]]
  atoms = pdbEntry.atoms
  for dimension in range(3):
    for index, coord in enumerate(pdbEntry.coords):
      mins[dimension] = min(
          mins[dimension], coord[dimension] -
          (pdb.radiiDefault[atoms[index][0]]*radScaleUse))
      maxs[dimension] = max(
          maxs[dimension], coord[dimension] +
          (pdb.radiiDefault[atoms[index][0]]*radScaleUse))
  difference = 0
  for dimension in range(3):
    difference = max(difference, maxs[dimension] - mins[dimension])
  length = difference + 2. * probe
  gridScale = 1. / float(gridSpacing)
  possibleGridSizes = srf.keys()
  possibleGridSizes.sort()
  for possibleGridSize in possibleGridSizes:
    percentFill = gridScale * 100 * length / (possibleGridSize - 1)
    if percentFill < 99:
      break   # keep these settings, they are good enough
    if possibleGridSizes[-1] == possibleGridSize:
      print "no grid size big enough, either make new version of tst or " + \
          "adjust grid size parameter"
      sys.exit(1)
  srfExecutable = pathTo + srf[possibleGridSize]
  if not os.path.exists(srfExecutable):
    print "the surface preparation executable does not exist at: ", \
        srfExecutable
    exit(1)
  if "tri" == whichSrf:
    execString = srfExecutable + "  " + \
        pdbFileName + " " + str(contour)   # run trisrf
  elif "mesh" == whichSrf:
    execString = srfExecutable + "  " + \
        pdbFileName + " " + str(probe) + " " + str(radScaleUse)  # run meshsrf
  #print percentFill, execString
  try:
    os.unlink("trisrf.tri")
  except OSError:
    pass  # this is okay, just making sure it is deleted
  trisrfProc = os.popen4(execString)
  if "tri" == whichSrf:  # pick which method to use, default
    trisrfProc[0].write(str(percentFill) + "\n33\n")
  elif "mesh" == whichSrf:  # alternate, better for tunnels
    trisrfProc[0].write(str(percentFill) + "\n")
  trisrfProc[0].flush()
  trisrfProc[0].close()
  finishedRunningSrf = trisrfProc[1].read()
  log = open(tstFileName + ".log", 'w')
  log.write(finishedRunningSrf)
  if "tri" == whichSrf:  # pick which method to use, default
    try:
      os.rename("trisrf.tri", triFileName)
      os.rename("trisrf.phi", phiFileName)
    except OSError:  # actual problem
      print "trisrf did not make .tri file, check logs"
      log.close()
      return False
  elif "mesh" == whichSrf:  # alternate, better for tunnels
    try:
      os.rename("meshsrfA.tri", triFileName)
      os.rename("meshsrfA.phi", phiFileName)
    except OSError:  # actual problem
      print "meshsrf did not make .tri file, check logs"
      log.close()
      return False
  if not os.path.exists(pathTo + gen):
    print "the surface generation executable does not exist at: " + pathTo + \
        gen
    exit(1)
  trigenProc = os.popen4(pathTo + gen + " " + triFileName + " " + tstFileName)
  trigenProc[0].flush()
  trigenProc[0].close()
  finishedRunningGen = trigenProc[1].read()
  try:
    os.unlink("trilinel.dat")
    os.unlink("triline.usr")
    os.unlink("trinext.dat")
    os.unlink("trisrf.pdb")
    os.unlink("trisrf.rec")
    os.unlink("trisrf.usr")
    os.unlink("fort.10")
    os.unlink("trigen.py")
    os.unlink("triline.py")
  except OSError:
    pass  # again, just cleaning up junk files
  try:
    os.unlink("mesh.pdb")
    os.unlink("meshline.usr")
    os.unlink("meshlinel.dat")
    os.unlink("meshtri.dat")
    os.unlink("fort.10")
    os.unlink("trigen.py")
    os.unlink("triline.py")
  except OSError:
    pass  # again, just cleaning up junk files
  log.write(finishedRunningGen)
  log.close()
  return True  # indicates success

#this is the main loop here, usually called from tstTravelDepth now
if -1 != string.find(sys.argv[0], "tstCreate"):
  if len(sys.argv) >= 3:
    pdbFileName = sys.argv[1]
    gridSpacing = sys.argv[2]
    if len(sys.argv) > 3:
      whichSurfaceCode = sys.argv[3]
      if len(sys.argv) > 4:
        pathTo = sys.argv[4]
        makeTst(pdbFileName, gridSpacing, pathTo, whichSrf=whichSurfaceCode)
      else:
        makeTst(pdbFileName, gridSpacing, whichSrf=whichSurfaceCode)
    else:
      makeTst(pdbFileName, gridSpacing)
  else:
    print "Usage: tstCreate.py pdbFile gridSpacing " + \
        "[tri|mesh] [pathToExecutables]"
