#!/usr/bin/env python

#compares 2 paths found/created... gives something like rmsd
#everything in source is matched to closest in target... then "rmsd" computed
#1.5.2007 added support for HOLE output .sph files

import pdb
import geometry
import sys
import string

def readCGOPath(filename):
  result = readCGOPathRadius(filename)
  return result[0]

def readCGOPathWithRadius(filename):
  result = readCGOPathRadius(filename)
  output = []
  for index in xrange(len(result[0])):
    output.append(list(result[0][index]))
    output[-1].insert(0, result[1][index])
  return output

def readCGOPathRadius(filename):
  raw = []
  if filename:
    file = open(filename, 'r')
    try:
      for line in file:
        raw.append(string.rstrip(line))
    except StopIteration:
      pass  # EOF
    file.close()
  data = []
  for rawLine in raw:
    for commasSplit in string.split(rawLine, ","):
      for spaceSplit in string.split(commasSplit, ' '):
        if len(spaceSplit) > 0:
          data.append(spaceSplit)
  sphereIndices = []
  for index, token in enumerate(data):
    if 'SPHERE' == token:
      sphereIndices.append(index)
  result = []
  radii = []
  for index in sphereIndices:
    result.append(
        (float(data[index + 1]), float(data[index + 2]),
            float(data[index + 3])))
    try:
      radii.append((float(data[index+4])))
    except ValueError:
      print index, data[index:index+5]
  return result, radii

def readSphPath(sphFileName):
  result = readSphPathRadius(sphFileName)
  return result[0]  # just path

def readSphPathRadius(sphFileName):
  result, radii = [], []
  sphData = pdb.pdbData(sphFileName)     # sph is basically PDB format
  seenZeroYet = False  # stupid hack for sph file
  for index, factors in enumerate(sphData.factors):  # if second factor
    resNum = sphData.resNums[index]
    if factors[1] > 0.0:                  # is non-zero, then part of path
      if resNum > 0:
        result.append(sphData.coords[index])
        radii.append(sphData.factors[index][0])
      elif resNum == 0 and not seenZeroYet:
        result.append(sphData.coords[index])
        radii.append(sphData.factors[index][0])
        seenZeroYet = True
  return result, radii

def readPointFile(pointFileName):
  raw = []
  if pointFileName:
    file = open(pointFileName, 'r')
    try:
      for line in file:
        raw.append(string.rstrip(line))
    except StopIteration:
      pass  # EOF
    file.close()
  result = []
  for line in raw:
    coords = string.split(line)
    if len(coords) == 3:
      xyz = []
      for coord in coords:
        xyz.append(float(coord))
      result.append(xyz)
  return result

def comparePaths(
    source=False, target=False, sourceDataManual=False,
    core=100., sourceRadiiManual=False):
  result = comparePathsCoverage(
      source, target, sourceDataManual, core, sourceRadiiManual)
  return result[0]  # just rmsd

def comparePathsCoverage(
    source=False, target=False,
    sourceDataManual=False, core=100., sourceRadiiManual=False):
  '''get both the rmsd and how much of the target was seen'''
  result = comparePathsManyMetrics(
      source, target, sourceDataManual, sourceRadiiManual, core)
  return result[:2]  # just first two

def comparePathsManyMetrics(
    source=False, target=False, sourceDataManual=False,
    sourceRadiiManual=False, core=100):
  '''do lots of different metrics to compare the paths'''
  #open and read in both first
  sourceData, sourceRadii, targetData, targetRadii = False, False, False, False
  if source:
    if -1 != string.find(source, "py"):
      sourceData, sourceRadii = readCGOPathRadius(source)
    elif -1 != string.find(source, "sph"):
      sourceData, sourceRadii = readSphPathRadius(source)
  else:
    sourceData = sourceDataManual
    sourceRadii = sourceRadiiManual
  if target:
    if -1 != string.find(target, "py"):
      targetData, targetRadii = readCGOPathRadius(target)
    elif -1 != string.find(target, "sph"):
      targetData, targetRadii = readSphPathRadius(target)
  #now have the data... now do the one-sided RMSD thing
  sumDistanceSquared, sumWeighted = 0.0, 0.0
  sumRadiiDiff = 0.0
  #figure out which % of the center (core) to use
  outside = (100. - core)/200.
  dataCounted = 0.
  targetMappedTo, targetMappedToIndices = [], []
  withinOne, withinRadius = 0, 0
  for index, sourceDatum in enumerate(sourceData):
    #make sure in 'core'
    if float(index)/float(len(sourceData)) >= outside and \
       float(index)/float(len(sourceData)) < (1. - outside):
      dataCounted += 1.
      #find match
      closest, distance = targetData[0], geometry.distL2(
          targetData[0], sourceDatum)
      tarRad = 1.
      for tarIndex, targetDatum in enumerate(targetData):
        thisDist = geometry.distL2(targetDatum, sourceDatum)
        if thisDist < distance:
          closest = targetDatum
          distance = thisDist
          tarRad = targetRadii[tarIndex]
      if distance < 1.:
        withinOne += 1
      if distance < tarRad:
        withinRadius += 1
      if sourceRadii:
        sumRadiiDiff += abs(tarRad - sourceRadii[index])**2.
      if closest not in targetMappedTo:
        targetMappedTo.append(closest)
        targetIndex = targetData.index(closest)
        targetMappedToIndices.append(targetIndex)
      sumDistanceSquared += distance**2.
      sumWeighted += distance**2.*(1./(tarRad+.0000000001))
  targetMappedToIndices.sort()
  if dataCounted > 0.:
    prmsd = (sumDistanceSquared/float(dataCounted))**0.5
    radiicomp = (sumRadiiDiff/float(dataCounted))**0.5
    coverage = float(len(targetMappedTo)) / float(len(targetData))
    span = float(
        targetMappedToIndices[-1] - targetMappedToIndices[0]+1.) / float(
        len(targetData))
    #percentage of length of target from first covered to last covered
    wrmsd = (sumWeighted/dataCounted)**0.5
    less1 = float(withinOne)/float(dataCounted)
    lessrad = float(withinRadius)/float(dataCounted)
  else:
    prmsd = "err"
    radiicomp = "err"
    coverage = "err"
    span = "err"
    wrmsd = "err"
    less1 = "err"
    lessrad = "err"
  return prmsd, coverage, span, wrmsd, less1, lessrad, radiicomp

#this is main
if -1 != string.find(sys.argv[0], "comparePaths"):
  if len(sys.argv) >= 3:
    source = sys.argv[1]
    target = sys.argv[2]
    core = 100.
    if len(sys.argv) >= 4:
      core = float(sys.argv[3])
    metrics = comparePathsManyMetrics(source, target, core=core)
    prmsd = metrics[0]
    wrmsd = metrics[3]
    print wrmsd
  else:
    print "Usage: comparePaths.py source target [%core]"
    print "where source and target are two python cgo files for paths"
    print "%core determines how much of the source to use, 100 is default"
