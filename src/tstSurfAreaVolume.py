#!/usr/bin/env python
#ryan g. coleman ryangc@mail.med.upenn.edu
#geometry is great

import sys
import string
import math
import tstdata
import cavity
import geometry

def calculateArea(allTris, triPoints, pointXYZ):
  '''calculates whole surface area based on triangles'''
  area = 0.0
  for tri in allTris:
    points = triPoints[tri-1][1:]
    area += geometry.calcTriArea(
        pointXYZ[points[0]-1][1:], pointXYZ[points[1]-1][1:],
        pointXYZ[points[2]-1][1:])
  return area

def calculateVolume(allTris, triPoints, pointXYZ):
  '''calculates whole area based on triangles and a random point (0,0,0).
  individual formula is based on http://mathworld.wolfram.com/Plane.html'''
  volume = 0.0
  areaTotalCheck = 0.0
  origin = (0., 0., 0.)  # could be any point, but this one is fine
  for tri in allTris:
    #calculate volume
    triPointXYZ = []  # 3 points will be here soon
    for point in triPoints[tri-1][1:]:
      triPointXYZ.append(pointXYZ[point-1][1:])
    area = geometry.calcTriAreaList(triPointXYZ)
    normal = geometry.getTriNormalList(triPointXYZ)
    centerTri = geometry.getAverage(triPointXYZ)
    planeD = geometry.calculatePlaneD(normal, centerTri)
    #print triPointXYZ, normal, area
    height = geometry.planeDistToOrigin(normal+[planeD])
    #this height is signed and we can just add it to volume.
    #volume is 1/3 area * height
    thisVol = area * height / 3.
    #print area, normal
    #print height, thisVol,
    volume += thisVol
    areaTotalCheck += area  # free to check other computation.
    #determine whether the point is on the same side as the normal. +/-
    #can do this but not necessary
    #side = geometry.checkPlaneSide(normal+[planeD], origin)
    #if (side and thisVol < 0.) or (not side and thisVol > 0.) or thisVol == 0.:
    #  print side, thisVol
  if volume > 0.:
    return volume
  else:
    return -volume  # since sometimes the origin is inside the shape

def calculateSphericity(area, volume):
  '''from wikipedia http://en.wikipedia.org/wiki/Sphericity
  sphericity = pi^(1/3)(6volume)^(2/3) / area'''
  return ((math.pi**(1./3.))*((6*volume)**(2./3.)))/area

def tstSurfVolAnalysis(tstFileName):
  tstD = tstdata.tstData(
      tstFileName, necessaryKeys=tstdata.tstData.necessaryKeysForSV)
  #could be less
  surfPoints, surfTris, cavPoints, cavTris = cavity.findBiggestDisjointSets(
      tstD.dict['POINT_XYZ'], tstD.dict['TRIANGLE_POINT'],
      tstD.dict['POINT_NEIGHBOR'])
  chVol = False
  if 'CONVEX_HULL_TRI_POINT_LIST' in tstD.dict:
    chVol = True
    chTris = [entry[0] for entry in tstD.dict['CONVEX_HULL_TRI_POINT_LIST']]
  #allTris = surfTris.union(cavTris) #combined
  #print len(tstD.dict['TRIANGLE_POINT']), len(surfTris), len(cavTris) #check
  allArea = calculateArea(
      surfTris, tstD.dict['TRIANGLE_POINT'], tstD.dict['POINT_XYZ'])
  allVol = calculateVolume(
      surfTris, tstD.dict['TRIANGLE_POINT'], tstD.dict['POINT_XYZ'])
  cavArea = calculateArea(
      cavTris, tstD.dict['TRIANGLE_POINT'], tstD.dict['POINT_XYZ'])
  cavVol = calculateVolume(
      cavTris, tstD.dict['TRIANGLE_POINT'], tstD.dict['POINT_XYZ'])
  chArea = False
  if chVol:
    chArea = calculateArea(
        chTris, tstD.dict['CONVEX_HULL_TRI_POINT_LIST'], tstD.dict['POINT_XYZ'])
    chVol = calculateVolume(
        chTris, tstD.dict['CONVEX_HULL_TRI_POINT_LIST'], tstD.dict['POINT_XYZ'])
  #allcavVol = calculateVolume(allTris, \
  #                     tstD.dict['TRIANGLE_POINT'], \
  #                     tstD.dict['POINT_XYZ'])
  #print allVol, cavVol, allcavVol  # again another check not necessary anymore
  return allArea, allVol, allArea + cavArea, allVol - cavVol, chArea, chVol

#this is where main is...
if -1 != string.find(sys.argv[0], "tstSurfAreaVolume.py"):
  if 1 < len(sys.argv):
    print "tst\tsurfaceArea\tvolume\tcavSurfArea\tcavVolume\tspher\t" + \
        "cavSpher\tchArea\tchVol"
    for tstfile in sys.argv[1:]:
      surfArea, volume, cavA, cavV, chA, chV = tstSurfVolAnalysis(tstfile)
      spher1 = geometry.calculateSphericity(surfArea, volume)
      spher2 = geometry.calculateSphericity(cavA, cavV)
      print tstfile, "\t", surfArea, "\t", volume, "\t", cavA, "\t",
      print cavV, "\t", spher1, "\t", spher2, "\t", chA, "\t", chV
  else:
    print "Usage: tstSurfAreaVolume.py tstfile [list of more tstfiles]"
