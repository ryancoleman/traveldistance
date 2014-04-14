#!/usr/bin/env python

import sys,string
import tstdata
import geometry
import statistics
import math #for pi

def tstEdgeCurvature(trianglePoint, pointXyz, pointTriangle, pointNeighbor):
  '''for each edge, calculate the angle between the triangles around it.
  calculate point curvature based on average of these for each point'''
  triXyz = {}
  for triPtList in trianglePoint:
    tri = triPtList[0]
    xyz = []
    for pt in triPtList[1:]:
      xyz.append(pointXyz[pt-1][1:])
    triXyz[tri] = xyz
  edgeAngle = {} #store edge angles as they are found so don't duplicate work
  pointMeanAngle = [] #once all edges found, calculate mean, store in tst format
  pointWeightedMeanAngle = [] #weight by edge length
  for pointNeighborList in pointNeighbor:
    mainPt = pointNeighborList[0]
    angles = []
    weightedAngles = []
    for otherPt in pointNeighborList[2:]: #pN[1] is count
      ptList = [mainPt, otherPt]
      ptList.sort()
      ptTuple = tuple(ptList) #canonicalized format
      edgeLength = geometry.distL2(pointXyz[mainPt-1][1:], \
                                   pointXyz[otherPt-1][1:])
      if ptTuple in edgeAngle: #already done
        angles.append(edgeAngle[ptTuple])
        weightedAngles.append(edgeAngle[ptTuple]*edgeLength)
      else: #have to compute it
        mainTris = set(pointTriangle[mainPt-1][2:])
        otherTris = set(pointTriangle[otherPt-1][2:])
        tris = list(mainTris.intersection(otherTris)) #will almost always be 2
        #for now assume only 2
        normalA = geometry.getTriNormalList(triXyz[tris[0]])
        normalB = geometry.getTriNormalList(triXyz[tris[1]])
        unsignedAngle = geometry.getAngle(normalA, normalB) #unsigned
        centerTriA = geometry.getAverage(triXyz[tris[0]])
        planeA = geometry.calculatePlaneD(normalA, centerTriA)
        ptsB = set(trianglePoint[tris[1]-1][1:])
        edgePts = set(ptList)
        otherB = pointXyz[list(ptsB.difference(edgePts))[0]-1][1:]
        side = geometry.checkPlaneSide(normalA+[planeA], otherB)
        if side:
          angle = - unsignedAngle * 180 / math.pi #concave negative
        else:
          angle = unsignedAngle * 180 / math.pi #convex positive
        edgeAngle[ptTuple] = angle
        angles.append(angle)
        weightedAngles.append(angle*edgeLength)
    pointMeanAngle.append([mainPt, statistics.computeMean(angles)])
    pointWeightedMeanAngle.append([mainPt, \
                                  statistics.computeMean(weightedAngles)])
  return edgeAngle, pointMeanAngle, pointWeightedMeanAngle

#this is main
if -1 != string.find(sys.argv[0], "tstCurvature.py"):
  for tstFileName in sys.argv[1:]:
    tstD = tstdata.tstData(tstFileName, \
                           necessaryKeys=tstdata.tstData.necessaryKeysForCurve)
    eA, pA, pWA = tstEdgeCurvature(tstD.dict['TRIANGLE_POINT'], \
                     tstD.dict['POINT_XYZ'], \
                     tstD.dict['POINT_TRIANGLE'], tstD.dict['POINT_NEIGHBOR'])
    '''
    #append curvature to tst file
    tstFile = open(tstFileName, 'a')
    tstdata.writeEntrySingleFloat(pWA, "POINT_CURVATURE_EDGE LIST", \
                            "END POINT_CURVATURE_EDGE", tstFile)
    tstFile.close()
    '''
    curves,absCurves = [],[]
    for pointWeightCurv in pWA:
      curves.append(pointWeightCurv[1])
      absCurves.append(abs(pointWeightCurv[1]))
    meanCurv = statistics.computeMean(curves)
    meanAbsCurv = statistics.computeMean(absCurves)
    curves,absCurves = [],[]
    for pointWeightCurv in eA.values():
      curves.append(pointWeightCurv)
      absCurves.append(abs(pointWeightCurv))
    meanCurvE = statistics.computeMean(curves)
    meanAbsCurvE = statistics.computeMean(absCurves)
    print tstFileName, meanCurv, meanAbsCurv, meanCurvE, meanAbsCurvE
