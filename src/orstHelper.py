#Ryan G. Coleman, Kim A. Sharp crystal.med.upenn.edu, ryan.g.coleman@gmail.com
#bunch of helper methods for orthogonal range searching in the context of depth

import rangesearch
import geometry
import operator
import sys

def getIntersectingPts(
    startPt, endPt, longEdge, shortAxis1, shortAxis2, orst, maxIntersect):
  '''helper function that finds the intersecting points on a line,
  knows how to perturb line if intersecting points are odd
  newly added range search tree must be passed in'''
  divideByZeroProblems = False
  try:  # to catch divide by zero issues
    intersectingPts = []
    #intersectingTris = []  # don't need except for debugging
    result = orst.rangeQuery(
        startPt[shortAxis1]-longEdge, startPt[shortAxis1]+longEdge,
        startPt[shortAxis2]-longEdge, startPt[shortAxis2]+longEdge)
    #result[x][3] contains tri tuple indices, if all three there, do...
    checkResults = []
    for resultEntry in result:   # put triangle index into new struct, 3 copies
      for triIndex in resultEntry[3]:   # will be there if whole tri is inside
        checkResults.append(triIndex)
    checkResults.sort(key=operator.itemgetter(10))
    lastPt, ptCount = -1, 0
    for point in checkResults:
      if False == maxIntersect or len(intersectingPts) < maxIntersect:
        if point == lastPt:
          ptCount += 1
        else:
          lastPt, ptCount = point, 1
        if ptCount == 3:
          triTuple = point
          #now check to make sure at least one point of the triangle is on each
          #side of the x/y shortaxis1/2 planes
          short1BelowCt, short1AboveCt, short2BelowCt, short2AboveCt = \
              0, 0, 0, 0
          for pointInTri in triTuple[:3]:
            if pointInTri[shortAxis1] < startPt[shortAxis1]:
              short1BelowCt += 1
            elif pointInTri[shortAxis1] > startPt[shortAxis1]:
              short1AboveCt += 1
            else:  # weird equality case, increment both
              short1BelowCt += 1
              short1AboveCt += 1
            if pointInTri[shortAxis2] < startPt[shortAxis2]:
              short2BelowCt += 1
            elif pointInTri[shortAxis2] > startPt[shortAxis2]:
              short2AboveCt += 1
            else:  # weird equality case, increment both
              short2BelowCt += 1
              short2AboveCt += 1
          #now do check, only do rest if necessary
          if short1BelowCt > 0 and short1AboveCt > 0 and \
             short2BelowCt > 0 and short2AboveCt > 0:
            triPts0 = triTuple[0]
            triPts1 = triTuple[1]
            triPts2 = triTuple[2]
            posPt, maxIt = False, 5000
            while False == posPt:
              posPt = geometry.linePlaneIntersectionNumeric(
                  triPts0, triPts1, triPts2, startPt, endPt)
              if False == posPt:
                triPts0, triPts1, triPts2 = geometry.perturbTriangle(
                    triPts0, triPts1, triPts2)
                maxIt -= 1
                if maxIt < 0:
                  print "perturbed 5000 times", triPts0, triPts1, triPts2,
                  print startPt, endPt, "giving up"
                  sys.exit(1)
            if posPt is not False:
              if geometry.intPointInsideTriTuple(triTuple, posPt):
                intersectingPts.append(posPt)
                #intersectingTris.append([triPts0, triPts1, triPts2])
    #print "inter", len(intersectingPts), len(intersectingTris)
    #tstdebug.debugTris([startPt,endPt],
    #intersectingTris,intersectingPts,"debug."+
    #str(x)+"."+str(y)+"."+str(z)+".line.py")
  except ZeroDivisionError:  # caused by intPt = triangle vertex
    divideByZeroProblems = True
  if divideByZeroProblems or len(intersectingPts) % 2 == 1:
    return False
  else:
    return intersectingPts

def buildOrst(
    triTuples, shortAxis1, shortAxis2, limitPoints, pointTriList, pointList):
  '''helper that builds the ortho range search tree from the triangle points'''
  inputData = []
  for pointTriNum in xrange(len(pointTriList)):
    inputData.append(False)  # initialize to empty
  for ptIndex, pointTri in enumerate(pointTriList):
    if pointTri[0] in limitPoints:
      triList = pointTri[2:]
      modifiedTriList = []
      for thisTri in triList:
        modifiedTriList.append(triTuples[thisTri])
      #disabled, use tstCheckTris for now
      #if pointTri[1] == 13:
      #  #print triList, modifiedTriList, "something wrong with tst file"
      pointXYZ = pointList[pointTri[0]-1]
      inputData[ptIndex] = [
          pointXYZ[1], pointXYZ[2], pointXYZ[3], modifiedTriList]
  newInputData = []
  for inputRow in inputData:
    if inputRow is not False:
      newInputData.append(inputRow)
  if len(newInputData) > 0:
    orst = rangesearch.orthoRangeSearchTree(
        newInputData, shortAxis1, shortAxis2)
    return orst
  else:
    return False

def decideInsideLong(
    emptyGrid, triTuples, longAxis, allPoints, pointTriList, pointList,
    triList, valueToSet, valueFromSet=False, maxIntersect=False):
  '''helper function that does proper stuff depending on longAxis value'''
  #don't bother figuring out what's inside the ms if already bad grid
  lenX = len(emptyGrid)
  lenY = len(emptyGrid[0])
  lenZ = len(emptyGrid[0][0])
  lens = (lenX, lenY, lenZ)
  shortAxis1, shortAxis2 = 0, 0
  if longAxis == 0:
    shortAxis1, shortAxis2 = 1, 2
  elif longAxis == 1:
    shortAxis1, shortAxis2 = 0, 2
  elif longAxis == 2:
    shortAxis1, shortAxis2 = 0, 1
  longEdge = 1.0000001 * geometry.getLongestEdge(triList, pointList, longAxis)
  #build orthogonal range search tree structure now, save in orst
  orst = buildOrst(
      triTuples, shortAxis1, shortAxis2, allPoints, pointTriList, pointList)
  if False != orst:
    for short1 in range(0, lens[shortAxis1]):
      #print " "
      for short2 in range(0, lens[shortAxis2]):
        startPt, endPt = [], []  # set these up in next switch
        if longAxis == 0:
          x, y, z = 0, short1, short2
          yCoord = emptyGrid[x][y][z][2]
          zCoord = emptyGrid[x][y][z][3]
          xCoordStart = emptyGrid[0][y][z][1]
          xCoordEnd = emptyGrid[lenX - 1][y][z][1]
          startPt = [xCoordStart, yCoord, zCoord]
          endPt = [xCoordEnd, yCoord, zCoord]
        elif longAxis == 1:
          x, y, z = short1, 0, short2
          xCoord = emptyGrid[x][y][z][1]
          zCoord = emptyGrid[x][y][z][3]
          yCoordStart = emptyGrid[x][0][z][2]
          yCoordEnd = emptyGrid[x][lenY - 1][z][2]
          startPt = [xCoord, yCoordStart, zCoord]
          endPt = [xCoord, yCoordEnd, zCoord]
        elif longAxis == 2:
          x, y, z = short1, short2, 0
          xCoord = emptyGrid[x][y][z][1]
          yCoord = emptyGrid[x][y][z][2]
          zCoordStart = emptyGrid[x][y][0][3]
          zCoordEnd = emptyGrid[x][y][lenZ - 1][3]
          startPt = [xCoord, yCoord, zCoordStart]
          endPt = [xCoord, yCoord, zCoordEnd]
        intersectingPts, maxIt, usedStart, usedEnd = False, 5000, startPt, endPt
        while False == intersectingPts:
          intersectingPts = getIntersectingPts(
              usedStart, usedEnd, longEdge, shortAxis1, shortAxis2, orst,
              maxIntersect)
          usedStart, usedEnd = geometry.perturbLine(
              longAxis, shortAxis1, shortAxis2, usedStart, usedEnd, maxIt)
          maxIt -= 1
          if maxIt < 0:
            print "had to perturb line 5000 times..."
            print usedStart, usedEnd
            sys.exit(1)
        if len(intersectingPts) > 0:
          #print len(intersectingPts),
          #check even-ness... perhaps perturb if odd
          if len(intersectingPts) % 2 == 1:
            pass
            #print "odd number of intersecting points!!"
            #perturb starting line, try again
          #need to sort based on ascending long-axis int pt
          #longInts mean long dimension intercepts
          longInts = [xb[longAxis] for xb in intersectingPts]  # make quick list
          longInts.sort()
          #then decide inside/outside ms points, put -2 in grid if inside
          lPlace, lInside = 0, False  # records which intercepts
                       # have been seen, and whether currently outside or inside
          for longCount in range(0, lens[longAxis]):
            x, y, z = -1, -1, -1
            if longAxis == 0:
              x, y, z = longCount, short1, short2
            elif longAxis == 1:
              x, y, z = short1, longCount, short2
            elif longAxis == 2:
              x, y, z = short1, short2, longCount
            gridPiece = emptyGrid[x][y][z]
            while lPlace < len(longInts) and \
                gridPiece[longAxis+1] > longInts[lPlace]:  # switch
              lPlace += 1
              lInside = not lInside
            if lInside:  # replace the encoding int
              if False == valueFromSet or valueFromSet == gridPiece[0]:
                newGridPiece = valueToSet, gridPiece[1], gridPiece[2], \
                    gridPiece[3]
                emptyGrid[x][y][z] = newGridPiece
  #emptyGrid has been modified, nothing returned...

def decideInside(
    emptyGrid, triTuples, allPoints, pointTriList, pointList, triList,
    valueToSet, valueFromSet=False, maxIntersect=False):
  '''moved this into a function for better readability, etc.
  just modifies emptyGrid, doesn't return anything
  calls helper function after determining which axis is longer'''
  #don't bother figuring out what's inside the ms if already bad grid
  lenX = len(emptyGrid)
  lenY = len(emptyGrid[0])
  lenZ = len(emptyGrid[0][0])
  #loop through shortest 2 first is most efficient
  if lenZ >= lenX and lenZ >= lenY:
    decideInsideLong(
        emptyGrid, triTuples, 2, allPoints, pointTriList, pointList,
        triList, valueToSet, valueFromSet, maxIntersect)
  elif lenY >= lenX and lenY >= lenZ:
    decideInsideLong(
        emptyGrid, triTuples, 1, allPoints, pointTriList, pointList,
        triList, valueToSet, valueFromSet, maxIntersect)
  elif lenX >= lenZ and lenX >= lenY:
    decideInsideLong(
        emptyGrid, triTuples, 0, allPoints, pointTriList, pointList,
        triList, valueToSet, valueFromSet, maxIntersect)
  #that's the end, emptyGrid has been modified, nothing returned...

def decideInsidePhi(
    phiData, triTuples, allPoints, pointTriList, pointList, triList,
    valueToSet=1, valueToUnset=0, maxIntersect=False):
  '''helper function that puts values into a phimap instead of a grid.
  cubic so doesn't matter which axis, just use x.
  ToSet is inside surface, ToUnset is outside the surface.
  the values in PhiData are destroyed and replaced with values.'''
  longEdge = 1.0000001 * geometry.getLongestEdge(triList, pointList, 0)
  #build orthogonal range search tree structure now, save in orst
  orst = buildOrst(
      triTuples, 1, 2, allPoints, pointTriList, pointList)
  if False != orst:
    for short1 in xrange(0, phiData.gridDimension):
      #print " "
      for short2 in xrange(0, phiData.gridDimension):
        startPt, endPt = [], []  # set these up next
        x, y, z = 0, short1, short2
        yCoord = phiData.getXYZ(x, y, z)[1]
        zCoord = phiData.getXYZ(x, y, z)[2]
        xCoordStart = phiData.getXYZ(x, y, z)[0]
        xCoordEnd = phiData.getXYZ(phiData.gridDimension - 1, y, z)[0]
        startPt = [xCoordStart, yCoord, zCoord]
        endPt = [xCoordEnd, yCoord, zCoord]
        intersectingPts, maxIt, usedStart, usedEnd = False, 5000, startPt, endPt
        while False == intersectingPts:
          intersectingPts = getIntersectingPts(
              usedStart, usedEnd, longEdge, 1, 2, orst, maxIntersect)
          usedStart, usedEnd = geometry.perturbLine(
              0, 1, 2, usedStart, usedEnd, maxIt)
          maxIt -= 1
          if maxIt < 0:
            print "had to perturb line 5000 times..."
            print usedStart, usedEnd
            sys.exit(1)
        if len(intersectingPts) > 0:
          #need to sort based on ascending long-axis int pt
          #longInts mean long dimension intercepts
          longInts = [xb[0] for xb in intersectingPts]  # make quick array
          longInts.sort()
          #then decide inside/outside ms points, put -2 in grid if inside
          lPlace, lInside = 0, False  # records which intercepts
                       # have been seen, and whether currently outside or inside
                       #outside is False!
          for longCount in xrange(0, phiData.gridDimension):
            x, y, z = longCount, short1, short2
            gridPiece = phiData.getXYZ(x, y, z)
            while lPlace < len(longInts) and \
                gridPiece[0] > longInts[lPlace]:  # switch
              lPlace += 1
              lInside = not lInside
            if lInside:  # replace the encoding int
              phiData.setValue(x, y, z, valueToSet)
            else:
              phiData.setValue(x, y, z, valueToUnset)
        else:  # no intersecting points
          for longCount in xrange(0, phiData.gridDimension):
            x, y, z = longCount, short1, short2
            phiData.setValue(x, y, z, valueToUnset)
  #print phiData.countValues()
  #phiData has been modified, no return value
