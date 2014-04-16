#ryangc@mail.med.upenn.edu ryan.g.coleman@gmail.com
#Ryan G Coleman, Kim Sharp http://crystal.med.upenn.edu

#file contains lots of path primitive methods

import pdb
import geometry
import statistics
import tstdata

def pathLength(path):
  '''computes the total length from point to point'''
  length = 0
  lastPathPt = path[0][1:4]  # init for loop
  for nextPathPtRad in path[1:]:
    nextPathPt = nextPathPtRad[1:4]
    length += geometry.distL2(nextPathPt, lastPathPt)
    lastPathPt = nextPathPt
  return length

def pathCrowFliesLength(path):
  '''computes the distance from the first to the last point'''
  firstPt = path[0][1:4]
  lastPt = path[-1][1:4]
  return geometry.distL2(firstPt, lastPt)

def computeWindingMetric(path):
  '''returns length/crowflies length, higher is windingier'''
  return pathLength(path) / (pathCrowFliesLength(path) + .00000000001)

def findAllMinima(path, eps=0.01):
  '''finds the minima, including one minima from a length that isn't
  strictly a minima (several low values together). returns location of minima'''
  firstInMinima = -1  # means not in a minima
  listIndices = []
  for index, pathPt in enumerate(path):
    thisRad = pathPt[0]
    if firstInMinima == -1 and index != 0:  # first doesn't count
      if thisRad < path[index-1][0]:
        firstInMinima = index
    elif index != 0:
      lastRad = path[firstInMinima][0]
      if lastRad == thisRad or abs(lastRad - thisRad) < eps:
        if thisRad < lastRad:  # without eps
          firstInMinima = index
      elif lastRad > thisRad:
        firstInMinima = index
      else:  # lastrad < thisRad
        listIndices.append(firstInMinima)
        firstInMinima = -1
  if firstInMinima != -1:
    listIndices.append(firstInMinima)
  return listIndices

def insideMaxRadius(path, indexOne, indexTwo):
  '''finds the max radius inbetween the two indices, inclusive'''
  maxInsideRadius = 0.
  for pathPt in path[indexOne:indexTwo+1]:
    if not maxInsideRadius or (pathPt[0] and pathPt[0] > maxInsideRadius):
      maxInsideRadius = pathPt[0]
  return maxInsideRadius

def insideTwoMinimaRadiusMax(path):
  '''finds 2 global minima, finds maxima between them'''
  listMinima = findAllMinima(path)
  if len(listMinima) == 0:
    minRad = pathMinRadius(path)
    return minRad, minRad, minRad
  elif len(listMinima) == 1:
    minRad = path[listMinima[0]][0]
    return minRad, minRad, minRad
  else:  # actually enough minima
    lowest = listMinima[0]
    secondLowest = listMinima[1]
    if path[lowest][0] > path[secondLowest][0]:
      lowest, secondLowest = secondLowest, lowest  # switch
    for otherIndex in listMinima[2:]:
      radius = path[otherIndex][0]
      if radius < path[lowest][0]:
        secondLowest, lowest = lowest, otherIndex
      elif radius < path[secondLowest][0]:
        secondLowest = otherIndex
    minimas = path[lowest][0], path[secondLowest][0]  # save for later
    if lowest > secondLowest:
      lowest, secondLowest = secondLowest, lowest
    insideMax = insideMaxRadius(path, lowest, secondLowest)
    return insideMax, minimas[0], minimas[1]

def pathMinimaCount(path):
  '''finds the number of minima along the path'''
  minimaIndices = findAllMinima(path)
  return len(minimaIndices)

def averageTheta(path):
  '''calculates the mean theta of ab to bc for each triple of points'''
  if len(path) <= 2:
    return [], 0  # no reason to count, not enough nodes
  thetas = []
  firstPt = path[0][1:4]  # first node
  secondPt = path[1][1:4]  # second node
  for node in path[2:]:
    firstVec = geometry.getVector(secondPt, firstPt)
    secondVec = geometry.getVector(node[1:4], secondPt)
    theta = geometry.getAngle(firstVec, secondVec)
    thetas.append(theta)
    firstPt = secondPt
    secondPt = node[1:4]
  averageTheta = statistics.computeMean(thetas)
  return thetas, averageTheta

def pathMaxDistance(nodePath, distanceName):
  '''finds the max of the distance like 'traveldepth'''
  maxDist = nodePath[0].distances[distanceName]
  for node in nodePath[1:]:  # no need to check first again
    if node.distances[distanceName] > maxDist:
      maxDist = node.distances[distanceName]
  return maxDist

def pathMinRadius(path):
  '''finds the overall minimum of the radius (stored in 0th element)'''
  minRadius = False
  for pathPt in path:
    if not minRadius or pathPt[0] < minRadius:
      minRadius = pathPt[0]
  return minRadius

def pathMaxInsideRadius(path):
  '''walks in from both ends to find minima, then finds maximum between them,
  returns max, min1, min2 even though might all be the same'''
  bottom, top = 0, len(path) - 1
  bottomVal, topVal = path[bottom][0], path[top][0]
  while bottom < len(path) and path[bottom+1][0] <= bottomVal:
    bottom += 1
    bottomVal = path[bottom][0]
  while top > bottom and path[top-1][0] <= topVal:
    top -= 1
    topVal = path[top][0]
  maxInsideRadius = 0.
  for pathPt in path[bottom:min(top + 1, len(path) - 1)]:
    if not maxInsideRadius or (pathPt[0] and pathPt[0] > maxInsideRadius):
      maxInsideRadius = pathPt[0]
  if not maxInsideRadius:
    maxInsideRadius = path[bottom][0]
  #print bottom, top, maxInsideRadius, path[bottom], path[bottom + 1]
  return maxInsideRadius, bottomVal, topVal  # max and the extremal mins

def outputNodesText(listNodes, outFileName):
  '''outptuts coordinates in text format'''
  fileTemp = open(outFileName, 'w')
  for node in iter(listNodes):
    for coord in node.xyz:
      fileTemp.write(str(coord) + " ")
    fileTemp.write("\n")
  fileTemp.close()

def getResiduesBetweenPoints(pair, pdbD):
  '''a pair of nodes, finds center and radius, returns residues within'''
  aXYZ = pair[0].getXYZ()
  bXYZ = pair[1].getXYZ()
  radius = geometry.distL2(aXYZ, bXYZ)/2.
  pointRad = [radius]
  for index in range(3):
    pointRad.append((aXYZ[index]+bXYZ[index])/2.)
  resList = getNearbyResidues([pointRad], pdbD)
  return resList

def getNearbyResidues(pointPath, pdbD, nearbyDistance=0.):
  '''returns a list of residues in the pdbD near the pointpath within rad+nd'''
  residuesNearPath = []
  for pathPt in pointPath:
    for index, coord in enumerate(pdbD.coords):
      distanceBetween = geometry.distL2(pathPt[1:4], coord)
      if distanceBetween < pathPt[0] + nearbyDistance:
        residueNumber = pdbD.resNums[index]
        chain = pdbD.chains[index]
        resChain = str(residueNumber) + str(chain)
        if resChain not in residuesNearPath:  # guarantee uniqueness
          residuesNearPath.append(resChain)
  residuesNearPath.sort()
  return residuesNearPath

def outputNearbyResidues(pointPath, outName, pdbRawData, nearbyDistance):
  '''outputs 2 different files. 1. residues near the minimum radius
  2. residues near the entire path.
  near is defined as radius + nearbyDistance'''
  pdbD = pdb.pdbData()
  pdbD.processLines(pdbRawData)
  minRadiusPoint = pointPath[0]
  for pathPt in pointPath:
    if pathPt[0] < minRadiusPoint[0]:
      minRadiusPoint = pathPt
  residuesNearMin = getNearbyResidues([minRadiusPoint], pdbD, nearbyDistance)
  pdbNearMin = pdbD.getListResiduesChains(residuesNearMin)
  pdbNearMin.write(outName + ".residues.pathmin.pdb")
  residuesNearPath = getNearbyResidues(pointPath, pdbD, nearbyDistance)
  pdbNearPath = pdbD.getListResiduesChains(residuesNearPath)
  pdbNearPath.write(outName + ".residues.path.pdb")
  #that's it... files written, return residue lists in case
  return residuesNearMin, residuesNearPath

def outputRadiiTxt(origPath, txtfile):
  distanceRadiusPairs = [(0, origPath[0][0])]  # first pair is 0, first radius
  lastPt, lastDist = origPath[0][1:], 0
  for pt in origPath[1:]:  # all but first
    newDistance = geometry.distL2(pt[1:], lastPt)
    newRadius = pt[0]
    lastPt = pt[1:]
    lastDist += newDistance
    distanceRadiusPairs.append((lastDist, newRadius))
  outputFile = open(txtfile, 'w')
  for distance, radius in distanceRadiusPairs:
    outputFile.write(str(distance) + ", " + str(radius) + " \n")
  outputFile.close()

def checkPath(path, loopPointsList, pointXYZ, xyzStart=1):
  '''checks to see if the path intersects any topological loop on the surf'''
  pathThrough = False, False  # return a tuple... false is for failure
  for loopPts in loopPointsList:
    triangles = tstdata.trianglinizeLoop(loopPts)
    numberIntersects = 0
    #do intersection checks... carefully
    for triangle in triangles:
      lastPathPt = path[0][xyzStart:xyzStart + 3]  # init for loop
      for nextPathPtRad in path[1:]:
        nextPathPt = nextPathPtRad[xyzStart:xyzStart+3]
        triPts0 = pointXYZ[triangle[0]-1][1:]
        triPts1 = pointXYZ[triangle[1]-1][1:]
        triPts2 = pointXYZ[triangle[2]-1][1:]
        posPt, maxIt = False, 5000
        while False == posPt:
          posPt = geometry.linePlaneIntersectionNumeric(
              triPts0, triPts1, triPts2, lastPathPt, nextPathPt)
          if False == posPt:
            triPts0, triPts1, triPts2 = geometry.perturbTriangle(
                triPts0, triPts1, triPts2)
            maxIt -= 1
            if maxIt < 0:
              print "had to perturb points 5000 times", triPts0, triPts1,
              print triPts2, lastPathPt, nextPathPt, "giving up"
              sys.exit(1)
        if posPt is not False:
          if geometry.distL2(lastPathPt, nextPathPt) >= \
              geometry.distL2(lastPathPt, posPt) and \
              geometry.distL2(lastPathPt, nextPathPt) >= \
              geometry.distL2(nextPathPt, posPt):
            if geometry.intPointInsideTri(triPts0, triPts1, triPts2, posPt):
              numberIntersects += 1
        lastPathPt = nextPathPt  # for next loop
    #print numberIntersects  # for debugging...
    if 1 == numberIntersects % 2:  # if intersects odd number of times
      pathThrough = triangles, loopPts
      break  # no need to do more checks, one is good enough
  return pathThrough  # in case caller wants to do something with it.

def checkPathPoints(path, loopPointsList, pointXYZ):
  '''checks to see if the path intersects any topological loop on the surf'''
  pathThrough = False, False  # return a tuple... false is for failure
  for loopPts in loopPointsList:
    triangles = tstdata.trianglinizeLoop(loopPts)
    numberIntersects = 0
    #do intersection checks... carefully
    for triangle in triangles:
      lastPathPt = path[0]  # init for loop
      for nextPathPt in path[1:]:
        triPts0 = pointXYZ[triangle[0]-1][1:]
        triPts1 = pointXYZ[triangle[1]-1][1:]
        triPts2 = pointXYZ[triangle[2]-1][1:]
        posPt, maxIt = False, 5000
        while False == posPt:
          posPt = geometry.linePlaneIntersectionNumeric(
              triPts0, triPts1, triPts2, lastPathPt, nextPathPt)
          if False == posPt:
            triPts0, triPts1, triPts2 = geometry.perturbTriangle(
                triPts0, triPts1, triPts2)
            maxIt -= 1
            if maxIt < 0:
              print "had to perturb points 5000 times",
              print triPts0, triPts1, triPts2,
              print lastPathPt, nextPathPt, "giving up"
              sys.exit(1)
        if posPt is not False:
          if geometry.distL2(lastPathPt, nextPathPt) >= \
              geometry.distL2(lastPathPt, posPt) and \
              geometry.distL2(lastPathPt, nextPathPt) >= \
              geometry.distL2(nextPathPt, posPt):
            if geometry.intPointInsideTri(triPts0, triPts1, triPts2, posPt):
              numberIntersects += 1
        lastPathPt = nextPathPt  # for next loop
    if 1 == numberIntersects % 2:  # if intersects odd number of times
      pathThrough = triangles, loopPts
      break  # no need to do more checks, one is good enough
  return pathThrough  # in case caller wants to do something with it.

def subtractPaths(pathA, pathB):
  '''removes all points in pathB from pathA, returns new path'''
  newPath = []
  for point in pathA:
    if point not in pathB:
      newPath.append(point)
  return newPath
