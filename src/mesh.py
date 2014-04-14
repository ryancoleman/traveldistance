#ryan g. coleman ryangc@mail.med.upenn.edu crystal.med.upenn.edu
#copyright Ryan Coleman, Kim Sharp 2007
#contains explicit mesh data structures and code

import math
import string
import unionfind2
import priodict
import geometry  # lots of geometric primitives
import tstdebug
import tm3  # tree/graph structure...
import statistics  # average and such
import pdb
import pca
useNumeric = True  # use numeric, if available
try:
  try:  # use numeric (older faster code)
    import Numeric
  except ImportError:  # fallback to numpy (newer slower code)
    import numpy.oldnumeric as Numeric
except ImportError:  # otherwise fallback
  useNumeric = False  # have to do worse linear dimesion estimate

def returnBoxStructure(curBox):
  '''returns a dict of indexes connected in orthogonal directions only,
  upper corner given as input, used only for mesh helping, so in this module.
  code is hard coded nastiness for speed, basically (bad excuse, code with loops
  would likely not be slower, but might be uglier. easier to debug perhaps).'''
  boxStruct = {}
  thisBox = (curBox[0], curBox[1], curBox[2])
  boxStruct[thisBox] = [
      (thisBox[0] + 1, thisBox[1], thisBox[2]),
      (thisBox[0], thisBox[1] + 1, thisBox[2]),
      (thisBox[0], thisBox[1], thisBox[2] + 1)]
  thisBox = (curBox[0] + 1, curBox[1], curBox[2])
  boxStruct[thisBox] = [
      (thisBox[0] - 1, thisBox[1], thisBox[2]),
      (thisBox[0], thisBox[1] + 1, thisBox[2]),
      (thisBox[0], thisBox[1], thisBox[2] + 1)]
  thisBox = (curBox[0], curBox[1] + 1, curBox[2])
  boxStruct[thisBox] = [
      (thisBox[0] + 1, thisBox[1], thisBox[2]),
      (thisBox[0], thisBox[1] - 1, thisBox[2]),
      (thisBox[0], thisBox[1], thisBox[2] + 1)]
  thisBox = (curBox[0], curBox[1], curBox[2] + 1)
  boxStruct[thisBox] = [
      (thisBox[0] + 1, thisBox[1], thisBox[2]),
      (thisBox[0], thisBox[1] + 1, thisBox[2]),
      (thisBox[0], thisBox[1], thisBox[2] - 1)]
  thisBox = (curBox[0] + 1, curBox[1] + 1, curBox[2])
  boxStruct[thisBox] = [
      (thisBox[0] - 1, thisBox[1], thisBox[2]),
      (thisBox[0], thisBox[1] - 1, thisBox[2]),
      (thisBox[0], thisBox[1], thisBox[2] + 1)]
  thisBox = (curBox[0] + 1, curBox[1], curBox[2] + 1)
  boxStruct[thisBox] = [
      (thisBox[0] - 1, thisBox[1], thisBox[2]),
      (thisBox[0], thisBox[1] + 1, thisBox[2]),
      (thisBox[0], thisBox[1], thisBox[2] - 1)]
  thisBox = (curBox[0], curBox[1] + 1, curBox[2] + 1)
  boxStruct[thisBox] = [
      (thisBox[0] + 1, thisBox[1], thisBox[2]),
      (thisBox[0], thisBox[1] - 1, thisBox[2]),
      (thisBox[0], thisBox[1], thisBox[2] - 1)]
  thisBox = (curBox[0] + 1, curBox[1] + 1, curBox[2] + 1)
  boxStruct[thisBox] = [
      (thisBox[0] - 1, thisBox[1], thisBox[2]),
      (thisBox[0], thisBox[1] - 1, thisBox[2]),
      (thisBox[0], thisBox[1], thisBox[2] - 1)]
  return boxStruct

def createValidFunctionSurfaceZero(listSurfRef):
  '''creates a validator function that returns true only if the lookup on the
  listSurfRef of the surfaceReferencePoint of the surfaceNode is 0. used to
  expand patch only as long as it is hydrophobic'''
  lookupDict = {}
  for items in listSurfRef:
    if items[1] == 0.:
      lookupDict[items[0]] = True
    else:
      lookupDict[items[0]] = False
  return lambda surfNode: lookupDict[surfNode]

class node(object):
  '''only made when you want to give a list of nodes to an outside class/func'''
  classification = False
  #outside(CH)=0, inside(MS)=1, between(CH+MS)=2, surf=3, cavity surface = 5

  def __init__(self, referenceIn, xyzIn, pathXyzIn, distIn):
    self.reference = referenceIn
    self.xyz = xyzIn
    self.pathXyz = pathXyzIn
    self.distances = distIn

class mesh(object):
  '''mesh for travel distance, unified data structure makes computing any
  travel distance easy'''

  def __init__(
      self, grid, pointXYZ, pointNeighbor, gridSize,
      outside, inside, between, allPoints=False, cavPoints=False):
    '''creates a mesh from a grid data structure and a tst data structure.
    if cavPoints then create 1 edge per cavpointlist which is the shortest edge
    possible that connects the cavity to the surface.
    pointsHydro is a point->-1, 0,+1 charge mapping, for calculating % apolar
    surface area of pockets'''
    self.gridNodes = []   # lists of nodes, all other dicts keyed off these
    self.surfaceNodes = []
    self.nodeDist = {}  # moving data storage out of node class into dicts
    self.nodeNeighbors = {}
    self.nodeClass = {}
    self.nodeXyz = {}
    self.nodeOffGridXyz = {}
    self.nodeNode = {}
    #want somewhere to store tracebacks information
    self.tracebacks = {}  # dictionary keyed on the name just like distances
    self.gridNodeInitialize(grid, outside, inside, between)
    self.surfNodeInitialize(grid, pointXYZ, gridSize, cavPoints)
    gridErrors, listBadPairs, badNodesSet = self.surfNodeConnect(
        pointXYZ, pointNeighbor)
    if gridErrors:
      print "warning, grid errors, try changing threshold:",
      print gridErrors, len(listBadPairs), len(badNodesSet)
      exit(1)
    self.gridNodeConnect(grid)
    if cavPoints:  # want to add edges into cavities
      #print "cavpoints", cavPoints
      self.addCavityEdges(cavPoints)
    #now all nodes have been added and neighbors calculated

  def setPtHydro(self, pointsHydro):
    '''
    should be a tst formatted record of the surface point to -1, 0, 1 hydro
    '''
    self.nodeHydro = {}
    if pointsHydro:
      for point, charge in pointsHydro:
        self.nodeHydro[point] = charge

  def setPtCurvature(self, pointsCurv):
    '''should be a tst formatted record of the surface point to curvature'''
    self.nodeCurv = {}
    if pointsCurv:
      for point, curv in pointsCurv:
        self.nodeCurv[point] = curv

  def gridNodeInitialize(self, grid, outside, inside, between):
    '''
    initialize the gridNodes, grid should be outside/inside/between already
    '''
    for indexX, rowX in enumerate(grid):
      for indexY, rowY in enumerate(rowX):
        for indexZ, entryZ in enumerate(rowY):
          indices = (indexX, indexY, indexZ)
          #indices is the reference back into grid
          value = entryZ[0]
          self.gridNodes.append(indices)
          self.nodeXyz[indices] = entryZ[1:]
          self.nodeOffGridXyz[indices] = entryZ[1:]
          self.nodeDist[indices] = {}
          self.nodeNeighbors[indices] = []
          if outside == value:
            self.nodeClass[indices] = 0
          elif inside == value:
            self.nodeClass[indices] = 1
          elif between == value:
            self.nodeClass[indices] = 2
          else:
            print "warning: invalid grid value"

  def surfNodeInitialize(self, grid, pointXYZ, gridSize, cavPoints):
    '''now initialize the surface nodes, connect to 2 gridnodes each, each other
    and connect back from said gridnodes.  if cavpoints then label each one as 5
    more to do later in another step'''
    if cavPoints:
      allCavPoints = set()
      for cavPointList in cavPoints:
        allCavPoints.update(cavPointList)
    realMins = grid[0][0][0][1:4]
    for point in pointXYZ:  # assuming no cavities
      ptIndex, ptXYZ = int(point[0]), point[1:4]
      self.surfaceNodes.append(ptIndex)
      self.nodeXyz[ptIndex] = ptXYZ
      self.nodeOffGridXyz[ptIndex] = ptXYZ
      self.nodeClass[ptIndex] = 3
      self.nodeDist[ptIndex] = {}
      self.nodeNeighbors[ptIndex] = []
      if cavPoints and ptIndex in allCavPoints:
        self.nodeClass[ptIndex] = 5
      maxDir, maxDist = -1, -1.
      for direction in xrange(3):  # find axis it is on
        distance = (ptXYZ[direction] - realMins[direction])
        distance = distance % gridSize  # now between 0 and gridsize
        distance = distance / gridSize  # now between 0 and 1
        if distance > 0.5:
          distance = 1. - distance
        if distance > maxDist:
          maxDir, maxDist = direction, distance
      gridIndices = [], []
      for direction in range(3):
        if direction == maxDir:
          gridIndices[0].append(int(math.floor(
              (ptXYZ[direction] - realMins[direction]) / gridSize)))
          gridIndices[1].append(int(math.ceil(
              (ptXYZ[direction] - realMins[direction]) / gridSize)))
        else:
          for gridIndex in range(2):
            gridIndices[gridIndex].append(int(round(
                (ptXYZ[direction] - realMins[direction]) / gridSize)))
      newGridIndices = []
      for gridIndex in range(2):  # convert to tuple not list
        newGridIndices.append(tuple(gridIndices[gridIndex]))
      for gridIndex in newGridIndices:
        gridPoint = grid[gridIndex[0]][gridIndex[1]][gridIndex[2]][1:4]
        distanceBetween = geometry.distL2(gridPoint, ptXYZ)
        #connect the surface nodes to the grid nodes and back
        self.connectNeighbors(ptIndex, gridIndex, distanceBetween)
        if self.nodeClass[gridIndex] == 0 or self.nodeClass[gridIndex] == 2:
          #outside or between
          dirVector = list(self.nodeXyz[gridIndex])
          offGridXYZ = []
          for coord in range(len(dirVector)):
            dirVector[coord] -= self.nodeXyz[ptIndex][coord]
            dirVector[coord] *= .1
            offGridXYZ.append(self.nodeXyz[ptIndex][coord] + dirVector[coord])
          #print self.gridNodes[gridIndex].xyz, curSurfNode.xyz, offGridXYZ
          self.nodeOffGridXyz[ptIndex] = tuple(offGridXYZ)  # gets set here

  def connectNeighbors(self, node1, node2, distanceBetween):
    '''helper method to connect 2 nodes'''
    self.connectNeighbor(node1, node2, distanceBetween)
    self.connectNeighbor(node2, node1, distanceBetween)  # reverse

  def connectNeighbor(self, node1, node2, distanceBetween):
    '''only does 1 way connection from 1 to 2'''
    if node1 not in self.nodeNeighbors:
      self.nodeNeighbors[node1] = []  # init list
    if (node2, distanceBetween) not in self.nodeNeighbors[node1]:
      self.nodeNeighbors[node1].append((node2, distanceBetween))

  def surfNodeConnect(self, pointXYZ, pointNeighbor):
    '''now connect up all the surface nodes, check for errors'''
    gridErrors, listBadPairs, badNodesSet = False, [], set()
    for point in pointXYZ:  # assuming no cavities
      ptIndex, ptXYZ = int(point[0]), tuple(point[1:4])
      neighbors = []
      for neighbor, neighDist in self.nodeNeighbors[ptIndex]:
        neighbors.append(neighbor)
      if self.nodeClass[neighbors[0]] == self.nodeClass[neighbors[1]]:
        gridErrors = True
        listBadPairs.append(neighbors)
        for neighbor in neighbors:
          badNodesSet.add(neighbor)
      #now connect surface nodes to each other
      for otherPoint in pointNeighbor[ptIndex-1][2:]:
        otherXYZ = pointXYZ[otherPoint-1][1:4]
        self.connectNeighbor(
            otherPoint, ptIndex, geometry.distL2(ptXYZ, otherXYZ))
    return gridErrors, listBadPairs, badNodesSet

  def gridNodeConnect(self, grid):
    '''now connect the gridNodes, go through each box..'''
    lenX, lenY, lenZ = len(grid), len(grid[0]), len(grid[0][0])
    for indexX in xrange(lenX - 1):  # don't do last one
      for indexY in xrange(lenY - 1):
        for indexZ in xrange(lenZ - 1):
          indices = (indexX, indexY, indexZ)
          #indices is the reference back into grid
          boxStruct = returnBoxStructure(indices)
          unionCenters = unionfind2.unionFind()  # initialize
          for oneCenter in boxStruct.keys():
            for otherCenter in boxStruct[oneCenter]:
              if self.nodeClass[oneCenter] == 1:  # oneNode.isInside():
                if self.nodeClass[otherCenter] == 1:  # otherNode.isInside():
                  unionCenters.union(oneCenter, otherCenter)
              elif self.nodeClass[oneCenter] in [0, 2]:  # isOutsideOrBetween()
                if self.nodeClass[otherCenter] in [0, 2]:  # isOutsideOrBetween
                  unionCenters.union(oneCenter, otherCenter)
          unionCenterLists = unionCenters.toLists()
          for oneSet in unionCenterLists:  # connect all within each set
            for indexOne, nodeOne in enumerate(oneSet):
              gridPtOne = grid[nodeOne[0]][nodeOne[1]][nodeOne[2]][1:4]
              for indexTwo, nodeTwo in enumerate(oneSet):
                if indexOne < indexTwo:
                  gridPtTwo = grid[nodeTwo[0]][nodeTwo[1]][nodeTwo[2]][1:4]
                  #have to add grid neighbors
                  distBetween = geometry.distL2(gridPtOne, gridPtTwo)
                  self.connectNeighbors(nodeOne, nodeTwo, distBetween)

  def addCavityEdges(self, cavPtLists):
    '''adds cavity edges which are the shortest travel in distance from one pt
    in each cavity to the surface.'''
    travelInName = "travelinforcavities"
    self.calculateTravelDistance(travelInName, [3], [1, 2, 5])
    #necessary for cavities
    for cavPtList in cavPtLists:
      minEnd, minStart, minDist = False, False, 100.**100.  # huge number
      for aPt in cavPtList:  # now go through whole set, find min
        traceSurfNodes = self.getFinalTracebackSet(set([aPt]), travelInName)
        for traceSurfNode in traceSurfNodes:
          travDist = self.nodeDist[aPt][travelInName]
          if travDist < minDist:
            minDist = travDist
            minStart = traceSurfNode
            minEnd = aPt
      #print "new cavity edge", minEnd, minStart, minDist  # debugging
      #add edges between the 2 surface nodes minEnd and minStart
      self.connectNeighbors(minEnd, minStart, minDist)
    #should delete travelInName stuff now, save memory
    del self.tracebacks[travelInName]
    for aNode in self.gridNodes + self.surfaceNodes:
      if travelInName in self.nodeDist[aNode]:
        del self.nodeDist[aNode][travelInName]

  def calculateNearbyAtoms(self, pdbD, nearbyDistance=0.):
    '''find all the atoms near each surface node'''
    outputList = []
    inputList = []
    for surfNode in self.surfaceNodes:
      inputList.append(self.nodeXyz[surfNode])
    nearbyDict = pdbD.getNearbyAtoms(inputList, nearbyDistance)
    for surfNode in self.surfaceNodes:
      tempList = [surfNode]
      tempList.extend(nearbyDict[tuple(self.nodeXyz[surfNode])])
      outputList.append(tempList)
    return outputList

  def calcTravelDistPair(
      self, validCodes, initialNode, endingNode,
      startDist=0, endMaxDist=False, limitSet=False):
    '''calculates the distance between 2 nodes, honoring validCodes, limitSet,
    endMaxDist and startDist, but no other options. Deletes distance dict after
    finding the distance, also doesn't compute the tracebacks at all'''
    name = str(validCodes) + str(initialNode)  # temporary name
    self.calculateTravelDistance(
        name, validCodes, validCodes, startDist=startDist,
        endMaxDist=endMaxDist, limitSet=limitSet, initialNodes=[initialNode],
        endingNodes=[endingNode], doTracebacks=False)
    returnDist = self.nodeDist[endingNode][name]
    if name in self.tracebacks:
      del self.tracebacks[name]  # don't delete if doesn't exist
    for aNode in self.gridNodes + self.surfaceNodes:
      if name in self.nodeDist[aNode]:
        del self.nodeDist[aNode][name]
    return returnDist

  def calcTravelDistPairLists(
      self, validCodes, initialNodes, endingNodes,
      startDist=0, endMaxDist=False, limitSet=False):
    '''calculates the distance between 2 lists of nodes, honoring validCodes,
    limitSet, endMaxDist and startDist, but no other options.
    Deletes distance dict after
    finding the distance, also doesn't compute the tracebacks at all'''
    name = str(validCodes) + str(initialNodes) + str(endingNodes)
    #temporary name
    returnDist = self.calculateTravelDistance(
        name, validCodes, validCodes, startDist=startDist,
        endMaxDist=endMaxDist, limitSet=limitSet, initialNodes=initialNodes,
        endingNodes=endingNodes, doTracebacks=False)
    if name in self.tracebacks:
      del self.tracebacks[name]  # don't delete if doesn't exist
    for aNode in self.gridNodes + self.surfaceNodes:
      if name in self.nodeDist[aNode]:
        del self.nodeDist[aNode][name]
    return returnDist

  def calcTravelDistOneToMany(
      self, validCodes, initialNodes, endingNodesList,
      startDist=0, endMaxDist=False, limitSet=False):
    '''calculates the distance between from a list of nodes to many other lists
    of nodes, honoring validCodes,
    limitSet, endMaxDist and startDist, but no other options.
    Deletes distance dict after
    finding the distance, also doesn't compute the tracebacks at all'''
    name = str(validCodes) + str(initialNodes)  # temporary name
    self.calculateTravelDistance(
        name, validCodes, validCodes, startDist=startDist,
        endMaxDist=endMaxDist, limitSet=limitSet, initialNodes=initialNodes,
        doTracebacks=False)
    if name in self.tracebacks:
      del self.tracebacks[name]  # don't delete if doesn't exist
    nodesListDist = {}
    for endingNodes in endingNodesList:
      shortestDist = 10000000000.
      for aNode in endingNodes:
        if name in self.nodeDist[aNode] and \
            self.nodeDist[aNode][name] < shortestDist:
          shortestDist = self.nodeDist[aNode][name]
      nodesListDist[endingNodes] = shortestDist
    for aNode in self.gridNodes + self.surfaceNodes:
      if name in self.nodeDist[aNode]:
        del self.nodeDist[aNode][name]
    return nodesListDist

  def calculateCircularPatch(
      self, validCodes, initialNode, startDist=0,
      endMaxArea=False, limitSet=False):
    '''calculates the distance from a node, honoring validCodes, limitSet,
    endMaxDist and startDist, but no other options. returns the circular
    or spherical patch constructed. doesn't do tracebacks.'''
    name = str(validCodes) + str(initialNode)  # temp name
    self.calculateTravelDistance(
        name, validCodes, validCodes, startDist=startDist,
        endMaxArea=endMaxArea, limitSet=limitSet, initialNodes=[initialNode],
        doTracebacks=False)
    patch = self.getOrder(name)
    surfAreaList = self.getSurfaceAreaCumulative(patch)
    return patch, surfAreaList

  def growCircularPatch(
      self, validCodes, initialNode, validFunc,
      startDist=0, endMaxArea=False, limitSet=False):
    '''grows a circular patch that meets the requirements of the valid function,
    returns the patch and the area'''
    name = str(validCodes) + str(initialNode)  # temporary name
    self.calculateTravelDistance(
        name, validCodes, validCodes, startDist=startDist,
        endMaxArea=endMaxArea, limitSet=limitSet, initialNodes=[initialNode],
        doTracebacks=False, newNodeTest=validFunc)
    patch = self.getOrder(name)
    surfArea = self.getSurfaceArea(patch)
    return patch, surfArea

  def growSphericalVolume(
      self, validCodes, initialNode, startDist=0, endMaxVol=False,
      limitSet=False):
    '''grows a spherical volume where each node
    can be inside or outside or between based on the valid
    codes parameter and the limitSet if desired. endmaxvol is when to quit'''
    name = str(validCodes) + str(initialNode)  # temporary name
    self.calculateTravelDistance(
        name, validCodes, validCodes, startDist=startDist, endMaxVol=endMaxVol,
        limitSet=limitSet, initialNodes=[initialNode], doTracebacks=False)
    patch = self.getOrder(name)
    vol = self.getVolume(patch)
    return patch, vol

  def groupByTracebacks(
      self, validCodes, initialNodes, startDist=0, limitSet=False):
    '''method to grow from a set of initial nodes, then group all points into
    groups based on which initial node they trace back to'''
    name = str(initialNodes)   # any name will do
    #does travel out from just the initialNodes, not all surface nodes.
    self.calculateTravelDistance(
        name, validCodes, validCodes, initialNodes=initialNodes,
        startDist=startDist, limitSet=limitSet, doTracebacks=True)
    #now calculation done, for each node, find the init node it traces back to
    groupDict = {}
    for startNode in initialNodes:
      groupDict[startNode] = []
    for aNode in self.gridNodes + self.surfaceNodes:
      if (self.nodeClass[aNode] in validCodes) and \
          (not limitSet or aNode in limitSet):  # takes care of limitset
        #now check to see what it traces back to
        startingNodes = self.getFinalTracebackSet(set([aNode]), name)
        for startNode in startingNodes:  # shouldn't be many duplicates
          groupDict[startNode].append(aNode)
    for startNode in initialNodes:
      groupDict[startNode] = tuple(groupDict[startNode])  # tuples are hashable
    return groupDict

  def calculateTravelDistance(
      self, name, initializationCodes, validCodes, startDist=0,
      endMaxDist=False, limitSet=False, initialNodes=None, endingNodes=None,
      doTracebacks=True, endMaxArea=False, endMaxVol=False,
      newNodeTest=lambda node: True):
    '''does all the work of travel distance calculations for any algorithm
    codes are outside(CH)=0, inside(MS)=1, between(CH+MS)=2, surf=3, cavsurf=5
    returns (and saves) traceback, other data left in nodes,
    restricts computation to nodes in the limitSet if it exists.
    if intialNodes and endingNodes are present, they are used,
    and initCodes ignored. endingnodes just waits to find 1, not all.
    doTracebacks is boolean, either keep track of them or don't (save time).
    newNodeTest is a function that should take a node, and return true if
    it is okay to continue or false if the travel distance calculation should
    stop, and if it is false then it does stop. default is always true.
    '''
    currentNodes = priodict.priorityDictionary()
    #currentNodes holds data on nodes left to process
    if doTracebacks:
      history = {}  # keyed on (node, dist) pairs, goes to prevBox
      self.tracebacks[name] = {}
      #tracebacks assembled during run, keys on ending grid point,
      # points back to list of direct ancestors
    if initialNodes is not None:  # i want to type if initialNodes are not None
      for initialNode in initialNodes:
        if newNodeTest(initialNode):
          self.nodeDist[initialNode][name] = startDist
          currentNodes[initialNode] = startDist
    else:  # initialize distance dictionary for this name for each node
      for aNode in self.gridNodes + self.surfaceNodes:
        if (self.nodeClass[aNode] in initializationCodes) and \
            (not limitSet or aNode in limitSet):  # takes care of limitset
            if newNodeTest(aNode):
              self.nodeDist[aNode][name] = startDist
              actuallyInit = False
              for neighNode, neighDist in self.nodeNeighbors[aNode]:
                if self.nodeClass[neighNode] in validCodes:
                  actuallyInit = True
                  break
              if actuallyInit:  # this way heap is not overloaded...
                currentNodes[aNode] = startDist
    lastTravelDistance = startDist
    lastTravelArea = 0.  # initialize but isn't updated unless needed
    lastTravelVol = 0.  # same as area, init now, not updated unless needed
    #initialization steps over
    while (not endMaxDist or endMaxDist >= lastTravelDistance) and \
        (not endMaxArea or endMaxArea >= lastTravelArea) and \
        (not endMaxVol or endMaxVol >= lastTravelVol) and \
        len(currentNodes) > 0:  # this is the main loop
      currentNode = currentNodes.smallest()
      if not newNodeTest(currentNode):
        break  # stop calculation now, this node is bad (for some reason)
      lastTravelDistance = currentNodes.pop(currentNodes.smallest())
      if endMaxArea:  # keeping track of the area as a stopping criterion
        lastTravelArea += self.nodeSurfArea[currentNode]
      if endMaxVol:  # keeping track of the area as a stopping criterion
        lastTravelVol += self.nodeVol[currentNode]
      if name not in self.nodeDist[currentNode] or \
          self.nodeDist[currentNode][name] >= lastTravelDistance:
        #update the dist, add neighbors to heap
        self.nodeDist[currentNode][name] = lastTravelDistance
        if doTracebacks:
          prevNodeList = history.pop((currentNode, lastTravelDistance), False)
          if prevNodeList and len(prevNodeList) > 0:
            self.tracebacks[name][currentNode] = [
                prevNodeList, currentNode, lastTravelDistance]
        if endingNodes is not None and currentNode in endingNodes:
          break  # stop calculation now, found one of the ending nodes
        for neighborNode, nbDist in self.nodeNeighbors[currentNode]:
          if self.nodeClass[neighborNode] in validCodes and \
              name not in self.nodeDist[neighborNode]:
            if not limitSet or neighborNode in limitSet:
              newDist = lastTravelDistance + nbDist
              if neighborNode not in currentNodes or \
                  newDist <= currentNodes[neighborNode]:
                currentNodes[neighborNode] = newDist  # updates prio dict
                if doTracebacks:
                  prevNodeList = history.get((neighborNode, newDist), [])  # old
                  prevNodeList.append(currentNode)  # now add to it
                  history[(neighborNode, newDist)] = prevNodeList  # now update
    #print lastTravelDistance
    return lastTravelDistance

  def clusterNodesByFunction(
      self,  validCodes, limitSet=False, nodeTest=lambda node: True):
    '''clusters all validnodes together based on the nodeTest function
    codes are outside(CH)=0, inside(MS)=1, between(CH+MS)=2, surf=3
    restricts computation to nodes in the limitSet if it exists.
    if 2 nodes pass the nodeTest function and are neighbors then they are
    joined together. the clusters are returned. this is for growphobic2'''
    currentNodes = set()  # holds data on nodes left to process
    nodeClusters = unionfind2.unionFind()  # initialize structure
    for aNode in self.gridNodes + self.surfaceNodes:
      if (self.nodeClass[aNode] in validCodes) and \
          (not limitSet or aNode in limitSet):  # takes care of limitset
        if nodeTest(aNode):
          currentNodes.add(aNode)  # list to process later
          nodeClusters.find(aNode)  # make each its own set for now
    #initialization steps over
    while len(currentNodes) > 0:  # this is the main loop
      currentNode = currentNodes.pop()  # already been checked by nodeTest
      for neighborNode, nbDist in self.nodeNeighbors[currentNode]:
        if self.nodeClass[neighborNode] in validCodes:  # complicated logic
          if not limitSet or neighborNode in limitSet:
            if nodeTest(neighborNode):
              nodeClusters.union(currentNode, neighborNode)  # join
    #it is really that simple
    return nodeClusters

  def getOrder(self, name):
    '''gets the order in which nodes were found according to the distance'''
    unorderedList = []
    for aNode in self.gridNodes + self.surfaceNodes:
      if name in self.nodeDist[aNode] and self.nodeDist[aNode][name]:
        #exists and notfalse
        unorderedList.append(aNode)
    unorderedList.sort(
        lambda x, y: cmp(self.nodeDist[x][name], self.nodeDist[y][name]))
    return unorderedList  # now actually ordered

  def getMaxName(self, name, iterableNodes):
    '''gets the max distance metric defined by name on iterablenodes'''
    maxi = 0.
    for aNode in iterableNodes:
      if name in self.nodeDist[aNode]:
        maxi = max(maxi, self.nodeDist[aNode][name])
    return maxi

  def getSumName(self, name, iterableNodes):
    '''sums the distance metric defined by name on iterablenodes'''
    sum = 0
    for aNode in iterableNodes:
      if name in self.nodeDist[aNode]:
        sum += self.nodeDist[aNode][name]
    return sum

  def findCountMouths(self, iterableNodes):
    '''find mouths, i.e. between nodes whose neighbor between nodes aren't all
    in the group. separate them into disjoint sets and return the number of sets
    and their sizes'''
    setNodes = set(iterableNodes)
    mouths = set()
    for aNode in iterableNodes:
      if aNode in self.nodeTag:
        if self.nodeTag[aNode] > 0:
          mouths.add(aNode)
    #return mouths #debug finding code
    mouthSets = unionfind2.unionFind()  # initialize
    for aNode in mouths:
      for neighborNode, neighDist in self.nodeNeighbors[aNode]:
        if neighborNode in mouths:
          mouthSets.union(aNode, neighborNode)
    mouthLists = mouthSets.toLists()
    #want to sort by length so biggest is first
    mouthLists.sort(lambda x, y: cmp(len(x), len(y)))
    mouthLists.reverse()
    return mouthLists

  def getMaxDimension(self, iterableNodes, notJustSurface=False):
    '''finds the max difference in each orthogonal dimension. returns max.
    this is linear not quadratic, and always a under-approximation'''
    mins = [100000., 100000., 100000.]
    maxs = [-100000., -100000., -100000.]
    for aNode in iterableNodes:
      if 3 == self.nodeClass[aNode] or 5 == self.nodeClass[aNode] or \
          notJustSurface:  # is surface
        xyz = self.nodeXyz[aNode]
        for direction in xrange(3):
          if xyz[direction] < mins[direction]:
            mins[direction] = xyz[direction]
          if xyz[direction] > maxs[direction]:
            maxs[direction] = xyz[direction]
    bestDist = 0.
    for direction in xrange(3):
      if maxs[direction]-mins[direction] > bestDist:
        bestDist = maxs[direction]-mins[direction]
    return bestDist

  def getFarthestDistance(self, iterableNodes, notJustSurface=False):
    '''looks at all pairs, finds the largest distance, using brute force,
    which appears to be faster than hard to implement sub-quadratic solutions.
    also made it only consider surface nodes since it was too slow.
    falls back to max dimension if points > 100'''
    if useNumeric:  # can do good estimate using PCA
      pointList = []
      for aNode in iterableNodes:
        if 3 == self.nodeClass[aNode] or 5 == self.nodeClass[aNode] or \
            notJustSurface:  # is surface/don't care
          pointList.append(self.nodeXyz[aNode])
      if len(pointList) > 0:
        pcaLong = pca.findLongestDimension(pointList)  # must be list
        return pcaLong
      else:
        return 0.
    else:
      if len(iterableNodes) > 20:  # do linear estimate
        return self.getMaxDimension(iterableNodes, notJustSurface)
      bestDist = 0.
      for aCount, aNode in enumerate(iterableNodes):
        if 3 == self.nodeClass[aNode] or 5 == self.nodeClass[aNode] or \
            notJustSurface:  # is surface
          for bCount, bNode in enumerate(iterableNodes):
            if aCount > bCount:
              if 3 == self.nodeClass[bNode] or 5 == self.nodeClass[bNode] or \
                  notJustSurface:  # is inside
                newDistSq = geometry.distL2Squared(
                    self.nodeXyz[aNode], self.nodeXyz[bNode])
                if newDistSq > bestDist:
                  bestDist = newDistSq
            elif aCount <= bCount:
              break  # no need to compare every pair twice
      return math.sqrt(bestDist)

  def getVolume(self, iterableNodes):
    '''gets the volume of a set of nodes'''
    vol = 0.
    for aNode in iterableNodes:
      vol += self.nodeVol[aNode]
    return vol

  def getApolarSurfaceArea(self, iterableSurfaceNodes):
    '''gets the surface area of just apolar points'''
    surfArea = 0.
    for aNode in iterableSurfaceNodes:
      try:
        if self.nodeHydro[aNode] == 0:
          surfArea += self.nodeSurfArea[aNode]
      except KeyError:
        pass  # this is if you pass in non-surface nodes, which could happen
    return surfArea

  def getCurvatures(self, iterableSurfaceNodes):
    '''gets the mean curvature and the absolute mean curvature (roughness)'''
    curvs, absCurvs = [], []
    for aNode in iterableSurfaceNodes:
      if aNode in self.nodeCurv:
        curv = self.nodeCurv[aNode]
        curvs.append(curv)
        absCurvs.append(abs(curv))
    return statistics.computeMean(curvs), statistics.computeMean(absCurvs)

  def getHydroSurfaceAreas(self, iterableSurfaceNodes):
    '''gets the surface area of each kind of surface point (neg, apolar, pos)'''
    surfAreas = [0., 0., 0.]  # -1, 0, +1 charge
    for aNode in iterableSurfaceNodes:
      if aNode in self.nodeHydro:
        for index, type in enumerate([-1, 0, 1]):
          if self.nodeHydro[aNode] == type:
            surfAreas[index] += self.nodeSurfArea[aNode]
    return surfAreas

  def getSurfaceArea(self, iterableSurfaceNodes):
    '''gets the surface area of a set of surface nodes'''
    surfArea = 0.
    for aNode in iterableSurfaceNodes:
      if aNode in self.nodeSurfArea:
        surfArea += self.nodeSurfArea[aNode]
    return surfArea

  def getSurfaceAreaCumulative(self, iterableSurfaceNodes):
    '''gets the surface area of a set of surface nodes, returns cumulat list'''
    surfArea, surfList = 0., []
    for aNode in iterableSurfaceNodes:
      surfArea += self.nodeSurfArea[aNode]
      surfList.append(surfArea)
    return surfList

  def setSurfaceArea(self, triPointList):
    self.nodeSurfArea = {}
    for aNode in self.surfaceNodes:  # initialize to 0
      self.nodeSurfArea[aNode] = 0.
    for triangle in triPointList:
      points = triangle[1:]
      xyzs = []
      for point in points:
        xyzs.append(self.nodeXyz[point])
      triArea = geometry.calcTriAreaList(xyzs)
      for point in points:
        self.nodeSurfArea[point] += (triArea/3.)  # increment by 1/3rd of tri

  def setVolume(self, gridSize):
    '''sets the volume by simple method where any mesh vertex gets the volume
    of the cube and any surface node gets 0 volume. not perfect but good enough
    for a first pass'''
    self.nodeVol = {}
    for aNode in self.surfaceNodes:  # set all these to 0
      self.nodeVol[aNode] = 0.
    for aNode in self.gridNodes:  # set all these to the volume of a cube
      self.nodeVol[aNode] = gridSize**3.0

  def getSurfaceTravelDistance(self, name):
    '''outputs a list of point numbers and travel distances'''
    pointNums = self.surfaceNodes
    pointNums.sort()
    pointTravelDistance = []
    for pointNum in pointNums:
      thisDist = self.nodeDist[pointNum][name]
      pointTravelDistance.append([pointNum, thisDist])
    return pointTravelDistance

  def getGridTravelDistance(self, grid, name, defaultValue=-1):
    '''copies the travel distance into new grid, doesn't modify passed grid,
    sets grids from False to defaultValue'''
    newGrid = []
    for indexX, rowX in enumerate(grid):
      newX = []
      for indexY, rowY in enumerate(rowX):
        newY = []
        for indexZ, entryZ in enumerate(rowY):
          indices = (indexX, indexY, indexZ)
          #indicies is the reference back into grid
          if name in self.nodeDist[indices]:
            newZ = self.nodeDist[indices][name]
            if not newZ:  # means there is no distance in the dictionary
              newZ = defaultValue
          else:
            newZ = defaultValue
          newY.append((newZ, entryZ[1], entryZ[2], entryZ[3]))
        newX.append(newY)
      newGrid.append(newX)
    return newGrid

  def getAtomNumbers(self, nodes, pointAtomList):
    '''for each node, if it is surface, find the atoms nearby,
    return line numss'''
    numbers = []
    for aNode in nodes:
      if 3 == self.nodeClass[aNode] or 5 == self.nodeClass[aNode]:  # surface
        #print aNode.surfacePointReference, len(pointAtomList),
        #print  pointAtomList[aNode.surfacePointReference-1][0]  # debug
        for number in pointAtomList[aNode-1][1:]:
          if number not in numbers:  # uniqueness
            numbers.append(number)
    return numbers

  def getResChainsNodes(self, numbers, pdbData):
    '''for each node, if it is surface, find the atoms nearby, find the residues
    and return the list of residues by number+chain id'''
    resChainList = pdbData.getResidueChainsFromNums(numbers)
    return pdb.turnListIntoString(resChainList)

  def getResNamesChainsNodes(self, numbers, pdbData):
    '''for each node, if it is surface, find the atoms nearby, find the residues
    and return the list of residues by number+chain id'''
    resChainList = pdbData.getResidueNamesChainsFromNums(numbers)
    return pdb.turnListIntoString(resChainList)

  def getAtomNamesChainsNodes(self, numbers, pdbData):
    '''for each node, if it is surface, find the atoms nearby, find the residues
    and return the list of residues by number+chain id'''
    resChainList = pdbData.getAtomsFromNums(numbers)
    return pdb.turnListIntoString(resChainList)

  def calculatePocketProperties(
      self, nodes, name, pdbData, pointAtomList, ligandNodes, doPCA=True):
    '''wrapper for all these functions since it is called multiple times'''
    surfArea = self.getSurfaceArea(nodes)
    negSA, apolarSA, posSA = self.getHydroSurfaceAreas(nodes)
    meanCurv, meanAbsCurv = self.getCurvatures(nodes)
    volume = self.getVolume(nodes)
    meanTD = self.getSumName(name, nodes) / float(len(nodes))
    maxTD = self.getMaxName(name, nodes)
    mouthLists = self.findCountMouths(nodes)
    numbers = self.getAtomNumbers(nodes, pointAtomList)
    nearbyRes = self.getResChainsNodes(numbers, pdbData)
    nearbyResNames = self.getResNamesChainsNodes(numbers, pdbData)
    nearbyAtoms = self.getAtomNamesChainsNodes(numbers, pdbData)
    if doPCA:  # normal case, turned off for speed if features unused
      pts = []
      for node in nodes:
        pts.append(self.nodeXyz[node])
      dimensions = pca.findDimensions(pts)
    else:  # fake dimensions data
      dimensions = [-1., -1., -1.]  # obviously faked
    biggestMouths = [0., 0.]  # going to use second size as well at some point
    mouthResNames = "+"  # empty is default
    if len(mouthLists) > 0:
      biggestMouths[0] = self.getVolume(mouthLists[0])
      if 0. == biggestMouths[0]:  # means the mouth is all surface
        biggestMouths[0] = self.getSurfaceArea(mouthLists[0])
      farthestPointMouth = self.getFarthestDistance(mouthLists[0])
      if 0. == farthestPointMouth:
        farthestPointMouth = self.getFarthestDistance(
            mouthLists[0], notJustSurface=True)
      numbers = self.getAtomNumbers(mouthLists[0], pointAtomList)
      mouthResNames = self.getResNamesChainsNodes(numbers, pdbData)
    else:
      farthestPointMouth = 0.  # no mouths
      mouthResNames = "+"
    #now want to find the pocket type, i.e. bottleneck, groove, tunnel, cleft
    pocketType = "n/a"  # default
    sizeInt, sizeUnion, jaccard, fraclig = 0., 0., 0., 0.
    if ligandNodes is not None:
      nodeSet = set(nodes)
      unionSet = nodeSet.union(ligandNodes)
      sizeUnion = len(unionSet)
      if sizeUnion > 0:
        intSet = nodeSet.intersection(ligandNodes)
        sizeInt = len(intSet)
        fraclig = float(len(intSet)) / float(len(ligandNodes))
        jaccard = float(sizeInt) / float(sizeUnion)
    #return a LOT of things
    return volume, surfArea, meanTD, mouthLists, nearbyRes, nearbyResNames, \
        nearbyAtoms, biggestMouths[0], farthestPointMouth, pocketType, \
        mouthResNames, maxTD, negSA, apolarSA, posSA, meanCurv, meanAbsCurv, \
        dimensions, sizeInt, sizeUnion, jaccard, Fraclig

  def pocketMapping(
      self, name, validCodes, pointAtomList, pdbData, groupName='group',
      outName="test", ligandNodes=None, doPCA=True):
    '''actual implementation of pocket mapping algorithm.
    1. sorts all between and surface points
    2. goes through all points, if neighbor in cluster add, if not, new cluster
    3. builds tree and connected leaf data structure.
    4. calculates pocket properties
    5. optional pocket property is the intersection/overlap of the ligand vol.
    (6. doPCA optionally turns on/off the pocket PCA feature computation)'''
    sortedNodes = []
    self.nodeTag = {}  # tag to number of 'valid' neighbors
                       # decrement all neighbors when a node is processed
                       # if tag is 0, not a mouth, if >0, is a mouth
    for aNode in self.gridNodes + self.surfaceNodes:
      if self.nodeClass[aNode] in validCodes and name in self.nodeDist[aNode]:
        sortedNodes.append(aNode)
        goodNeighborCount = 0
        for neighborNode, neighDist in self.nodeNeighbors[aNode]:
          if name in self.nodeDist[neighborNode]:
          #if self.nodeClass[neighborNode] != 1:  # not inside
            goodNeighborCount += 1
        self.nodeTag[aNode] = goodNeighborCount
    sortedNodes.sort(self.__distCompare(name))  # sorts based on distance
    sortedNodes.reverse()  # want largest=biggest travel depth first
    groups = unionfind2.unionFindAttach()
    clusters = unionfind2.unionFindAttach()
    #actually clusters, eventually gets to 1 @ ch
    localMaxima = []
    borders = []
    groupCount = 0
    tm3treeNames = [
        "NAME", "Surface Area", "Volume", "threshold", "meanTD", "maxTD",
        "height", "mean height", "mouths", "Area of Biggest Mouth",
        "Diameter of Biggest Mouth", "Negative Surface Area",
        "Fraction Negative Surface Area", "Apolar Surface Area",
        "Fraction Apolar Surface Area", "Positive Surface Area",
        "Fraction Positive Surface Area", "mean Curvature",
        "mean absolute Curvature", "longest dimension", "middle dimension",
        "short dimension", "Type of Pocket", "Residue Number List",
        "Residue Name List", "Atom Name List", "Mouth Residue Name List"]
    if ligandNodes is not None:
      tm3treeNames.extend([
          "Ligand Intersection", "Ligand Overlap", "Ligand Jaccard",
          "Fraction Ligand Volume"])
    tree = tm3.tmTree(tm3treeNames)
    lastNode = None
    surfNodeToLeaf = {}
    for count in xrange(len(sortedNodes)):
      thisNode = sortedNodes[count]
      thisThresh = self.nodeDist[thisNode][name]
      #print thisNode, thisThresh
      lastNode = thisNode
      #groups.find(thisNode)  # insert always, since new
      #clusters.find(thisNode)
      alreadySeenNeighbors = []
      for neighborNode, neighDist in self.nodeNeighbors[thisNode]:
        if groups.check(neighborNode):  # means neighbor was seen and is valid
          addIt = True  # want only unique groups
          for alreadySeenNeighbor in alreadySeenNeighbors:
            if not clusters.different(alreadySeenNeighbor, neighborNode):
              addIt = False
              break
          if addIt:
            alreadySeenNeighbors.append(neighborNode)
      if len(alreadySeenNeighbors) == 0:
        localMaxima.append(thisNode)
        groupCount += 1  # getting lost somewhere, not all groups accounted
        groups.find(thisNode, set([str(groupCount)]))  # attach name
        clusters.find(thisNode, set([(str(groupCount), ())]))  # attach name
        #print clusters.getAttached(thisNode)
      elif len(alreadySeenNeighbors) == 1:  # new part of old group
        groups.find(thisNode)
        groups.union(thisNode, alreadySeenNeighbors[0])
        clusters.union(thisNode, alreadySeenNeighbors[0])
      elif len(alreadySeenNeighbors) > 1:
        groups.find(thisNode)  # insert always, since new
        groups.union(thisNode, alreadySeenNeighbors[0])  # arbitrary choice
        childrenNames, childrenNodes, chMeshNodes = [], [], []
        for alreadySeenNeighbor in alreadySeenNeighbors:
          chMeshNodes.append(clusters.getList(alreadySeenNeighbor))
          attachedData = clusters.getAttached(alreadySeenNeighbor)
          #print attachedData
          for attachedDatum in attachedData:
            childrenNames.append(attachedDatum[0])
            childrenNodes.append(attachedDatum[1])
        #print childrenNames
        #for each child, we have to get all the data and make the new node
        newChildNodes = []
        for count, aChild in enumerate(childrenNames):
          for aNode in chMeshNodes[count]:  # add to surfnodetoleaf now
            if (3 == self.nodeClass[aNode] or 5 == self.nodeClass[aNode]) and \
                aNode not in surfNodeToLeaf:  # no overwrite, surface
              surfNodeToLeaf[aNode] = aChild
          volume, surfArea, meanTD, mouthLists, nearbyRes, nearbyResNames, \
              nearbyAtoms, biggestMouth, farthestPointMouth, pocketType, \
              mouthResNames, maxTD, negSA, apolarSA, posSA, \
              meanCurv, meanAbsCurv, dimensions, \
              sizeInt, sizeUnion, jaccard, fraclig = \
              self.calculatePocketProperties(
                  chMeshNodes[count], name, pdbData, pointAtomList,
                  ligandNodes, doPCA=doPCA)
          try:
            fractionNeg = negSA/surfArea
            fractionApolar = apolarSA/surfArea
            fractionPos = posSA/surfArea
          except ZeroDivisionError:  # surfarea is 0
            fractionNeg = 0.
            fractionApolar = 0.
            fractionPos = 0.
          tm3nodeProps = [
              aChild, surfArea, volume, thisThresh, meanTD, maxTD,
              maxTD - thisThresh, meanTD-thisThresh, len(mouthLists),
              float(biggestMouth), farthestPointMouth, negSA, fractionNeg,
              apolarSA, fractionApolar, posSA, fractionPos, meanCurv,
              meanAbsCurv, dimensions[0], dimensions[1], dimensions[2],
              pocketType, nearbyRes, nearbyResNames, nearbyAtoms, mouthResNames]
          if ligandNodes is not None:
            tm3nodeProps.extend([sizeInt, sizeUnion, jaccard, fraclig])
          newChildNodes.append(tree.addNode(
              tm3nodeProps, childrenNodes[count]).getId())
        for alreadySeenNeighbor in alreadySeenNeighbors:
          clusters.union(thisNode, alreadySeenNeighbor)
        groupCount += 1
        clusters.clearAttached(thisNode)
        clusters.find(thisNode, set([(str(groupCount), tuple(newChildNodes))]))
        #print clusters.getAttached(thisNode)
        if 3 == self.nodeClass[thisNode] or 5 == self.nodeClass[thisNode]:
          borders.append(thisNode)
      #do neighbor groups leaf attachment now
      thisId = list(groups.getAttached(thisNode))[0]  # only one
      for neighbor, neighDist in self.nodeNeighbors[thisNode]:
        if neighbor in self.nodeTag:
          self.nodeTag[neighbor] -= 1
    #print len(clusters.toLists())   # debugging
    #set the root and all that needs to be done to finalize the tree
    attachedData = clusters.getAttached(lastNode)
    childrenNodes = []
    for attachedDatum in attachedData:
      newName = attachedDatum[0]
      childrenNodes.extend(list(attachedDatum[1]))
    allNodes = clusters.getList(lastNode)
    #these are disabled and should only be used to check and see if cavities
    #are working correctly.... next 4 lines
    #clustersLeft = [] #want to process cavities and insert them into the tree
    #for clusterList in clusters.toLists():
    #  if allNodes[0] not in clusterList:  # as good as doing equality checking
    #    clustersLeft.append(clusterList)
    #print len(clustersLeft) #cavity debugging
    volume, surfArea, meanTD, mouthLists, nearbyRes, nearbyResNames, \
        nearbyAtoms, biggestMouth, farthestPointMouth, pocketType, \
        mouthResNames, maxTD, negSA, apolarSA, posSA, \
        meanCurv, meanAbsCurv, dimensions, sizeInt, sizeUnion, jaccard, \
        fraclig = self.calculatePocketProperties(
            allNodes, name, pdbData, pointAtomList, ligandNodes, doPCA=doPCA)
    thisThresh = self.nodeDist[lastNode][name]
    try:
      fractionNeg = negSA/surfArea
      fractionApolar = apolarSA/surfArea
      fractionPos = posSA/surfArea
    except ZeroDivisionError:  # surfarea is 0
      fractionNeg = 0.
      fractionApolar = 0.
      fractionPos = 0.
    tm3nodeProps = [
        newName, surfArea, volume, thisThresh, meanTD, maxTD,
        maxTD - thisThresh, meanTD-thisThresh, len(mouthLists),
        float(biggestMouth), farthestPointMouth, negSA, fractionNeg,
        apolarSA, fractionApolar, posSA, fractionPos, meanCurv, meanAbsCurv,
        dimensions[0], dimensions[1], dimensions[2], pocketType,
        nearbyRes, nearbyResNames, nearbyAtoms, mouthResNames]
    if ligandNodes is not None:
      tm3nodeProps.extend([sizeInt, sizeUnion, jaccard, fraclig])
    rootNode = tree.addNode(tm3nodeProps, childrenNodes)
    for aNode in allNodes:
      if (3 == self.nodeClass[aNode] or 5 == self.nodeClass[aNode]) and \
          aNode not in surfNodeToLeaf:  # surface
        surfNodeToLeaf[aNode] = newName
    tree.setRoot(rootNode)
    for listCount, nodeList in enumerate(groups.toLists()):
      listColor = listCount
      if listCount % 2 == 1:
        listColor += 1000
      for aNode in nodeList:
        if aNode not in borders:
          self.nodeDist[aNode][groupName] = listColor
        else:
          self.nodeDist[aNode][groupName] = -1000
    return localMaxima, borders, tree, surfNodeToLeaf

  def getMinMaxMeanNodes(self, iterableNodes, distanceName):
    '''gets the min, max, and average of a certain distance over the nodes'''
    minD, maxD = None, None
    sumD = 0.
    for node in iterableNodes:
      if distanceName in self.nodeDist[node]:
        distance = self.nodeDist[node][distanceName]
      else:
        distance = 0.  # means it was far outside convex hull
      if minD is None or distance < minD:
        minD = distance
      if maxD is None or distance > maxD:
        maxD = distance
      sumD += distance
    if len(iterableNodes) > 0:
      meanD = sumD / float(len(iterableNodes))
    else:
      meanD = 0.
    return minD, maxD, meanD

  def getWithinNodesNoInside(self, xyzRadiusList):
    '''
    for each ball, finds all nodes not inside the MS that are within the ball
    returns the union of all such nodes
    '''
    returnSet = set()
    for aNodeKey in self.gridNodes + self.surfaceNodes:
      aNode = self.nodeXyz[aNodeKey]
      if self.nodeClass[aNodeKey] != 1:  # is not inside
        for coordRad in xyzRadiusList:
          distanceTo = geometry.distL2(coordRad[:3], aNode)
          if distanceTo < coordRad[3]:
            returnSet.add(aNodeKey)
            break  # don't need to check any more vdw radii
    return returnSet

  def getTracebackSet(self, nodeSet, name):
    '''finds all the nodes that traceback from nodes in this set'''
    returnSet = set()
    listToProcess = nodeSet.copy()
    while len(listToProcess) > 0:
      currentNode = listToProcess.pop()
      if currentNode in self.tracebacks[name]:
        listBack, curNode, dist = self.tracebacks[name][currentNode]
        for nextNode in listBack:
          if nextNode not in listToProcess and nextNode not in returnSet:
            listToProcess.add(nextNode)
      returnSet.add(currentNode)
    return returnSet

  def getFinalTracebackSet(self, nodeSet, name):
    '''finds all the initial nodes that traceback from nodes in this set'''
    returnSet = set()
    listToProcess = nodeSet.copy()
    while len(listToProcess) > 0:
      currentNode = listToProcess.pop()
      if currentNode in self.tracebacks[name]:
        listBack, curNode, dist = self.tracebacks[name][currentNode]
        for nextNode in listBack:
          if nextNode not in listToProcess and nextNode not in returnSet:
            listToProcess.add(nextNode)
      else:  # has no ancestors, must be starting node
        returnSet.add(currentNode)
    return returnSet

  def getBetweenNodes(self):
    '''returns a list of just the between nodes'''
    returnList = []
    for aNode in self.gridNodes:
      if self.nodeClass[aNode] == 2:  # somewhere in the between
        returnList.append(aNode)
    return returnList

  def getSurfaceNodes(self):
    '''returns a list of just the surface nodes'''
    return self.surfaceNodes

  def getSurfaceNodesSorted(self, name):
    '''returns list sorted on some distance (usually travel)'''
    surfNodes = self.getSurfaceNodes()
    surfNodes.sort(self.__distCompare(name))  # sorts based on distance
    surfNodes.reverse()  # want high to low
    return surfNodes

  def countNodes(self):
    return len(self.gridNodes) + len(self.surfaceNodes)

  def __numericCompare(self, x, y):
    '''helper function to do comparisons'''
    if x > y:
      return 1
    elif x == y:
      return 0
    else:  # x < y
      return -1

  def __distCompare(self, name):
    '''helper comparison method based on the name distance of passed in nodes'''
    returnFunct = lambda x, y: self.__numericCompare(
        self.nodeDist[x][name], self.nodeDist[y][name])
    return returnFunct  # yes returns a function since sort() only takes 2 args

  def __shortenPath(self, path):
    '''finds shortcuts in the path'''
    shortPath = []
    partialPath = path[:]  # beginning to end of path
    while len(partialPath) > 0:
      thisNode = partialPath.pop()  # remove end
      shortPath.insert(0, thisNode)  # do this regardless
      if len(partialPath) > 0:
        bestShortcut, shortcutIndex = False, len(partialPath)
        for neighbor, neighDist in self.nodeNeighbors[thisNode]:
          if neighbor != partialPath[-1]:  # one must be, but we ignore it
            if neighbor in partialPath:  # found a shortcut
              lastIndex = partialPath.index(neighbor)
              if lastIndex < shortcutIndex:
                shortcutIndex, bestShortcut = lastIndex, neighbor
        if bestShortcut:
          partialPath = partialPath[:shortcutIndex+1]
    return shortPath

  def __findPathOutBranchBound(self, startNode, centerName, excludePts):
    '''uses branch and bound to find a path from start to outside convexhull'''
    visitedNodes = set(excludePts)  # all the things seen or not allowed
    tree = {}  # empty tree, goes up
    nodesLeft, valuesLeft = [startNode], [self.nodeDist[startNode][centerName]]
    while len(nodesLeft) > 0:
      thisIndex = valuesLeft.index(max(valuesLeft))  # expand higher first
      thisNode = nodesLeft.pop(thisIndex)
      thisValue = valuesLeft.pop(thisIndex)
      visitedNodes.add(thisNode)
      for neighbor, neighDist in self.nodeNeighbors[thisNode]:
        if neighbor not in visitedNodes:  # don't want loops
          if neighbor not in nodesLeft:  # already on the stack
            if self.nodeClass[neighbor] != 1:
              #don't care, don't go through inside
              if self.nodeClass[neighbor] == 0:  # then we're done, outside
                path = [thisNode, neighbor]
                value = min(
                    self.nodeDist[thisNode][centerName],
                    self.nodeDist[neighbor][centerName])
                while path[0] != startNode:
                  path.insert(0, tree[path[0]])
                  value = min(value, self.nodeDist[path[0]][centerName])
                shortPath = self.__shortenPath(path)  # shorten path if possible
                #print len(path), len(shortPath)
                return shortPath, value
              else:  # add to the list of things to check
                nodesLeft.append(neighbor)
                valuesLeft.append(self.nodeDist[neighbor][centerName])
                tree[neighbor] = thisNode
    #failure...
    return False, False

  def getNeighborsNotInside(self, iterableNodes, excludeSet=False):
    '''returns neighbors of iterableNodes not in excludeSet or iterableNodes,
    also does not include neighbors on the inside of the surface'''
    neighborNodes = set()
    for node in iter(iterableNodes):
      for neighborNode, neighDist in self.nodeNeighbors[node]:
        if neighborNode not in iterableNodes and (
            (not excludeSet) or neighborNode not in excludeSet):
          if self.nodeClass[neighborNode] != 1:  # .isInside():
            neighborNodes.add(neighborNode)
    return neighborNodes

  def getNeighborsOutsideOrBetween(self, iterableNodes, excludeSet=False):
    '''returns neighbors of iterableNodes not in excludeSet or iterableNodes,
    also does not include neighbors on the inside of the surface'''
    neighborNodes = set()
    for node in iter(iterableNodes):
      for neighborNode, neighDist in self.nodeNeighbors[node]:
        if neighborNode not in iterableNodes and (
            (not excludeSet) or neighborNode not in excludeSet):
          if 0 == self.nodeClass[neighborNode] or \
              2 == self.nodeClass[neighborNode]:  # outside or between
            neighborNodes.add(neighborNode)
    return neighborNodes

  def getNeighborsSurface(self, iterableNodes, excludeSet=False):
    '''returns neighbors of iterableNodes not in excludeSet or iterableNodes,
    only returns surface nodes'''
    neighborNodes = set()
    for node in iter(iterableNodes):
      for neighborNode, neighDist in self.nodeNeighbors[node]:
        if neighborNode not in iterableNodes and (
            (not excludeSet) or neighborNode not in excludeSet):
          if self.nodeClass[neighborNode] in [3, 5]:
            neighborNodes.add(neighborNode)
    return neighborNodes

  def getNeighborsOutsideOrBetween(self, iterableNodes, excludeSet=False):
    '''returns neighbors of iterableNodes not in excludeSet or iterableNodes,
    only returns nodes outside or between'''
    neighborNodes = set()
    for node in iter(iterableNodes):
      for neighborNode, neighDist in self.nodeNeighbors[node]:
        if neighborNode not in iterableNodes and (
            (not excludeSet) or neighborNode not in excludeSet):
          if self.nodeClass[neighborNode] in [0, 2]:
            neighborNodes.add(neighborNode)
    return neighborNodes

  def containsOutside(self, iterableNodes):
    '''returns true if any are outside'''
    for node in iter(iterableNodes):
      if self.nodeClass[node] == 0:
        return True
    return False

  def getJustSurfaceNodes(self, iterableNodes):
    '''returns a set of just the surface nodes in iterablenodes'''
    surfaceNodes = set()
    for node in iter(iterableNodes):
      if self.nodeClass[node] in [3, 5]:
        surfaceNodes.add(node)
    return surfaceNodes

  def findDisjointUnions(self, iterableNodes):
    '''takes an iterable list of nodes, breaks into disjoint sets based on
    neighbors and returns a list of these sets'''
    discs = unionfind2.unionFind()
    for node in iter(iterableNodes):
      discs.find(node)  # make sure it is added even if no neighbors
      for neighborNode, neighDist in self.nodeNeighbors[node]:
        if neighborNode in iterableNodes:  # otherwise don't care
          discs.union(node, neighborNode)
    return discs.toLists()

  def __addPointNeighbors(self, pointNeighbors, pointOne, pointTwo):
    '''adds to the tree, which is a two-way connected dictionary'''
    if pointOne not in pointNeighbors:            # initialize list
      pointNeighbors[pointOne] = []
    if pointTwo not in pointNeighbors[pointOne]:  # no duplicates
      pointNeighbors[pointOne].append(pointTwo)
    if pointTwo not in pointNeighbors:            # initialize list
      pointNeighbors[pointTwo] = []
    if pointOne not in pointNeighbors[pointTwo]:  # no duplicates
      pointNeighbors[pointTwo].append(pointOne)

  def __calculateMinPointToDist(self, point, iterableSet):
    '''iterates over whole set, comparing to point, returns smallest'''
    minDist = 10000000000.
    for otherPt in iter(iterableSet):
      minDist = min(minDist, geometry.distL2Squared(
          self.nodeXyz[point], self.nodeXyz[otherPt]))
    return math.sqrt(minDist)

  def __decideWhichCloser(self, thisPt, loopPts):
    '''checks both loops, sees which the point is closer to'''
    minDists, lower = [False, False], 1
    for index in xrange(len(loopPts)):
      minDists[index] = self.__calculateMinPointToDist(thisPt, loopPts[index])
    if minDists[0] < minDists[1]:
      lower = 0
    return lower, minDists[lower]  # in case you want dist too

  def __divideTwoCloseGroups(self, loopPtsOne, loopPtsTwo):
    '''this is where the plugs are made. algorithm splits all points into 2
    sets based on which loop it is closer to. plugs are adjacent to each other
    and adjacent (eventually) to a loop point.'''
    loopPts = (loopPtsOne, loopPtsTwo)
    allPts, triedPts = set(), set()
    allPts.update(self.getNeighborsOutsideOrBetween(loopPtsOne))
    allPts.update(self.getNeighborsOutsideOrBetween(loopPtsTwo))
    pointGroups = (set(loopPtsOne), set(loopPtsTwo))
    for pointGroup in pointGroups:
      allPts.difference_update(pointGroup)
    removePts = set()
    for point in allPts:
      if self.nodeClass[point] not in [0, 2]:  # not out or between
        removePts.add(point)
    allPts.difference_update(removePts)
    #don't want neighboring points on the surface
    #only points within maxDist of any loop point are considered, this ensures
    # a complete plug across the whole opening
    while len(allPts) > 0:
      thisPt = allPts.pop()
      triedPts.add(thisPt)
      lower, distClose = self.__decideWhichCloser(thisPt, loopPts)
      other = (lower + 1) % 2
      neighborOpposite = False
      neighborsNotInside = self.getNeighborsOutsideOrBetween([thisPt])
      for neighborPt in neighborsNotInside:
        if neighborPt in pointGroups[other]:
          neighborOpposite = True
          break
      if not neighborOpposite:  # yet
        for neighborPt in neighborsNotInside:
          if neighborPt in allPts or neighborPt in triedPts:
            #otherwise we don't know if it is a neighbor of a neighbor yet
            if neighborPt not in pointGroups[lower]:
              neighborLower, neighborDist = self.__decideWhichCloser(
                  neighborPt, loopPts)
              if neighborLower == other:
                #have to check and make sure this neighbor is adjacent to
                # other plug
                neighborNeighbors = self.getNeighborsNotInside([neighborPt])
                for neighborNeighborCheck in neighborNeighbors:
                  if neighborNeighborCheck in pointGroups[other]:
                    #has to be adjacent to the correct plug
                    neighborOpposite = True
                    break  # out of 4 loop 3 lines previous
                if neighborOpposite:
                  break  # out of 4 loop 12 lines previous
      if neighborOpposite:  # now take action
        pointGroups[lower].add(thisPt)
        for neighborPt in neighborsNotInside:
          if neighborPt not in pointGroups[0] and \
              neighborPt not in pointGroups[1]:
            if self.nodeClass[neighborPt] in [0, 2]:  # out or betw
              allPts.add(neighborPt)
    #print len(pointGroups[0]), len(pointGroups[1])
    return pointGroups

  def findOnePlug(
      self, indexLoop, regLoopPtList, debugOut, outFileName, centerName):
    '''finds a plug for one regular loop'''
    ringNeighbors = {}  # dictionary keyed on disc to lists of discs
    ringToDisc = {}  # mapping from rings to discs
    ringToExclude = {}
    ringToCenterName = {}
    possibleHoleStarts = []
    rings = [False, False]  # two lists of sets
    discs = [False, False]  # two lists of sets
    for indexSide, regLoopPtListSide in enumerate(regLoopPtList):
      ringThisSide = set(regLoopPtListSide)
      rings[indexSide] = frozenset(ringThisSide)
    for side in xrange(len(discs)):
      other = (side + 1) % 2
      if 1 == other:  # only these some things once
        groups = self.__divideTwoCloseGroups(rings[side], rings[other])
        wholeDiscSet = set()  # cumulative over matching pairs of discs
        for indexGroup, group in enumerate(groups):
          wholeDiscSet.update(group)
          discThisGroup = frozenset(group)
          ringToDisc[rings[indexGroup]] = discThisGroup
          discs[indexGroup] = discThisGroup
          if debugOut:
            pts = []  # here to end of loop is debugging
            for node in iter(group):
              pts.append(self.nodeXyz[node])
            if indexGroup == 0:  # debugging
              tstdebug.pointDebug(
                  pts, filename=outFileName + ".plug." + str(indexLoop) +
                  "." + str(indexGroup) + ".py", mainColor=(.1, .9, .1))
            else:
              tstdebug.pointDebug(
                  pts, filename=outFileName + ".plug." + str(indexLoop) +
                  "." + str(indexGroup) + ".py", mainColor=(.9, .1, .1))
        possibleHoleStarts.append(
            self.findMaxima(centerName, wholeDiscSet)[0][0])
        ringToExclude[rings[side]] = ringToDisc[rings[other]]
        #don't go through other side
        ringToExclude[rings[other]] = ringToDisc[rings[side]]
        #don't go through other side
      ringNeighbors[rings[side]] = [rings[other]]  # list of sets
    return ringNeighbors, ringToDisc, ringToExclude, ringToCenterName, \
        possibleHoleStarts

  def findOnePlugMaxima(
      self, ringCount, ring, ringToDiscAll, ringToExcludeAll,
      debugOut, outFileName, centerName):
    '''finds and returns the maxima'''
    disc = ringToDiscAll[ring]
    excludeSet = ringToExcludeAll[ring]
    #print len(ring), len(disc),
    maximaDisc, maximaDist = self.findMaxima(centerName, disc)
    if debugOut:
      pts = []  # here to end of loop is debugging
      for node in iter(maximaDisc):
        pts.append(self.nodeXyz[node])
        #print node.classification
      tstdebug.pointDebug(
          pts, filename=outFileName + ".plug.localmaxima." + str(ringCount) +
          ".py", mainColor=(.1, .3, .9))
    return maximaDisc, maximaDist

  def growOneMaxima(
      self, ringToMaximaAll, ring, maximaDisc, centerName,
      ringToExcludeAll, pointNeighbors, outsidePoints,
      ringNeighborsAll, ringToDiscAll):
    '''grows out from one maxima, adds best path, doesn't return since stuff
    modified in place (pointNeighbors, ringToMaximaAll, outsidePoints'''
    excludeSet = ringToExcludeAll[ring]
    ringToMaximaAll[ring], paths, bestPath = [], [], 0
    for index, maximaDiscTry in enumerate(maximaDisc):
      path, value = self.__findPathOutBranchBound(
          maximaDiscTry, centerName, excludeSet)
      #print len(path)
      if path and len(path) > 1:
        paths.append(path)
        if len(path) < len(paths[bestPath]):
          bestPath = len(paths) - 1
    if len(paths) > 0:
      ringToMaximaAll[ring].append(paths[bestPath][0])   # the maxima used
      for pathIndex in xrange(len(paths[bestPath])-1):
        #does pairs, don't do last one
        otherIndex = pathIndex + 1
        self.__addPointNeighbors(
            pointNeighbors,
            paths[bestPath][pathIndex], paths[bestPath][otherIndex])
      if paths[bestPath][-1] not in outsidePoints:  # want unique points
        outsidePoints.append(paths[bestPath][-1])
    else:  # len(paths) == 0:  # no paths found!!!!
      ringToMaximaAll[ring] = [maximaDisc[0]]  # arbitrary but necessary
      #print indexMin, maximaDisc, ring, paths
    #print outsidePoints[-1].isOutside()  # checks... shouldn't be necessary
    #connect maxima (starting points)
    biggestMaxima = False
    if ring in ringNeighborsAll:  # do after only both have maxima
      otherRing = ringNeighborsAll[ring][0]
      if otherRing in ringToMaximaAll:
        discs = [ringToDiscAll[ring], ringToDiscAll[ring]]
        #each maxima on each side needs to map to all on opposite side
        for maximaOneSide in ringToMaximaAll[ring]:
          for maximaOtherSide in ringToMaximaAll[otherRing]:
            maximas = [maximaOneSide, maximaOtherSide]
            if self.nodeDist[maximas[0]][centerName] > \
                self.nodeDist[maximas[1]][centerName]:
              biggestMaxima = maximas[0]
            else:
              biggestMaxima = maximas[1]
            #print maximas,
            connPath, value = self.__connectTwoMaxima(
                centerName, discs, maximas)
            #print value
            for pathIndex in xrange(len(connPath) - 1):
              #does pairs, don't do last one
              otherIndex = pathIndex + 1
              self.__addPointNeighbors(
                  pointNeighbors, connPath[pathIndex], connPath[otherIndex])
    #return almost nothing, modified in place
    return biggestMaxima

  def findMaxima(self, name, iterableNodes):
    '''finds the maxima among the iterableNodes according to name'''
    maxSet, maxVal = set(), -1.
    for node in iter(iterableNodes):
      if self.nodeDist[node][name] == maxVal:
        maxSet.add(node)
      elif self.nodeDist[node][name] > maxVal:
        maxSet = set([node])  # start over
        maxVal = self.nodeDist[node][name]
    #change set to list here...
    returnMaxList = []
    while len(maxSet) > 0:
      returnMaxList.append(maxSet.pop())
    return tuple(returnMaxList), maxVal

  def __connectTwoMaxima(self, centerName, discs, maximas):
    '''connects two maxima together, discs and maximas are tuples.
    uses branch and bound to find the best path with the lowest dist on it.
    doesn't allow path to wander very far from straightline path, nodes should
    be close and is better than restricting to discs which are sometimes
    not completely connected with each other.'''
    visitedNodes = set()  # all the things seen
    tree = {}  # empty tree, goes up
    maxDistToCheck = geometry.distL2(
        self.nodeXyz[maximas[0]], self.nodeXyz[maximas[1]]) * 1.5
    #could be lower but be safe
    #print maxDistToCheck
    nodesLeft, valuesLeft = [
        maximas[0]], [self.nodeDist[maximas[0]][centerName]]
    while len(nodesLeft) > 0:
      thisIndex = valuesLeft.index(max(valuesLeft))  # expand higher first
      thisNode = nodesLeft.pop(thisIndex)
      thisValue = valuesLeft.pop(thisIndex)
      visitedNodes.add(thisNode)
      for neighbor, neighborDist in self.nodeNeighbors[thisNode]:
        if not discs or neighbor in discs[0] or neighbor in discs[1]:
          if neighbor not in visitedNodes:  # don't want loops
            if neighbor not in nodesLeft:  # already on the stack
              if neighbor == maximas[1]:  # then we're done
                path = [thisNode, neighbor]
                value = min(
                    self.nodeDist[thisNode][centerName],
                    self.nodeDist[neighbor][centerName])
                while path[0] != maximas[0]:
                  path.insert(0, tree[path[0]])
                  value = min(value, self.nodeDist[path[0]][centerName])
                return path, value
              elif 1 != self.nodeClass[neighbor]:  # not inside
                #add to the list of things to check...if dist okay
                thisDist = geometry.distL2(
                    self.nodeXyz[maximas[1]], self.nodeXyz[neighbor])
                if thisDist < maxDistToCheck:
                  nodesLeft.append(neighbor)
                  valuesLeft.append(self.nodeDist[neighbor][centerName])
                  tree[neighbor] = thisNode
                else:  # so we don't keep trying this one
                  visitedNodes.add(neighbor)
    #if we get here it means we've tried and tried but can't get there from here
    #possibly retry without restricting to discs
    if discs:
      #print "trying again"
      return self.__connectTwoMaxima(centerName, False, maximas)
    else:  # really give up
      #print "giving up, giving in"
      return False, False

  def __findPathPointToPoint(
      self, pointStart, pointEnd, pointNeighbors, centerName):
    '''uses depth first search to find all non-loop paths from start to end'''
    paths = [[pointStart]]
    completePaths = []
    #seenPts = set([pointStart])
    while len(paths) > 0:
      curPath = paths.pop()
      curEnd = curPath[-1]
      neighbors = pointNeighbors[curEnd]
      neighbors.sort(self.__distCompare(centerName))  # sorts based on distance
      #neighbors is sorted low to high now
      for neighbor in neighbors:
        if neighbor == pointEnd:  # found one path
          completePath = curPath[:]
          completePath.append(neighbor)
          completePaths.append(completePath)  # and continue
        elif neighbor not in curPath:  # no loops
          newPath = curPath[:]
          #seenPts.add(neighbor)
          newPath.append(neighbor)
          paths.append(newPath)    # highs put on last... so depth first
    return completePaths  # return the list of paths found

  def __findPathPointToMultiple(
      self, pointStart, pointsEnd, pointNeighbors, centerName, plugPoints):
    '''uses depth first search to find all non-loop paths from start to ends'''
    paths = [[pointStart]]
    completePaths = []
    #seenPts = set([pointStart])
    while len(paths) > 0:
      curPath = paths.pop()
      curEnd = curPath[-1]
      neighbors = pointNeighbors[curEnd]
      neighbors.sort(self.__distCompare(centerName))  # sorts based on distance
      #neighbors is sorted low to high now
      for neighbor in neighbors:
        if neighbor in pointsEnd:  # found one path
          if neighbor != pointStart or len(curPath) > 2:
            completePath = curPath[:]
            completePath.append(neighbor)
            completePaths.append(completePath)  # and continue
        elif neighbor not in curPath:  # no cycles
          newPath = curPath[:]
          newPath.append(neighbor)
          paths.append(newPath)     # highs put on last... so depth first
    return completePaths  # return the list of paths found

  def __getPathPlugs(self, plugPoints, pathPoints):
    '''helper to find the plugs each path passes through'''
    plugs = "-"
    for indexPlug, plugPoint in enumerate(plugPoints):
      for pathPoint in pathPoints:
        if plugPoint.reference == pathPoint.reference:
          plugs += str(indexPlug) + "-"
          break
    return plugs

  def getPaths(self, centerName, pointNeighbors, outsidePoints, plugPoints):
    '''finds all paths from each outside point to all others,
       also tracks which plugs each path passes through,
       reports shortest path unique wrt begin, end, plugs'''
    paths = []
    #force paths to be unique wrt start, end, and plugPoints passed through
    for indexOne, pointOne in enumerate(outsidePoints):
      otherPoints = outsidePoints[indexOne:]  # include starting point
      if len(otherPoints) > 0:
        pathEndDict = {}
        allPaths = self.__findPathPointToMultiple(
            pointOne, otherPoints, pointNeighbors, centerName, plugPoints)
        if allPaths and len(allPaths) > 0:  # otherwise no paths found
          for path in allPaths:
            if path and len(path) > 0:  # paranoia
              #print path
              indexTwo = outsidePoints.index(path[-1])
              if indexTwo not in pathEndDict:
                pathEndDict[indexTwo] = []
              pathEndDict[indexTwo].append(path)
        for indexTwo in pathEndDict.keys():
          possPaths = pathEndDict[indexTwo]
          if 1 == len(possPaths):
            nodePath = self.makeNodes(possPaths[0])
            plugs = self.__getPathPlugs(plugPoints, nodePath)
            paths.append((indexOne, indexTwo, plugs, nodePath))
          elif 1 < len(possPaths):
            pathPlugDict = {}
            for path in possPaths:
              pathNodes = self.makeNodes(path)
              plugs = self.__getPathPlugs(plugPoints, pathNodes)
              if plugs not in pathPlugDict:
                pathPlugDict[plugs] = []
              pathPlugDict[plugs].append(path)
            for plugs in pathPlugDict.keys():
              possPlugPaths = pathPlugDict[plugs]
              if 1 == len(possPlugPaths):
                pathNodes = self.makeNodes(possPlugPaths[0])
                paths.append((indexOne, indexTwo, plugs, pathNodes))
              elif 1 < len(possPlugPaths):
                #pick min-length path
                minLengthPathIndex = 0
                for pathIndexHere, path in enumerate(possPlugPaths):
                  if len(path) < len(possPlugPaths[minLengthPathIndex]):
                    minLengthPathIndex = pathIndexHere
                smallPath = possPlugPaths[minLengthPathIndex]
                pathNodes = self.makeNodes(smallPath)
                paths.append((indexOne, indexTwo, plugs, pathNodes))
    return paths

  def makeNodesDict(self, nodeRefDict):
    '''makes a dictionary of nodes to nodes from a dict of refs to refs'''
    nodeDict = {}
    for refKey, refValueList in nodeRefDict.iteritems():
      if refKey not in self.nodeNode:
        nodeTemp = node(
            refKey, self.nodeXyz[refKey],
            self.nodeOffGridXyz[refKey], self.nodeDist[refKey])
        self.nodeNode[refKey] = nodeTemp
      nodeKey = self.nodeNode[refKey]
      nodeList = []
      for nodeRef in refValueList:
        if nodeRef not in self.nodeNode:
          self.nodeNode[nodeRef] = node(
              nodeRef, self.nodeXyz[nodeRef],
              self.nodeOffGridXyz[nodeRef], self.nodeDist[nodeRef])
        nodeList.append(self.nodeNode[nodeRef])
      nodeDict[nodeKey] = nodeList
    return nodeDict

  def makeNodes(self, nodeRefList):
    '''makes a list of node objects from a list of references.'''
    nodeList = []
    for nodeRef in nodeRefList:
      if nodeRef not in self.nodeNode:
        self.nodeNode[nodeRef] = node(
            nodeRef, self.nodeXyz[nodeRef],
            self.nodeOffGridXyz[nodeRef], self.nodeDist[nodeRef])
      nodeList.append(self.nodeNode[nodeRef])
    return nodeList

  def debugMesh(self, filename):
    for code in range(4):  # all 4 codes
      fileTemp = open(filename + "." + str(code) + ".py", 'w')
      fileTemp.write(
          "from pymol.cgo import * \nfrom pymol import cmd\nsurfObj = [ ")
      for aNode in self.gridNodes.values() + self.surfaceNodes.values():
        if code == self.nodeClass[aNode]:
          fileTemp.write("COLOR, ")
          fileTemp.write("0.1, ")
          fileTemp.write("0.95, ")
          fileTemp.write("0.95, ")
          fileTemp.write("SPHERE, ")
          for coord in self.nodeXyz[aNode]:
            fileTemp.write(str(coord) + ", ")
          fileTemp.write("0.2, ")
      fileTemp.write(" ]\n")
      fileTemp.write(
          "cmd.load_cgo(surfObj,'grid" + str(code) + filename + "')\n")
      fileTemp.close()

class meshFromSpheres(mesh):
  '''a class that creates a new mesh from a given list of spheres,
  can be made to be a VDW surface or solvent accessible surface, etc.'''

  def __init__(self, grid, mins, maxs, gridSize, spheres):
    '''creates a mesh that lies on a previously specified grid'''
    self.gridNodes = []   # dictionaries of nodes, keyed on gridIndex or ptIndex
    self.surfaceNodes = []

  def __init__(self, spheres, gridSpacing=1.0, outputName="debug.spheres.py"):
    '''creates a mesh that lies on a new grid with given spacing'''
    self.surfaceNodes = []
    self.surfacePoints = []
    self.nodeDist = {}  # moving data storage out of node class into dicts
    self.nodeNeighbors = {}
    self.nodeClass = {}
    self.nodeXyz = {}
    self.nodeOffGridXyz = {}
    self.nodeNode = {}
    mins, maxs = geometry.findMinsMaxsSpheres(spheres)  # call geometry method
    #add 1.5*gridSpacing to min and max, ensure gridpoints outside surface
    for coord in range(3):
      mins[coord] = mins[coord] - gridSpacing * 1.5
      maxs[coord] = maxs[coord] + gridSpacing * 1.5
    #add to maxs to bring all dimension differences to a multiple of spacing
    for coord in range(3):
      difference = maxs[coord] - mins[coord]
      maxs[coord] += difference % gridSpacing
    self.__createEmptyGridNodes(mins, maxs, gridSpacing)
    #gridNodes has now been modified, mins maxs and spacing are now saved
    self.__createSurfaceNodes(spheres)
    tempPts = []
    for node in self.surfaceNodes:
      tempPts.append(self.nodeXyz[node])
    tstdebug.pointDebug(
        tempPts, filename=outputName, mainColor=(.9, .4, .1),
        radius=gridSpacing/10.)

  def __createEmptyGridNodes(self, mins, maxs, spacing, classification=0):
    '''makes a grid full of nodes, set to default outside'''
    self.gridNodes = []   # dictionaries of nodes, keyed on gridIndex or ptIndex
    lengths = []
    for coord in range(3):  # get the length in each dimension
      lengths.append(1+int((maxs[coord]-mins[coord])/spacing))
    for xIndex in xrange(lengths[0]):
      for yIndex in xrange(lengths[1]):
        for zIndex in xrange(lengths[2]):
          indices = (xIndex, yIndex, zIndex)
          point = []
          for coord in xrange(3):
            point.append(indices[coord] * spacing + mins[coord])
          self.gridNodes.append(indices)
          self.nodeXyz[indices] = point
          self.nodeOffGridXyz[indices] = point
          self.nodeDist[indices] = {}
          self.nodeNeighbors[indices] = []
          self.nodeClass[indices] = 0
    self.lengths = lengths  # save this too
    self.mins = mins  # save these three things, they will no longer change
    self.maxs = maxs
    self.spacing = spacing
    #no need to return, as gridNodes have been modified

  def __updateLineBreaks(self, linebreaks, coords):
    '''linebreaks is a list of in/out vertices, coords is a single pair,
    update and return new linebreaks'''
    if 0 == len(linebreaks):
      return coords  # easy case
    else:  # now harder case, find place to put
      newLB = []
      lineBreakIndex = 0
      while lineBreakIndex < len(linebreaks) and \
          linebreaks[lineBreakIndex] < coords[0]:  # find first first
        newLB.append(linebreaks[lineBreakIndex])
        lineBreakIndex += 1
      if lineBreakIndex == len(linebreaks):
        #exhausted original list, tag to end
        newLB.extend(coords)
        return newLB
      #otherwise keep going,
      if 1 == lineBreakIndex % 2:  # means we ignore the first, superseded
        pass
      else:  # new, so insert
        newLB.append(coords[0])
      while lineBreakIndex < len(linebreaks) and \
          linebreaks[lineBreakIndex] < coords[1]:  # find second now
        lineBreakIndex += 1  # not adding to lb
      if lineBreakIndex == len(linebreaks):
        #exhausted original list, tag to end
        newLB.append(coords[1])
        return newLB
      if 1 == len(newLB) % 2 and 0 == lineBreakIndex % 2:
        newLB.append(coords[1])
      while lineBreakIndex < len(linebreaks):
        newLB.append(linebreaks[lineBreakIndex])
        lineBreakIndex += 1
      return newLB  # either way return the new one

  def __createSurfNodesOneDir(self, indexOne, indexTwo, spheres, axis):
    '''creates surface nodes along one line in the axis direction'''
    if 2 == axis:
      minIndices = (indexOne, indexTwo, 0)
      maxIndices = (indexOne, indexTwo, self.lengths[axis] - 1)
    elif 1 == axis:
      minIndices = (indexOne, 0, indexTwo)
      maxIndices = (indexOne, self.lengths[axis] - 1, indexTwo)
    elif 0 == axis:
      minIndices = (0, indexOne, indexTwo)
      maxIndices = (self.lengths[axis] - 1, indexOne, indexTwo)
    extremeNodes = (minIndices, maxIndices)
    extrema = []
    for node in extremeNodes:
      extrema.append(self.nodeXyz[node])
    #algorithm works by dividing line into inside/outside junctions
    #for now assume spheres haven't been checked to make sure they intersect
    lineBreaks = []
    for sphere in spheres:
      points = geometry.lineSphereIntersection(extrema[0], extrema[1], sphere)
      if points:  # failure is indicated by False
        axisCoords = []
        for point in points:
          axisCoords.append(point[axis])
        newLB = self.__updateLineBreaks(lineBreaks, axisCoords)
        #print lineBreaks, axisCoords, newLB #debugging
        #print len(lineBreaks)
        lineBreaks = newLB
    #linebreaks containts axis coords for surface nodes
    lineBreakIndex, lastGridNode = 0, False
    for axisIndex in xrange(self.lengths[axis]-1):
      if 2 == axis:
        curIndices = (indexOne, indexTwo, axisIndex)
      elif 1 == axis:
        curIndices = (indexOne, axisIndex, indexTwo)
      elif 0 == axis:
        curIndices = (axisIndex, indexOne, indexTwo)
      curGridNode = curIndices
      if lineBreakIndex < len(lineBreaks) and \
          self.nodeXyz[curGridNode][axis] > lineBreaks[lineBreakIndex]:
        #here we make a surface node and connect it to the gridnodes
        pointIndex = len(self.surfacePoints) + 1
        if 2 == axis:
          surfacePoint = [
              self.nodeXyz[curGridNode][0], self.nodeXyz[curGridNode][1],
              lineBreaks[lineBreakIndex]]
        elif 1 == axis:
          surfacePoint = [
              self.nodeXyz[curGridNode][0], lineBreaks[lineBreakIndex],
              self.nodeXyz[curGridNode][2]]
        elif 0 == axis:
          surfacePoint = [
              lineBreaks[lineBreakIndex], self.nodeXyz[curGridNode][1],
              self.nodeXyz[curGridNode][2]]
        pointListRef = [pointIndex]
        pointListRef.extend(surfacePoint)
        self.surfacePoints.append(pointListRef)
        newSurfNode = surfaceNode(pointIndex, surfacePoint)
        self.surfaceNodes[pointIndex] = newSurfNode
        lastDist = geometry.distL2(self.nodeXyz[lastGridNode], surfacePoint)
        curDist = geometry.distL2(self.nodeXyz[curGridNode], surfacePoint)
        self.connectNeighbors(curGridNode, newSurfNode, curDist)
        self.connectNeighbors(lastGridNode, newSurfNode, lastDist)
        if 0 == lineBreakIndex % 2:  # odd
          self.nodeClass[curGridNode] = 1
        lineBreakIndex += 1
      elif lastGridNode:
        thisDist = geometry.distL2(
            self.nodeXyz[lastGridNode], self.nodeXyz[curGridNode])
        self.connectNeighbors(curGridNode, lastGridNode, thisDist)
        self.nodeClass[curGridNode] = self.nodeClass[lastGridNode]
      lastGridNode = curGridNode

  def __createSurfaceNodes(self, spheres):
    '''creates surface nodes, sets grid nodes to interior if inside sphere,
    connects surface nodes to grid nodes and grid nodes to each other'''
    for eachX in xrange(self.lengths[0]):
      for eachY in xrange(self.lengths[1]):
        self.__createSurfNodesOneDir(eachX, eachY, spheres, 2)
    for eachX in xrange(self.lengths[0]):
      for eachZ in xrange(self.lengths[2]):
        self.__createSurfNodesOneDir(eachX, eachZ, spheres, 1)
    for eachY in xrange(self.lengths[1]):
      for eachZ in xrange(self.lengths[2]):
        self.__createSurfNodesOneDir(eachY, eachZ, spheres, 0)

class meshJustSurface(mesh):
  '''mesh for just the surface, no grid ever read in'''

  def __init__(self, pointXYZ, pointNeighbor):
    '''creates a mesh from just the surface'''
    self.gridNodes = []    # this will be empty
    self.surfaceNodes = []
    self.nodeDist = {}  # moving data storage out of node class into dicts
    self.nodeNeighbors = {}
    self.nodeClass = {}
    self.nodeXyz = {}
    self.nodeOffGridXyz = {}
    self.nodeNode = {}
    #want somewhere to store tracebacks information
    self.tracebacks = {}  # dictionary keyed on the name just like distances
    self.surfNodeInitialize(pointXYZ)
    self.surfNodeConnect(pointXYZ, pointNeighbor)
    #now all nodes have been added and neighbors calculated

  def surfNodeInitialize(self, pointXYZ):
    '''now initialize the surface nodes'''
    for point in pointXYZ:
      ptIndex, ptXYZ = int(point[0]), point[1:4]
      self.surfaceNodes.append(ptIndex)
      self.nodeXyz[ptIndex] = ptXYZ
      self.nodeOffGridXyz[ptIndex] = ptXYZ
      self.nodeClass[ptIndex] = 3
      self.nodeDist[ptIndex] = {}
      self.nodeNeighbors[ptIndex] = []

  def surfNodeConnect(self, pointXYZ, pointNeighbor):
    '''now connect up all the surface nodes'''
    for point in pointXYZ:  # assuming no cavities
      ptIndex, ptXYZ = int(point[0]), tuple(point[1:4])
      #now connect surface nodes to each other
      for otherPoint in pointNeighbor[ptIndex-1][2:]:
        otherXYZ = pointXYZ[otherPoint-1][1:4]
        self.connectNeighbor(
            ptIndex, otherPoint, geometry.distL2(ptXYZ, otherXYZ))
