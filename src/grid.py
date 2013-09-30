#ryangc@mail.med.upenn.edu ryan.g.coleman@gmail.com
#Ryan G Coleman, Kim Sharp http://crystal.med.upenn.edu

#contains lots of grid primitives, there is not a grid class
import math
from geometry import dist, distL2

def getIndices(mins, gridSize, pt):
  '''helper function to find the box a point is in'''
  xIndex = int(math.floor((pt[0]-mins[0])/gridSize))
  yIndex = int(math.floor((pt[1]-mins[1])/gridSize))
  zIndex = int(math.floor((pt[2]-mins[2])/gridSize))
  #print xIndex,yIndex,zIndex,mins,pt,maxs
  return xIndex,yIndex,zIndex   

def assignAtomDepths(gridD, gridSize, mins, maxs, pdbD, minVal=0.):
  atomDepths = []
  for coord in pdbD.coords:
    gridIndex = getIndices(mins, gridSize, coord)
    #print gridIndex, coord, len(gridD), len(gridD[0]), len(gridD[0][0])
    atomDepths.append(max(minVal, \
                   gridD[gridIndex[0]][gridIndex[1]][gridIndex[2]][0]))
  return atomDepths  

def fillBoxesSets(grid):
  outsideBoxes,insideBoxes = set(), set()
  #useful variables to set
  lenX = len(grid)
  lenY = len(grid[0])
  lenZ = len(grid[0][0])
  for x in xrange(0,lenX):
    for y in xrange(0,lenY):
      for z in xrange(0,lenZ):
        thisBox = grid[x][y][z]
        if thisBox[0] == -1:
          outsideBoxes.update([(x,y,z)])
        elif thisBox[0] == 0:
          insideBoxes.update([(x,y,z)])
          #re-encode box to be large so can tell when 'real' value is present
          newBox = lenX+lenY+lenZ,thisBox[1],thisBox[2],thisBox[3]
          grid[x][y][z] = newBox
        elif thisBox[0] == -2:
          insideBoxes.update([(x,y,z)])
  return outsideBoxes,insideBoxes
  #grid modified in place as well

def getMaximaOnly(grid, maxGreater=0):
  '''returns a copy of the grid, with everything but the maxima removed
  should be run after finalizing grid distances so all non-negative'''
  lens = [len(grid),len(grid[0]),len(grid[0][0])]
  newGrid = []
  for indexX,rowX in enumerate(grid):
    newX = []
    for indexY,rowY in enumerate(rowX):
      newY = []
      for indexZ,entryZ in enumerate(rowY):
        thisBox = grid[indexX][indexY][indexZ][0]
        maxima = maxGreater
        requiredBoxes = 26
        for adjBox in getAllAdjacentBoxes((indexX,indexY,indexZ), \
                                          lens[0], lens[1], lens[2]):
          requiredBoxes -= 1
          if grid[adjBox[0]][adjBox[1]][adjBox[2]][0] >= thisBox:
            maxima -= 1
        if maxima >= 0 and requiredBoxes == 0:
          newY.append(entryZ)
        else:
          newY.append((0,entryZ[1],entryZ[2],entryZ[3]))
      newX.append(newY)
    newGrid.append(newX)
  return newGrid
        
def getMaximaRanking(grid):
  '''returns a copy of the grid, with the value being the number of neighbors
  not greater than this one
  should be run after finalizing grid distances so all non-negative'''
  lens = [len(grid),len(grid[0]),len(grid[0][0])]
  newGrid = []
  for indexX,rowX in enumerate(grid): 
    newX = []
    for indexY,rowY in enumerate(rowX):
      newY = []
      for indexZ,entryZ in enumerate(rowY):
        thisBox = grid[indexX][indexY][indexZ][0]
        maxima = 0
        requiredBoxes = 26
        for adjBox in getAllAdjacentBoxes((indexX,indexY,indexZ), \
                                          lens[0], lens[1], lens[2]):
          requiredBoxes -= 1
          if grid[adjBox[0]][adjBox[1]][adjBox[2]][0] <= thisBox: #ties ties
            maxima += 1   
        if maxima >= 0 and requiredBoxes == 0 and thisBox > 0:
          newY.append((maxima,entryZ[1],entryZ[2],entryZ[3]))
        else:
          newY.append((0,entryZ[1],entryZ[2],entryZ[3]))
      newX.append(newY)
    newGrid.append(newX)
  return newGrid
        
def calculateOffsets(grid1, grid2, gridSize):
  '''figures out the difference in offsets between two grids
  (with the same spacing)'''
  differences = [] 
  for index,value in enumerate(grid1[0][0][0][1:]):
    differences.append(value-grid2[0][0][0][index+1])
  offsets = [int((x/gridSize)) for x in differences]
  #print grid1[0][0][0][1:], grid2[0][0][0][1:],differences, offsets
  return offsets
       
#helper function that calculates L1 distance from end-point to end-point
def calcEdgeGridDist(pt1, pt2, mins, maxs, gridSize, metric='L1'):
  return dist(getIndices(mins,gridSize,pt1), \
              getIndices(mins,gridSize,pt2),metric=metric)

def getAllAdjacentBoxes(curBox, lenX, lenY, lenZ):
  returnVec = []
  returnVec.extend(getAdjacentBoxes(curBox, lenX, lenY, lenZ))
  returnVec.extend(getSideAdjacentBoxes(curBox, lenX, lenY, lenZ))
  returnVec.extend(getCornerAdjacentBoxes(curBox, lenX, lenY, lenZ))
  return returnVec

def getAllAdjacentBoxesOnce(curBox, lenX, lenY, lenZ, extraEdges=False):
  #when called on all curBoxes, only returns each pair once
  returnVec = getAllAdjacentBoxes(curBox, lenX, lenY, lenZ)
  #add distance info
  newReturnVec = []
  for box in returnVec:
    newReturnVec.append((box, distL2(curBox,box)))
  if extraEdges and extraEdges.has_key(curBox):
    for adjBox,adjDist,adjGridDist in extraEdges[curBox]: #tuple unpack
      newReturnVec.append((adjBox,adjDist))
  newVec = [box for box in newReturnVec if curBox[0:2] <= box[0][0:2]]
  return newVec

#helper function, gets up to 6 adjacent boxes
def getAdjacentBoxes(curBox, lenX, lenY, lenZ):
  adjVec = []
  for x in [curBox[0]-1,curBox[0]+1]:
    if x >= 0 and x < lenX:
      adjVec.append((x,curBox[1],curBox[2]))
  for y in [curBox[1]-1,curBox[1]+1]:
    if y >= 0 and y < lenY:
      adjVec.append((curBox[0],y,curBox[2]))
  for z in [curBox[2]-1,curBox[2]+1]:
    if z >= 0 and z < lenZ:
      adjVec.append((curBox[0],curBox[1],z))
  return adjVec

#helper function, gets side adjacent boxes (not directly connected to center)
def getSideAdjacentBoxes(curBox, lenX, lenY, lenZ):
  adjVec = []
  for x in [curBox[0]-1,curBox[0]+1]:
    if x >= 0 and x < lenX:
      for y in [curBox[1]-1,curBox[1]+1]:
        if y >= 0 and y < lenY:
          adjVec.append((x,y,curBox[2]))
      for z in [curBox[2]-1,curBox[2]+1]:
        if z >= 0 and z < lenZ:
          adjVec.append((x,curBox[1],z))
  for y in [curBox[1]-1,curBox[1]+1]:
    if y >= 0 and y < lenY:
      for z in [curBox[2]-1,curBox[2]+1]:
        if z >= 0 and z < lenZ:
          adjVec.append((curBox[0],y,z))
  return adjVec

#helper function, gets 'corner' boxes
def getCornerAdjacentBoxes(curBox, lenX, lenY, lenZ):
  adjVec = []
  for x in [curBox[0]-1,curBox[0]+1]:
    if x >= 0 and x < lenX:
      for y in [curBox[1]-1,curBox[1]+1]:
        if y >= 0 and y < lenY:
          for z in [curBox[2]-1,curBox[2]+1]:
            if z >= 0 and z < lenZ:
              adjVec.append((x,y,z))
  return adjVec

def makeGridFromPhi(phiData, threshold=6.0, inside=-2.0, outside=-1.0):
  newGrid = []
  mins, maxs = phiData.getMinsMaxs()
  gap = 1./phiData.scale
  for x in xrange(phiData.gridDimension):
    newX = []
    for y in xrange(phiData.gridDimension):
      newY = []
      for z in xrange(phiData.gridDimension):
        value = phiData.phiArray[z*(phiData.gridDimension**2) + \
                                   y*phiData.gridDimension + x]
        where = 0.0
        if False == threshold:
          where = value
        else:
          if value < threshold:
            where = outside
          else:
            where = inside
        newTuple = where, mins[0]+(x*gap), mins[1]+(y*gap), mins[2]+(z*gap)
        newY.append(newTuple) #start inside ch, easier to check outsideness
      newX.append(newY)
    newGrid.append(newX)
  return newGrid

#helper routine, makes tuples of encoding, centerX, centerY, centerZ
def makeNewEmptyGrid(mins, maxs, gap, value=0):
  newGrid = []
  for x in xrange(int(math.ceil((maxs[0]-mins[0])/gap))):
    newX = []
    for y in xrange(int(math.ceil((maxs[1]-mins[1])/gap))):
      newY = []
      for z in xrange(int(math.ceil((maxs[2]-mins[2])/gap))):
        newTuple = value, mins[0]+0.5*gap+(x*gap), mins[1]+0.5*gap+(y*gap), \
                                   mins[2]+0.5*gap+(z*gap)
        newY.append(newTuple) #start inside ch, easier to check outsideness
      newX.append(newY)
    newGrid.append(newX)
  return newGrid
  
def findPointMinsMaxs(phiData, pointXYZ, pointList):
  minsPts = pointXYZ[0][1:]
  maxsPts = pointXYZ[0][1:]
  for point in pointList:
    xyz = pointXYZ[point-1][1:]
    for coord in range(3):
      minsPts[coord] = min(minsPts[coord], xyz[coord])   
      maxsPts[coord] = max(maxsPts[coord], xyz[coord])   
  mins, maxs = phiData.getMinsMaxs()
  gap = 1./phiData.scale
  newMins = list(getIndices(mins, gap, minsPts))
  newMaxs = list(getIndices(mins, gap, maxsPts)) #so they initialize to the pts
  return newMins, newMaxs

def makeTrimmedGridFromPhi(phiData, pointXYZ, pointList, \
                           threshold=6.0, inside=-2.0, outside=-1.0, border=2):
  '''makes a trimmed grid from the phi data, hopefully faster than 
  makeGridFromPhi then trimGrid, does not support threshold=False'''
  mins, maxs = phiData.getMinsMaxs()
  gap = 1./phiData.scale
  newMins, newMaxs = findPointMinsMaxs(phiData, pointXYZ, pointList)
  for x in xrange(phiData.gridDimension):
    for y in xrange(phiData.gridDimension):
      for z in xrange(phiData.gridDimension):
        if x < newMins[0] or x > newMaxs[0] or \
           y < newMins[1] or y > newMaxs[1] or \
           z < newMins[2] or z > newMaxs[2]:
          value = phiData.phiArray[z*(phiData.gridDimension**2) + \
                                              y*phiData.gridDimension + x]
          if value >= threshold: #inside the surface
            newMins[0] = min(x, newMins[0])
            newMins[1] = min(y, newMins[1])
            newMins[2] = min(z, newMins[2])
            newMaxs[0] = max(x, newMaxs[0])
            newMaxs[1] = max(y, newMaxs[1])
            newMaxs[2] = max(z, newMaxs[2])
  #add border, careful about current border
  newMins = [max(0,x-border) for x in newMins]
  newMaxs = [min(phiData.gridDimension,x+border) for x in newMaxs]
  newGrid = []
  for x in xrange(phiData.gridDimension):
    newX = []
    for y in xrange(phiData.gridDimension):
      newY = []
      for z in xrange(phiData.gridDimension):
        indices = [x,y,z]
        good = True
        for coord in xrange(3):
          if indices[coord] < newMins[coord] or indices[coord] >=newMaxs[coord]:
            good = False
        if good:
          value = phiData.phiArray[z*(phiData.gridDimension**2) + \
                                           y*phiData.gridDimension + x]
          where = 0.0
          if value < threshold:
            where = outside
          else:
            where = inside
          newTuple = where, mins[0]+(x*gap), mins[1]+(y*gap), mins[2]+(z*gap)
          newY.append(newTuple) #start inside ch, easier to check outsideness
      if len(newY) > 0:
        newX.append(newY)
    if len(newX) > 0:
      newGrid.append(newX)
  lens = [len(newGrid),len(newGrid[0]),len(newGrid[0][0])]
  newMinVals = [x-gap/2. for x in newGrid[0][0][0][1:]]
  newMaxVals = [x+gap/2. for x in newGrid[-1][-1][-1][1:]]
  return newGrid,newMinVals,newMaxVals

def trimGrid(grid, gridSize, pointXYZ, pointList, \
             inside=-2.0, border=1):
  '''trims grid so boundary on each coordinate is only 1 grid cube
  returns grid, mins, maxs of new grid'''
  minsPts = pointXYZ[0][1:]
  maxsPts = pointXYZ[0][1:]
  for point in pointList:
    xyz = pointXYZ[point-1][1:]
    for coord in range(3):
      minsPts[coord] = min(minsPts[coord], xyz[coord])   
      maxsPts[coord] = max(maxsPts[coord], xyz[coord])   
  lens = [len(grid),len(grid[0]),len(grid[0][0])]
  newMins, newMaxs = [10000000.0,10000000.0,100000000.0],[-1000000.,-1000000.,-1000000]
  for indexX,rowX in enumerate(grid):
    for indexY,rowY in enumerate(rowX):
      for indexZ,entryZ in enumerate(rowY):
        indices = [indexX,indexY,indexZ]
        if inside == entryZ[0]:
          for coord in xrange(3): #this makes sure the interior points are inside new grid
            if indices[coord] < newMins[coord]:
              newMins[coord] = indices[coord]
            if indices[coord] > newMaxs[coord]:
              newMaxs[coord] = indices[coord]
        for coord in xrange(3): #this makes sure the points are inside new grid
          if entryZ[coord+1] > minsPts[coord]-gridSize/2.:
            if indices[coord] < newMins[coord]:
              newMins[coord] = indices[coord]
          if entryZ[coord+1] < maxsPts[coord]+gridSize/2.:
            if indices[coord] > newMaxs[coord]:
              newMaxs[coord] = indices[coord]
  #add border, careful about current border
  newMins = [max(0,x-border) for x in newMins]
  newMaxs = [min(lens[coord],newMaxs[coord]+border) for coord in xrange(3)]
  newGrid = []
  for indexX,rowX in enumerate(grid):
    newX = []
    for indexY,rowY in enumerate(rowX):
      newY = []
      for indexZ,entryZ in enumerate(rowY):
        indices = [indexX,indexY,indexZ]
        good = True
        for coord in xrange(3):
          if indices[coord] < newMins[coord] or indices[coord] >=newMaxs[coord]:
            good = False
        if good:
          newY.append(entryZ)
      if len(newY) > 0:
        newX.append(newY)
    if len(newX) > 0:
      newGrid.append(newX)
  newMinVals = [x-gridSize/2. for x in grid[newMins[0]][newMins[1]][newMins[2]][1:]]
  newMaxVals = [x+gridSize/2. for x in grid[newMaxs[0]-1][newMaxs[1]-1][newMaxs[2]-1][1:]]
  return newGrid,newMinVals,newMaxVals

def copyGrid(grid):
  #returns a copy of the grid
  newGrid = []
  for indexX,rowX in enumerate(grid):
    newX = []
    for indexY,rowY in enumerate(rowX):
      newY = []
      for indexZ,entryZ in enumerate(rowY):
        newY.append(entryZ)
      newX.append(newY)
    newGrid.append(newX)
  return newGrid
        
def resetGrid(grid, value=0.):
  '''returns a copy of the grid where all values set to some number'''
  for indexX,rowX in enumerate(grid):
    for indexY,rowY in enumerate(rowX):
      for indexZ,entryZ in enumerate(rowY):
        newEntry = value,entryZ[1],entryZ[2],entryZ[3]
        grid[indexX][indexY][indexZ] = newEntry
        
#change -2-travel dist to travel dist
def finalizeGridTravelDist(grid, gridSize):
  maxTD = 0.0
  for indexX,rowX in enumerate(grid):   
    for indexY,rowY in enumerate(rowX):
      for indexZ,entryZ in enumerate(rowY):
        if entryZ[0] == -1:
          newEntry = 0,entryZ[1],entryZ[2],entryZ[3]
          grid[indexX][indexY][indexZ] = newEntry
        elif entryZ[0] < -2:
          newEntry = ((-2-entryZ[0])*gridSize),entryZ[1],entryZ[2],entryZ[3]
          grid[indexX][indexY][indexZ] = newEntry
          maxTD = max(maxTD, newEntry[0])
        elif entryZ[0] > 0:
          newEntry = entryZ[0]*gridSize,entryZ[1],entryZ[2],entryZ[3]
          grid[indexX][indexY][indexZ] = newEntry
          maxTD = max(maxTD, newEntry[0])
  #grid modified in place, return maximum travel distance, useful sometimes
  return maxTD

def findPointsInCube(index, mins, maxs, gridSize, allPoints, points):
  '''returns a vector of all the indices of points in cube'''
  returnVec = []
  for pointIndex in allPoints:
    xyz = points[pointIndex-1][1:]
    cube = getIndices(mins, gridSize, xyz)
    if cube[0] == index[0] and \
       cube[1] == index[1] and \
       cube[2] == index[2]:
      returnVec.append(pointIndex)
  return returnVec

def findLongSurfEdges(pointList, pointNeighborList, gridSize, mins, maxs):
  '''returns the extra edges, i.e. a dictionary of all surface edges between  
  grid cubes with their euclidean distance between grid cube centers'''
  extraEdges = {} #empty dictionary
  surfaceEdgeBoxes = {}
  for pointNeighbors in pointNeighborList:
    pointStart = pointList[pointNeighbors[0]-1]
    startIndex = getIndices(mins, gridSize, pointStart[1:])
    surfaceEdgeBoxes[startIndex] = True #set them all to true, has_key is check
    endList = []
    for neighbors in pointNeighbors[2:]: #pN[1] is # of neighbors
      pointEnd = pointList[neighbors-1]
      endIndex = getIndices(mins, gridSize, pointEnd[1:])  
      gridLength = calcEdgeGridDist(pointStart[1:], pointEnd[1:], \
                                    mins, maxs, gridSize, metric='LINF')
      realLength = calcEdgeGridDist(pointStart[1:], pointEnd[1:], \
                                    mins, maxs, gridSize, metric='L2')
      if gridLength > 0: #no reason to add the ones within a grid cube
        if (endIndex, realLength, gridLength) not in endList: #no duplicates
          endList.append((endIndex, realLength, gridLength))
          #also want to add boxes connecting to surfaceEdgeBoxes dict??
    if len(endList) > 0:
      if startIndex not in extraEdges:
        extraEdges[startIndex] = endList
      else: #append
        newList = extraEdges[startIndex] + endList
        extraEdges[startIndex] = newList
      #print startIndex, endList
  return extraEdges, surfaceEdgeBoxes

def assignPointsValues(pointList, gridD, gridSize, mins, maxs, allPoints=False):
  '''helper function, finds which grid cube each surface point is in, determines
  depth value to assign
  record each points value based on which grid box it is in
  encoded value is -(val)-3 if < -2, 0 and -1 both map to 0, pos map to 1+num'''
  pointTravelDist = []
  for point in pointList:
    pointXYZ = point[1:]
    #compute x,y,z indices into grid     (use mins, maxs, gridSize)
    xIndex, yIndex, zIndex = getIndices(mins, gridSize, pointXYZ)
    if allPoints and point[0] not in allPoints:
      #print point[0], " not added, set to -1" #debugging fix loop
      pointTravelDist.append([point[0], -1])
    elif gridD[xIndex][yIndex][zIndex][0] > 0:
      pointTravelDist.append([point[0], \
                             (gridD[xIndex][yIndex][zIndex][0])*gridSize])
    elif gridD[xIndex][yIndex][zIndex][0] >= -1:
      pointTravelDist.append([point[0], 0])
    elif gridD[xIndex][yIndex][zIndex][0] == -2:  #problem
      #print point[0], " not added, set to -2" #debugging fix loop
      pointTravelDist.append([point[0], -2])
    else:
      pointTravelDist.append([point[0], \
                             (-2-gridD[xIndex][yIndex][zIndex][0])*gridSize])
  return pointTravelDist
