#Ryan G. Coleman, Kim A. Sharp, http://crystal.med.upenn.edu
#These two functions used to be the main travel distance loop,
#but have been replaced

from priodict import priorityDictionary  # david eppstein's heap dictionary
import grid  # moved primitive grid functions

def figureOutEdge(
    curBox, curDist, oldEncode, currentBoxes, history, adjBox, gapDist,
    surfEdge=False):
  '''this is the inner loop of travel distance, where an edge is processed.'''
  if (oldEncode >= 0 and oldEncode > curDist + gapDist) and \
     (adjBox not in currentBoxes or curDist + gapDist < currentBoxes[adjBox]):
    currentBoxes[adjBox] = curDist + gapDist
    prevBoxList = history.get((adjBox, curDist + gapDist), [])
    prevBoxList.append(curBox[0])
    history[(adjBox, curDist + gapDist)] = prevBoxList
  if not surfEdge:
    if (oldEncode < -1 and curBox[1] >= 0) and (
        oldEncode == -2 or oldEncode < -2-(gapDist)-curDist) and (
            adjBox not in currentBoxes or
            curDist + gapDist < currentBoxes[adjBox]):
      currentBoxes[adjBox] = curDist + gapDist
      prevBoxList = history.get((adjBox, curDist + gapDist), [])
      prevBoxList.append(curBox[0])
      history[(adjBox, curDist + gapDist)] = prevBoxList
  else:  # surfEdge is True
    if (oldEncode < -1) and (
        oldEncode == -2 or oldEncode < -2-(gapDist)-curDist) and (
            adjBox not in currentBoxes or
            curDist + gapDist < currentBoxes[adjBox]):
      currentBoxes[adjBox] = curDist + gapDist
      prevBoxList = history.get((adjBox, curDist + gapDist), [])
      prevBoxList.append(curBox[0])
      history[(adjBox, curDist + gapDist)] = prevBoxList
  #no return value, currentBoxes and history are modified

def calcTravelDist(
    gridD, pointList, gridSize, mins, maxs, allPoints,  extraEdges,
    surfaceEdgeBoxes,  filename="default.tst", maxDistanceCutoff=False,
    volumePoints=False):
  '''old method (v1.0) for doing travel distance.
  grid encoding is -1 = starting grid cubes
                    0 = valid to spread to
                   -2 = invalid to spread to'''
  #useful variables to set
  lenX, lenY, lenZ = len(gridD), len(gridD[0]), len(gridD[0][0])
  #get these in sets now, makes later steps easier
  outsideBoxes, insideBoxes = grid.fillBoxesSets(gridD)
  currentBoxes = priorityDictionary()  # heap as dictionary
  history = {}  # plain dictionary, keyed on (box,dist) pairs, goes to prevBox
  tracebacks = {}  # assembled during run, keyed on ending grid point
  while len(outsideBoxes) > 0:  # can do all these without checking
    currentBoxes[outsideBoxes.pop()] = 0  # 0 is starting grid distance
    #no history for these
  pointTravelDist = []  # this is set when loop is done
  currentTravelDist = 0.0
  while (not maxDistanceCutoff or currentTravelDist <= maxDistanceCutoff) \
      and len(currentBoxes) > 0:
    curBox = [currentBoxes.smallest(), 0]  # dummy second element
    curBox[1] = currentBoxes.pop(currentBoxes.smallest())
    prevBoxList = history.pop((curBox[0], curBox[1]), False)  # get history info
    if prevBoxList and len(prevBoxList) > 0:
      tracebacks[curBox[0]] = [prevBoxList, curBox[0], curBox[1]]
    oldCurEncode = gridD[curBox[0][0]][curBox[0][1]][curBox[0][2]][0]
    curInside = False
    if oldCurEncode < -1:  # encoding still keeps track of inside/outside
      curInside = True
      curBox[1] = -curBox[1]-2  # transform to what code expects
    currentTravelDist = curBox[1]
    #print len(currentBoxes), curBox   # debugging
    #don't update if a lower depth has already been processed
    if (curBox[1] == 0 or (curBox[1] > 0 and oldCurEncode > curBox[1])
       or oldCurEncode == -2 or (curBox[1] < 0 and oldCurEncode < curBox[1])) \
       and (oldCurEncode > -2 or curBox[0] in surfaceEdgeBoxes):
      if curBox[1] != 0:  # only do these if not inside
        #update the current box we are examining
        oldBox = gridD[curBox[0][0]][curBox[0][1]][curBox[0][2]]
        newBox = curBox[1], oldBox[1], oldBox[2], oldBox[3]
        #print "writing ", curBox
        gridD[curBox[0][0]][curBox[0][1]][curBox[0][2]] = newBox
      curDist = curBox[1]
      if curDist < 0:  # adjust
        curDist = -curDist - 2  # make positive, subtract 2
      #go to up to 6 adjacent boxes, only good returned
      for adjBox in grid.getAdjacentBoxes(curBox[0], lenX, lenY, lenZ):
        oldEncode = gridD[adjBox[0]][adjBox[1]][adjBox[2]][0]
        figureOutEdge(
            curBox, curDist, oldEncode, currentBoxes, history, adjBox, 1)
      #also do side boxes
      for adjBox in grid.getSideAdjacentBoxes(curBox[0], lenX, lenY, lenZ):
        oldEncode = gridD[adjBox[0]][adjBox[1]][adjBox[2]][0]
        figureOutEdge(
            curBox, curDist, oldEncode, currentBoxes, history, adjBox, 2.**.5)
      #also do corner boxes
      for adjBox in grid.getCornerAdjacentBoxes(curBox[0], lenX, lenY, lenZ):
        oldEncode = gridD[adjBox[0]][adjBox[1]][adjBox[2]][0]
        figureOutEdge(
            curBox, curDist, oldEncode, currentBoxes, history, adjBox, 3.**.5)
      if curBox[0] in extraEdges:
        for adjBox, adjDist, adjGridDist in extraEdges[curBox[0]]:
          #print curBox,adjBox,adjDist,adjGridDist, lenX, lenY, lenZ
          oldEncode = gridD[adjBox[0]][adjBox[1]][adjBox[2]][0]
          figureOutEdge(
              curBox, curDist, oldEncode, currentBoxes,
              history, adjBox, adjDist, True)
  #check completely done stuff, make sure all boxes have value
  #first case, a box outside the MS and inside the CH did not get a value
  # propagated to it
  pointTravelDist = grid.assignPointsValues(
      pointList, gridD, gridSize, mins, maxs, allPoints)
  if volumePoints:
    volumePointDepths = grid.assignPointsValues(
        volumePoints, gridD, gridSize, mins, maxs)
  if volumePoints:
    return pointTravelDist, tracebacks, volumePointDepths
  else:
    return pointTravelDist, tracebacks, False
