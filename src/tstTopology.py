#ryan g. coleman ryangc@mail.med.upenn.edu crystal.med.upenn.edu
#Copyright 2006, Ryan G Coleman, Kim A Sharp

import tstdata #data structure
import sys #for exiting
from priodict import priorityDictionary
import geometry
import tstdebug # for debugging
import mesh
import math

def fillInHolesAndGrow(triPointsOrig, pointTriListOrig, pointXyzOrig, \
                       normXyzOrig, numHandles, outFileName, debugOut, \
                       meshData, centerName):
  #copy all these since they will be modified
  loopPtsSave,loopTrisSave = [],[]
  triPoints = triPointsOrig[:] #so we don't mess up original
  pointTriList = pointTriListOrig[:] #so we don't mess up original
  pointXyz = pointXyzOrig #actually not modified, points never added or deleted
  normXyz = normXyzOrig[:] #might change these if i get around to it
  newTris = set() #list of the new tri indices added to triPoints
  returnLoops,returnPoints = [],[]
  trisRemoved = set()
  ringNeighborsAll = {}
  ringToDiscAll = {}
  ringToExcludeAll = {}
  ringToCenterNameAll = {}
  possibleHoleStartsAll = []
  ringToPotentialMaximaAll = {}
  for countHandleRemoved in range(numHandles):
    #print countHandleRemoved, "removed handles of", numHandles
    #run again... prepare for next loop
    #print countHandles(triPoints, pointTriList) #only necessary when debugging
    regLoopOkay,altOffset = False,0
    offsetOrder = __makeOffsetOrder(len(triPoints))
    while not regLoopOkay:
      #print altOffset, #necessary for debugging only
      regLoops, regPts, smallestList, regLoopOkay, loopTris, loopPts \
                     = __getListRegularLoop(triPoints, pointTriList, pointXyz, \
                  offsetOrder[altOffset], debugOut, countHandleRemoved, \
                  outFileName, normXyz, trisRemoved, newTris)
      #print regLoopOkay, smallestList
      if 0 == countHandleRemoved:
        loopTrisSave = loopTris
        loopPtsSave = loopPts
      if len(smallestList) > 0:
        regLoopOkay, smallest = False, False
        smallestIndex = 0
        while not regLoopOkay and len(smallestList) > smallestIndex:
          smallest = smallestList[smallestIndex]
          #print len(regPts[smallest][0]) + len(regPts[smallest][1])
          regLoopOkay = __checkRegularLoop(regLoops, regPts, smallest, \
                                           trisRemoved, newTris)
          #print regLoopOkay, smallest
          if regLoopOkay: #do more checks
            #additional check of whether plug maxima is inside convex hull here
            ringNeighbors, ringToDisc, ringToExclude, ringToCenterName, \
                 possHoleStarts = meshData.findOnePlug(countHandleRemoved, \
                       regPts[smallest], debugOut, outFileName, centerName)
            ringToPotentialMaxima = {}
            for ringCount,ring in enumerate(ringToDisc.keys()):
              maximaDisc, maximaDist = meshData.findOnePlugMaxima(\
                    countHandleRemoved*2+ringCount, ring, ringToDisc, \
                    ringToExclude, debugOut, outFileName, centerName)
              ringToPotentialMaxima[ring] = maximaDisc
              for potentialMaxima in maximaDisc:
                if meshData.nodeClass[potentialMaxima] == 0:
                  #print "rejecting based on plug outside convex hull", countHandleRemoved
                  regLoopOkay = False
          if not regLoopOkay: #try another
            smallestIndex += 1
        if regLoopOkay: #passed all tests, so we are keeping this reg loop
          ringNeighborsAll.update(ringNeighbors)
          ringToDiscAll.update(ringToDisc)
          ringToExcludeAll.update(ringToExclude)
          ringToCenterNameAll.update(ringToCenterName)
          #possibleHoleStartsAll.extend(possHoleStarts) #not done like this anymore
          ringToPotentialMaximaAll.update(ringToPotentialMaxima)
      altOffset += 1
      if altOffset > 30: #hacky kludge to get results
        #print "couldn't find any more appropriate loops"
        break
    #found one to change... add it to list to be returned in ascending order of size
    #print regKinds[smallest]
    for index in xrange(len(returnLoops)+1):
      if index == len(returnLoops):
        returnLoops.append(regLoops[smallest])
        returnPoints.append(regPts[smallest]) #stops anyway after this iteration
      elif len(returnLoops[index]) >= len(regLoops[smallest]):
        returnLoops.insert(index, regLoops[smallest])
        returnPoints.insert(index, regPts[smallest])
        break #stop the for loop, we added it
    #sometimes a triangle is made of 3 points adjacent in the loop... this is
    #a big problem later (makes edges with more than 2 triangles)
    #so we detect here and remove them with the loop to maintain topological correctness
    trisToRemove = regLoops[smallest][:] #copy
    loopsToRemove = [regPts[smallest][0][:],regPts[smallest][1][:]] #copy both
    #print len(triPoints), len(pointTriList), len(trisRemoved), len(newTris)
    __removeLoop(loopsToRemove, trisToRemove, pointTriList, triPoints, \
       trisRemoved, newTris, debugOut, outFileName, pointXyz,countHandleRemoved)
    #print len(triPoints), len(pointTriList), len(trisRemoved), len(newTris)
    #now triPoints and pointTriList have been modified and are consistent
  returnLoops.reverse() #actually want reverse order
  returnPoints.reverse()
  #now we do the bits that used to be in centergrow, i.e. grow tree and return
  pointNeighbors, outsidePoints, ringToMaximaAll = {}, [], {}
  for ringCount, ring in enumerate(ringToDiscAll.keys()):
    biggestMaxima = meshData.growOneMaxima(ringToMaximaAll, ring, \
         ringToPotentialMaximaAll[ring], centerName, ringToExcludeAll, \
         pointNeighbors, outsidePoints, ringNeighborsAll, ringToDiscAll)
    if biggestMaxima:
      possibleHoleStartsAll.append(biggestMaxima)
  possHoleStarts = meshData.makeNodes(possibleHoleStartsAll)
  pointNeighborsOut = meshData.makeNodesDict(pointNeighbors)
  outsidePointsOut = meshData.makeNodes(outsidePoints)
  return loopTrisSave, loopPtsSave, returnLoops, returnPoints, \
             pointNeighbors, pointNeighborsOut, \
             outsidePoints, outsidePointsOut, possHoleStarts

def __makeOffsetOrder(length):
  originalList = [x for x in xrange(length)]
  orderList, counter  = [], 0
  while len(originalList) > 0:
    if counter >= len(originalList):
      counter = counter % len(originalList)
    orderList.append(originalList.pop(counter))
    counter += 100
  return orderList

def __getListRegularLoop(triPoints, pointTriList, pointXyz, altOffset, \
      debugOut, countHandleRemoved, outFileName, normXyz, trisRemoved, newTris):
  '''gets one regular loop, possibly bad'''
  regLoopOkay = False
  loops = getLoops(triPoints, pointTriList, pointXyz, altOffset)
  if debugOut:
    tstdebug.debugTriangleList(loops, triPoints, pointXyz, outFileName + ".orig.loops."+str(countHandleRemoved)+".py")
  #print len(loops), (numHandles - countHandleRemoved) * 2 #no longer necessary, checks implemented
  loopPoints,loopsNew,regLoops,regPts,regKinds = regularizeLoops( \
                       triPoints, pointTriList, pointXyz, normXyz, loops)
  if debugOut:
    tstdebug.debugTriangleList(loopsNew, triPoints, pointXyz, outFileName + ".loops."+str(countHandleRemoved)+".py")
    tstdebug.debugTriangleList(regLoops, triPoints, pointXyz, outFileName+".regular.loops."+str(countHandleRemoved)+".py")
  loopTrisSave = loopsNew
  loopPtsSave = [] #might need reset
  for oneLoopPoints in loopPoints:
    loopPtsSave.append(oneLoopPoints[0]) #only need 1 of the 2
  #pick the smallest regular loop that goes around the hole and fill it in
  smallest = [] #now a list of possibilities sorted in ascending order of length
  for countRegLoop in xrange(len(regLoops)):
    regKind = regKinds[countRegLoop]
    if regKind[1] > regKind[0]: #means we are keeping it
      if len(smallest) == 0:
        smallest.append(countRegLoop)
      else:
        insertPoint = 0
        thisLen = len(regPts[countRegLoop][0]) + len(regPts[countRegLoop][1])
        while thisLen > len(regPts[smallest[insertPoint]][0]) + \
                           len(regPts[smallest[insertPoint]][1]):
          insertPoint += 1
          if insertPoint >= len(smallest):
            break
        smallest.insert(insertPoint, countRegLoop)
  if len(smallest) == 0: #still invalid...give up
    regLoopOkay = False
  return regLoops, regPts, smallest, regLoopOkay, loopTrisSave, loopPtsSave

def __checkRegularLoop(regLoops, regPts, oneToCheck, trisRemoved, newTris):
  #make sure doesn't contain invalid triangles
  regLoopOkay = True
  for tri in regLoops[oneToCheck]:
    if tri in trisRemoved or tri in newTris:
      #print "problem with loop", tri in trisRemoved, tri in newTris
      #print regLoops[oneToCheck]
      regLoopOkay = False #do loop again, otherwise done
  #this checks for loops in the boundary... redo if found
  if regLoopOkay: #no point in checking if already bad
    loopsToCheck = [regPts[oneToCheck][0][:],regPts[oneToCheck][1][:]] #copy both
    for loopToCheck in loopsToCheck:
      loopToCheck.sort()  #all duplicates adjacent
      last = -1 #impossible number
      for data in loopToCheck:
        if last == data:
          regLoopOkay = False # bad
        last = data
  return regLoopOkay

def __removeLoop(loopsToRemove, trisToRemove, pointTriList, triPoints, \
                 trisRemoved, newTris, debugOut, outFileName, pointXyz, \
                 countHandleRemoved):
  '''remove a loop from the structure, all data in first 4 parameters is
  modified in place, no returns'''
  if loopsToRemove:
    for checkLoop in loopsToRemove:
      needsChecked = True
      while len(checkLoop) > 2 and needsChecked:
        needsChecked = False #will be reset if checks need done again
        secondLastPt = checkLoop[-2]
        lastPt = checkLoop[-1]
        for thisPt in checkLoop:
          for triToCheck in pointTriList[thisPt-1][2:]:
            points = triPoints[triToCheck-1][1:]
            if lastPt in points and secondLastPt in points:
              #problem found... add this stuff to repair it
              needsChecked = True #redo this loop
              trisToRemove.append(triToCheck)
              #print lastPt, secondLastPt, checkLoop, points, thisPt
              checkLoop.remove(lastPt) #basically adding tri to loop and removing one point
              #print "added tri and removed pt", triToCheck, lastPt
          if needsChecked:
            break #start over with the checks
          secondLastPt = lastPt #update for next run
          lastPt = thisPt
    trisRemoved.update(trisToRemove)
    #now we have to fill in the hole.... tricky stuff.
    if debugOut:
      tstdebug.debugTriangleList([list(trisToRemove)], triPoints, pointXyz, outFileName+".removeTris."+str(countHandleRemoved)+".py")
    for delTri in trisToRemove: #get removed from triPoints
      for delTriPoint in triPoints[delTri-1][1:]:
        #these points can have the tri removed
        oldPTL = pointTriList[delTriPoint-1]
        oldEnd = oldPTL[2:]
        oldEnd.remove(delTri)
        newPTL = [oldPTL[0],oldPTL[1] - 1] + oldEnd
        pointTriList[delTriPoint-1] = newPTL
      triPoints[delTri-1] = [False] #can't fully delete because of indexing
    #now time to add new ones
    for loopIndex,vertexLoop in enumerate(loopsToRemove):
      trisToAdd = __trianglinizeLoopNoDupes(vertexLoop, \
                                                  triPoints, pointTriList)
      if debugOut:
        tstdebug.debugTrianglesNotOrig(trisToAdd, pointXyz, outFileName + ".addTris."+str(countHandleRemoved)+"."+str(loopIndex)+".py",ptList=list(vertexLoop))
      for addTriQuestionable in trisToAdd:
        #have to make this tri counterclockwise...
        addTri = makeClockwise(addTriQuestionable,triPoints,pointTriList)
        newTriIndex = len(triPoints) + 1
        newTris.add(newTriIndex)
        newTriRecord = [newTriIndex, addTri[0], addTri[1], addTri[2]]
        triPoints.append(newTriRecord)
        for addTriPoint in addTri:
          oldPTL = pointTriList[addTriPoint-1]
          oldEnd = oldPTL[2:]
          oldEnd.append(newTriIndex)
          newPTL = [oldPTL[0], oldPTL[1] + 1] + oldEnd
          pointTriList[addTriPoint-1] = newPTL
  #no returns as data has  been modified in place

def fillInHoles(triPointsOri, pointTriListOri, pointXyzOri, normXyzOri, \
                numHandles, outFileName, debugOut):
  """recursively find and fill in holes, returning regularized loops"""
  #copy all these since they will be modified
  loopPtsSave,loopTrisSave = [],[]
  triPoints = triPointsOri[:]
  pointTriList = pointTriListOri[:]
  pointXyz = pointXyzOri #actually not modified, points never added or deleted
  normXyz = normXyzOri[:] #might change these if i get around to it
  newTris = set() #list of the new tri indices added to triPoints
  returnLoops,returnPoints = [],[]
  trisRemoved = set()
  for countHandleRemoved in range(numHandles):
    #print countHandleRemoved, "removed handles of", numHandles
    #run again... prepare for next loop
    #print countHandles(triPoints, pointTriList) #only necessary when debugging
    regLoopOkay,altOffset = False,0
    while not regLoopOkay:
      #print altOffset, #necessary for debugging only
      regLoops, regPts, smallestList, regLoopOkay, loopTris, loopPts \
                     = __getListRegularLoop(triPoints, pointTriList, pointXyz, \
                  altOffset, debugOut, countHandleRemoved, outFileName, \
                  normXyz, trisRemoved, newTris)
      if 0 == countHandleRemoved:
        loopTrisSave = loopTris
        loopPtsSave = loopPts
      altOffset += 1
      smallestIndex, regLoopOkay = 0, False
      smallest = False
      while not regLoopOkay:
        #first round of checks
        smallest = smallestList[smallestIndex]
        regLoopOkay = __checkRegularLoop(regLoops, regPts, smallest, \
                                         trisRemoved, newTris)
        if not regLoopOkay:
          smallestIndex += 1
          if smallestIndex >= len(smallestList):
            break #out of while loop
    #found one to change... add it to list to be returned in ascending order of size
    #print regKinds[smallest]
    for index in range(len(returnLoops)+1):
      if index == len(returnLoops):
        returnLoops.append(regLoops[smallest])
        returnPoints.append(regPts[smallest]) #stops anyway after this iteration
      elif len(returnLoops[index]) >= len(regLoops[smallest]):
        returnLoops.insert(index, regLoops[smallest])
        returnPoints.insert(index, regPts[smallest])
        break #stop the for loop, we added it
    #sometimes a triangle is made of 3 points adjacent in the loop... this is
    #a big problem later (makes edges with more than 2 triangles)
    #so we detect here and remove them with the loop to maintain topological correctness
    trisToRemove = regLoops[smallest][:] #copy
    loopsToRemove = [regPts[smallest][0][:],regPts[smallest][1][:]] #copy both
    #print len(triPoints), len(pointTriList), len(trisRemoved), len(newTris)
    __removeLoop(loopsToRemove, trisToRemove, pointTriList, triPoints, \
       trisRemoved, newTris, debugOut, outFileName, pointXyz,countHandleRemoved)
    #print len(triPoints), len(pointTriList), len(trisRemoved), len(newTris)
    #now triPoints and pointTriList have been modified and are consistent
  returnLoops.reverse() #actually want reverse order
  returnPoints.reverse()
  return loopTrisSave, loopPtsSave, returnLoops, returnPoints

def __trianglinizeLoopNoDupes(vertexLoop, triPoints, pointTriList):
  allTris = [] #needed for findtri.. routine
  for triIndex in triPoints:
    if triIndex[0] != False:
      allTris.append(triIndex[0])
  modVertLoop = vertexLoop[:]
  everythingOkay = False
  while not everythingOkay:
    everythingOkay = True
    possibleTrisToAdd = tstdata.trianglinizeLoop(modVertLoop)
    for tri in possibleTrisToAdd:
      for edge in makeAllEdges(tri):
        otherTris = findAllTrisAcrossEdge(edge, allTris,\
                                             pointTriList,triPoints)
        if len(otherTris) > 1: #2 or greater is bad, means more than 2 tris per edge
          everythingOkay = False
    if not everythingOkay:
      temp = modVertLoop.pop(0)
      modVertLoop.append(temp)
      if modVertLoop == vertexLoop:
        return False #utter failure
  return possibleTrisToAdd


def fillInHolesTst(tstd, numHandles, outFileName, debugOut):
  """just a wrapper and easier way to call fillInHoles"""
  pointTriList = tstd.dict['POINT_TRIANGLE']
  triPoints = tstd.dict['TRIANGLE_POINT']
  pointXyz = tstd.dict['POINT_XYZ']
  normXyz = tstd.dict['NORM_XYZ']
  return fillInHoles(triPoints, pointTriList, pointXyz, normXyz, numHandles, \
                     outFileName, debugOut)

def expandDisc(triPoints, pointTriList, pointXyz, alternateStartOffset=0, \
               verboseDebug=False, verboseDebugFileName="temp.py"):
  '''entirely based on kim sharp's idea to expand a triangle, keeping boundary
  a simple disc, which leads to a collection of strips which are topologically
  important. if no handles, then expands until only a single triangle left'''
  #randomly pick the first triangle to start
  first = len(triPoints) - alternateStartOffset  #start at end, avoids new added tris better
  while triPoints[first-1] == [False]:
    first -= 1 #in case deleted during fill in hole procedure
    if first <= 0:
      print "serious error in finding appropriate loops"
      sys.exit(1)
  #print first,
  discTris, availableTris = set([first]), priorityDictionary()
  availTrisToEdge = {}
  clockwiseBoundaryPoints = triPoints[first-1][1:]
  unseenTris = set()
  for triPointRec in triPoints:
    if triPointRec != [False]:
      unseenTris.add(triPointRec[0])
  unseenTris.remove(first)
  if verboseDebug:
    uniqueStr = str(triPoints.count([False])) + "."
    numberTri = len(discTris)
    strNum = '%06d' % numberTri
    fractionSeen = float(numberTri)/float(len(triPoints))
    triColor = (1-fractionSeen,fractionSeen,0.)
    tstdebug.debugTriangles(discTris, triPoints, pointXyz, \
                strNum+"."+uniqueStr+verboseDebugFileName, \
                ptList=clockwiseBoundaryPoints, \
                triColor=triColor)
  for edge in makeAllEdges(triPoints[first-1][1:]):
    adjTri = findTriAcrossEdge(edge,unseenTris,pointTriList,triPoints)
    if -1 != adjTri:
      lengthToTri = lengthAcrossTris(first, adjTri, edge, triPoints, pointXyz)
      unseenTris.remove(adjTri)
      availableTris[adjTri] = lengthToTri #okay for now, later be careful
      availTrisToEdge[adjTri] = edge
  stripTris = set() #the ones that can't be added yet or at all
  #stopping criteria is that we've tried to add all possible tris on border
  curDistance = 0
  while len(availableTris) > 0:
    #print len(clockwiseBoundaryPoints), len(discTris), len(unseenTris), len(availableTris) #debug
    nextPossibleTri = availableTris.smallest()
    nextEdge = availTrisToEdge[nextPossibleTri]
    del availTrisToEdge[nextPossibleTri] #remove clutter
    curDistance = availableTris[nextPossibleTri]
    del availableTris[nextPossibleTri]
    possPoints = triPoints[nextPossibleTri-1][1:]
    if __checkBoundaryCounterEdge(clockwiseBoundaryPoints, nextEdge):
      nextEdge.reverse() #necessary to get edges pointed correctly
    otherEdges = makeOtherEdges(possPoints, nextEdge) #always need this
    #print nextEdge, otherEdges, clockwiseBoundaryPoints #debugging
    otherPoint = otherEdges[0][1] # or [1][0], they are the same
    try:
      found = clockwiseBoundaryPoints.index(otherPoint)
    except ValueError:
      found = -1
    safeToAdd = False
    numberEdgesInBoundary = 0
    if -1 == found: #other point not in boundary, safe to add, easy case
      safeToAdd = True
    if not safeToAdd: #still possibly okay, as long as exactly one edge is in
      for otherEdge in otherEdges:
        if __checkBoundaryEdge(clockwiseBoundaryPoints, otherEdge):
          numberEdgesInBoundary += 1
      if 1 == numberEdgesInBoundary:  #exactly 1
        safeToAdd = True
    #now it is either safe or not, proceed with updates and prepare for next iteration
    if not safeToAdd:
      #really we don't know if these are okay unless we try to add all others adjacent to it
      stripTris.add(nextPossibleTri) #possibly bad... don't know for sure
    else: #add it
      discTris.add(nextPossibleTri)
      if verboseDebug:
        numberTri = len(discTris)
        strNum = '%06d' % numberTri
        fractionSeen = float(numberTri)/float(len(triPoints))
        triColor = (1-fractionSeen,fractionSeen,0.)
        tstdebug.debugTriangles(discTris, triPoints, pointXyz, \
                    strNum+"."+uniqueStr+verboseDebugFileName, \
                    ptList=clockwiseBoundaryPoints, \
                    triColor=triColor)
      #unmark up to 2 striptris... as long as they aren't in disc already
      for otherEdge in otherEdges:
        unmarkTri = findTriAcrossEdge(otherEdge, \
                stripTris, pointTriList, triPoints)
        if -1 == unmarkTri:
          unmarkTri = findTriAcrossEdge(otherEdge, \
                  unseenTris, pointTriList, triPoints)
        if -1 == unmarkTri:
          unmarkTri = findTriAcrossEdge(otherEdge, \
                  availableTris.keys(), pointTriList, triPoints)
        if unmarkTri > -1: #else it is in the disc already
          combDistance = lengthAcrossTris(unmarkTri, nextPossibleTri, \
                                          otherEdge, triPoints, pointXyz)
          combDistance += curDistance
          if unmarkTri not in availableTris or \
                          combDistance < availableTris[unmarkTri]:
            availableTris[unmarkTri] = combDistance #can only go on here 3 times
            availTrisToEdge[unmarkTri] = otherEdge
          if unmarkTri in stripTris:
            stripTris.remove(unmarkTri) #only 1 of these will work
          elif unmarkTri in unseenTris:
            unseenTris.remove(unmarkTri)#only 1 of these will work
      #also update boundary points, 2 cases
      if 0 == numberEdgesInBoundary: #insertion case
        #add the other point in between the two
        try:
          insertion = clockwiseBoundaryPoints.index(nextEdge[0])
        except ValueError:
          print "serious error in tst data structure 2"
          print nextEdge, clockwiseBoundaryPoints
          tstdebug.debugTriangleList([list(availableTris.keys())], triPoints, pointXyz, "debug.availableTris.py")
          tstdebug.debugTriangleList([list(stripTris)], triPoints, pointXyz, "debug.stripTris.py")
          tstdebug.debugTriangleList([list(discTris)], triPoints, pointXyz, "debug.discTris.py",ptListList=[clockwiseBoundaryPoints])
          tstdebug.debugTriangleList([[nextPossibleTri]], triPoints, pointXyz, "debug.nextTri.py")
          sys.exit(1)
        clockwiseBoundaryPoints.insert(insertion+1, otherPoint)
      else: #deletion case, triangle takes out something
        try:
          deletion = clockwiseBoundaryPoints.index(otherPoint)
        except ValueError:
          print "serious error in tst data structure 3"
          sys.exit(1)
        #delete the one before or after the otherpoint that is in the tri
        tryDelete1 = (deletion - 1 + len(clockwiseBoundaryPoints)) % len(clockwiseBoundaryPoints)
        tryDelete1Value = clockwiseBoundaryPoints[tryDelete1]
        tryDelete2 = (deletion + 1 + len(clockwiseBoundaryPoints)) % len(clockwiseBoundaryPoints)
        tryDelete2Value = clockwiseBoundaryPoints[tryDelete2]
        try:
          possPoints.index(tryDelete1Value)
          delete1 = True
        except ValueError:
          delete1 = False
        try:
          possPoints.index(tryDelete2Value)
          delete2 = True
        except ValueError:
          delete2 = False
        if delete1:
          clockwiseBoundaryPoints.remove(tryDelete1Value)
        elif delete2:
          clockwiseBoundaryPoints.remove(tryDelete2Value)
        else:
          print "serious error in tst data structure 4"
          sys.exit(1)
  stripTris.update(unseenTris)
  return list(stripTris), clockwiseBoundaryPoints, list(discTris)

def __checkUniquePoints(pointListList):
  pointSet = set()
  for pointList in pointListList:
    for point in pointList:
      if point in pointSet:
        #print "point duplicated", point
        return False
      else:
        pointSet.add(point)
  return True  #made it here, should be okay

def expandLoop(loop, boundaries, triPoints, pointTriList, pointXyz, \
               verboseDebug=False, verboseDebugFileName="temp.py"):
  '''takes a loop, expands both sides, when they touch, return the tracebacks'''
  loopTris, availableTris = set(loop), \
                            [priorityDictionary(), priorityDictionary()]
  availableTriEdges = [{},{}] #keeps track of the edges crossed to get to tris
  if not __checkUniquePoints(boundaries):
    return False #don't try to do loops that have non-unique boundary points
  clockwiseBoundaryPointList = boundaries #length 2 list. checked clockwiseness
  unseenTris = set()
  for tri in range(len(triPoints)):
    unseenTris.add(tri+1) #index off by one
  tracebacks = {}
  if verboseDebug:
    uniqueStr = str(triPoints.count([False])) + "." + str(len(loop)) + "."
    numberTri = len(loopTris)
    strNum = '%06d' % numberTri
    fractionSeen = float(numberTri)/float(len(triPoints))
    triColor = (1-fractionSeen,fractionSeen,0.)
    tstdebug.debugTrianglesPtList(loopTris, triPoints, pointXyz, \
                strNum+"."+uniqueStr+verboseDebugFileName, \
                ptListList=clockwiseBoundaryPointList, \
                triColor=triColor)
  for loopTri in loop:
    tracebacks[loopTri] = [False] #false denotes end keyed from original loop
    unseenTris.remove(loopTri)
  for side,clockwiseBoundaryPts in enumerate(clockwiseBoundaryPointList):
    last = clockwiseBoundaryPts[-1]
    for current in clockwiseBoundaryPts:
      edge = [last,current]
      last = current
      adjTri = findTriAcrossEdge(edge,unseenTris,pointTriList,triPoints)
      if -1 != adjTri:
        unseenTris.remove(adjTri)
        other = (side + 1) % 2
        prevTri = findTriAcrossEdge(edge,loopTris,pointTriList,triPoints)
        if adjTri not in tracebacks:
          tracebacks[adjTri] = [] #initialize
        tracebacks[adjTri].append(prevTri)
        if  adjTri in availableTris[other]: #already have crossover
          #print "crossover1"
          return __crossoverReportLoop(adjTri, prevTri, tracebacks, loop, triPoints, pointTriList, pointXyz)
        combinedLength = lengthAcrossTris(prevTri,adjTri,edge,triPoints, \
                                          pointXyz)
        if adjTri not in availableTris[side] or \
                  availableTris[side][adjTri] > combinedLength:
          if adjTri in availableTris[side]:
            #print adjTri in availableTris[side], combinedLength
            del availableTris[side][adjTri]
          availableTris[side][adjTri] = combinedLength #would add current if not 0
          availableTriEdges[side][adjTri] = edge #otherwise something faster already seen
  side = 0 #side is set to which side is being processed, other is set to the other obviously
  stripTris = set() #the ones that can't be added yet or at all
  #all these need to be set before entering loop
  #side, stripTris, unseentris, clockwiseBoundaryPointList, loopTris, availableTris
  #stopping criteria is that we've tried to add all possible tris on border
  while len(availableTris[0]) > 0 or len(availableTris[1]) > 0:
    #print len(availableTris[0]), len(availableTris[1])
    side,smallest,other,otherValue = 0,1000000000.,1,1000000000.
    if len(availableTris[0]) > 0:
      smallest = availableTris[0][availableTris[0].smallest()]
    if len(availableTris[1]) > 0:
      otherValue = availableTris[1][availableTris[1].smallest()]
    if smallest > otherValue:
      side,smallest,other = 1,otherValue,0
    #now we've picked the side with the smallest one to expand on
    #print len(availableTris[0]), len(availableTris[1]), len(unseenTris),smallest
    nextPossibleTri = availableTris[side].smallest()
    nextEdge = availableTriEdges[side][nextPossibleTri]
    possPoints = triPoints[nextPossibleTri-1][1:]
    otherEdges = makeOtherEdges(possPoints, nextEdge) #always need this
    #print nextEdge, otherEdges, "checking edges"
    otherPoint = otherEdges[0][1] # or [1][0], they are the same
    if otherPoint in nextEdge:
      print nextEdge, otherEdges, "other edges problem"
      sys.exit(1)
    if otherPoint in clockwiseBoundaryPointList[other]:
      countOtherEdges = 0
      for otherEdge in otherEdges:
        if __checkBoundaryEdge(clockwiseBoundaryPointList[other], otherEdge):
          otherTri = findTriAcrossEdge(otherEdge, loopTris, pointTriList, triPoints)
          countOtherEdges += 1
      #print countOtherEdges
      if countOtherEdges > 0:
        prevTri = findTriAcrossEdge(nextEdge, loopTris, pointTriList, triPoints)
        #print "crossover2", countOtherEdges, prevTri, nextPossibleTri, otherTri
        #tstdebug.debugTriangleList([tracebacks.keys()], triPoints, pointXyz, "temp.built.up.to.this.py")
        return __crossoverReportLoop(nextPossibleTri, otherTri, tracebacks, loop, triPoints, pointTriList, pointXyz)
    try:
      found = clockwiseBoundaryPointList[side].index(otherPoint)
    except ValueError:
      found = -1
    safeToAdd = False
    if -1 == found: #other point not in boundary, safe to add, easy case
      safeToAdd = True
      numberEdgesInBoundary = 0
    if not safeToAdd: #still possibly okay, as long as exactly one edge is in
      numberEdgesInBoundary = 0
      for otherEdge in otherEdges:
        if __checkBoundaryEdge(clockwiseBoundaryPointList[side], otherEdge):
          numberEdgesInBoundary += 1
      if 1 == numberEdgesInBoundary:  #exactly 1
        safeToAdd = True
    #now it is either safe or not, proceed with updates and prepare for next iteration
    if not safeToAdd:
      #really we don't know if these are okay unless we try to add all others adjacent to it
      del availableTris[side][nextPossibleTri]
      del availableTriEdges[side][nextPossibleTri]
      #del tracebacks[nextPossibleTri] #remove since not added through here
      stripTris.add(nextPossibleTri) #possibly bad... don't know for sure
    else: #add it
      oldMetric = availableTris[side][nextPossibleTri]
      del availableTris[side][nextPossibleTri]
      del availableTriEdges[side][nextPossibleTri]
      loopTris.add(nextPossibleTri)
      #unmark 2 striptris
      for otherEdge in otherEdges:
        unmarkTri = findTriAcrossEdge(otherEdge, \
                  stripTris, pointTriList, triPoints)
        if -1 == unmarkTri:
          unmarkTri = findTriAcrossEdge(otherEdge, \
                  unseenTris, pointTriList, triPoints)
        if unmarkTri > -1:
          combinedLength = lengthAcrossTris(unmarkTri, nextPossibleTri, \
                          otherEdge, triPoints, pointXyz)
          if unmarkTri not in availableTris[side] or \
               availableTris[side][unmarkTri] > oldMetric + combinedLength:
            if unmarkTri in availableTris[side]:
              #print unmarkTri in availableTris[side], oldMetric + combinedLength
              del availableTris[side][unmarkTri]
            availableTris[side][unmarkTri] = oldMetric + combinedLength #can only go on here 3 times
            otherEdgeTemp = otherEdge[:]
            otherEdgeTemp.reverse()
            availableTriEdges[side][unmarkTri] = otherEdgeTemp #otherwise better way already seen
          if unmarkTri in availableTris[other]:
            #print "crossover3"
            return __crossoverReportLoop(unmarkTri, nextPossibleTri, tracebacks, loop, triPoints, pointTriList, pointXyz)
          if unmarkTri not in tracebacks:
            tracebacks[unmarkTri] = []
          tracebacks[unmarkTri].append(nextPossibleTri)
          if unmarkTri in stripTris:
            stripTris.remove(unmarkTri) #only 1 of these will work
          elif unmarkTri in unseenTris:
            unseenTris.remove(unmarkTri)#only 1 of these will work
          else:
            print "serious error in tst data structure 1 loop"
            tstdebug.debugTriangleList([list(availableTris[side].keys())], triPoints, pointXyz, "debug.availableTris.py")
            tstdebug.debugTriangleList([list(stripTris)], triPoints, pointXyz, "debug.stripTris.py")
            tstdebug.debugTriangleList([[unmarkTri]], triPoints, pointXyz, "debug.unmarkTri.py")
            sys.exit(1)
      #also update boundary points, 2 cases
      if 0 == numberEdgesInBoundary: #insertion case
        #add the other point in between the two
        try:
          insertion = clockwiseBoundaryPointList[side].index(nextEdge[0])
        except ValueError:
          print nextEdge, clockwiseBoundaryPointList[side]
          print triPoints[nextPossibleTri-1]
          print "serious error in tst data structure 2 loop"
          tstdebug.debugTriangleList([list(availableTris[side].keys())], triPoints, pointXyz, "debug.availableTris.py")
          tstdebug.debugTriangleList([list(stripTris)], triPoints, pointXyz, "debug.stripTris.py")
          tstdebug.debugTriangleList([[unmarkTri]], triPoints, pointXyz, "debug.unmarkTri.py")
          sys.exit(1)
        clockwiseBoundaryPointList[side].insert(insertion+1, otherPoint)
      else: #deletion case, triangle takes out something
        try:
          deletion = clockwiseBoundaryPointList[side].index(otherPoint)
        except ValueError:
          print "serious error in tst data structure 3 loop"
          sys.exit(1)
        #delete the one before or after the otherpoint that is in the tri
        tryDelete1 = (deletion - 1 + len(clockwiseBoundaryPointList[side])) % len(clockwiseBoundaryPointList[side])
        tryDelete1Value = clockwiseBoundaryPointList[side][tryDelete1]
        tryDelete2 = (deletion + 1 + len(clockwiseBoundaryPointList[side])) % len(clockwiseBoundaryPointList[side])
        tryDelete2Value = clockwiseBoundaryPointList[side][tryDelete2]
        try:
          possPoints.index(tryDelete1Value)
          delete1 = True
        except ValueError:
          delete1 = False
        try:
          possPoints.index(tryDelete2Value)
          delete2 = True
        except ValueError:
          delete2 = False
        if delete1:
          clockwiseBoundaryPointList[side].remove(tryDelete1Value)
        elif delete2:
          clockwiseBoundaryPointList[side].remove(tryDelete2Value)
        else:
          print "serious error in tst data structure 4 loop"
          sys.exit(1)
      if verboseDebug:
        uniqueStr = str(triPoints.count([False])) + "." + str(len(loop)) + "."
        numberTri = len(loopTris)
        strNum = '%06d' % numberTri
        fractionSeen = float(numberTri)/float(len(triPoints))
        triColor = (1-fractionSeen,fractionSeen,0.)
        tstdebug.debugTrianglesPtList(loopTris, triPoints, pointXyz, \
                   strNum+"."+uniqueStr+verboseDebugFileName, \
                   ptListList=clockwiseBoundaryPointList, \
                   triColor=triColor)
    #do one from the other side
    side = other
  #this is actually failure...
  return False

#easier way to call if you want
def expandDiscTst(tstd):
  return expandDisc(tstd.dict['TRIANGLE_POINT'], tstd.dict['POINT_TRIANGLE'], tstd.dict['POINT_XYZ'])

def getLoopsTst(tstd):
  pointTriList = tstd.dict['POINT_TRIANGLE']
  triPoints = tstd.dict['TRIANGLE_POINT']
  return getLoops(triPoints, pointTriList, tstd.dict['POINT_XYZ'])

def getLoops(triPoints, pointTriList, pointXyz,  \
                                      alternateStartOffset=0):
  #print "getLoops start", len(triPoints), len(pointTriList), alternateStartOffset
  stripTris, boundaryPts, discTris = expandDisc(triPoints, pointTriList, \
                                     pointXyz,  alternateStartOffset)
  #print len(stripTris), len(boundaryPts), len(discTris),
  seenTris = set() #stores the ones we've seen already
  curTri = stripTris.pop()
  seenTris.add(curTri)
  tree,loopConnections = makeTreeBranch(curTri, seenTris, stripTris, \
                                        triPoints, pointTriList)
  #print len(tree), len(loopConnections), curTri
  #now make loops of triangles...
  connectionTargets = set()
  for targets in loopConnections.itervalues():
    for target in targets:
      #print "target", target
      connectionTargets.add(target)
  listHalfLoops = buildListHalfLoops(tree,connectionTargets)
  #insert into loopConnections...
  for key,values in loopConnections.iteritems():
    newList = []
    for value in values:
      newList.append(listHalfLoops[value])
    loopConnections[key] = newList
  loops = makeLoops(tree,loopConnections)
  return loops

def getLoopsAndVerticesTst(tstd):
  pointTriList = tstd.dict['POINT_TRIANGLE']
  triPoints = tstd.dict['TRIANGLE_POINT']
  pointXyz = tstd.dict['POINT_XYZ']
  return getLoopsAndVertices(triPoints, pointTriList, pointXyz)

def getLoopsAndVertices(triPoints, pointTriList, pointXyz):
  loops = getLoops(triPoints, pointTriList, pointXyz)
  vertexLoops = []
  #tstdebug.debugTriangleList(loops, triPoints, pointXyz, "temp.original.loops.py")
  for index,loop in enumerate(loops):
    #tstdebug.debugTriangleList([loop], triPoints, pointXyz, "temp.original."+str(index)+".loop.py")
    vertexLoops.append(getVertexBoundary(loop, pointTriList, triPoints))
  return loops,vertexLoops

def regularizeLoopsTst(tstd, loops):
  '''run after getLoopsAndVertices, uses loops constructed, finds new ones
  that are hopefully regular or minimal around necks of handles'''
  pointTriList = tstd.dict['POINT_TRIANGLE']
  triPoints = tstd.dict['TRIANGLE_POINT']
  pointXyz = tstd.dict['POINT_XYZ']
  normXyz = tstd.dict['NORM_XYZ']
  return regularizeLoops(triPoints, pointTriList, pointXyz, normXyz, loops)

def regularizeLoops(triPoints, pointTriList, pointXyz, normXyz, loops):
  '''run after getLoopsAndVertices, uses loops constructed, finds new ones
  that are hopefully regular or minimal around necks of handles'''
  vertexLoopsSave,loopsSave = [],[]
  for index,loop in enumerate(loops):
    boundaries = getVertexBoundary(loop, pointTriList, triPoints, \
                                                  pointXyz, returnBoth=True)
    if boundaries == False: #special failure case
      pass
    else:
      vertexLoopsSave.append(boundaries)
      loopsSave.append(loop)
  loopsSave2 = []
  vertexLoopsSave2 = []
  regularLoops, regularLoopPts = [],[]
  for index in range(len(loopsSave)):
    regularLoop = expandLoop(loopsSave[index], vertexLoopsSave[index], \
                                          triPoints, pointTriList, pointXyz)
    if regularLoop: #more error checking involved here
      loopBoundaryPts = getVertexBoundary(regularLoop, pointTriList, \
                                 triPoints, pointXyz, returnBoth=True)
      if loopBoundaryPts == False: #special failure case
        pass
      else:
        regularLoopPts.append(loopBoundaryPts)
        regularLoops.append(regularLoop)
        loopsSave2.append(loopsSave[index])
        vertexLoopsSave2.append(vertexLoopsSave[index])
  kinds = [] #either "around handle" or "around hole"
  for index in range(len(loopsSave2)):
    kind = classifyRegularLoop(regularLoops[index], \
                            regularLoopPts[index], triPoints, pointXyz)
    kinds.append(kind)
  return vertexLoopsSave2, loopsSave2, regularLoops, regularLoopPts, kinds

def classifyRegularLoop(loopTris, loopPts, triPoints, pointXyz):
  '''now based on triangle normals pointing towards or away from center'''
  geometricCenter = averagePoints(loopPts[0] + loopPts[1], pointXyz)
  countSame,countDiff = 0,0
  for loopTri in loopTris:
    triPointXYZ = []
    for point in triPoints[loopTri-1][1:]:
      triPointXYZ.append(pointXyz[point-1][1:])
    normal = geometry.getTriNormal(triPointXYZ[0],triPointXYZ[1],triPointXYZ[2])
    centerTri = geometry.getAverage(triPointXYZ)
    vector = []
    for index in range(3):
      vector.append(centerTri[index] - geometricCenter[index])
    #make sure vector and normal are nonzero, both can happen, don't count it
    if (vector[0] == 0. and vector[1] == 0. and vector[2] == 0.) or \
         (normal[0] == 0. and normal[1] == 0. and normal[2] == 0.):
      pass #this condition happens when the center of the tri is near the
           #average of the points in the loop or when the tri is really small
    else:
      angle = geometry.getAngle(vector, normal)
      if angle < math.pi/2.:
        countSame += 1
      else:
        countDiff += 1
  #print countSame,countDiff
  return [countSame, countDiff]

def eulerCharacteristic(triPoints, pointTriList):
  vertices = 0
  faces = 0
  for triRec in triPoints:
    if triRec != [False]:
      faces += 1
  edges = 0
  for pointTriListRec in pointTriList:
    edges += pointTriListRec[1]
    if pointTriListRec[1] > 0:
      vertices += 1
  edges = edges / 2
  euler = vertices - edges + faces
  #print vertices, edges, faces, euler
  return euler

#assumes single piece
def countHandles(triPoints, pointTriList):
  euler = eulerCharacteristic(triPoints, pointTriList)
  handles = (euler-2)/-2
  return handles

#"private" helper functions follow

def averagePoints(ptList, pointXyz):
  geometricCenter,count = [0.,0.,0.],0
  for point in ptList:
    xyz = pointXyz[point-1][1:]
    count += 1
    for index in range(3):
      geometricCenter[index] = geometricCenter[index] + xyz[index]
  for index in range(3):
    geometricCenter[index] /= float(count)
  return geometricCenter

#takes a strip of triangles in order, returns one of the vertex boundaries
#returns clockwise loops
def getVertexBoundary(loop, pointTriList, triPoints, pointXyz=False, returnBoth=False):
  loopSet = set(loop)
  boundaryLists = [[],[]]
  for triangle in loop:
    edges = makeAllEdges(triPoints[triangle-1][1:])
    nonAdjacentEdge = False
    for edge in edges:
      adjTri = findTriAcrossEdge(edge,loopSet,pointTriList,triPoints,notThisTris=set([triangle]))
      if -1 == adjTri:
        if nonAdjacentEdge != False:
          print "serious error in loop structure 1"
          print loop, triangle
          tstdebug.debugTriangleList([loop], triPoints, pointXyz, "crash.loop.py")
          tstdebug.debugTriangleList([[triangle]], triPoints, pointXyz, "crash.tri.py")
          sys.exit(1)
        nonAdjacentEdge = edge
    if nonAdjacentEdge == False:
      #if pointXyz != False:
        #print triangle
        #tstdebug.debugTriangleList([[triangle]], triPoints, pointXyz, "crash.tri.py")
      #print "serious error in loop structure 7, no pointXyz needed to debug"
      return False #forces re-try
    added = False
    for boundaryList in boundaryLists:
      for point in nonAdjacentEdge:
        if len(boundaryList) > 0 and \
           boundaryList[-1] == point: #checking against last
          nonAdjacentEdge.remove(point)
          boundaryList.extend(nonAdjacentEdge)
          added = True
    if not added: #check list length, may need to reverse
      for boundaryList in boundaryLists:
        if 2 == len(boundaryList):
          for point in nonAdjacentEdge:
            if len(boundaryList) > 0 and boundaryList[0] == point:
              boundaryList.reverse()
              nonAdjacentEdge.remove(point)
              boundaryList.extend(nonAdjacentEdge)
              added = True
    if not added: #hopefully empty place to start new one
      for blIndex,boundaryList in enumerate(boundaryLists):
        if 0 == len(boundaryList) and not added:
          boundaryList.extend(nonAdjacentEdge)
          added = True
    if not added:
      print "serious error in loop structure 2"
      sys.exit(1)
  #trim off the repeated ends
  for boundaryList in boundaryLists:
    if len(boundaryList) > 0: #this only happens once, because the beg/end
      boundaryList.pop()      #are always the same by construction
  #special failure case... one border empty on certain occasions,
  #return false and let caller handle it
  if len(boundaryLists[0]) == 0 or len(boundaryLists[1]) == 0:
    return False
  #otherwise proceed...
  #okay now make sure both clockwise
  clockwiseSet = [False, False]
  for triangle in loop:
    edges = makeAllEdges(triPoints[triangle-1][1:])
    for edge in edges:
      for blIndex,boundaryList in enumerate(boundaryLists):
        if not clockwiseSet[blIndex]:
          if edge[0] in boundaryList:
            found = boundaryList.index(edge[0])
            next = (found + 1) % len(boundaryList)
            previous = (found - 1) % len(boundaryList)
            if boundaryList[next] == edge[1]:
              clockwiseSet[blIndex] = True
            elif boundaryList[previous] == edge[1]: #have to reverse
              boundaryList.reverse()
              clockwiseSet[blIndex] = True
  if not clockwiseSet[0] or not clockwiseSet[1]:
    print "serious error, can't make boundary lists clockwise"
    print clockwiseSet, boundaryLists, loop
    blah = [] #crash on purpose
    print blah[-1]
  if returnBoth: #return both, obviously
    return boundaryLists
  else:  #just return the smaller list
    if len(boundaryLists[0]) <= len(boundaryLists[1]):
      return boundaryLists[0]
    else:
      return boundaryLists[1]

#traverses tree(again), uses halfloops already make to construct new loops
#get rid of recursion, replace with stack
def makeLoops(tree, loopConnections):
  stackTreeTraversed = [(tree, [])] #starting stack
  returnList = []
  while len(stackTreeTraversed) > 0: #until stack is empty
    tree, traversed = stackTreeTraversed.pop(0)
    #the 1 == len(tree) case is executed no matter what
    traversed.insert(0, tree[0])
    if tree[0] in loopConnections:
      for target in loopConnections[tree[0]]:
        newLoop = target
        endNewLoop = traversed[:]
        #trim ends now, i.e. find the common ancestor (might be root)
        stillTrimming, index = True, 0
        while stillTrimming and index < len(newLoop): #walk down path
          node = newLoop[index]
          if endNewLoop[-1] == node: #if the first and last are the same
            endNewLoop.pop()         #we can eliminate from the loop
            index += 1
          else:
            stillTrimming = False
        finalLoop = newLoop[index-1:]
        finalLoop.extend(endNewLoop)
        returnList.append(finalLoop)
    if 2 == len(tree): #trunk (not a branch or leaf)
      stackTreeTraversed.insert(0, (tree[1],traversed[:]))
    elif 3 == len(tree): #branch
      for branch in tree[1:3]:
        stackTreeTraversed.insert(0, (branch,traversed[:]))
    #end of cases
  #idea is clean up 'double loops'
  #print len(returnList), "length of loop list"
  newReturnList = [] #created when there is a connector across the branches
  for loop in returnList: #down to the correct connector
    for source in loopConnections.keys():
      for targetList in loopConnections[source]: #once for each loop connector
        if len(targetList) > 0: #shouldn't happen but check anyway
          target = targetList[-1]
          newLoop = loop
          sourceIndex,targetIndex = False, False
          try:
            sourceIndex = newLoop.index(source)
            targetIndex = newLoop.index(target)
          except ValueError:
            pass #can just go back and try next connector...this one isn't here
          if sourceIndex and targetIndex:
            if sourceIndex == targetIndex + 1: #normal connector case...
              pass #nothing needs done
            else:    #we have a shortcut to make
              if sourceIndex > targetIndex:
                newLoop = loop[targetIndex:sourceIndex+1] #cut out extra stuff.
              else: #target > source
                newLoop = loop[:sourceIndex+1]
                newLoop.extend(loop[targetIndex:])
        loop = newLoop
    newReturnList.append(loop) #may have been modified
  returnList = newReturnList
  return returnList

#traverses tree, build list of half loops from root to each target, postorder
# traversal, this version doesn't use recursion (recursion depth problems)
def buildListHalfLoops(tree,connectionTargets):
  stackTreesTraversed = [(tree,[])] #starting data
  halfLoops = {} #dict keyed on members of connection targets
  while len(stackTreesTraversed) > 0:
    curTree,curTraversed = stackTreesTraversed.pop(0) #pop off the front
    #length is 1, 2, or 3... all cases now
    if 1 == len(curTree): #leaf
      if curTree[0] in connectionTargets:
        curTraversed.append(curTree[0])
        halfLoops[curTree[0]] = curTraversed[:] #dictionary keyed on connectionTarget
      #don't add anything to stack...
    elif 2 == len(curTree): #trunk (no branching)
      curTraversed.append(curTree[0])
      if curTree[0] in connectionTargets:
        halfLoops[curTree[0]] = curTraversed[:]
      #now add rest of trunk to stack
      stackTreesTraversed.insert(0,(curTree[1],curTraversed[:]))
    elif 3 == len(curTree): #branch
      curTraversed.append(curTree[0])
      if curTree[0] in connectionTargets:
        halfLoops[curTree[0]] = curTraversed[:]
      for branch in curTree[1:3]: #both branches
        stackTreesTraversed.insert(0,(branch,curTraversed[:]))
  return halfLoops #no more in stack, so return dictionary built

# tree making function  ---changed from recursion to iteration
def makeTreeBranch(firstTri, seenTris, unseenTris, triPoints, pointTriList, \
                   lastTri=False):
  loopConnections = {}
  #print "mtb", firstTri, len(seenTris), len(unseenTris), len(triPoints), len(pointTriList), len(loopConnections), lastTri,
  tree = ['replace']
  stackToProcess = [[firstTri, lastTri, seenTris, unseenTris, tree, 0]] #init
  while len(stackToProcess) > 0:
    firstTri, lastTri, seenTris, unseenTris, treeInsert, point = stackToProcess.pop()
    if 0 == len(unseenTris):
      treeInsert[point] = [firstTri] #that's it, no branching
    if False == lastTri: #initialize tree... first step
      pass
    #find adjacent triangles...
    edges = makeAllEdges(triPoints[firstTri-1][1:])
    adjacentTris = []
    for edge in edges:
      adjTri = findTriAcrossEdge(edge, unseenTris,pointTriList,triPoints)
      if -1 != adjTri:
        adjacentTris.append(adjTri)
      else:
        connectionTri = findTriAcrossEdge(edge, seenTris, pointTriList, \
                            triPoints, notThisTris=set([lastTri,firstTri]))
        if -1 != connectionTri: #record the connection in loopConnections now
          if connectionTri not in loopConnections or \
               firstTri not in loopConnections[connectionTri]:
            if firstTri not in loopConnections:
              loopConnections[firstTri] = [connectionTri]
            else:
              loopConnections[firstTri].append(connectionTri)
    #now decide what to do
    if 0 == len(adjacentTris):
      treeInsert[point] = [firstTri] #that's it, no branching
    elif 1 == len(adjacentTris): #going down a single branch.. normal
      nextTri = adjacentTris[0]
      unseenTris.remove(nextTri)
      seenTris.add(nextTri)
      treeInsert[point] = [firstTri, 'replace']
      stackToProcess.insert(0, [nextTri,firstTri,seenTris,unseenTris, \
                                treeInsert[point], 1]) #the 'replace'
    elif 2 == len(adjacentTris): #splitting, either at a 3-way or just starting
      for nextTri in adjacentTris:
        unseenTris.remove(nextTri)
        seenTris.add(nextTri) #both these must be done before recursive calls
      treeInsert[point] = [firstTri, 'replace1', 'replace2']
      for index,nextTri in enumerate(adjacentTris):
        stackToProcess.insert(0, [nextTri,firstTri,seenTris,unseenTris, \
                                  treeInsert[point], index+1]) #the 'replace'
    elif 3 == len(adjacentTris): #very bad... messy for data structure ,restart
      newStartTri = unseenTris.pop()
      unseenTris.extend(seenTris)
      seenTris.clear()
      return makeTreeBranch(newStartTri, seenTris, unseenTris,
                            triPoints, pointTriList)
      #just give up and start over... don't want degree 3 vertices
  return tree[0], loopConnections #dumb but necessary to return tree[0] not tree

def __checkBoundaryEdge(clockwiseBoundaryPoints, edge):
  '''returns true if either call is true'''
  if __checkBoundaryCounterEdge(clockwiseBoundaryPoints, edge):
    return True
  elif __checkBoundarySameEdge(clockwiseBoundaryPoints, edge):
    return True
  else:
    return False

def __checkBoundaryCounterEdge(clockwiseBoundaryPoints, counterEdge):
  #edge is counterclockwise, list is clockwise, return true or false
  try:
    found = clockwiseBoundaryPoints.index(counterEdge[0])
  except ValueError:
    return False #definitely not, one point not in list
  check = (found - 1 + len(clockwiseBoundaryPoints)) % len(clockwiseBoundaryPoints)
  if clockwiseBoundaryPoints[check] == counterEdge[1]:
    return True #the boundary contains the edge
  else:
    return False #or it does not

def __checkBoundarySameEdge(clockwiseBoundaryPoints, sameEdge):
  #edge is counterclockwise, list is clockwise, return true or false
  try:
    found = clockwiseBoundaryPoints.index(sameEdge[1])
  except ValueError:
    return False #definitely not, one point not in list
  check = (found - 1 + len(clockwiseBoundaryPoints)) % len(clockwiseBoundaryPoints)
  if clockwiseBoundaryPoints[check] == sameEdge[0]:
    return True #the boundary contains the edge
  else:
    return False #or it does not

#takes truple and tuple and returns other 2 possible tuples
def makeOtherEdges(possPoints, nextEdge):
  startIndex = possPoints.index(nextEdge[0])
  nextIndex = (startIndex + 1) % 3
  finalIndex = (startIndex + 2) % 3
  if nextEdge[1] != possPoints[finalIndex]: #edge needs reversed
    nextEdgeTemp = nextEdge[:]
    nextEdgeTemp.reverse()
    return makeOtherEdges(possPoints, nextEdgeTemp)
  return [[possPoints[finalIndex],possPoints[nextIndex]], \
          [possPoints[nextIndex], possPoints[startIndex]]]

#clockwise due to ordering of possPoints (from triPoints)
def makeAllEdges(possPoints):
  return [[possPoints[0],possPoints[1]], \
          [possPoints[1],possPoints[2]], \
          [possPoints[2], possPoints[0]]]

#validTris is a set
def findTriAcrossEdge(edge, validTris,pointTriList,triPoints, notThisTris=set()):
  nextPossibleTri = findAllTrisAcrossEdge(edge,validTris,\
                                             pointTriList,triPoints,notThisTris)
  if [] == nextPossibleTri:
    return -1
  else:
    #if len(nextPossibleTri) > 1:
    #  print nextPossibleTri, edge
    return nextPossibleTri[0]

def findAllTrisAcrossEdge(edge, validTris, pointTriList, triPoints,\
                                                   notThisTris=set()):
  nextPossibleTri = [] #default value, indicates failure
  for possibleTri in pointTriList[edge[0]-1][2:]:
    if possibleTri in validTris:
      if possibleTri not in notThisTris:
        possPoints = triPoints[possibleTri-1][1:]
        for possPoint in possPoints:
          if possPoint == edge[1]:
            nextPossibleTri.append(possibleTri)
  return nextPossibleTri

def __crossoverReportLoop(crossoverTri, otherTri, tracebacks, oldLoop, \
                          triPoints, pointTriList, pointXyz):
  '''makes loops for expandLoop procedure when crossover happens'''
  #print crossoverTri, otherTri,
  loop = [crossoverTri]
  currentTri = crossoverTri
  while False != currentTri:
    if tracebacks.has_key(currentTri):
      #if len(tracebacks[currentTri]) > 1:
      #  #print tracebacks[currentTri], currentTri
      currentTri = tracebacks[currentTri][-1] #last is best
      if currentTri != False:
        loop.append(currentTri)
    else:
      currentTri = False
  loop.insert(0, otherTri)
  currentTri = otherTri
  while False != currentTri:
    if tracebacks.has_key(currentTri):
      #if len(tracebacks[currentTri]) > 1:
      #  #print tracebacks[currentTri], currentTri
      currentTri = tracebacks[currentTri][-1] #last is best
      if currentTri != False:
        loop.insert(0,currentTri)
    else:
      currentTri = False
  #check to make sure ends match... remove one
  if loop[0] == loop[-1]:
    loop.pop()
    return loop
  else: #have to add some from the original loop
    start = oldLoop.index(loop[-1])
    end = oldLoop.index(loop[0])
    if abs(start-end) < abs(len(oldLoop)-(max(start,end))) + abs(min(start,end)):
      if end < start: #count down instead of up
        for index in range(end+1, start):
          loop.insert(0,oldLoop[index])
      else: #count up like normal
        for index in range(start+1,end):
          loop.append(oldLoop[index])
    else: #wrap-around
      if end < start: #count up from start then up to end
        for index in range(start+1, len(oldLoop)):
          loop.append(oldLoop[index])
        for index in range(end):
          loop.append(oldLoop[index])
      else: #count up from end then up to start
        for index in range(end+1, len(oldLoop)):
          loop.insert(0, oldLoop[index])
        for index in range(start):
          loop.insert(0, oldLoop[index])
    loop = __possiblyTrimLoop(loop, oldLoop, triPoints, pointTriList, pointXyz)
  #print oldLoop, loop
  return loop

def __getNextLoop(index, loopLength):
   return (index + 1) % loopLength

def __getPrevLoop(index, loopLength):
    return (index - 1 + loopLength) % loopLength

def  __possiblyTrimLoop(loop, oldLoop, triPoints, pointTriList, pointXyz):
  '''trims and returns loop, cutting out pieces of oldloop if they connect
  in 3-way overlaps'''
  returnLoop = [] #copy instead of break
  loopSet = set(loop)
  oldLoopSet = set(oldLoop)
  for index, triangle in enumerate(loop):
    prev = __getPrevLoop(index, len(loop))
    next = __getNextLoop(index, len(loop))
    edges = makeAllEdges(triPoints[triangle-1][1:])
    for edge in edges:
      adjTri = findTriAcrossEdge(edge,loopSet,pointTriList,triPoints,notThisTris=set([triangle]))
      if adjTri != loop[next] and adjTri != loop[prev]:
        if adjTri in oldLoopSet:
          adjIndex = loop.index(adjTri) #find where it is
          #print index, next, prev, adjIndex, loop, oldLoop #debugging
          #have to find which way to cut:
          while len(returnLoop) == 0:
            if loop[next] in oldLoopSet:
              if adjIndex > index:
                for place,tri in enumerate(loop):
                  if place <= index or adjIndex <= place:
                    returnLoop.append(tri)
              else: #adjIndex < index
                for place,tri in enumerate(loop):
                  if place >= adjIndex and index >= place:
                    returnLoop.append(tri)
            elif loop[prev] in oldLoopSet:
              if adjIndex > index:
                for place,tri in enumerate(loop):
                  if place >= index and adjIndex >= place:
                    returnLoop.append(tri)
              else: #adjIndex < index
                for place,tri in enumerate(loop):
                  if place <= adjIndex or index <= place:
                    returnLoop.append(tri)
            else: #advance next and prev pointers
              prev = __getPrevLoop(prev, len(loop))
              next = __getNextLoop(next, len(loop))
          if len(returnLoop) > 1:
            #returnloop should be okay... but possible we have 2 places wrong.. so
            #call this loop again...
            #print "problem found and fixed once"
            #print returnLoop
            return  __possiblyTrimLoop(returnLoop, oldLoop, triPoints, pointTriList, pointXyz)
  return loop #no problems found, return original

def lengthEdges(adjTri, triPoints, pointXyz, excludeEdges=[]):
  edges = makeAllEdges(triPoints[adjTri-1][1:])
  for edge in excludeEdges:
    try:
      edges.remove(edge)
    except ValueError:
      copyEdge = edge[:]
      copyEdge.reverse()
      edges.remove(copyEdge)
  #now get combined length of these two edges... that plus old is metric
  combinedLength = 0
  for leftEdge in edges:
    points = []
    for point in leftEdge:
      points.append(pointXyz[point-1][1:])
    combinedLength += geometry.distL2(points[0], points[1])
  return combinedLength

def lengthAcrossTris(triOne, triTwo, adjEdge, triPoints, pointXyz):
  """measures from one point to the other, across midpoint of adjacent edge"""
  edges = [makeAllEdges(triPoints[triOne-1][1:]), \
           makeAllEdges(triPoints[triTwo-1][1:])]
  farPoints = []
  for tris in edges:
    pts = [pointXyz[tris[0][0]-1][1:], \
           pointXyz[tris[1][0]-1][1:], \
           pointXyz[tris[2][0]-1][1:]   ]
    triCenter = (pts[0][0] + pts[1][0] + pts[2][0])/3., \
                (pts[0][1] + pts[1][1] + pts[2][1])/3., \
                (pts[0][2] + pts[1][2] + pts[2][2])/3.
    farPoints.append(triCenter) #should only be one left
  edgePoints = []
  for point in adjEdge:
    edgePoints.append(pointXyz[point-1][1:])
  edgeCenter = (edgePoints[0][0] + edgePoints[1][0])/2., \
               (edgePoints[0][1] + edgePoints[1][1])/2., \
               (edgePoints[0][2] + edgePoints[1][2])/2.
  return geometry.distL2(edgeCenter, farPoints[0]) + \
            geometry.distL2(edgeCenter, farPoints[1])

def makeClockwise(triQuestionable, triPoints, pointTriList):
  #print triQuestionable,
  validTris = set()
  for tri in triPoints:
    validTris.add(tri[0])
  edge = [triQuestionable[0], triQuestionable[1]]
  adjacentTri = findTriAcrossEdge(edge, validTris,pointTriList,triPoints)
  adjTriPts = triPoints[adjacentTri-1][1:]
  #print adjTriPts,
  #stupid cases now...
  if adjTriPts.index(edge[0]) + 1 == adjTriPts.index(edge[1]):
    #print "reversing"
    triQuestionable.reverse()
  elif adjTriPts.index(edge[0]) - 2 == adjTriPts.index(edge[1]):
    #print "reversing"
    triQuestionable.reverse()
  #print triQuestionable
  return triQuestionable
