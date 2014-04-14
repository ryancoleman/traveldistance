#ryan g. coleman ryangc@mail.med.upenn.edu
#this file is a module to do 2D orthogonal range searching
#by the simplest algorithm with the following properties:
# space requirements for tree O(n log n)
# time requirements for building tree O(n log n)
# query time  O(log^2 n + k) where k is the number of points in the range
#also allows any data to be stored hanging off the points in the data structure
#written in something like 12 hours march 23rd through 25th 2005
#3.2008 upgraded to python , including sets.Set() -> set() conversion

import operator # for itemgetter

#main class in this module
class orthoRangeSearchTree(object):
  """builds a 2D orthogonal range search tree and allows queries from it.

  Constructor takes points in the plane and builds a 2D range search tree
  in O(n log n) time and space, where n is the number of points passed in.

  Queries in the form of rectangles return all points in the rectangle in
  O(log^2 n + k) where k is the number of returned points.

  (note on complexity: uses's python tuned up quicksort, which in worst-case
  quadratic, but almost always (n log n))

  This algorithm uses the simplest form of the algorithm for this, with the
  correction allowing multiple points with the same x or y coordinates.

  Specially adapted to contain any amount of data tied to each point, just
  put it in the 'point' list passed in.

  Algorithms adapted from Chapter 5 of:

  Computational Geometry: Algorithms and Applications
  Second Edition
  Mark de Berg, Otfried Schwarzkopf TU Eindhoven (the Netherlands)
  Marc van Kreveld, Mark Overmars, Utrecht University (the Netherlands)

  2nd rev. ed. 2000. 367 pages, 370 fig.
  Hardcover DM 59
  ISBN: 3-540-65620-0
  """

  __orst = [] #orthogonal range search tree data structure saved here
  __inputData = [] #data about points, etc.

  ###following 3 methods are helpers for building binary search arrays
  #takes an integer, gives the right child in a binsearch array
  def __right(self, x):
    return int(2+int(2*x))

  #takes an integer, gives the left child in a binsearch array
  def __left(self, x):
    return int(1+int(2*x))

  #finds parent in a binsearch array
  #returns -1 if root
  def __par(self, x):
    if x % 2 == 0: #even
      return int(int(x/2)-1)
    else:
      return int(int((x+1)/2)-1)

  #finds sibling in a binsearch array
  #returns -1 if root
  def __sibling(self, x):
    if x % 2 == 0: #even
      return int(x-1)
    else:
      return int(x+1)

  #finds the median of a number where number is list length
  #returns 0 for x=0
  def __median(self, x):
    if x % 2 == 0: #even
      return int((x-2)/2)
    else:
      return int((x-1)/2)

  #private helper method used to make binary search trees in arrays
  #leaves inputList alone, copies what is necessary and returns it
  #inputList should contain sublists of length 2,where index 0 is the
  #sorted element, and index 1 is carried along (pointer to real data)
  def __makeBinTreeArray(self, inputList):
    #inputlist must be sorted on the 0th element of its sublist
    bsta = inputList[:] #Binary Search Tree Array
    #correct for non power of 2 in array by shifting stuff around
    order = 0
    while 2**order < len(bsta):
      order += 1
    if 2**order != len(bsta):
      diff = 2**order - len(bsta)
      bsta = inputList[-diff:] + inputList[:len(inputList)-diff]
    #corrected list built now
    for index in xrange(len(inputList)-1): #add empty records for internal nodes
      bsta.insert(0, [False])
    #now do the bottom-up build thing
    current = len(bsta)-1 #have to iterate backwards
    while current > 0: #don't need to do root
      if current % 2 == 1: #can skip following for right child
        #print "current is ", current #debugging
        rightMost = current
        rightVal = [bsta[rightMost][0]]
        while rightMost < len(bsta):
          rightVal = [bsta[rightMost][0]]
          #print rightMost, rightVal #debugging
          rightMost = self.__right(rightMost)
        bsta[self.__par(current)] = rightVal
      current -= 1 #decrement
    #print bsta #debug
    return bsta

  #tree format is [x, associated y-Tree, left, right]
  #2D tree is made top down (second level trees made bottom up)
  #top level tree is sorted by x
  def __makeTopLevel(self,inputData,xSorted,ySorted,initOrder):
    #print "debug recursive start " + str(len(xSorted))
    #initOrder = x, initOrder+1 = y
    #effectively make a 1-d range tree
    thisBSTA = self.__makeBinTreeArray(ySorted)
    if len(ySorted) == 1 and len(xSorted) == 1: #leaf node case
      thisEntry = [xSorted[0], thisBSTA, False, False]
      return thisEntry
    else: #split and recurse
      medianX = self.__median(len(xSorted))
      xSorLeft, xSorRight, ySorLeft, ySorRight = [],[],[],[] #new lists to make
      for index, xEntry in enumerate(xSorted):
        #find associated yEntry
        #print inputData[xEntry[1]], inputData[xEntry[1]][initOrder+1] #debug
        yEntry = inputData[xEntry[1]][initOrder+1]
        if index <= medianX:
          xSorLeft.append(xEntry)
          ySorLeft.append(yEntry)
        else:
          xSorRight.append(xEntry)
          ySorRight.append(yEntry)
      #print xSorLeft, xSorRight, ySorLeft, ySorRight        #debug
      #x arrays are sorted, y aren't
      ySorLeft.sort()
      ySorRight.sort()
      leftNode = self.__makeTopLevel(inputData,xSorLeft,ySorLeft,initOrder)
      rightNode = self.__makeTopLevel(inputData,xSorRight,ySorRight,initOrder)
      thisEntry = [xSorted[medianX], thisBSTA, leftNode, rightNode]
      return thisEntry
      #recursive calls complete data structure building

  #main init (constructor for the class)
  def __init__(self,inputData,xAxis=0,yAxis=1):
    """takes a planar set of points and constructs the range tree.

    a list of list is passed in as inputData with at least 2 elements in the
    sublists, a x and y coordinate.  Which is x and which is y can be specified
    by the xAxis and yAxis parameters.  The rest of the data in each sublist is
    kept at the point, and is all returned when queried.

    O(n log n) time and space, where n is the number of points passed in.

    (note on complexity: uses's python tuned up quicksort, which in worst-case
    quadratic, but almost always (n log n))

    To resolve a complete ordering, if points have same x and y coordinates,
    the passed-in order is used.

    Arguments:
    inputData -- no default, list is modified/destroyed/don't mess with it!
    xAxis -- default 0
    yAxis -- default 1
    """

    self.__orst = [] #clear this
    #first, must define complete ordering of points, sort by both x and y
    #record initial ordering
    self.__inputData = inputData #keep reference, assume we can fiddle with
                                 #passed in data and caller won't mess it up
    #save input ordering of points...
    for index,point in enumerate(inputData):
      point.append(index)
    initOrder = len(inputData[0])
    #sort on x, copy out x-coordinate and initial index
    inputData.sort(key=operator.itemgetter(xAxis)) #lambda x,y:cmp(x[xAxis],y[xAxis]))
    xSorted = [(xxx[xAxis],xxx[initOrder-1]) for xxx in inputData]
    for index,point in enumerate(inputData):
      point.append(xSorted[index])    #need reference back into this
    #sort on y, copy out y-coordinate and initial index
    inputData.sort(key=operator.itemgetter(yAxis)) #lambda x,y:cmp(x[yAxis],y[yAxis]))
    ySorted = [(xxx[yAxis],xxx[initOrder-1]) for xxx in inputData]
    for index,point in enumerate(inputData):
      point.append(ySorted[index])    #need reference back into this
    #resort back to original order
    inputData.sort(key=operator.itemgetter(initOrder-1)) #lambda x,y:cmp(x[initOrder-1],y[initOrder-1]))

    #2D tree is made top down (second level trees made bottom up)
    self.__orst = self.__makeTopLevel(inputData,xSorted,ySorted,initOrder)
    #data has been saved in __orst and __inputData for later
    #print self.__orst


  #private helper that finds the node where left and right split
  #or the leaf if they don't split, works for both 1d and 2d
  #returns the subtree of orst where split occurs
  def __findSplitNodeTopLevel(self,orst,leftLimit,rightLimit):
    curTree = orst
    curVal = curTree[0][0]
    while curTree[2] != False and curTree[3] != False and \
        (rightLimit <= curVal or leftLimit > curVal): #stop at leaf or diverge
      if rightLimit <= curVal:
        #print "go left", curTree
        curTree = curTree[2] #left subtree
      else:
        #print "go right", curTree
        curTree = curTree[3] #right subtree
      curVal = curTree[0][0]
    #print "split subtree is", curTree
    return curTree

  #private helper that finds the node where left and right split
  #or the leaf if they don't split, works for both 1d and 2d
  #returns the index into bsta where the split occurs
  def __findSplitNode(self,bsta,leftLimit,rightLimit):
    current = 0 #root
    curVal = bsta[current][0]
    while self.__left(current) < len(bsta) and \
        (rightLimit <= curVal or leftLimit > curVal):
      #stop once at leaf node, or at split
      #if in here, go whichever way down tree
      if rightLimit <= curVal:
        #print "go left", current, leftLimit, rightLimit, curVal
        current = self.__left(current)
      else:
        #print "go right", current, leftLimit, rightLimit, curVal
        current = self.__right(current)
      curVal = bsta[current][0]
    #print "split node is ", current, " for ", leftLimit, curVal,rightLimit, bsta
    return current

  #private helper does reporting of 1d subtrees
  def __oneDReportSubtree(self,y1,y2,bsta,appendResult,yTop):
    #print "1D subtree report "+ str(y1)+" "+str(y2)+" "+str(yTop)+" "+str(bsta)
    if yTop < len(bsta): #don't try to access missing nodes
      if len(bsta[yTop]) == 2: #at a leaf
        if y1 <= bsta[yTop][0] and bsta[yTop][0] <= y2:
          #print "adding", self.__inputData[bsta[yTop][1]]
          appendResult.append(self.__inputData[bsta[yTop][1]]) #done
      else: #not a leaf--recurse
        self.__oneDReportSubtree(y1,y2,bsta,appendResult,self.__left(yTop))
        self.__oneDReportSubtree(y1,y2,bsta,appendResult,self.__right(yTop))

  #private helper that does 1d-queries on passed in bsta
  #y1 and y2 are range, appendResult is array of results being built
  def __oneDRangeQuery(self,y1,y2,bsta,appendResult):
    #print "1D query " + str(y1) + " " + str(y2) + " " + str(bsta)
    botSplit = self.__findSplitNode(bsta,y1,y2)
    if len(bsta[botSplit]) == 2: #denotes a leaf
      #print "a leaf"
      if y1 <= bsta[botSplit][0] and bsta[botSplit][0] <= y2:
        #print "adding ", self.__inputData[bsta[botSplit][1]]
        appendResult.append(self.__inputData[bsta[botSplit][1]]) #done
    else:
      #print "not a leaf"
      #print "LEFT Y"      #LEFT
      current = self.__left(botSplit)
      while len(bsta[current]) != 2: #un-leaf
        if y1 <= bsta[current][0]:
          self.__oneDReportSubtree(
                 y1,y2,bsta,appendResult,self.__right(current))
          current = self.__left(current)
        else:
          current = self.__right(current)
      #clean up for current leaf
      if y1 <= bsta[current][0] and bsta[current][0] <= y2:
        #print "adding ", self.__inputData[bsta[current][1]]
        appendResult.append(self.__inputData[bsta[current][1]]) #done
      #print "RIGHT Y"      #RIGHT
      current = self.__right(botSplit)
      while len(bsta[current]) != 2: #un-leaf
        if y2 > bsta[current][0]:
          self.__oneDReportSubtree(
                 y1,y2,bsta,appendResult,self.__left(current))
          current = self.__right(current)
        else:
          current = self.__left(current)
      #clean up for current leaf
      if y1 <= bsta[current][0] and bsta[current][0] <= y2:
        #print "adding ", self.__inputData[bsta[current][1]]
        appendResult.append(self.__inputData[bsta[current][1]]) #done
      #DONE


  def rangeQuery(self,x1,x2,y1,y2):
    """returns all the points in the rectangle defined by those corners.

    Queries in the form of rectangles return all points in the rectangle in
    O(log^2 n + k) where k is the number of returned points.

    Range can be made inclusive by adding epsilon (small positive) to x2 and y2

    Arguments:
    x1 -- x-coordinate number 1 (left)
    x2 -- x-coordinate number 2 (right)
    y1 -- y-coordinate number 1 (down)
    y2 -- y-coordinate number 2 (up)
    """

    #reminder __orst contains tree, __inputData contains data
    #print str(self.__orst) + " end of data structure"
    topSplit = self.__findSplitNodeTopLevel(self.__orst,x1,x2) #tree returned
    posReturn = []
    if topSplit[2] != False and topSplit[3] != False: #leaf check
      #print "not a leaf",  topSplit
      #print "LEFT X"      #LEFT
      current = topSplit[2]
      while current[2] != False and current[3] != False: #un-leaves
        if x1 <= current[0][0]:
          self.__oneDRangeQuery(y1, y2, current[3][1],posReturn)
          current = current[2]
        else:
          current = current[3]
      #clean up leaf
      #print "clean up leaf", current
      if x1 <= current[0][0] and current[0][0] <= x2:
        self.__oneDRangeQuery(y1,y2,current[1],posReturn)
      #print "RIGHT Y"      #RIGHT
      current = topSplit[3]
      while current[2] != False and current[3] != False: #un-leaves
        if x2 > current[0][0]:
          self.__oneDRangeQuery(y1,y2, current[2][1],posReturn)
          current = current[3]
        else:
          current = current[2]
      #clean up leaf
      #print "clean up leaf", current
      if x1 <= current[0][0] and current[0][0] <= x2:
        self.__oneDRangeQuery(y1,y2,current[1],posReturn)
    else: #only one possible return, maybe none
      #print "a leaf--only chance"
      if x1 <= topSplit[0][0] and \
         topSplit[0][0] <= x2:
        self.__oneDRangeQuery(y1,y2,topSplit[1],posReturn)
    return posReturn

  def printOrst(self):
    """prints the internal structure of the orthogonal range search tree."""
    print self.__orst

  def printData(self):
    """prints the internal structure of the stored data."""
    print self.__inputData

  #helper for next function
  def __getNode(self, node):
    returnN = []
    returnN.append(node[0][0])
    if node[2] != False:
      returnN.append(self.__getNode(node[2]))
    else:
      returnN.append('$')
    if node[3] != False:
      returnN.append(self.__getNode(node[3]))
    else:
      returnN.append('$')
    return returnN


  def printOrstTopLevelTree(self):
    """prints the nodes in the top level (x-coord) of the tree structure."""
    print self.__getNode(self.__orst)

  #for testing/debugging purposes
  def naiveRangeQuery(self,x1,x2,y1,y2):
    """returns (hopefully) same as rangeQuery(.) but in O(n) time as check."""

    returnVec = []
    for xxx in self.__inputData:
      lXXX = len(xxx)
      xVal = xxx[lXXX-2][0]
      yVal = xxx[lXXX-1][0]
      if x1 <= xVal and xVal <= x2 and y1 <= yVal and yVal <= y2:
        returnVec.append(xxx)
    return returnVec

#following is code for examples/testing
#not real unit testing at all.. can be deleted
def generateRandomPoints(size=500):
  """generates random points, returned in vector format."""
  import random
  returnVec = []
  for index in xrange(size):
    x = round(random.uniform(-500,500),2)
    y = round(random.uniform(-500,500),2)
    returnVec.append([x,y])
  return returnVec

def generateRandomRange():
  """returns an xMin,xMax,yMin,yMax quadruple."""
  import random
  x1 = round(random.uniform(-1000,1000),3)
  x2 = round(random.uniform(-1000,1000),3)
  y1 = round(random.uniform(-1000,1000),3)
  y2 = round(random.uniform(-1000,1000),3)
  if x1 > x2:
    temp = x1
    x1 = x2
    x2 = temp
  if y1 > y2:
    temp = y1
    y1 = y2
    y2 = temp
  return (x1,x2,y1,y2)

def runTests(count=1, rangeCount=10, sizeTree=500):
  """runs count number of tests of building and rangeCount number of query checks.

  Generates random x,y pairs, builds trees, runs both kinds of queries and makes
  sure the results match.

  Reports differences.
  """

  print "prints a . for every check passed and a newline for every tree checked"
  for index in xrange(count):
    pairs = generateRandomPoints(sizeTree)
    orst = orthoRangeSearchTree(pairs)
    print "tree built"
    for queryIndex in xrange(rangeCount):
      testQuery = generateRandomRange()
      real = orst.rangeQuery(testQuery[0],testQuery[1],testQuery[2],testQuery[3])
      dumb = orst.naiveRangeQuery(\
                             testQuery[0],testQuery[1],testQuery[2],testQuery[3])
      realXY = [(xx[0],xx[1]) for xx in real]
      dumbXY = [(xx[0],xx[1]) for xx in dumb]
      realSet = set(realXY)
      dumbSet = set(dumbXY)
      if len(realSet.symmetric_difference(dumbSet)) != 0:
        print len(realXY), len(dumbXY)
        print realSet.symmetric_difference(dumbSet)
        print testQuery
        return orst, testQuery #break out for now, allow debugging
      else:
        print ".",
    print " "
  return 0,0
