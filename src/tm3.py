#!/usr/bin/env python2.5

#this file is for storing and writing 'tm3' files.
#http://www.cs.umd.edu/hcil/treemap/doc4.1/create_TM3_file.html
#for use with the treemap program(s)
#basically it is a table of values about things
#and a hierarchy of those things
#that writes/reads from a tab-delimited format
#also includes lots of methods to search and compare these structures

import string, sys, operator #python libraries
import statistics  #averages and such
import dot         #dot files for reading/writing to graphviz or aisee

#these types are defined here for later
floatType = type(1.)
intType = type(1)
strType = type("1")
#dateType not supported as of now

def conversion(stringType, data):
  '''a quick method to convert types'''
  if stringType == 'STRING':
    return str(data)
  elif stringType == 'FLOAT':
    return float(data)
  elif stringType == 'INTEGER':
    return int(data)

def getLineMst(connections):
  '''takes a set of connections that form a line, puts them in order, returns'''
  connsLimit2 = {}
  for connection in connections:
    node1fu,node2fu,score,extra1,extra2 = connection #unpack
    node1 = int(string.split(node1fu, "f")[0][1:])
    node2 = int(string.split(node2fu, "f")[0][1:])
    for node, otherNode in ((node1,node2),(node2,node1)):
      if node not in connsLimit2:
        connsLimit2[node] = []
      connsLimit2[node].append(otherNode)
  #find start
  lowestSingle = None   
  for head,tails in connsLimit2.iteritems():
    if len(tails) == 1:
      if lowestSingle is None:
        lowestSingle = head
      elif head < lowestSingle:
        lowestSingle = head
  outList = [lowestSingle]
  secondLast = lowestSingle
  lastSeen = connsLimit2[lowestSingle][0]
  #outList.append(lastSeen)
  done = False
  while not done:
    tails = connsLimit2[lastSeen]
    if len(tails) == 1:
      done = True #end of list so can quit
    else: #determine which of theses is the one you haven't already seen
      if tails[0] == secondLast:
        newSeen = tails[1]
      else:
        newSeen = tails[0]
    if not done:
      outList.append(lastSeen)
      secondLast = lastSeen #updates for next round
      lastSeen = newSeen
  return outList

def getStandardColumnsMeanStddev():
  '''returns a dict that has been pre-calculated based on the cast/surfnet
  data set, which seem to be pretty reasonable values for the mean and stddev'''
  columnToMean = {1: 2021.5149160870699, 2: 4611.729596728087, 
                  6: 5.1688493192518639, 7: 1.0553398650040864, 
                  8: 1.397754830189109, 9: 1787.604724497015, 
                  10: 18.882273257259961, 17: -7.1678825213817818, 
                  19: 18.9604617953798, 20: 15.147237413589277, 
                  21: 12.084412899466979} 
  colToStddev = {1: 3893.0666553174933, 2: 10372.377624807914, 
                 6: 8.1288279132363161, 7: 1.5324445735835186, 
                 8: 2.8725940273764614, 9: 3608.5463294631663, 
                 10: 26.678589504726215, 17: 8.5735418064033091, 
                 19: 26.833793052868906, 20: 22.42064979940854, 
                 21: 19.025192702415332}
  return columnToMean, colToStddev

def calcColumnsMeanStddev(columnList, tmDataList):
  '''returns a dict of column number to mean and another dict to stddev.'''
  columnsToMean = {}
  columnsToStddev = {}
  for column in columnList:
    colData = []
    for tmData in tmDataList:
      colData.extend(tmData.getListColumn(column))
    colAvg = statistics.computeAverage(colData)
    colStddev = statistics.computeStdDev(colData, colAvg)
    columnsToMean[column] = colAvg
    columnsToStddev[column] = colStddev
  return columnsToMean, columnsToStddev

def findSelfScore(tmData, columnList, resColNums, colToMean, colToStddev):
  '''does all against all comparison of nodes looking for similar ones
  looks only within 1 structure. returns  mean of 1/distance
  only consider non-overlapping other pockets '''
  matchMatrix = {} #dict of dicts of scores, keyed on nodes
  nodeList  = tmData.tree.keys()
  scores = []
  for tmNodeCount, node1 in enumerate(nodeList):
    matchMatrix[node1] = {} #init sub-dict
    for tmNode2Count, node2 in enumerate(nodeList):
      if tmNodeCount < tmNode2Count:
        resOverlap = dot.compareColumnsResidues(node1, node2, resColNums)
        if resOverlap == 100.: #if there is NO overlap within same protein
          score = dot.compareColumns(node1, node2, columnList, \
                                     colToMean, colToStddev)
          matchMatrix[node1][node2] = score
          scores.append(score)
  selfScore = {}
  for node in nodeList:
    lessScore, countScore = 0, 0
    for otherNode in nodeList:
      if node in matchMatrix and otherNode in matchMatrix[node]:
        if matchMatrix[node][otherNode] == 0.:
          lessScore += 1./0.00000000000000001
        else:
          lessScore += 1./matchMatrix[node][otherNode]
        countScore += 1
      elif otherNode in matchMatrix and node in matchMatrix[otherNode]:
        if matchMatrix[otherNode][node] == 0.:
          lessScore += 1./0.00000000000000001
        else:
          lessScore += 1./matchMatrix[otherNode][node]
        countScore += 1
    try:
      selfScore[node] = lessScore/float(countScore)
    except ZeroDivisionError:
      selfScore[node] = 0.
  return selfScore

def findSimilarNodes(tmDataList, columnListNames,  \
                     maxThr, maxConnectionCount, skipConnCount, \
                     resCols, sizeColName, sameStruct=False, sizeMin=-1., \
                     sizeMax=10000000000, mst=False, selfColListNames=None, \
                     doSelfScore=True, justNodes=None, lineMst=False, \
                     lineMstEnds=False, refinePockets=False, possNodes=None, \
                     dotResidues=None, printThings=True, calcColMeansStd=True, \
                     alpha=1.,clusterOutput=False):
  '''does all against all comparison of nodes looking for similar ones
  does lots of other stuff (need to update this documentation)'''
  justNodesValues = None
  if justNodes is not None:
    justNodesValues = justNodes.values()
  if selfColListNames is None:
    selfColListNames = columnListNames
  columnList = tmDataList[0].titlesToColumns(columnListNames)
  selfColList = tmDataList[0].titlesToColumns(selfColListNames)
  resColNums = tmDataList[0].titlesToColumns(resCols)
  sizeCol = tmDataList[0].titleToColumn(sizeColName)
  if calcColMeansStd:
    colToMean, colToStddev = calcColumnsMeanStddev(columnList, tmDataList)
  else:
    colToMean, colToStddev = getStandardColumnsMeanStddev()
  selfScores = {} #make combined list
  if doSelfScore:
    for tmData in tmDataList:
      selfSc = findSelfScore(tmData, selfColList, resColNums, colToMean, \
             colToStddev)
      selfScores.update(selfSc)
  #print selfScores.values()
  if len(tmDataList) == 1: #can't compare if only 1 struct
    #want to output nodes + self scores in sorted order for debugging of 
    #self scoring system
    treeName = tmDataList[0].inputFileName #only one
    selfScoreList = list(selfScores.iteritems())
    selfScoreList.sort(key=operator.itemgetter(1)) #sort by score
    if printThings:
      for node, selfScore in selfScoreList:
        if node.attributes[sizeCol] > sizeMin and \
                          node.attributes[sizeCol] < sizeMax:
          print dot.outputDrawStr(treeName, node) + "; " + str(selfScore)
    return None,None,None #no matches to return
  else: #more than 1 tree, do comparisons
    dotData = dot.dot(tmDataList)
    matchList = dotData.computeSearchConnections( \
                               maxThr, columnList, colToMean, colToStddev, \
                               resColNums, sizeCol, selfScores, sameStruct, \
                               sizeMin=sizeMin, sizeMax=sizeMax, \
                               doSelfScore=doSelfScore, \
                               justNodes=justNodesValues, alpha=alpha)
    if refinePockets:
      #iteration done here
      notDone = True
      eliminatedNodes = []
      while notDone:
        changeNode, examinedNode = dotData.refinePockets(columnList, colToMean,\
                               colToStddev, resColNums, sizeCol, selfScores, \
                               sizeMin=sizeMin, sizeMax=sizeMax, \
                               doSelfScore=doSelfScore, justNodes=justNodes, 
                               matchList=matchList, possNodes=possNodes, \
                               notTheseNodes=eliminatedNodes)
        if changeNode: #only recompute if there was a change
          if justNodes is not None:
            justNodesValues = justNodes.values()
          matchList = dotData.computeSearchConnections( \
                               maxThr, columnList, colToMean, colToStddev, \
                               resColNums, sizeCol, selfScores, sameStruct, \
                               sizeMin=sizeMin, sizeMax=sizeMax, \
                               doSelfScore=doSelfScore, \
                               justNodes=justNodesValues, alpha=alpha)
          eliminatedNodes = [] #reset since there was a change
        else: #worst wasn't changed this round
          eliminatedNodes.append(examinedNode)
          if examinedNode is None:
            notDone = False 
      if dotResidues is not None:
        dotResidues.setBestNodes(justNodes)
    connections = 0
    while connections <= maxConnectionCount:
      dotData.addSearchConnections(maxThr, maxConnCount=connections, mst=mst)
      #clusters = dotData.getClustersConnections()
      #print len(clusters)
      if not mst:
        if printThings:
          dotData.tempRemoveKeepersWriteGdl("search."+ \
                     string.zfill(connections,6) + ".gdl", \
                     justNodes=justNodesValues)
      elif mst: #change filenames, don't write edges
        if printThings:
          dotData.tempRemoveKeepersWriteGdl("mst."+ \
                   string.zfill(connections,4) + ".gdl", edges=False, \
                   justNodes=justNodesValues)
      if connections + skipConnCount > maxConnectionCount: #last time only
        if clusterOutput: #do but record the num of clusters as they are added
          dotData.addSearchConnections(maxThr, maxConnCount=connections*25, \
                 mst=mst, clusterOutput=clusterOutput)
        if lineMst: #do all again 
          dotData.addSearchConnections(maxThr, maxConnCount=connections, \
                                     lineMst=True)
          if printThings:
            dotData.tempRemoveKeepersWriteGdl("linemst."+ \
                   string.zfill(connections,4) + ".gdl", edges=False, \
                   justNodes=justNodesValues)
          mslList = getLineMst(dotData.connections)
          mslScore = statistics.listOrderCorrectness(mslList)
          if dotResidues is not None:
            dotResidues.mslScore = mslScore
          if printThings:
            print "linemst",
            for item in mslList:
              print item,
            print " " 
            print "linemstscore, ", mslScore
        if lineMstEnds: #do one more time
          #want to do same as linemst but give hints as to the endpoints 
          startNode = justNodes[tmDataList[0]]
          endNode = justNodes[tmDataList[-1]]
          dotData.addSearchConnections(maxThr, maxConnCount=connections, \
                            lineMst=True, startNode=startNode, endNode=endNode)
          if printThings:
            dotData.tempRemoveKeepersWriteGdl("linemstends."+ \
                   string.zfill(connections,4) + ".gdl", edges=False, \
                   justNodes=justNodesValues)
          mslList = getLineMst(dotData.connections)
          if printThings:
            print "linemstends", 
            for item in mslList:
              print item,
            print " " 
      connections += skipConnCount
    matchList.sort(lambda x,y: cmp(x[4],y[4])) #sort by score. best first
    if printThings:
      for match in matchList:
        print "%5.2f ;" % match[4],
        print dot.outputDrawStr(match[0].inputFileName, match[2]) + ";",
        print dot.outputDrawStr(match[1].inputFileName, match[3]) + ";"
      #print  " "
    return matchList

class tmNode(object):
  '''holds information about each node'''

  def __init__(self, attributes, tree=None): #assume user passes matching lists of att's
    '''puts the list of attributes in the new node'''
    self.attributes = attributes
    self.tree = tree #None otherwise a pointer to the tree this node is in

  def __repr__(self):
    '''prints out tab delimited with extra tab at end'''
    retStr = ""
    for attribute in self.attributes:
      retStr += str(attribute) + "\t"
    retStr += "\t"
    return retStr

  def typeRow(self):
    '''prints out INTEGER FLOAT or STRING for each attribute with tabs'''
    retStr = ""
    for attribute in self.attributes:
      if type(attribute) == floatType:
        retStr += "FLOAT\t"
      elif type(attribute) == intType:
        retStr += "INTEGER\t"
      elif type(attribute) == strType:
        retStr += "STRING\t"
      else:
        raise TypeError("Unknown type, must be integer, float or string")
    return retStr

  def getId(self):
    '''returns the first attribute, hopefully unique'''
    return str(self.attributes[0])

class tmTree(object):
  '''holds info about the hierarchy of nodes, the titles for each attribute,
  and is able to export + import .tm3 files'''

  def __init__(self, attributeTitles):
    '''just store the titles'''
    self.attributeTitles = attributeTitles
    self.tree = {} #dict stores children information for nodes
    self.idNode = {} #stores id->node map
    self.neighbors = {} #stores leaf neighbor info for tnv file
    self.parent = {} #built from tree when necessary

  def titleRow(self):
    '''export titles'''
    retStr = ""
    for attTitle in self.attributeTitles:
      retStr += attTitle + "\t"
    return retStr

  def addNode(self, attributes, childrenList=None):
    '''makes a node, if no children it is a leaf, otherwise connect it, return
    node for caller'''
    newNode = tmNode(attributes, self)
    self.idNode[newNode.getId()] = newNode
    #print self.idNode
    self.tree[newNode] = []
    if childrenList is not None:
      for childId in childrenList:
        self.tree[newNode].append(self.idNode[childId])
    return newNode

  def addChild(self, parentNode, childNode):
    '''using parent, adds the child to its tree. creates tree if necessary'''
    if parentNode not in self.tree:
      self.tree[parentNode] = []
    self.tree[parentNode].append(childNode) #no return necessary

  def setRoot(self, theNode):
    '''sets the root, important for output'''
    self.root = theNode

  def buildParent(self):
    '''uses self.tree to build self.parent'''
    for parent, children in self.tree.iteritems():
      for child in children:
        self.parent[child] = parent

  def isLeaf(self, aNodeId):
    '''true if this has no children'''
    children = self.tree[self.idNode[aNodeId]]
    if 0 == len(children):
      return True
    else:
      return False

  def isSingleChild(self, aNodeId):
    '''true if this has exactly 1 child'''
    if self.root == self.idNode[aNodeId]:
      return False #special case
    children = self.tree[self.idNode[aNodeId]]
    if 1 == len(children):
      return True
    else:
      return False

  def mergeSingleChild(self, aNodeId):
    '''checks to make sure node has a single child. deletes this node, making 
    the child and parent attach to each other'''
    if self.isSingleChild(aNodeId): #check otherwise don't do it
      children = self.tree[self.idNode[aNodeId]]
      del self.tree[self.idNode[aNodeId]]
      parent = self.parent[self.idNode[aNodeId]]
      del self.parent[self.idNode[aNodeId]]
      self.tree[parent].remove(self.idNode[aNodeId])
      self.tree[parent].extend(children)
      self.parent[children[0]] = parent
      del self.idNode[aNodeId] #delete it from this as well

  def deleteLeaf(self, aNodeId):
    '''removes the leaf with that id from the tree and the idnode dict'''
    children = self.tree[self.idNode[aNodeId]]
    if len(children) == 0:
      thisParent = self.parent[self.idNode[aNodeId]]
      del self.parent[self.idNode[aNodeId]] #now delete entry in parent
      oldChildren = self.tree[thisParent][:] #copy this
      oldChildren.remove(self.idNode[aNodeId]) #should always work
      self.tree[thisParent] = oldChildren
      del self.idNode[aNodeId] #delete it from this as well

  def bottomUpTraversal(self):
    '''returns a list of the nodes in bottom up order'''
    stackList = [self.root]
    outputList = []
    while len(stackList) > 0:
      current = stackList.pop()
      outputList.insert(0, current)
      children = self.tree[current]
      for child in children:
        stackList.append(child)
    return outputList

  def getListColumn(self, column):
    '''gets a list of all the values of all the nodes for one column'''
    retList = []
    for node in self.idNode.values():
      retList.append(node.attributes[column])
    return retList

  def connectLeafNodes(self, nodeA, nodeB):
    '''connects the leaf nodes for .tnv file generation'''
    if nodeA != nodeB:
      if nodeA not in self.neighbors:
        self.neighbors[nodeA] = []
      if nodeB not in self.neighbors[nodeA]:
        self.neighbors[nodeA].append(nodeB)
      if nodeB not in self.neighbors:
        self.neighbors[nodeB] = []
      if nodeA not in self.neighbors[nodeB]:
        self.neighbors[nodeB].append(nodeA)

  def getLeafToGroup(self):
    '''returns a dict of each leaf to all groups above it in the tree.
    actually returns a dict of each node to all nodes above it. since some 
    leaves are actually branch points.'''
    leafToGroup = {}
    stackList = [[self.root]]
    while len(stackList) > 0:
      current = stackList.pop()
      children = self.tree[current[-1]]  #first find the children
      for child in children:             #build the stack up
        newList = current[:] #copy
        newList.append(child)
        stackList.append(newList)
      leafToGroup[int(current[-1].getId())] = []      #dict constructed here
      for group in current:
        leafToGroup[int(current[-1].getId())].append(int(group.getId()))
    return leafToGroup

  def write(self, fileName):
    '''exports the tm3 nodes and tree as a tm3 file'''
    outFile = open(fileName, 'w')
    outFile.write(self.titleRow() + "\n")
    outFile.write(self.root.typeRow() + "\n") #header done
    stackList = [[self.root]]
    while len(stackList) > 0:
      current = stackList.pop()
      hierarchy = ""
      for node in current:
        hierarchy += node.getId() + "\t"
      children = self.tree[current[-1]]
      for child in children:
        newList = current[:] #copy
        newList.append(child)
        stackList.append(newList)
      outFile.write(str(current[-1]) + hierarchy + "\n")
    outFile.close()

  def writeTNV(self, fileName):
    '''exports the tm3 leaf nodes as a .tnv file'''
    outFile = open(fileName, 'w')
    outFile.write(self.attributeTitles[0] + "\n")
    for aNode in self.neighbors.keys():
      outFile.write(aNode + "\t")
      for neighborNode in self.neighbors[aNode]:
        outFile.write(neighborNode + "\t")
      outFile.write("\n")
    outFile.close()
  
  def titlesToColumns(self, titles):
    '''returns list of columns that the named titles are in'''
    listOut = []
    for title in titles:
      listOut.append(self.titleToColumn(title))
    return listOut

  def titleToColumn(self, title):
    '''returns the column that the named title is in'''
    return self.attributeTitles.index(title)

  def getAllNodes(self):
    '''returns a list of all the nodes in no order whatsoever'''
    return self.idNode.values()

  def compareResidueNodes(self, tmNode1, tmNode2, \
                          columnName="Residue Name List"):
    '''takes two nodes, finds the overlap between their lining residues'''
    resColumn = self.titleToColumn(columnName)
    resIds1 = set(string.split(tmNode1.attributes[resColumn],"+"))
    resIds2 = set(string.split(tmNode2.attributes[resColumn],"+"))
    try:
      resIds1.remove('') #how did these get here? stupid string methods
      resIds2.remove('') 
    except KeyError:
      pass #ignore if not there
    union = float(len(resIds1.union(resIds2)))
    if union > 0.: #otherwise automatically no match
      #print resIds1, resIds2, union , len(resIds1.intersection(resIds2))
      thisScore = (float(len(resIds1.intersection(resIds2)))/union) 
      return thisScore
    else:
      return 0. #no overlap (union), so score is automatically 0

  def convertResSetToIdentities(self, oldSet):
    '''turns a set of residues with numbers to a set with counts. e.g.
    ALA127, ALA128, GLY78 becomes ALAL0, ALA1, GLY0'''
    newSet = set()
    for item in oldSet:
      newItemPrefix = item[:3]
      count = 0
      while newItemPrefix + str(count) in newSet:
        count += 1 #increment until we're adding the unique residue number
      newSet.add(newItemPrefix + str(count)) #new unique one
    #print oldSet, newSet #debug the input/output
    return newSet

  def compareResidueIdentityMultipleNodes(self, tmNodeList, \
                          columnName="Residue Name List"):
    '''computes the overlap/tanimoto score over an entire list of nodes'''
    if 1 == len(tmNodeList):
      return 1. #max score since only 1 node
    resColumn = self.titleToColumn(columnName)
    resIds = []
    for tmNode in tmNodeList:
      resIds1 = set(string.split(tmNode.attributes[resColumn],"+"))
      try:
        resIds1.remove('') #how did these get here? stupid string methods
      except KeyError:
        pass #ignore if not there
      #now convert to ALA1, ALA2, etc
      newResIds1 = self.convertResSetToIdentities(resIds1)    
      resIds.append(newResIds1)
    unionSet = resIds[0].copy()
    for otherSet in resIds[1:]:
      unionSet.update(otherSet) #update is union_update renamed
    union = float(len(unionSet))
    if union > 0.: #otherwise automatically no match
      intersectSet = resIds[0].copy()
      for otherSet in resIds[1:]:
        intersectSet.intersection_update(otherSet)
      thisScore = (float(len(intersectSet))/union) 
      return thisScore
    else:
      return 0. #no overlap (union), so score is automatically 0

  def compareResidueIdentityNodes(self, tmNode1, tmNode2, \
                          columnName="Residue Name List"):
    '''takes two nodes, finds the overlap between their lining residues.
    only care about residue type not number'''
    resColumn = self.titleToColumn(columnName)
    resIds1 = set(string.split(tmNode1.attributes[resColumn],"+"))
    resIds2 = set(string.split(tmNode2.attributes[resColumn],"+"))
    try:
      resIds1.remove('') #how did these get here? stupid string methods
      resIds2.remove('') 
    except KeyError:
      pass #ignore if not there
    #now convert to ALA1, ALA2, etc
    newResIds1 = self.convertResSetToIdentities(resIds1)    
    newResIds2 = self.convertResSetToIdentities(resIds2)    
    union = float(len(newResIds1.union(newResIds2)))
    if union > 0.: #otherwise automatically no match
      thisScore = (float(len(newResIds1.intersection(newResIds2)))/union) 
      return thisScore
    else:
      return 0. #no overlap (union), so score is automatically 0

  def compareResidueList(self, residueList, tmNode, \
                         columnName="Residue Name List", \
                         lessThanRes=True, penalty=2):
    '''compares the column named to the residue list in +RES1+RES2 format'''
    resColumn = self.titleToColumn(columnName)
    resIds1 = set(string.split(residueList, "+"))
    resIds2 = set(string.split(tmNode.attributes[resColumn],"+"))
    resPocketCount = len(resIds2)
    try:
      resIds1.remove('') #how did these get here? stupid string methods
      resIds2.remove('') 
    except KeyError:
      pass #maybe not there after all
    union = float(len(resIds1.union(resIds2)))
    if union > 0.: #otherwise automatically no match
      #print resIds1, resIds2, union , len(resIds1.intersection(resIds2))
      if not lessThanRes:
        thisScore = 100-(float(len(resIds1.intersection(resIds2)))/union) * 100.
        return thisScore, resPocketCount 
      else: #means want to return a score based on having just list, not more
        lList = float(len(resIds1))
        pList = float(len(resIds2))
        intersect = float(len(resIds1.intersection(resIds2)))
        return 100 - ((intersect/lList) * 100.) + \
                         (((pList-intersect)/pList)*100.*penalty), resPocketCount
    else:
      return 100., resPocketCount #no intersection, worst score possible

  def findBestResidueListMatch(self, resList, colName="Residue Name List", \
                               lessThanRes=False, penalty=1.):
    '''go through all the tmnodes in this tree, find the best match'''
    bestMatch, bestScore, bestSize  = False, 100.*(penalty+10.), 100000.
    someOverlapNodes = []
    for aNode in self.getAllNodes():
      thisScore, thisSize = self.compareResidueList(resList, aNode, colName, \
                                          lessThanRes, penalty)
      if thisScore < bestScore or \
                (thisScore == bestScore and thisSize < bestSize):
        bestScore = thisScore
        bestMatch = aNode
        bestSize = thisSize
      if thisScore < 50. and thisSize > 0.:
        someOverlapNodes.append(aNode)
    if 0 == len(someOverlapNodes): #should have at least 1
      someOverlapNodes.append(bestMatch)
    return bestMatch, bestScore, someOverlapNodes

  def getMaxTravelDepth(self, maxTDcolName='maxTD'):
    '''looks at all nodes, finds the max maxTD, returns it'''
    colNumber = self.titleToColumn(maxTDcolName)
    maxMaxTD = 0.
    for aNode in self.getAllNodes():
      maxMaxTD = max(maxMaxTD, aNode.attributes[colNumber])
    return maxMaxTD

class tmTreeFromFile(tmTree):
  '''a subclass that can be read from a file instead of made from scratch'''

  def __init__(self, fileName):
    self.inputFileName = fileName
    treeFile = open(fileName, 'r')
    try:
      idTree = {} #temporary tree, used to make self.tree later
      attributeLine = treeFile.readline()            
      attributeTitles = string.split(string.strip(attributeLine), '\t')
      tmTree.__init__(self, attributeTitles)    
      types = string.split(treeFile.readline())
      for line in treeFile: #the rest of the lines
        tokens =  string.split(line)
        attributes = []
        for count,type in enumerate(types):
          attributes.append(conversion(type, tokens[count]))
        newNode = tmNode(attributes, self)
        thisId = newNode.getId()
        if thisId not in idTree:
          idTree[thisId] = [] #leaves have no children
        self.idNode[thisId] = newNode
        rootId = tokens[len(types)]
        parentId = tokens[-2]    
        if rootId == thisId:  #this is the root
          self.root = newNode
        else:
          if parentId not in idTree:
            idTree[parentId] = []
          idTree[parentId].insert(0, thisId) #this sorting keeps output the same
    except StopIteration:
      pass #eof signal
    finally:
      treeFile.close() #close the file
      #now make actual tree from idTree
      for parent, children in idTree.iteritems():
        childList = []
        for child in children:
          childList.append(self.idNode[child])
        self.tree[self.idNode[parent]] = childList
      self.buildParent() #also build the parent list, necessary for node dels
    
#main is only run for testing import of tm3 from commandline        
if -1 != string.find(sys.argv[0], "tm3.py"):
  for filename in sys.argv[1:]:
    tmData = tmTreeFromFile(filename)
    #tmData.write("temp.tm3") #can test output, should be identical
    #for node in tmData.bottomUpTraversal():
    #  if node in tmData.tree: #may have already been merged
    #    if tmData.isLeaf(node.getId()):
    #      print node.getId(),
    #print

