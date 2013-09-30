#!/usr/bin/env python2.5

import tm3, string, sys, statistics, operator
import munkreskuhn
import unionfind2
import geometry
import collections #for defaultdict

colorDict = {0:'white',      
   1:'blue',             11:'darkyellow',     21:'lightcyan',
   2:'red',              12:'darkmagenta',        22:'lilac',       
   3:'green',            13:'darkcyan',           23:'turquoise',
   4:'yellow',           14:'gold',               24:'aquamarine',
   5:'magenta',          15:'lightgrey',          25:'khaki',
   6:'cyan',             16:'lightblue',         26:'purple',
   7:'darkgrey',         17:'lightred',          27:'yellowgreen',
   8:'darkblue',         18:'lightgreen',      28:'pink',
   9:'darkred',      19:'lightyellow',      29:'orange',
   10:'darkgreen',      20:'lightmagenta',       30:'orchid',
                                          31:'black'}

borderStyleDict={0:'solid', 1:'dashed', 2:'dotted', 3:'double', 4:'triple'}

def mapToColorScale(value, maxVal=1., colorScale=[8,1,16,0,17,2,9]):
  '''the value should be between 0 and maxVal, map to the colorscale'''
  increments = float(len(colorScale))
  increment = maxVal/increments
  which = min(int(value/increment),len(colorScale)-1)
  #print value, colorScale[which], colorDict[colorScale[which]]
  return colorDict[colorScale[which]]

def makeTstStr(orig):
  '''truncates everything after tst'''
  return orig[:string.rfind(orig,".tree")]

def outputDrawStr(treeName, node):
  '''prints a string to copy/paste into pymol to display the group'''
  return  "tstOpenDraw " + makeTstStr(treeName) +  \
           ",group="+ node.getId()

class dot(object):
  '''this class takes some tm3 files, produces a dot file of each tree and
  the relationships between them. dot files are standard files used with 
  graphviz to produce nice looking output'''
  
  class subgraph(object):
    '''holds the data for one tree'''
    
    def __init__(self, label, inTree, namePrefix, root):
      '''initialize the data structure'''
      self.label = label
      self.tree = {} #where data is stored
      self.prefix = namePrefix
      self.nodeToName = {}
      self.nameToNode = {}
      self.rootName = self.nodeToNameFun(namePrefix, root)
      self.keepers = {}
      self.keepers[self.rootName] = True
      for parent,children in inTree.iteritems():
        self.nodeToName[parent] = self.nodeToNameFun(namePrefix, parent)
        self.nameToNode[self.nodeToName[parent]] = parent
        newChildren = []
        for child in children:
          nodeName = self.nodeToNameFun(namePrefix, child)
          self.nameToNode[nodeName] = child
          newChildren.append(nodeName)
          if nodeName not in self.keepers:
            self.keepers[nodeName] = False
        if len(newChildren) > 0:
          self.tree[namePrefix+"_"+parent.getId()] = newChildren
          if len(newChildren) > 1: #disabled but used to keep all branch pts
            self.keepers[self.nodeToName[parent]] = True #if multi child keep

    def resetKeepers(self):
      '''reset the keepers to initial state'''
      self.keepers = {}
      for parent,children in self.tree.iteritems():
        self.keepers[parent] = False
        for child in children:
          self.keepers[child] = False
      for parent,children in self.tree.iteritems():
        if len(children) > 1: #keep all branch pts
          self.keepers[parent] = True
      self.keepers[self.rootName] = True

    def nodeToNameFun(self, namePrefix, node):
      '''method that outputs the name/label for a given node'''
      return namePrefix + "_" + node.getId() 

    def __repr__(self):
      '''returns a string containing all relevant data for dot file'''
      retStr = "subgraph cluster" + self.prefix + " {\n"
      retStr += "  node [shape=box,height=0,margin=0,0];\n"
      for parent,children in self.tree.iteritems():
        for child in children:
          retStr += "  " + parent + " -- " + child + ";\n"
      retStr += "  label = \""+self.label+"\";\n}\n"
      return retStr
   
    def gdlStr(self, colornumber=0, edges=True, edgeColor="lightgrey", \
               justNodes=None, colorMap=None):
      '''returns a string containing all relevant data for gdl file.
      either each tree is colored uniquely. or nodes are colored by residue 
      overlap to connected nodes and the borders are colored by tree.'''
      if colorMap is None:
        color = str(colorDict[colornumber % 32])
        bordercolor = 'black'
        borderstyle = 'solid'
      else:
        bordercolor = str(colorDict[colornumber % 32])
        borderstyle = str(borderStyleDict[(colornumber/32) % 5])
      retStr = " graph: { title:\"" + self.prefix + "\"\n"
      #retStr += "layoutalgorithm : normal\n"
      retStr += "edge.arrowstyle: line\n"
      namesDrawn = {}
      for parent,children in self.tree.iteritems():
        parNode = self.nameToNode[parent]
        if parent not in namesDrawn:
          drawStr = outputDrawStr(self.label, self.nameToNode[parent])
          parId = parNode.getId()
          if parent == self.rootName and \
                     (justNodes is None or parNode in justNodes):
            if colorMap is not None:
              color = mapToColorScale(colorMap[parent])
            retStr += " node: { title:\"" + parent + "\" info2:\"" + parent + \
               "\" info1:\"" +drawStr + ";\" label:\"" + parId + \
               "\" color: " + color +" shape: hexagon bordercolor: " + \
               bordercolor + " borderstyle: " + borderstyle + " }\n"
            namesDrawn[parent] = True
          elif justNodes is None or parNode in justNodes:
            if colorMap is not None:
              color = mapToColorScale(colorMap[parent])
            retStr += " node: { title:\"" + parent + "\" info2:\"" + parent + \
               "\" info1:\"" +drawStr + ";\" label:\"" + parId + \
               "\" color: " + color +" bordercolor: " + bordercolor + \
               " borderstyle: " + borderstyle + " }\n"
            namesDrawn[parent] = True
        for child in children:
          if child not in namesDrawn:
            drawStr = outputDrawStr(self.label, self.nameToNode[child])
            childNode = self.nameToNode[child]
            childId = childNode.getId()
            if justNodes is None or childNode in justNodes:
              if colorMap is not None:
                color = mapToColorScale(colorMap[child])
              retStr += " node: { title:\"" + child + "\" info2:\"" + child + \
                "\" info1:\"" +drawStr + ";\" label:\"" + childId + \
                "\" color: " + color +" bordercolor: "+ bordercolor + \
                " borderstyle: " + borderstyle + " }\n"
              namesDrawn[child] = True
          if edges and (justNodes is None):
            retStr += " edge: { source: \""+parent+"\" target: \"" +child+ \
                    "\" color: "+edgeColor+" thickness: 1 }\n"
      retStr += " }\n"
      return retStr

    def removeNonKeepers(self):
      '''removes the nodes not marked as keepers'''
      stack = [self.rootName]
      while len(stack) > 0:
        current = stack.pop()
        children = []
        if current in self.tree: #no children, no problems
          children = self.tree[current]
        newChildren = children[:] #copy
        childStack = children[:]
        while len(childStack) > 0:
          child = childStack.pop()
          if self.keepers[child]:
            stack.append(child)
          else:
            grandChildren = []
            if child in self.tree:
              grandChildren = self.tree[child]
              del self.tree[child]
            newChildren.remove(child)
            newChildren.extend(grandChildren)
            childStack.extend(grandChildren)
        self.tree[current] = newChildren

  def __init__(self, tmDataList):
    '''takes a list of tmData, constructs individual subgraphs for each. '''
    self.tmToSubgraph = {}
    self.tmDataList = tmDataList
    self.namesToNode = {}
    digits = len(str(len(tmDataList))) #how many digits necessary to represent
    for tmDataCount,tmData in enumerate(tmDataList):
      outTmStr = string.zfill(tmDataCount, digits) #so you can sort
      newSubgraph = self.subgraph(tmData.inputFileName, tmData.tree, \
                 "s"+outTmStr+"f"+tmData.inputFileName[:4], tmData.root)
      self.tmToSubgraph[tmData] = newSubgraph
      self.namesToNode.update(newSubgraph.nameToNode) #put into big dict
    #that's all, built basic dot structure

  def addConnections(self, columnListNames, threshold, remove=True, \
                     resIds=False):
    '''adds either all connections or sequential connections, removes nodes 
    with no matches (if desired)'''
    self.threshold = threshold # save it for later
    self.connections = []
    columnList = self.tmDataList[0].titlesToColumns(columnListNames)
    #just compare 1 with 2, 2 with 3, etc..
    tmDataLast = self.tmDataList[0]
    for tmDataCount,tmData in enumerate(self.tmDataList[1:]):
      subGraph1 = self.tmToSubgraph[tmDataLast]
      subGraph2 = self.tmToSubgraph[tmData]
      tempConns = matchTrees(tmDataLast, tmData, columnList, threshold,resIds)
      for connection in tempConns:
        if connection[2] <= self.threshold:
          newConn = [ subGraph1.nodeToName[connection[0]], \
                      subGraph2.nodeToName[connection[1]], connection[2] ]
          subGraph1.keepers[newConn[0]] = True
          subGraph2.keepers[newConn[1]] = True
          self.connections.append(newConn)
      tmDataLast = tmData #update for next round
    if remove: #sometimes may not want to do this
      for subGraph in self.tmToSubgraph.values():
        subGraph.removeNonKeepers()

  def computeSearchConnections(self, threshold, columnList, columnsToMean, \
                           columnsToStddev, \
                           resColNums, sizeCol, selfScores, sameStruct, \
                           sizeMin=-1,sizeMax=10000000000, doSelfScore=True, \
                           justNodes=None, returnMatrix=False, alpha=1.):
    '''computes the entire matrix of possible connections amongst nodes, 
    lots of options that need documenting'''
    self.threshold = threshold # save it for later
    matchList = []
    if returnMatrix:
      matchDict = {}
    for tmDataCount, tmData in enumerate(self.tmDataList):
      for tmDataCount2, tmData2 in enumerate(self.tmDataList):
        if (sameStruct and tmDataCount == tmDataCount2) or  \
                           tmDataCount < tmDataCount2:
          for tmNodeCount, node1 in enumerate(tmData.tree.keys()):
            if justNodes is None or node1 in justNodes:
              for tmNode2Count, node2 in enumerate(tmData2.tree.keys()):
                if (tmDataCount<tmDataCount2 or tmNodeCount<tmNode2Count) and \
                   node1.attributes[sizeCol] < sizeMax and \
                   node2.attributes[sizeCol] < sizeMax and \
                   node1.attributes[sizeCol] > sizeMin and \
                   node2.attributes[sizeCol] > sizeMin and \
                   (justNodes is None or node2 in justNodes):
                  #actually want to compare these 2 nodes
                  score = compareColumns(node1, node2, columnList, \
                                         columnsToMean, columnsToStddev)
                  okayWithin = True
                  if score < threshold:
                    if tmDataCount == tmDataCount2: #same protein struct/tree
                      score2 = compareColumnsResidues(node1, node2, resColNums)
                      if score2 < 100.: #if there is overlap within same protein
                        okayWithin = False
                    if okayWithin: #no res overlap
                      #print score, selfScores[node1], selfScores[node2] 
                      totalScore = score #normal distance function
                      if doSelfScore: #do this instead
                        totalScore = -1./score
                        totalScore += alpha*(selfScores[node1] + selfScores[node2])
                      matchList.append((tmData,tmData2,node1,node2,totalScore))
                      if returnMatrix:
                        if node1 not in matchDict:
                          matchDict[node1] = {}
                        matchDict[node1][node2] = totalScore
    self.matchList = matchList
    if not returnMatrix:
      return matchList
    else: #return the matrix, turn dict into matrix
      tooBigNodes = []
      for tmData in self.tmDataList:
        for node1 in tmData.tree.keys():
          if node1.attributes[sizeCol] > sizeMax:
            tooBigNodes.append(node1)
      rows = matchDict.keys()
      colsSet = set()
      for eachRow in matchDict.values():
        colsSet.update(eachRow.keys())
      cols = list(colsSet)
      matchMatrix = []
      for row in rows:
        newRow = []
        for col in cols:
          newRow.append(matchDict[row][col])
        matchMatrix.append(newRow)
      return rows, cols, matchMatrix, tooBigNodes

  def refinePockets(self, columnList, columnsToMean, columnsToStddev, \
                           resColNums, sizeCol, selfScores, \
                           sizeMin=-1,sizeMax=10000000000, doSelfScore=True, \
                           justNodes=None, matchList=None, possNodes=None, 
                           notTheseNodes=None):
    '''finds the pocket with the highest mean difference, tries to find a better
    pocket for that tree, iterate until optimum?'''
    #go through matchlist, max lists of dists for each node
    if notTheseNodes is None:
      notTheseNodes = []
    justNodesToScores = {}
    nodeToTree = {}
    for aMatch in matchList:
      tmData1,tmData2,node1,node2,score = aMatch #unpack
      for aNode in (node1,node2):
        try:
          justNodesToScores[aNode].append(score) 
        except KeyError:
          justNodesToScores[aNode] = [score] 
      nodeToTree[node1] = tmData1
      nodeToTree[node2] = tmData2
    justNodesMeans = []
    for justNode in justNodesToScores.keys():
      if justNode not in notTheseNodes:
        meanScore = statistics.computeAverage(justNodesToScores[justNode])
        justNodesMeans.append((justNode, meanScore))
    if len(justNodesMeans) == 0: #nothing to do
      return False, None
    justNodesMeans.sort(key=operator.itemgetter(1)) #sort to find biggest
    worstNode, worstScore = justNodesMeans[-1] #unpack the one to fix
    worstTree = nodeToTree[worstNode]
    newPossNodes = possNodes[worstTree]
    possNodeScores = {}
    for possNode in newPossNodes:
      if possNode.attributes[sizeCol] < sizeMax and \
                                         possNode.attributes[sizeCol] > sizeMin:
        possNodeScores[possNode] = []
        for tmData in self.tmDataList:
          if tmData != worstTree: #don't do this tree against itself
            justNode = justNodes[tmData]
            score = compareColumns(justNode, possNode, columnList, \
                                         columnsToMean, columnsToStddev)
            totalScore = score
            if doSelfScore:
              totalScore += min(selfScores[node1], selfScores[node2])
            possNodeScores[possNode].append(totalScore)
    possNodeMean = []
    for possNode in newPossNodes:
      if possNode.attributes[sizeCol] < sizeMax and \
                                         possNode.attributes[sizeCol] > sizeMin:
        meanScore = statistics.computeAverage(possNodeScores[possNode])
        possNodeMean.append((possNode, meanScore))
    possNodeMean.sort(key=operator.itemgetter(1))
    newBestNode = possNodeMean[0][0]
    if worstNode != newBestNode:
      justNodes[worstTree] = newBestNode
      return True, worstNode #indicates change
    else:
      return False, worstNode #the node was already the best

  def getClustersConnections(self):
    '''using the data in self.connections, count the number of clusters.'''
    clusters = unionfind2.unionFind()
    for connection in self.connections:
      node1, node2 =  connection[0], connection[1] #unpack
      clusters.union(node1, node2) #really is that simple
    return clusters.toLists()

  def addSearchConnections(self, totalThreshold, remove=False, mst=False, \
                           maxConnCount=100000000000, lineMst=False, \
                           startNode=None, endNode=None, clusterOutput=False):
    '''adds the connections to self.connections if they meet the requirements'''
    tempConns = self.matchList[:] #copy and destroy possibly
    self.connections = []
    if clusterOutput:
      clusters = unionfind2.unionFind()
      overlapFunction = self.tmDataList[0].compareResidueIdentityMultipleNodes
      overlapCache = {}
      treeCountCache = {}
    for tmData in self.tmDataList:
      self.tmToSubgraph[tmData].resetKeepers()     
    if mst: #init this data structure
      mstUF = unionfind2.unionFind()
    if lineMst: #init this data structure
      mstUF = unionfind2.unionFind()
      connsLimit2 = {}
      if startNode is not None and endNode is not None: #limit endpoints
        connsLimit2[startNode] = [endNode]  #for lineMstEnds given hints
    tempConns.sort(key=operator.itemgetter(4)) #best first
    for aMatch in tempConns:
      tmData,tmData2,node1,node2,totalScore = aMatch #unpack
      mstOkay = (not mst) and (not lineMst) #iff both false, everything is fine 
      if mst: #do checks for mst
        if mstUF.different(node1, node2): #calls find on node1+2 to init them
          mstOkay = True
          mstUF.union(node1, node2)
      elif lineMst: #if okay to mst might check for linemst
        if mstUF.different(node1, node2): #calls find on node1+2 to init them
          if (node1 not in connsLimit2 or len(connsLimit2[node1]) == 1) and \
                     (node2 not in connsLimit2 or len(connsLimit2[node2]) == 1):
            #only now we know it is completely okay
            if node1 not in connsLimit2:
              connsLimit2[node1] = [] 
            connsLimit2[node1].append(node2)
            if node2 not in connsLimit2:
              connsLimit2[node2] = [] 
            connsLimit2[node2].append(node1)
            mstOkay = True
            mstUF.union(node1, node2)
      if mstOkay: #means either everything is fine or not mst
        if totalScore < totalThreshold and len(self.connections) < maxConnCount:
          subGraph1 = self.tmToSubgraph[tmData]
          subGraph2 = self.tmToSubgraph[tmData2]
          newConn = [ subGraph1.nodeToName[node1], \
                      subGraph2.nodeToName[node2], totalScore, \
                      node1, node2 ]
          subGraph1.keepers[newConn[0]] = True
          subGraph2.keepers[newConn[1]] = True
          self.connections.append(newConn)
          if clusterOutput:
            clusters.union(node1, node2) #really is that simple
            clustLists = clusters.toLists()
            for aCluster in clustLists:
              aCluster.sort()
              tupleCluster = tuple(aCluster)
              if tupleCluster not in overlapCache:
                aOverlap = overlapFunction(aCluster)
                overlapCache[tupleCluster] = aOverlap
              else:
                aOverlap = overlapCache[tupleCluster]
              if tupleCluster not in treeCountCache:
                treeSet = set()
                for node in aCluster:
                  treeSet.add(node.tree)
                treeCount = len(treeSet)
                treeCountCache[tupleCluster] = treeCount
              else:
                treeCount = treeCountCache[tupleCluster]
              if aOverlap >= 0.0 or len(clustLists) < 5:
                print len(self.connections), len(clustLists), 
                print "len:", len(aCluster), "over:", aOverlap, "count:", treeCount,
                print outputDrawStr(aCluster[0].tree.inputFileName, aCluster[0])
    if remove: #sometimes may not want to do this here
      for subGraph in self.tmToSubgraph.values():
        subGraph.removeNonKeepers()

  def decideWidth(self, value, threshold, minV=0.5, maxV=4):
    '''decides how wide the line for the connections should be, lower=thicker'''
    closeToMax = (threshold - value) / threshold
    return max(min(minV, minV + closeToMax * (maxV-minV)), maxV)

  def write(self, outputFile):
    '''writes the dot file'''
    outFile = open(outputFile, 'w')
    outFile.write("graph G {\n")
    #first write an invisible node and connect it to each subgraph's root
    outFile.write("head [style=invisible];\n")
    outFile.write("node [shape=box,height=0,margin=0,0];\n")
    for count, subgraph in enumerate(self.tmToSubgraph.values()):
      outFile.write("head -- " + subgraph.rootName + " [style=invisible, len=1];\n")
      if count < len(self.tmToSubgraph) - 1: #isn't last
        outFile.write(subgraph.rootName + "space [style=invisible];\n")
        outFile.write("head -- " + subgraph.rootName+"space [style=invisible];\n")
        outFile.write(subgraph.rootName + "space2 [style=invisible];\n")
        outFile.write("head -- " + subgraph.rootName+"space2 [style=invisible];\n")
    for subgraph in self.tmToSubgraph.values():
      outFile.write(str(subgraph))
    #now write connections between trees
    for connection in self.connections:
      width = self.decideWidth(connection[2], self.threshold)
      label = " " +str(int(connection[2])) #no need for all those decimal points
      outFile.write(connection[0] + " -- " + connection[1] + " [label = \"" + \
                    label + "\" penwidth=" + str(width) + "];\n")
    #emphasize nodes with matches
    for connection in self.connections:
      outFile.write(connection[0] + "[shape=egg,height=0,margin=0,0];\n")
      outFile.write(connection[1] + "[shape=egg,height=0,margin=0,0];\n")
    outFile.write("}\n")
    outFile.close()

  def tempRemoveKeepersWriteGdl(self, outputFile, force=True, edges=True, \
                                justNodes=None):
    '''backs up the tree structure, removes non keepers, writes gdl, restore
    tree structure so it can be used again'''
    backupTrees = {}
    for subGraph in self.tmToSubgraph.values():
      backupTrees[subGraph] = subGraph.tree.copy()
      subGraph.removeNonKeepers()
    self.writeGdl(outputFile, force, edges=edges, justNodes=justNodes, \
                  colorByResidueTanimoto=True)
    #self.writeGdl(outputFile + ".notreeedges.gdl", force=force, edges=False)
    for subGraph in self.tmToSubgraph.values():
      subGraph.tree = backupTrees[subGraph]

  def writeGdl(self, outputFile, force=True, edges=True, justNodes=None, 
               colorByResidueTanimoto=False):
    '''writes a GDL file. if the colorByResidueTanimoto flag is true, then 
    for each node, all connections are examined and the average tanimoto score
    is used to color the node instead of the color of the tree'''
    if colorByResidueTanimoto:
      otherNodesDict = {} #key on name, then node, then other nodes
      for conn in self.connections:
        if conn[0] not in otherNodesDict:
          otherNodesDict[conn[0]] = [conn[3]]
        otherNodesDict[conn[0]].append(conn[4])
        if conn[1] not in otherNodesDict:
          otherNodesDict[conn[1]] = [conn[4]]
        otherNodesDict[conn[1]].append(conn[3])
      nodeToMeanTanimoto = collections.defaultdict(float) #defaults to 0 overlap
      compareResFunction = self.tmDataList[0].compareResidueIdentityNodes
      for oneNodeName in otherNodesDict.iterkeys():
        totalTanimoto = 0.
        oneNode = otherNodesDict[oneNodeName][0]
        otherNodes = otherNodesDict[oneNodeName][1:]
        for otherNode in otherNodes:
          totalTanimoto += compareResFunction(oneNode, otherNode)
        nodeToMeanTanimoto[oneNodeName] = totalTanimoto/float(len(otherNodes))
        #if len(otherNodes) < 2: #require at least 5 connections 
        #  nodeToMeanTanimoto[oneNodeName] = 0. 
      #if len(nodeToMeanTanimoto) > 0:
      #  print max(nodeToMeanTanimoto.values())
      #else:
      #  print 0.
    outFile = open(outputFile, 'w')
    outFile.write("graph : {\n")
    outFile.write("edge.arrowstyle: none\n")
    if force:
      outFile.write("layoutalgorithm : forcedir\n")
      outFile.write("attraction   : 40\n")
      outFile.write("repulsion    : 99\n")
      outFile.write("gravity      : 0.0\n")
      outFile.write("fdmax        : 500\n")
      outFile.write("tempmax      : 254\n")
      outFile.write("temptreshold : 3\n")
      outFile.write("tempscheme   : 3\n")
      outFile.write("tempfactor   : 1.08\n")
      outFile.write("randomfactor : 100\n")
    if justNodes is None:
      outFile.write("node: {title: \"head\" invisible: yes}\n")
      for count, tmTree in enumerate(self.tmDataList):
        subgraph = self.tmToSubgraph[tmTree]
        outFile.write("edge: {source: \"head\" target: \"" + \
                      subgraph.rootName + "\" linestyle:invisible }\n")
    for sgnum, tmTree in enumerate(self.tmDataList):
      subgraph = self.tmToSubgraph[tmTree]
      if colorByResidueTanimoto:
        outFile.write(subgraph.gdlStr(sgnum, edges=edges, justNodes=justNodes, \
                      colorMap=nodeToMeanTanimoto))
      else: #normal coloring by tree
        outFile.write(subgraph.gdlStr(sgnum, edges=edges, justNodes=justNodes))
    #now write connections between trees
    for connection in self.connections:
      width = self.decideWidth(connection[2], self.threshold)
      #label = " " +str(int(connection[2])) #no need for all those decimal points
      #outFile.write(connection[0] + " -- " + connection[1] + " [label = \"" + \
      #              label + "\" penwidth=" + str(width) + "];\n")
      outFile.write("edge: { source: \"" + connection[0] + "\" target: \"" + \
                    connection[1] + "\" thickness: " + str(width) + \
                    " priority: " + str(width) + "}\n")
    outFile.write("}\n")
    outFile.close()

  def findChains(self):
    '''find and print chains that go from tree to tree'''
    print "tstMultiOpen ",
    for tm3tree in self.tmDataList:
      end = tm3tree.inputFileName.find(".tree")
      print tm3tree.inputFileName[:end] + ", ",
    print " "
    longestChain = len(self.tmToSubgraph)
    connDict = {}
    for connection in self.connections: #rely on 0->1 ordering being consistent
      if connection[0] not in connDict:
        connDict[connection[0]] = []
      connDict[connection[0]].append(connection[1])    
    for key in connDict.keys(): #start compression or putting into bigger list
      listChains = connDict[key]
      if len(listChains) == 1: #haven't seen before and potential to grow
        while listChains[-1] in connDict:
          listChains.extend(connDict[listChains[-1]])
      connDict[key] = listChains
    realChains = []
    for key,listChains in connDict.iteritems():
      if len(listChains) + 1 == longestChain:
        realChains.append([key])
        realChains[-1].extend(listChains)
    #use data is self.namesToNode to print out volume changes?
    printCols = self.tmToSubgraph.keys()[0].titlesToColumns(["Volume","meanTD"])
    for chain in realChains:
      print "start :",
      for item in chain:
        actualNode = self.namesToNode[item]
        for printCol in printCols:
          print str(round(float(actualNode.attributes[printCol]), 1)) + " ",
        print " -- ",
      print " : " 
      print "tstMultiDraw groups=[",
      for item in chain[:-1]:
        actualNode = self.namesToNode[item]
        print actualNode.getId()+ ", ",            
      actualNode = self.namesToNode[chain[-1]]
      print actualNode.getId(),            
      print "]" 

class dotResidues(dot):
  '''specialized class for finding residue overlaps'''

  def __init__(self, tmDataList, residueList, threshold=50, \
               colNameIn="Residue Name List", penalty=1., printHelpful=True):
    '''adds the many tmdata, finds pocket with best residue overlap in each,
    connect them and output many things'''
    self.tmToSubgraph = {}
    self.tmDataList = tmDataList
    self.threshold = threshold # save it for later
    self.namesToNode = {}
    for tmDataCount,tmData in enumerate(tmDataList):
      newSubgraph = self.subgraph(tmData.inputFileName, tmData.tree, \
                 "s"+str(tmDataCount)+"f"+tmData.inputFileName[:4], tmData.root)
      self.tmToSubgraph[tmData] = newSubgraph
      self.namesToNode.update(newSubgraph.nameToNode) #put into big dict
    self.connections = []
    self.treeToBestNode = {}
    self.treeToPossNodes = {}
    self.bestList = []
    self.bestScores = []
    for tmDataCount,tmData in enumerate(tmDataList):
      bestNode, bestScore, possNodes = tmData.findBestResidueListMatch(\
                   residueList, colName=colNameIn, penalty=penalty)
      self.treeToBestNode[tmData] = bestNode
      self.treeToPossNodes[tmData] = possNodes
      self.bestList.append(bestNode)
      self.bestScores.append(bestScore)
    tmDataLast = tmDataList[0]
    for tmDataCount,tmData in enumerate(tmDataList[1:]):
      subGraph1 = self.tmToSubgraph[tmDataLast]
      subGraph2 = self.tmToSubgraph[tmData]
      best1 = self.treeToBestNode[tmDataLast]
      best2 = self.treeToBestNode[tmData]
      newConn = [ subGraph1.nodeToName[best1], \
                      subGraph2.nodeToName[best2], 0. ] #fake size
      subGraph1.keepers[newConn[0]] = True
      subGraph2.keepers[newConn[1]] = True
      self.connections.append(newConn)
      tmDataLast = tmData #update for next round
    if printHelpful:
      self.printHelpful()

  def setBestNodes(self, newBestNodes):
    '''resets the best nodes to the new dictionary'''
    self.treeToBestNode = newBestNodes
    self.bestList = []
    for tmData in self.tmDataList:
      self.bestList.append(newBestNodes[tmData])

  def printHelpful(self, prefix="temp"):
    '''prints helpful output for debugging visually with pymol'''
    print "tstMultiOpen ",
    for tm3tree in self.tmDataList[:-1]:
      end = tm3tree.inputFileName.find(".tree")
      print tm3tree.inputFileName[:end] + ", ",
    for tm3tree in [self.tmDataList[-1]]:
      end = tm3tree.inputFileName.find(".tree")
      print tm3tree.inputFileName[:end] 
    #printCols = self.tmToSubgraph.keys()[0].titlesToColumns(["Volume","meanTD"])
    #print "start :",
    #for item in self.bestList:
    #  actualNode = item
    #  for printCol in printCols:
    #    print str(round(float(actualNode.attributes[printCol]), 1)) + " ",
    #  print " -- ",
    #print " : "
    print "tstMultiDraw groups=[",   
    for item in self.bestList[:-1]:
      actualNode = item
      print actualNode.getId()+ ", ",
    actualNode = self.bestList[-1]
    print actualNode.getId(),
    print "]"
    outFile = open(prefix + ".tab.txt", 'w')
    outFile.write("tm3file\t")
    outFile.write(self.tmDataList[0].titleRow() + "\n")
    for tmDataCount,tmData in enumerate(self.tmDataList):
       bestNode = self.treeToBestNode[tmData]
       outFile.write(tmData.inputFileName + "\t")
       outFile.write(str(bestNode) + "\n")
    outFile.close()
    self.writeAnimation(outFileName=prefix + ".animation.py")
    self.writeAnimation(bb=True, outFileName=prefix + ".animation_bb.py")

  def writeAnimation(self, outFileName="temp.animation.py", \
                     width=3, bb=False, maxTD=True):
    '''uses data in tmDataList and bestList to output an animation script'''
    if maxTD:
      maxMaxTD = 0.
      for tm3tree in self.tmDataList:
        thisMax = tm3tree.getMaxTravelDepth()
        maxMaxTD = max(thisMax, maxMaxTD)
      outMaxString = ",maxval="+str(maxMaxTD)
    outFile = open(outFileName, 'w')
    outFile.write("from pymol import cmd\n")
    for index, tm3tree in enumerate(self.tmDataList):
      end = tm3tree.inputFileName.find(".tree")
      tstName = tm3tree.inputFileName[:end] 
      pdbName = string.replace(tstName, ".cav.tst", ".pdb") #this should work
      pdbName = string.replace(pdbName, ".nocav.tst", ".pdb") #or this
      pocketId = self.bestList[index].getId()
      outFile.write("cmd.do(\"")
      outFile.write("delete all")
      outFile.write("\")\ncmd.do(\"")
      outFile.write("load " + pdbName)
      outFile.write("\")\ncmd.do(\"")
      outFile.write("color grey")
      outFile.write("\")\ncmd.do(\"")
      if not bb:
        outFile.write("show sticks")
      else:
        outFile.write("hide all")
        outFile.write("\")\ncmd.do(\"")
        outFile.write("show ribbon")
      outFile.write("\")\ncmd.do(\"")
      outFile.write("tstOpenDraw " + tstName + ", group=" + str(pocketId) )
      if maxTD:
        outFile.write(outMaxString)
      outFile.write("\")\ncmd.do(\"")
      outFile.write("ray")
      outFile.write("\")\ncmd.do(\"")
      if not bb:
        outFile.write("png temp." + str(index).zfill(width) + "." + tstName + \
                               ".png")
      else:
        outFile.write("png temp.bb." + str(index).zfill(width) + "." + tstName \
                             + ".png")
      outFile.write("\")\n")
    outFile.close()

class dotRoot(dotResidues):
  '''specialized class for finding residue overlaps'''

  def __init__(self, tmDataList, printHelpful=True):
    '''adds the many tmdata, finds pocket with best residue overlap in each,
    connect them and output many things'''
    self.tmToSubgraph = {}
    self.tmDataList = tmDataList
    self.namesToNode = {}
    for tmDataCount,tmData in enumerate(tmDataList):
      newSubgraph = self.subgraph(tmData.inputFileName, tmData.tree, \
                 "s"+str(tmDataCount)+"f"+tmData.inputFileName[:4], tmData.root)
      self.tmToSubgraph[tmData] = newSubgraph
      self.namesToNode.update(newSubgraph.nameToNode) #put into big dict
    self.connections = []
    self.treeToBestNode = {}
    self.treeToPossNodes = {}
    self.bestList = []
    for tmDataCount,tmData in enumerate(tmDataList):
      self.treeToBestNode[tmData] = tmData.root #just root
      self.treeToPossNodes[tmData] = tmData.root
      self.bestList.append(tmData.root)
    tmDataLast = tmDataList[0]
    for tmDataCount,tmData in enumerate(tmDataList[1:]):
      subGraph1 = self.tmToSubgraph[tmDataLast]
      subGraph2 = self.tmToSubgraph[tmData]
      best1 = self.treeToBestNode[tmDataLast]
      best2 = self.treeToBestNode[tmData]
      newConn = [ subGraph1.nodeToName[best1], \
                      subGraph2.nodeToName[best2], 0. ] #fake size
      subGraph1.keepers[newConn[0]] = True
      subGraph2.keepers[newConn[1]] = True
      self.connections.append(newConn)
      tmDataLast = tmData #update for next round
    if printHelpful:
      self.printHelpful()

def compareColumns(tmNode1, tmNode2, columnList, columnsToMean, \
                   columnsToStddev):
  '''compares 2 nodes,  report dist between vectors of z-scores
  '''
  zScores1, zScores2 = [],[]
  for col in columnList:
    colMean = columnsToMean[col]
    colStddev = columnsToStddev[col]
    zScore1 = (tmNode1.attributes[col]-colMean)/colStddev
    zScore2 = (tmNode2.attributes[col]-colMean)/colStddev
    zScores1.append(zScore1)
    zScores2.append(zScore2)
  #distZ = geometry.distL2(zScores1, zScores2)  
  distZ = geometry.dist(zScores1, zScores2, metric='L1')   #L1 seems to be best
  return distZ

def compareColumnsResidues(tmNode1, tmNode2, columnList):
  '''
  if resIds is true, the column is a string of resnums+chain, then do int/union.
  '''
  resIds1 = set(string.split(tmNode1.attributes[columnList[0]],"+"))
  resIds2 = set(string.split(tmNode2.attributes[columnList[0]],"+"))
  resIds1.remove('') #how did these get here? stupid string methods
  resIds2.remove('')
  union = float(len(resIds1.union(resIds2)))
  if union > 0.: #otherwise automatically no match
    #print resIds1, resIds2, union , len(resIds1.intersection(resIds2))
    return 100 - (float(len(resIds1.intersection(resIds2)))/union) * 100. #want complete match to be equal to 0
  else:
    return 100 #no intersection, worst score possible

def matchTrees(tmData1, tmData2, columnList, threshold, resIds=False):
  retList = []
  if not resIds:
    columnsToMean, columnsToStddev = tm3.calcColumnsMeanStddev(columnList, \
                                                           [tmData1, tmData2])
  revData = False
  if len(tmData1.tree.keys()) < len(tmData2.tree.keys()):
    tmData1, tmData2 = tmData2, tmData1
    revData = True
  cost = []
  names1, names2 = [],[]
  for node2 in tmData2.tree.keys():
    names2.append(node2)
  for node1 in tmData1.tree.keys():
    names1.append(node1)
    costRow = []
    for node2 in tmData2.tree.keys():
      if not resIds:
        score = compareColumns(node1, node2, columnList, \
                  columnsToMean, columnsToStddev)
      else:
        score = compareColumnsResidues(node1, node2, columnList)
      costRow.append(score)
    cost.append(costRow)
  matches = munkreskuhn.assignAndReturnMatches(cost)
  returnConns = []
  for match in matches:
    if revData:
      returnConns.append((names2[match[1]], names1[match[0]], match[2]))
    else:
      returnConns.append((names1[match[0]], names2[match[1]], match[2]))
  return returnConns

if -1 != string.find(sys.argv[0], "dot.py"):
  pass #actually do nothing
