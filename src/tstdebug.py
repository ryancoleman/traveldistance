
#ryan.g.coleman@gmail.com ryangc@mail.med.upenn.edu
#Ryan Coleman, Kim Sharp, univ. of penn., http://crystal.med.upenn.edu
#contains various debugging functions, usually these write a pymol script

from grid import getIndices

def debugTravelGrid(grid, filename, maxTD=5.0):
  fileTemp = open(filename, 'w')
  fileTemp.write("from pymol.cgo import * \nfrom pymol import cmd\nsurfObj = [ ")
  for oneX in grid:
    for oneY in oneX:
      for onePoint in oneY:
        if onePoint[0] != -1 and onePoint[0] != -2:
          fileTemp.write("COLOR, ")
          if onePoint[0] == -2:
            fileTemp.write("0.1, ")
            fileTemp.write("0.1, ")
            fileTemp.write("0.1, ")
          else:
            depth = onePoint[0]
            if depth < 0:
              depth = -2-depth
            fileTemp.write(str(1.-(depth/maxTD)))
            fileTemp.write(", ")
            fileTemp.write(str((depth/maxTD)))
            fileTemp.write(", ")
            fileTemp.write("0.01, ")
          fileTemp.write("SPHERE, ")
          for coord in onePoint[1:]:
            fileTemp.write(str(coord) + ", ")
          fileTemp.write("0.2, \n")
  fileTemp.write(" ]\n")
  fileTemp.write("cmd.load_cgo(surfObj,'grid" +filename + "')\n")
  fileTemp.close()

def debugTravelSurfGrid(grid, filename, extraEdges, mins,maxs,gridSize, maxTD=5.):
  fileTemp = open(filename, 'w')
  fileTemp.write("from pymol.cgo import * \nfrom pymol import cmd\nsurfObj = [ ")
  for oneX in grid:
    for oneY in oneX:
      for onePoint in oneY:
        pointIndex = grid.getIndices(mins,gridSize,onePoint[1:])
        if extraEdges.has_key(pointIndex):
          fileTemp.write("COLOR, ")
          if onePoint[0] == -2:
            fileTemp.write("0.1, ")
            fileTemp.write("0.1, ")
            fileTemp.write("0.1, ")
          else:
            depth = onePoint[0]
            if depth < 0:
              depth = -2-depth
            fileTemp.write(str(1.-(depth/maxTD)))
            fileTemp.write(", ")
            fileTemp.write(str((depth/maxTD)))
            fileTemp.write(", ")
            fileTemp.write("0.01, ")
          fileTemp.write("SPHERE, ")
          for coord in onePoint[1:]:
            fileTemp.write(str(coord) + ", ")
          fileTemp.write("0.2, \n")
  fileTemp.write(" ]\n")
  fileTemp.write("cmd.load_cgo(surfObj,'grid" +filename + "')\n")
  fileTemp.close()

def debugGridCountVals(grid):
  output = {}
  for indexX,rowX in enumerate(grid):
    for indexY,rowY in enumerate(rowX):
      for indexZ,entryZ in enumerate(rowY):
        if entryZ[0] in output:
          output[entryZ[0]] += 1
        else:
          output[entryZ[0]] = 1
  vals = output.keys()
  vals.sort()
  for val in vals:
    print str(val) + ": " + str(output[val]) + ", ",
  print " "

#debugging in macpymol, shows the grid
def debugGrid(grid, filename):
  fileTemp = open(filename, 'w')
  fileTemp.write("from pymol.cgo import * \nfrom pymol import cmd\nsurfObj = [ ")
  for oneX in grid:
    for oneY in oneX:
      for onePoint in oneY:
        fileTemp.write("COLOR, ")
        if onePoint[0] == -1:
          fileTemp.write("0.1, ")
          fileTemp.write("0.1, ")
          fileTemp.write("0.95, ")
        elif onePoint[0] == 0:
          fileTemp.write("0.1, ")
          fileTemp.write("0.95, ")
          fileTemp.write("0.1, ")
        elif onePoint[0] == -2:
          fileTemp.write("0.95, ")
          fileTemp.write("0.1, ")
          fileTemp.write("0.1, ")
        fileTemp.write("SPHERE, ")
        for coord in onePoint[1:]:
          fileTemp.write(str(coord) + ", ")
        fileTemp.write("0.2, ")
  fileTemp.write(" ]\n")
  fileTemp.write("cmd.load_cgo(surfObj,'grid" +filename + "')\n")
  fileTemp.close()

def debugGridContours(grid, filename):
  fileTemp = open(filename, 'w')
  fileTemp.write("from pymol.cgo import * \nfrom pymol import cmd\nsurfObj = [ ")
  for oneX in grid:
    for oneY in oneX:
      for onePoint in oneY:
        if onePoint[0] > 6.0:
          fileTemp.write("COLOR, ")
          fileTemp.write("0.1, ")
          fileTemp.write("0.1, ")
          fileTemp.write("0.95, ")
          fileTemp.write("SPHERE, ")
          for coord in onePoint[1:]:
            fileTemp.write(str(coord) + ", ")
          fileTemp.write("0.2, ")
  fileTemp.write(" ]\n")
  fileTemp.write("cmd.load_cgo(surfObj,'grid" +filename + "')\n")
  fileTemp.close()

#debugging in macpymol, shows the grid
def debugGridNoBlue(grid, filename):
  fileTemp = open(filename, 'w')
  fileTemp.write("from pymol.cgo import * \nfrom pymol import cmd\nsurfObj = [ ")
  for oneX in grid:
    for oneY in oneX:
      for onePoint in oneY:
        if onePoint[0] != -1:
          fileTemp.write("COLOR, ")
          if onePoint[0] == 0:
            fileTemp.write("0.1, ")
            fileTemp.write("0.95, ")
            fileTemp.write("0.1, ")
          elif onePoint[0] == -2:
            fileTemp.write("0.95, ")
            fileTemp.write("0.1, ")
            fileTemp.write("0.1, ")
          fileTemp.write("SPHERE, ")
          for coord in onePoint[1:]:
            fileTemp.write(str(coord) + ", ")
          fileTemp.write("0.2, ")
  fileTemp.write(" ]\n")
  fileTemp.write("cmd.load_cgo(surfObj,'grid" +filename + "')\n")
  fileTemp.close()

#debugging in macpymol, shows the grid
def debugGridJustGreen(grid, filename):
  fileTemp = open(filename, 'w')
  fileTemp.write("from pymol.cgo import * \nfrom pymol import cmd\nsurfObj = [ ")
  for oneX in grid:
    for oneY in oneX:
      for onePoint in oneY:
        if onePoint[0] == 0:
          fileTemp.write("COLOR, ")
          if onePoint[0] == 0:
            fileTemp.write("0.1, ")
            fileTemp.write("0.95, ")
            fileTemp.write("0.1, ")
          fileTemp.write("SPHERE, ")
          for coord in onePoint[1:]:
            fileTemp.write(str(coord) + ", ")
          fileTemp.write("0.2, ")
  fileTemp.write(" ]\n")
  fileTemp.write("cmd.load_cgo(surfObj,'grid" +filename + "')\n")
  fileTemp.close()

#debugging in macpymol
def debugTris(line, triList, ptList, fileName):
  surfObj = "from pymol.cgo import * \nfrom pymol import cmd\nsurfObj = [ "
  for point in ptList:
    surfObj += "COLOR, "
    surfObj += "0.1, "
    surfObj += "0.1, "
    surfObj += "0.95, "
    surfObj += "SPHERE, "
    for coord in point:
      surfObj += str(coord) + ", "
    surfObj += "0.1, "
  surfObj += " BEGIN, TRIANGLES, "
  for tri in triList:
    for point in tri:
      surfObj += "COLOR, "
      surfObj += "1.0, "
      surfObj += "0.0, "
      surfObj += "0.5, "
      surfObj += "VERTEX, "
      for coord in point:
        surfObj += str(coord) + ", "
  surfObj += "END, BEGIN, LINES, "
  for point in line:
    surfObj += "COLOR, "
    surfObj += "1.0, "
    surfObj += "0.5, "
    surfObj += "0.0, "
    surfObj += "VERTEX, "
    for coord in point:
      surfObj += str(coord) + ", "
  surfObj += "END ]\n"
  surfObj += "cmd.load_cgo(surfObj,'surf" +str(line[0][0]) + "." + str(line[0][1]) + "')\n"
  fileOut = open(fileName, 'w')
  fileOut.write(surfObj)
  fileOut.close()

def debugTriangles(triList, allTris, ptXYZ, fileName, ptList=False, \
                   triColor=(1.,0.,.5)):
  surfObj = "from pymol.cgo import * \nfrom pymol import cmd\nsurfObj = [ "
  surfObj += " BEGIN, TRIANGLES, "
  for tri in triList:
    for point in allTris[tri-1][1:]:
      surfObj += "COLOR, "
      surfObj += str(triColor[0]) + ", "
      surfObj += str(triColor[1]) + ", "
      surfObj += str(triColor[2]) + ", "
      surfObj += "VERTEX, "
      for coord in ptXYZ[point-1][1:]:
        surfObj += str(coord) + ", "
  surfObj += "END ]\n"
  surfObj += "cmd.load_cgo(surfObj,'tris" + fileName + "')\n"
  surfObj += "\nsurfObj2 = [ "
  surfObj += "LINEWIDTH, 3,  "
  surfObj += "BEGIN, LINES, "
  lastPt = ptList[len(ptList)-1]
  for point in ptList:
    surfObj += "COLOR, "
    surfObj += "0.0, "
    surfObj += "0.5, "
    surfObj += "0.9, "
    surfObj += "VERTEX, "
    for coord in ptXYZ[lastPt-1][1:]:
      surfObj += str(coord) + ", "
    surfObj += "VERTEX, "
    for coord in ptXYZ[point-1][1:]:
      surfObj += str(coord) + ", "
    lastPt = point
  surfObj += "END ]\n"
  surfObj += "cmd.load_cgo(surfObj2,'lines" + fileName + "')\n"
  fileOut = open(fileName, 'w')
  fileOut.write(surfObj)
  fileOut.close()

def debugTrianglesPtList(triList, allTris, ptXYZ, fileName, ptListList=False, \
                   triColor=(1.,0.,.5)):
  surfObj = "from pymol.cgo import * \nfrom pymol import cmd\nsurfObj = [ "
  surfObj += " BEGIN, TRIANGLES, "
  for tri in triList:
    for point in allTris[tri-1][1:]:
      surfObj += "COLOR, "
      surfObj += str(triColor[0]) + ", "
      surfObj += str(triColor[1]) + ", "
      surfObj += str(triColor[2]) + ", "
      surfObj += "VERTEX, "
      for coord in ptXYZ[point-1][1:]:
        surfObj += str(coord) + ", "
  surfObj += "END ]\n"
  surfObj += "cmd.load_cgo(surfObj,'tris" + fileName + "')\n"
  for ptListIndex, ptList in enumerate(ptListList):
    surfObj += "\nsurfObjL"+ str(ptListIndex) + " = [ "
    surfObj += "LINEWIDTH, 3,  "
    surfObj += "BEGIN, LINES, "
    lastPt = ptList[len(ptList)-1]
    for point in ptList:
      surfObj += "COLOR, "
      surfObj += "0.0, "
      surfObj += "0.3, "
      surfObj += "0.9, "
      surfObj += "VERTEX, "
      for coord in ptXYZ[lastPt-1][1:]:
        surfObj += str(coord) + ", "
      surfObj += "VERTEX, "
      for coord in ptXYZ[point-1][1:]:
        surfObj += str(coord) + ", "
      lastPt = point
    surfObj += "END ]\n"
    surfObj += "cmd.load_cgo(surfObjL" + str(ptListIndex) + ",'lines" +str(ptListIndex)+ fileName + "')\n"
  fileOut = open(fileName, 'w')
  fileOut.write(surfObj)
  fileOut.close()

def debugTrianglesNotOrig(triList, ptXYZ, fileName, ptList=False):
  surfObj = "from pymol.cgo import * \nfrom pymol import cmd\nsurfObj = [ "
  surfObj += " BEGIN, TRIANGLES, "
  for tri in triList:
    for point in tri:
      surfObj += "COLOR, "
      surfObj += "1.0, "
      surfObj += "0.0, "
      surfObj += "0.5, "
      surfObj += "VERTEX, "
      for coord in ptXYZ[point-1][1:]:
        surfObj += str(coord) + ", "
  surfObj += "END ]\n"
  surfObj += "cmd.load_cgo(surfObj,'tris" + fileName + "')\n"
  surfObj += "\nsurfObj2 = [ "
  surfObj += "BEGIN, LINES, "
  lastPt = ptList[len(ptList)-1]
  for point in ptList:
    surfObj += "COLOR, "
    surfObj += "0.0, "
    surfObj += "0.5, "
    surfObj += "0.9, "
    surfObj += "VERTEX, "
    for coord in ptXYZ[lastPt-1][1:]:
      surfObj += str(coord) + ", "
    surfObj += "VERTEX, "
    for coord in ptXYZ[point-1][1:]:
      surfObj += str(coord) + ", "
    lastPt = point
  surfObj += "END ]\n"
  surfObj += "cmd.load_cgo(surfObj2,'lines" + fileName + "')\n"
  fileOut = open(fileName, 'w')
  fileOut.write(surfObj)
  fileOut.close()

def debugTriangleList(triListList, allTris, ptXYZ, fileName, ptListList=False):
  surfObj = "from pymol.cgo import * \nfrom pymol import cmd\n"
  top = float(len(triListList))
  for index, triList in enumerate(triListList):
    surfObj += "surfObj" + str(index) + " = [  BEGIN, TRIANGLES, "
    #topTri = float(len(triList))
    for indexTri,tri in enumerate(triList):
      for point in allTris[tri-1][1:]:
        surfObj += "COLOR, "
        surfObj += str(1.-float(float(index)/top))   #modify color somewhat
        surfObj += ", "
        surfObj += str(float(float(index)/top))   #modify color somewhat
        surfObj += ", "
        #surfObj += str(float(float(indexTri)/topTri))   #modify color somewhat
        surfObj += str(0.9)
        surfObj += ", "
        surfObj += "VERTEX, "
        for coord in ptXYZ[point-1][1:]:
          surfObj += str(coord) + ", "
    surfObj += "END ]\n"
  if ptListList:
    for index,ptList in enumerate(ptListList):
      surfObj += "\nsurfObjX" + str(index) + " = [ "
      surfObj += "BEGIN, LINES, "
      lastPt = ptList[len(ptList)-1]
      for point in ptList:
        surfObj += "COLOR, "
        surfObj += "0.0, "
        surfObj += "0.5, "
        surfObj += "0.9, "
        surfObj += "VERTEX, "
        for coord in ptXYZ[lastPt-1][1:]:
          surfObj += str(coord) + ", "
        surfObj += "VERTEX, "
        for coord in ptXYZ[point-1][1:]:
          surfObj += str(coord) + ", "
        lastPt = point
      surfObj += "END ]\n"
  for index, triList in enumerate(triListList):
    surfObj += "cmd.load_cgo(surfObj" + str(index) + ",'tris" + str(index) + "." + fileName + "')\n"
    if ptListList:
      if index < len(ptListList): #should be equal but just in case...
        surfObj += "cmd.load_cgo(surfObjX" + str(index) + " ,'lines" + str(index) + "." + fileName + "')\n"
  fileOut = open(fileName, 'w')
  fileOut.write(surfObj)
  fileOut.close()

def debugTracebacks(grid, tracebacks, filename, minDist=0.0):
  fileTemp = open(filename, 'w')
  fileTemp.write("from pymol.cgo import * \nfrom pymol import cmd\nsurfObj = [ ")
  for startList,end,dist in tracebacks.values():
    if dist >= minDist:
      for start in startList:
        fileTemp.write("CYLINDER, ")
        for gridPt in (start,end):
          for coord in grid[gridPt[0]][gridPt[1]][gridPt[2]][1:]:
            fileTemp.write(str(coord) + ", ")
        #radius
        fileTemp.write("0.1, ")
        #color
        fileTemp.write("1.0, ")
        fileTemp.write("0.1, ")
        fileTemp.write("0.1, ")
        fileTemp.write("0.5, ")
        fileTemp.write("0.6, ")
        fileTemp.write("0.6, ")
  fileTemp.write(" ]\n")
  fileTemp.write("cmd.load_cgo(surfObj,'" + str(minDist) + "." +  filename + "')\n")
  fileTemp.close()

def debugOneTraceback(grid, tracebacks, filename, ending):
  leftToDo = [ending] #queue
  fileTemp = open(filename, 'w')
  fileTemp.write("from pymol.cgo import * \nfrom pymol import cmd\nsurfObj = [ ")
  while len(leftToDo) > 0:
    #print leftToDo
    currentEnd = leftToDo.pop(0)
    if tracebacks.has_key(currentEnd):
      startList,end,dist = tracebacks[currentEnd]
      for start in startList:
        leftToDo.append(start)
        fileTemp.write("CYLINDER, ")
        for gridPt in (start,end):
          for coord in grid[gridPt[0]][gridPt[1]][gridPt[2]][1:]:
            fileTemp.write(str(coord) + ", ")
        #radius
        fileTemp.write("0.1, ")
        #color
        fileTemp.write("1.0, ")
        fileTemp.write("0.1, ")
        fileTemp.write("0.1, ")
        fileTemp.write("0.4, ")
        fileTemp.write("0.6, ")
        fileTemp.write("0.6, ")
  fileTemp.write(" ]\n")
  fileTemp.write("cmd.load_cgo(surfObj,'" + filename + "')\n")
  fileTemp.close()

def findDeepestTraceback(tracebacks):
  maxDist,maxEnd  = 0.0, False
  for startList,end,dist in tracebacks.values():
    if dist > maxDist:
      maxDist, maxEnd = dist,end
  return maxEnd

def debugDeepestTraceback(grid, tracebacks, filename):
  debugOneTraceback(grid,tracebacks,filename,findDeepestTraceback(tracebacks))

def debugSetGridSpheres(outputSet, gridSize, filename, endPoints=False, radius=False,mainColor=(0.0,.01,.9)):
  endPointsSet = set()
  if endPoints:
    for singleEndPoint in endPoints:
      endPointsSet.add(singleEndPoint[1:])
  #outputSet should be a list of grid'points' where 1:4 is the xyz coords
  #makes spheres at each grid center
  fileTemp = open(filename, 'w')
  fileTemp.write("from pymol.cgo import * \nfrom pymol import cmd\nsurfObj = [ ")
  for point in outputSet:
    fileTemp.write("COLOR, ")
    if tuple(point[1:]) in endPointsSet: #color it red
      fileTemp.write("0.95, ")
      fileTemp.write("0.1, ")
      fileTemp.write("0.1, ")
    else:
      for color in mainColor:
        fileTemp.write(str(color) + ", ")
    fileTemp.write("SPHERE, ")
    for coord in point[1:]:
      fileTemp.write(str(coord) + ", ")
    if radius:
      fileTemp.write(str(point[0]) + ", \n")
    else:
      fileTemp.write(str(gridSize/2.) + ", \n")
  fileTemp.write(" ]\n")
  fileTemp.write("cmd.load_cgo(surfObj,'" +filename + "')\n")
  fileTemp.close()

def outputGridClusters(grid, unionfindStructure, triGrid, actualSet=False):
  #writes a grid where everything is 0, the biggest set is +2, the ones with tris are -2
  #the ones outside the convex hull are marked as +4 or -4
  unionLists = unionfindStructure.toLists()
  while 1 > len(unionLists): #probably never happens but better than crashing
    unionLists.append([])
  unionLists.sort(lambda x,y: cmp(len(x), len(y)))
  #biggest 2 are last now
  bigSets = (set(), set())
  for index,bigList in enumerate(unionLists[-2:]):
    for value in bigList:
      bigSets[index].add(value[1:])
  if actualSet:
    bigSets = (bigSets[1], set())
    for gridPt in actualSet:
      bigSets[1].add(gridPt[1:])
  newGrid = []
  for indexX,rowX in enumerate(grid):
    newX = []
    for indexY,rowY in enumerate(rowX):
      newY = []
      for indexZ,entryZ in enumerate(rowY):
        value = 0
        if entryZ[1:] in bigSets[1]:
          value = 2
          if len(triGrid[indexX][indexY][indexZ][0]) > 0:
            value = 4
        elif entryZ[1:] in bigSets[0]:
          value = -2
          if len(triGrid[indexX][indexY][indexZ][0]) > 0:
            value = -4
        newEntry = value,entryZ[1],entryZ[2],entryZ[3]
        newY.append(newEntry)
      newX.append(newY)
    newGrid.append(newX)
  return newGrid

def pointDebug(pointList, radius=0.5, filename="temp.py", mainColor=(0.0,.01,.9)):
  fileTemp = open(filename, 'w')
  fileTemp.write("from pymol.cgo import * \nfrom pymol import cmd\nsurfObj = [ ")
  for point in pointList:
    fileTemp.write("COLOR, ")
    for color in mainColor:
      fileTemp.write(str(color) + ", ")
    fileTemp.write("SPHERE, ")
    for coord in point:
      fileTemp.write(str(coord) + ", ")
    fileTemp.write(str(radius) + ", \n")
  fileTemp.write(" ]\n")
  fileTemp.write("cmd.load_cgo(surfObj,'" +filename + "')\n")
  fileTemp.close()

colorGrads = {'flat':[1.0,1.0,1.0, 1.0,1.0,1.0, 1.0,1.0,1.0],
              'black':[0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0],
              'mono':[1.0,1.0,1.0, 0.5,0.5,0.5, 0.0,0.0,0.0],
              'gwg':[0.5,0.5,0.5, 1.0,1.0,1.0, 0.0,1.0,0.0],
              'rwb':[1.0,0.0,0.0, 1.0,1.0,1.0, 0.0,0.0,1.0],
              'bwr':[0.0,0.0,1.0, 1.0,1.0,1.0, 1.0,0.0,0.0],
              'rgb':[1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0],
              'bgr':[0.0,0.0,1.0, 0.0,1.0,0.0, 1.0,0.0,0.0],
              'octant':[1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0, 0.5,0.5,0.0,
                        0.5,0.0,0.5, 0.0,0.5,0.5, 1.0,1.0,1.0, 0.1,0.1,0.1]}

def decideColor(value, min, med, max, gradient,cycle=False):
  '''tstdisp can't be imported outside of pymol, so this function is copied.
  yes i realize this is bad. no real good way around it'''
  if cycle: #new case, make more like topo map...
    #ignore mediumValueColor
    #use the colors in the gradient as whole colors, switch based on mod value
    whichCol = int((value-min)/4) % (len(gradient)/3)
    red = gradient[int(whichCol*3+0)]
    gre = gradient[int(whichCol*3+1)]
    blu = gradient[int(whichCol*3+2)]
    return [red,gre,blu]
  elif len(gradient) == 9:
    if min > value:
      return [gradient[0],gradient[1],gradient[2]]
    elif min <= value and value < med:
      alpha = (value-min)/(med-min)
      red = (1.-alpha)*gradient[0] + (alpha)*gradient[3]
      gre = (1.-alpha)*gradient[1] + (alpha)*gradient[4]
      blu = (1.-alpha)*gradient[2] + (alpha)*gradient[5]
      return [red,gre,blu]
    elif med <= value and value <= max:
      alpha = (value-med)/(max-med)
      red = (1.-alpha)*gradient[3] + (alpha)*gradient[6]
      gre = (1.-alpha)*gradient[4] + (alpha)*gradient[7]
      blu = (1.-alpha)*gradient[5] + (alpha)*gradient[8]
      return [red,gre,blu]
    elif value > max:
      return [gradient[6],gradient[7],gradient[8]]
    else:
      print "something funny with colors"
      return [0.0,0.0,0.0]
  elif len(gradient) == 24:
    valueOct = int(value)-1
    #print valueOct
    red = gradient[valueOct*3]
    gre = gradient[valueOct*3+1]
    blu = gradient[valueOct*3+2]
    return [red,gre,blu]
