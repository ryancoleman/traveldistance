#!/usr/bin/env python2.5
#script to read .tst format and make surfaces in pymol
#ryan coleman

#from pymol.cgo import *
#from pymol import cmd
import math,  string
#for running commandline instead of through pymol
import sys, os
import tstdata
  
def tstConvexHull(execKHull, tstFileName, execQHull=False, runKHull=False):
  tstD = tstdata.tstData(tstFileName) #read the file into the data structure
  tempFile = open(tstFileName + ".tempQHullFile", 'w')
  tempFile.write("3\n") #dimension
  tempFile.write(str(len(tstD.dict['POINT_XYZ'])) + "\n") #number of points
  for points in tstD.dict['POINT_XYZ']:
    tempFile.write(str(points[1]) + " ") 
    tempFile.write(str(points[2]) + " ") 
    tempFile.write(str(points[3]) + "\n")
  tempFile.close()
  #call to qhull
  if runKHull:
    if os.path.exists(execKHull):
      os.popen("cat " + tempFile.name+ " | " + execKHull + " TO " + \
               tempFile.name+ ".output") #run khull
    else:
      print "khull does not exist or is not in the correct place:", execKHull
      exit(1)
  else:
    if os.path.exists(execQHull):
      os.popen("cat " + tempFile.name+ " | " + execQHull + " i n TO " + \
               tempFile.name+ ".output") #run qhull
    else:
      print "qhull does not exist or is not in the correct place:", execQHull
      exit(1)
  tempOutputFile = open(tempFile.name + ".output", 'r')
  totalFaces =tempOutputFile.readline() #first line contains number of faces
  numberFaces = 0 #actually triangles now
  convexHullFaces = [] #actually triangles now
  triNumToFaceNum = {} #dict related triangle index (1-indexed) to face index (0-indexed)
  pointsInTriList = []
  for count in xrange(len(tstD.dict['POINT_XYZ'])):
    pointsInTriList.append([0])
  for count in xrange(int(totalFaces)):
    line = tempOutputFile.readline()
    #read in the faces, usually triangles but can be any polygon
    tokens = string.split(line)
    #i bet these are 0-indexed, so make them 1-indexed
    faceIndices = [int(x) + 1 for x in tokens]
    for numberTri in xrange(len(faceIndices) - 2): #turn polygons into tris
      triIndices = [faceIndices[0]]
      for oneIndex in faceIndices[numberTri+1:numberTri+3]:
        triIndices.append(oneIndex)
      numberFaces += 1 #1-index the list
      triNumToFaceNum[numberFaces] = count
      triIndices.insert(0, numberFaces)
      convexHullFaces.append(triIndices)    
      for points in triIndices[1:]:
        pointsInTriList[points-1][0] += 1
        pointsInTriList[points-1].append(triIndices[0])
  dimension = tempOutputFile.readline() #dimension+1  of the normals (x, y, z, offset)
  totalFaces = tempOutputFile.readline() #repeat of # of faces
  outputLines = tempOutputFile.readlines()
  triangleNormals = []
  triangleNormalsCount = 0
  for triangle in convexHullFaces:
    faceNormalTokens = string.split(outputLines[triNumToFaceNum[int(triangle[0])]])
    faceNormal = [float(x) for x in faceNormalTokens]
    triangleNormalsCount += 1
    faceNormal.insert(0, triangleNormalsCount)
    triangleNormals.append(faceNormal)
  tempOutputFile.close()
  chNormals = []
  for index, pointsAndFaces in enumerate(pointsInTriList):
    normal = [0.0,0.0,0.0]
    if pointsAndFaces[0] > 0:
      normalTotal = [0.0,0.0,0.0]
      for faces in pointsAndFaces[1:]:
        for coord in range(3):
          normalTotal[coord] += triangleNormals[faces-1][coord+1]
      normal = [x/float(pointsAndFaces[0]) for x in normalTotal]
    normal.insert(0, index+1)  
    chNormals.append(normal)
  #add index data to pointsInTriList
  for index,pointInTri in enumerate(pointsInTriList):
    pointInTri.insert(0,index+1)
  #add to data
  tstD.dict['CONVEX_HULL_TRI_POINT_LIST'] = convexHullFaces
  tstD.dict['CONVEX_HULL_TRI_NORM_LIST'] = triangleNormals
  tstD.dict['NORM_XYZ_CONVEX_HULL'] = chNormals
  tstD.dict['CONVEX_HULL_POINT_TRI_LIST'] = pointsInTriList
  #add to file
  tstFile = open(tstFileName, 'a')
  tstFile.write("CONVEX_HULL_TRI_POINT_LIST\n")
  for line in convexHullFaces:
    lineOut = ""
    for count in xrange(len(line)):
      lineOut += "%8d " % line[count]
    tstFile.write(lineOut)
    tstFile.write("\n")
  tstFile.write("END CONVEX_HULL_TRI_POINT_LIST\n")
  tstFile.write("CONVEX_HULL_TRI_NORM_LIST\n")
  for line in triangleNormals:
    lineOut = "%8d" % line[0]
    for count in xrange(1,len(line)):
      lineOut += "%+15.8f " % line[count]
    noPlusLine = string.replace(lineOut, "+", " ")
    tstFile.write(noPlusLine)
    tstFile.write("\n")
  tstFile.write("END CONVEX_HULL_TRI_NORM_LIST\n")
  tstFile.write("NORM_XYZ_CONVEX_HULL\n")
  for line in chNormals:
    lineOut = "%8d" % line[0]
    for count in xrange(1,len(line)):
      lineOut += "%+15.8f " % line[count]
    noPlusLine = string.replace(lineOut, "+", " ")
    tstFile.write(noPlusLine)
    tstFile.write("\n")
  tstFile.write("END NORM_XYZ_CONVEX_HULL\n")
  tstFile.write("CONVEX_HULL_POINT_TRI_LIST\n")
  for line in pointsInTriList:
    lineOut = ""
    for count in xrange(len(line)):
      lineOut += "%8d " % line[count]
    tstFile.write(lineOut)
    tstFile.write("\n")
  tstFile.write("END CONVEX_HULL_POINT_TRI_LIST\n")
  tstFile.close()

#this is where main will go
if -1 != string.find(sys.argv[0], "tstConvexHull.py"): 
  if len(sys.argv) > 2: #else do nothing, read in as module
    for tstFile in sys.argv[2:]:
      print "running convex hull for", tstFile
      tstConvexHull(False, tstFile, sys.argv[1])
  else:
    print "Usage: tstConvexHull qhull-executable tstFile [list of additional tst files]"
