#module contains tstData class
#ryan g coleman ryangc@mail.med.upenn.edu

import string
import bz2
import geometry

class tstData(object):
  '''this class holds all the data from a tst data file in a dictionary
  it automatically reads in any standard formatted record
  (and ignores parentheses)'''

  #one example of necessaryKeys that can be used
  necessaryKeysForMesh = [
      'CONVEX_HULL_TRI_POINT_LIST', 'CONVEX_HULL_POINT_TRI_LIST', 'POINT_XYZ',
      'TRIANGLE_POINT', 'POINT_NEIGHBOR']
  necessaryKeysForPocket = [
      'CONVEX_HULL_TRI_POINT_LIST', 'CONVEX_HULL_POINT_TRI_LIST', 'POINT_XYZ',
      'TRIANGLE_POINT', 'POINT_NEIGHBOR', 'POINT_TRIANGLE', 'PDB_RECORD',
      'POINT_PDB_RECORD']
  necessaryKeysForSV = [
      'POINT_XYZ', 'TRIANGLE_POINT', 'NORM_XYZ', 'POINT_NEIGHBOR',
      'CONVEX_HULL_TRI_POINT_LIST']
  necessaryKeysForHoles = [
      'CONVEX_HULL_TRI_POINT_LIST', 'CONVEX_HULL_POINT_TRI_LIST', 'POINT_XYZ',
      'NORM_XYZ', 'TRIANGLE_POINT', 'POINT_TRIANGLE', 'PDB_RECORD',
      'POINT_PDB_RECORD', 'TRIANGLE_PDB_RECORD', 'POINT_NEIGHBOR']
  necessaryKeysForDisp = [
      'TRIANGLE_POINT', 'POINT_XYZ', 'NORM_XYZ', 'CONVEX_HULL_TRI_POINT_LIST',
      'NORM_XYZ_CONVEX_HULL', 'CURVATURE_XYZ', 'PROPERTY_XYZ', 'POTENTIAL_XYZ',
      'POINT_OCTANT', 'DEPTH_TRAVEL_DIST', 'DEPTH_EUCLID_DIST',
      'DEPTH_TRAVEL_DIST_EXTREMA', 'DEPTH_TRAVEL_DIST_TOPO', 'CHARGE_XYZ',
      'HYDROPHOBIC_XYZ', 'POINT_PDB_RECORD']
  necessaryKeysForCurve = [
      'POINT_XYZ', 'TRIANGLE_POINT', 'POINT_TRIANGLE', 'POINT_NEIGHBOR']

  def __init__(self, filename, necessaryKeys=False, drawCountTmp=0):
    '''necessaryKeys can be set to only read in certain records from tst'''
    self.dict = {}  # this dictionary will hold a key:value pair
    #where the key is the record name
    #and the value is the list (possibly multiply dimensioned)
    self.drawCount = drawCountTmp  # counts number of draws, names each surface
    #consecutively, so user can turn them on/off
    if filename:  # otherwise we don't do anything
      if filename.endswith("bz2"):
        bzTstFile = bz2.BZ2File(filename, 'r')
        records, recordNames = self.readFile(bzTstFile, necessaryKeys)
        bzTstFile.close()
      else:  # assume normal tst file
        tstFile = open(filename, 'r')
        records, recordNames = self.readFile(tstFile, necessaryKeys)
        tstFile.close()
      #puts data into dictionary for easy look-up
      for count in xrange(len(records)):
        self.dict[recordNames[count]] = records[count]
      self.__neighbors = False

  def readFile(self, tstFile, necessaryKeys=False):
    recordNames = []
    records = []
    curRecName = ''
    curRec = []
    recordTypes = {}  # maps curRecName to int, float, or string
    try:
      for line in tstFile:
        if curRecName == '':  # look for new one
          if len(line) > 0:
            tokens = string.split(line)
            for token in tokens:
              if len(token) > 0:
                curRecName = token
                break  # the first positive length token is the name
                if curRecName in recordTypes:  # old, so delete
                  del recordTypes[curRecName]
        else:  # we are inside a record, or almost done with one
          if string.find(line, "END", 0, 5) != -1:  # found the end
            if not necessaryKeys or curRecName in necessaryKeys:
              recordNames.append(curRecName)
              records.append(curRec)
            curRecName = ''
            curRec = []
          else:  # add this line to the current record
            if not necessaryKeys or curRecName in necessaryKeys:
              #try to avoid some issues with floats that are too long by
              #inserting spaces before minus '-' signs...
              newLine = string.replace(line, '-', ' -')
              tokens = string.split(newLine)
              if curRecName in recordTypes:
                if recordTypes[curRecName] == 'int':
                  ints = [int(x) for x in tokens]
                  if len(ints) == 1:
                    curRec.append(ints[0])
                  elif len(ints) > 1:
                    curRec.append(ints)
                elif recordTypes[curRecName] == 'float':
                  floats = [float(x) for x in tokens]
                  if len(floats) == 1:
                    curRec.append(floats[0])
                  elif len(floats) > 1:
                    curRec.append(floats)
                elif recordTypes[curRecName] == 'string':
                  curRec.append(line[:-1])  # remove then
              else:
                try:   # first try to convert to ints
                  ints = [int(x) for x in tokens]
                  if len(ints) == 1:
                    curRec.append(ints[0])
                  elif len(ints) > 1:
                    curRec.append(ints)
                  recordTypes[curRecName] = 'int'
                except ValueError:
                  try:   # first try to convert to floats
                    floats = [float(x) for x in tokens]
                    if len(floats) == 1:
                      curRec.append(floats[0])
                    elif len(floats) > 1:
                      curRec.append(floats)
                    recordTypes[curRecName] = 'float'
                  except ValueError:  # if completely failed, use strings
                    curRec.append(line[:-1])  # remove then
                    recordTypes[curRecName] = 'string'
    except StopIteration:
      pass  # read EOF,okay to quit, if file is bad might lose last record
    return records, recordNames

  def eulerCharacteristic(self):
    vertices = len(self.dict['POINT_PDB_RECORD'])
    faces = len(self.dict['TRIANGLE_PDB_RECORD'])
    edges = 0
    for pointTriList in self.dict['POINT_TRIANGLE']:
      edges += pointTriList[1]
    edges = edges / 2
    euler = vertices - edges + faces
    #print vertices, edges, faces, euler
    return euler

  def countHandles(self):
    '''
    assumes single piece
    '''
    euler = self.eulerCharacteristic()
    handles = (euler-2)/-2
    return handles

  def buildNeighbors(self, allPoints):
    '''builds a neighbor dictionary'''
    if (not self.__neighbors) or 0 == len(self.__neighbors):  # remake
      self.__neighbors = {}
      for pointStart in allPoints:
        neighborList = self.dict['POINT_NEIGHBOR'][pointStart-1]
        startXYZ = self.dict['POINT_XYZ'][pointStart-1][1:]
        tempList = []
        for neighborPoint in neighborList[2:]:  # first 2 are p#, order
          endXYZ = self.dict['POINT_XYZ'][neighborPoint-1][1:]
          distance = geometry.distL2(startXYZ, endXYZ)
          tempList.append([neighborPoint, distance])
        self.__neighbors[pointStart] = tempList
    return self.__neighbors

def writeEntryIntegers(data, title, endtitle, file):
  '''
  file should be open already, data is the data,
  title and endtitle are the headers
  '''
  file.write(title + "\n")
  for line in data:
    lineOut = ""
    for entry in line:
      lineOut += "%8d" % entry
    file.write(lineOut + "\n")
  file.write(endtitle + "\n")

def writeEntryFloats(data, title, endtitle, file):
  '''
  file should be open already, data is the data,
  title and endtitle are the headers
  '''
  file.write(title + "\n")
  for line in data:
    lineOut = "%8d" % line[0]
    for entry in line[1:]:
      lineOut += "%+9.4f" % entry
    noPlusLine = string.replace(lineOut, "+", " ")
    file.write(noPlusLine + "\n")
  file.write(endtitle + "\n")

def writeEntrySingleFloat(data, title, endtitle, file):
  '''
  file should be open already, data is the data,
  title and endtitle are the headers
  '''
  file.write(title + "\n")
  for line in data:
    lineOut = "%8d" % line[0]
    for entry in line[1:]:
      lineOut += "%+15.6f" % entry
    noPlusLine = string.replace(lineOut, "+", " ")
    file.write(noPlusLine + "\n")
  file.write(endtitle + "\n")

def writeEntryString(data, title, endtitle, file):
  file.write(title + "\n")
  for line in data:
    if line.endswith("\n"):
      file.write(line)
    else:
      file.write(line + "\n")
  file.write(endtitle + "\n")

def trianglinizeLoop(loopPoints):
  tris = []
  if len(loopPoints) >= 2:  # otherwise can't do anything
    startPt = loopPoints[0]
    lastPt = loopPoints[1]
    for nextPt in loopPoints[2:]:
      tris.append([startPt, lastPt, nextPt])
      lastPt = nextPt
  return tris

class tstDataWritable(tstData):

  def __init__(self, filename, necessaryKeys=False, drawCountTmp=0):
    '''necessaryKeys can be set to only read in certain records from tst'''
    self.dict = {}  # this dictionary will hold a key:value pair
    #where the key is the record name
    #and the value is the list (possibly multiply dimensioned)
    self.drawCount = drawCountTmp  # counts number of draws, names each surface
    #consecutively, so user can turn them on/off
    if necessaryKeys:  # need to add ones we need to write the data back out
      necessaryKeys.extend([
          'TRIANGLE_NEIGHBOR', 'POINT_XYZ', 'TRIANGLE_POINT',
          'POINT_TRIANGLE', 'POINT_NEIGHBOR', 'NORM_XYZ', 'POINT_PDB_RECORD',
          'TRIANGLE_PDB_RECORD', 'PDB_RECORD', 'CURVATURE_XYZ',
          'PROPERTY_XYZ'])
    if filename:
      if filename.endswith("bz2"):
        bzTstFile = bz2.BZ2File(filename, 'r')
        records, recordNames = self.readFile(bzTstFile, necessaryKeys)
        bzTstFile.close()
      else:  # assume normal tst file
        tstFile = open(filename, 'r')
        records, recordNames = self.readFile(tstFile, necessaryKeys)
        tstFile.close()
      #puts data into dictionary for easy look-up
      for count in xrange(len(records)):
        self.dict[recordNames[count]] = records[count]
      self.__neighbors = False

  def write(self, filename):
    '''
    writes completely new file out
    this code *should* be compatible
    needs to be updated if necessaryKeys was used to read file
    '''
    tstFile = open(filename, 'w')
    writeEntryIntegers(
        self.dict['TRIANGLE_NEIGHBOR'],
        "TRIANGLE_NEIGHBOR LIST (CLOCKWISE)",
        "END TRIANGLE_NEIGHBOR LIST", tstFile)
    writeEntryFloats(
        self.dict['POINT_XYZ'],
        "POINT_XYZ LIST",
        "END POINT_XYZ LIST", tstFile)
    writeEntryIntegers(
        self.dict['TRIANGLE_POINT'],
        "TRIANGLE_POINT LIST (CLOCKWISE really, this was checked)",
        "END TRIANGLE_POINT LIST", tstFile)
    writeEntryIntegers(
        self.dict['POINT_TRIANGLE'],
        "POINT_TRIANGLE LIST (CLOCKWISE)",
        "END POINT_TRIANGLE LIST", tstFile)
    writeEntryIntegers(
        self.dict['POINT_NEIGHBOR'],
        "POINT_NEIGHBOR LIST (CLOCKWISE)",
        "END POINT_NEIGHBOR LIST (CLOCKWISE)", tstFile)
    writeEntryFloats(
        self.dict['NORM_XYZ'],
        "NORM_XYZ LIST",
        "END NORM_XYZ LIST", tstFile)
    writeEntryIntegers(
        self.dict['POINT_PDB_RECORD'],
        "POINT_PDB_RECORD",
        "END POINT_PDB_RECORD", tstFile)
    writeEntryIntegers(
        self.dict['TRIANGLE_PDB_RECORD'],
        "TRIANGLE_PDB_RECORD",
        "END TRIANGLE_PDB_RECORD", tstFile)
    writeEntryString(
        self.dict['PDB_RECORD'],
        "PDB_RECORD", "END PDB_RECORD", tstFile)
    writeEntrySingleFloat(
        self.dict['CURVATURE_XYZ'], "CURVATURE_XYZ",
        "END CURVATURE_XYZ", tstFile)
    writeEntrySingleFloat(
        self.dict['PROPERTY_XYZ'], "PROPERTY_XYZ", "END PROPERTY_XYZ", tstFile)
    tstFile.close()
