#script to read .tst format and make surfaces in pymol
#ryan g. coleman ryangc ATSYMBOL mail.med.upenn.edu

from pymol.cgo import *
from pymol import cmd
import math
import string

#this class holds all the data from a tst data file in a dictionary
#it automatically reads in any standard formatted record
#(and ignores parentheses)
class tstData(object):
  #this holds any key from tst files necessary for display. if you add a new
  #key, add it here so it can be read in. done to save memory
  #by not reading keys that are never used.
  necessaryKeysForDisp = [
      'TRIANGLE_POINT', 'POINT_XYZ', 'NORM_XYZ', 'CONVEX_HULL_TRI_POINT_LIST',
      'NORM_XYZ_CONVEX_HULL', 'CURVATURE_XYZ', 'PROPERTY_XYZ', 'POTENTIAL_XYZ',
      'POINT_OCTANT', 'DEPTH_TRAVEL_DIST', 'DEPTH_EUCLID_DIST',
      'DEPTH_TRAVEL_DIST_EXTREMA', 'DEPTH_TRAVEL_DIST_TOPO', 'CHARGE_XYZ',
      'HYDROPHOBIC_XYZ', 'POINT_PDB_RECORD', 'PDB_RECORD',
      'POINT_CURVATURE_EDGE', 'POINT_LEAF', 'LEAF_GROUP', 'POINT_NEARBY_ATOM']

  def __init__(self, filename, necessaryKeys=necessaryKeysForDisp):
    if filename:
      tstFile = open(filename, 'r')
      self.dict = {}
      records, recordNames = self.readFile(tstFile, necessaryKeys)
      tstFile.close()
      #puts data into dictionary for easy look-up
      for count in range(len(records)):
        self.dict[recordNames[count]] = records[count]

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
                  curRec.append(line[:-1])  # remove the \n
              else:
                try:  # first try to convert to ints
                  ints = [int(x) for x in tokens]
                  if len(ints) == 1:
                    curRec.append(ints[0])
                  elif len(ints) > 1:
                    curRec.append(ints)
                  recordTypes[curRecName] = 'int'
                except ValueError:
                  try:  # first try to convert to floats
                    floats = [float(count) for count in tokens]
                    if len(floats) == 1:
                      curRec.append(floats[0])
                    elif len(floats) > 1:
                      curRec.append(floats)
                    recordTypes[curRecName] = 'float'
                  except ValueError:  # if completely failed, use strings
                    curRec.append(line[:-1])  # remove the \n
                    recordTypes[curRecName] = 'string'
    except StopIteration:
      pass  # read EOF, okay to quit, if file is bad might lose last record
    return records, recordNames

class displayObj(object):
  '''contains all methods and local variables used to display things'''
  tstDatas = []  # global arrays basically
  tstNames = []
  pointGroupCache = [False]
  drawCount = 0
  lastGroup = None

  def tstOpen(self, filename):
    '''
    this function is the opening to pymol to open a tst file, read in the data
    overwrites old data, only allow 1 file at a time--perhaps change later
    '''
    if len(self.tstDatas) > 0:
      self.tstDatas[0] = tstData(filename)
      self.tstNames[0] = filename[:4]
    else:
      self.tstDatas.append(tstData(filename))
      self.tstNames.append(filename[:4])
    self.pointGroupCache = [False]

  def tstMultiOpen(self, filenames):
    '''uses arguments to open all the files passed in'''
    tstDatasLocal, tstNamesLocal = [], []  # reset no matter what
    for filename in filenames:
      tstDatasLocal.append(tstData(filename))
      tstNamesLocal.append(filename[:4])
    pointGroupCacheLocal = [False for count in range(len(filenames))]
    #diff for each
    self.tstDatas = tstDatasLocal
    self.tstNames = tstNamesLocal
    self.pointGroupCache = pointGroupCacheLocal
    #actually open and plan to use all tstDatas files

  #follows is all color stuff
  def decideColor(self, value, min, med, max, gradient, cycle=False):
    #print value, min, med, max, gradient  # useful debugging routine
    value = value[0]
    if cycle:  # new case, make more like topo map...
      #ignore mediumValueColor
      #use the colors in the gradient as whole colors,
      #switch based on mod value
      whichCol = int((value-min)/4) % (len(gradient)/3)
      red = gradient[int(whichCol*3 + 0)]
      gre = gradient[int(whichCol*3 + 1)]
      blu = gradient[int(whichCol*3 + 2)]
      return [red, gre, blu]
    elif len(gradient) == 9:
      if gradient[0] == gradient[3] and gradient[3] == gradient[6] and \
          gradient[1] == gradient[4] and gradient[4] == gradient[7] and \
          gradient[2] == gradient[5] and gradient[5] == gradient[8]:
        #all the same
        return [gradient[0], gradient[1], gradient[2]]  # just output first
      if min > value:
        return [gradient[0], gradient[1], gradient[2]]
      elif min <= value and value < med:
        alpha = (value-min)/(med-min)
        red = (1.-alpha)*gradient[0] + (alpha)*gradient[3]
        gre = (1.-alpha)*gradient[1] + (alpha)*gradient[4]
        blu = (1.-alpha)*gradient[2] + (alpha)*gradient[5]
        return [red, gre, blu]
      elif med <= value and value <= max:
        alpha = (value-med)/(max-med)
        red = (1.-alpha)*gradient[3] + (alpha)*gradient[6]
        gre = (1.-alpha)*gradient[4] + (alpha)*gradient[7]
        blu = (1.-alpha)*gradient[5] + (alpha)*gradient[8]
        return [red, gre, blu]
      elif value > max:
        return [gradient[6], gradient[7], gradient[8]]
      else:
        print "something funny with colors"
        return [0.0, 0.0, 0.0]
    elif len(gradient) == 24:
      valueOct = int(value) - 1
      #print valueOct
      red = gradient[valueOct * 3]
      gre = gradient[valueOct * 3 + 1]
      blu = gradient[valueOct * 3 + 2]
      return [red, gre, blu]

  colorModes = {
      'octant': ['octant', 'octant'],
      'curve': ['curve', 'gwg'],
      'curve2': ['curve2', 'gwg'],
      'elec': ['pot', 'rwb'],
      'prop': ['prop', 'rwb'],
      'red': ['prop', 'red'],
      'green': ['prop', 'green'],
      'blue': ['prop', 'blue'],
      'magenta': ['prop', 'magenta'],
      'yellow': ['prop', 'yellow'],
      'cyan': ['prop', 'cyan'],
      'flat': ['prop', 'flat'],
      'depth': ['depth', 'rgb'],
      'depthmono': ['depth', 'mono'],
      'depthmonorev': ['depth', 'monorev'],
      'charge': ['charge', 'rwb', 0.],
      'hydro': ['hydro', 'rwb', 0.],
      'euclid': ['euclid', 'rgb']}
  surfaces = {
      'ms': ['TRIANGLE_POINT', 'POINT_XYZ', 'NORM_XYZ'],
      'ch': ['CONVEX_HULL_TRI_POINT_LIST', 'POINT_XYZ', 'NORM_XYZ_CONVEX_HULL']}
  colorVals = {
      'curve': 'CURVATURE_XYZ',
      'curve2': 'POINT_CURVATURE_EDGE',
      'prop': 'PROPERTY_XYZ',
      'pot': 'POTENTIAL_XYZ',
      'octant': 'POINT_OCTANT',
      'depth': 'DEPTH_TRAVEL_DIST',
      'charge': 'CHARGE_XYZ',
      'hydro': 'HYDROPHOBIC_XYZ',
      'euclid': 'DEPTH_EUCLID_DIST',
      'depth3ext': 'DEPTH_TRAVEL_DIST_EXTREMA',
      'depth3topo': 'DEPTH_TRAVEL_DIST_TOPO'}
  colorGrads = {
      'flat': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
      'black': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
      'mono': [1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0],
      'monorev': [0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0],
      'gwg': [0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0],
      'rwb': [1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0],
      'bwr': [0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0],
      'rgb': [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
      'gbw': [0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
      'land': [0.0, 1.0, 0.0, 0.54, 0.27, 0.07, 1.0, 1.0, 1.0],
      'bgr': [0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0],
      'red': [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
      'green': [0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0],
      'blue': [0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0],
      'yellow': [1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0],
      'magenta': [1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0],
      'cyan': [0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0],
      'octant': [
          1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.5, 0.5, 0.0,
          0.5, 0.0, 0.5, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0, 0.1, 0.1, 0.1]}
  surfTypes = ['tris', 'edges', 'pts', 'pdb']
  pointAtomMap = 'POINT_PDB_RECORD'

  def printBadChoice(self):
    print "Choose valid color magic mode ---"
    print "Valid mode values are ",
    for val in self.colorModes.keys():
      print val + " ",
    print "\ntstDraw mode=blah"
    print "\nOR\nChoose valid color value and gradient---"
    print "Valid color values are ",
    for val in self.colorVals.keys():
      print val + " ",
    print "\nValid color gradients are ",
    for gra in self.colorGrads.keys():
      print gra + " ",
    print "\ntstDraw color=blah, grad=blah"
    print "Valid surface types are ",
    for su in self.surfTypes:
      print su + " ",
    print "\nValid surfaces to draw are ",
    for su in self.surfaces:
      print su + " ",
    print "\ntstDraw ... surf=blah which=blah"

  def displayPdb(self, tstD, tstName, group=False, pointGroupCacheTemp=None):
    '''creates an actual pdb object of the protein or just the part near the
    surface/group of interest'''
    if not group:  # just display the whole thing
      pdbText = "\n".join(tstD.dict['PDB_RECORD'])
      cmd.read_pdbstr(pdbText, tstName + "." + str(self.drawCount) + ".pdb")
    else:
      pointNearbyAtomList = tstD.dict['POINT_NEARBY_ATOM']
      import sets
      atomSet = sets.Set()  # old since pymol is based on 2.3 forever
      for pointNearbyAtom in pointNearbyAtomList:
        try:
          point = pointNearbyAtom[0]
          groups = pointGroupCacheTemp[point]
          if group in groups:
            atomSet.update(pointNearbyAtom[1:])
        except TypeError:
          pass  # means a point had no nearby atoms, which can happen given dist
        #print point, groups, atomSet
      pdbStr = ""
      for atomLine in atomSet:
        pdbStr += tstD.dict['PDB_RECORD'][int(atomLine)-1] + "\n"
      cmd.read_pdbstr(
          pdbStr, tstName + "." + str(self.drawCount) + ".pdb." + str(group))

  def tstMultiDraw(
      self, mode='', color='depth', grad='rgb', surf='tris', which='ms',
      minval='', medval='', maxval='', cycle=False, masklow='', maskhigh='',
      groups=None):
    '''calls tstdraw multiple times. groups should be a list of equal length
    to tstDatas and have None wherever you don't want to draw that one.
    or don't pass in groups and everything gets draw fully'''
    if groups is None:
      groups = []
    else:
      newG = []
      for token in string.split(groups[1:-1], ","):
        try:
          newG.append(int(token))
        except ValueError:
          newG.append(False)
      groups = newG
    for whichDraw in range(len(self.tstDatas)):
      if len(groups) == len(self.tstDatas):
        #draw groups differently, otherwise ignore
        if groups[whichDraw]:
          self.tstDraw(
              mode=mode, color=color, grad=grad, surf=surf, which=which,
              minval=minval, medval=medval, maxval=maxval, cycle=cycle,
              masklow=masklow, maskhigh=maskhigh, group=groups[whichDraw],
              drawWhichData=whichDraw)
      else:
        self.tstDraw(
            mode=mode, color=color, grad=grad, surf=surf, which=which,
            minval=minval, medval=medval, maxval=maxval, cycle=cycle,
            masklow=masklow, maskhigh=maskhigh, drawWhichData=whichDraw)

  def tstDrawParent(self, increment=1, drawWhichData=0):
    '''using self.lastGroup, finds the parent of that node and draws that'''
    increment = int(increment)
    if len(self.tstDatas) == 0:
      print "use tstOpen ----.tst before tstDraw\n"
      return 1
    else:
      tstD = self.tstDatas[drawWhichData]
    if 'LEAF_GROUP' not in tstD.dict:
      print "tst file does not have group/pocket data in it. run pocketmap"
      return 1
    if self.lastGroup is None:
      print "no previous group, use tstDrawDeepest or tstDraw group=id"
      return 1
    lastGroupInt = int(self.lastGroup)
    ancestors = tstD.dict['LEAF_GROUP'][lastGroupInt - 1]  # off by 1
    if increment < (len(ancestors)-2):  # actually an ancestor to show
      parent = ancestors[-1-increment]  # -1 is the group, -2 is the parent
    else:  # just return root, went too far up
      parent = ancestors[1]
    print "tstDraw group=" + str(parent)
    self.tstDraw(group=parent, drawWhichData=drawWhichData)

  def tstDrawDeepest(self, drawWhichData=0):
    '''finds the deepest leaf (lowest numbered leaf) in the tst file, should be
    leaf 1 but don't want to assume'''
    if len(self.tstDatas) == 0:
      print "use tstOpen ----.tst before tstDraw\n"
      return 1
    else:
      tstD = self.tstDatas[drawWhichData]
    deepest = 1e10000
    if 'LEAF_GROUP' not in tstD.dict:
      print "tst file does not have group/pocket data in it. run pocketmap"
      return 1
    for leafGroup in tstD.dict['LEAF_GROUP']:
      if leafGroup[0] < deepest:
        deepest = leafGroup[0]
    print "tstDraw group=" + str(deepest)
    self.tstDraw(group=deepest, drawWhichData=drawWhichData)

  def tstDraw(
      self, mode='', color='depth', grad='rgb', surf='tris', which='ms',
      minval='', medval='', maxval='', cycle=False, masklow='', maskhigh='',
      group=None, drawWhichData=0):
    '''this function allows user to draw the surface'''
    if len(self.tstDatas) == 0:
      print "use tstOpen ----.tst before tstDraw\n"
      return 1
    else:
      tstD = self.tstDatas[drawWhichData]
      tstName = self.tstNames[drawWhichData]
    if len(mode) > 0:
      if mode in self.colorModes.keys():
        color = self.colorModes[mode][0]
        grad = self.colorModes[mode][1]
        if len(self.colorModes[mode]) > 2:
          medval = self.colorModes[mode][2]
      else:
        self.printBadChoice()
        return 1
    if surf not in self.surfTypes:
      self.printBadChoice()
      return 1
    if group is not None:
      group = int(group)
      if not self.pointGroupCache[drawWhichData]:
        self.pointGroupCache[drawWhichData] = {}  # make a cache
        for pointLeaf in tstD.dict['POINT_LEAF']:
          if len(pointLeaf) > 1:
            point = pointLeaf[0]
            groupsTemp = tstD.dict['LEAF_GROUP'][pointLeaf[1] - 1]
            try:
              self.pointGroupCache[drawWhichData][point] = groupsTemp[1:]
            except TypeError:
              self.pointGroupCache[drawWhichData][point] = False
              #XXX problem here fix sometime
    if surf == 'pdb':
      if group is not None:
        self.displayPdb(
            tstD, tstName, group, self.pointGroupCache[drawWhichData])
      else:
        self.displayPdb(tstD, tstName)
      self.drawCount = self.drawCount + 1  # increment unique drawCount number
      return 1  # quit now, don't do the rest of this stuff
    surfObj = []
    if surf == 'tris':
      surfObj = [BEGIN, TRIANGLES]
    elif surf == 'edges':
      surfObj = [BEGIN, LINES]
    elif surf == 'pts':
      surfObj = [BEGIN, POINTS]
    donePts = []  # only used for points
    #fetch the data we need
    if which not in self.surfaces.keys():
      self.printBadChoice()
      return 1
    triPoints = tstD.dict[self.surfaces[which][0]]
    #so we can get convex hull, etc.
    pointXYZ = tstD.dict[self.surfaces[which][1]]
    pointNorm = tstD.dict[self.surfaces[which][2]]
    #hard, variable part is color setup
    localColorVals = self.colorVals.copy()
    if color not in self.colorVals.keys() and color in tstD.dict.keys():
      for newColor in tstD.dict.keys():
        localColorVals[newColor] = newColor
    if grad not in self.colorGrads.keys() or color not in localColorVals.keys():
      self.printBadChoice()
      return 1
    pointColorValue = tstD.dict[localColorVals[color]]
    #all of these should look up a single value
    #if the sizes don't match, assume mapping atom property to surface....
    if len(pointColorValue) != len(pointXYZ):
      newPCV = []
      pointAtom = tstD.dict[self.pointAtomMap]  # changes between the 2
      for pointIndex in pointXYZ:
        atomRef = pointAtom[pointIndex[0] - 1][1]  # 1-indexed
        atomPCV = pointColorValue[atomRef - 1][1]  # also 1-indexed
        newPCV.append([pointIndex[0], atomPCV])
      pointColorValue = newPCV  # now with atom property mapped to surface
    #allow full user control of min/med/max values if they want it
    #otherwise figure them out/guess
    minColorValue = 0.0
    if minval == '':
      #just a note, I'm not sure why these are using reduce & lambda
      #this code runs in pymol which is using a very old python
      #perhaps it had (or I thought it had) to do this rather than something
      #simpler.
      minColorValue = reduce(
          lambda xVal, yVal: min(xVal, yVal),
          [xVal[1] for xVal in pointColorValue])
    else:
      minColorValue = float(minval)
    maxColorvalue = 0.0
    if maxval == '':
      maxColorValue = reduce(
          lambda xVal, yVal: max(xVal, yVal),
          [xVal[1] for xVal in pointColorValue])
    else:
      maxColorValue = float(maxval)
    mediumColorValue = 0.0
    if medval == '':
      mediumColorValue = (minColorValue + maxColorValue) / 2.0
    else:
      mediumColorValue = float(medval)
    if cycle:  # is True, want to make something more like a topo map
      minColorValue = int(math.floor(minColorValue))
      maxColorValue = int(math.ceil(maxColorValue))
      mediumColorValue = int(mediumColorValue)
    if grad != 'octant':
      for outVal in [minColorValue, mediumColorValue, maxColorValue]:
        print str(outVal) + " maps to " + str(self.decideColor(
            [outVal], minColorValue, mediumColorValue,
            maxColorValue, self.colorGrads[grad], cycle))
    if masklow == '':
      masklow = reduce(
          lambda xVal, yVal: min(xVal, yVal),
          [xVal[1] for xVal in pointColorValue])
    masklow = float(masklow)
    if maskhigh == '':
      maskhigh = reduce(
          lambda xVal, yVal: max(xVal, yVal),
          [xVal[1] for xVal in pointColorValue])
    maskhigh = float(maskhigh)
    for oneTriPoints in triPoints:
      if surf == 'tris' or surf == 'pts':
        withinMask = True
        for point in oneTriPoints[1:]:
          for onePointColorVal in pointColorValue[point-1][1:]:
            if float(onePointColorVal) < masklow or \
                float(onePointColorVal) > maskhigh:
              withinMask = False
        if withinMask and group is not None:
          #needs checked more, now for group
          for point in oneTriPoints[1:]:
            groupsTemp = self.pointGroupCache[drawWhichData][point]
            if not (group in groupsTemp):
              withinMask = False
        if withinMask:
          for point in oneTriPoints[1:]:
            #all but the first, which is the tri#
            if surf == 'tris' or point not in donePts:
              if surf == 'pts':
                donePts.append(point)
              surfObj.append(COLOR)
              colorRGB = [1.0, 1.0, 1.0]
              if not color == 'flat':
                colorRGB = self.decideColor(
                    pointColorValue[point - 1][1:],
                    minColorValue, mediumColorValue, maxColorValue,
                    self.colorGrads[grad], cycle)
              for col in colorRGB:
                surfObj.append(col)
              surfObj.append(NORMAL)
              for norm in pointNorm[point-1][1:]:
                surfObj.append(norm)
              surfObj.append(VERTEX)
              for coordinate in pointXYZ[point-1][1:]:
                surfObj.append(coordinate)
      elif surf == 'edges':
        points = oneTriPoints[1:]
        pointsSh = [points[1], points[2], points[0]]
        for one, two in zip(points, pointsSh):
          if one < two:
            for point in [one, two]:
              surfObj.append(COLOR)
              colorRGB = self.decideColor(
                  pointColorValue[point-1][1:], minColorValue,
                  mediumColorValue, maxColorValue, self.colorGrads[grad], cycle)
              for col in colorRGB:
                surfObj.append(col)
              surfObj.append(NORMAL)
              for norm in pointNorm[point-1][1:]:
                surfObj.append(norm)
              surfObj.append(VERTEX)
              for coordinate in pointXYZ[point-1][1:]:
                surfObj.append(coordinate)
    surfObj.append(END)
    if group is None:
      cmd.load_cgo(surfObj, str(tstName) + '.s.' + str(self.drawCount))
    else:
      cmd.load_cgo(
          surfObj, str(tstName) + '.' + str(group) + '.s.' +
          str(self.drawCount))
      self.lastGroup = group
    self.drawCount = self.drawCount + 1  # increment unique drawCount number

  #new function to draw lots of surfaces at once for multiple files
  #only draws ten, use offset to get more
  def tstOpenDraw(
      self, fileMask, mode='', offset=0, minval='', medval='', maxval='',
      group=None):
    minvalIn, medvalIn, maxvalIn = minval, medval, maxval
    offset = int(offset)
    print fileMask
    from glob import glob
    files = glob(fileMask)
    for filename in files[offset:offset+10]:
      self.tstOpen(filename)
      self.tstDraw(
          mode, minval=minvalIn, medval=medvalIn, maxval=maxvalIn, group=group)

#this is really really dumb but i can't get it to work otherwise
displayInstance = displayObj()

def tstOpen(filename):
  displayInstance.tstOpen(filename)

def tstMultiOpen(*filenames):
  displayInstance.tstMultiOpen(filenames)

def tstOpenDraw(
    fileMask, mode='', offset=0, minval='', medval='', maxval='', group=None):
  displayInstance.tstOpenDraw(
      fileMask=fileMask, mode=mode, offset=offset,
      minval=minval, medval=medval, maxval=maxval, group=group)

def tstDraw(
    mode='', color='depth', grad='rgb', surf='tris', which='ms', minval='',
    medval='', maxval='', cycle=False, masklow='', maskhigh='', group=None):
  displayInstance.tstDraw(
      mode=mode, color=color, grad=grad, surf=surf, which=which,
      minval=minval, medval=medval, maxval=maxval, cycle=cycle,
      masklow=masklow, maskhigh=maskhigh, group=group)

def tstDrawDeepest():
  displayInstance.tstDrawDeepest()

def tstDrawParent(increment=1):
  displayInstance.tstDrawParent(increment=increment)

def tstMultiDraw(
    mode='', color='depth', grad='rgb', surf='tris', which='ms', minval='',
    medval='', maxval='', cycle=False, masklow='', maskhigh='', groups=None):
  displayInstance.tstMultiDraw(
      mode=mode, color=color, grad=grad, surf=surf, which=which,
      minval=minval, medval=medval, maxval=maxval, cycle=cycle,
      masklow=masklow, maskhigh=maskhigh, groups=groups)

cmd.extend("tstOpen", tstOpen)
cmd.extend("tstDraw", tstDraw)
cmd.extend("tstDrawDeepest", tstDrawDeepest)
cmd.extend("tstDrawParent", tstDrawParent)
cmd.extend("tstOpenDraw", tstOpenDraw)
cmd.extend("tstMultiOpen", tstMultiOpen)
cmd.extend("tstMultiDraw", tstMultiDraw)
