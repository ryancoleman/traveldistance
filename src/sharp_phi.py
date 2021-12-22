#!/usr/bin/env python
# ryan g. coleman, ryangc ATSYMBOL mail.med.upenn.edu
# kim sharp lab
# phi.py enables read/write of binary phi-maps from delphi
# usage is import, then pass in filename to object creation

import struct
import array
import sys
import string
import os
import math
#import gzip, bz2 #for compressed file reading (not enabled yet)

# format follows
#       character*20 toplabel
#       character*10 head,character*60 title
#       real*4 phi(65,65,65)
#       character*16 botlabel
#       real*4 scale, oldmid(3)

class phi(object):
  gridSizes = 193, 129, 65  # maybe should be parameter? #257 removed for now

  def __init__(self, phiFileName=False, is64=False):
    '''reads the phi file from disk'''
    self.oldmid = [0., 0., 0.]
    self.__minsmaxs = None
    self.__boundaries = None
    if phiFileName:  # otherwise just creating an empty phi map for writing
      for gridSize in phi.gridSizes:
        try:
          phiFile = open(phiFileName, 'rb')  # b is for binary, r is for read
          tempArray = array.array('f')
          junk = struct.unpack('4s', phiFile.read(4))
          (check,) = struct.unpack('4s', phiFile.read(4))
          #print "junk, check:", junk, check
          if check == " now":
            #print "32bit phimap"
            pass
          else:
            print "64bit phimap, will downconvert"
            is64 = True
          if not is64:
            (temptop,) = struct.unpack('16s', phiFile.read(16))
            self.toplabel = check + temptop
          else:
            (temptop,) = struct.unpack('20s', phiFile.read(20))
            self.toplabel = temptop
          #print "toplabel:", self.toplabel
          junk = struct.unpack('8s', phiFile.read(8))
          if is64:
            junk = struct.unpack('8s', phiFile.read(8))
          (self.head,) = struct.unpack('10s', phiFile.read(10))
          #print "head:", self.head
          (self.title,) = struct.unpack('60s', phiFile.read(60))
          #print "title:", self.title
          junk = struct.unpack('8s', phiFile.read(8))
          if is64:
            junk = struct.unpack('8s', phiFile.read(8))
          #next line raises error if grid too big
          #GxGxG -> packed into an array xyz order samplePhi = array.array('f')
          tempArray.fromfile(phiFile, gridSize**3)
          junk = struct.unpack('8s', phiFile.read(8))
          if is64:
            junk = struct.unpack('8s', phiFile.read(8))
          self.gridDimension = gridSize
          self.phiArray = tempArray
          break  # read successfully, just go on and read the last bits
        except EOFError:
          phiFile.close()
      (self.botlabel,) = struct.unpack('16s', phiFile.read(16))
      #print "botlabel:", self.botlabel
      junk = struct.unpack('8s', phiFile.read(8))
      if is64:
        junk = struct.unpack('8s', phiFile.read(8))
      (self.scale, self.oldmid[0], self.oldmid[1], self.oldmid[2],) = \
          struct.unpack('ffff', phiFile.read(16))
      #print "scale, oldmid:", self.scale, self.oldmid
      junk = struct.unpack('4s', phiFile.read(4))
      phiFile.close()

  def write(self, phiFileName=False):
    '''write data to member data structure manually,
    then call this to write to file
    the pad lines reproduce the binary padding of an original
    fortran formatted phi file'''
    if phiFileName:  # do nothing if no filename given
      phiFile = open(phiFileName, 'wb')  # b may be unnecessary, have to check
      phiFile.write(struct.pack('4b', 20, 0, 0, 0))  # pad
      phiFile.write(struct.pack('20s', self.toplabel))
      phiFile.write(struct.pack('8b', 20, 0, 0, 0, 70, 0, 0, 0))  # pad
      phiFile.write(struct.pack('10s', self.head))
      phiFile.write(struct.pack('60s', self.title))
      phiFile.write(struct.pack('8b', 70, 0, 0, 0, 4, -61, 16, 0))  # pad
      self.phiArray.tofile(phiFile)  # array
      phiFile.write(struct.pack('8b', 4, -61, 16, 0, 16, 0, 0, 0))  # pad
      phiFile.write(struct.pack('16s', self.botlabel))
      phiFile.write(struct.pack('8b', 16, 0, 0, 0, 16, 0, 0, 0))  # pad
      phiFile.write(
          struct.pack(
              'ffff', self.scale, self.oldmid[0], self.oldmid[1],
              self.oldmid[2]))
      phiFile.write(struct.pack('4b', 16, 0, 0, 0))  # pad
      phiFile.close()

  def getMinsMaxs(self):
    '''finds the positions of the extreme grid corners'''
    if self.__minsmaxs is None:
      mins, maxs = [], []
      for center in self.oldmid:
        mins.append(center-((self.gridDimension-1.)/(2.*self.scale)))
        maxs.append(center+((self.gridDimension-1.)/(2.*self.scale)))
      self.__minsmaxs = mins, maxs
    return self.__minsmaxs

  def getMinMaxValues(self):
    '''finds the minimum and maximum value'''
    return min(self.phiArray), max(self.phiArray)

  def getMaxValues(self):
    '''just the max'''
    return max(self.phiArray)

  def countValues(self):
    '''counts the occurence of each value'''
    counts = {}
    for value in self.phiArray:
      if value in counts:
        counts[value] += 1
      else:
        counts[value] = 1
    return counts

  def histogramValues(self, width=1.):
    '''makes a basic histogram'''
    ends = self.getMinMaxValues()
    bars = int(math.ceil((ends[1]-ends[0])/width)+1)
    counts = [0 for xVal in xrange(bars)]  # just make a list of 0s of right len
    for value in self.phiArray:
      counts[int(math.floor((value-ends[0])/width))] += 1
    return counts

  def getXYZ(self, xInd, yInd, zInd):
    '''returns the xyz coordinate of the center of the box'''
    mins, maxs = self.getMinsMaxs()
    gap = 1./self.scale
    return mins[0]+(xInd*gap), mins[1]+(yInd*gap), mins[2]+(zInd*gap)

  def getValue(self, xInd, yInd, zInd):
    '''for a given set of indices, return the value in the array'''
    index = int(zInd*(self.gridDimension**2.) + yInd*self.gridDimension + xInd)
    return self.phiArray[index]

  def setValue(self, xInd, yInd, zInd, value):
    '''puts the value into the phi array'''
    index = int(zInd*(self.gridDimension**2.) + yInd*self.gridDimension + xInd)
    self.phiArray[index] = value

  def transform(self, threshold=6.0, inside=-2.0, outside=-1.0):
    '''for every value in the array, change it to inside or outside,
    destructively overwrites old values'''
    for index in xrange(len(self.phiArray)):
      value = self.phiArray[index]
      if value < threshold:
        where = outside
      else:
        where = inside
      self.phiArray[index] = where

  def findBoundaries(
      self, inside=-2.0, border=2, pointXYZ=None, pointList=None):
    '''finds the extreme x, y, z positions that enclose all inside positions'''
    if self.__boundaries is None:  # need to calculate it
      if pointXYZ is not None:
        self.__boundaries = self.findPointMinsMaxs(pointXYZ, pointList)
      else:
        self.__boundaries = [
            self.gridDimension, self.gridDimension, self.gridDimension], [
            0, 0, 0]
      for x in xrange(self.gridDimension):
        for y in xrange(self.gridDimension):
          for z in xrange(self.gridDimension):
            if x < self.__boundaries[0][0] or x > self.__boundaries[1][0] or \
                y < self.__boundaries[0][1] or y > self.__boundaries[1][1] or \
                z < self.__boundaries[0][2] or z > self.__boundaries[1][2]:
              value = self.getValue(x, y, z)
              if value == inside:
                indices = (x, y, z)
                for coord in range(3):
                  self.__boundaries[0][coord] = min(
                      self.__boundaries[0][coord], indices[coord])
                  self.__boundaries[1][coord] = max(
                      self.__boundaries[1][coord], indices[coord])
      for coord in range(3):
        self.__boundaries[0][coord] = max(0, self.__boundaries[0][coord]-border)
        self.__boundaries[1][coord] = min(
            self.gridDimension, self.__boundaries[1][coord]+border)
    return self.__boundaries

  def getBoundaryLengths(self, inside=-2.0, border=2):
    '''calls findBoundaries if necessary, returns the lengths (max-min)'''
    if self.__boundaries is None:  # need to calculate it
      self.findBoundaries(inside, border)
    lengths = [
        self.__boundaries[1][0] - self.__boundaries[0][0],
        self.__boundaries[1][1] - self.__boundaries[0][1],
        self.__boundaries[1][2] - self.__boundaries[0][2]]
    return lengths

  def createFromGrid(
      self, grid, gridSize, defaultValue=0.0, toplabel="",
      head="", title="", botlabel="", lowestGridSize=65):
    '''does grid->phi data structure conversion'''
    self.toplabel = toplabel[:20]  # easy stuff first
    self.head = head[:10]
    self.title = title[:60]
    self.botlabel = botlabel[:16]
    lens = [len(grid), len(grid[0]), len(grid[0][0])]
    #have to expand to valid gridSize
    newGridSize = 0
    for possibleGridSize in self.gridSizes:
      good = True
      if possibleGridSize < lowestGridSize:
        good = False
      for oneLength in lens:
        if oneLength > possibleGridSize:
          good = False
      if good:
        newGridSize = possibleGridSize
    self.gridDimension = newGridSize
    #now take care of the grid
    self.phiArray = array.array('f')
    for z in xrange(self.gridDimension):
      for y in xrange(self.gridDimension):
        for x in xrange(self.gridDimension):
          if x < lens[0] and y < lens[1] and z < lens[2]:
            self.phiArray.append(grid[x][y][z][0])
          else:  # outside real grid
            self.phiArray.append(defaultValue)
    #scale and oldmid are all that is left
    self.scale = 1./gridSize
    for coord in range(3):
      self.oldmid[coord] = grid[0][0][0][coord+1] - \
          (gridSize/2.) + (self.gridDimension/self.scale)/2.
    #data should be ready for writing now

  def findPointMinsMaxs(self, pointXYZ, pointList):
    minsPts = pointXYZ[0][1:]
    maxsPts = pointXYZ[0][1:]
    for point in pointList:
      xyz = pointXYZ[point-1][1:]
      for coord in range(3):
        minsPts[coord] = min(minsPts[coord], xyz[coord])
        maxsPts[coord] = max(maxsPts[coord], xyz[coord])
    newMins = list(self.getIndices(minsPts))
    newMaxs = list(self.getIndices(maxsPts))  # so they initialize to pts
    return newMins, newMaxs

  def getIndices(self, pt):
    '''helper function to find the box a point is in'''
    mins, maxs = self.getMinsMaxs()
    gridSize = 1./self.scale
    xIndex = int(math.floor((pt[0]-mins[0])/gridSize))
    yIndex = int(math.floor((pt[1]-mins[1])/gridSize))
    zIndex = int(math.floor((pt[2]-mins[2])/gridSize))
    #print xIndex, yIndex, zIndex, mins, pt, maxs
    return xIndex, yIndex, zIndex

if -1 != string.find(sys.argv[0], "sharp_phi.py"):
  phiData = phi(sys.argv[1])
  print phiData.countValues()
  print phiData.scale
  print phiData.oldmid
