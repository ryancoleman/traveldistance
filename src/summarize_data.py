#!/usr/bin/env python

#ryan g. coleman ryan.g.coleman ATSYMBOL gmail.com ryangc ATSYMBOL mail.med.upenn.edu
#kim sharp lab http://crystal.med.upenn.edu

import string
import sys
import statistics
import analyze_packing

if -1 != string.find(sys.argv[0], "summarize_data.py"):
  try:
    summaryData = analyze_packing.readSummaryFile(sys.argv[1])
    listPdbs = [
        analyze_packing.readList(fileName) for fileName in sys.argv[2:4]]
    for count in xrange(len(listPdbs[0])):
      for lindex in range(2):  # put both lists on same line
        print listPdbs[lindex][count], "\t",
        for name in summaryData.keys():
          print summaryData[name][listPdbs[lindex][count]], "\t",
      print " "  # end of line
  except IndexError:
    print "summarize_data.py summary.txt file1.txt file2.txt"
    print "uses information in summary.txt to compare matching files in 1 and 2"
    sys.exit(1)
