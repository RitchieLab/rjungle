#!/usr/bin/env python

import sys

if len(sys.argv) < 2:
	print "usage: %s <rjpedfile>" % sys.argv[0]
	sys.exit(1)
rjpedPath = sys.argv[1]

print "reading MAF files ..."
markerChr = {}
markerPos = {}
markerDupe = {}
for c in range(1,23):
	with open('chr%d.maf' % c) as mafFile:
		for line in mafFile:
			words = line.split()
			if words[1] in markerDupe:
				markerDupe[words[1]] += 1
			elif words[1] in markerChr:
				del markerChr[words[1]]
				del markerPos[words[1]]
				markerDupe[words[1]] = 1
			else:
				markerChr[words[1]] = c
				markerPos[words[1]] = int(words[2])
		#foreach line
	#with mafFile
#foreach chr
print "... OK: %d markers (%d duplicates dropped)" % (len(markerChr),len(markerDupe))

if rjpedPath.endswith('.rjped'):
	mapPath = rjpedPath[:-6] + '.map'
else:
	mapPath = rjpedPath + '.map'
print "generating map file '%s' ..." % (mapPath)
with open(rjpedPath, 'rU') as rjpedFile:
	markerList = rjpedFile.next().split()[1:]
#with rjpedFile
l = 0
with open(mapPath, 'w') as mapFile:
	for marker in markerList:
		if marker in markerChr:
			l += 1
			mapFile.write("%d %s %d\n" % (markerChr[marker],marker,markerPos[marker]))
	#foreach marker
#with mapFile
print "... OK: %d lines" % l
