#!/usr/bin/env python

import sys

if len(sys.argv) < 2:
	print "usage: %s <rjpedfile>" % sys.argv[0]
	sys.exit(1)
rjpedPath = sys.argv[1]

print "reading MAF files ..."
markerChr = {}
markerPos = {}
markerAll = {}
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
				del markerAll[words[1]]
				markerDupe[words[1]] = 1
			else:
				markerChr[words[1]] = c
				markerPos[words[1]] = int(words[2])
				markerAll[words[1]] = {
					'-1' : " 0 0",
					'0' : (" %s %s" % (words[4],words[4])),
					'1' : (" %s %s" % (words[3],words[4])),
					'2' : (" %s %s" % (words[3],words[3]))
				}
		#foreach line
	#with mafFile
#foreach chr
print "... OK: %d markers (%d duplicates dropped)" % (len(markerChr),len(markerDupe))

print "reading phenotype files ..."
sampleSex = {}
samplePhe = {}
with open('CAPandPrince_IDSexHDL.txt') as pheFile:
	header = pheFile.next().strip()
	if header != "ID\tSex\tmHDLC_SVEV":
		print "ERROR: unrecognized file header: %s" % header
		sys.exit(1)
	for line in pheFile:
		words = line.rstrip(" \r\n").split("\t")
		if words[0] in sampleSex:
			print "ERROR: duplicate sample %s" % words[0]
			sys.exit(1)
		sampleSex[words[0]] = words[1]
		samplePhe[words[0]] = words[2]
	#foreach line
#with pheFile
print "... OK: %d samples" % len(sampleSex)

if rjpedPath.endswith('.rjped'):
	pedPath = rjpedPath[:-6] + '.ped'
else:
	pedPath = rjpedPath + '.ped'
print "converting '%s' to '%s' ..." % (rjpedPath,pedPath)
with open(rjpedPath, 'rU') as rjpedFile:
	with open(pedPath, 'w') as pedFile:
		markerList = rjpedFile.next().split()[1:]
		l = 1
		for line in rjpedFile:
			line = line.rstrip(" \r\n")
			l += 1
			s1 = line.find(' ')
			sample = line[0:s1]
			if sample in sampleSex:
				pedFile.write("%s %s 0 0 %s %s" % (sample,sample,sampleSex[sample],samplePhe[sample]))
				c = 1
				for marker in markerList:
					c += 1
					if s1 < 0:
						print "ERROR: not enough input columns at line %d column %d" % (l,c)
						sys.exit(1)
					s0 = s1 + 1
					s1 = line.find(' ',s0)
					if marker in markerAll:
						if s1 < 0:
							pedFile.write(markerAll[marker][line[s0:]])
						else:
							pedFile.write(markerAll[marker][line[s0:s1]])
				pedFile.write("\n")
			#if sample ok
		#foreach input line
	#with pedFile
#with rjpedFile
print "... OK: %d lines" % l
