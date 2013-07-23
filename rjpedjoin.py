#!/usr/bin/env python

print "identifying common samples ..."
sampleCount = {}
for c in xrange(1,23):
	with open('chr%d.samples' % c, 'rU') as sampleFile:
		# skip header
		sampleFile.next()
		for line in sampleFile:
			sample = line.strip()
			if sample in sampleCount:
				sampleCount[sample] += 1
			else:
				sampleCount[sample] = 1
sampleSet = set()
for sample in sampleCount:
	if sampleCount[sample] == 22:
		sampleSet.add(sample)
print "... OK: %d samples" % len(sampleSet)

print "joining rjped files ..."
chrFile = [ open('all.rjped','w') ]
chrLine = [ None ]
for c in xrange(1,23):
	chrFile.append( open('chr%d.rjped' % c, 'rU') )
	chrLine.append( chrFile[-1].next().rstrip("\r\n") )
	# join the headers
	if c == 1:
		chrFile[0].write(chrLine[c])
	else:
		s = chrLine[c].find(' ')
		chrFile[0].write(chrLine[c][s:])
#foreach chr
chrFile[0].write("\n")
lines = 0
try:
	while True:
		lines += 1
		# from each file, read lines until we get a sample from the common set
		for c in xrange(1,23):
			while True:
				line = chrFile[c].next().rstrip(" \r\n")
				s = line.find(' ')
				if line[0:s] in sampleSet:
					if c == 1:
						sample = line[0:s]
					elif sample != line[0:s]:
						print "ERROR: sample order mismatch; chr%d '%s' vs chr%d '%s'" % (1,sample,c,line[0:s])
						sys.exit(1)
					chrLine[c] = line
					break
		# add the lines to the output
		chrFile[0].write(sample)
		for c in xrange(1,23):
			s = chrLine[c].find(' ')
			chrFile[0].write(chrLine[c][s:])
		#foreach chr
		chrFile[0].write("\n")
	#forever until StopIteration
except StopIteration:
	pass
chrFile[0].close()
for c in xrange(1,23):
	n = 0
	try:
		while True:
			chrFile[c].next()
			n += 1
	except StopIteration:
		pass
	chrFile[c].close()
	if n > 0:
		print "    WARNING: chr%d.rjped has %d leftover lines" % (c,n)
print "... OK: %d lines" % lines
