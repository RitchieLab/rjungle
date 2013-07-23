#!/usr/bin/env python

import sys
import os


def ftail(filePtr, numLines=10, lineEnd='\n', bufferSize=1024):
	numLines,bufferSize = int(numLines),int(bufferSize)
	filePtr.seek(0, 2)
	fSize = filePtr.tell()
	lines = 0
	block = -1
	data = []
	# read blocks from the end until we have enough lines
	while lines < numLines and fSize > 0:
		if (fSize > bufferSize):
			filePtr.seek(block * bufferSize, 2)
			chunk = filePtr.read(bufferSize)
		else:
			filePtr.seek(0, 0)
			chunk = filePtr.read(fSize)
		fSize -= bufferSize
		if block == -1 and chunk.endswith(lineEnd):
			chunk = chunk[:-len(lineEnd)]
		lines += chunk.count(lineEnd)
		block -= 1
		data.append(chunk)
	# reverse the blocks in-place
	for d in range(0, len(data) // 2):
		data[d],data[-d-1] = data[-d-1],data[d]
	# rejoin the blocks and then resplit on line boundaries
	return ''.join(data).split(lineEnd)[-numLines:]
#ftail()


def go(basePath, topStr, varList):
	# validate top
	topList = topStr.split(',')
	for n,top in enumerate(topList):
		if not str(top).isdigit():
			sys.stderr.write("ERROR: invalid top %s, expected integer" % top)
			exit(2)
		topList[n] = int(top)
	
	# validate variables and initialize storage
	varImp = {}
	varTop = {}
	for token in varList:
		for var in token.split(','):
			varImp[var] = 0.0
		varTop[token] = {}
		for top in topList:
			varTop[token][top] = 0
	
	# process each importance file
	pMax = 0
	for p in range(1, 101):
		pVarImp = {}
		filePath = os.path.join(basePath, str(p), 'out.importance')
		with open(filePath,'r') as impFile:
			# verify file header
			impFile.seek(0, 0)
			header = impFile.next()
			if not header.startswith("iteration id varname "):
				sys.stderr.write("ERROR: unknown file header for %s: %s\n" % (filePath, header))
				continue
			# identify last iteration for this run
			lastLine = ftail(impFile, 1, bufferSize=128)[-1]
			lastIter = int(lastLine.split()[0])
			# read all variable importances for the last iteration
			impFile.seek(0, 0)
			impFile.next()
			for line in impFile:
				words = line.split()
				if int(words[0]) == lastIter:
					pVarImp[words[2]] = float(words[3])
		#with impFile
		
		# tally importance scores
		for var in varImp:
			if var not in pVarImp:
				sys.stderr.write("WARNING: no importance data for var %s in %s\n" % (var, filePath))
			else:
				varImp[var] += pVarImp[var]
		
		# for each top, identify importance cutoff and tally hits
		impList = sorted(pVarImp.values(), reverse=True)
		for top in topList:
			if top > len(impList):
				sys.stderr.write("WARNING: top %d greater than var count %d in %s\n" % (top, len(impList), filePath))
			else:
				cut = impList[top - 1]
				for token in varList:
					hit = True
					for var in token.split(','):
						if var not in pVarImp:
							sys.stderr.write("WARNING: no importance data for var %s in %s\n" % (var, filePath))
							hit = False
							break
						elif pVarImp[var] < cut:
							hit = False
							break
					if hit:
						varTop[token][top] += 1
			#if top valid
		#foreach top
		
		pMax += 1
	#for p=1..100
	
	# average importance scores
	for var in varImp:
		varImp[var] /= 100.0
	for token in varList:
		if token not in varImp:
			varImp[token] = 0.0
			tokVars = token.split(',')
			for var in tokVars:
				varImp[token] += varImp[var]
			varImp[token] /= len(tokVars)
	
	# print results
	sys.stdout.write("\tavg")
	for top in topList:
		sys.stdout.write("\ttop%d" % top)
	sys.stdout.write("\n")
	for token in varList:
		sys.stdout.write("%s\t%f" % (token,varImp[token]))
		for top in topList:
			sys.stdout.write("\t%d" % varTop[token][top])
		sys.stdout.write("\n")
#go()


if __name__ == "__main__":
	if len(sys.argv) <= 3:
		name = os.path.basename(sys.argv[0])
		sys.stderr.write("""
usage: %s <basepath> <top> <var> [var] ...
  <basepath> = path prefix in which '1/' '2/' etc. are found
  <top> = #(s) for which to tally occurences of variables in the top #
  <var> = variables to tally and summarize (i.e. 'RL0.42' or 'C7')
All numbered directories inside <basepath> will be scanned for the RJ output
file 'out.importance'. Those files will be analyzed, computing average
importance scores for all named <var>s as well as tallies of occurrences of
<var>s in the <top>s. You may specify multiple values for <top> and each <var>
using commas.

Example:
  %s /path/ 10,20 C1 C2 C1,C2
will generate a 3x3 table:
         avg  top10  top20
  C1
  C2
  C1,C2
where the avg of C1,C2 is the average of the individual averages, and the
topX of C1,C2 is the number of occurrences of C1 and C2 together in the topX.
""" % (name,name))
		exit(1)
	go(sys.argv[1], sys.argv[2], sys.argv[3:])
