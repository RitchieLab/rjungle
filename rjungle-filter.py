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


def go(source, result, top, prefix):
	# process result file
	sys.stderr.write("reading %s for %d best vars matching '%s' ...\n" % (result,top,prefix))
	prefix = prefix.upper()
	varSet = set()
	with open(result,'r') as impFile:
		# verify file header
		header = impFile.next()
		if not header.startswith("iteration id varname "):
			sys.stderr.write("ERROR: unknown results file header: %s\n" % (header))
			exit(1)
		# identify best variables
		for line in ftail(impFile, top):
			var = line.split()[2].upper()
			if prefix.startswith('!'):
				if not var.startswith(prefix[1:]):
					varSet.add(var)
			elif var.startswith(prefix):
				varSet.add(var)
	#with impFile
	sys.stderr.write("... OK: %d matches\n" % len(varSet))
	
	# process source file
	sys.stderr.write("filtering %s ...\n" % source)
	with open(source,'r') as srcFile:
		# identify the column headers we want
		header = srcFile.next()
		colSet = set()
		for c,col in enumerate(header.split()):
			if c == 0:
				sys.stdout.write(col)
			elif col.upper() in varSet:
				sys.stdout.write(" %s" % col)
				colSet.add(c)
		sys.stdout.write("\n")
		# filter each subsequent row to matching columns only
		for line in srcFile:
			for c,val in enumerate(line.split()):
				if c == 0:
					sys.stdout.write(val)
				elif c in colSet:
					sys.stdout.write(" %s" % val)
			sys.stdout.write("\n")
	#with srcFile
	sys.stderr.write("... OK\n")
#go()


if __name__ == "__main__":
	if len(sys.argv) <= 3:
		name = os.path.basename(sys.argv[0])
		sys.stderr.write("""
usage: %s <sourcefile> <resultfile> <top> [prefix]
  <sourcefile> = source data file to filter
  <resultfile> = rjungle importance result file
  <top> = # of variables to filter to
  <prefix> = optional variable prefix to further filter on
The <top> best variables from the <resultfile> will be read and used to
filter the columns of the <sourcefile>.  If a <prefix> is specified, then
among the <top> filtered variables, only those whose names begin with the
<prefix> will be included.  If the prefix begins with the character '!',
then it will be applied in reverse (i.e. '!x' will match variables that do
not begin with 'x').

Example:
  %s mydata.txt rjungle.importance 10 SNP
will generate a filtered copy of 'mydata.txt', including only columns that
match the 10 best variables in 'rjungle.importance' whose names also begin
with 'SNP'.
""" % (name,name))
		exit(1)
	source = sys.argv[1]
	result = sys.argv[2]
	top = int(sys.argv[3])
	prefix = sys.argv[4] if len(sys.argv) > 4 else ''
	go(source, result, top, prefix)
