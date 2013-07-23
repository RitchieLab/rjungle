#!/usr/bin/env python

import sys
import argparse

#TODO: compressed input/output

def ped2rjped(pedPath, mapPath, phePath, pheFormat, allPaths, quiet, out, msg):
	msg.write("reading .map file ...")
	msg.flush()
	markerList = []
	markerMajor = {}
	numDupe = 0
	with open(mapPath,'rU') as mapFile:
		# no header
		for line in mapFile:
			marker = line.split()[1]
			markerList.append(marker)
			if marker in markerMajor:
				numDupe += 1
			markerMajor[marker] = None
	#with mapFile
	msg.write(" OK: %d markers (%d duplicates)\n" % (len(markerList),numDupe))
	
	if phePath:
		msg.write("reading phenotype file ...")
		msg.flush()
		samplePheno = {}
		pheHeader = pheFormat[0]
		pheSample = pheFormat[1] - 1
		phePheno = pheFormat[2] - 1
		numDupe = 0
		numShort = 0
		with open(phePath,'rU') as pheFile:
			for line in pheFile:
				if pheHeader > 0:
					pheHeader -= 1
					continue
				words = line.split()
				if len(words) <= pheSample or len(words) <= phePheno:
					numShort += 1
				elif words[pheSample] in samplePheno:
					numDupe += 1
				else:
					samplePheno[words[pheSample]] = words[phePheno]
			#foreach line
		#with pheFile
		msg.write(" OK: %d phenotypes (%d duplicates, %d lines too short)\n" % (len(samplePheno),numDupe,numShort))
	#if phePath
	
	if allPaths:
		msg.write("reading allele files ...")
		msg.flush()
		numOK = 0
		numDupe = 0
		numExtra = 0
		for allPath in allPaths:
			with open(allPath,'rU') as allFile:
				for line in allFile:
					words = line.split()
					marker = words[1]
					major = words[4]
					if marker not in markerMajor:
						numExtra += 1
					elif markerMajor[marker]:
						numDupe += 1
					else:
						numOK += 1
						markerMajor[marker] = major
				#foreach line
			#with allFile
		#foreach allPath
		msg.write(" OK: %d markers (%d duplicates, %d extras, %d missing)\n" % (numOK,numDupe,numExtra,(len(markerMajor)-numOK)))
	else:
		msg.write("tallying alleles in .ped file ...")
		msg.flush()
		allelesList = [ {} for marker in markerList ]
		with open(pedPath,'rU') as pedFile:
			l = 0
			for line in pedFile:
				l += 1
				# skip 6 columns
				line.rstrip(" \r\n")
				s0 = line.find(' ')
				s0 = line.find(' ',s0+1)
				s0 = line.find(' ',s0+1)
				s0 = line.find(' ',s0+1)
				s0 = line.find(' ',s0+1)
				s1 = line.find(' ',s0+1)
				# tally alleles, two per marker
				m = 0
				first = True
				while m < len(allelesList):
					if s1 < 0:
						msg.write("ERROR: expected %d columns at line %d, but only found %d\n" % (6+2*len(markerList),l,(6 if first else 7)+2*m))
						sys.exit(1)
					# identify allele
					s0 = s1 + 1
					s1 = line.find(' ',s0)
					if s1 < 0:
						a = line[s0:]
					else:
						a = line[s0:s1]
					# tally allele
					if a in allelesList[m]:
						allelesList[m][a] += 1
					else:
						allelesList[m][a] = 1
					# every second allele, tally towards the next marker
					if not first:
						m += 1
					first = not first
				#foreach allele
				if s1 >= 0:
					msg.write("ERROR: expected %d columns at line %d, but found %d+\n" % (6+2*len(markerList),l,(6 if first else 7)+2*m))
					sys.exit(1)
			#foreach line
		#with pedFile
		msg.write(" OK: %d samples\n" % (l,))
		
		msg.write("analyzing allele tallies ...")
		msg.flush()
		m = 0
		numOK = 0
		for marker in markerList:
			total = 0
			major = None
			for a in sorted(allelesList[m]):
				if a != '0':
					major = a
					total += int(allelesList[m][a])
			# no confidence threshold, 50.1% is major
			if major:
				numOK += 1
				markerMajor[marker] = major
			m += 1
		#foreach marker
		msg.write(" OK: %d markers\n" % (numOK,))
	#if/else allPaths
	
	msg.write("converting .ped file to .rjped ...")
	msg.flush()
	out.write("FID IID PAT MAT SEX PHENOTYPE ")
	out.write(" ".join(markerList))
	out.write("\n")
	numMissing = 0
	with open(pedPath,'rU') as pedFile:
		l = 0
		for line in pedFile:
			l += 1
			# copy first 5-6 columns
			line.rstrip(" \r\n")
			s0 = line.find(' ')
			sample = line[0:s0]
			s0 = line.find(' ',s0+1)
			s0 = line.find(' ',s0+1)
			s0 = line.find(' ',s0+1)
			s0 = line.find(' ',s0+1)
			s1 = line.find(' ',s0+1)
			# drop samples with no phenotype
			if phePath:
				if sample not in samplePheno:
					numMissing += 1
					continue
				out.write(line[0:s0+1])
				out.write(samplePheno[sample])
			else:
				out.write(line[0:s1])
			# convert alleles, two per marker
			m = 0
			for marker in markerList:
				if s1 < 0:
					msg.write("ERROR: expected %d columns at line %d, but only found %d\n" % (6+2*len(markerList),l,6+2*m))
					sys.exit(1)
				# identify alleles
				s0 = s1 + 1
				s1 = line.find(' ',s0)
				if s1 < 0:
					msg.write("ERROR: expected %d columns at line %d, but only found %d\n" % (6+2*len(markerList),l,7+2*m))
					sys.exit(1)
				a1 = line[s0:s1]
				s0 = s1 + 1
				s1 = line.find(' ',s0)
				if s1 < 0:
					a2 = line[s0:]
				else:
					a2 = line[s0:s1]
				# convert
				if a1 == a2:
					if a1 == '0':
						out.write(' -1') # 0 0
					elif a1 == markerMajor[marker]:
						out.write(' 0') # M M
					else:
						out.write(' 2') # m m
				elif a1 == '0':
					if a2 == markerMajor[marker]:
						out.write(' 0') # 0 M
					else:
						out.write(' 1') # 0 m
				elif a2 == '0':
					if a1 == markerMajor[marker]:
						out.write(' 0') # M 0
					else:
						out.write(' 1') # m 0
				elif a1 == markerMajor[marker] or a2 == markerMajor[marker]:
					out.write(' 1') # M m , m M
				else:
					out.write(' 2') # m m
				m += 1
			#foreach marker
			if s1 >= 0:
				msg.write("ERROR: expected %d columns at line %d, but found %d+\n" % (6+2*len(markerList),l,6+2*m))
				sys.exit(1)
			out.write("\n")
		#foreach line
	#with pedFile
	msg.write(" OK: %d lines (%d dropped)\n" % (l,numMissing))
#ped2rjped()


if __name__ == "__main__":
	versMaj,versMin,versRev,versDate = 0,1,0,'2012-02-17'
	version = "%d.%d.%d (%s)" % (versMaj, versMin, versRev, versDate)
	
	# define usage
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description="ped2rjped version %s" % version,
		epilog="""
example: %(prog)s -p my.ped -m my.map -t 1 3 2 my.phenotypes -a my.chr*.maf -o my.rjped
"""
	)
	parser.add_argument('-p', '--ped', action='store', type=str, metavar='<file>', required=True,
		help="input .ped file (required)"
	)
	parser.add_argument('-m', '--map', action='store', type=str, metavar='<file>', required=True,
		help="input .map file (required)"
	)
	parser.add_argument('-a', '--alleles', action='append', type=str, nargs='*', metavar='<file>',
		help="input .maf file(s) (default: determine minor alleles from .ped file; slower)"
	)
	parser.add_argument('-t', '--phenotype', action='store', type=str, metavar='<file>',
		help="alternate phenotype file (default: use .ped phenotypes)"
	)
	parser.add_argument('-T', '--phenotype-format', action='store', type=int, nargs=3, metavar='<num>', default=[1,1,2],
		help="alternate phenotype file format, in 3 values: number of header lines; column containing sample identifiers; column containing phenotypes (default: 1 1 2)"
	)
	parser.add_argument('-o', '--output', action='store', type=str, metavar='<file>',
		help="output .rjped file (default: stdout)"
	)
	parser.add_argument('-q', '--quiet', action='store_true',
		help="don't print progress messages"
	)
	parser.add_argument('--version', action='version', version=version)
	
	# parse arguments and run
	args = parser.parse_args()
	allPaths = None
	if args.alleles:
		allPaths = []
		for alleles in args.alleles:
			allPaths.extend(alleles)
	
	with (open(args.output,'w') if (args.output and args.output != '-') else sys.stdout) as out:
		msg = (sys.stdout if (args.output and args.output != '-') else sys.stderr)
		ped2rjped(args.ped, args.map, args.phenotype, args.phenotype_format, allPaths, args.quiet, out, msg)
#__main__
