#!/usr/bin/env python

import sys
import argparse
import math

#from zfile import zopen


def bimbam2ped(samplePath, genoPath, cutoff, maf, sampleThreshold, markerThreshold, outputPath, sampleDropPath, markerDropPath, quiet, resume):
	msg = (sys.stderr if ((not outputPath) or (outputPath == '-')) else sys.stdout)
	markers = markerOK = samples = sampleOK = resumeHeader = resumeLines = None
	
	# process sample file (no header!)
	if not quiet:
		msg.write("reading sample file ...")
		msg.flush()
	samples = []
	with open(samplePath,'rU') as sampleFile:
		for line in sampleFile:
			samples.append(line.strip())
	if not quiet:
		msg.write(" OK: %d samples\n" % len(samples))
	
	if resume:
		raise Exception("NOT IMPLEMENTED") #TODO: have to read every marker from .geno file, can't just read non-dropped markers from the output header
		
		# read the sampleDrop file
		while (not sampleOK) and sampleDropPath: #while so we can break
			try:
				if not quiet:
					msg.write("reading previous sample-drop file ...")
					msg.flush()
				with open(sampleDropPath,'rU') as sampleDropFile:
					header = sampleDropFile.next().strip()
					if header != "sample_id drop_reason":
						msg.write("WARNING: unrecognized sample-drop file header: %s\n" % header)
						break
					sampleDropSet = set()
					for line in sampleDropFile:
						sample = line.split(None,1)[0]
						sampleDropSet.add(sample)
				#with sampleDropFile
				sampleOK = []
				for sample in samples:
					sampleOK.append(-1 if sample in sampleDropSet else 1)
				if not quiet:
					msg.write(" OK: %d samples were dropped\n" % len(sampleDropSet))
			except IOError:
				if not quiet:
					msg.write(" ERROR: file could not be read\n")
				sampleOK = None
			break
		
		# read the previous output file if any, to count lines and get the markers out of the header
		while outputPath: #while so we can break
			try:
				if not quiet:
					msg.write("reading previous output file ...")
					msg.flush()
				with open(outputPath,'rU') as resumeFile:
					resumeHeader = resumeFile.next().strip()
					if not resumeHeader.startswith("sample_id "):
						msg.write("WARNING: unrecognized output file header\n")
						break
					markers = resumeHeader.split()[1:]
					markerOK = [1 for marker in markers]
					resumeLines = 0
					for line in resumeFile:
						resumeLines += 1
				#with resumeFile
				if not quiet:
					msg.write(" OK: %d lines were written\n" % resumeLines)
			except IOError:
				if not quiet:
					msg.write(" ERROR: file could not be read\n")
				markers = markerOK = resumeHeader = resumeLines = None
			break
		
		# if we didn't get the markerOK list from the output header, use the geno and markerDrop files instead
		while (not markers) and markerDropPath: #while so we can break
			try:
				if not quiet:
					msg.write("reading previous marker-drop file ...")
					msg.flush()
				with open(markerDropPath) as markerDropFile:
					header = markerDropFile.next().strip()
					if header != "marker drop_reason":
						msg.write("WARNING: unrecognized marker-drop file header: %s\n" % header)
						break
					markerDropSet = set()
					for line in markerDropFile:
						marker = line.split(None,1)[0]
						markerDropSet.add(marker)
					if not quiet:
						msg.write(" OK: %d markers were dropped\n" % len(markerDropSet))
					with open(genoPath,'rU') as genoFile:
						if not quiet:
							msg.write("reading bimbam genotype file ...")
							msg.flush()
						markers = []
						markerOK = []
						for line in genoFile:
							marker = line.split(None,1)[0]
							markers.append(marker)
							markerOK.append(-1 if marker in markerDropSet else 1)
						if not quiet:
							msg.write(" OK: %d markers\n" % len(markers))
					#with genoFile
				#with markerDropFile
			except IOError:
				if not quiet:
					msg.write(" ERROR: file could not be read\n")
				markers = markerOK = None
			break
		
		# if we're missing anything, abort the resume
		if (not samples) or (not sampleOK) or (not markers) or (not markerOK):
			msg.write("ERROR: could not identify dropped samples and markers, unable to resume\n")
			sys.exit(1)
		
		# compute the stats
		if not quiet:
			msg.write("resumed statistics:\n")
		numMarkerOK = markerOK.count(1)
		numMarkerDrop = len(markers) - numMarkerOK
		numSampleOK = sampleOK.count(1)
		numSampleDrop = len(samples) - numSampleOK
		if not quiet:
			msg.write("  kept %d samples, %d markers\n" % (numSampleOK,numMarkerOK))
			msg.write("  dropped %d samples, %d markers\n" % (numSampleDrop,numMarkerDrop))
		
	else:
		
		# run through the bimbam/geno file, tallying up hits and misses for each row/col
		if not quiet:
			msg.write("analyzing bimbam genotype file ...")
			msg.flush()
		numInputCols = 3 + (2 * len(samples))
		sampleRange = range(0,len(samples))
		sampleMarkersOK = [0 for s in sampleRange]
		markers = []
		markerSamplesOK = []
		markerMinorAlleles = []
		markerGenotype = [-1 for s in sampleRange]
		markerOK = []
		with open(genoPath,'rU') as genoFile:
			for line in genoFile:
				line = line.split()
				if len(line) != numInputCols:
					msg.write("ERROR: expected %d columns at line %d, found %d!\n" % (numInputCols,len(markerOK)+1,len(line)))
					sys.exit(1)
				markers.append(line[0])
				markerSamplesOK.append(0)
				markerMinorAlleles.append(0)
				for s in sampleRange:
					f2 = float(line[3+2*s])
					f1 = float(line[4+2*s])
					if f2 >= cutoff:
						markerGenotype[s] = 2
					elif f1 >= cutoff:
						markerGenotype[s] = 1
					elif (1.0 - f2 - f1) >= cutoff:
						markerGenotype[s] = 0
					else:
						markerGenotype[s] = -1
					if markerGenotype[s] >= 0:
						markerSamplesOK[-1] += 1
						markerMinorAlleles[-1] += markerGenotype[s]
				if (float(markerSamplesOK[-1]) / float(len(samples))) < markerThreshold:
					markerOK.append(-1)
				elif (float(markerMinorAlleles[-1]) / (2 * float(markerSamplesOK[-1]))) < maf:
					markerOK.append(-2)
				else:
					markerOK.append(1)
					for s in sampleRange:
						if markerGenotype[s] >= 0:
							sampleMarkersOK[s] += 1
			if not quiet:
				msg.write(" OK: %d markers\n" % len(markers))
		#with bimbam/geno file
		
		# compute stats for the main run(s)
		if not quiet:
			msg.write("generating statistics (%1.2f cutoff, %1.2f maf, %1.2f s-threshold, %1.2f m-threshold) ..." % (cutoff,maf,sampleThreshold,markerThreshold))
			msg.flush()
		numMarkerOK = markerOK.count(1)
		numMarkerDropGeno = markerOK.count(-1)
		numMarkerDropMAF = markerOK.count(-2)
		if (numMarkerOK + numMarkerDropGeno + numMarkerDropMAF) != len(markers):
			msg.write("ERROR: marker miscount!  %d + %d + %d != %d\n" % (numMarkerOK,numMarkerDropGeno,numMarkerDropMAF,len(markers)))
			sys.exit(1)
		sampleOK = []
		for s in sampleRange:
			if (float(sampleMarkersOK[s]) / float(numMarkerOK)) < sampleThreshold:
				sampleOK.append(-1)
			else:
				sampleOK.append(1)
		numSampleOK = sampleOK.count(1)
		numSampleDrop = sampleOK.count(-1)
		if (numSampleOK + numSampleDrop) != len(samples):
			msg.write("ERROR: sample miscount!  %d + %d != %d\n" % (numSampleOK,numSampleDrop,len(samples)))
			sys.exit(1)
		if sampleDropPath:
			with open(sampleDropPath,'w') as sampleDropFile:
				sampleDropFile.write("sample_id drop_reason\n")
				for s in range(0,len(samples)):
					if sampleOK[s] == -1:
						sampleDropFile.write("%s \"classified %d of %d markers (%1.6f)\"\n" % (samples[s],sampleMarkersOK[s],numMarkerOK,(float(sampleMarkersOK[s])/float(numMarkerOK))))
					elif sampleOK[s] != 1:
						sampleDropFile.write("%s ?\n" % samples[s])
		if markerDropPath:
			with open(markerDropPath,'w') as markerDropFile:
				markerDropFile.write("marker drop_reason\n")
				for m in range(0,len(markers)):
					if markerOK[m] == -1:
						markerDropFile.write("%s \"classified %d of %d samples (%1.6f)\"\n" % (markers[m],markerSamplesOK[m],len(samples),(float(markerSamplesOK[m])/float(len(samples)))))
					elif markerOK[m] == -2:
						markerDropFile.write("%s \"minor allele frequency %1.6f\"\n" % (markers[m],(float(markerMinorAlleles[m])/(2*float(markerSamplesOK[m])))))
					elif markerOK[m] != 1:
						markerDropFile.write("%s ?\n" % markers[m])
		if not quiet:
			msg.write(" OK: kept %d samples, %d markers\n" % (numSampleOK,numMarkerOK))
			msg.write("  dropped %d samples under %1.2f threshold\n" % (numSampleDrop,sampleThreshold))
			msg.write("  dropped %d markers under %1.2f threshold\n" % (numMarkerDropGeno,markerThreshold))
			msg.write("  dropped %d markers under %1.2f MAF\n" % (numMarkerDropMAF,maf))
		
	#if resume
	
	# The input has a row per marker and two columns per sample (person),
	# but the output needs a row per sample (person) and one column per
	# marker.  Since the file might be too big to read and transpose
	# entirely in memory, we'll have to re-read it several times, converting
	# X input columns into output rows each iteration.  In practice we need
	# about 8 bytes per sample-marker (1-char string plus array overhead),
	# so we'll process X columns per pass such that we need ~1gb memory.
	samplesPerPass = min(len(samples), int(math.floor((1024 * 1024 * 1024) / (numMarkerOK * 8))))
	numPasses = int(math.ceil(float(len(samples)) / float(samplesPerPass)))
	
	# process bimbam/geno file
	if not quiet:
		msg.write("processing bimbam genotype file (%d samples per pass over %d passes):\n" % (samplesPerPass,numPasses))
	with (sys.stdout if ((not outputPath) or (outputPath == '-')) else open(outputPath,'a' if resume else 'w')) as out:
		startingPass = 0
		if resumeLines:
			if (resumeLines % samplesPerPass) != 0:
				msg.write("ERROR: previous output file contains %d lines = %1.2f passes\n" % (resumeLines, 1.0*resumeLines/samplesPerPass))
				sys.exit(1)
			startingPass = resumeLines / samplesPerPass
		if not resumeHeader:
			out.write("sample_id")
			for m in range(0,len(markers)):
				if markerOK[m] > 0:
					out.write(" %s" % markers[m])
			out.write("\n")
		for p in range(startingPass,numPasses):
			# initialize output lines
			s0 = p * samplesPerPass
			s1 = min(len(samples), s0 + samplesPerPass)
			sampleRange = range(s0,s1)
			output = [ ([samples[s]] if sampleOK[s] > 0 else None) for s in sampleRange ]
			if not quiet:
				msg.write("  pass %d of %d (%d samples) ..." % (p+1,numPasses,len(output)))
				msg.flush()
			with open(genoPath,'rU') as genoFile:
				m = 0
				for line in genoFile:
					if markerOK[m] > 0:
						line = line.split()
						for s in sampleRange:
							if sampleOK[s] > 0:
								f2 = float(line[3+2*s])
								f1 = float(line[4+2*s])
								if f2 >= cutoff:
									output[s-s0].append("2")
								elif f1 >= cutoff:
									output[s-s0].append("1")
								elif (1.0 - f2 - f1) >= cutoff:
									output[s-s0].append("0")
								else:
									output[s-s0].append("-1")
							#if col's sample is ok
						#foreach col/sample
					#if line's marker is ok
					m += 1
				#foreach line/marker
			#with bimbam/geno file
			for s in sampleRange:
				if sampleOK[s] > 0:
					output[s-s0].append("\n")
					out.write(" ".join(output[s-s0]))
			if not quiet:
				msg.write(" OK\n")
		#foreach pass
	#with outputFile
#bimbam2ped()


if __name__ == "__main__":
	versMaj,versMin,versRev,versDate = 0,1,0,'2012-02-13'
	version = "%d.%d.%d (%s)" % (versMaj, versMin, versRev, versDate)
	
	# define usage
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description="bimbam2ped version %s" % version,
		epilog="""
example: %(prog)s -s samples.txt -b input.bimbam -c .85 -t 0.95 -o output.ped
"""
	)
	parser.add_argument('-s', '--samples', action='store', type=str, metavar='file', required=True,
		help="sample id list (in order of bimbam columns)"
	)
	parser.add_argument('-b', '--bimbam', action='store', type=str, metavar='file', required=True,
		help="input bimbam file"
	)
	parser.add_argument('-c', '--cutoff', action='store', type=float, metavar='float', default=0.85,
		help="minimum confidence score for classification, below which the data point will be coded as missing (default: .85)"
	)
	parser.add_argument('-m', '--maf', action='store', type=float, metavar='float', default=0.05,
		help="minimum minor allele frequency, below which the marker will be dropped (default: .05)"
	)
	parser.add_argument('-t', '--sample-threshold', action='store', type=float, metavar='float', default=0.95,
		help="minimum portion of classified markers, below which an entire sample will be dropped (default: .95)"
	)
	parser.add_argument('-T', '--marker-threshold', action='store', type=float, metavar='float', default=0.95,
		help="minimum portion of classified samples, below which an entire marker will be dropped (default: .95)"
	)
	parser.add_argument('-o', '--output', action='store', type=str, metavar='file',
		help="output file for genotype classifications (default: stdout)"
	)
	parser.add_argument('-S', '--sample-drop', action='store', type=str, metavar='file',
		help="output file to list dropped samples (default: none)"
	)
	parser.add_argument('-M', '--marker-drop', action='store', type=str, metavar='file',
		help="output file to list dropped markers (default: none)"
	)
	#parser.add_argument('-r', '--resume', action='store_true',
	#	help="try to resume an interrupted run (requires sample-drop and marker-drop files)"
	#)
	parser.add_argument('-q', '--quiet', action='store_true',
		help="don't print progress messages"
	)
	parser.add_argument('--version', action='version', version=version)
	
	# parse arguments
	args = parser.parse_args()
	
	# run join
	bimbam2ped(args.samples, args.bimbam, args.cutoff, args.maf, args.sample_threshold, args.marker_threshold, args.output, args.sample_drop, args.marker_drop, args.quiet, args.resume)
#__main__
