#!/usr/bin/env python

import sys
import argparse
import math
import time

#from zfile import zopen

def bimbam2dose(samplePath, genoPath, phenoPath, minMAF, minSampleConf, minMarkerConf, sampleOutPath, markerOutPath, outputPath, quiet):
	msg = (sys.stderr if ((not outputPath) or (outputPath == '-')) else sys.stdout)
	
	# process sample file (no header!)
	if not quiet:
		msg.write("reading sample file ...")
		msg.flush()
	t0 = time.time()
	samples = []
	with open(samplePath,'rU') as sampleFile:
		for line in sampleFile:
			samples.append(line.strip())
	numSamples = len(samples)
	t1 = time.time()
	if not quiet:
		msg.write(" OK: %d samples in %1.1fs\n" % (numSamples,t1-t0))
	
	# process phenotype file (if any)
	samplePheno = phenotypes = None
	if phenoPath:
		if not quiet:
			msg.write("reading phenotype file...")
			msg.flush()
		t0 = time.time()
		samplePheno = {}
		phenotypes = []
		with open(phenoPath,'rU') as phenoFile:
			for line in phenoFile:
				words = line.rstrip().split()
				sample = words[0]
				pheno = words[1]
				
				samplePheno[sample] = pheno
				phenotypes.append(pheno)
			#foreach line/phenotype
		#with phenoFile
		numSampleNoPheno = len(set(samples) - set(samplePheno.keys()))
		t1 = time.time()
		if not quiet:
			msg.write(" OK: %d phenotypes in %1.1fs\n" % (len(phenotypes),t1-t0))
			if numSampleNoPheno > 0:
				msg.write("WARNING: no phenotype for %d of %d samples (%1.1f%%)\n" % (numSampleNoPheno,numSamples,100.0*numSampleNoPheno/numSamples))
	#if phenoPath
	
	# analyze MAF and average sample/marker confidence in the genotype file
	if not quiet:
		msg.write("analyzing genotype file ...")
		msg.flush()
	t0 = time.time()
	numGenoCols = 3 + (2 * numSamples)
	sampleRange = [ s for s in xrange(0,numSamples) ]
	markers = []
	markerMin = []
	markerMaj = []
	markerConf = []
	markerMAF = []
	sampleConf = [ 0.0 for s in sampleRange ]
	with open(genoPath,'rU') as genoFile:
		l = 0
		for line in genoFile:
			l += 1
			words = line.rstrip().split()
			if len(words) != numGenoCols:
				msg.write("ERROR: expected %d columns at line %d, found %d!\n" % (numGenoCols,l,len(words)))
				sys.exit(1)
			markers.append(words[0])
			markerMin.append(words[1])
			markerMaj.append(words[2])
			markerConf.append(0.0)
			markerMAF.append(0.0)
			m = len(markers) - 1
			for s in sampleRange:
				g2 = float(words[3+2*s])
				g1 = float(words[4+2*s])
				g0 = 1.0 - g2 - g1
				conf = max(g0,g1,g2)
				sampleConf[s] += conf
				markerConf[m] += conf
				markerMAF[m] += (2.0*g2 + 1.0*g1)
			markerMAF[m] /= 2.0*numSamples
		#foreach line/marker
	#with genoFile
	numMarkers = len(markers)
	markerRange = [ m for m in xrange(0,numMarkers) ]
	t1 = time.time()
	if not quiet:
		msg.write(" OK: %d markers in %1.1fs\n" % (numMarkers,t1-t0))
	
	# compute sample stats
	if not quiet:
		msg.write("generating sample statistics (%s minimum confidence) ..." % (minSampleConf,))
		msg.flush()
	t0 = time.time()
	for s in sampleRange:
		sampleConf[s] /= float(numMarkers)
	sampleOK = [ (sampleConf[s] >= minSampleConf) for s in sampleRange ]
	numSampleOK = sampleOK.count(True)
	numSampleDrop = sampleOK.count(False)
	if sampleOutPath:
		with open(sampleOutPath,'w') as sampleOutFile:
			sampleOutFile.write("sample_id\tavg_confidence\n")
			for s in sampleRange:
				sampleOutFile.write("%s\t%f\n" % (samples[s],sampleConf[s]))
		#with sampleOutFile
	t1 = time.time()
	if not quiet:
		msg.write(" OK: %d of %d samples dropped (%1.1f%%) in %1.1fs\n" % (numSampleDrop,numSamples,100.0*numSampleDrop/numSamples,t1-t0))
	
	# compute marker stats
	if not quiet:
		msg.write("generating marker statistics (%s minimum confidence, %s minimum MAF) ..." % (minMarkerConf,minMAF))
		msg.flush()
	t0 = time.time()
	for m in markerRange:
		markerConf[m] /= float(numSamples)
	markerOK = [ (markerConf[m] >= minMarkerConf and markerMAF[m] >= minMAF) for m in markerRange ]
	numMarkerOK = markerOK.count(True)
	numMarkerDrop = markerOK.count(False)
	if markerOutPath:
		with open(markerOutPath,'w') as markerOutFile:
			markerOutFile.write("marker_id\tavg_confidence\tMAF\n")
			for m in markerRange:
				markerOutFile.write("%s\t%f\t%f\n" % (markers[m],markerConf[m],markerMAF[m]))
		#with markerOutFile
	t1 = time.time()
	if not quiet:
		msg.write(" OK: %d of %d markers dropped (%1.1f%%) in %1.1fs\n" % (numMarkerDrop,numMarkers,100.0*numMarkerDrop/numMarkers,t1-t0))
	
	if outputPath != False:
		# The input has a row per marker and two columns per sample (person),
		# but the output needs a row per sample (person) and one column per
		# marker.  Since the file might be too big to read and transpose
		# entirely in memory, we'll have to re-read it several times, converting
		# X input columns into output rows each iteration.  In practice we need
		# about 72 bytes per sample-marker, so we'll process X columns per pass
		# such that we need ~1gb memory.
		samplesPerPass = min(numSamples, int(math.floor((1024 * 1024 * 1024) / (numMarkerOK * 72))))
		numPasses = int(math.ceil(float(numSamples) / float(samplesPerPass)))
		samplesPerPass = int(math.ceil(float(numSamples) / float(numPasses)))
		
		# process genotype file
		if not quiet:
			msg.write("processing genotype file (%d samples per pass over %d passes):\n" % (samplesPerPass,numPasses))
		t0 = time.time()
		with (sys.stdout if ((not outputPath) or (outputPath == '-')) else open(outputPath,'w')) as out:
			out.write("sample_id")
			if samplePheno:
				out.write(" PHENOTYPE")
			for m in markerRange:
				if markerOK[m]:
					out.write(" %s_%s" % (markers[m],markerMin[m]))
			out.write("\n")
			for p in xrange(0,numPasses):
				if not quiet:
					msg.write("  pass %d of %d ..." % (p+1,numPasses))
					msg.flush()
				t0p = time.time()
				s0 = p * samplesPerPass
				s1 = min(numSamples, s0 + samplesPerPass)
				passSampleRange = range(s0,s1)
				if samplePheno:
					output = [ ([samples[s],samplePheno[samples[s]] if samples[s] in samplePheno else "-9"] if sampleOK[s] else None) for s in passSampleRange ]
				else:
					output = [ ([samples[s]] if sampleOK[s] else None) for s in passSampleRange ]
				with open(genoPath,'rU') as genoFile:
					m = 0
					for line in genoFile:
						if markerOK[m]:
							words = line.rstrip().split()
							for s in passSampleRange:
								if sampleOK[s]:
									g2 = float(words[3+2*s])
									g1 = float(words[4+2*s])
									g0 = 1.0 - g2 - g1
									output[s-s0].append("%1.2f" % (2.0*g2 + 1.0*g1))
								#if col's sample is ok
							#foreach col/sample
						#if line's marker is ok
						m += 1
					#foreach line/marker
				#with genoFile
				for s in passSampleRange:
					if sampleOK[s]:
						out.write(" ".join(output[s-s0]))
						out.write("\n")
				t1p = time.time()
				if not quiet:
					msg.write(" OK: %d samples in %1.1fs\n" % (s1-s0,t1p-t0p))
			#foreach pass
		#with output
		t1 = time.time()
		msg.write("conversion complete in %1.1fs\n" % (t1-t0))
	#if convert
#bimbam2dose()


if __name__ == "__main__":
	versName,versMaj,versMin,versRev,versDate = 'bimbam2dose',0,1,0,'2012-03-29'
	version = "%d.%d.%d (%s)" % (versMaj, versMin, versRev, versDate)
	
	# define usage
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description="%s version %s" % (versName, version),
		epilog="""
example: %(prog)s -s samples.txt -g input.geno -p phenotypes.txt -m 0.1 -t 0.5 -T 0.5 -S out.samples -M out.markers -o out.dose
"""
	)
	parser.add_argument('-s', '--samples', action='store', type=str, metavar='file', required=True,
		help="sample list file (in order of genotype file columns)"
	)
	parser.add_argument('-g', '--genotypes', action='store', type=str, metavar='file', required=True,
		help="BIMBAM genotype file"
	)
	parser.add_argument('-p', '--phenotypes', action='store', type=str, metavar='file',
		help="phenotype file (default: none)"
	)
	parser.add_argument('-m', '--maf', action='store', type=float, metavar='float', default=0.0,
		help="minimum minor allele frequency, below which the marker will be dropped (default: 0)"
	)
	parser.add_argument('-t', '--sample-confidence', action='store', type=float, metavar='float', default=0.0,
		help="minimum average sample confidence, below which the sample will be dropped (default: 0)"
	)
	parser.add_argument('-T', '--marker-confidence', action='store', type=float, metavar='float', default=0.0,
		help="minimum average marker confidence, below which the marker will be dropped (default: 0)"
	)
	parser.add_argument('-S', '--sample-output', action='store', type=str, metavar='file',
		help="output file for average sample confidence scores (default: none)"
	)
	parser.add_argument('-M', '--marker-output', action='store', type=str, metavar='file',
		help="output file for average marker confidence scores (default: none)"
	)
	parser.add_argument('-a', '--analyze-only', action='store_true',
		help="only analyze the genotype file for average confidence scores, don't perform conversion to dose file"
	)
	parser.add_argument('-o', '--output', action='store', type=str, metavar='file',
		help="output file for dosages (default: stdout)"
	)
	parser.add_argument('-q', '--quiet', action='store_true',
		help="don't print progress messages"
	)
	parser.add_argument('--version', action='version', version=version)
	
	# parse arguments
	args = parser.parse_args()
	outputPath = False if args.analyze_only else args.output
	
	# run
	bimbam2dose(
			args.samples, args.genotypes, args.phenotypes,
			args.maf,
			args.sample_confidence, args.marker_confidence,
			args.sample_output, args.marker_output,
			outputPath, args.quiet
	)
#__main__
