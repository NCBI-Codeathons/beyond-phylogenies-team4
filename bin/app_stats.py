#!/usr/bin/env python

import sys
import csv
import glob
import os.path

import flacodb as fdb

SPIKE = [21563, 25384]
SPIKELEN = 3821
DEPTHS = [20, 50, 100, 200]

def frac_depth(filename, depth):
    nbases = 0
    ngood = 0

    with open(filename, "r") as f:
        c = csv.reader(f, delimiter="\t")
        for row in c:
            if row[0] == "MN908947.3":
                start = int(row[1])
                end = int(row[2])
                size = end - start
                d = int(row[3])
                nbases += size
                if d >= depth:
                    ngood += size
    if nbases:
        return 100.0*ngood/nbases
    else:
        return 0.0

def get_consensus_stats(consfile):
    seq = ""
    longestn = 0
    with open(consfile, "r") as f:
        f.readline()            # skip header
        for line in f:
            seq += line.rstrip("\r\n").upper()
    spike_n = seq[SPIKE[0]:SPIKE[1]].count("N")

    # Find length of longest stretch of Ns
    instretch = False
    stretchlen = 0
    for b in seq:
        if instretch:
            if b == "N":
                stretchlen += 1
            else:
                longestn = max(longestn, stretchlen)
                stretchlen = 0
                instretch = False
        else:
            if b == "N":
                stretchlen = 1
                instretch = True
                
    return ( 100 - 100.0 * spike_n / SPIKELEN, stretchlen )

def main(run, sample, sampledir, alnstats, hidepth):
    DB = fdb.FLACOdb()

    alnstats = [int(x) for x in alnstats.split(":")]
    if alnstats[0]:
        covidpct = 100.0*alnstats[2]/alnstats[0]
    else:
        covidpct = 0.0

    data = [
        ("ITotalReads", alnstats[0]),
        ("IHsReads", alnstats[1]),
        ("ICovReads", alnstats[2]),
        ("FCovPct", covidpct)]

    coveragefiles = glob.glob(sampledir + "/*.pathogen-coverage-report.tsv")
    if coveragefiles:
        other = ""
        with open(coveragefiles[0], "r") as f:
            for line in f:
                covdata = line.rstrip("\r\n").split("\t")
                if covdata[0] == "SARS-CoV-2":
                    avgdepth = int(covdata[1]) / 30000.0
                    data.append( ("FAvgDepth", avgdepth) )
                    data.append( ("ITotalCoverage", covdata[1]) )
                    data.append( ("IMedianDepth", covdata[4]) )
                else:
                    other += ":".join(covdata) + "|"
        data.append( ("TOtherCov", other) )

    bedfiles = glob.glob(sampledir + "/*.pathogen_full_res.bed")
    if bedfiles:
        bed = bedfiles[0]
        for dep in DEPTHS:
            fracdepth = frac_depth(bedfiles[0], dep)
            data.append( ("FHiDepthFrac" + str(dep), fracdepth) )

    covfiles = glob.glob(sampledir + "/consensus/*.fasta")
    if covfiles:
        consensus_stats = get_consensus_stats(covfiles[0])
        data.append( ("FSpikeCov", consensus_stats[0]) )
        data.append( ("IMaxNStretch", consensus_stats[1]) )

    pngfiles = glob.glob(sampledir + "/*_plotcov.png")
    if pngfiles:
        data.append( ("TPlotCov", sample + "." + run + ".plotcov.png") ) # collect-stats renames the PNG file

    DB.storeValues(sample, run, data)

#    sys.stdout.write(sample + "\t" + "\t".join([str(x) for x in alnstats]) + "\t" + covidpct + "\t" + "\t".join(covdata[1:]) + "\t{}\t{:.2f}\n".format(avgdepth, fracdepth))

if __name__ == "__main__":
    run = sys.argv[1]
    sample = sys.argv[2]
    sampledir = sys.argv[3]
    alnstats = sys.argv[4]
    if len(sys.argv) > 5:
        hidepth = sys.argv[5]
    else:
        hidepth = 250
    main(run, sample, sampledir, alnstats, hidepth)

