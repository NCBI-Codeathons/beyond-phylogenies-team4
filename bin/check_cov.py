#!/usr/bin/env python

import sys
import csv
import glob

import flacodb as fdb

SPIKE = [21563, 25384]
SPIKELEN = 3821

def frac_depth(filename, depth, spike=False):
    """Returns the fraction of bases with depth higher than `depth' in BED file `filename'."""
    nbases = 0
    ngood = 0
    spikegood = 0

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
                    if SPIKE[0] <= start < end <= SPIKE[1]:
                        spikegood += size
    if nbases:
        if spike:
            return 100.0*spikegood/SPIKELEN
        else:
            return 100.0*ngood/nbases
    else:
        return 0.0

def main(levels, samplesfiles):
    DB = fdb.FLACOdb()
    lvls = [int(l) for l in levels.split(",")]
    sys.stdout.write("Sample\tAvgDepth\tMedianDepth\t" + "\t" .join(["{}X".format(l) for l in lvls]) + "\t" + "\t".join(["spike({}X)".format(l) for l in lvls]) + "\n")
    DB.opendb()
    for smpfile in samplesfiles:
        run = smpfile.split("-")[0]
        with open(smpfile, "r") as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                sample = row[0]
                sampledir = "AppResults/{}/Files/".format(row[1])
                coveragefiles = glob.glob(sampledir + "/*.pathogen_full_res.bed")
                if coveragefiles:
                    smpdata = DB.getSampleData(sample, run)
                    filename = coveragefiles[0]
                    vals =   ["{:.2f}%".format(frac_depth(filename, lvl)) for lvl in lvls]
                    spvals = ["{:.2f}%".format(frac_depth(filename, lvl, spike=True)) for lvl in lvls]
                    avgd = smpdata['AvgDepth'] if 'AvgDepth' in smpdata else 0.0
                    mdep = smpdata['MedianDepth'] if 'MedianDepth' in smpdata else 0
                    sys.stdout.write("{}\t{:.1f}\t{}\t{}\t{}\n".format(sample, avgd, mdep, "\t".join(vals), "\t".join(spvals)))
    DB.closedb()

if __name__ == "__main__":
    args = sys.argv[1:]
    main(args[0], args[1:])

