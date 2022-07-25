#!/usr/bin/env python

import sys
import csv
from Bio import SeqIO

from flacodb import FLACOdb

class Sample(object):
    name = ""
    sampling_date = ""
    city = ""
    
    def __init__(self, name, db):
        self.name = name
        query = "SELECT Location, SamplingDate FROM Samples WHERE SampleName=?;"
        row = db.execute(query, name).fetchone()
        if row:
            self.city = row[0]
            self.sampling_date = row[1]
        else:
            sys.stderr.write("Error: no metadata for sample `{}'.\n".format(name))

    def seqname(self):
        if self.name.startswith("Haiti"):
            return "hCoV-19/Haiti/UF-{}/{}|{}|{}|NorthAmerica".format(self.name, self.sampling_date[:4], self.name, self.sampling_date)
        else:
            return "hCoV-19/USA/FL-{}/{}|{}|{}|NorthAmerica".format(self.name, self.sampling_date[:4], self.name, self.sampling_date)

def makeFasta(infiles):
    DB = FLACOdb()
    DB.opendb()
    try:
        for infile in infiles:
            sample, filename = infile.split(":")
            smp = Sample(sample, DB)
            if smp:
                sys.stdout.write(">{}\n".format(smp.seqname()))
                with open(filename, "r") as f:
                    f.readline()
                    for line in f:
                        sys.stdout.write(line)
    finally:
        DB.closedb()

def makeOneFasta(filename):
    DB = FLACOdb()
    DB.opendb()
    try:
        with open(filename, "rt") as f:
            for line in f:
                if line[0] == '>':
                    sample = line.rstrip("\r\n")[1:]
                    smp = Sample(sample, DB)
                    if smp:
                        sys.stdout.write(">{}\n".format(smp.seqname()))
                    else:
                        sys.stdout.write(line)
                else:
                    sys.stdout.write(line)
    finally:
        DB.closedb()

if __name__ == "__main__":
    args = sys.argv[1:]
    if args[0] == '-f':
        makeOneFasta(args[1])
    else:
        makeFasta(args)
