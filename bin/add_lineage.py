#!/usr/bin/env python

import sys
import csv

def readLineages(filename):
    data = {}
    with open(filename, "r") as f:
        f.readline()
        f.readline()
        c = csv.reader(f, delimiter='\t')
        for row in c:
            if "|" in row[0]:
                smp = row[0].split("|")[1]
                data[smp] = row[1]
    return data

def addLineage(filename, calls):
    with open(filename, "r") as f:
        sys.stdout.write(f.readline().strip() + "\tPangoLineage\n")
        c = csv.reader(f, delimiter='\t')
        for row in c:
            smp = row[0]
            call = calls[smp] if smp in calls else "N/A"
            sys.stdout.write("\t".join(row) + "\t" + call + "\n")

if __name__ == "__main__":
    calls = readLineages(sys.argv[1])
    addLineage(sys.argv[2], calls)
