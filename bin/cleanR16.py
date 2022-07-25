#!/usr/bin/env python

import sys, csv

__doc__ = """Remove samples from R16 that were redone in R18."""

def readToRemove(filename):
    with open(filename, "r") as f:
        return f.read().strip().split("\n")

def remove(filename, toremove):
    with open(filename, "r") as f:
        c = csv.reader(f, delimiter='\t')
        for row in c:
            if row[0] not in toremove:
                sys.stdout.write("\t".join(row) + "\n")

if __name__ == "__main__":
    remove(sys.argv[1], readToRemove(sys.argv[2]))
