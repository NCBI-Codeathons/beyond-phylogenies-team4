#!/usr/bin/env python

import sys
import csv

import flacodb as fdb

"""
taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,version,pangolin_version,pangoLEARN_version,pango_version,status,note
hCoV-19/USA/FL-Miami-500/2021|Miami-500|2021-05-06|NorthAmerica,
P.1,
0.0,
0.9960763399489543,
Gamma (P.1-like),
0.875000,
0.000000,
PLEARN-v1.2.12,3.0.5,
2021-06-05,
v1.2.12,
passed_qc,
scorpio call: Alt alleles 14; Ref alleles 0; Amb alleles 1
"""

FIELDS = {'lineage': -1,
          'ambiguity_score': -1,
          'scorpio_call': -1,
          'pangolin_version': -1,
          'pangoLEARN_version': -1,
          'status': -1,
          'note': -1}

def findFields(hdr):
    idx = 0
    for f in hdr:
        if f in FIELDS:
            FIELDS[f] = idx
        idx += 1

def main(run, filename):
    DB = fdb.FLACOdb()
    DB.opendb()
    try:
        with open(filename, "r") as f:
            hdr = f.readline().rstrip("\r\n").split(",")
            findFields(hdr)
            c = csv.reader(f, delimiter=',')
            for line in c:
                if "|" in line[0]:
                    parts = line[0].split("|")
                    sample = parts[1]
                    ascore = line[FIELDS['ambiguity_score']]
                    if ascore:
                        prob = 1.0 - float(ascore)
                    else:
                        prob = 0.0
                    data = [ ("TPangoLineage", line[FIELDS['lineage']]),
                             ("FPangoProb", prob),
                             ("TScorpioCall", line[FIELDS['scorpio_call']]),
                             ("TPangoVersion", line[FIELDS['pangoLEARN_version']]),
                             ("TPangoStatus", line[FIELDS['status']]),
                             ("TPangoNote", line[FIELDS['note']]) ]
                    DB.storeValues(sample, run, data)
    finally:
        DB.closedb()
        
if __name__ == "__main__":
    args = sys.argv[1:]
    main(args[0], args[1])
