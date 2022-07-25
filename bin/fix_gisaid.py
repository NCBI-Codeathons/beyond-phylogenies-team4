#!/usr/bin/env python

import sys
import csv

def main(filename):
    with open(filename, "r") as f:
        c = csv.reader(f, delimiter='\t')
        for row in c:
            name = row[0]
            gisaid = row[1]
            if "Miami" in name:
                newname= name.replace("Miami", "MSL")
            elif "Saliva" in name:
                newname= name.replace("Saliva", "SAL")
            elif "saliva" in name:
                newname= name.replace("saliva", "SAL")
            else:
                newname = name
            sys.stdout.write(f"{gisaid}\t{name}\t{newname}\n")

if __name__ == "__main__":
    main(sys.argv[1])
