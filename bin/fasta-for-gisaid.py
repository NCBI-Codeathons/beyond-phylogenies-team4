#!/usr/bin/env python

import sys

def main(infile, outfile):
    with open(outfile, "w") as out, open(infile, "r") as f:
        for line in f:
            if line[0] == '>':
                parts = line.split("|")
                out.write(parts[0] + "\n")
                sys.stdout.write(parts[1] + "\n")
            else:
                out.write(line)

if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) >= 2:
        infile = args[0]
        outfile = args[1]
        if infile == "-":
            infile = "/dev/stdin"
        if outfile == "-":
            outfile = "/dev/stdout"
        main(infile, outfile)
    else:
        sys.stderr.write("Usage: fasta-for-gisaid.py input.fa output.fa\n")
