#!/usr/bin/env python3

import sys

from datetime import date, timedelta

from FlacoUtils import dateDays

class GisaidFilt(object):
    filename = ""
    outfile = "/dev/stdout"
    nameWanted = None
    nameNotWanted = None
    wantedIDs = []
    unwantedIDs = []
    fromDate = None
    toDate = None

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == "-o":
                self.outfile = a
                prev = ""
            elif prev == "-n":
                self.nameWanted = a
                prev = ""
            elif prev == "-N":
                self.nameNotWanted = a
                prev = ""
            elif prev == "-a":
                self.fromDate = dateDays(a)
                prev = ""
            elif prev == "-b":
                self.toDate = dateDays(a)
                prev = ""
            elif prev == "-w":
                with open(a, "r") as f:
                    self.wantedIDs = f.read().split("\n")
                prev = ""
            elif prev == "-W":
                with open(a, "r") as f:
                    self.unwantedIDs = f.read().split("\n")
                prev = ""
            elif a in ["-o", "-n", "-N", "-a", "-b", "-w", "-W"]:
                prev = a
            else:
                self.filename = a
        return self.filename

    def usage(self):
        sys.stdout.write("""gisaidfilt.py - Filter GISAID fasta files

Usage: gisaidfilt.py [options] filename.fa

Options:

  -o O | Write output to file O (default: standard output)
  -n N | Output sequences whose name contains N
  -N N | Output sequences whose name does not contain N
  -w W | Output sequences whose ID is contained in file W
  -W W | Output sequences whose name is not contained in file W
  -a A | Do not output sequences before date A
  -b B | Do not output sequences after date B

Dates for the -a and -b options should be specified in ISO format (YYYY-MM-DD) optionally
followed by +N or _N to specify a number of days before or after that date. For example,
to only return sequences within 28 days of Jan 1st 2021, use:

  -a 2021-01-01_28 -b 2021-01-01+28

""")

    def run(self):
        good = False
        self.badDates = 0
        nin = nout = 0

        if self.nameWanted:
            sys.stderr.write("Removing sequences not containing `{}'\n".format(self.nameWanted))
        if self.nameNotWanted:
            sys.stderr.write("Removing sequences containing `{}'\n".format(self.nameNotWanted))
        if self.fromDate:
            sys.stderr.write("Removing sequences before {}\n".format(date.isoformat(self.fromDate)))
        if self.toDate:
            sys.stderr.write("Removing sequences after {}\n".format(date.isoformat(self.toDate)))

        with open(self.outfile, "w") as out:
            with open(self.filename, "r") as f:
                for line in f:
                    if line[0] == '>':
                        nin += 1
                        good = self.checkHeader(line.rstrip("\r\n"))
                        if good:
                            out.write(line)
                            nout += 1
                        sys.stderr.write("\r{} / {} written.".format(nout, nin))
                    elif good:
                        out.write(line)
        sys.stderr.write("\r{} / {} written.\n".format(nout, nin))
        if self.badDates:
            sys.stderr.write("{} bad dates ignored.\n".format(self.badDates))

    def checkHeader(self, hdr):
        parts = hdr.split("|")
        name = parts[0]
        seqid = parts[1]
        if self.nameWanted and self.nameWanted not in name:
            return False
        if self.nameNotWanted and self.nameNotWanted in name:
            return False
        if self.wantedIDs and seqid not in self.wantedIDs:
            return False
        if self.unwantedIDs and seqid in self.unwantedIDs:
            return False
        try:
            thisdate = date.fromisoformat(parts[2])
        except ValueError:
            #sys.stderr.write("Bad date {} - ignoring.\n".format(parts[2]))
            self.badDates += 1
            return False
        if self.fromDate and thisdate < self.fromDate:
            return False
        if self.toDate and thisdate > self.toDate:
            return False
        return True

if __name__ == "__main__":
    args = sys.argv[1:]
    G = GisaidFilt()
    if "-h" in args:
        G.usage()
    elif G.parseArgs(args):
        G.run()
    else:
        G.usage()
