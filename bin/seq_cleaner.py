#!/usr/bin/env python

import sys
import csv
from Bio import Seq, SeqIO

SPIKE = [21563, 25384]
SPIKELEN = 3821

class Cleaner(object):
    infile = None
    outfile = "/dev/stdout"
    badfile = None
    badstream = None
    reportfile = None
    min_length = 0
    n_frac = 0
    n_stretch = None
    remgaps = False
    wantedfile = None
    wanted = None
    musthavedate = False
    partialmatch = []
    exclude = False             # If True, strings in partialmatch are excluded
    remdup = False

    def parseArgs(self, args):
        if "-h" in args or "--help" in args:
            return False
        self.partialmatch = []
        prev = ""
        for a in args:
            if prev == "-l":
                self.min_length = int(a)
                prev = ""
            elif prev == "-f":
                self.n_frac = float(a)
                prev = ""
            elif prev == "-b":
                self.badfile = a
                prev = ""
            elif prev == "-r":
                self.reportfile = a
                prev = ""
            elif prev == "-n":
                self.n_stretch = int(a)
                prev = ""
            elif prev == "-w":
                self.wantedfile = a
                prev = ""
            elif prev == "-p":
                self.partialmatch = [a]
                prev = ""
            elif prev == "-x":
                self.partialmatch = [a]
                self.exclude = True
                prev = ""
            elif prev == "-P":
                with open(a, "r") as f:
                    for line in f:
                        self.partialmatch.append(line.strip())
                prev = ""
            elif a in ["-l", "-f", "-b", "-r", "-n", "-w", "-p", "-P", "-x"]:
                prev = a
            elif a == "-g":
                self.remgaps = True
            elif a == "-d":
                self.musthavedate = True
            elif a == "-u":
                self.remdup = True
            elif self.infile is None:
                self.infile = a
            else:
                self.outfile = a
        return self.infile

    def usage(self):
        sys.stdout.write("""Usage: seq_cleaner.py [options] input.fasta [output.fasta]

Clean sequences in the input FASTA file on the basis of minimum length, and maximum N frequency.

Options:

  -l L | Only output sequences with length greater then or equal to L.
  -f F | Only output sequences with a fraction of N smaller than F.
  -w W | Only output sequences whose name is listed in file W.
  -p P | Output sequences that contain string P anywhere in the header.
  -x X | Output sequences that do NOT contain string X in the header.
  -P P | Like -p, but reads strings to match from file P.
  -d   | Only output sequences with a valid date in the header.
  -u   | Remove duplicate sequences (base on seq name, only keeps first).
  -g   | Remove all gaps from sequences.
  -b B | Write "bad" sequences to file B.
  -r R | Write report to file R.

""")
        return False

    def writeBad(self, seq_record):
        if self.badstream:
            SeqIO.write(seq_record, self.badstream, "fasta")

    def readWanted(self):
        if self.wantedfile:
            self.wanted = set()
            with open(self.wantedfile, "r") as f:
                c = csv.reader(f, delimiter='\t')
                for row in c:
                    if row and row[0][0] != '#':
                        self.wanted.add(row[0])
            sys.stderr.write("{} wanted sequences.\n".format(len(self.wanted)))

    def clean(self):
        nin = nout = nshort = nmissing = 0
        seqdata = []
        seen = set()

        if self.badfile:
            self.badstream = open(self.badfile, "w")
        self.readWanted()
        try:
            with open(self.outfile, "w") as out:
                for seq_record in SeqIO.parse(self.infile, "fasta"):
                    nin += 1
                    header = seq_record.name
                    if self.remdup and header in seen:
                        continue
                    else:
                        seen.add(header)
                    if self.partialmatch:
                        if self.exclude:
                            good = True
                            for p in self.partialmatch:
                                if p in header:
                                    good = False
                                    break
                        else:
                            good = False
                            for p in self.partialmatch:
                                if p in header:
                                    good = True
                                    break
                        if not good:
                            continue
                    header_fields = header.split("|")
                    if len(header_fields) > 2:
                        name = header_fields[1]
                        if self.wanted and name not in self.wanted:
                            continue
                        date = header_fields[2]
                        if self.musthavedate and date == '':
                            continue
                    seq = str(seq_record.seq).upper()
                    if self.remgaps:
                        seq = seq.replace("-", "")
                        seq_record.seq = Seq.Seq(seq)
                    sl = len(seq)
                    if sl == 0:
                        sys.stderr.write("Warning: sequence {} is empty.\n".format(seq_record.name))
                        continue
                    nc = seq.count("N") + seq.count("-")
                    nspike = seq[SPIKE[0]:SPIKE[1]].count("N") + seq[SPIKE[0]:SPIKE[1]].count("-")
                    nf = 100 - 100.0*nc/sl
                    sf = 100 - 100.0*nspike/SPIKELEN
                    status = "NO"

                    if sl < self.min_length:
                        nshort += 1
                        self.writeBad(seq_record)
                    elif nf < self.n_frac:
                        nmissing += 1
                        self.writeBad(seq_record)
                    else:
                        if self.n_stretch:
                            self.remove_nstretch(seq_record)
                        SeqIO.write(seq_record, out, "fasta")
                        nout += 1
                        status = "OK"
                    seqdata.append([seq_record.name, sl, "{:.1f}".format(nf), "{:.1f}".format(sf), status])

        finally:
            if self.badstream:
                self.badstream.close()

        if self.reportfile:
            with open(self.reportfile, "w") as rep:
                rep.write("#Input sequences:\t{}\n".format(nin))
                rep.write("#Too short:\t{}\n".format(nshort))
                rep.write("#More than {}% Ns:\t{}\n".format(int(100 * self.n_frac), nmissing))
                rep.write("#Good sequences:\t{}\n".format(nout))
                rep.write("#\n")
                rep.write("#Sequence\tLen\tNfrac\tSpike\tStatus\n")
                for row in seqdata:
                    rep.write("{}\t{}\t{}\t{}\t{}\n".format(*row))
        sys.stderr.write("{} sequences read, {} written ({} filtered).\n".format(nin, nout, nin-nout))
            
    def remove_nstretch(self, seq_record):
        cleanbases = []
        instretch = False
        lenstretch = 0
        for b in seq_record.seq:
            if instretch:
                if b == "N":
                    lenstretch += 1
                else:
                    if lenstretch <= self.n_stretch:
                        cleanbases.append("N"*lenstretch)
                    lenstretch = 0
                    instretch = False
            else:
                if b == "N":
                    instretch = True
                    lenstretch = 1
                else:
                    cleanbases.append(b)
        seq_record.seq = Seq.Seq("".join(cleanbases))

if __name__ == "__main__":
    args = sys.argv[1:]
    C = Cleaner()
    if C.parseArgs(args):
        C.clean()
    else:
        C.usage()
