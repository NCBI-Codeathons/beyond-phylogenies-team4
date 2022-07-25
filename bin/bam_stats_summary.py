#!/usr/bin/env python

import sys
import flacodb as fdb

def perc(x, t):
    if t == 0:
        return "\t{}\t0%".format(x)
    else:
        return "\t{}\t{:.2f}%".format(x, 100.0*x/t)

def main(short=False):
    totalreads = 0
    counts = {'hs': 0, 'EBV': 0}

    for line in sys.stdin:
        r = line.split("\t")
        ref = r[0]
        nreads = int(r[2])
        totalreads += nreads

        if "decoy" in ref:
            pass
        elif ref == "*":
            pass
        elif ref == "chrEBV":
            counts['EBV'] += nreads
        elif ref[:3] == 'chr':
            counts['hs'] += nreads
        else:
            counts[ref] = nreads

    if short:
        sys.stdout.write("{}:{}:{}\n".format(totalreads, counts['hs'], counts['MN908947.3'] + counts['NC_004718.3']))
    else:
        sys.stdout.write("""#Group\tReads\tPctReads
Total reads\t{}\t
Human{}
SARS-COV-2{}
EBV{}
""".format(totalreads, perc(counts['hs'], totalreads), perc(counts['MN908947.3'] + counts['NC_004718.3'], totalreads), perc(counts['EBV'], totalreads)))
        for k in sorted(counts.keys()):
            if k in ['hs', 'EBV', 'MN908947.3', 'NC_004718.3']:
                pass
            else:
                sys.stdout.write("{}{}\n".format(k, perc(counts[k], totalreads)))

if len(sys.argv) == 1:
    main()
else:
    main(True)

