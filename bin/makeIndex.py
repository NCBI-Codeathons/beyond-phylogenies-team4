#!/usr/bin/env python

import sys
from os import getenv
from collections import defaultdict
import flacodb as fdb

def samplesByMonth(F):
    counts = defaultdict(int)
    for row in F.execute("SELECT SamplingDate FROM Samples WHERE SamplingDate <> '';").fetchall():
        m = row[0][:7]
        counts[m] += 1
    months = sorted(counts.keys())
    return [ [m, counts[m]] for m in months]

def samplesByGroup(F):
    counts = defaultdict(int)
    for row in F.execute("SELECT SampleName FROM Samples;").fetchall():
        m = row[0].split("-")[0]
        counts[m] += 1
    nvax = int(F.execute("SELECT count(*) FROM Samples WHERE Vaccinated='Y';").fetchone()[0])
    counts["Vaccinated"] = nvax

    # Ad-hoc fixes to group counts
    counts["Miami"] += counts["MSL"]
    for k in counts.keys():
        if k.startswith("Haiti") and k != "Haiti":
            counts["Haiti"] += counts[k]
    return counts

def writeStats(F, runcounts):
    totsamples = sum(runcounts.values())
    months = samplesByMonth(F)
    n1 = int(F.execute("select count(*) from sampledata where Field='FCovPct' and Fval>20;").fetchone()[0])
    n2 = int(F.execute("select count(*) from sampledata where Field='TPangoLineage';").fetchone()[0])
    sys.stdout.write("""<UL Class="w3-ul w3-border">
<LI><B>Runs:</B> {}</LI>
<LI><B>Total samples:</B> {:,}</LI>
<LI><B>Samples with SARS-CoV-2 reads:</B> {:,} ({}%)</LI>
<LI><B>Samples with lineage call:</B> {:,} ({}%)</LI>
<LI><B>Samples by month:</B><PRE>
""".format(len(runcounts), totsamples, n1, int(100.0*n1/totsamples), n2, int(100.0*n2/n1)))
    for mo in months:
        sys.stdout.write("&nbsp;&nbsp;{} {:6}\n".format(mo[0], mo[1]))
    sys.stdout.write("</PRE></UL>")

def runstable(runsfile, runcounts):
    with open(runsfile, "r") as f:
        runs = f.read().strip().split("\n")
    nr = 0
    for run in runs:
        desc = "" if nr == 0 else "Run #{}".format(nr)
        sys.stdout.write("""	  <tr>
	    <td>{0}</td>
            <TD class='w3-right-align'>{2}</Td>
	    <TD><A href='{0}-stats/{0}-stats.html'>View</A> / <A href="{0}-stats/{0}-stats.txt" download>Data</A></TD>
	    <TD><A href='{0}-stats/{0}-cov.html'>View</A> / <A href="{0}-stats/{0}-cov.txt" download>Data</A></TD>
	    <TD><I>{1}</I></TD>
""".format(run, desc, runcounts[run]))
        nr += 1

def runlineages(runsfile, runcounts):
    with open(runsfile, "r") as f:
        runs = f.read().strip().split("\n")
    nr = 0
    for run in runs:
        sys.stdout.write("""	  <tr>
	    <td>{0}</td>
            <TD class='w3-right-align'>{1}</Td>
	    <TD><A href='{0}-stats/{0}-results.html'>View</A> / <A href="{0}-stats/{0}-results.txt" download>Data</A></TD>
	    <TD><A href='{0}-stats/{0}-lintable.html'>View</A> / <A href="{0}-stats/{0}-lintable.txt" download>Data</A></TD>
""".format(run, runcounts[run]))
        nr += 1

def groupstable(groupsfile, groupcounts):
    with open(groupsfile, "r") as f:
        groups = f.read().strip().split("\n")
    for grp in groups:
        grp = grp.split(",")[0]
        sys.stdout.write("""	  <tr>
	    <td>{0}</td>
            <td class='w3-right-align'>{1}</td>
	    <TD><A href='group-results/{0}-stats.html'>View</A> / <A href="group-results/{0}-stats.txt" download>Data</A></TD>
	    <TD><A href='group-results/{0}-cov.html'>View</A> / <A href="group-results/{0}-cov.txt" download>Data</A></TD>
	    <TD><I></I></TD>
""".format(grp, groupcounts[grp]))

def grouplineages(groupsfile, groupcounts):
    with open(groupsfile, "r") as f:
        groups = f.read().strip().split("\n")
    nr = 0
    for group in groups:
        group = group.split(",")[0]
        sys.stdout.write("""	  <tr>
	    <td>{0}</td>
            <TD class='w3-right-align'>{1}</Td>
	    <TD><A href='group-results/{0}-results.html'>View</A> / <A href="group-results/{0}-results.txt" download>Data</A></TD>
	    <TD><A href='group-results/{0}-lintable.html'>View</A> / <A href="group-results/{0}-lintable.txt" download>Data</A></TD>
""".format(group, groupcounts[group]))
        nr += 1

def main(template, runsfile, groupsfile):
    F = fdb.FLACOdb()
    F.opendb()
    try:
        runcounts = F.getRunCounts()
        groupcounts = samplesByGroup(F)
        with open(template, "r") as f:
            for line in f:
                if line.startswith("##Runs"):
                    runstable(runsfile, runcounts)
                elif line.startswith("##Groups"):
                    groupstable(groupsfile, groupcounts)
                elif line.startswith("##Stats"):
                    writeStats(F, runcounts)
                elif line.startswith("##LinRuns"):
                    runlineages(runsfile, runcounts)
                elif line.startswith("##LinGroups"):
                    grouplineages(groupsfile, groupcounts)
                else:
                    sys.stdout.write(line)
    finally:
        F.closedb()

if __name__ == "__main__":
    args = sys.argv[1:]
    if args:
        template = getenv("FLACO_TEMPLATES") + "/indexTemplate.html"
        main(template, args[0], args[1])
    else:
        sys.stdout.write("Usage: makeIndex.py RUNS GROUPS\n")
