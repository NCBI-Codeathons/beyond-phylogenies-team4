#!/usr/bin/env python3

from os import getenv
import sys
import csv
import os.path
import sqlite3 as sql
import FlacoUtils as fu

TABLES = [
    """DROP TABLE Samples;""",
    """CREATE TABLE Samples (
    TubeID int,
    SampleNumber int,
    SampleName text primary key,
    Patient text,
    SampleType text,
    Timepoint int,
    DaysAfterBaseline int,
    Location text,
    GenomeStatus text,
    Lineage text,
    Notes text,
    DiagnosisLab text,
    Vaccinated char(1) default 'N',
    VaccineAdministered text,
    SamplingDate text,
    City text,
    County text,
    ZipCode text,
    Sex char(1),
    Age int,
    DiseaseOutcome text,
    CoMorbidityFactors text,
    Travel text,
    ViralLoad float,
    LibraryPrep text,
    CT float
);""",

  """DROP TABLE SampleData;""",
  """CREATE TABLE SampleData (
    SampleName text,
    Run text,
    Field text,
    Ival int,
    Fval float,
    Tval text);""",
    """CREATE INDEX SampleData_Sample on SampleData(SampleName);""",
    """CREATE INDEX SampleData_Run on SampleData(Run);""",
    """CREATE INDEX SampleData_Field on SampleData(Field);""",

  """DROP TABLE SampleGISAID;""",
  """CREATE TABLE SampleGISAID (
     SampleName text,
     GISAID text);""",
  """CREATE INDEX SampleGISAID_Sample on SampleGISAID(SampleName);""",
  """CREATE INDEX SampleGISAID_GISAID on SampleGISAID(GISAID);""",

  """DROP TABLE IF EXISTS Runs;""",
  """CREATE TABLE Runs (
     Name text,
     RunDate text,
     RunNumber int,
     Flowcell text);""",
    ]

METADATA_FIELDS = ["Tube ID", "Sample Number ICBR", "Sample Name",
                   "Patient", "Sample Type", "Timepoint (longitudinal)",
                   "Days after baseline", "Location", "Genome status",
                   "Lineage", "Notes", "Lab that Diagnosed Infection",
                   "Vaccinated? (Y/N)", "Vaccine administered, if any",
                   "Sampling Date", "city of residence", "county of residence", "Patient zip code",
                   "Sex", "Age", "disease outcome (if known)", "co-morbidity factors (if known)",
                   "history of travel in the previous two weeks before diagnosis (if known)",
                   "Viral load or Information", "Manual or Automated Library prep?"]

QUERABLE_METADATA = ["SamplingDate", "Vaccinated", "Sex", "Age", "ViralLoad", "CT"]


FIELDSMAP = {"Tube ID": "TubeID",
             "Sample Number ICBR": "SampleNumber",
             "Sample Name": "SampleName",
             "Patient": "Patient",
             "Sample Type": "SampleType",
             "Timepoint (longitudinal)": "Timepoint",
             "Days after baseline": "DaysAfterBaseline",
             "Location": "Location",
             "Genome status": "GenomeStatus",
             "Lineage": "Lineage",
             "Notes": "Notes",
             "Lab that Diagnosed Infection": "DiagnosisLab",
             "Vaccinated? (Y/N)": "Vaccinated",
             "Vaccine administered, if any": "VaccineAdministered",
             "Sampling Date": "SamplingDate",
             "city of residence": "City",
             "county of residence": "County",
             "Patient zip code": "ZipCode",
             "Sex": "Sex",
             "Age": "Age",
             "disease outcome (if known)": "DiseaseOutcome",
             "co-morbidity factors (if known)": "CoMorbidityFactors",
             "history of travel in the previous two weeks before diagnosis (if known)": "Travel",
             "Viral load or Information": "ViralLoad",
             "Manual or Automated Library prep?": "LibraryPrep"}

### Query language

def smaller(a, b):
    return a < b

def larger(a, b):
    return a > b

def smallerq(a, b):
    return a <= b

def largerq(a, b):
    return a >= b

def equal(a, b):
    return a == b

def different(a, b):
    return a != b

def contains(a, b):
    return b.lower() in a.lower()

CLAUSEMAP = [[">=", largerq],
             ["<=", smallerq],
             ["==", equal],
             ["=", equal],
             ["!=", different],
             ["<", smaller],
             [">", larger],
             ["%", contains]]

def safeFloat(f):
    try:
        return float(f)
    except ValueError:
        return f

class QueryTerm(object):
    """An object that represents clauses of the form 'variable op value', where
op can be one of: >, <, >=, <=, ==, !=."""

    varname = ""
    value = None
    operator = None
    required = False

    def __init__(self, clause):
        if clause[0] == '+':
            self.required = True
            clause = clause[1:]
        for cm in CLAUSEMAP:
            p = clause.find(cm[0])
            if p > 0:
                left = clause[:p]
                right = safeFloat(clause[p+len(cm[0]):])
                self.varname = left.strip()
                self.value = right
                self.operator = cm[1]
                break

        # If there's no operator, assume we simply want to
        # to display the var's value
        if not self.operator:
            self.varname = clause
            self.value = "***nosuchvalue***"
            self.operator = different

    def satisfied(self, datum):
        if datum:
            if self.operator:
                return self.operator(datum, self.value)
            else:
                return True

        # All tests return False on missing data
        return False

### Main Class

class FLACOdb(object):
    dbfile = None
    outfile = "/dev/stdout"
    wantedfile = None
    fmt = "csv"                 # or "html"
    args = []
    cmdopts = []                # Command options
    wanted = []                 # List of wanted samples
    wantedruns = []             # Only retrieve samples in specified run(s)
    wantedgroups = []           # Only retrieve samples in specified group(s)
    unwantedgroups = []
    vax = False                 # True if we want vaccinated only
    novax = False               # True if we want unvaccinated only
    fromDate = None             # Minimum sampling date
    toDate = None               # Maximum sampling date
    multiRun = False            # Output data for all runs a sample belongs to

    _conn = None                # DB connection
    _lvl = 0                    # For nested opendb() calls

    def __init__(self):
        self.dbfile = getenv("FLACO_DB")
        self.args = []
        
    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == "-db":
                self.dbfile = a
                prev = ""
            elif prev == "-o":
                self.outfile = a
                prev = ""
            elif prev == "-f":
                if a in ["csv", "html", "w3"]:
                    self.fmt = a
                prev = ""
            elif prev == "-g":
                self.wantedgroups = a.split(",")
                prev = ""
            elif prev == "-x":
                self.unwantedgroups = a.split(",")
                prev = ""
            elif prev == "-w":
                for name in a.split(","):
                    self.wanted.append(name)
                prev = ""
            elif prev == "-W":
                self.wantedfile = a
                prev = ""
            elif prev == "-r":
                wr = a.split(",")
                for r in wr:
                    self.wantedruns.append(r)
                prev = ""
            elif prev == "-a":
                self.fromDate = str(fu.dateDays(a))
                prev = ""
            elif prev == "-b":
                self.toDate = str(fu.dateDays(a))
                prev = ""
            elif a in ["-o", "-db", "-f", "-w", "-W", "-r", "-a", "-b", "-g", "-x"]:
                prev = a
            elif a == "-vax":
                self.vax = True
            elif a == "-novax":
                self.novax = True
            elif a == "-m":
                self.multiRun = True
            elif a[0] == "-":
                self.cmdopts.append(a)
            else:
                self.args.append(a)
        return self.args

    def parseGroups(self, a):
        grps = a.split(",")
        for grp in grps:
            if grp[0] == '-':
                self.unwantedgroups.append(grp[1:])
            else:
                self.wantedgroups.append(grp)

    def opendb(self):
        if self._conn is None:
            self._conn = sql.connect(self.dbfile)
            self._conn.row_factory = sql.Row
        self._lvl += 1

    def closedb(self):
        self._lvl += -1
        if self._lvl == 0:
            self._conn.commit()
            self._conn.close()
            self._conn = None

    def execute(self, query, *args):
        return self._conn.execute(query, args)
        
    def execute1(self, query, args):
        return self._conn.execute(query, args)
        
    def initialize(self):
        self.opendb()
        for tab in TABLES:
            self.execute(tab)
        self.closedb()

    ## Data operations

    def getFilteredSamples(self):
        """Returns a list of sample names using the constraints in wanted, vax/novax, from/toDate."""
        names = []
        result = []

        self.opendb()
        try:
            if self.wanted:
                names = self.wanted
            elif self.wantedruns:
                for r in self.wantedruns:
                    names += self.getRunSamples(r)
                    
            for row in self.execute("SELECT SampleName, Vaccinated, SamplingDate FROM Samples;"):
                if names and row[0] not in names:
                    continue
                if self.vax and row[1] != 'Y':
                    continue
                if self.novax and row[1] != 'N':
                    continue
                if self.fromDate and (row[2] is None or row[2] < self.fromDate):
                    continue
                if self.toDate and (row[2] is None or row[2] > self.toDate):
                    continue
                grp = row[0].split("-")[0]
                if grp in self.unwantedgroups:
                    continue
                if self.wantedgroups:
                    if grp not in self.wantedgroups:
                        continue
                result.append(row[0])
        finally:
            self.closedb()
        return result
            
    def getValue(self, sample, run, field):
        """Return the value of `field' for `sample' in `run'."""
        tag = field[0]
        query = "SELECT {}val FROM SampleData WHERE SampleName=? and Run=? and Field=?;".format(tag)
        return self.execute(query, sample, run, field).fetchone()

    def storeValue(self, sample, run, field, value):
        """Store or update the value of `field' for `sample' in `run', setting it to `value'."""
        tag = field[0]
        row = self.getValue(sample, run, field)
        if row:
            query = "UPDATE SampleData SET {}val=? WHERE SampleName=? and Run=? and Field=?;".format(tag)
            self.execute(query, value, sample, run, field)
        else:
            query = "INSERT INTO SampleData(SampleName, Run, Field, {}val) VALUES (?, ?, ?, ?);".format(tag)
            self.execute(query, sample, run, field, value)

    def storeValues(self, sample, run, pairs):
        """Store multiple values at once for `sample' in `run'."""
        self.opendb()
        try:
            if self.checkSampleKnown(sample):
                for pair in pairs:
                    self.storeValue(sample, run, pair[0], pair[1])
            else:
                sys.stderr.write("Warning: unknown sample `{}'.\n".format(sample))
        finally:
            self.closedb()

    def test(self):
        self.opendb()
        self.storeValue("Smp1", "Run1", "Iuno", 10)
        self.storeValue("Smp1", "Run1", "Tdue", "abc")
        self.storeValue("Smp1", "Run1", "Iuno", 20)
        self.storeValue("Smp2", "Run1", "Iuno", 30)
        self.storeValue("Smp1", "Run2", "Iuno", 40)
        self.storeValue("Smp4", "Run2", "Iuno", 50)
        self.storeValue("Smp5", "Run2", "Iuno", 60)
        print(self.checkSampleKnown("Path-Saliva-100"))
        print(self.checkSampleKnown("Smp10"))
        self.closedb()

    def checkSampleKnown(self, sample):
        self.opendb()
        try:
            return self.execute("SELECT SampleName FROM Samples WHERE SampleName=?", sample).fetchone()
        finally:
            self.closedb()
    
    def getSampleRuns(self, sample):
        """Return all the runs that this sample appears in."""
        runs = []
        self.opendb()
        try:
            for row in self.execute("SELECT DISTINCT Run FROM SampleData WHERE SampleName=?", sample):
                runs.append(row[0])
        finally:
            self.closedb()
        return runs

    def getRunSamples(self, run):
        """Return all the samples in this run."""
        samples = []
        self.opendb()
        try:
            for row in self.execute("SELECT DISTINCT SampleName FROM SampleData WHERE Run=?", run):
                samples.append(row[0])
        finally:
            self.closedb()
        return samples

    def getRunCounts(self):
        """Return a dictionary of run names and number of samples they contain."""
        runcounts = {}
        self.opendb()
        try:
            for row in self.execute("SELECT DISTINCT Run FROM SampleData;"):
                run = row[0]
                cnt = self.execute("SELECT COUNT(DISTINCT SampleName) FROM SampleData WHERE Run=?", run).fetchone()
                runcounts[run] = int(cnt[0])
        finally:
            self.closedb()
        return runcounts

    def getSampleData(self, sample, run):
        """Return everything known about `sample' in `run' as a dictionary."""
        data = {'SampleName': sample, 'Run': run}
        metadata_query = "SELECT {} FROM Samples WHERE SampleName=?".format(",".join(QUERABLE_METADATA))
        self.opendb()
        try:
            mdrow = self.execute(metadata_query, sample).fetchone()
            for qf in QUERABLE_METADATA:
                data[qf] = mdrow[qf]
            for row in self.execute("SELECT Field, Ival, Fval, Tval FROM SampleData WHERE SampleName=? and Run=?;", sample, run):
                tag = row[0][0]
                field = row[0][1:]
                if tag == "I":
                    data[field] = int(row[1])
                elif tag == "F":
                    data[field] = float(row[2])
                elif tag == "T":
                    data[field] = row[3]
        finally:
            self.closedb()
        return data

    def getSampleLineage(self, sample, run):
        self.opendb()
        try:
            return self.execute("""SELECT a.SampleName, a.SamplingDate, b.Tval 
FROM samples a, sampledata b 
WHERE a.samplename=? AND b.run=? AND a.samplename = b.samplename AND b.Field='TPangoLineage' AND b.Tval != ''
ORDER BY a.SamplingDate DESC;""", sample, run).fetchone()
        finally:
            self.closedb()

    ## Import

    def importSamples(self):
        filename = args[1]
        fieldNames = []
        positions = []
        dateCol = None
        
        sys.stderr.write("*** Parsing file {}\n".format(filename))
        with open(filename, "r", encoding="iso-8859-1") as f:
            c = csv.reader(f, delimiter='\t')
            hdr = c.__next__()
            pos = 0
            for colname in hdr:
                if colname in FIELDSMAP:
                    fieldNames.append(FIELDSMAP[colname])
                    positions.append(pos)
                    if FIELDSMAP[colname] == "SamplingDate":
                        dateCol = pos
                pos += 1
            sys.stderr.write("*** {} valid columns detected\n".format(len(positions)))
            query = "INSERT OR REPLACE INTO Samples ({}) VALUES ({});".format(", ".join(fieldNames), ", ".join(["?"]*len(positions)))
            nrows = 0
            self.opendb()
            try:
                for row in c:
                    if dateCol:
                        row[dateCol] = fu.formatDate(row[dateCol])
                    rvals = [ row[idx] for idx in positions ]
                    self.execute1(query, rvals)
                    nrows += 1
            finally:
                self.closedb()
            sys.stderr.write("*** {} rows imported.\n".format(nrows))

    def getSamples(self):
        """Return all samples."""
        samples = self.getFilteredSamples()
        with open(self.outfile, "w") as out:
            for smp in samples:
                out.write(smp + "\n")
        sys.stderr.write("*** {} sample names written.\n".format(len(samples)))

    def getRunData(self):
        """Show all runs, the samples belonging to one or more runs."""
        self.opendb()
        try:
            with open(self.outfile, "w") as out:
                if len(self.args) == 1: # show all runs?
                    counts = self.getRunCounts()
                    out.write("#RunNumber\tRun\tDate\tFlowcell\tNumSamples\n")
                    for row in self.execute("""SELECT RunNumber, Name, RunDate, Flowcell FROM Runs ORDER BY RunDate;""").fetchall():
                        out.write("{}\t{}\t{}\t{}\t{}\n".format(row["RunNumber"], row["Name"], row["RunDate"], row["Flowcell"], counts[row["Name"]]))
                else:
                    out.write("#Run\tSample\n")
                    for run in self.args[1:]:
                        samples = self.getRunSamples(run)
                        for smp in samples:
                            out.write("{}\t{}\n".format(run, smp))
        finally:
            self.closedb()

    def getMetadata(self):
        """Return metadata for all samples."""
        self.opendb()
        try:
            with open(self.outfile, "w") as out:
                FIELDS = [ FIELDSMAP[m] for m in METADATA_FIELDS ]
                if self.fmt == "html":
                    W = fu.TSVWriterHTML(out, FIELDS)
                    W.title = "Sample Metadata"
                    W.formats = [""]*len(FIELDS)
                else:
                    W = fu.TSVWriter(out, FIELDS)
                W.header()
                samplecol = FIELDS.index("SampleName")
                query = "SELECT {} FROM Samples ORDER BY TubeID;".format(", ".join(FIELDS))
                samplenames = self.getFilteredSamples()
                n = 0
                for row in self.execute(query).fetchall():
                    if row[samplecol] in samplenames:
                        W.writeRow(row)
                        n += 1
        finally:
            self.closedb()
        sys.stderr.write("*** {} metadata entries written.\n".format(n))

    def getSampleStats(self):
        """Print sample data for one or more samples."""

        FIELDS = ['Run', 'SampleName', 'SamplingDate', 'TotalReads', 'HsReads', 'CovReads', 'CovPct']
        FORMAT = ["S",   "I",          "S",            "D",          "D",       "D",        "P"]

        n = 0
        with open(self.outfile, "w") as out:
            if self.fmt == "html":
                W = fu.StatsWriterHTML(out, FIELDS)
                W.title = ", ".join(self.args[1:]) + " sample statistics"
                W.formats = FORMAT
            elif self.fmt == "w3":
                W = fu.StatsWriterW3(out, FIELDS, FORMAT)
                W.toggles = True
                W.template = getenv("FLACO_TEMPLATES") + "/statsTemplate.html"
                W.title = ", ".join(self.args[1:]) + " - sample statistics"
                W.formats = FORMAT
            else:
                W = fu.StatsWriter(out, FIELDS)
            W.header()
            self.opendb()
            try:
                for smp in self.getFilteredSamples():
                    for run in self.getSampleRuns(smp):
                        if not self.wantedruns or run in self.wantedruns:
                            data = self.getSampleData(smp, run)
                            fu.floatize(data, "CovPct")
                            W.writeRow(data)
                            n += 1

                W.close()
            finally:
                self.closedb()
        sys.stderr.write(""" *** Stats for {} samples written ({}).\n""".format(n, self.fmt))

    def getCoverageStats(self):
        """Print sample data for one or more samples."""

        FIELDS = ['Run', 'SampleName', 'SamplingDate', 'TotalCoverage', 'MedianDepth', 'AvgDepth', 'HiDepthFrac20', 'HiDepthFrac100', 'HiDepthFrac200', "MaxNStretch", "SpikeCov"]
        FORMAT = ["S",   "I",          "S",            "D",             "D",           "R",        "P",             "P",              "P",              "D",           "P"]

        n = 0
        with open(self.outfile, "w") as out:
            if self.fmt == "html":
                W = fu.StatsWriterHTML(out, FIELDS)
                W.title = ", ".join(self.args[1:]) + " coverage statistics"
                W.formats = FORMAT
            elif self.fmt == "w3":
                W = fu.StatsWriterW3(out, FIELDS, FORMAT)
                W.toggles = True
                W.template = getenv("FLACO_TEMPLATES") + "/statsTemplate.html"
                W.title = ", ".join(self.args[1:]) + " - coverage statistics"
                W.formats = FORMAT
            else:
                W = fu.StatsWriter(out, FIELDS)
            W.header()
            self.opendb()
            try:
                for smp in self.getFilteredSamples():
                    for run in self.getSampleRuns(smp):
                        if not self.wantedruns or run in self.wantedruns:
                            data = self.getSampleData(smp, run)
                            fu.floatize(data, 'AvgDepth')
                            fu.floatize(data, 'HiDepthFrac20')
                            fu.floatize(data, 'HiDepthFrac100')
                            fu.floatize(data, 'HiDepthFrac200')
                            fu.floatize(data, "SpikeCov")
                            W.writeRow(data)
                            n += 1

                W.close()
            finally:
                self.closedb()
        sys.stderr.write(""" *** Coverage stats for {} samples written ({}).\n""".format(n, self.fmt))

    def getSampleResults(self):
        
        FIELDS = ['Run', 'SampleName', 'SamplingDate', 'PangoLineage', 'ScorpioCall', 'PangoProb', 'PangoVersion', 'PangoNote']
        FORMAT = ["S",   "I",          "S",            "S",            "S",           "F",         "S",            ""]

        n = 0
        with open(self.outfile, "w") as out:
            if self.fmt == "html":
                W = fu.ResultsWriterHTML(out, FIELDS)
                W.title = ", ".join(self.args[1:]) + " - lineage results"
                W.formats = FORMAT
            elif self.fmt == "w3":
                W = fu.StatsWriterW3(out, FIELDS, FORMAT)
                W.template = getenv("FLACO_TEMPLATES") + "/resultsTemplate.html"
                W.title = ", ".join(self.args[1:]) + " lineage results"
                W.toggles = False
                W.formats = FORMAT

            else:
                W = fu.Writer(out, FIELDS)
            W.header()
            self.opendb()
            try:
                for smp in self.getFilteredSamples():
                    for run in self.getSampleRuns(smp):
                        if not self.wantedruns or run in self.wantedruns:
                            data = self.getSampleData(smp, run)
                            if 'PangoLineage' in data:
                                fu.floatize(data, "PangoProb")
                                W.writeRow(data)
                                n += 1
                W.close()
            finally:
                self.closedb()
        sys.stderr.write(""" *** Results for {} samples written ({}).\n""".format(n, self.fmt))

    def lineagesTable(self):
        lineages = list(set(fu.LINEAGES.values())) # need to remove duplicates
        FIELDS = ["Date", "N"] + lineages + ["Other"]
        FORMAT = ["S", "D"] + ["P"]*len(lineages) + ["P"] 

        n = 0
        with open(self.outfile, "w") as out:
            if self.fmt == "w3":
                W = fu.StatsWriterW3(out, FIELDS, FORMAT)
                W.template = getenv("FLACO_TEMPLATES") + "/linTemplate.html"
                W.title = ", ".join(self.args[1:]) + " lineages over time"
            else:
                W = fu.Writer(out, FIELDS)

            W.header()
            result = []
            self.opendb()
            try:
                for smp in self.getFilteredSamples():
                    for run in self.getSampleRuns(smp):
                        if not self.wantedruns or run in self.wantedruns:
                            lg = self.getSampleLineage(smp, run)
                            if lg and lg["SamplingDate"]:
                                result.append([fu.dateToDays(lg["SamplingDate"]), lg["SamplingDate"], lg["Tval"]])
            finally:
                self.closedb()
            if result:
                result.sort(key=lambda r: r[0], reverse=True)
                for row in fu.lineagesInWindow(result, lineages):
                    W.writeRow(row)
                    n += 1
            W.close()
        sys.stderr.write(""" *** Lineages for {} time periods written ({}).\n""".format(n, self.fmt))

    def readWanted(self):
        if not self.wantedfile:
            return
        if not os.path.isfile(self.wantedfile):
            sys.stderr.write("Warning: file {} does not exist.\n".format(self.wantedfile))
            return
        with open(self.wantedfile, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    self.wanted.append(line.split("\t")[0])
    
    def runQuery(self):
        terms = []
        nr = 0
        count = False

        if "-l" in self.cmdopts:
            return self.showFields()

        if "-c" in self.cmdopts:
            count = True

        # Parse query terms
        for clause in self.args[1:]:
            term = QueryTerm(clause)
            if term.operator:
                terms.append(term)

        # Execute query
        with open(self.outfile, "w") as out:
            if not count:
                out.write("SampleName\tRun")
                for term in terms:
                    out.write("\t" + term.varname)
                out.write("\n")
            self.opendb()
            try:
                for smp in self.getFilteredSamples():
                    for run in self.getSampleRuns(smp):
                        if not self.wantedruns or run in self.wantedruns:
                            data = self.getSampleData(smp, run)
                            if self.querySatisfied(terms, data):
                                nr += 1
                                if not count:
                                    out.write("{}\t{}".format(data["SampleName"], data["Run"]))
                                    for term in terms:
                                        out.write("\t{}".format(data[term.varname] if term.varname in data else "???"))
                                    out.write("\n")
            except BrokenPipeError:
                pass
            finally:
                self.closedb()
            if count:
                out.write("{}\n".format(nr))
        sys.stderr.write("*** {} samples returned.\n".format(nr))

    def showFields(self):
        typenames = {'T': "[Text] ",
                     'I': "[Int]  ",
                     'F': "[Float]"}
        sys.stdout.write("# Data fields:\n")
        self.opendb()
        try:
            for row in self.execute("SELECT DISTINCT Field FROM SampleData;"):
                if row['Field']:
                    field = row['Field']
                    tag = field[0]
                    name = field[1:]
                    sys.stdout.write("{} {}\n".format(typenames[tag], name))
        finally:
            self.closedb()
        sys.stdout.write("\n# Metadata fields:\n")
        for qf in QUERABLE_METADATA:
            sys.stdout.write("        {}\n".format(qf))
        sys.stdout.write("\n")

    def querySatisfied(self, terms, data):
        for term in terms:
            if term.varname in data:
                if not term.satisfied(data[term.varname]):
                    return False
            else:
                return not term.required
        return True

    # GISAID submission

    GISAID_PARAMS = {'submitter': 'ariva',
                     'orig_lab': 'unknown',
                     'orig_lab_address': 'unknown',
                     'laboratory': '"Salemi Lab, University of Florida"',
                     'lab_address': '"Emerging Pathogens Institute, 2055 Mowry Rd, Gainesville, FL 32610"',
                     'author': 'M. Salemi'}

    GISAID_HEADER = "submitter,fn,covv_virus_name,covv_type,covv_passage,covv_collection_date,covv_location,covv_add_location,covv_host,covv_add_host_info,covv_sampling_strategy,covv_gender,covv_patient_age,covv_patient_status,covv_specimen,covv_outbreak,covv_last_vaccinated,covv_treatment,covv_seq_technology,covv_assembly_method,covv_coverage,covv_orig_lab,covv_orig_lab_addr,covv_provider_sample_id,covv_subm_lab,covv_subm_lab_addr,covv_subm_sample_id,covv_authors,covv_comment,comment_type"

    def read_config(self, filename):
        with open(filename, "r") as f:
            for line in f:
                if not line:
                    continue
                if line[0] == "#":
                    continue
                p = line.find(":")
                if p > 0:
                    key = line[:p]
                    val = line[p+1:].strip()
                    self.GISAID_PARAMS[key] = val

    def getSampleInfo(self, fasta):
        result = []
        outfasta = os.path.splitext(fasta)[0] + ".gisaid.fa"
        sys.stderr.write("Fasta file for submission: {}\n".format(outfasta))
        n = 0
        with open(fasta, "r") as f, open(outfasta, "w") as out:
            for line in f:
                if not line:
                    continue
                if line[0] == '>':
                    hdr = line[1:].split("|")
                    result.append( (hdr[0], hdr[1]) )
                    out.write(">" + hdr[0] + "\n")
                    n += 1
                else:
                    out.write(line)
        sys.stderr.write("{} sequences prepared for upload to GISAID.\n".format(n))
        return outfasta, result

    def writeGISAID(self):
        if len(args) > 2:
            self.read_config(self.args[2])
        subfasta, sampleInfo = self.getSampleInfo(self.args[1])
        with open(self.outfile, "w") as out:
            sys.stdout.write(self.GISAID_HEADER + "\n")
            query = "SELECT SamplingDate, Location, Sex, Age, SampleType FROM Samples WHERE SampleName = ?"
            query2 = "SELECT Fval FROM SampleData WHERE SampleName = ? AND Field='FAvgDepth';"
            self.opendb()
            try:
                for (virus, sample) in sampleInfo:
                    data = self.execute(query, sample).fetchone()
                    cov = self.execute(query2, sample).fetchone()
                    if cov:
                        cov = str(int(cov[0]))
                    else:
                        cov = "unknown"
                    samplingdate = data[0][:7]
                    location = data[1]
                    if data[2] == "F":
                        gender = "Female"
                    elif data[2] == "M":
                        gender = "Male"
                    else:
                        gender = "unknown"
                    age = data[3] or "unknown"
                    sys.stdout.write(",".join([
                        self.GISAID_PARAMS['submitter'],
                        subfasta,
                        virus,
                        "betacoronavirus",
                        "Original",
                        samplingdate,
                        "North America / USA / Florida",
                        "",
                        "Human",
                        "",
                        "",
                        gender,
                        str(age),
                        "unknown",  # Patient status
                        data[4],    # Specimen source
                        "",
                        "",
                        "",
                        "Illumina NovaSeq",
                        "unknown",  # assembly method
                        cov,
                        self.GISAID_PARAMS['orig_lab'],
                        self.GISAID_PARAMS['orig_lab_address'],
                        "unknown",         # original sample ID
                        self.GISAID_PARAMS['laboratory'],
                        self.GISAID_PARAMS['lab_address'],
                        sample,
                        self.GISAID_PARAMS['author'],
                        "",
                        ""]) + "\n")

            finally:
                self.closedb()

    def gisaid_import(self):
        if len(self.args) < 3:
            return self.usage(self.args[0])
        idsfile = self.args[2]
        q1 = "SELECT GISAID FROM SampleGISAID WHERE SampleName=?"
        q2 = "UPDATE SampleGISAID SET GISAID=? WHERE SampleName=?"
        q3 = "INSERT INTO SampleGISAID(SampleName, GISAID) VALUES (?, ?)"
        n = 0
        self.opendb()
        try:
            with open(idsfile, "r") as f:
                for line in f:
                    parts = line.strip().split("; ")
                    if len(parts) == 2 and parts[1][:7] == "EPI_ISL":
                        n += 1
                        smpid = parts[0].split("/")[2]
                        if smpid.startswith("FL-"):
                            smpid = smpid[3:]
                        #sys.stderr.write(smpid + " " + parts[1] + "\n")
                        found = self.execute(q1, smpid).fetchone()
                        if found:
                            self.execute(q2, parts[1], smpid)
                        else:
                            self.execute(q3, smpid, parts[1])
        finally:
            self.closedb()
        sys.stderr.write("{} GISAID IDs imported.\n".format(n))

    def gisaid_ids(self, names):
        q = "SELECT SampleName, GISAID FROM SampleGISAID;"
        q2 = "SELECT GISAID FROM SampleGISAID WHERE SampleName=?;"
        self.opendb()
        try:
            with open(self.outfile, "w") as out:
                if names:
                    for name in names:
                        gisaid = self.execute(q2, name).fetchone() or "???"
                        out.write("{}\t{}\n".format(name, gisaid[0]))
                else:
                    FIELDS = ["SampleName", "GISAID"]
                    if self.fmt == "html":
                        W = fu.TSVWriterHTML(out, FIELDS)
                        W.title = "Sample Metadata"
                        W.formats = [""]*len(FIELDS)
                    else:
                        W = fu.TSVWriter(out, FIELDS)
                    W.header()

                    samplenames = self.getFilteredSamples()
                    for row in self.execute(q).fetchall():
                        if row[0] in samplenames:
                            W.writeRow(row)
        finally:
            self.closedb()

    def gisaidMain(self):
        if len(self.args) == 1:
            return self.usage(self.args[0])
        if args[1] == "ids":
            self.gisaid_ids(args[2:])
        elif args[1] == "import":
            self.gisaid_import()
        else:
            self.writeGISAID()

    def addControl(self):
        self.opendb()
        try:
            for smpname in self.args[1:]:
                self.execute("INSERT INTO Samples(SampleName, SamplingDate) VALUES (?, '2020-01-01');", smpname)
        finally:
            self.closedb()

    def run(self):
        cmd = self.args[0]
        if "-h" in args or "--help" in args:
            return self.usage(cmd)
        self.readWanted()
        try:
            if cmd == "init":
                self.initialize()
            elif cmd == "import":
                self.importSamples()
            elif cmd == "samples":
                self.getSamples()
            elif cmd[:3] == "run":
                self.getRunData()
            elif cmd[:4] == "meta":
                self.getMetadata()
            elif cmd[:5] == "stats":
                self.getSampleStats()
            elif cmd[:3] == "cov":
                self.getCoverageStats()
            elif cmd == "results":
                self.getSampleResults()
            elif cmd == "query":
                self.runQuery()
            elif cmd == "gisaid":
                self.gisaidMain()
            elif cmd == "lintable":
                self.lineagesTable()
            elif cmd == "addctl":
                self.addControl()
        except BrokenPipeError:
            pass

    def usage(self, what=None):
        sys.stdout.write("flacodb.py - the FLACO database.\n\n")
        if what is None:
            sys.stdout.write("""Usage: flacodb.py command [options]

where command is one of: run, import, metadata, stats, cov, results, query, gisaid. 
Use "-h command" to get help on a command. 

General options (that apply to all commands except 'run'):

  -db DB       | Use database file DB. Normally there is no need to specify this.
  -o O         | write output to file O (default: standard output).
  -f csv/html  | Choose output format (default: csv).
  -r R1,R2,... | Only output samples from the specified run(s).
  -g G1,G2,... | Only output samples from the specified group(s) (first part of sample name).
  -x G1,G2,... | Exclude samples from the specified group(s) (first part of sample name).
  -w S1,S2,... | Only output samples with the specified name(s).
  -W W         | Like -w, but reads sample names from file W (one per line, first column).
  -a A         | Only output samples with collection date equal to A or later.
  -b B         | Only output samples with collection date equal to B or earlier.
  -vax         | Only output vaccinated patients.
  -novax       | Only output unvaccinated patients.

Dates for the -a and -b options should be written in ISO format: YYYY-MM-DD. They
can be optionally followed by +N or _N to indicate N days before or after the 
specified date.

""")
        elif what == "run":
            sys.stdout.write("""Usage: flacodb.py run [runs...]

With no arguments, show all runs in the database in tab-delimited format, 
with the following columns:

  Run number
  Run name
  Run date
  Flowcell
  Number of samples

When run names are specified, writes a tab-delimited file with two columns:

    run    sample

for all samples in the specified runs.
""")
        elif what == "import":
            sys.stdout.write("""Usage: flacodb.py import metadata.txt

Read metadata for new samples from a tab-delimited file. The first row in the 
file should be a header containing column names.

""")

        elif what == "gisaid":
            sys.stdout.write("""Usage: flacodb.py gisaid sequences.fa [paramsfile]
       flacodb.py gisaid import ids-file
       flacodb.py gisaid ids [names...]

The first form of this command prepares a fasta file for submission to GISAID, writing 
the corresponding metadata. Sequences from the input fasta file are written to a new file 
with extension .gisaid.fa, leaving only the first component of the sequence names. Metadata 
are written to standard output, or to the file specified with the -o argument.

If present, paramsfile can be used to override some of the metadata fields on a 
per-project basis. It should be a text file with one "field: value" pair on each line.
The following fields can be specified in this file:

  submitter: [gisaid account name]
  orig_lab: [lab of origin for samples]
  orig_lab_address: [address of orig_lab]
  laboratory: [submitting laboratory]
  lab_address: [address of submitting laboratory]
  author: [people involved in submission]

Values containing commas should be enclosed in double quotes.


The second form reads a file created by submission to GISAID containing the EPI_ISL
identifier associated with each submitted sequence, storing them into the database.
""")


if __name__ == "__main__":
    args = sys.argv[1:]
    F = FLACOdb()
    if F.parseArgs(args):
        F.run()
    else:
        F.usage()

