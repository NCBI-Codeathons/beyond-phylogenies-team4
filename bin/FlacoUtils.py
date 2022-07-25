### Utils

from datetime import date, timedelta, datetime
from collections import defaultdict

LINEAGES = {'B.1.1.7': "B.1.1.7 (Alpha)",
            'B.1.617.2': "B.1.617.2 / AY (Delta)",
            'AY.': "B.1.617.2 / AY (Delta)",
            'BA.1': "BA.1 (Omicron)",
            'BA.2': "BA.2 (Omicron)",
            'BA.2.': "BA.2.* (Other BA.2)",
            'B.1': "B.1.* (other B.1)",
            'P.1': "P.1 (Gamma)"}

def dateDays(s):
    """Parse an ISO format date spec YYYY-MM-DD optionally followed by +N or _N,
where N is a number of days. Returns the date object after adding or subtracting
the specified days."""
    if "+" in s:
        parts = s.split("+")
        dt = date.fromisoformat(parts[0])
        delta = timedelta(days=int(parts[1]))
        return dt + delta
    elif "_" in s:
        parts = s.split("_")
        dt = date.fromisoformat(parts[0])
        delta = timedelta(days=int(parts[1]))
        return dt - delta
    else:
        return date.fromisoformat(s)

def formatDate(d):
    if d.count("/") == 2:
        parts = d.split("/")
        if len(parts[2]) == 2:
            parts[2] = "20" + parts[2]
        return "{}-{:02d}-{:02d}".format(parts[2], int(parts[0]), int(parts[1]))
    return d

def dateToDays(d):
    """Convert a date in ISO format to a number of days."""
    year = int(d[:4])
    month = int(d[5:7])
    day = int(d[8:])
    return date(year, month, day).toordinal()

def daysToDate(n):
    d = date.fromordinal(n)
    return d.isoformat()

def dget(dict, key, default=""):
    if key in dict:
        return dict[key]
    else:
        return default

def floatize(dict, key):
    """If key is in dict, replace its value with a floating point representation with two decimals."""
    if key in dict:
        dict[key] = "{:.2f}".format(dict[key])

### Function to build lineages table over sliding window

def lineagesInWindow(data, lineages):
    step = 7
    winsize = 28
    first = data[0][0]
    last  = data[-1][0]
    wstart = first
    wend   = wstart - winsize

    result = []

    while True:
        counts = {'Other': 0}
        for l in lineages:
            counts[l] = 0
        for d in data:
            if d[0] < wend:
                break
            if d[0] < wstart:
                addLineage(d[2], counts)
        nlin = sum(counts.values())
        if nlin:
            row = {'Date': daysToDate(wstart), "N": nlin}
            for l in lineages + ['Other']:
                row[l] = 100.0*counts[l] / nlin
            result.append(row)
        wstart = wstart - step
        wend = wstart - winsize
        if wstart < last:
            break
    return result

def addLineage(lin, counts):
    for x in LINEAGES:
        if lin.startswith(x):
            counts[LINEAGES[x]] += 1
            return
    counts['Other'] += 1

### Writers

class Writer(object):
    stream = None
    headers = None
    formats = None
    className = ""              # for HTML output

    def __init__(self, stream, headers, formats=None):
        self.stream = stream
        self.headers = headers
        self.formats = formats

    def header(self):
        self.stream.write("\t".join(self.headers) + "\n")

    def writeRow(self, row):
        self.writeRowTabDelimited(row)
        self.processRow(row)

    def writeRowTabDelimited(self, row):
        values = [ dget(row, f) for f in self.headers ]
        self.stream.write("\t".join([str(x) for x in values]) + "\n")

    def processRow(self, row):
        pass

    def close(self):
        pass

    def cell(self, value, fmt):
        if value == '':
            return "<TD class='{}'>&nbsp;</TD>".format(self.className)
        elif fmt == "D":
            return "<TD class='{}' align='right'>{:,}</TD>".format(self.className, int(value))
        elif fmt == "I":
            return "<TD class='{}'><I>{}</I></TD>".format(self.className, value)
        elif fmt == "B":
            return "<TD class='{}'><B>{}</B></TD>".format(self.className, value)
        elif fmt == "R":
            return "<TD class='{}' align='right'>{}</TD>".format(self.className, value)
        elif fmt == "P":
            return "<TD class='{}' align='right'>{}%</TD>".format(self.className, value)
        else:
            return "<TD class='{}'>{}</TD>".format(self.className, value)

class TSVWriter(Writer):

    def writeRowTabDelimited(self, row):
        self.stream.write("\t".join([str(x) for x in row]) + "\n")

class StatsWriter(Writer):
    nsamples = 0
    sumTotalReads = 0
    sumCovReads = 0

    def processRow(self, row):
        self.nsamples += 1
        self.sumTotalReads += int(dget(row, 'TotalReads', 0))
        self.sumCovReads += int(dget(row, 'CovReads', 0))

    def close(self):
        self.stream.write("\n")
        self.stream.write("Total reads:\t{}\n".format(self.sumTotalReads))
        self.stream.write("Total Cov reads:\t{}\n".format(self.sumCovReads))

class StatsWriterHTML(StatsWriter):
    title = ""
    className = "stats"

    css = """
BODY {
  margin: 0px;
  padding: 0px;
}
TABLE.stats {
  border: 2px solid black;
  border-collapse: collapse;
}
TH.stats {
  padding: 4px;
  background: #EEEEEE;
}
TD.stats {
  padding: 4px;
  border-top: 1px solid black;
}
.wide {
  width:90%;
}
.toggle {
  display: none;
}
TABLE.banner {
  border: 0px;
  border-spacing: 0px;
  border-collapse: collapse;
}
.banner {
  padding: 4px;
  background: #3333FF;
  color: white;
  font-weight: bold;
  font-family: sans-serif;
}
"""

    js = """
function openRow(id) {
  row = document.getElementById(id);
  if (row.style.display == "none") {
    row.style.display = "table-row";
  } else {
    row.style.display = "none";
  }
}
"""

    def header(self):
        self.stream.write("""<!DOCTYPE html>
<HTML>
  <HEAD>
    <TITLE>{}</TITLE>
    <STYLE>
{}
    </STYLE>
    <SCRIPT>
{}
    </SCRIPT>
  </HEAD>
  <BODY>
    <TABLE width='100%' class='banner'>
      <TR>
        <TD class='banner'>University of Florida</TD>
        <TD class='banner' align='center'>Salemi Laboratory</TD>
        <TD class='banner' align='center'>ICBR Bioinformatics Core</TD>
        <TD class='banner' align='right'>{}</TD>
      </TR>
    </TABLE>
    <CENTER>
      <H2>{}</H2>

      <TABLE class="{} wide">
        <TR><TH class="{}">&nbsp;</TH><TH class="{}">&nbsp;</TH>""".format(self.title, self.css, self.js, 
                                                                           datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                                                           self.title, self.className, self.className, self.className))
        
        for hdr in self.headers:
            self.stream.write("""<TH class="{}">{}</TH>""".format(self.className, hdr))
        self.stream.write("""</TR>\n""")

    def writeRow(self, row):
        self.stream.write("""<TR><TD class='{}' align='middle' onclick='openRow("row{}c");'><B>C</B></TD><TD class='{}' align='middle' onclick='openRow("row{}v");'><B>V</B></TD>""".format(self.className, self.nsamples, self.className, self.nsamples))
        for i in range(len(self.headers)):
            val = dget(row, self.headers[i])
            fmt = self.formats[i]
            self.stream.write(self.cell(val, fmt))
        self.stream.write("</TR>\n")
        if "PlotCov" in row:
            self.stream.write("""<TR id="row{}c" style="display: none" class="toggle"><TD colspan="14" align='center'><IMG src='{}' /></TD></TR>\n""".format(
                self.nsamples, row["PlotCov"]))
        if "OtherCov" in row:
            self.stream.write("""<TR id="row{}v" style="display: none" class="toggle"><TD colspan="14"><BLOCKQUOTE>{}</BLOCKQUOTE></TD></TR>\n""".format(self.nsamples, self.toTable(row["OtherCov"])))

        self.processRow(row)

    def toTable(self, data):
        tbl = """<TABLE class="{}">\n""".format(self.className)
        rows = data.split("|")
        row0 = rows[0].split(":")
        tbl += """<TR class="{}"><TH class="{}">{}</TH><TH class="{}">{}</TH><TH class="{}">{}</TH><TH class="{}">{}</TH></TR>\n""".format(
            self.className, self.className, row0[0], self.className, row0[1], self.className, row0[3], self.className, row0[4])
        for row in rows[1:]:
            row = row.split(":")
            if len(row) >= 5:
                tbl += """<TR class="{}"><TD class="{}"><I>{}</I></TD><TD class="{}" align="right">{}</TD><TD class="{}" align="right">{:.2f}</TD><TD class="{}" align="right">{:.2f}</TD></TR>\n""".format(self.className, self.className, row[0], self.className, row[1], self.className, float(row[3]), self.className, float(row[4]))
        tbl += "</TABLE>\n"
        return tbl

    def close(self):
        self.stream.write("""</TABLE>
</CENTER>
</BODY>
</HTML>
""")

class StatsWriterW3(StatsWriterHTML):
    template = ""
    toggles = False
    rows = []

    def header(self):
        self.rows = []

    def writeRow(self, data):
        self.rows.append(data)

    def close(self):
        with open(self.template, "r") as f:
            for line in f:
                if line.startswith("##Title"):
                    self.stream.write(self.title)
                elif line.startswith("##Header"):
                    if self.toggles:
                        self.stream.write("<th></th>")
                        col = 1
                    else:
                        col = 0
                    for i in range(len(self.headers)):
                        fld = self.headers[i]
                        fmt = self.formats[i]
                        kind = "kind='N'" if fmt in "FDRP" else ""
                        self.stream.write("""<th dir='A' {} onclick="sortTable('maintable', {});">{}</th>\n""".format(kind, col, fld))
                        col += 1
                elif line.startswith("##Rows"):
                    for r in self.rows:
                        self.outputRow(r)
                elif line.startswith("##Linplot"):
                    self.writeLinplot()
                elif line.startswith("##Pieplot"):
                    self.writePieplot()
                else:
                    self.stream.write(line)

    def cell(self, value, fmt):
        if value == '':
            return "<TD>&nbsp;</TD>"
        elif fmt == "D":
            return "<TD class='w3-right-align' val='{}'>{:,}</TD>".format(int(value), int(value))
        elif fmt == "I":
            return "<TD val='{}'><I>{}</I></TD>".format(value, value)
        elif fmt == "B":
            return "<TD val='{}'><B>{}</B></TD>".format(value, value)
        elif fmt == "R":
            return "<TD class='w3-right-align' val='{}'>{}</TD>".format(value, value)
        elif fmt == "P":
            return "<TD class='w3-right-align' val='{}'>{:.2f}%</TD>".format(value, float(value))
        else:
            return "<TD val='{}'>{}</TD>".format(value, value)

    def outputRow(self, row):
        ncols = len(self.headers) + 2
        self.stream.write("<TR>")
        if self.toggles:
            if "PlotCov" in row:
                self.stream.write("""<TD class='w3-center' onclick='openModal("{}");'><B>C</B></TD>""".format(row["PlotCov"]))
            else:
                self.stream.write("""<TD>&nbsp;</TD>""")
        for i in range(len(self.headers)):
            val = dget(row, self.headers[i])
            fmt = self.formats[i]
            self.stream.write(self.cell(val, fmt))
        self.stream.write("</TR>\n")
#        if self.toggles:
#            if "PlotCov" in row:
#                self.stream.write("""<TR id="row{}c" style="display:none"><TD colspan="{}" class='w3-center'><IMG src='{}' /></TD></TR>\n""".format(
#                    self.nsamples, ncols, row["PlotCov"]))
#            if "OtherCov" in row:
#                self.stream.write("""<TR id="row{}v" style="display:none"><TD colspan="{}"><BLOCKQUOTE>{}</BLOCKQUOTE></TD></TR>\n""".format(
#                    self.nsamples, ncols, self.toTable(row["OtherCov"])))
        self.processRow(row)

    def writeLinplot(self):
        ntr = 1
        traces = []
        dates = ", ".join(["'" + r["Date"] + "'" for r in self.rows])
        for L in self.headers[2:]:
            trace = "trace" + str(ntr)
            traces.append(trace)
            self.stream.write("""var {} = {{
  x: [{}],
  y: [{}],
  name: '{}',
  type: 'bar'
}};
""".format(trace, dates, ", ".join([str(r[L]) for r in self.rows]), LINEAGES[L] if L in LINEAGES else L))
            ntr += 1

        self.stream.write("""var data = [{}];\n""".format(", ".join(traces)))
        self.stream.write("""var layout = {barmode: 'stack'};
Plotly.newPlot('linplot', data, layout);
""")

    def lineagesTable(self):
        lincounts = defaultdict(int)
        for row in self.rows:
            if "PangoLineage" in row:
                lincounts[row["PangoLineage"]] += 1
        table = []
        nr = 0
        for k in lincounts:
            if k != "None":
                table.append([k, lincounts[k]])
                nr += lincounts[k]
        table.sort(key=lambda p: p[1], reverse=True)
        for row in table:
            pct = 100.0*row[1]/nr
            self.stream.write("""<TR>
  <TD>{}</TD>
  <TD class="w3-right-align">{}</TD>
  <TD class="w3-right-align">{:.2f}%</TD>
  <TD><PROGRESS value="{}" max="100">{}</PROGRESS></TD>
</TR>""".format(
                row[0], row[1], pct, int(pct), pct))

    def writePieplot(self):
        lineages = list(LINEAGES.values())
        counts = {'Other': 0}
        for l in lineages:
            counts[l] = 0
        total = 0
        for row in self.rows:
            if "PangoLineage" in row:
                addLineage(row["PangoLineage"], counts)
                total += 1
        if total:
            values = [ "{:.2f}".format(100.0 * counts[x] / total) for x in lineages + ["Other"] ]
            labels = [ "'" + x + "'" for x in lineages ] + ["'Other'"]

            self.stream.write("""var data = [{{
  values: [{}],
  labels: [{}],
  type: 'pie'
}}];
Plotly.newPlot('pieplot', data);
""".format(", ".join(values), ", ".join(labels)))

class TSVWriterHTML(StatsWriterHTML):

    def writeRow(self, row):
        self.stream.write("""<TR>""")
        for i in range(len(self.headers)):
            val = row[i]
            fmt = self.formats[i]
            self.stream.write(self.cell(val, fmt))
        self.stream.write("</TR>\n")

class ResultsWriterHTML(StatsWriterHTML):
    lincounts = defaultdict(int)

    def header(self):
        self.stream.write("""<!DOCTYPE html>
<HTML>
  <HEAD>
    <TITLE>{}</TITLE>
    <STYLE>
{}
    </STYLE>
    <SCRIPT>
{}
    </SCRIPT>
  </HEAD>
  <BODY>
    <TABLE width='100%' class='banner'>
      <TR>
        <TD class='banner'>University of Florida</TD>
        <TD class='banner' align='center'>Salemi Laboratory</TD>
        <TD class='banner' align='center'>ICBR Bioinformatics Core</TD>
        <TD class='banner' align='right'>{}</TD>
      </TR>
    </TABLE>
    <CENTER>
      <H2>{}</H2>

      <TABLE class="{} wide">
        <TR>""".format(self.title, self.css, self.js, 
                       datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                       self.title, self.className))
        
        for hdr in self.headers:
            self.stream.write("""<TH class="{}">{}</TH>""".format(self.className, hdr))
        self.stream.write("""</TR>\n""")

    def writeRow(self, row):
        self.stream.write("""<TR>""")
        for i in range(len(self.headers)):
            val = dget(row, self.headers[i])
            fmt = self.formats[i]
            self.stream.write(self.cell(val, fmt))
        self.stream.write("</TR>\n")
        self.processRow(row)

    def processRow(self, row):
        if 'PangoLineage' in row:
            l = row['PangoLineage']
            self.lincounts[l] += 1

    def close(self):
        table = []
        nr = 0
        for k in self.lincounts:
            table.append([k, self.lincounts[k]])
            nr += self.lincounts[k]
        table.sort(key=lambda p: p[1], reverse=True)
        self.stream.write("""</TABLE><BR><BR><TABLE class='{}'>
<TR class='{}'><TH class='{}'>Lineage</TH><TH class='{}'>Count</TH><TH class='{}'>Percent</TH></TR>
""".format(self.className, self.className, self.className, self.className, self.className))
        for row in table:
            self.stream.write("<TR class='{}'><TD class='{}'>{}</TD><TD class='{}' align='right'>{}</TD><TD class='{}' align='right'>{:.2f}%</TD></TR>\n".format(
                self.className, self.className, row[0], self.className, row[1], self.className, 100.0*row[1]/nr))
        self.stream.write("""</TABLE><BR>
</CENTER>
</BODY>
</HTML>
""")

