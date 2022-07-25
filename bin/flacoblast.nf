#!/usr/bin/env nextflow

params.cmd = "missing"
params.dbfile = "/dev/null"
params.blastcmd = "/dev/null"

/* makedb */

Channel
  .fromPath(params.dbfile)
  .splitText()
  .set { db }

process Makedb {
  module "ncbi_blast"
  memory "5G"
  time "2:00:00"

  when:
  params.cmd == "makedb"

  input:
  val onedbfile from db

  """
  #!/bin/bash
  makeblastdb -dbtype nucl -in ${launchDir}/${onedbfile}
  """
}

/* blast */

Channel
  .fromPath(params.blastcmd)
  .splitCsv(sep: '\t')
  .set { blastrows }

process Blast {
  module "ncbi_blast"
  memory "10G"
  time "24:00:00"

  when:
  params.cmd == "blast"

  input:
  val row from blastrows

  """
  #!/bin/bash
  DIR=${launchDir}/${row[0]}
  DB=${launchDir}/${row[1]}
  for fa in ${row[2]};
  do
    NAME=\${fa%.*}
    OUT=${launchDir}/\${NAME}.blast.csv
    blastn -db \$DB -query ${launchDir}/\$fa -outfmt 6 > \$OUT
  done
  """
}
