require(ape)
require(parallel)
require(foreach)
require(dplyr)

fa <- read.FASTA("updated_2021-08-16.aln")

aln_data <- data.frame(
  SampleName=gsub("hCoV-19/USA/FL-(.+)\\/.+", "\\1", names(fa)),
  SamplingDate=gsub("hCoV-19/USA/FL-.+/\\d{4}\\|.+\\|(\\d{4}-\\d{2}-\\d{2})\\|.+", "\\1", names(fa)),
  ID=names(fa)
)


dup_names <- select(read.delim("DUPLICATED-SAMPLES.txt", sep='\t', header=F, col.names=c("ROW", "SampleName", "SamplingDate")), -ROW) %>%
  filter(grepl("b", SampleName))

i=1
while(i<= nrow(aln_data)) {
  invisible(foreach(j = 1:nrow(dup_names)) %do% {
    if(isTRUE(grepl(aln_data$SampleName[i], dup_names$SampleName[j], ignore.case=TRUE) &
       aln_data$SamplingDate[i] == dup_names$SamplingDate[j])) {
      aln_data$ID_2[i] = gsub("(hCoV-19/USA/FL-)(.+)(\\/.+)", "\\1\\2b\\3", aln_data$ID[i])
    } else { aln_data$ID_2[i] <- NA}
  })
i=i+1
}

aln_data <- filter(aln_data, !is.na(ID_2))

i=1
while(i <= length(fa)) {
  invisible(foreach(j = 1:nrow(aln_data)) %do% {
    if(isTRUE(aln_data$ID[j] == names(fa)[i])) {
      names(fa)[i] = aln_data$ID_2[j] 
      } else {names(fa)[i] = names(fa)[i]}
  })
i=i+1
}

tree <- read.tree("uncondensed-final-tree.nh")

i=1
while(i <= length(tree$tip.label)) {
  invisible(foreach(j = 1:nrow(aln_data)) %do% {
    if(isTRUE(aln_data$ID[j] == tree$tip.label[i])) {
      tree$tip.label[i] = aln_data$ID_2[j]
      } else {tree$tip.label[i] = tree$tip.label[i]}
  })
i=i+1
}


write.FASTA(fa, "updated_2021-08-16_2.aln")
write.tree(tree, "uncondensed-final-tree2.nh")