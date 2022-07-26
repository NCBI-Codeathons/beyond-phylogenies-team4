rm(list=ls())

require(parallel)
require(purrr)
require(dplyr)
require(tidyverse)
require(ape)


BPS <- purrr::map_dfr(list.files(pattern="*.tsv", full.names = TRUE),
                      ~read_tsv(.x) %>% mutate(file = sub(".tsv$", "", basename(.x))))

#BPS <- filter(BPS, `#Sample` != "MN908947.3")

# Find min(BPS) and filter samples, excluding files not corresponding to min(BPS)
#min_BPS <- BPS %>%
#  group_by(`#Sample`) %>%
#  mutate(min_BPS = min(`Parsimony score`)) %>%
#  filter(`Parsimony score` == min_BPS) 


#If duplicated across clusters, we're in trouble...
#if (any(duplicated(min_BPS$`#Sample`))) {
#  message("Samples with equivalent placement across multiple clusters detected. Need to give a look...")
#} else {
  BPS <- group_by(BPS, file) %>%
    group_split()
  
 # new_path <- paste0("samples_", Sys.Date())
 # fa <- read.FASTA(file=paste0("../", new_path, "/",  new_path, ".fasta") )
  yyyymm <- gsub(".+(\\d{4}-\\d{2}).+", "\\1", BPS[[1]]$file[1])
  file=list.files(path=paste0("../", yyyymm), pattern="_masked.aln$", full.names=T)
  fa <- read.FASTA(file)

  sub_fa <- list()
  for (i in seq_along(BPS)) {
    sub_fa[[i]] <- fa[names(fa) %in% BPS[[i]]$`#Sample`]
    write.FASTA(sub_fa[[i]], 
                file=paste0(gsub("(.+)_parsimony-scores", "\\1", BPS[[i]]$file[1]), "_updated.fasta"))}
  
