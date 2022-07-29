#voc test script
#use if you want a random sample of lineages to test for filtering process

n = 20
voc_temp = data.frame(lineage = sample(metadata$lineage, n),
                      WHO = rep("unknown", n))

write_delim(voc_temp, "voc_test.tab", delim = "\t")

rm(n, voc_temp)