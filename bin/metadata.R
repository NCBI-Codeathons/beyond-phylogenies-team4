rm(list=ls())

require(ape)
require(optparse)
require(parallel)
require(dplyr)

`%notin%` <- Negate(`%in%`)
numCores=detectCores()

option_list = list(
make_option(c("-t", "--tree"), type="character", default=list.files(pattern=".nh"),
              help="new tree file [default= usher output]", metavar="character"),
make_option(c("-m", "--metadata"), type="character", default=list.files(pattern="metadata"),
              help="metadata file [default= contains 'metadata']", metavar="character"),
make_option(c("-n", "--columnName"), type="character", default="SampleName",
              help="metadata column containing sample/sequence ID [default='SampleName']", metavar="character") ,
make_option(c("-l", "--lineages"), type="character", default=list.files(pattern="lineages"),
              help="lineages file [default= pangolin output contains 'lineages']", metavar="character")
          
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


metadata = read.delim(opt$metadata, sep='\t', header=T, stringsAsFactors = F)
lineages = read.delim(opt$lineages, sep=',', header=T)
lineages = dplyr::rename(lineages, ID=taxon)


metadata = mutate(metadata, ID=metadata[,colnames(metadata) %in% opt$columnName]) %>%
  mutate(ID = paste(ID, gsub("(\\d{4}).+", "\\1", SamplingDate), sep="/"))
         
if(isTRUE(any(metadata$ID %notin% lineages$ID))) {
  for (i in 1:nrow(metadata)) {
    metadata$ID[i] = lineages$ID[grepl(metadata$ID[i], lineages$ID, fixed=T)]
  }
}
  
metadata = left_join(lineages, metadata, by="ID")
  
write.table(metadata, file=paste0("updated_metadata_", Sys.Date(), ".tab"), quote=F, row.names=F, sep='\t')










