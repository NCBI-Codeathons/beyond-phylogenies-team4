rm(list=ls())

require(ape)
require(optparse)
require(parallel)
require(dplyr)

`%notin%` <- Negate(`%in%`)
numCores=detectCores()

option_list = list(
make_option(c("-o", "--origin"), type="character", default="USA",
              help="country of origin of metadata [default=USA]", metavar="character"),
make_option(c("-m", "--metadata"), type="character", default=list.files(pattern="metadata"),
              help="metadata file [default= contains 'metadata']", metavar="character"),
make_option(c("-n", "--columnName"), type="character", default="SampleName",
              help="metadata column containing sample/sequence ID [default='SampleName']", metavar="character") ,
make_option(c("-l", "--lineages"), type="character", default=list.files(pattern="lineages"),
              help="lineages file [default= pangolin output contains 'lineages']", metavar="character")
          
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


metadata = read.delim(opt$metadata, sep='\t', header=T, stringsAsFactors = F) %>%
  mutate(Country = opt$origin,
         Run = gsub(".+\\(\\d{4}-\\d{2}-\\d{2}).+", "\\1", opt$metadata))

lineages = read.delim(opt$lineages, sep=',', header=T)
lineages = dplyr::rename(lineages, ID=taxon) %>%
  mutate(DATE = gsub(".+\\|(\\d{4}-\\d{2}-\\d{2})\\|.+", "\\1", ID)) %>%
  mutate(DATE=ifelse(grepl("00", DATE), gsub("00", "01", DATE), DATE)) %>%
  mutate(DATE=as.Date(DATE)) %>%
  mutate(Country = gsub("hCoV-19/([A-Za-z]+)/.+", "\\1", ID)) %>%
  dplyr::select(ID, DATE, lineage)

metadata = mutate(metadata, ID=metadata[,colnames(metadata) %in% opt$columnName]) %>%
  mutate(ID = paste(ID, gsub("(\\d{4}).+", "\\1", SamplingDate), sep="/"))
         
if(isTRUE(any(metadata$ID %notin% lineages$ID))) {
  for (i in 1:nrow(metadata)) {
    metadata$ID[i] = ifelse(length(lineages$ID[grepl(metadata$ID[i], lineages$ID, fixed=T)])>0,
                            lineages$ID[grepl(metadata$ID[i], lineages$ID, fixed=T)], NA)
  }
}
  
for (i in seq_along(colnames(metadata))) {
  if (tryCatch({
    isTRUE(any(grepl("date$", colnames(metadata)[i], ignore.case = T)))}, error = function(e) stderr() )) {
    colnames(metadata)[i] <- "DATE"
    metadata$DATE=as.Date(metadata$DATE)
  } else {
    mutate(metadata,
           DATE = as.Date(gsub(".+\\|(\\d{4}-\\d{2}-\\d{2})\\|.+", "\\1", ID)))
  }
}
  
  metadata = left_join(lineages, metadata, by=c("ID", "DATE")) %>%
    dplyr::select(ID, DATE, everything())



updated_meta = list.files(path="../", pattern="updated_metadata")

if(isTRUE(length(updated_meta)>0)) {
	dates = as.Date(gsub(".+(\\d{4}-\\d{2}-\\d{2}).+", "\\1", updated_meta))
	latest_date = max(dates)
	latest_meta = paste0("../", updated_meta[grepl(latest_date, updated_meta)])
	latest_meta = read.delim(latest_meta, sep='\t', header=T) %>%
	  mutate(DATE=as.Date(DATE))
	
	for(i in 1:ncol(metadata)) {
	  for (j in 1:ncol(latest_meta))
	  if(colnames(metadata)[i] == colnames(latest_meta)[j] &
	     class(metadata[,i]) != class(latest_meta[,j])) {
	    class(metadata[,i]) = class(latest_meta[,j])
	  }
	}
	metadata = full_join(latest_meta, metadata)
	}
  
write.table(metadata, file=paste0("../updated_metadata_", Sys.Date(), ".tab"), quote=F, row.names=F, sep='\t')










