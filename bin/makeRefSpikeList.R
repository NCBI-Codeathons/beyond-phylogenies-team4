rm(list=ls())

require(ape)
require(dplyr)
require(optparse)
require(parallel)
require(foreach)
`%notin%` <- Negate(`%in%`)

today=gsub("-", "", Sys.Date())

option_list = list(
  make_option(c("-r", "--ref"), type="character",  
  default=paste0("/blue/salemi/brittany.rife/cov/phyloFLACO/ref_lineages/ref_lin_aln/ref_lineages_", today, ".fasta"),
  help="alignment of representative lineage sequences", metavar="character"),
  
  make_option(c("-l", "--reflin"), type="character", 
  default="/blue/salemi/brittany.rife/cov/phyloFLACO/ref_lineages/pango-designation/lineages.csv",
  help="list of lineages associated with representative lineage sequences", metavar="character"),
  
  make_option(c("-a", "--add"), type="character",
              default="Y",
              help="Specify (Y/N) whether to add to existing reference lineages")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


print("Reading in reference alignment and list.")

ref <- read.dna(opt$ref, format="fasta", as.matrix=T)
ref_list <- read.csv(opt$reflin, col.names=c("partialID", "Lineage"))

print("Converting partial names to full GISAID names.")
names_ref <- data.frame(ID=rownames(ref))
names_ref$partialID = gsub("hCoV-19/(.+/\\d{4})\\|.+", "\\1", names_ref$ID)

ref_list <- merge(names_ref, ref_list, by="partialID", all.x=TRUE, all.y=FALSE)

if (isTRUE(any(grepl("MN908947", names_ref$ID))==1) &
    any(grepl("MN908947", ref_list$ID))==0){
  ref_list <- rbind(ref_list, data.frame(
    partialID = names_ref$partialID[grepl("MN908947", names_ref$ID)],
    ID = names_ref$ID[grepl("MN908947", names_ref$ID)],
    Lineage = "Wuhan2"
  ))
}

print("Extracting spike protein sequences.")
## Look at just spike
spike_ref <- trans(ref[,21563:25384])

print("Grouping Spike protein reference sequences based on lineage and reformatting fasta files.")
r_reform <- list()
invisible(foreach (i=1:nrow(spike_ref)) %do% {
  r_reform[[i]] <- as.character(spike_ref[i,])[1,]
})
r_reform <- data.frame(do.call(rbind, r_reform))
rownames(r_reform) <- rownames(spike_ref)

r_reform <- r_reform[order(match(rownames(r_reform),ref_list$ID)),]


invisible(foreach (i=1:nrow(r_reform)) %do% {
  r_reform$lineage[i] <- ref_list$Lineage[ref_list$ID==rownames(r_reform)[i]]
})

r_reform2 <- group_split(r_reform, lineage, .keep=F) 
class(r_reform2) <- "list"
names(r_reform2) <- unique(r_reform$lineage)



if(opt$add=="N"){
print("Saving new reference list in RDS format")
saveRDS(r_reform2, "R_list_reference_spike.RDS")
} else {
  if(opt$add=="Y") {
    print("Retrieving original reference list for the addition of new lineages")
    r_reform_original <- readRDS("/blue/salemi/brittany.rife/cov/phyloFLACO/ref_lineages/R_list_reference_spike.RDS")
    
 # Merge matching names
    r_reform_original2 <- r_reform_original
       invisible(foreach(i=seq_along(r_reform_original)) %do% {
      invisible(foreach(j=seq_along(r_reform2)) %do% {
        if(names(r_reform2)[j] == names(r_reform_original)[i]) {
          r_reform_original2[[i]] <- merge(r_reform_original[[i]], r_reform2[[j]], all=T)
        } else {
          r_reform_original2[[i]] <- r_reform_original2[[i]]
        } # End if statement
      }) # End loop along r_reform2
       }) # End loop along r_reform_original
       
       new_lin <- names(r_reform2)[names(r_reform2) %notin% names(r_reform_original)]
       new_lin <- r_reform2[new_lin]
       
       r_reform_original2 <- append(r_reform_original2, new_lin)
       saveRDS(r_reform2_original2, "R_list_reference_spike.RDS")
        
    } # End if statement
    } # End if-else statement

