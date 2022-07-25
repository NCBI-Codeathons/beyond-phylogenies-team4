rm(list=ls())

.packages = c("ape", "tidyr", "dplyr", "optparse", "parallel", "ggplot2", "grid", "lubridate")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# .inst <- .packages %in% installed.packages()
# if(length(.packages[!.inst]) > 0) BiocManager::install(.packages[!.inst], type = "source", checkBuilt = TRUE)

# Load packages into session 
invisible(lapply(.packages, library, warn.conflicts=FALSE, character.only=TRUE))

`%notin%` <- Negate(`%in%`)
numCores = try(Sys.getenv("SLURM_CPUS_ON_NODE"))
if (numCores == "") {
  numCores = detectCores()
}


option_list = list(
  make_option(c("-q", "--query"), type="character", 
              help="alignment to be queried against representative lineage sequences", 
              metavar="character"),
  make_option(c("-m", "--meta"), type="character", 
              #default=list.files(pattern="metadata.csv")[[1]],
              help="metadata file with from FLACODB or Panglin output", 
              metavar="character"),
  make_option(c("-s", "--subs"), type="character",
              default="spike_muts.txt",
              help="file with list of mutations",
              metavar="character"),
  make_option(c("-d", "--min_date"), type="character",
              default="2021-01-01",
              help="date filter",
              metavar="character")
  #  make_option(c("-f", "--ref"), type="character", 
  #              default="/blue/salemi/brittany.rife/cov/phyloFLACO/ref_lineages/R_list_reference_spike.RDS",
  #              help="RDS suptype list of Spike ref aln (default path is /blue/salemi/brittany.rife/cov/phyloFLACO/ref_lineages/R_list_reference_spike.RDS)", 
  #              metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("Reading in query alignment file. Please make sure your query sequences have been aligned to the MN9098947 reference sequence!")

fa <- read.dna(opt$query, format="fasta", as.matrix=T)

#Pull out reference
ref <- fa[which(grepl("MN908947", rownames(fa))),]
# If reference strain found in query, remove


print("Reading in query metadata")

message("Searching for metadata file...")
if(endsWith(opt$meta, "csv")) {
  metadata <- read.csv(opt$meta, header=T, stringsAsFactors = F)
} else {
  if(endsWith(opt$meta, "tab") | endsWith(opt$meta, "txt")) {
    metadata <- read.table(opt$meta, sep='\t', header=T, stringsAsFactors = F)
  }
}



## Extract fasta headers
names_q <- data.frame(ID=rownames(fa))

## Allowing for use of DYNAMITE-style metadata file or simple output from Pangolin (needs revising)
print("Determining metadata file format.")


##The following is needed because flacodb only provides a partial name for the sequences, but the sequence headers are in GISAID format!
if("SampleName" %in% colnames(metadata)) {
  names_q$partialID = gsub("hCoV-19/USA/FL-(.+)/\\d{4}\\|.+", "\\1", names_q$ID)
  print("Checking that all sequences are accounted for in metadata...")
  if (isTRUE(all(metadata$SampleName %in% names_q$partialID))) {
    print("All sequences are accounted for.")
  } else {
    print("Not all sequences are accounted for. Still proceeding, but please check that metadata is correct.")
    metadata <-  metadata[which(metadata$SampleName %in% names_q$partialID),] # Need to add this because not all consensus sequences are being retrieved for some reason.
      }
  metadata <- merge(metadata, names_q, by.x="SampleName", by.y="partialID") %>%
    rename(ID=SampleName)
  
  print("FLACODB metadata file detected.")
} else {
  metadata <- merge(names_q, metadata, by.x="ID", by.y="taxon") 
  print("Pangolin output file detected.")
}
#Get rid of duplicate sequences
# If reference strain found in query, remove
fa <- fa[which(!grepl("MN908947", rownames(fa))),]
metadata <- metadata[which(!grepl("MN908947", metadata$ID)),]

#Get rid of duplicate sequences
fa <- fa[which(!duplicated(rownames(fa))),]
metadata <- metadata[which(!duplicated(metadata$ID)),]




# Filter metadata on specific date
metadata <-   mutate(metadata, Date = as.Date(gsub(".+\\|(\\d{4}-\\d{2}-\\d{2})\\|.+", "\\1", ID))) %>%
  filter(Date >= opt$min_date)

metadata <- mutate(metadata,
                   PangoLineage = if_else(grepl("B.1.1.7|Q", lineage), "Alpha",
                                 if_else(grepl("B.1.351", lineage), "Beta",
                                         if_else(grepl("P.1", lineage), "Gamma",
                                                 if_else(grepl("B.1.429|B.1.427", lineage), "Epsilon",
                                                         if_else(grepl("B.1.525", lineage), "Eta",
                                                                 if_else(grepl("B.1.526", lineage), "Iota",
                                                                         if_else(grepl("B.1.617.1", lineage), "Kappa",
                                                                                 if_else(grepl("B.1.617.3", lineage), "B.1.617.3",
                                                                                         if_else(grepl("B.1.621|B.1.621.1", lineage), "Mu",
                                                                                                 if_else(grepl("P.2", lineage), "Zeta",
                                                                                                         if_else(grepl("B.1.617.2|AY", lineage), "Delta",
                                                                                                                 if_else(grepl("B.1.1.529|BA", lineage), "Omicron",
                                                                                                                         if_else(grepl("None", lineage), "Unknown",
                                                                                                                         "Other"))))))))))))))
#  mutate(Variant = if_else(WHO %in% VBM | WHO %in% VOC, WHO, "Other")) 

## Filter out non-Alachua sequences
metadata <- metadata %>% 
  mutate(location=gsub("hCoV-19/USA/FL-([A-Za-z]+)-.+", "\\1", ID)) %>%
  mutate(location=if_else(grepl("shands|path|STP|UF|SED", location, ignore.case=T), "Alachua", "FL")) %>%
  filter(location=="Alachua")


# May need this if good dates not obtainable:
fa <- fa[rownames(fa) %in% metadata$ID,]




print("Creating query Spike alignment")

spike_fa <- trans(fa[,21563:25384])


print("Reformatting query Spike alignment for cross-referencing.")
# Now query fasta
q_reform <- mclapply (1:nrow(spike_fa), function(i) {
  as.character(spike_fa[i,])[1,]
}, mc.cores=numCores)
q_reform <- data.frame(do.call(rbind, q_reform))
rownames(q_reform) <- rownames(spike_fa)

q_reform <- q_reform[order(match(rownames(q_reform),metadata$ID)),]

muts <- read.delim(opt$subs, col.names = "name") %>%
  mutate(wt = gsub("([A-Z]+)\\d+.+", "\\1", name),
         site = as.numeric(gsub("[A-Z]+(\\d+).+", "\\1", name)),
         mut = gsub("[A-Z]+\\d+(.+)", "\\1", name))

MOI <- mclapply (1:nrow(metadata), function(qseq) {
#   
#   if (isFALSE(any(q_reform[,452] == "L" | q_reform[,452] == "R"))) {
#     print("Neither L/R found at position 452. Please check alignment.")
#     stop() }
#   if (isFALSE(any(q_reform[,452] == "L" | q_reform[,452] == "Q"))) {
#     print("Neither L/Q found at position 452. Please check alignment.")
#     stop() }
#   if (isFALSE(any(q_reform[,484] == "E" | q_reform[,484] == "K"))) {
#     print("Neither E/K found at position 484. Please check alignment.")
#     stop() }
#   if (isFALSE(any(q_reform[,490] == "F" | q_reform[,490] == "S"))) {
#     print("Neither F/S found at position 490. Please check alignment.")
#     stop() }
#   if (isFALSE(any(q_reform[,501] == "N" | q_reform[,501] == "Y"))) {
#     print("Neither N/Y found at position 501. Please check alignment.")
#     stop() }
#   
  result <- lapply(1:nrow(muts), function(i) {
     if(isTRUE(q_reform[qseq,muts$site[i]] %in% c(muts$wt[i], "X"))) {
       return("")
       } else {
         return(q_reform[qseq,muts$site[i]])
       }
    })
  result <- as.data.frame(do.call(cbind, result))
#   
#   #For binary
#   result <- lapply(1:nrow(muts), function(i) {
#     if(isTRUE(q_reform[qseq,muts$site[i]] %in% c(muts$wt[i], "X"))) {
#       return(0)
#     } else {
#       return(1)
#     }
#   })
#   result <- as.data.frame(do.call(cbind, result))
#   
  colnames(result) = muts$name
  result$ID=metadata$ID[qseq]
  return(result)
}, mc.cores=numCores)

MOI <- do.call(rbind, MOI)
MOI <- left_join(MOI, select(metadata, ID, PangoLineage), by="ID") %>%
  mutate(Date=gsub(".+(\\d{4}-\\d{2}-\\d{2}).+", "\\1", ID)) %>%
  select(ID, Date, everything())
write.csv(MOI, file=paste0("Mutations_of_Interest_", Sys.Date(), ".csv"), quote=F, row.names=F)
# 
# binned <- gather(omicron_muts, mutation, present, muts$name[1]:last(muts$name)) %>%
#   mutate(PangoLineage = if_else(grepl("AY.*", PangoLineage), "Delta", PangoLineage),
#          PangoLineage = if_else(PangoLineage=="B.1.617.2", "Delta", PangoLineage)) %>%
#   filter(PangoLineage!="None") %>%
#   mutate(Date=as.Date(Date)) 
# 
# 
# df <- mutate(binned, day = as.numeric(1),
#              week = isoweek(Date),
#              month = lubridate::month(Date, label = TRUE),
#              year = lubridate::year(Date),
#              weekdate = as.Date(paste(year,week,day), format="%Y %U %u")) %>%
#   mutate(weekdate = if_else(week==53 & year==2020, as.Date(paste(year, 52, day), format="%Y %U %u"), 
#                             if_else(week==53 & year==2021, as.Date(paste(2021, 1, day), format="%Y %U %u"), weekdate))) %>%
#   group_by(PangoLineage, mutation, weekdate) %>%
#   arrange(weekdate) %>%
#   summarize(n=sum(present)) %>%
#   group_by(PangoLineage, mutation) %>%
#   mutate(cum_n=cumsum(n)) %>%
#   distinct()
# 
# 
# major <- group_by(df, PangoLineage, mutation) %>%
#   summarize(max_cum=max(cum_n)) %>%
#   filter(max_cum >= 100)
# 
# major <- unique(major$PangoLineage)
# 
# df2 <- filter(df, PangoLineage %in% major)
# 
# 
# ggplot(df2) +
#   geom_bar(aes(x=weekdate, y=cum_n, fill=PangoLineage), position="stack", stat="identity", width=7) +
#   theme_minimal() +
#   scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
#   scale_x_date(date_breaks = "1 month", date_labels = "%b")+
#   facet_wrap(~mutation, ncol=2, scales = "free_y") +
#   ylab("Number of genomes") +
#   xlab("Time") +
#   #  scale_fill_viridis_d() +
#   scale_fill_brewer(palette="Set3")
# 
# ggsave(plot=last_plot(), "Omicron_mutations_time.pdf", width=11, height=8.5, units="in")
# 




#print("Reading in list of reference Spike alingments (by lineage).")
#r_reform <- readRDS(opt$ref)
r_reform=list(MN908947=trans(ref[,21563:25384]))
r_reform[[1]] <- as.character(r_reform[[1]])[1,]

print("Finding mutational differences.")
#not_found <- data.frame(ID=NA, Lineage=NA)
mut_diffs <- mclapply(1:nrow(metadata), function(qseq) {

  #lin = which(names(r_reform) == metadata$PangoLineage[qseq])
  lin=1
  vec <- vector()
  vec[1] <- metadata$ID[qseq]
  vec[2] <- metadata$PangoLineage[qseq]
  
  #if (isTRUE(length(lin) != 0)) {

    muts <- lapply(seq_along(q_reform[qseq,]), function(aa) {
       

      if (isTRUE(any(q_reform[qseq,aa] == r_reform[[lin]][aa]) | # Needs to be changed back to [,aa] if pangolin lineages used again
                 q_reform[qseq,aa] == "X")){
        muts[[aa]] <- ""
      } else {
        #Don't need to report all reference lineage amino acids if all the same!
#        if (isTRUE(length(unique(unlist(r_reform[[lin]][,aa]))) == 1)) {
#          mut_diffs[[qseq]][aa+1] <- paste0(metadata$PangoLineage[qseq],
#                                          "-",
#                                          unique(unlist(r_reform[[lin]][,aa])),
#                                          aa,
#                                          q_reform[qseq,aa])
#        } else {
        muts[[aa]] <- paste0(paste(unlist(r_reform[[lin]][aa]), collapse="/"), # Needs to be changed back to [,aa] if pangolin lineages used again
                                          aa,
                                          q_reform[qseq,aa])
 #       }# End if-else statement
      } # End if-else statement
    }) # End for loop
    muts <- as.vector(t(do.call(rbind, muts)))
    
#  } else {
#    mut_diffs[[qseq]] <- c(metadata$ID[qseq], rep(NA, ncol(q_reform[qseq,])))
#    not_found <- rbind(not_found, data.frame(ID=metadata$ID[qseq], Lineage=metadata$PangoLineage[qseq]))
#  } # End if-else statement
    return(c(vec, muts))
}, mc.cores=numCores) # End for loop
mut_diffs_df <- as.data.frame(do.call(rbind, mut_diffs)) 
#mut_diffs_df <- Filter(function(x) length(unique(x[x!=""]))>1, mut_diffs_df)


#If lineages not found in reference, print message and save file
#if(isTRUE(nrow(not_found)>1)) {
#  print(paste("The following lineages were present in the query but not in the reference list - ",
#              unique(not_found$Lineage), ".Need to update from pangolin github."))
#  write.csv(not_found[-1,], paste0("lineages_not_represented_", Sys.Date(), ".csv"), quote=F, row.names=F)
#}



#Create frequencies table for mutation differences
message("Calculating frequencies for mutations.")


#Get rid of all sequences for which lineages were not represented (in a different list)
mut_diffs_df <- mut_diffs_df[rowSums(is.na(mut_diffs_df)) != ncol(mut_diffs_df), ]
colnames(mut_diffs_df) <- c("ID", "Lineage", rep(1:(ncol(mut_diffs_df)-2)))

#print("Saving report of mutational changes (does not include sequences with non-represented lineages)")

#mut_diffs_df <- lapply(1:ncol(mut_diffs_df), function(i) {
#  lapply(1:nrow(mut_diffs_df), function(j) {
#    if(grepl("X", mut_diffs_df[j,i])) {
#    mut_diffs_df[j,i] <- ""}
#  })
#})

write.csv(mut_diffs_df, paste0("Mutation_changes_sites_", Sys.Date(), ".csv"), quote=F, row.names=F)

mut_diffs_dt <- gather(mut_diffs_df, Site, Mutation, "1":"1274") %>%
  filter(Mutation!= "") %>%
  group_by(Lineage, Site) %>%
  mutate(n=length(ID)) %>%
  group_by(Lineage, Site, Mutation) %>%
  summarize(Frequency=length(ID)/n) %>%
  distinct() %>%
#  tidyr::separate(Mutation, c("Lineage", "Mutation"), sep="-") %>%
  filter(!grepl("X", Mutation)) %>%
  mutate(Site=as.numeric(Site)) %>%
  arrange(Site) 

write.csv(mut_diffs_dt, file=paste0("Mutation_changes_frequency_", Sys.Date(), ".csv"), quote=F, row.names=F)

## If you would like to focus on specific lineages, just modify the code below:

# mut_diffs_dt <- filter(mut_diffs_dt, 
#                        Lineage=="P.1" | Lineage=="B.1.1.7")


## Plotting ############################################################################################################
print("Plotting frequencies along Spike gene")

require(ggplot2)

# Domain information obtained from http//www.ncbi.nlm.nih.gov/protein/YP_009724390.1 on the Wuhan-1 strain
print("Plotting mutations.")
domains=data.frame(x1=c(0, 13, 319, 543, 918, 1162, 1233),
                   x2=c(680, 304, 541, 1208, 983, 1203, 1273),
                   y1=c(0, 0, 0, 0, 0, 0, 0),
                   y2=c(Inf, 1, 1, Inf, 1, 1, 1),
                   Domain=c('S1 subunit', 'NTD', 'RBD',
                            'S1/S2 cleavage site + \nS2 subunit',
                            'HR1', 'HR2', 'CoV_S2_C'))
domains$Domain <- factor(domains$Domain, levels=c('S1 subunit', 'NTD', 'RBD',
                                                  'S1/S2 cleavage site + \nS2 subunit',
                                                     'HR1', 'HR2', 'CoV_S2_C'))



#note <- textGrob("***For mutations at Freq<1%,\n see Frequency table", gp=gpar(fontsize=13, fontface="bold", color="red"))

mut_diffs_dt %>%
  filter(., Frequency<1) %>%
  mutate(Mutation = if_else(Frequency<=0.01, "***", Mutation)) %>%
  drop_na() %>%
  ggplot() +
  geom_rect(data=domains, aes(xmin=x1-0.5, xmax=x2+0.5, ymin=y1, ymax=y2, fill=Domain), color="white", alpha=0.8) +
  geom_bar(aes(x=Site, y=Frequency, group=Mutation), stat="identity", position="dodge", color="black", width=1) +
  geom_text(aes(x=Site, y=Frequency, group=Mutation, label=Mutation),  vjust=0.01,  hjust=0.01, angle=45) +
  facet_wrap(~Lineage, drop=T, ncol=1) +
  theme_minimal() +
  scale_y_continuous(limits=c(0,1)) +
  theme(panel.grid = element_blank(), text = element_text(size=20), axis.text.x = element_text(angle=45)) +
  xlab("Amino acid position") +
  coord_cartesian(ylim=c(0,1), clip="off") +
#  scale_fill_brewer(palette = "Set1")
  scale_fill_brewer(type = "qual")
#  annotation_custom(note, xmin=1200,xmax=2000,ymin=0.10,ymax=0.10)

n <- length(unique(mut_diffs_dt$Lineage))
ggsave(plot=last_plot(), file=paste0("Lineage_mutations_", Sys.Date(), ".pdf"), width=11, height=3*n, units="in")






# Test ##############################################################################################################

#r_list <- data.frame(ID=paste0("Ref", seq(1:4)), lineage=c("B", "B", "B.1", "B.1"))
#q_list <- data.frame(ID=paste0("Q", seq(1:4)), lineage=c("B", "B", "B.1", "B.1"))

#q <- trans(rDNAbin(nrow=4, ncol=30000, base.freq = rep(0.25, 4), prefix = "Q"))
#r <- trans(rDNAbin(nrow=4, ncol=30000, base.freq = rep(0.25, 4), prefix = "Ref"))


# Reformat reference fasta in to list for each lineage

# r_reform <- list()
# for (i in 1:nrow(r)) {
#   r_reform[[i]] <- as.character(r[i,])[1,]
# }
# r_reform <- data.frame(do.call(rbind, r_reform))
# rownames(r_reform) <- rownames(r)
# 
# r_reform <- r_reform[order(match(rownames(r_reform),r_list$ID)),]
# 
# for (i in 1:nrow(r_reform)) {
#   r_reform$lineage[i] <- r_list$lineage[r_list$ID==rownames(r_reform)[i]]
# }
# 
# r_reform <- group_split(r_reform, lineage, keep=F) %>%
#   setNames(unique(r_reform$lineage))
# 
# # Now query fasta
# q_reform <- list()
# for (i in 1:nrow(q)) {
#   q_reform[[i]] <- as.character(q[i,])[1,]
# }
# q_reform <- data.frame(do.call(rbind, q_reform))
# rownames(q_reform) <- rownames(q)
# 
# q_reform <- q_reform[order(match(rownames(q_reform),q_list$ID)),]
# 
# 
# 
# 
# mut_diffs <- list()
# for (qseq in 1:nrow(q_list)) {
#   lin = which(names(r_reform) == q_list$lineage[qseq])
#   
#   mut_diffs[[qseq]] <- vector(length=ncol(q_reform[qseq,]))
#   
#   for (aa in seq_along(q_reform[qseq,])) {
#     
#     if (isTRUE(any(q_reform[qseq,aa] == r_reform[[lin]][,aa]) & q_reform[qseq,aa] != "X")){
#       mut_diffs[[qseq]][aa] <- ""
#     } else {
#       mut_diffs[[qseq]][aa] <- paste0(r_list$lineage[qseq], 
#                                   "-", 
#                                   paste(unlist(r_reform[[lin]][,aa]), collapse="/"), 
#                                   aa, 
#                                   q_reform[qseq,aa])
#     }
#   }
# }
# 
# names(mut_diffs) = q_list$ID
# 
# 
# mut_diffs2 <- do.call(rbind, mut_diffs)
# mut_diffs_dt <- table(mut_diffs2)



