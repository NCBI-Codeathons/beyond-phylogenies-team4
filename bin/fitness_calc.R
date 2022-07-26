rm(list=ls())

require(ape)
require(parallel)
require(lubridate)
require(dplyr)
require(optparse)
require(rlist)
require(tidyr)


numCores=detectCores()
`%notin%` <- Negate(`%in%`)


option_list = list(
  make_option(c("-f", "--fasta"), type="character", default=list.files(pattern="*.fasta")[[1]], 
              help="new fasta file [default= .fasta extension]", metavar="character"),
  make_option(c("-p", "--path"), type="character", default="/blue/salemi/brittany.rife/cov/test/",
              help="path to metadata file [default= /blue/salemi/brittany.rife/cov/test/]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Read in all trees in folder and assign names
message("Reading in trees...")
filenames <- list.files(pattern=".tree", full.names=F)
trees <- lapply(filenames, function(t) {
  t2 <- read.tree(t)
  t2$edge.length[is.nan(t2$edge.length)] <- 0
  return(t2)
}) 
names(trees) <- gsub("(.)\\.tree", "\\1", filenames)

# Read in the fasta with all recently sampled sequences and extract headers (ignoring reference sequence)
message("Reading in most recent fasta...")
fa <- read.FASTA(opt$fasta) 
taxids <- data.frame(ID=names(fa))
taxids <- filter(taxids, ID!="NC_045512.2")

# Read in metadata file from specified path
message("Reading in most recent metadata file")
#metadata <- read.csv(file=paste0(opt$path,list.files(pattern="metadata")[[1]]), header=T, stringsAsFactors = F)
metadata <- read.csv(file=paste0("../", list.files(path="../", pattern="metadata")[[1]]), header=T, stringsAsFactors = F)
#metadata <- read.csv("/Users/macbook/Dropbox (UFL)/COVID/in_house/Starting_tree/test_using_partial/metadata.csv")

# Subdivide metadata information according to trees
message("Subdividing metadata according to trees...")
seq_data <- lapply(trees, function(t) {
  df <- data.frame(ID=taxids[taxids$ID %in% t$tip.label,]) %>%
    left_join(., select(metadata, ID, country, DATE))
  return(df)
}) 

## Determine fitness for background data only ##################################################################
get_bg_data <- function() {
  bg_new_seq_data <- seq_data[grep("background", names(seq_data))][[1]]
  if (nrow(bg_new_seq_data)>0) {
    bg_new_seq_data$origin <- "added"
    countries <- bg_new_seq_data$country
    bt <- trees[grep("background", names(trees))][[1]]
    assign("bt", bt, envir = globalenv())
    bg_old_seq_data <- data.frame(ID=bt$tip.label[bt$tip.label %notin% taxids$ID])
    bg_old_seq_data <- left_join(bg_old_seq_data, select(metadata, ID, country, DATE)) %>%
      filter(country %in% unique(countries))
    bg_old_seq_data$origin <- "initial"
    
    bg_combined <- list()
    for (i in 1:nrow(bg_new_seq_data)) {
      bg_combined[[i]] <- rbind(bg_old_seq_data[bg_old_seq_data$country==bg_new_seq_data$country[i],], 
                                bg_new_seq_data[i,])
    }
    
    bt_pruned <- vector("list", length = length(bg_combined))
    bt_added <- vector("list", length = length(bg_combined))
    
    for (i in seq_along(countries)) {
      # Remove all non-related countries from original data and all new sequence data
      bt_pruned[[i]] <- drop.tip(bt, c(metadata$ID[metadata$country!=countries[i]],
                                       taxids$ID)) 
      # Remove all non-related countries from original data and the added sequence for each list
      bt_added[[i]] <- drop.tip(bt, c(metadata$ID[metadata$country!=countries[i]],
                                      taxids$ID[taxids$ID != last(bg_combined[[i]]$ID)]))
    }
    assign("bt_pruned", bt_pruned, envir = globalenv())
    assign("bt_added", bt_added, envir = globalenv())
  } else {
    bg_combined <- NULL
  }
  
  return(bg_combined)
}
bg_data <- get_bg_data()



if(!is.null(bg_data)) {
  # Calculate time fitness score
  message("Calculating the temporal fitness...")
  time_fitness <- function(group) {
    n.samples = length(group$ID)
    if(n.samples>1){
      first.date = min(decimal_date(as.Date(group$DATE)))
      last.date = max(decimal_date(as.Date(group$DATE)))
      ideal.time = seq(from = first.date,
                       to = last.date,
                       by =  (last.date-first.date)/(n.samples-1)) #n.samples is the number of genomes
      
      ### worst possible time spread, all sequences in one time point
      max.time.spread = as.numeric(sum(abs(ideal.time-rep(first.date, n.samples))))
      
      timeframe = decimal_date(as.Date(group$DATE))
      
      # how good this sample is in terms of time spread?
      time.spread = as.numeric(sum(abs(sort(timeframe)-ideal.time)))
      
      # fitness score with rescaling
      tem.diversity = 1 - (time.spread/max.time.spread) +10^(-64)
    } else {
      tem.diversity=0
    }
    
    #There you go.  tem.diversity is a [0, 1] bound value, the higher the better.
    return(tem.diversity)
  }
  
  F_score <- list(pre=NA, post=NA)
  F_score$pre <- lapply(bg_data, function(x) {
    x = x[x$origin=="initial",]
    result <- time_fitness(x)
    return(result)
  })
  F_score$post <- mclapply(bg_data, function(x) {
    result <- time_fitness(x)
    return(result)
  }, mc.cores=numCores)
  F_score <- data.frame(pre=unlist(F_score$pre), post=unlist(F_score$post))
  
  for (i in 1:nrow(F_score)) {
    if(F_score$pre[i]==0 & F_score$post[i]==0) {
      F_score$post[i] = 1
    }
  }
  
  
  # Calculate diversity fitness score
  message("Calculating the genetic fitness...")
  calculate_p_dist <- function(pruned_tree, added_tree) {
    d_score_pre <- cophenetic(pruned_tree) 
    d_score_pre <- melt(d_score_pre)[melt(upper.tri(d_score_pre))$value,] %>%
      mutate(Var="pre")
    d_score_post <- cophenetic(added_tree)
    d_score_post <- melt(d_score_post)[melt(upper.tri(d_score_post))$value,] %>%
      mutate(Var="post")
    
    d_diff <- rbind(select(d_score_pre, Var, value), select(d_score_post, Var, value))
    p <- wilcox.test(value~Var, data=d_diff)
    return(p$p.value)
  }
  p_dist <- mcmapply(function(x,y) {
    tryCatch(calculate_p_dist(x, y), 
             error=function(e) e=0)
  }, bt_pruned, bt_added, mc.cores=numCores)
  

## Create discard pile for sequences that do not add temporal signal and genetic diversity
  message("Creating discard pile...")
  discard <- NA
  for (i in 1:nrow(F_score)) {
    if(isTRUE(F_score$post[i] <= F_score$pre[i] | 
              p_dist[i] > 0.05)) {
      discard[i] <- last(bg_data[[i]]$ID)
    }
  }
  discard <- discard[!is.na(discard)]
  
  ## Output fasta for discard pile and downstream tree reconstruction (needs to have .aln extension)
  message("Saving discard pile sequences to fasta for tree reconstruction (.aln extension)")
  discard_fa <- fa[names(fa) %in% discard]
 
    
  write.FASTA(discard_fa, file=paste0("../discard_", Sys.Date(), "/discard_", Sys.Date(), ".aln"))
  
  
  # Update existing background tree
  message("Updating existing background tree...")
  bt <- drop.tip(bt, discard)
  path_to_bg_tree <- paste0("background_", Sys.Date(), "_updated")
  write.tree(bt, file=paste0("../", path_to_bg_tree, "/", path_to_bg_tree, ".tree"))
  
  ## Update existing background fasta
  message("Updating existing background fasta file...")
#  old_bg_fa <- read.FASTA(paste0("../", path_to_bg_tree, "/", path_to_bg_tree, ".fasta"))
  new_bg_fa <- fa[names(fa) %in% bt$tip.label]
  write.FASTA(new_bg_fa, paste0("../", path_to_bg_tree, "/", path_to_bg_tree, ".fasta"))
  
} else{
  message("All samples were assigned to clusters.")
}



### Determine if added sequences part of cluster or clade ###############################################
message("Parsing cluster data...")
get_cluster_data <- function() {
  
  cluster_new_seq_data <- seq_data[grep("^c", names(seq_data))]
  cluster_new_seq_data <- list.filter(cluster_new_seq_data, nrow(.) > 0)
  
  if (length(cluster_new_seq_data)>0) {
    cluster_new_seq_data <- mclapply(cluster_new_seq_data, function(c) {
      c$origin <- "added"
      return(c)
    }, mc.cores=numCores)
    
    
    cluster_trees <- trees[grep("^c", names(trees))]
    
    cluster_old_seq_data <- list()
    for (i in seq_along(cluster_new_seq_data)) {
      cluster_old_seq_data[[i]] <- data.frame(ID=cluster_trees[[i]]$tip.label[cluster_trees[[i]]$tip.label %notin% taxids$ID])
      cluster_old_seq_data[[i]] <- left_join(cluster_old_seq_data[[i]],
                                             select(metadata, ID, country, DATE))
    }
    cluster_old_seq_data <- mclapply(cluster_old_seq_data, function(c) {
      c$origin <- "initial"
      return(c)
    }, mc.cores=numCores)
    
    clust_combined <- list()
    cluster_trees2 <- vector("list", length=length(cluster_new_seq_data))
    
    for (i in seq_along(cluster_new_seq_data)) {
      clust_combined[[i]] <- vector("list", length=nrow(cluster_new_seq_data[[i]]))
      for (j in 1:nrow(cluster_new_seq_data[[i]])) {
        clust_combined[[i]][[j]] <- rbind(cluster_old_seq_data[[i]], cluster_new_seq_data[[i]][j,])
        cluster_trees2[[i]][[j]] <- drop.tip(cluster_trees[[i]], taxids$ID[taxids$ID != cluster_new_seq_data[[i]]$ID[j]])
      }
      names(clust_combined)[i] <- gsub("(c\\d+).+", "\\1", names(cluster_trees)[[i]])
      names(cluster_trees2)[i] <- gsub("(c\\d+).+", "\\1", names(cluster_trees)[[i]])
    }
    assign("cluster_trees2", cluster_trees2, envir=globalenv())
  } else {
    clust_combined <- NULL
  }
  return(clust_combined)
}
clust_data <- get_cluster_data()

if(!is.null(clust_data)) {
  # Read in branch length threshold
  message("Using branch length cutoff from initial cluster identification...")
  branch_length_limit <- read.table("../branch_length_limit.txt", header=F)
  
  clades <- list()
  for (i in seq_along(cluster_trees2)) {
    clades[[i]] <- mclapply(cluster_trees2[[i]], function(x) {
      x_tbl <- tidytree::as_tibble(x) 
      return(x_tbl)
    }, mc.cores=numCores)
    names(clades)[i] <- names(cluster_trees2)[i]
  }  
  
  # Use branchwise algorithm to determine if new sequence should be considered part of cluster or not.
  message("Calculating mean branch lengths...")
  bl_score <- list()
  for (i in seq_along(clades)) {
    bl_score[[i]] <- mclapply(clades[[i]], function(x) {
      score <- NA
      mean_bl <- mean(x$branch.length)
      if (mean_bl <= branch_length_limit ) {
        score = 1
      } else { 
        score=0 }
      return(data.frame(score=score, ID=x$label[x$label %in% taxids$ID]))
    }, mc.cores=numCores)
    names(bl_score)[i] <- names(clades)[i]
  }
  
  ## Subset bl_score according to positive cluster association
  message("Determining if added cluster sequences actually part of clusters...")
  bl_score <- mclapply(bl_score, function(x) {
    x <- list.filter(x, score==1)
    return(do.call(rbind, x))
  }, mc.cores=numCores)
  bl_score <- plyr::compact(bl_score)
  
  if(isTRUE(length(bl_score) >0)) {
    # If part of cluster, update traits file
    message("Updating cluster traits file...")
    cluster_traits <- read.csv(paste0("../", list.files(path="../", pattern="trait*")[[1]]), stringsAsFactors = F) %>%
      group_by(cluster_id) %>%
      group_split()
    
    cluster_traits <- lapply(cluster_traits, spread, field, trait) 
    
    for (i in seq_along(bl_score)) {
      for (j in seq_along(cluster_traits)) {
        if (cluster_traits[[j]]$cluster_id[1] == names(bl_score)[i]) {
          cluster_traits[[j]] <- left_join(cluster_traits[[j]], metadata[metadata$ID %in% bl_score[[i]]$ID,])
        }
      }
    }
    
    # Convert cluster_traits back to long format
    cluster_traits <- lapply(cluster_traits,  function(x) {
      gather(x, field, trait, -cluster_id, -ID, -DATE, -birth_origin, factor_key=TRUE)
    })
    cluster_traits <- do.call(rbind, cluster_traits) 
    
    # Just replace old file for now.
    write.csv(cluster_traits, file=paste0("../", list.files(path="../", pattern="trait*")[[1]]))
  } else{
    message("Added samples do not contribute to cluster activity...")
  }
} else {
  message("No samples were added to cluster subtrees.")
}

