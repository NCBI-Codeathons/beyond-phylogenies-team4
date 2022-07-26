rm(list=ls())
#setwd(getSrcDirectory()[1]) #When working on the cluster
#dirname(rstudioapi::getActiveDocumentContext()$path) # If working in Rstudio
rt0 <- Sys.time()

# List of packages for session
.packages <-  c("optparse", "remotes", "phytools", "data.tree", 
                "tidytree", "lubridate", "rlist", "familyR", "tidyverse", 
                "ggtree", "treeio", "parallel", "geiger", "tibble", "treedater")  # May need to incorporate code for familyR (https://rdrr.io/github/emillykkejensen/familyR/src/R/get_children.R) i fno longer supported.

## Will need to remove install section if using on cluster ###################################
# Install CRAN packages (if not already installed)
# Load packages into session 
lapply(.packages, require, character.only=TRUE)

option_list = list(
   make_option(c("-t", "--tree"), type="character", default="", 
              help="tree file name [default= .nwk extension]", metavar="character"),
  make_option(c("-s", "--seqLen"), type="numeric", default=30000, 
              help="sequence length [default= 10000]", metavar="numeric"),
  make_option(c("-c", "--cluster"), type="character", default="b", 
              help="choice of cluster algorithm from c (Phylopart's cladewise) or b (DYNAMITE's branchwise) [default= dynamite]", metavar="character"),
  make_option(c("-l", "--leaves"), type="character", default="addLeaves", 
              help="choice of transformation to tree from bifurcating or addLeaves [default=addLeaves]", metavar="character"),
  make_option(c("-p", "--path"), type="character", default="../",
              help="new file directory [default= <main directory>]", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$tree)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Additional options
`%notin%` <- Negate(`%in%`) # Just plain useful
# The option below is useful when dealing with dates of internal nodes downstream
options(digits=15)

numCores = try(Sys.getenv("SLURM_CPUS_ON_NODE"))
if (numCores == "") {
  numCores = detectCores()
}

write("Checking for tree file...")
checkFortree <- function(tree_file) {
  print("Make sure tree_file is scaled in substitutions/site and not in units of time.")
  # Check tree_file format -- either Newick or Nexus -- using suffices
  if (endsWith(tree_file, "nwk") || 
      endsWith(tree_file,"newick") || 
      endsWith(tree_file,"tree") ||
      endsWith(tree_file,"treefile" )){
    print(paste("Newick file detected:", tree_file, sep=" "), stderr())
    sub_tree <- read.tree(tree_file)
  } else if (endsWith(tree_file, "nex") || endsWith(tree_file,"nexus") || endsWith(tree_file, "nxs")){
    print("Nexus file detected")
    sub_tree <- read.tree(tree_file) ##Note: will not read in bracketed annotations!
    return(tree_file)
  }  else {
    # Neither Newick nor Nexus identified -- stop
    print("Cannot identify tree_file format. Use nexus (nex,nxs,nexus), newick (nwk,newick), our BEAST output (tre, tree_files) formats.")
  }# end conditional statement
  #assign("sub_tree_file", sub_tree_file, envir=globalenv())
  return(sub_tree)
}# end checkfortree_file function
sub_tree <- checkFortree(opt$tree)

metadata <- read.csv(file=paste0("../",list.files(path="../", pattern="metadata")[[1]]), 
                              header=T, stringsAsFactors = F) %>%
  dplyr::filter(., ID %in% sub_tree$tip.label)

checkForDates <- function(metadata) {
  for (i in seq_along(colnames(metadata))) {
    if (tryCatch({
      isTRUE(any(grepl("^id", colnames(metadata)[i], ignore.case = T)))}, error = function(e) stderr() )) {
      colnames(metadata)[i] <- "ID"
    } #end first ef statement
    if (tryCatch({
      isTRUE(any(grepl("date$", colnames(metadata)[i], ignore.case = T)))}, error = function(e) stderr() )) {
      colnames(metadata)[i] <- "DATE"
    } # End second if statement
  } # End for loop
  
  assign("metadata", metadata, envir = globalenv())
  createSTS <- function(metadata) {
    date_range <- seq(as.Date("01/01/1900", "%d/%m/%Y"), as.Date(Sys.Date(), "%d/%m/%Y"), by="day")
    sts <- metadata$DATE
    if (isTRUE(class(sts) == "integer")) {
      sts <- as.Date(ISOdate(sts))
    } else if (isTRUE(class(sts) == "numeric")) {
      sts <- as.Date(date_decimal(sts))
    } else if (isTRUE(class(sts) == "character")) {
      sts <- as.Date(sts, tryFormats = c("%Y-%m-%d", "%Y/%m/%d", "%m/%d/%Y", "%m-%d-%Y"))
      if (isTRUE(min(sts, na.rm=T) < min(date_range))) { # Attempt to correct if in wrong format, which we will know if the first date is before 1900s
        sts <- as.Date(metadata$DATE, tryFormats = c("%d/%m/%Y", "%d-%m-%Y"))
      }
    } # if-else statement
    metadata$DATE <- sts
    assign("metadata", metadata, envir = globalenv()) # Need dates for metadata table as well for downstream analysis
    sts <- setNames(sts, metadata$ID)
    sts <- sts[names(sts) %in% sub_tree$tip.label]
    return(sts)
  }
  sts <- createSTS(metadata)
  
  ## Checkpoint
  #  if (isTRUE(max(as.Date(ISOdate(sts, 1, 1)))> Sys.Date())) {
  if (isTRUE(max(sts) > Sys.Date())) {
    
    write(paste("The following sequences likely have incorret dates:",
                return(sts[sts==max(sts)]),
                "Please start DYNAMITE again, placing a .tab metadata file with consistent date information in current folder.",
                sep=" "))
    stop()
  } # End checkpoint if statement
  ## Will need sts in downstream analyses
  return(sts)
} # End checkForDates function
sts <- checkForDates(metadata)

DateTree <- function(sub_tree, seqLen) {
  # result <- lsd2(inputTree=sub_tree,
  #                inputDate=decimal_date(sts),
  #                estimateRoot="as",
  #                constraint=TRUE,
  #                variance=1,
  #                ZscoreOutlier = 3,
  #                outFile = "lsd2_results",
  #                seqLen = seqLen,
  #                nullblen=-1
  #                )
  # assign("sub_tree", result$newickTree, envir=globalenv())
  # assign("time_tree", result$dateNexusTreeFile@phylo, envir=globalenv())
  # # assign("time_tree_data", result$dateNexusTreeFile, envir=globalenv()) ## Something wrong with node numbers in this tree, so grab from time_tree
  # assign("tmrca", result$tMRCA, envir=globalenv())
  sub_tree <- sub_tree
  result <- dater(sub_tree, decimal_date(sts), s=seqLen, ncpu=numCores, omega0=8E-04)
  assign("sub_tree", result$intree, envir=globalenv())
  assign("time_tree", as.phylo(as_tibble(result)), envir=globalenv())
  # assign("time_tree_data", result$dateNexusTreeFile, envir=globalenv()) ## Something wrong with node numbers in this tree, so grab from time_tree
  assign("tmrca", result$timeOfMRCA, envir=globalenv())
  
}# End DateTree function
DateTree(sub_tree, opt$seqLen)

# Function to specify most recent sampling date (mrsd)
findMRSD <- function(time_tree) {
  date.mrsd <- max(sts[names(sts) %in% time_tree$tip.label])
  mrsd <- decimal_date(date.mrsd)
  
  gg_tree <- ggtree(time_tree, mrsd=date.mrsd) + theme_tree2()
  assign("gg_tree", gg_tree, envir = globalenv())
  node_dates <- gg_tree$data ## Could also do nodeHeights here....
  assign("node_dates", node_dates, envir=globalenv())
  assign("num.mrsd", mrsd, envir=globalenv())
  assign("date.mrsd", date.mrsd, envir=globalenv())
  print(paste0("The updated most recent sampling date is ", date.mrsd))
} # End findMRSD function
findMRSD(time_tree)

write("Defining well-supported clades within the tree using node labels...")
define.clades <- function(sub_tree) {
  ## Grab only subtrees that are supported by bootstrap values >90
  ## We may need to change this in case people have other support values
  ## Note that subtrees() function in ape renumbers nodes, so useless here, since at the end we wish to recombine the information
  write("Creating list of all subtrees...This may take a while.")
  family_tree <- tidytree::as_tibble(sub_tree)
  ## Need to relabel columns so that "parent" and "node" are "from" and "to" for familyR::get_children function
  colnames(family_tree)[1:2] <- c("from", "to")
  ## The dataframe needs to be transformed back into a data.table for future analyses
  family_tree <- data.table::as.data.table(family_tree)
  assign("family_tree", family_tree, envir = globalenv())
  
  write("Extract well-supported subtrees.")
  supported_nodes <- list()
  if(max(as.numeric(sub_tree$node.label), na.rm=T) >1) {
    supported_nodes <- suppressWarnings(dplyr::filter(family_tree, as.numeric(label) > 90))
  } else {
    supported_nodes <- suppressWarnings(dplyr::filter(family_tree, as.numeric(label) > 0.90))
  }
  assign("supported_nodes", supported_nodes, envir = globalenv())  
  
  ## Creating a list for all sub_trees using the familyR package (get_children function)
  clades <- list()
  for (node in unique(supported_nodes$to)) {
    clades[[node]] <- familyR::get_children(family_tree, node)
  } # End loop along family_tree_supported
  ## This function introduces several null items in the list, which can be removed by the following:
  clades<-plyr::compact(clades)
  ## Merge information across 'nodes' and 'edges' dataframes within each nested list
  for (i in seq_along(clades)) {
    clades[[i]] <- merge(clades[[i]]$edges,clades[[i]]$nodes, by.x = "to", by.y = "id")
  } # End loop along clades
  ## Restructure list of clades for easy visualization and manipulation and remove
  ## Zeroth level (contains root branch length) from first clade (full tree) if exists
  clades <- mclapply(clades, function(x) {
    dplyr::select(x, from, everything()) %>%
      arrange(level) %>%
      filter(!level==0)
  }, mc.cores=numCores) %>%
    list.filter(from[1] != rootnode(sub_tree)) # Whole tree is included because support is alwasy 100% for root, so discard# End loop along clades
  
  # Save to global environment or merge next function.
  assign("clades", clades, envir=globalenv())
  
} # End defineClades function
define.clades(sub_tree)

merge.nested.clust <- function(clusters) {
  copy <- clusters
  result <- list()
  unwanted <- list()
  for (ct in seq_along(clusters)) {
    for (cc in seq_along(copy)) {
      if (isTRUE(all(clusters[[ct]]$to %in% copy[[cc]]$to) &
                 length(clusters[[ct]]$to) != length(copy[[cc]]$to))) {
        unwanted[[ct]] <- clusters[[ct]]
      } else{NULL}
    } # End loop along copy
  } # End loop along true
  result <- setdiff(clusters, unwanted)
  for (j in seq_along(result)) {
    names(result)[[j]] <- paste0("c", j)
  }
  return(result)
}  
merge.overlap.clust <- function(clusters) {
  copy <- clusters
  unwanted <- list()
  result <- list()
  for (ct in seq_along(clusters)) {
    for (cc in seq_along(copy)) {
      if (isTRUE(sum(copy[[cc]]$to %in% clusters[[ct]]$to) > 0.1*length(copy[[cc]]$to)) &
          isTRUE(names(copy)[[cc]] != names(clusters)[[ct]])) {
        unwanted[[cc]] <- copy[[cc]]
        clusters[[ct]] <- full_join(copy[[cc]], clusters[[ct]], by=c("from", "to", "branch.length", "label"))
      } else{clusters[[ct]] <- clusters[[ct]]}
    } # End loop along copy
  } # End loop along true
  result <- setdiff(clusters, Filter(Negate(function(x) is.null(unlist(x))), unwanted)) %>%
    lapply(., function(x){
      dplyr::select(x, from, to, branch.length, label) %>%
        dplyr::arrange(from,to) 
    }) %>%
    unique()
  
  return(result)
}  # In case you want to remove this and consider only fully nested clusters
branchLengthLimit <- read.table("../branch_length_limit.txt", header=F)[1,1]
branchWise <- function(tree, branch_length_limit, make_tree) {
  pickClust <- function(clade){
    
    ## Need to add mean_bl column  to original clade list
    clade$mean_bl <- rep(Inf, nrow(clade))
    
    
    ## Set initial values
    current_level <- as.numeric(2)
    current_mean_bl <- as.numeric(-Inf)
    
    
    ## Initiate subclade using first two edges connected to root of clade (level=1)
    sub_clade <- filter(clade, level==1,
                        branch.length <= branch_length_limit)
    
    
    
    while(isTRUE(current_mean_bl <= branch_length_limit)){
      
      # Create a vector of nodes sampled from the subsequent level
      # Each iteration chooses amongst nodes that are connected to the current sub-clade:
      
      next_level_nodes <- filter(clade, level == current_level,
                                 from %in% sub_clade$to)
      
      # A shortlist of possible enlargements of the sub-clade is kept to be able
      # to compare each potential enlargement of the sub-clade and always keep the enlargement
      # if the mean branch length is under the limit
      #
      # The shortlist is enlarged by vertices that are:
      #  1) adjacent to the most recent added node(s)
      #  2) not already IN the sub_clade
      new_node <- dplyr::setdiff(next_level_nodes, sub_clade)
      sub_clade <- rbind(sub_clade,new_node, fill=T)
      
      
      # The branch length is NOT calculated by the branch length of an individual
      # edges leading to nodes in the shortlist BUT on the mean of the nodes in the previous level
      # and added node.
      for (x in 1:nrow(sub_clade)) {
        if (isTRUE(sub_clade$level[x] < current_level &
                   sub_clade$mean_bl[x] != Inf)) {
          sub_clade$mean_bl[x] <- sub_clade$mean_bl[x]
        } else{
          sub_clade$mean_bl[x] <- sum(c(sub_clade$branch.length[x],
                                        subset(sub_clade$branch.length, sub_clade$level < sub_clade$level[x])))/
            length(c(sub_clade$branch.length[x],
                     subset(sub_clade$branch.length, sub_clade$level < sub_clade$level[x])))
        }
      }
      
      # Identify nodes with mean branch length > current mean branch length limit
      # and remove from shortlist
      unwanted_edges_at_current_level <- filter(sub_clade, level==current_level, mean_bl >= branch_length_limit)
      sub_clade <- setdiff(sub_clade, unwanted_edges_at_current_level)
      
      ## Redefine current level and mean branch length, taking into account whether
      ## or not all nodes belonging to the current level have been filtered out
      if (max(sub_clade$level)==current_level) {
        #current_mean_bl <- min(sub_clade$mean_bl[level=current_level])
        current_mean_bl <- mean(sub_clade$branch.length)
        current_level <- current_level+1
      } else {
        current_mean_bl = Inf
        current_level <- current_level+1
      } # End ifelse statment
    } # End tree traversal (while loop)
    if (nrow(sub_clade) == 0) {
      return(NULL)
    }else {
      return(sub_clade)
    }
  } # End pickClust function
  addLeaves <- function(tree, clusters) {
    sub_tree <- as_tibble(tree)
    x <- clusters  # Copy clusters list for ease
    for (c in seq_along(x)) {
      x[[c]] <- x[[c]] %>% dplyr::rename(parent=from, node=to) # comparison of clusters and subtree is based on same column names
      for (node in 2:nrow(x[[c]])) {
        if (length(x[[c]]$parent[x[[c]]$parent==x[[c]]$parent[node]]) == 1) { # For each row in a cluster, find parent nodes
          ## that only have one edge (need two for downstream analysis)
          p <- x[[c]]$parent[node]
          n <- x[[c]]$node[node]
          l <- as.data.frame(filter(sub_tree, parent == p & node != n)) ## Find the other edge in sub_tree
          l$branch.length =
            x[[c]]$branch.length[x[[c]]$parent == l$parent &
                                   x[[c]]$node != l$node] ## Replace the current branch.length with the branch.length of the
          ## other edge in the current list
          x[[c]] <- rbind(x[[c]], l) ## Add the new edge to the current cluster
        } # End if statement ## If no lonely edges, let the cluster list be itself
      } # End loop along row
    } # End loop along clusters
    return(x)
  } # End function; this added articial leaves to create bifurcating tree, but changed to the bifurcate f(x) above, which just involes dropping tips.
  
  #  branch_length_limit <- findThreshold(sub_tree)
  
  clusters <- mclapply(clades, pickClust, mc.cores=numCores) %>%
    compact() %>%
    merge.nested.clust() %>%
    merge.overlap.clust() %>%
    merge.nested.clust()
  ## Remove singleton nodes (non-bifurcating branches) by using bifurcate() function or extracting entire clade.
  
  clusters <- addLeaves(tree, clusters)
  clusters <- mclapply(clusters, as_tibble, mc.cores=numCores)
  clusters <-  list.filter(clusters, sum(label %in% tree$tip.label) >= 5)
  
  return(clusters)
}

message("Picking clusters according to initially defined branch length limit")
clusters <- branchWise(sub_tree, branchLengthLimit, make_tree=opt$leaves)
new_name <- gsub("(.+)_updated.+", "\\1", opt$path)

if (isTRUE(length(clusters) > 0)) {
  names(clusters) <- paste0(new_name, "_", names(clusters))
  
  
  message("Writing cluster (and background) trees to separate files...")
  #Have no idea what cluster numbers already in folder, so need origin and date information in file name
  
  clusterPhylo <- function(clusters, sub_tree, time_tree) {
    
    clusters_time_tree <- mclapply(clusters, function(x) {
      tips <- x$label[x$label %in% time_tree$tip.label]
      full_clade <- extract.clade(time_tree, findMRCA(time_tree, tips))
      cluster <- keep.tip(full_clade, tips) # Need to make sure not full clade, but only tips present in cluster
      return(cluster)
    }, mc.cores=numCores)
    assign("clusters_time_tree", clusters_time_tree, envir = globalenv() )
    for (i in seq_along(clusters_time_tree)) {
      # write.tree(clusters_time_tree[[i]], file=paste0("../", new_name, "_", names(clusters_sub_tree)[i], ".tree"))
    }
    
    clusters_sub_tree <- mclapply(clusters, function(x) {
      tips <- x$label[x$label %in% sub_tree$tip.label]
      full_clade <- extract.clade(sub_tree, findMRCA(sub_tree, tips))
      cluster <- keep.tip(full_clade, tips) # Need to make sure not full clade, but only tips present in cluster
      return(cluster)
    }, mc.cores=numCores)
    assign("clusters_sub_tree", clusters_sub_tree, envir = globalenv() )
    for (i in seq_along(clusters_sub_tree)) {
      write.tree(clusters_sub_tree[[i]], file=paste0("../", new_name, "_", names(clusters_sub_tree)[i], ".tree"))
    }
  }
  clusterPhylo(clusters, sub_tree, time_tree)
  
  unclusteredPhylo <- function(clusters_sub_tree, sub_tree, time_tree) {
    cluster_tips <- unlist(mclapply(clusters_sub_tree, function(x) {
      x$tip.label
    }, mc.cores=numCores))
    background_subtree <- drop.tip(sub_tree, cluster_tips)
    assign("background_subtree", background_subtree, envir = globalenv() )
    write.tree(background_subtree, file=opt$tree)
    background_timetree <- drop.tip(time_tree, cluster_tips)
    #   write.tree(background_timetree, file=paste0("background_timetree_", Sys.Date(), ".tree"))
  }
  unclusteredPhylo(clusters_sub_tree, sub_tree, time_tree)
  
  write("Performing quick data manipulation... Hold on to your butts!")
  dataManip <- function(clusters) {
    
    ## Need to convert clusters to list if only a single cluster found (will be in dataframe format, rather than list)
    
    if (class(clusters) == "data.frame") {
      clusters <- list(clusters) 
    } else if (class(clusters) == "list") {
      clusters <- clusters
    } else {warning("class of clusters data unknown, error in dataManip() function or above")}
    
    metadata2 <- gather(metadata, field, trait, -ID, -DATE, -GISAID, factor_key=TRUE)
    cluster_data <- list()
    
    if (isTRUE(exists("asr", envir = globalenv()))) { #Only use if ASR incorporated above
      asr2 <- bind_rows(asr, .id = "field")
      asr2$to <- as.integer(asr2$to)
      
      anc.state.data <- mclapply(clusters, function(x) {
        anc.field <- dplyr::select(asr2, field, to) %>%
          dplyr::filter(to == x$parent[1])
        anc.trait <- dplyr::select(asr2,trait, to) %>%
          dplyr::filter(to == x$parent[1])
        x <- merge(anc.field, anc.trait, by="to")
        return(x)
      }, mc.cores = numCores)
      
      
      for (i in seq_along(clusters)) {
        cluster_data[[i]] <- merge(clusters[[i]], metadata2, by.x = "label", by.y="ID", all=T) %>%
          dplyr::select(parent, everything())
        cluster_data[[i]] <- merge(cluster_data[[i]], anc.state.data[[i]], by.x=c("parent", "field", "trait"), by.y=c("to", "field", "trait"), all=T)
      } 
    } else {
      
      for (i in seq_along(clusters)) {
        clusters[[i]]$ID <- clusters[[i]]$label
        cluster_data[[i]] <- left_join(clusters[[i]], metadata2, by="ID") %>%
          dplyr::select(parent, node, ID, DATE, GISAID, field, trait)
      }
      
      
    }
    names(cluster_data) <- names(clusters)
    return(cluster_data)
    
  } # End dataManip() function
  cluster_data <- dataManip(clusters)
  
  ## Now DYNAMITE determines if clusters are related by connecting the children of each cluster to the root of remaining clusters
  connectClust <- function(sub_tree, cluster_data) {
    
    cluster_data <- cluster_data
    dup_cluster_data <- cluster_data
    full_tree <- as_tibble(sub_tree)
    
    for (i in seq_along(cluster_data)) {
      
      for (j in seq_along(dup_cluster_data)) {
        dup_cluster_data[[j]]$birth_origin <- NA
        #      dup_cluster_data[[j]]$birth_origin <- NA
        ## If the parent of origin node of one cluster (found by referencing full tree) is a child node in one of the other clusters... 
        if (isTRUE(names(cluster_data)[[i]] != names(dup_cluster_data)[[j]] &
                   (tidytree::parent(full_tree, dup_cluster_data[[j]]$parent[1])$node %in% cluster_data[[i]]$node |
                    ## Or if the parent of the parent of the origin node of one cluster (found by referencing full tree) is a child node in one of the other clusters... (meaning clusters are allowed to be separated by up to two branches)
                    tidytree::parent(full_tree, dup_cluster_data[[j]]$parent[1])$parent %in% cluster_data[[i]]$node))) {
          #        dup_cluster_data[[j]]$birth_origin == paste0("cluster_", cluster_data[[i]]$parent[1] )
          dup_cluster_data[[j]]$birth_origin = paste0(names(cluster_data)[[i]])
        } #End if statement
      } # End for loop along duplicate list
    } # End for loop along original list
    return(dup_cluster_data)
  } #End connectClust function
  cluster_data <- connectClust(sub_tree, cluster_data) 
  
  
  cluster_data <- dplyr::bind_rows(cluster_data, .id = "cluster_id") %>%
    dplyr::filter(!is.na(DATE))
  
  ## Add current date to cluster data
  cluster_data$date_added <- Sys.Date()
  
  # message("Annotating trees with cluster assignments...")
  # 
  # TreeAnno <- function(tree, clusters) {
  #   clusters.tbl <- dplyr::bind_rows(clusters, .id="cluster_id")
  #   tree.tbl <- as_tibble(tree)
  #   tree.tbl$cluster_id <- "Background"
  #   for (i in 1:nrow(clusters.tbl)) {
  #     for (j in 1:nrow(tree.tbl)) {
  #       if (isTRUE(tree.tbl$parent[j] == clusters.tbl$parent[i])) {
  #         tree.tbl$cluster_id[j] = clusters.tbl$cluster_id[i]
  #       } # End if statement
  #     } # end loop along tree
  #   } # End loop along cluster_data
  #   class(tree.tbl) = c("tbl_tree", class(tree.tbl))
  #   t2 <- as.treedata(tree.tbl)
  #   return(t2)
  # }
  # 
  # annotated_subtree <- TreeAnno(sub_tree, clusters)
  # annotated_timetree <- TreeAnno(time_tree, clusters)
  
  
  cluster_traits <- read.csv(paste0("../", list.files(path="../", pattern="trait*")[[1]]), stringsAsFactors = F) %>%
    
    cluster_traits <- rbind(cluster_traits, select(cluster_data, -parent, -node))    
  
  write("Data are now being exported as 'cluster_info_<tree>.RDS' and 'dynamite_<tree>.tree.'")
  write.csv(cluster_traits, file=paste0("../", list.files(path="../", pattern="trait*")[[1]]))
  # write.beast(annotated_subtree, paste0("dynamite_subtree_", opt$tree, ".tree")) #### NEED THIS OUTPUT####################################
  # write.beast(annotated_timetree, paste0("dynamite_timetree_", opt$tree, ".tree")) #### NEED THIS OUTPUT####################################
  
  rt1 <- Sys.time()
  rt1-rt0
  
  ## Output separate fastas for each tree
  message("Updating individual fastas...")
  
  # Read in full fasta
  fa <- read.FASTA(paste0("../samples_",Sys.Date(), "/samples_", Sys.Date(), ".fasta"))
  
  
  cluster_fa_list <- mclapply(clusters_sub_tree, function(x) {
    fa[names(fa) %in% x$tip.label]
  })
  
  for (i in seq_along(cluster_fa_list)) {
    write.FASTA(cluster_fa_list[[i]], paste0("../", new_name, "_", names(cluster_fa_list)[i], ".fasta"))
  }
  
  background_fa <- fa[names(fa) %in% background_subtree$tip.label]
  write.FASTA(background_fa, paste0(opt$path, new_name, ".fasta"))
  
}






