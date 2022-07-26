require(ape)
require(optparse)

option_list = list(
  make_option(c("-t", "--tree"), type="character", default="uncondensed-final-tree.nh", 
              help=" tree file [default= .nh]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

tree <- read.tree(opt$tree)
tree <- multi2di(tree, random = TRUE, equiprob = TRUE)
tree$edge.length[is.nan(tree$edge.length)] <- 0

write.tree(tree, file=(paste0(opt$tree, ".tree")))
