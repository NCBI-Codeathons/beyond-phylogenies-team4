#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#HEADER CODE ####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


rm(list=ls())

if("optparse" %in% installed.packages()) {
 library(optparse)
} else {
  install.packages(optparse)
}

#potentially add as argument or determine naming convention

#specify location of files
.option_list = list(
  make_option(c("-m", "--metadata"), type="character",
              default=list.files(pattern =  "updated_metadata"), #example_files would nbeed to be changes to something permanent
              help="metadata\ file name [default= updated metadata file with tab extension]", metavar="character"),
  make_option(c("-v", "--voc"), type="character", 
              default=list.files(pattern="voc*"),
              help=" user input file with list of lineages of concern coinciding with metadata file  [default= file beginning with 'voc']", metavar="character")
                  )

.opt_parser = OptionParser(option_list=.option_list)
opt = parse_args(.opt_parser)


#could use pacman to change to include require
.packages = c("dplyr", "ape", "parallel", "foreach", "lubridate", "RColorBrewer", "ggpubr",
              "stringr", "tidyr", "readr") #new

lapply(.packages, library, warn.conflicts=FALSE, character.only=TRUE)

`%notin%` <- Negate(`%in%`)
numCores <- detectCores()

date=Sys.Date()
# date="2022-02-11"

# END OF HEADER CODE ####

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# READ DATA ####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#allows metadata to be tab delimitted or csv
if(str_detect(opt$metadata, pattern = ".tab$")) { #if file ends in csv read.csv else read.delim
  metadata = read_delim(opt$metadata, delim = "\t" , na = c(" ", "", "NA"), show_col_types = FALSE)
} else {
  metadata <- read.csv(opt$metadata, header=T) 
}

#FOR VOC
temp = unique(metadata$lineage)
  
#if no files with voc pattern use metadata with all lineages else use user input in voc.tab file
if(length(opt$voc) == 0) { 
  VOC = data.frame(lineage = temp,
                   WHO = rep("unknown", length(temp)))
} else {
  VOC = read_delim(opt$voc, delim = "\t" , na = c(" ", "", "NA"), show_col_types = FALSE)
}

#NEED TO BE ADDED AS INPUT OR STATIC DIRECTORY IS DETERMINED!!!
keep_col = c("City", "Age", "Sex", "Vaccinated", "Location")

traits = read.csv("example_files/trait_distributions_2021-02-10.csv") %>%
  filter(field %in% keep_col)

rm(temp)

# END OF READ DATA ####

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#PREPROCESSING CHECKS #### 
#clean columns of any special characters and pattern match and rename columns
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#remove any whitespace
colnames(metadata) = trimws(colnames(metadata))

#change to all title
colnames(metadata) = str_to_lower(colnames(metadata))

#replace all special characters with "_"
colnames(metadata) = str_replace(string = colnames(metadata), 
                                pattern = "[[:punct:]]",  #make all special characters
                                replacement = "_")

#COMMENTED OUT SECTION NOT WORKING####
#could get function to overwrite column names

# df = metadata
# new_name = "test"
# pattern = "^date$"
# rm(df, new_name, pattern)
# 
# col_fun = function(df, new_name, pattern) {
# 
#   i = which(stringr::str_detect(colnames(df), regex(pattern, ignore_case = T)))
# 
#   if(length(i) == 0) {
#     print(paste0("You are missing a ", new_name, " column for your samples or it doesn't match list of plauseable
#         matches. Please include a ", new_name, " column."))
#   } else {
#     colnames(df)[i] = new_name
#   }
#   
# }
# 
# col_fun(df = metadata,
#         new_name = "test",
#         pattern = "^date$")

# END OF FUNCTION####



#duplicate check
if(length(unique(metadata$id)) < nrow(metadata))


#DATE
#finds colnames with pattern date and replaces it with consistent format
new_name = "date"

i = which(stringr::str_detect(colnames(metadata), regex("^date$", ignore_case = T)))

  if(length(i) == 0) {
    print(paste0("You are missing a ", new_name, " column for your samples or it doesn't match list of plauseable
        matches. Please include a ", new_name, " column."))
  } else {
    colnames(metadata)[i] = new_name
  }


#ID
#finds colnames with pattern ID and replaces it with consistent format
new_name = "id"

i = which(stringr::str_detect(colnames(metadata), regex("^id", ignore_case = T)))

if(length(i) == 0) {
  print(paste0("You are missing a ", new_name, " column for your samples or it doesn't match list of plauseable
        matches. Please include a ", new_name, " column."))
} else {
  colnames(metadata)[i] = new_name
}

rm(i)



#Read in user input VOC 

if(length(list.files(pattern = "voc*")) > 0) { #does voc file exist?
  VOC = read.delim(opt$voc, sep='\t', header=T) #if true write in the file
} else { #if not keep all lineages and add WHO unknown
  
  print("if you wish to filter for specific variants with WHO lineages, create a voc.tab file. You can find the name from the following url 
         https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-classifications.html?CDC_AA_refVal=https%3A%2F%2Fwww.cdc.gov%2Fcoronavirus%2F2019-ncov%2Fvariants%2Fvariant-info.html")
  VOC = data.frame(lineage = unique(metadata$lineage),
             WHO = rep("unknown", length(unique(metadata$lineage)))
             )
}

# END OF PREPROCESSING CHECKS ####

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DATA CLEAN ####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

filter_col = c("location", "lineage", "date", "cluster", "WHO") #columns to drop if any are missing

temp = traits %>% distinct(cluster_id, ID)

## Will need to change "Location" column to include each state as separate location but each country collectively as one location  
metadata2 <- metadata %>%
  left_join(temp, by = c("id" = "ID")) %>%
  left_join(VOC, by = "lineage") %>% #will filter out lineages that dont match
  mutate(new_date = floor_date(date, "week", 1) + 7 , #rounds to the previous 
         age_group = cut(age, seq(0,100,25), 
                                c("age 0-25", "age 26-50", "age 51-75", 
                                  "76-100"))
         )
metadata2[is.na(metadata2)] <- "unknown"

write.csv(metadata2, 'metadata.csv')



#Below will need to be modified for new metadata2 file with all geographical locations. And .m1 above will ned to be revised to include geographical location info
locales = group_by(metadata2, location) %>%
  group_split()

locale_nm = unique(metadata2$location)

m1 <- lapply(locales, function(x) {
  x <- group_by(x, new_date) %>%
    mutate(total_date =  length(lineage)) %>%
    ungroup() %>%
    group_by(new_date, lineage) %>%
    mutate(total_prev=length(lineage)) %>%
    summarize(perc=total_prev/total_date) %>%
    distinct()
  return(x)
})

counts <- lapply(locales, function(x) {
  group_by(x, new_date) %>%
    summarize(total_date =  length(lineage)) %>%
    distinct()
})


#END OF DATA CLEAN ####

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PLOT PREANALYSIS ####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#create weekly rectangle shading objects
unique_dates <- unique(sort(metadata2$new_date))

rect_left <- seq(min(unique_dates), max(unique_dates), 7)
rect_left <- rect_left[seq(1, length(rect_left), 2)]
rectangles <- data.frame(
  xmin = rect_left,
  xmax = rect_left + 7,
  ymin = 0,
  ymax = Inf
)


#generate colors for plots
lineage_list = unique(metadata2$lineage)
col_list <- colorRampPalette(brewer.pal(8, "Pastel2"))(length(lineage_list))

# Set grouping for color assignment in plot
lineage_cols = data.frame(lineage = lineage_list,
                          col = col_list)

lim_x <- c( min(metadata2$date), max(metadata2$date))

# END OF PLOT PREPROCESSING ####

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LINEAGE PLOTTING ####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#create empty list for the for loop
plot_list <- list()
scale_factor <- list() 


for(x in seq_along(m1)) {

scale_factor[[x]] <- max(counts[[x]]$total_date)

#filter out lineages that aren't in the location while keeping the color scheme consistent
plot_col = lineage_cols %>% filter(lineage %in% m1[[x]]$lineage)
  
plot_list[[x]] <- 
ggplot(data=m1[[x]]) +
  #was using geom_area instead before, but not as nice
  geom_col(aes(x = new_date, fill = lineage, y = perc), position="stack", color="white") +
  xlab("Time") +
  scale_fill_manual(values= plot_col$col , name="Lineage", breaks= plot_col$lineage) + 
  theme_minimal() + coord_cartesian(xlim = lim_x) +
  scale_x_date(date_breaks="1 month", date_labels = "%b") +
  geom_line(data = counts[[x]], aes(x=new_date, y=total_date/scale_factor[[x]]), linetype="dashed") +
  scale_y_continuous(name = "% of Samples", labels=scales::percent,
                     sec.axis = sec_axis(trans=~.*scale_factor[[x]], name="Number of Samples", breaks=seq(0,scale_factor[[x]],200)),
                     limits=c(0,1.20),
                     breaks=seq(0,1,0.25)) +
  geom_rect(data=rectangles, aes(xmin=xmin-3.5, xmax=xmax-3.5, ymin=ymin, ymax=ymax), 
            fill='gray80', alpha=0.25) +
  annotate(geom="text", x=rectangles$xmin, y=1.125, size = 2, angle=90, hjust=.5, vjust=.5, label=rect_left) +
  ggtitle(locale_nm[[x]]) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x.bottom = element_text(size = 16),
        axis.text.y.left = element_text(size = 12),
        axis.title.y.left = element_text(size = 16),
        axis.title.y.right = element_text(size = 16))

}

#WRITE PLOT

if(dir.exists("plot_results")) {
  dir = "plot_results/"
} else {
  dir.create("plot_results")
  dir = "plot_results/"
}


for(x in seq_along(m1)) {
  
  ggsave(plot= plot_list[[x]], file=paste0(dir, locale_nm[[x]], "_","lineage_plot_", date, ".pdf"), width=11, height=8, units="in")
}



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# METADATA PLOTTING ####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

cluster_nm = unique(metadata2$cluster_id)
  
cluster_df = metadata2 %>% 
  group_by(cluster_id,new_date, location, age_group, vaccinated) %>%
  summarize(n = n()) %>%
  ungroup()


cluster_df = cluster_df%>%
  group_by(cluster_id) %>%
  group_split()


cluster_list = list()

for(x in seq_along(cluster_df)){
  
cluster_list[[x]] <-
  ggplot(cluster_df[[x]], aes(x = new_date, y = n, fill = vaccinated, color = vaccinated)) +
    geom_col(alpha=0.8) +
    facet_grid(rows = vars(age_group),
               cols = vars(location)) +
    theme_minimal() +
    ggtitle(cluster_nm[[x]]) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
}

for(x in seq_along(cluster_df)) {
  
  ggsave(plot= cluster_list[[x]], file=paste0(dir, cluster_nm[[x]], "_","cluster_traits_", date, ".pdf"), width=11, height=8, units="in")
}

