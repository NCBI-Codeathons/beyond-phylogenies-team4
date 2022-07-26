
rm(list=ls())

if("optparse" %in% installed.packages()) {
 library(optparse)
} else {
  install.packages(optparse)
}

.option_list = list(
  make_option(c("-m", "--metadata"), type="character", default=list.files(pattern="updated_metadata*.tab"),
              help="metadata file name [default= updated metadata file with tab extension]", metavar="character"),
make_option(c("-v", "--voc"), type="character", default=list.files(pattern="voc*"),
              help="file with list of variants of concern coinciding with metadata file  [default= file beginning with 'voc']", metavar="character")
)

.opt_parser = OptionParser(option_list=.option_list)
opt = parse_args(.opt_parser)


.packages = c("dplyr", "ape", "parallel", "foreach", "lubridate", "RColorBrewer", "ggpubr")

lapply(.packages, library, warn.conflicts=FALSE, character.only=TRUE)

`%notin%` <- Negate(`%in%`)
numCores <- detectCores()

date=Sys.Date()
# date="2022-02-11"

# Read in metadata
metadata <- read.delim(opt$metadata, sep='\t', header=T) %>%
  filter(!is.na(lineage) & lineage != "None")

for (i in seq_along(colnames(metadata))) {
  if (tryCatch({
    isTRUE(any(grepl("^id", colnames(metadata)[i], ignore.case = T)))}, error = function(e) stderr() )) {
    colnames(metadata)[i] <- "ID"
  } #end first if statement
  if (tryCatch({
    isTRUE(any(grepl("date$", colnames(metadata)[i], ignore.case = T)))}, error = function(e) stderr() )) {
    colnames(metadata)[i] <- "DATE"
  } else {
    mutate(metadata,
    Date = as.Date(gsub(".+\\|(\\d{4}-\\d{2}-\\d{2})\\|.+", "\\1", ID)))
}
  #end second if statement
}# End for loop

metadata <- metadata[!duplicated(metadata$ID),]

                  
## Define VOC (https://www.cdc.gov/coronavirus/2019-ncov/Variants/Variant-info.html) and group CA Variants (currently no Variants of high consequence)


VOC = read.delim(opt$voc, sep='\t', header=T)

## Was very specific here about VOC in metadata, but since uploading file, will need to add a bit of code to 
## merge file with metadata (and careful with nested lineages!


mutate(Variant = if_else(WHO %in% VBM | WHO %in% VOC, WHO, "Other")) 

# Set grouping for color assignment in plot
n1 <- length(VBM)
n2 <- length(VOC)

expanded_lineage_cols1 <- colorRampPalette(brewer.pal(8, "Pastel2"))(n1)
#expanded_lineage_cols2 <- colorRampPalette(brewer.pal(8, "Dark2"))(n2)
expanded_lineage_cols2 <- c("orange", "red", "darkred", "purple")


col <- setNames(c(expanded_lineage_cols1, expanded_lineage_cols2, "black"), c(VBM, VOC, "Other"))

lim_x <- c( min(metadata$Date), max(metadata$Date) )

week_calendar <- data.frame(date=seq(lim_x[1], lim_x[2], by="day")) %>%
  mutate(week=isoweek(date),year=year(date)) %>%
  group_by(year,week) %>%
  summarise(weekdate=min(date))

week_calendar <- week_calendar[-nrow(week_calendar),]


## Will need to change "Location" column to include each state as separate location but each country collectively as one location  
metadata <- metadata %>%
  mutate(day = as.numeric(1),
         week = isoweek(Date),
         month = lubridate::month(Date, label = TRUE),
         year = lubridate::year(Date),
         new_date = as.Date(paste(year,week,day), format="%Y %U %u"),
         new_date = if_else(week==53, as.Date(paste(year-1, 52, day), format="%Y %U %u"), new_date),
         new_date = if_else(week==52 & year==2022, as.Date(paste(year, 1, day), format="%Y %U %u"), new_date))
write.csv(metadata, 'metadata.csv')


#Below will need to be modified for new metadata file with all geographical locations. And .m1 above will ned to be revised to include geographical location info
locales = group_by(metadata, location) %>%
  group_split()

m1 <- lapply(locales, function(x) {
  x <- group_by(x, new_date) %>%
    mutate(total_date =  length(taxon)) %>%
    ungroup() %>%
    group_by(new_date, Variant) %>%
    mutate(total_prev=length(taxon)) %>%
    summarize(perc=total_prev/total_date) %>%
    distinct()
  return(x)
})

counts <- lapply(locales, function(x) {
  group_by(x, new_date) %>%
    summarize(total_date =  length(taxon)) %>%
    distinct()
})


unique_dates <- unique(sort(metadata$new_date))

rect_left <- seq(min(unique_dates), max(unique_dates), 7)
rect_left <- rect_left[seq(1, length(rect_left), 2)]
rectangles <- data.frame(
  xmin = rect_left,
  xmax = rect_left + 7,
  ymin = 0,
  ymax = Inf
)


plot_list <- list()
scale_factor <- list() 

#invisible(foreach(x=seq_along(m1)) %do% {
  x=1
  scale_factor[[x]] <- max(counts[[x]]$total_date)
  
#  plot_list[[x]] <- ggplot() +
  ggplot(data=m1[[x]]) +
    #was using geom_area instead before, but not as nice
    geom_col(data=m1[[x]], aes(x = new_date, fill = Variant, y = perc), position="stack", color="white") +
    xlab("Time") +
    scale_fill_manual(values=col, name="Lineage", breaks=names(col)) + 
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
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.text.x = element_text(size = 12, angle = 90),
          axis.title.x.bottom = element_text(size = 16),
          axis.text.y.left = element_text(size = 12),
          axis.title.y.left = element_text(size = 16),
          axis.title.y.right = element_text(size = 16))
  
#})

# p1 <- ggarrange(plotlist=plot_list,
#                 ncol=1, 
#                 labels = c("Alachua", "Rest of Florida", "Miami"),
#                 common.legend=T,
#                 legend="right")

ggsave(plot=last_plot(), file=paste0("Alachua_Lineages_", date, ".pdf"), width=11, height=8, units="in")

x=2
scale_factor[[x]] <- max(counts[[x]]$total_date)

#  plot_list[[x]] <- ggplot() +
ggplot(data=m1[[x]]) +
  #was using geom_area instead before, but not as nice
  geom_col(data=m1[[x]], aes(x = new_date, fill = Variant, y = perc), position="stack", color="white") +
  xlab("Time") +
  scale_fill_manual(values=col, name="Lineage", breaks=names(col)) + 
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
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x.bottom = element_text(size = 16),
        axis.text.y.left = element_text(size = 12),
        axis.title.y.left = element_text(size = 16),
        axis.title.y.right = element_text(size = 16))

ggsave(plot=last_plot(), file=paste0("Florida_Lineages_", date, ".pdf"), width=11, height=8, units="in")

x=3
scale_factor[[x]] <- max(counts[[x]]$total_date)

#  plot_list[[x]] <- ggplot() +
ggplot(data=m1[[x]]) +
  #was using geom_area instead before, but not as nice
  geom_col(data=m1[[x]], aes(x = new_date, fill = Variant, y = perc), position="stack", color="white") +
  xlab("Time") +
  scale_fill_manual(values=col, name="Lineage", breaks=names(col)) + 
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
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.title.x.bottom = element_text(size = 16),
        axis.text.y.left = element_text(size = 12),
        axis.title.y.left = element_text(size = 16),
        axis.title.y.right = element_text(size = 16))

ggsave(plot=last_plot(), file=paste0("Miami_Lineages_", date, ".pdf"), width=11, height=8, units="in")

#Report other Variants
unique(metadata$lineage[metadata$Date>="2021-12-01" & metadata$Variant %notin% c(VOC, VBM)])
table(metadata$lineage[metadata$Date>="2021-12-01" & metadata$Variant %notin% c(VOC, VBM)])

metadata[metadata$Date>="2021-12-01" & metadata$Variant %notin% c(VOC, VBM) & metadata$location=="Alachua",] %>%
  group_by(lineage) %>%
  summarize(n=length(taxon)) %>%

  ggplot(aes(x=2, y=n, fill=lineage)) +
  geom_col() +
  coord_polar(theta="y", start=0) +
  xlim(0.5,2.5) +
  theme_void() # remove background, grid, numeric labels

ggsave(plot=last_plot(), paste0("Non-VOCs_", date, ".pdf"), width=4, height=4, units="in")
length(metadata$taxon[metadata$WHO=="Omicron"])


