
rm(list=ls())
.packages = c("dplyr", "ape", "parallel", "foreach", "lubridate", "RColorBrewer", "ggpubr")

lapply(.packages, library, warn.conflicts=FALSE, character.only=TRUE)

`%notin%` <- Negate(`%in%`)
numCores <- detectCores()

date=Sys.Date()
# date="2022-02-11"

# Read in metadata
metadata <- read.csv(paste0("updated_lineages_", date, ".csv"), header=T) %>%
  filter(!is.na(lineage) & lineage != "None") %>%
  filter(!grepl("Haiti", taxon)) %>%
#  mutate(location=if_else(location=="Global","FL", location)) %>%
  select(taxon, lineage) %>%
  mutate(Date = as.Date(gsub(".+\\|(\\d{4}-\\d{2}-\\d{2})\\|.+", "\\1", taxon), 
                        format="%Y-%m-%d")) %>%
  mutate(location=gsub("hCoV-19/USA/FL-([A-Za-z]+)-.+", "\\1", taxon)) %>%
  mutate(location=if_else(grepl("shands|path|STP|UF|SED", location, ignore.case=T), "Alachua", 
  	if_else(grepl("miami|msl", location, ignore.case=T), "Miami", "FL")))
  

metadata <- metadata[!duplicated(metadata$taxon),]

                  
## Define VOC (https://www.cdc.gov/coronavirus/2019-ncov/Variants/Variant-info.html) and group CA Variants (currently no Variants of high consequence)

#VBM <- c("B.1.1.7/Q", "B.1.351", "P.1", "B.1.427/B.1.429", "B.1.525", "B.1.526", "B.1.617.1", "B.1.617.3", "B.1.621", "P2")
#VOC <- c( "B.1.617.2/AY", "B.1.1.529")

VBM <- c("20E (EU1)","Alpha", "Beta", "Gamma", "Epsilon", "Eta", "Iota", "Kappa", 
  "B.1.617.3", "Mu", "Zeta")
VOC <- c( "Delta", "Omicron (Original)", "Omicron (BA.1)", "Omicron (BA.2)")
metadata <- mutate(metadata,
                   WHO = if_else(grepl("B.1.177|B.1.177.73", lineage), "20E (EU1)",
                                  if_else(grepl("B.1.1.7|Q", lineage), "Alpha",
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
                                                                                                                           if_else(grepl("B.1.1.529", lineage), "Omicron (Original)",
                                                                                                                                   if_else(grepl("BA.1", lineage), "Omicron (BA.1)",
                                                                                                                                           if_else(grepl("BA.2", lineage), "Omicron (BA.2)",
                                                                                                                                   "Other")))))))))))))))) %>%
  mutate(Variant = if_else(WHO %in% VBM | WHO %in% VOC, WHO, "Other")) 

# Set grouping for color assignment in plot
n1 <- length(VBM)
n2 <- length(VOC)

expanded_lineage_cols1 <- colorRampPalette(brewer.pal(8, "Pastel2"))(n1)
#expanded_lineage_cols2 <- colorRampPalette(brewer.pal(8, "Dark2"))(n2)
expanded_lineage_cols2 <- c("orange", "red", "darkred", "purple")


col <- setNames(c(expanded_lineage_cols1, expanded_lineage_cols2, "black"), c(VBM, VOC, "Other"))

## Filtering base on just 2021
metadata <- filter(metadata, Date >= "2021-01-01")

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


