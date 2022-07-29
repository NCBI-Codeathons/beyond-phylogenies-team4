#GENERATE VARIABLES FOR TESTING

pacman::p_load(dplyr, lubridate, stringr, readr)

dir = "example_files/updated_metadata_2021-03-17.tab"
  
temp = read_delim(dir, delim = "\t" , na = c(" ", "", "NA"), show_col_types = FALSE)

#if there are duplicate columns delete the 2nd one
i = which(stringr::str_detect(colnames(temp), regex("lineage", ignore_case = T)))

if(length(i) >1) {
  temp =  temp[,-i[[2]]]
} 

#cluster_list = c("c1", "c2", "background")
florida_counties = c("Citrus", "Miami-Dade", "Gainsville")
vac_list = c("Y", "N")
sex_list = c("M", "F")

keep_col = c("ID", "DATE", "lineage", "location","age", "sex", "vaccinated")

temp = temp %>%
  select(-Location) %>% 
  mutate(location = sample(florida_counties, nrow(temp), replace = T),
      #   cluster = sample(cluster_list, nrow(temp), replace = T),
         age = round(rnorm(nrow(temp), mean = 40, sd =  10)),
         vaccinated = sample(vac_list, nrow(temp), replace = T),
         sex = sample(sex_list, nrow(temp), replace = T)
        ) %>%
  select(all_of(keep_col))
      
#sprinkle some NAs in there
s = round(runif(50, 1, max(nrow(temp))))
s = unique(s)
temp[s,"age"] = NA

s = round(runif(50, 1, max(nrow(temp))))
s = unique(s)
temp[s,"location"] = NA

s = round(runif(50, 1, max(nrow(temp))))
s = unique(s)
temp[s,"vaccinated"] = NA

s = round(runif(50, 1, max(nrow(temp))))
s = unique(s)
temp[s,"sex"] = NA

write_delim(temp, "updated_metadata_test.tab", delim = "\t")

rm(florida_counties, 
   vac_list, sex_list, temp, dir, i, keep_col)