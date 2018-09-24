library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(Hmisc)



#data from: https://portal.lternet.edu/nis/mapbrowse?packageid=knb-lter-ntl.262.2

#Density here is a somewhat specialized data type, in that it is based off counts from 
#subsamples, but then scaled to the full sample and rounded to the nearest thousand.
#This makes it very difficult to treat it as either a count type or a continous
#zero-inflated variable. Here it likely makes sense to treat it basically as a Tweedie
#distribution, with some confusing caveats...


abbr_first_word = function(x) {
  split_names = str_split_fixed(x,pattern = " ", n = 2) #split names at a space, into two words
  first_letter = paste(str_sub(split_names[,1],start = 1,end = 1),".",sep = "") #Turn the first word into a single letter
  paste(first_letter,split_names[,2],sep =  " ") #Re-paste the words together
}

abbr_first_word = function(x) {
  split_names = str_split_fixed(x,pattern = " ", n = 2) #split names at a space, into two words
  first_letter = paste(str_sub(split_names[,1],start = 1,end = 1),".",sep = "") #Turn the first word into a single letter
  paste(first_letter,split_names[,2],sep =  " ") #Re-paste the words together
}

zooplankton_data = read.csv("data/madisonlakeszoopoldnet.csv",stringsAsFactors = F)%>%
  filter(!taxon  %in% c("DAPHNIA MENDOTAE JUVENILE", #excluding juvenile stages
                        "DAPHNIA PULICARIA JUVENILE",
                        "DAPHNIA RETROCURVA JUVENILE",
                        "DIAPHANOSOMA BIRGEI JUVENILE",
                        "COPEPOD NAUPLII"),
         density !=43 #filtering out the single observation below 1000; very likely an error (i.e. not scaled properly)
         )%>%
  mutate(date = lubridate::as_date(sampledate), #adjusting labels
         day  = yday(date),
         taxon = tolower(taxon),
         taxon = Hmisc::capitalize(taxon),
         lake  = case_when(lakeid=="ME"~"Mendota",
                           lakeid=="MO"~"Menona",
                           lakeid=="KE"~"Kegonsa",
                           lakeid=="WA"~"Waubesa")
         )%>%
  rename(year = year4)%>%
  dplyr::select(taxon,lake, year, day, density)%>% #filtering out extraneous columns
  group_by(taxon)%>%
  filter(n()>400)%>% #filtering to only keep the most common species
  ungroup()%>%
  mutate(taxon = case_when(taxon=="Calanoida copepodites"~"Calanoid copepods",
                           taxon=="Cyclopoida copepodites"~"Cyclopoid copepods",
                           TRUE~abbr_first_word(taxon)))%>%
  complete(nesting(day, year,lake), 
           taxon, 
           fill = list(density=NA))%>%
  mutate(present = !is.na(density))%>%
  group_by(taxon)%>%
  mutate(
    density = ifelse(is.na(density), 0,density),
    density_adj = density + 1000,
    min_density = min(density[density>0]))%>%
  group_by(taxon,year,lake)%>%
  filter(any(present))%>%
  mutate(density_scaled = scale(log10(density+min_density)))%>%
  ungroup()%>%
  select(-present)

write.csv(zooplankton_data, "data/zooplankton_example.csv",row.names = FALSE)