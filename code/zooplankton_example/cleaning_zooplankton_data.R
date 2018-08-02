library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(Hmisc)


#data from: https://portal.lternet.edu/nis/mapbrowse?packageid=knb-lter-ntl.262.2

#Density here is a somewhat specialized data type, in that it is based off counts from 
#subsamples, but then scaled to the full sample and rounded to the nearest thousand.
#This makes it very difficult to treat it as either a count type or a continous
#zero-inflated variable. I've chosen to deal with this by replacing any zero
#value with the minumum observed density, then log-transforming and scaling within 
#species and year. 

zooplankton_data = read.csv("data/madisonlakeszoopoldnet.csv",stringsAsFactors = F)%>%
  filter(!taxon  %in% c("DAPHNIA MENDOTAE JUVENILE", #excluding juvenile stages
                        "DAPHNIA PULICARIA JUVENILE",
                        "DAPHNIA RETROCURVA JUVENILE",
                        "DIAPHANOSOMA BIRGEI JUVENILE",
                        "COPEPOD NAUPLII"))%>%
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
  mutate(density = ifelse(is.na(density), min(density,na.rm = T),density))%>%
  group_by(taxon,year,lake)%>%
  filter(any(present))%>%
  mutate(density_scaled = scale(log10(density)))%>%
  ungroup()%>%
  select(-present)

write.csv(zooplankton_data, "data/zooplankton_example.csv",row.names = FALSE)