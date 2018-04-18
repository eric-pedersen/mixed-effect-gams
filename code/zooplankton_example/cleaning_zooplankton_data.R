library(dplyr)
library(tidyr)
library(lubridate)
library(Hmisc)

#data from: https://portal.lternet.edu/nis/mapbrowse?packageid=knb-lter-ntl.262.2

zooplankton_data = read.csv("data/madisonlakeszoopoldnet.csv",stringsAsFactors = F)%>%
  filter(lakeid == "ME",
         !taxon  %in% c("DAPHNIA MENDOTAE JUVENILE", 
                        "DAPHNIA PULICARIA JUVENILE",
                        "COPEPOD NAUPLII"))%>%
  mutate(date = lubridate::as_date(sampledate),
         day  = yday(date),
         taxon = tolower(taxon),
         taxon = Hmisc::capitalize(taxon)
         )%>%
  rename(year = year4)%>%
  dplyr::select(taxon, year, day, density, avg_length)%>%
  group_by(taxon)%>%
  filter(n()>150)%>%
  ungroup()%>%
  complete(nesting(day, year), 
           taxon, 
           fill = list(density=NA, avg_length=NA))%>%
  mutate(present = as.numeric(!is.na(avg_length)))%>%
  group_by(taxon)%>%
  mutate(density = ifelse(is.na(density), min(density,na.rm = T),density))%>%
  group_by(taxon,year)%>%
  filter(any(present))%>%
  mutate(density_scaled = scale(log10(density)))%>%
  ungroup()

write.csv(zooplankton_data, "data/zooplankton_example.csv")