# we will assume birds range from latitude 0 to latitude 60,
# with birds varying in in the peak times they reach max and min latitude, 
# and in time they spend migrating vs. at breeding/feeding grounds 
# modelled as differences in the shape of the migration curve

library(dplyr)
library(tidyr)
set_seed = 10
species_data = data_frame(species = factor(paste("sp", 1:8, sep="")),
                          phase_shift = rnorm(8,0,6),
                          curve_shape = rnorm(8,-.25,0.4),
                          curve_width = 3
                          )

bird_move = crossing(latitude = seq(5,60,length=10),
                     species = species_data$species,
                     week = seq(4,52, by=6), n_indiv = 100)%>%
  left_join(species_data)%>%
  mutate(peak_lat = -cos((week+phase_shift)*(2*pi/52)),
         peak_lat = peak_lat*abs(peak_lat)^curve_shape,
         peak_lat = (peak_lat+1)*30,
         avg_abund = dnorm(latitude,peak_lat, curve_width))%>%
  group_by(species, week) %>%
  mutate(count = rmultinom(1,n_indiv[1],prob = avg_abund)[,1])

write.csv(bird_move,file = "data/bird_move.csv",row.names = F)
                          