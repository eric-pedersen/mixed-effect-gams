# we will assume birds range from latitude 0 to latitude 60,
# with birds varying in in the peak times they reach max and min latitude, 
# and in time they spend migrating vs. at breeding/feeding grounds 
# modelled as differences in the shape of the migration curve

library(dplyr)
library(tidyr)
set.seed(12)
n_sp = 6
species_data = data_frame(species = factor(paste("sp", 1:n_sp, sep="")),
                          phase_shift = sort(rnorm(n_sp,0,5)),
                          curve_shape = sort(rnorm(n_sp,-.1,0.4)),
                          curve_width =6
                          )

bird_move = crossing(latitude = seq(5,60,length=10),
                     species = species_data$species,
                     week = seq(0,52, by=4), 
                     n_indiv = 100,
                     obs_prob = 0.2)%>%
  left_join(species_data)%>%
  mutate(peak_lat = -cos((week+phase_shift)*(2*pi/52)),
         peak_lat = peak_lat*abs(peak_lat)^curve_shape,
         peak_lat = (peak_lat+1)*30,
         avg_abund = dnorm(latitude,peak_lat, curve_width))%>%
  group_by(species, week) %>%
  mutate(count = rmultinom(1,n_indiv[1],prob = avg_abund)[,1],
         count = rbinom(n(),count, obs_prob))%>%
  select(latitude, week, species, count)

write.csv(bird_move,file = "data/bird_move.csv",row.names = F)
                          