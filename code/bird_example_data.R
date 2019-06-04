# we will assume birds range from latitude 0 to latitude 60,
# with birds varying in in the peak times they reach max and min latitude, 
# and in time they spend migrating vs. at breeding/feeding grounds 
# modelled as differences in the shape of the migration curve

library(dplyr)
library(tidyr)

#Check for the R version number, and if greater than 3.6, switch to using the
#old random number generator, to ensure replicability

if(getRversion()>= 3.6) RNGversion("3.5.0")

set.seed(12)
n_sp = 6
species_data = data_frame(species = factor(paste("sp", 1:n_sp, sep="")),
                          phase_shift = sort(rnorm(n_sp,0,5)),
                          curve_shape = sort(rnorm(n_sp,-.1,0.4)),
                          curve_width =6
                          )

bird_move_simulations = crossing(latitude = seq(5,60,length=10),
                     species = species_data$species,
                     week = seq(0,52, by=4), 
                     n_indiv = 100,
                     obs_prob = 0.2)%>%
  left_join(species_data)%>%
  mutate(#This creates a phase-shifted version of a single cosine curve for each
         #species. This determines the latitude that, in any given week, 
         #the most birds of that species are present in
         peak_lat = -cos((week+phase_shift)*(2*pi/52)),
         #This changes the shape of the migration peak; flatter curves 
         #(= smaller shape values) correspond to more time spend in high and 
         #low latitudes (and thus a faster transition between summer and winter
         #grounds
         peak_lat = peak_lat*abs(peak_lat)^curve_shape,
         #scale peak latitude so it ranges from 0 to 30
         peak_lat = (peak_lat+1)*30,
         #transforms peak latitude for each week into a normal distribution
         #around that peak value
         avg_abund = dnorm(latitude,peak_lat, curve_width))%>%
  group_by(species, week) %>%
  #This samples from the expected distribution of birds in each site by first
  #distributing the tagged birds across their range with a multinomial, then
  #using a binomial distribution to choose a fraction of those
  mutate(count = rmultinom(1,n_indiv[1],prob = avg_abund)[,1],
         count = rbinom(n(),count, obs_prob))

#Extract the implied global function, by averagin across species realizations
bird_move_global = bird_move_simulations %>%
  group_by(week, species)%>%
  mutate(avg_abund = n_indiv*obs_prob*avg_abund/sum(avg_abund))%>%
  group_by(latitude, week)%>%
  summarize(global_function = mean(log(avg_abund)),
            n_indiv = n_indiv[1],
            obs_prob = obs_prob[1])%>%
  ungroup()%>%
  mutate(global_scaled_function = global_function-mean(global_function))
  
bird_move = bird_move_simulations %>%
  select(latitude, week, species, count)

write.csv(bird_move,file = "data/bird_move.csv",row.names = F)
write.csv(bird_move_global,file = "data/bird_move_global.csv",row.names = F)   
