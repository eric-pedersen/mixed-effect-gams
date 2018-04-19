library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)
source("../code/functions.R")

get_n_pen  = function(model) {
  family = model$family[[1]]
  if(family %in% c("Gamma","gaussian")){
    capture.output({out_val = nrow(gam.vcomp(model))-1})
  }else{
    capture.output({out_val = nrow(gam.vcomp(model))})
  }
  return(out_val)
}

get_n_coef = function(model) length(coef(model))

get_n_iter = function(model) model$outer.info$iter
get_n_out_iter = function(model) model$iter
  
bird_move = read.csv("../data/bird_move.csv")

comp_resources = crossing(model_number = c("1","2","3","4","5"),
                       data_source = factor(c("CO2","bird_move"),
                                     levels = c("CO2","bird_move")),
                       time = 0, n_smooths = 0,
                       n_coef = 0)

CO2$Plant_uo = factor(CO2$Plant, levels = levels(CO2$Plant), ordered = F)


comp_resources[1,"time"] = system.time(CO2_mod1 <- bam(log(uptake) ~ s(log(conc),k=5,m=2, bs="tp")+s(Plant_uo, k =12,  bs="re"), 
                                                    data= CO2,method="REML",
                                                    control = list(keepData=TRUE)))[3]

comp_resources[2,"time"] = system.time(bird_mod1 <- bam(count ~ te(week,latitude, bs= c("cc", "tp"), k=c(10,10)), 
                                                     data= bird_move, method="REML", family= poisson,
                                                     control = list(keepData=TRUE)))[3]

comp_resources[3,"time"] = system.time(CO2_mod2 <- bam(log(uptake) ~ s(log(conc),k=5,m=2, bs="tp")+
                                                       s(log(conc), Plant_uo, k=5, bs="fs",m=1), 
                                                     data= CO2,method="REML",
                                                     control = list(keepData=TRUE)))[3]


comp_resources[4,"time"] = system.time(bird_mod2 <- bam(count ~ te(week,latitude, bs= c("cc", "tp"), 
                                                                 k=c(10,10),m=c(2,2))+
                                                        te(week,latitude,species, bs= c("cc", "tp","re"), 
                                                           k=c(10,10,6),m = c(1,1,1)), 
                                                      data= bird_move, method="REML", family= poisson,
                                                      control = list(keepData=TRUE)))[3]


comp_resources[5,"time"] = system.time(
  CO2_mod3 <- bam(log(uptake) ~ s(log(conc),k=5,m=2, bs="tp")+
                    s(log(conc),by= Plant_uo, k =5,  bs="ts",m=1)+
                    s(Plant_uo,bs="re",k=12), 
                  data= CO2,method="REML",
                  control = list(keepData=TRUE)))[3]



comp_resources[6,"time"] = system.time(
  bird_mod3 <- bam(count ~ te(week,latitude, bs= c("cc", "tp"), 
                              k=c(10,10),m=c(2,2)) +
                     te(week,latitude, bs= c("cc", "tp"), 
                        k=c(10,10),m=c(1,1),by= species), 
                   data= bird_move, method="REML", family= poisson,
                   control = list(keepData=TRUE)))[3]


comp_resources[7,"time"] = system.time(
  CO2_mod4 <- bam(log(uptake) ~ s(log(conc), Plant_uo, k=5,  bs="fs",m=2), 
                  data= CO2,method="REML",
                  control = list(keepData=TRUE)))[3]


comp_resources[8,"time"] = system.time(
  bird_mod4 <- bam(count ~ te(week,latitude,species, bs= c("cc", "tp","re"), 
                              k=c(10,10,6),m = 2), 
                   data= bird_move, method="REML", family= poisson,
                   control = list(keepData=TRUE))
)[3]


comp_resources[9,"time"] = system.time(
  CO2_mod5 <- bam(log(uptake) ~ s(log(conc),by= Plant_uo, k =5,  bs="tp",m=2)+
                    s(Plant_uo,bs="re",k=12), data= CO2,method="REML",
                  control = list(keepData=TRUE))
  
)[3]

comp_resources[10,"time"] = system.time(
  bird_mod5 <- bam(count ~ te(week,latitude,by=species, bs= c("cc", "tp"), 
                              k=c(10,10),m = 2), 
                   data= bird_move, method="REML", family= poisson,
                   control = list(keepData=TRUE))
)[3]


comp_resources$model = list(CO2_mod1, bird_mod1, CO2_mod2, bird_mod2,
                            CO2_mod3, bird_mod3,CO2_mod4, bird_mod4, 
                            CO2_mod5, bird_mod5)

comp_resources = comp_resources %>%
  group_by(model_number, data_source)%>%
  mutate(n_smooths = get_n_pen(model[[1]]),
         n_coef = get_n_coef(model[[1]]),
         n_iter = get_n_iter(model[[1]]),
         n_iter_out = get_n_out_iter(model[[1]]))



