library(mgcv)
library(MASS)
library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)

source("code/functions.R")

#starting parameters  ####  
n_groups = 12 # number of groups
n_replicates = 20

total_amp = 1
noise  = total_amp/3 #variance of random noise around main function

# The length-scale parameter for the main and individual-level functions. Points
# with x-values that are more distant from each other than the length scale
# will be only weakly correlated
main_func_scale= 0.1
indiv_func_scale_base = 0.1

sim_data = rbind(
  expand.grid(trial = "varying main amp",rep = 1:n_replicates,
              n_data = c(50,100,200),
              main_func_amp = c(0,0.25,0.5,0.75,1.0),
              indiv_scale_diff = c(1,2)),
  expand.grid(trial = "varying group scale",rep = 1:n_replicates,
              n_data = c(50,100,200),
              main_func_amp = 0.5,
              indiv_scale_diff = c(1,2,4,8)
              ))
sim_data = data.frame(sim_data)
sim_data$seed = 1:nrow(sim_data)

aic_data = mutate(sim_data,
                  model_1 = 0,model_2a = 0,model_2b=0,model_3 =0,
                  model_4a= 0,model_4b = 0,model_5 =0)

for(i in 1:nrow(sim_data)){
  set.seed(sim_data$seed[i])
  n_data = sim_data$n_data[i]
  x = seq(0,1, length=n_data)
  indiv_scale_diff = sim_data$indiv_scale_diff[i]
  # variability of the main and individual level functions. Equal to the variance
  # of individual points drawn from the function at large distances from one another.
  # set so that the variance of the main and individual level-functions will sum
  # to a fixed value
  main_func_amp = sim_data$main_func_amp[i]
  indiv_func_amp = total_amp-main_func_amp #var

  indiv_func_scale = indiv_func_scale_base*indiv_scale_diff^rep(c(-1,1),each=n_groups/2)
  
  main_func = generate_smooth_func(x,n_funcs = 1,main_func_scale,main_func_amp)
  indiv_func = matrix(0, nrow=n_data, ncol=n_groups)
  for(j in 1:n_groups){
    indiv_func[,j] = generate_smooth_func(x,1,indiv_func_scale[j],
                                          indiv_func_amp)
  }
  full_func = indiv_func
  
  #adding individual and global functions together
  for(j in 1:n_groups){
    full_func[,j] = full_func[,j] + main_func[,1]
  }
  colnames(full_func) = paste("G",1:n_groups,sep="")
  
  full_data = full_func %>% 
    as.data.frame(.) %>%
    mutate(x = x, global_func = main_func[,1],
           indiv = 1:n_data)%>%
    gather(group,func_val, -x, -global_func, -indiv) %>%
    mutate(y= rnorm(n(),func_val,sqrt(noise)),
           indiv = paste(group, indiv, sep="_"))
  
  
  model_1 = bam(y~s(x,k=15),data=full_data,method="fREML", select=T)
  model_2a = update(model_1,formula. = y~s(x,k=15)+s(x,group,bs="fs",m=1,k=15))
  model_2b = update(model_1,formula. = y~ti(x,k=15)+ti(x,group,bs=c("tp","re"),
                                                       k=c(15,n_groups))+
                      ti(group,bs="re",k=n_groups))
  model_3 = update(model_1,formula. = y~s(x,k=15)+s(x,by=group,m=1,k=15)+group)
  model_4a = update(model_1,formula. = y~s(x,group,bs="fs",k=15))
  model_4b = update(model_1,formula. = y~te(x,group,bs=c("tp","re"),k=c(15,n_groups))+
                      te(group,bs="re",k=n_groups))
  model_5 = update(model_1,formula. = y~s(x,by=group,k=15)+group)
  
  aic_data$model_1[i] = AIC(model_1)
  aic_data$model_2a[i] = AIC(model_2a)
  aic_data$model_2b[i] = AIC(model_2b)
  aic_data$model_3[i] = AIC(model_3)
  aic_data$model_4a[i] = AIC(model_4a)
  aic_data$model_4b[i] = AIC(model_4b)
  aic_data$model_5[i] = AIC(model_5)
  print(i)
}
#generating main function and group-level functions ####
#Functions are genated using a Gaussian process with 

