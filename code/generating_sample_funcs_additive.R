library(mgcv)
library(MASS)
library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)

source("code/functions.R")
  

#starting parameters  ####  
n_data = 60 #number of data points per group
n_groups = 25 # number of groups
holdout_frac= 0.33 #fraction of data held out from training models

total_amp = 1
noise  = total_amp/2 #variance of random noise around main function

# variability of the main and individual level functions. Equal to the variance
# of individual points drawn from the function at large distances from one another.
# set so that the variance of the main and individual level-functions will sum
# to a fixed value
main_func_amp = 0.6 
indiv_func_amp = total_amp-main_func_amp #var

# The length-scale parameter for the main and individual-level functions. Points
# with x-values that are more distant from each other than the length scale
# will be only weakly correlated
main_func_scale= 0.2
indiv_func_scale = 0.1


x = seq(0,1,length=n_data)

#generating main function and group-level functions ####
#Functions are genated using a Gaussian process with 
main_func = generate_smooth_func(x,n_funcs = 1,main_func_scale,main_func_amp)
indiv_func = generate_smooth_func(x,n_funcs = n_groups,indiv_func_scale,
                                  indiv_func_amp)


#Plots the global and individual-level functions
matplot(x, indiv_func,type='l',ylim = range(cbind(main_func,indiv_func)))
points(x, main_func[,1],type="l",lwd= 2)

full_func = indiv_func

#adding individual and global functions together
for(i in 1:n_groups){
  full_func[,i] = full_func[,i] + main_func[,1]
}
colnames(full_func) = paste("G",1:n_groups,sep="")


#splitting data into full data, and training sets ####
full_data = full_func %>% 
  as.data.frame(.) %>%
  mutate(x = x, global_func = main_func[,1],
         indiv = 1:n_data)%>%
  gather(group,func_val, -x, -global_func, -indiv) %>%
  mutate(y= rnorm(n(),func_val,sqrt(noise)),
         indiv = paste(group, indiv, sep="_"))
  
# first training data: random removal of a fraction of the data
test_random_data = full_data %>%
  group_by(group)%>%
  filter(indiv %in% sample(indiv,n_data*(1-holdout_frac)))
  
#second training data: removal of a continuous block of each function 
test_block_data = full_data %>%
  group_by(group)%>%
  mutate(hold_start = runif(1, 0,1))%>%
  filter(!between((x-hold_start)%%1,0, holdout_frac))

full_data = full_data %>%
  mutate(in_random_holdout =!indiv%in%test_random_data$indiv,
         in_block_holdout =!indiv%in%test_block_data$indiv)



#data frame specifying model labels for later
model_labels = data.frame(original = c("global","fixef",
                                       "fs_basic", "fs_full"),
                          new = c("global smooth", "fixed effect smooths",
                                  "random effects smooths",
                                  "global smooth + random interactions"))
model_labels$new = factor(model_labels$new, levels=model_labels$new)


#Fitting randomly chosen training data ####
#model fit with only group-level intercepts and a global smoother
global_random_fit = gam(y~s(x,bs="tp", k=15)+group, 
                         data=test_random_data)

#model fit with a unique, unconnected smoother and mean for each group
fixef_random_fit  = gam(y~s(x,bs="tp",by=group, k=15)+group, 
                         data=test_random_data)

#model fit with random effect smoother between groups, integrating interaction and main terms into a single smoother
fs_basic_random_fit  = gam(y~s(x,group, bs=c("fs"),k=15), 
                             data=test_random_data)

#model fit with main smoothers and a seperate interaction random effect term
fs_full_random_fit  = gam(y~s(x, bs="tp",k=15) +
                               s(x,group,bs=c("fs"),m=1,
                                  k=15),
                        data=test_random_data,select=T)

#Testing models on randomly chosen testing data####
fit_random_data = full_data 
fit_random_data$global = as.vector(predict(global_random_fit,full_data))
fit_random_data$fixef = as.vector(predict(fixef_random_fit,full_data))
fit_random_data$fs_basic = as.vector(predict(fs_basic_random_fit,full_data))
fit_random_data$fs_full = as.vector(predict(fs_full_random_fit,full_data))
fit_random_data= fit_random_data %>%
  dplyr::select(global:fs_full, group,in_random_holdout,x,y)%>%
  gather(model, fit,global, fixef,fs_basic,fs_full )%>%
  mutate(model = model_labels$new[match(model, model_labels$original)])


#Fitting block-chosen training data ####
global_block_fit = gam(y~s(x,bs="tp", k=15)+group, 
                        data=test_block_data)
fixef_block_fit  = gam(y~s(x,bs="tp",by=group, k=15)+group, 
                        data=test_block_data)
fs_basic_block_fit  = gam(y~s(x,group,bs=c("fs"), 
                                   k=15),
                              data=test_block_data)
fs_full_block_fit  = gam(y~s(x,group,bs="fs",
                                  k=15,m=1)+s(x, bs="tp",k=15),
                            select=T,data=test_block_data)


#Testing models on block-chosen testing data####
fit_block_data = full_data 
fit_block_data$global = as.vector(predict(global_block_fit,full_data))
fit_block_data$fixef = as.vector(predict(fixef_block_fit,full_data))
fit_block_data$fs_basic = as.vector(predict(fs_basic_block_fit,full_data))
fit_block_data$fs_full = as.vector(predict(fs_full_block_fit,full_data))
fit_block_data= fit_block_data %>%
  dplyr::select(global:fs_full, group,in_block_holdout,x,y)%>%
  gather(model, fit,global, fixef,fs_basic,fs_full ) %>%
  mutate(model = model_labels$new[match(model, model_labels$original)])


# Calculating root-mean square errors for training, testing, and total data ####
fit_block_summary = fit_block_data %>%
  group_by(group,model) %>%
  summarise(rmse_holdout  = calc_rmse(fit[in_block_holdout],
                                    y[in_block_holdout]),
            rmse_insample = calc_rmse(fit[!in_block_holdout],
                                      y[!in_block_holdout]),
            rmse_total  = calc_rmse(fit,y))

fit_random_summary = fit_random_data %>%
  group_by(group,model) %>%
  summarise(rmse_holdout  = calc_rmse(fit[in_random_holdout],
                                      y[in_random_holdout]),
            rmse_insample = calc_rmse(fit[!in_random_holdout],
                                      y[!in_random_holdout]),
            rmse_total  = calc_rmse(fit,y))


print("random removal rmse")
print(fit_random_summary %>% group_by(model)%>%
        summarise(holdout =round(mean(rmse_holdout),1),
                  insample=round(mean(rmse_insample),1),
                  total = round(mean(rmse_total),1)))

print("block removal rmse")
print(fit_block_summary %>% group_by(model)%>%
        summarise(holdout =round(mean(rmse_holdout),1),
                  insample=round(mean(rmse_insample),1),
                  total = round(mean(rmse_total),1)))


#Model fit plots ####
random_fit_plot = ggplot(aes(x=x,y=y), data=full_data)+
  geom_point(color=ifelse(full_data$in_random_holdout,
                          "grey","black"))+
  facet_wrap(~group)+
  theme_bw()+
  geom_line(aes(y= fit,color=model), data=fit_random_data,
            size=1)+
  scale_color_brewer(palette="Set1")+ 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position="bottom")

block_fit_plot = ggplot(aes(x=x,y=y), data=full_data)+
  geom_point(color=ifelse(full_data$in_block_holdout,
                          "grey","black"))+
  geom_line(aes(y=func_val))+
  facet_wrap(~group)+
  theme_bw()+
  geom_line(aes(y= fit,color=model), data=fit_block_data,
            size=1)+
  scale_color_brewer(palette="Set1")+ 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position="bottom")



print(plot_grid(random_fit_plot, block_fit_plot,ncol = 2,align = "h"))
