library(mgcv)
library(MASS)
library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

source("code/functions.R")

#starting parameters  ####  
n_groups = 12 # number of groups
n_replicates = 20
n_holdout_data = 25

total_amp = 1
noise  = total_amp/3 #variance of random noise around main function

# The length-scale parameter for the main and individual-level functions. Points
# with x-values that are more distant from each other than the length scale
# will be only weakly correlated
main_func_scale= 0.1
indiv_func_scale_base = 0.1

sim_data = rbind(
  expand.grid(trial = "varying main amp",rep = 1:n_replicates,
              n_data = c(25,50,100,200),
              main_func_amp = c(0,0.25,0.5,0.75,1.0),
              indiv_scale_diff = c(1,2)),
  expand.grid(trial = "varying group scale",rep = 1:n_replicates,
              n_data = c(25,50,100,200),
              main_func_amp = c(0,0.5),
              indiv_scale_diff = c(1,2,4,8)
              ))
sim_data = data.frame(sim_data)
sim_data$seed = 1:nrow(sim_data)

aic_data = mutate(sim_data,
                  model_1 = 0,model_2a = 0,model_2b=0,model_3 =0,
                  model_4a= 0,model_4b = 0,model_5 =0)

rmse_data = aic_data

for(i in 1:nrow(sim_data)){
  set.seed(sim_data$seed[i])
  n_data = sim_data$n_data[i]
  x = seq(0,1, length=n_data+n_holdout_data)
  indiv_scale_diff = sim_data$indiv_scale_diff[i]
  # variability of the main and individual level functions. Equal to the variance
  # of individual points drawn from the function at large distances from one another.
  # set so that the variance of the main and individual level-functions will sum
  # to a fixed value
  main_func_amp = sim_data$main_func_amp[i]
  indiv_func_amp = total_amp-main_func_amp #var

  indiv_func_scale = indiv_func_scale_base*sqrt(indiv_scale_diff)^rep(c(-1,1),
                                                                      each=n_groups/2)
  
  main_func = generate_smooth_func(x,n_funcs = 1,main_func_scale,main_func_amp)
  indiv_func = matrix(0, nrow=n_data +n_holdout_data, ncol=n_groups)
  for(j in 1:n_groups){
    indiv_func[,j] = generate_smooth_func(x,1,indiv_func_scale[j],
                                          indiv_func_amp)
  }
  
  
  #adding individual and global functions together
  full_func = indiv_func
  for(j in 1:n_groups){
    full_func[,j] = full_func[,j] + main_func[,1]
  }
  colnames(full_func) = paste("G",1:n_groups,sep="")
  
  full_data = full_func %>% 
    as.data.frame(.) %>%
    mutate(x = x, global_func = main_func[,1],
           indiv = 1:(n_data+n_holdout_data))%>%
    gather(group,func_val, -x, -global_func, -indiv) %>%
    mutate(y= rnorm(n(),func_val,sqrt(noise)),
           indiv = paste(group, indiv, sep="_"))
  fit_x = sample(x,n_data)
  fit_data = full_data[full_data$x%in%fit_x,] 
  test_data = full_data[!full_data$x%in%fit_x,] 
  
  model_1 = bam(y~s(x,k=15),data=fit_data,method="fREML", select=T)
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
  
  
  rmse_data$model_1[i] = calc_rmse(test_data$y,predict(model_1,test_data))
  rmse_data$model_2a[i] = calc_rmse(test_data$y,predict(model_2a,test_data))
  rmse_data$model_2b[i] = calc_rmse(test_data$y,predict(model_2b,test_data))
  rmse_data$model_3[i] = calc_rmse(test_data$y,predict(model_3,test_data))
  rmse_data$model_4a[i] = calc_rmse(test_data$y,predict(model_4a,test_data))
  rmse_data$model_4b[i] = calc_rmse(test_data$y,predict(model_4b,test_data))
  rmse_data$model_5[i] = calc_rmse(test_data$y,predict(model_5,test_data))
  print(i)
}


#scale aic to compare all models against the best fit one
model_list = c("model_1","model_2a","model_2b","model_3","model_4a","model_4b","model_5")
model_list_simplified = c("model 1", "model 2", "model 3","model 4", "model 5")
aic_data[,model_list]= t(apply(aic_data[,model_list],MARGIN = 1,function(x)x-min(x)))

aic_data$best_model = apply(aic_data[,model_list],1, 
                            function(x) model_list[which(x==min(x))])

aic_data$best_model_simplified = str_replace(aic_data$best_model, 
                                             "(model)(_)([0-9])([abc]*)","\\1 \\3")

                                             
rmse_data$best_model = apply(rmse_data[,model_list],1, 
                            function(x) model_list[which(x==min(x))][1])

rmse_data$best_model_simplified = str_replace(rmse_data$best_model, 
                                             "(model)(_)([0-9])([abc]*)","\\1 \\3")

aic_summary_data = aic_data %>%
  group_by(trial,main_func_amp,indiv_scale_diff, n_data,best_model_simplified,
           best_model_simplified)%>%
  summarize(n_wins = n())%>%
  group_by(trial)%>%
  mutate(main_amp_lab = paste("var. expl. by\nglobal func: ",main_func_amp,sep=""),
         indiv_scale_lab= paste("ratio of scales\nbetween groups: ",
                                indiv_scale_diff,sep=""))

rmse_summary_data = rmse_data %>%
  group_by(trial,main_func_amp,indiv_scale_diff, n_data,best_model_simplified)%>%
  summarize(n_wins = n())%>%
  group_by(trial)%>%
  mutate(main_amp_lab = paste("var. expl. by\nglobal func: ",main_func_amp,sep=""),
         indiv_scale_lab= paste("ratio of scales\nbetween groups: ",
                                indiv_scale_diff,sep=""))
  

fit_plot_data = gather(aic_data, key = model, value = AIC,model_1:model_5)
fit_plot_data = left_join(fit_plot_data, 
                          gather(rmse_data, 
                                 key = model, value = RMSE,model_1:model_5)) %>%
  group_by(seed)%>%
  mutate(RMSE_scaled = exp(log(RMSE)-min(log(RMSE))),
         AIC_rank = rank(AIC),
         RMSE_rank = rank(RMSE))%>%
  mutate(main_amp_lab = paste("var. expl. by\nglobal func: ",main_func_amp,sep=""),
         indiv_scale_lab= paste("ratio of scales\nbetween groups: ",
                                indiv_scale_diff,sep=""))



color_list = c("model 1" = "white","model 2" = "#a6cee3","model 3" = "#b2df8a",
                "model 4"= "#1f78b4","model 5" ="#33a02c")

model_labels =  c("Model 1:\nGlobal trend","Model 2:\nGlobal trend\nshared smooths" ,
                  "Model 3:\nGlobal trend\nseperate smooths",
                  "Model 4:\nNo trend\nshared smooths","Model 5:\nno trend\nseperate smooths")


plot_base = list(geom_bar(stat="identity",color="black"),
                 scale_x_log10("number of observations", breaks=c(25,50,100,200)),
                 scale_fill_manual("", values=color_list,
                                   labels = model_labels),
                 scale_y_continuous("frequency model is\nchosen via AIC",
                                    limits=c(0,1),expand = c(0,0)),
                 theme_bw(),
                 theme(panel.margin.y = unit(1, "lines"),
                       legend.position="bottom")
                 )

aic_compare_main_plot = ggplot(aes(x=n_data, y=n_wins/n_replicates,
                                   fill=best_model_simplified),
                               data=filter(aic_summary_data,
                                           trial=="varying main amp")) +
  facet_grid(indiv_scale_lab~main_amp_lab)+
  plot_base


aic_compare_group_plot = ggplot(aes(x=n_data, y=n_wins/n_replicates,
                                    fill=best_model_simplified),
                               data=filter(aic_summary_data,
                                           trial=="varying group scale")) +
  facet_grid(main_amp_lab~indiv_scale_lab)+
  plot_base


rmse_compare_main_plot = ggplot(aes(x=n_data, y=n_wins/n_replicates,
                                   fill=best_model_simplified),
                               data=filter(rmse_summary_data,
                                           trial=="varying main amp")) +
  facet_grid(indiv_scale_lab~main_amp_lab)+
  plot_base+
  scale_y_continuous("frequency model is\nchosen via RMSE",
                     limits=c(0,1),expand = c(0,0))


rmse_compare_group_plot = ggplot(aes(x=n_data, y=n_wins/n_replicates,
                                    fill=best_model_simplified),
                                data=filter(rmse_summary_data,
                                            trial=="varying group scale")) +
  facet_grid(main_amp_lab~indiv_scale_lab)+
  plot_base+
  labs(y= "frequency model is chosen via RMSE")+
  scale_y_continuous("frequency model is\nchosen via RMSE",
                     limits=c(0,1),expand = c(0,0))



ggsave("figures/comparing AIC fits plot 1.png",aic_compare_main_plot,
       width=8, height=4, units="in",dpi = 300)


ggsave("figures/comparing AIC fits plot 2.png",aic_compare_group_plot,
       width=8, height=4, units="in",dpi = 300)


ggsave("figures/comparing RMSE fits plot 1.png",rmse_compare_main_plot,
       width=8, height=4, units="in",dpi = 300)


ggsave("figures/comparing RMSE fits plot 2.png",rmse_compare_group_plot,
       width=8, height=4, units="in",dpi = 300)
