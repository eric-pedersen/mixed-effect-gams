library(mgcv)
library(MASS)
library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)

source("code/functions.R")

#### Global parameters
set.seed(7)
n_groups =  4
n_points = 200
x = seq(0,1,length=n_points)
groups= letters[1:n_groups]

total_amp = 1
main_func_amp = 0.5
main_func_scale= 0.3
sep_scales = main_func_scale*2^(seq(-2,2,length=n_groups))

global_trend = c("Shared trend", "No shared trend" )
global_trend = factor(global_trend, levels =global_trend)
indiv_trend = c("No group-level trends",
                "Group-level trends\nsimilar smoothness",
                "Group-level trends\ndifferent smoothness")

models= as.data.frame(expand.grid(global_trend=global_trend, indiv_trend=indiv_trend))
models$model_num = c("model 1",NA,"model 2","model 4","model 3","model 5")

global_func = generate_smooth_func(x,n_funcs = 1,main_func_scale,total_amp)

indiv_func_single_smooth = generate_smooth_func(x,n_funcs = n_groups,
                                                main_func_scale,total_amp)
indiv_func_multi_smooth = matrix(0, ncol=n_groups,nrow=n_points)
for(i in 1:n_groups){
  indiv_func_multi_smooth[,i] = generate_smooth_func(x,n_funcs = 1,sep_scales[i],
                                                     total_amp)
}
colnames(indiv_func_multi_smooth) = groups
colnames(indiv_func_single_smooth) = groups


plotting_data = as.data.frame(expand.grid(x=x,groups=groups, 
                                          global_trend=global_trend, indiv_trend= indiv_trend))
plotting_data$global_function = 0
plotting_data$indiv_function = 0
plotting_data$total_function = 0

for(i in global_trend){
  for(j in indiv_trend){
    for(k in groups){
      
      indices = with(plotting_data, which(groups==k&global_trend==i&indiv_trend==j))

      if(j==indiv_trend[1]){
        plotting_data[indices,"indiv_function"] = 0
      }else if(j==indiv_trend[2]){
        plotting_data[indices,"indiv_function"] = indiv_func_single_smooth[,k]
      }else{
        plotting_data[indices,"indiv_function"] = indiv_func_multi_smooth[,k]
      }
      if(i==global_trend[2]){
        plotting_data[indices,"global_function"]= 0
      }else{
        plotting_data[indices,"global_function"]= global_func
        if(j!="No group-level trends"){
          plotting_data[indices,"global_function"] = main_func_amp*plotting_data[indices,"global_function"] 
          plotting_data[indices,"indiv_function"] = (1-main_func_amp)*plotting_data[indices,"indiv_function"] 
        }
      }
    }
  }
  
}
plotting_data$total_function = plotting_data$indiv_function+plotting_data$global_function

multismooth_plot= ggplot(data=plotting_data, aes(x=x, y=total_function))+
  geom_line(colour="darkgrey",aes(group=groups),size=1)+
  facet_grid(indiv_trend~global_trend)+
  scale_color_brewer(palette = "Set1")+
  geom_line(size=1, aes(y=global_function),linetype=2)+
  scale_y_continuous("f(x)")+
  scale_x_continuous(breaks=c(0,0.5, 1), labels=c("0","0.5","1"))+
  theme_bw()+
  geom_label(data=models,aes(label=model_num),x=0.5,y=2)+
  theme(legend.position="none",strip.background = element_blank(),
        text=element_text(size=14))
  


print(multismooth_plot)

ggsave("figures/fig.2 - alternate models of functional variability.png", multismooth_plot,
       width=6, height= 8, units = "in",dpi = 300)