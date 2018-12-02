library(mgcv)
library(MASS)
library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)


generate_smooth_func = function(x,n_funcs=1,length_scale=0.25,
                                amp = 1){
  # generates data from a Gaussian process with a squared exponential kernel
  # x is the input data, n_funcs is the number of functions to sample
  # length_scale determines how rapidly the function oscillates, and
  # amp determines the function's overall variance (over ranges of x >> length_scale,
  # var(func) will go to amp.).
  
  stopifnot(is.numeric(x))
  stopifnot(is.vector(x))
  x_dist = as.matrix(dist(x))
  x_cov = amp*exp(-x_dist^2/(2*length_scale^2))
  funcs = rmvn(n_funcs,rep(0, length=length(x)), x_cov)
  if(is.null(dim(funcs))){
    funcs = matrix(funcs, nrow=1)
  }
  return(t(funcs))
}


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

global_trend = c("Shared (Global) trend", "No shared trend" )
global_trend = factor(global_trend, levels =global_trend)
indiv_trend = c("No group-level trends",
                "Group-level trends\nsimilar smoothness\n(Shared penalty)",
                "Group-level trends\ndifferent smoothness\n(Individual penalties)")

models= as.data.frame(expand.grid(global_trend=global_trend, indiv_trend=indiv_trend))
models$model_num = c("model G",NA,"model GS","model S","model GI","model I")

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
  geom_label(data=models,aes(label=model_num),x=0.5,y=2)+
  theme_bw()+
  theme(legend.position="none",strip.background = element_blank(),
        text=element_text(size=14),
        panel.grid = element_blank())
  


print(multismooth_plot)

ggsave("figures/alternate_models.png", multismooth_plot,
       width=6, height= 8, units = "in",dpi = 300)