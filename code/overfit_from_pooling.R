library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(1)

calc_2nd_deriv = function(x,y){
  deriv_val = (lag(y) + lead(y) - 2*y)/(x-lag(x))^2
  deriv_val
}


freq_vals = c(1/4,1/2,1,2,4)
dat = crossing(x = seq(0,2*pi,length=150),freq = freq_vals)%>%
  mutate(y = sin(freq*x) +rnorm(n(), 0, 0.2),
         grp = paste("frequency = ",freq,sep= ""), 
         grp = factor(grp,  levels = paste("frequency = ",freq_vals,sep= "")))

mod1 = bam(y~s(x,k=30,grp, bs="fs"), data=dat)
mod2 = bam(y~s(x,k=30,by=grp)+s(grp,bs="re"), data=dat)




overfit_predict_data = crossing(x = seq(0,2*pi,length=500), freq = freq_vals)%>%
  mutate(grp = paste("frequency = ",freq,sep= ""), 
         grp = factor(grp,  levels = paste("frequency = ",freq_vals,sep= "")),
         y = sin(freq*x))%>%
  mutate(fit1 = as.numeric(predict(mod1,newdata = .,type = "response")),
         fit2 = as.numeric(predict(mod2,newdata = .,type="response")))

overfit_predict_data_long = overfit_predict_data %>%
  gather(model, value, y, fit1, fit2)%>%
  mutate(model = recode(model, y = "true value",fit1 = "model 4 fit",
                        fit2 = "model 5 fit"),
         model = factor(model, levels=  c("true value","model 4 fit", 
                                          "model 5 fit")))

deriv_est_data = overfit_predict_data%>%
  group_by(grp)%>%
  arrange(grp, x)%>%
  mutate(fit1_deriv = calc_2nd_deriv(x,fit1),
         fit2_deriv = calc_2nd_deriv(x,fit2))%>%
  summarize(freq = freq[1], fit1_int = sum(fit1_deriv^2*(x-lag(x)),na.rm = T),
            fit2_int = sum(fit2_deriv^2*(x-lag(x)),na.rm = T))%>%
  ungroup()%>%
  mutate(sqr_2nd_deriv = -freq^3*(sin(4*pi*freq)-4*pi*freq)/4)%>%
  gather(key=model,value = obs_sqr_deriv,fit1_int,fit2_int)%>%
  mutate(model = factor(ifelse(model=="fit1_int", "model 4 fit",
                               "model 5 fit"),
                        levels = c("model 4 fit", 
                                   "model 5 fit")))


deriv_plot =  ggplot(data=deriv_est_data, aes(x=sqr_2nd_deriv, y= obs_sqr_deriv,color= model))+
  geom_point()+
  scale_y_log10("integral of squared second\nderivative for fitted curves")+
  scale_x_log10("integral of squared second\nderivative for true curve")+
  scale_color_brewer(name=NULL,palette= "Set1")+
  geom_abline(color="black")+
  theme_bw()+
  theme(legend.position = "bottom")

fit_colors = c("black",RColorBrewer::brewer.pal(3, "Set1")[1:2])

overfit_vis_plot = ggplot(data=overfit_predict_data_long,aes(x=x,y= value,color=model))+
  geom_line()+
  scale_color_manual(values=fit_colors)+
  facet_grid(.~grp)+
  theme_bw()+
  theme(legend.position = "bottom",panel.grid = element_blank())
