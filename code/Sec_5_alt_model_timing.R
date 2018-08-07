library(tidyr)
library(dplyr)
library(mgcv)
library(ggplot2)
library(gamm4)

set.seed = 1 #ensures that each new model parameter set is an extension of the old one

n_x = 20
x = seq(-2,2, length=n_x)
n_steps = 7




fit_timing_data = data_frame(n_groups = 2^(1:n_steps),
                             gam=rep(0,length=n_steps),
                             `bam (discrete = FALSE)`= 0,
                             `bam (discrete = TRUE)` = 0,
                             gamm = 0, gamm4 = 0)

fac_all = paste("g", 1:max(fit_timing_data$n_groups),sep = "")


model_coefs_all = data_frame(fac=fac_all)%>%
  mutate(int = rnorm(n(), 0,0.1),
         x2  = rnorm(n(),0,0.2),
         logit_slope = rnorm(n(),0, 0.2))

for(i in 1:n_steps){

  n_g =  fit_timing_data$n_groups[i]
  
  fac_current = fac_all[1:n_g]
  fac_current = factor(fac_current, levels=  unique(fac_current))

  model_coefs = model_coefs_all%>%
    filter(fac %in% fac_current)%>%
    mutate(fac = factor(fac, levels= unique(fac)))
  
  set.seed = 1 #ensures that each new data set is an extension of the old one
  
  model_data = crossing(fac=fac_current, x=x)%>%
    left_join(model_coefs)%>%
    mutate(base_func  = dnorm(x)*10,
           indiv_func = int + x^2*x2 + 2*(exp(x*logit_slope)/(1+exp(x*logit_slope))-0.5),
           y = base_func + indiv_func + rnorm(n()))
  
  fit_timing_data$gam[i] = system.time(gam(y~s(x,k=10, bs="cp") + s(x,fac, k=10, bs="fs", xt=list(bs="cp"), m=2),
                             data= model_data, method="REML"))[3]
  
  fit_timing_data$`bam (discrete = FALSE)`[i] = system.time(bam(y~s(x,k=10, bs="cp") + s(x,fac, k=10, bs="fs", xt=list(bs="cp"),m=2),
                             data= model_data, discrete = FALSE))[3]
  
  fit_timing_data$`bam (discrete = TRUE)`[i] = system.time(bam(y~s(x,k=10, bs="cp") + s(x,fac, k=10, bs="fs", xt=list(bs="cp"),m=2),
                                                                data= model_data,discrete=TRUE))[3]
  
  
  fit_timing_data$gamm[i] = system.time(gamm(y~s(x,k=10, bs="cp") + s(x,fac, k=10, bs="fs", xt=list(bs="cp"),m=2),
                               data= model_data))[3]
  
  
  fit_timing_data$gamm4[i] = system.time(gamm4(y~s(x,k=10, bs="cp") + s(x,fac, k=10, bs="fs", xt=list(bs="cp"),m=2),
                                 data= model_data))[3]
}


fit_timing_long = fit_timing_data %>% 
  gather(model, timing, gam,`bam (discrete = FALSE)`, 
         `bam (discrete = TRUE)`, gamm, gamm4)%>%
  mutate(model =factor(model, levels = c("gam",
                                         "bam (discrete = FALSE)",
                                         "bam (discrete = TRUE)",
                                         "gamm", 
                                         "gamm4")))


timing_plot = ggplot(aes(n_groups, timing, color=model, linetype= model), 
                     data=fit_timing_long)+
  geom_line()+
  geom_point()+
  scale_color_manual(values = c("black", "#1b9e77","#1b9e77", "#d95f02", "#7570b3"))+
  scale_linetype_manual(values =c(1,1,2,1,1))+
  scale_y_log10("run time (seconds)", breaks = c(0.1,1,10,100), labels = c("0.1", "1","10", "100"))+
  scale_x_log10("number of groups", breaks = c(2,8,32,128))+
  
  theme_bw()+
  guides(color = guide_legend(nrow = 2,byrow = TRUE))+
  theme(panel.grid.minor  = element_blank(),panel.grid.major.x = element_blank(),
        legend.position = "bottom")

ggsave("figures/alternate_model_timing_plot.png",timing_plot,height=4,width=6,units="in", dpi=400)
