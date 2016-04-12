library(dplyr)
library(mgcv)
library(ggplot2)
library(cowplot)
n= 100
dat = data_frame(x =seq(-3,3,length=n),
                 y= rnorm(n, cos(x*pi),0.2),
                 in_sample = ifelse(abs(x)<1,1,0),
                 extrapolate_val = factor(sign(x)*as.numeric(abs(x)>1)))

dat$y = ifelse(dat$in_sample, dat$y,0)

model_m1 = gam(y~s(x,m=1,k=20), weights = in_sample,data=dat)
model_m2 = gam(y~s(x,m=2,k=20), weights = in_sample,data=dat)

predict_m1 = predict(model_m1,dat,se=T)
predict_m2 = predict(model_m2,dat,se=T)

predict_data = rbind(dat, dat)
predict_data$model = rep(c("tprs m=1","tprs m=2"),each=n)
predict_data$fit = as.numeric(c(predict_m1$fit,predict_m2$fit))
predict_data$se.fit = as.numeric(c(predict_m1$se.fit,predict_m2$se.fit))

in_sample_plot = ggplot(aes(x=x, y=fit), data=filter(predict_data, in_sample==1))+
  geom_line()+
  geom_rug(sides="b")+
  geom_point(aes(y=y))+
  geom_ribbon(aes(ymin =fit-2*se.fit, ymax=fit+2*se.fit),
              alpha=0.25)+
  theme_bw()+
  facet_grid(model~.)

out_sample_plot = ggplot(aes(x=x, y=fit,color=extrapolate_val
                             ,fill =extrapolate_val), data=predict_data)+
  geom_line(group=1)+
  geom_rug(sides="b")+
  scale_color_manual(values = c("darkgrey","black","darkgrey"))+
  scale_fill_manual(values = c("darkgrey","black","darkgrey"))+
  geom_ribbon(aes(ymin =fit-2*se.fit, ymax=fit+2*se.fit),
              alpha=0.33,color=NA)+
  facet_grid(model~.,scales="free_y")+
  theme_bw()+
  theme(legend.position ="none")


print(in_sample_plot)
print(out_sample_plot)


ggsave("figures/extrapolating_with_dummypoints.png",
       plot_grid(in_sample_plot,out_sample_plot,align = "h",nrow = 1),
       width = 12,height=6,units="in",dpi=300)


