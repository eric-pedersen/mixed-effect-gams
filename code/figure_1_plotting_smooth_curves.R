library(mgcv)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
k = 6

plotting_data = data.frame(x = seq(0,1,length=100))
cs_basis = smoothCon(s(x,bs="cs",k=k+1), data=plotting_data,
                     knots=NULL,absorb.cons=TRUE)[[1]]
tp_basis = smoothCon(s(x,bs="tp",k=k+1), data=plotting_data,
                     knots=NULL,absorb.cons=TRUE)[[1]]

#### Extract basis functions #### 
cs_basis_funcs = as.data.frame(cs_basis$X)
names(cs_basis_funcs) = paste("F",1:k, sep="")
cs_basis_funcs$x = plotting_data$x
cs_basis_funcs$model = "Cubic spline"

tp_basis_funcs = as.data.frame(tp_basis$X)
names(tp_basis_funcs) = paste("F",1:k, sep="")
tp_basis_funcs$x = plotting_data$x
tp_basis_funcs$model = "Thin plate spline"

spline_basis_funcs = rbind(cs_basis_funcs,tp_basis_funcs)
spline_basis_funcs = gather(spline_basis_funcs,func,value, -x,-model )



##### Extract penalty matrices ####
cs_basis_P = as.data.frame(cs_basis$S)
cs_basis_P = cs_basis_P/max(cs_basis_P)
names(cs_basis_P) = paste("F",1:k, sep="")
cs_basis_P$basis_y = factor(paste("F",1:k, sep=""), 
                            levels= rev(paste("F",1:k, sep="")))
cs_basis_P$model = "Cubic spline"

tp_basis_P = as.data.frame(tp_basis$S)
tp_basis_P = tp_basis_P/max(tp_basis_P)
names(tp_basis_P) = paste("F",1:k, sep="")
tp_basis_P$basis_y = factor(paste("F",1:k, sep=""), 
                            levels= rev(paste("F",1:k, sep="")))
tp_basis_P$model = "Thin plate spline"

spline_basis_penalties = rbind(cs_basis_P ,tp_basis_P )
spline_basis_penalties = gather(spline_basis_penalties,basis_x,value,
                                -basis_y,-model )




basis_func_plot = ggplot(aes(x=x,y=value),data=spline_basis_funcs)+
  geom_line()+
  scale_x_continuous(breaks=seq(0,1,length=5),
                     labels=c("0","0.25","0.5","0.75","1"))+
  facet_grid(model~func,switch="y",scales = "free_y")+
  theme_bw()+
  theme(strip.background = element_blank())

basis_penalty_plot = ggplot(aes(x=basis_x,y=basis_y,fill=value),
                            data=spline_basis_penalties)+
  geom_tile(color="black")+
  facet_grid(model~.)+
  theme_bw()+
  scale_fill_gradient2("penalty")+
  theme(strip.background = element_blank())+
  labs(x="", y="")+
  coord_fixed()+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme(axis.ticks=element_blank(),
        panel.grid=element_blank())

print(basis_func_plot)
print(basis_penalty_plot)
  
full_plot = plot_grid(basis_func_plot,basis_penalty_plot,ncol=2,
                      rel_widths = c(2,1),labels= c("A","B"))

ggsave("figures/fig.1 - example spline basis and penalties.png", full_plot, 
       width=12, height= 4, units = "in",dpi = 300)

