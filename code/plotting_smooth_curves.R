library(mgcv)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)

k = 6
plotting_data = data.frame(x = seq(0,1,length=100))

tp_basis = smoothCon(s(x,bs="tp",k=k), data=plotting_data,
                     knots=NULL,absorb.cons=FALSE)[[1]]

#### Extract basis functions #### 
tp_basis_funcr = as.data.frame(tp_basis$X)
names(tp_basis_funcr) = paste("F",1:k, sep="")
tp_basis_funcr$x = plotting_data$x
tp_basis_funcr$model = "Thin plate spline"

spline_basis_funcr = gather(tp_basis_funcr,func,value, -x,-model )



##### Extract penalty matrices ####

tp_basis_P = as.data.frame(tp_basis$S[[1]])
tp_basis_P = tp_basis_P/max(tp_basis_P)
names(tp_basis_P) = paste("F",1:k, sep="")
tp_basis_P$basis_y = factor(paste("F",1:k, sep=""), 
                            levels= rev(paste("F",1:k, sep="")))
tp_basis_P$model = "Thin plate spline"
spline_basis_penalties = gather(tp_basis_P ,basis_x,value,
                                -basis_y,-model )

spline_basis_penalties$penalty_type = "Smoothness penalty"


#### Creating plots for smoothness and null-space penalties
basis_func_plot = ggplot(aes(x=x,y=value),data=spline_basis_funcr)+
  geom_line()+
  scale_x_continuous(breaks=seq(0,1,length=3),
                     labels=c("0","0.5","1"))+
  facet_wrap(~func)+
  theme_bw()+
  theme(strip.background = element_blank(), panel.grid = element_blank())

basis_penalty_plot = ggplot(aes(x=basis_x,y=basis_y,fill=value),
                            data=spline_basis_penalties)+
  geom_tile(color="black")+
  theme_bw()+
  scale_fill_gradient2("penalty",high = "#b2182b",low="#2166ac",midpoint = 0 )+
  theme(strip.background = element_blank())+
  labs(x="", y="")+
  coord_fixed()+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme(axis.ticks=element_blank(),
        panel.grid=element_blank())

  
full_plot = plot_grid(basis_func_plot,basis_penalty_plot,ncol=2,
                      rel_widths = c(1,1),labels= c("A","B"))

ggsave("figures/example spline basis and penalties.png", full_plot, 
       width=8, height= 4, units = "in",dpi = 300)

