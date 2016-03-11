library(mgcv)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
make_null_space_penalty = function(x){
  stopifnot(is.matrix(x))
  eigen_trans = eigen(x)
  vals = eigen_trans$values
  vecs = eigen_trans$vectors
  null_vec = ifelse(vals<1e-5,1,0)
  null_pen_matrix = vecs%*%diag(null_vec)%*%solve(vecs)
  null_pen_matrix = null_pen_matrix/abs(max(null_pen_matrix))
  return(null_pen_matrix)
  
}

k = 6
plotting_data = data.frame(x = seq(0,1,length=100))
plotting_data$y = rnorm(100, sin(2*pi*plotting_data$x),0.2)
cr_basis = smoothCon(s(x,bs="cr",k=k+1), data=plotting_data,
                     knots=NULL,absorb.cons=TRUE)[[1]]
tp_basis = smoothCon(s(x,bs="tp",k=k+1), data=plotting_data,
                     knots=NULL,absorb.cons=TRUE)[[1]]

#### Extract basis functions #### 
cr_basis_funcr = as.data.frame(cr_basis$X)
names(cr_basis_funcr) = paste("F",1:k, sep="")
cr_basis_funcr$x = plotting_data$x
cr_basis_funcr$model = "Cubic spline"

tp_basis_funcr = as.data.frame(tp_basis$X)
names(tp_basis_funcr) = paste("F",1:k, sep="")
tp_basis_funcr$x = plotting_data$x
tp_basis_funcr$model = "Thin plate spline"

spline_basis_funcr = rbind(cr_basis_funcr,tp_basis_funcr)
spline_basis_funcr = gather(spline_basis_funcr,func,value, -x,-model )



##### Extract penalty matrices ####
cr_basis_P = as.data.frame(cr_basis$S[[1]])
cr_basis_P = cr_basis_P/max(cr_basis_P)
names(cr_basis_P) = paste("F",1:k, sep="")
cr_basis_P$basis_y = factor(paste("F",1:k, sep=""), 
                            levels= rev(paste("F",1:k, sep="")))
cr_basis_P$model = "Cubic spline"

tp_basis_P = as.data.frame(tp_basis$S[[1]])
tp_basis_P = tp_basis_P/max(tp_basis_P)
names(tp_basis_P) = paste("F",1:k, sep="")
tp_basis_P$basis_y = factor(paste("F",1:k, sep=""), 
                            levels= rev(paste("F",1:k, sep="")))
tp_basis_P$model = "Thin plate spline"

spline_basis_penalties = rbind(cr_basis_P ,tp_basis_P )
spline_basis_penalties = gather(spline_basis_penalties,basis_x,value,
                                -basis_y,-model )

spline_basis_penalties$penalty_type = "Smoothness penalty"

##### Create null penalty matrices ####
cr_basis_NP = as.data.frame(make_null_space_penalty(cr_basis$S[[1]]))
names(cr_basis_NP) = paste("F",1:k, sep="")
cr_basis_NP$basis_y = factor(paste("F",1:k, sep=""), 
                            levels= rev(paste("F",1:k, sep="")))
cr_basis_NP$model = "Cubic spline"

tp_basis_NP = as.data.frame(make_null_space_penalty(tp_basis$S[[1]]))
names(tp_basis_NP) = paste("F",1:k, sep="")
tp_basis_NP$basis_y = factor(paste("F",1:k, sep=""), 
                            levels= rev(paste("F",1:k, sep="")))
tp_basis_NP$model = "Thin plate spline"

spline_null_penalties = rbind(cr_basis_NP ,tp_basis_NP )
spline_null_penalties = gather(spline_null_penalties,basis_x,value,
                                -basis_y,-model )


spline_null_penalties$penalty_type = "Null space penalty"

spline_penalties = rbind(spline_basis_penalties,spline_null_penalties)
spline_penalties$penalty = factor(spline_penalties$penalty, 
                                  levels = c("Smoothness penalty", "Null space penalty"))


#### Creating plots for smoothness and null-space penalties
basis_func_plot = ggplot(aes(x=x,y=value),data=spline_basis_funcr)+
  geom_line()+
  scale_x_continuous(breaks=seq(0,1,length=5),
                     labels=c("0","0.25","0.5","0.75","1"))+
  facet_grid(model~func,switch="y",scales = "free_y")+
  theme_bw()+
  theme(strip.background = element_blank())

basis_penalty_plot = ggplot(aes(x=basis_x,y=basis_y,fill=value),
                            data=spline_penalties)+
  geom_tile(color="black")+
  facet_grid(model~penalty)+
  theme_bw()+
  scale_fill_gradient2("penalty",high = "#b2182b",low="#2166ac",midpoint = 0 )+
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
                      rel_widths = c(1.5,1),labels= c("A","B"))

ggsave("figures/fig.1 - example spline basis and penalties.png", full_plot, 
       width=14, height= 4, units = "in",dpi = 300)

