library(MASS)
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


##### Creating a random draw from the function ####
set.seed(6) #ensure we can 
coef_sample = rmvn(n = 1, mu = rep(0, times = k),V = 4*ginv(tp_basis$S[[1]]))
coef_sample[(k-1):k] = c(-0.1,0.12) #randomly draw the null space terms from a normal distribution
random_basis = as.matrix(tp_basis_funcr[,1:6])%*%diag(coef_sample) 
random_basis = as.data.frame(random_basis)
names(random_basis) = names(tp_basis_funcr)[1:6]

tp_example_curve = tp_basis_funcr %>%
  select(-matches("F[1-9]")) %>% #remove the old basis functions
  bind_cols(random_basis) %>%#append the new basis functions multiplied by a random sample of coefficients
  mutate(Ftotal = rowSums(select(., matches("F[1-9]"))))%>%
  gather(key = `basis function`, value = value, matches("F[1-9]"))

  
bs_func_labels = tp_example_curve %>%
  group_by(`basis function`)%>%
  summarise(value = value[x==max(x)],
            x     = max(x))%>%
  mutate(coef = round(coef_sample,2),
         basis_label = paste("paste(",`basis function`,"%*%", coef, ")", sep=""))




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


basis_sample_plot = ggplot(data= tp_example_curve, aes(x, value))+
  geom_line(aes(y = Ftotal),size=2)+
  geom_line(aes(group = `basis function`,color= `basis function`),size=0.5)+ 
  geom_text(data= bs_func_labels,
            aes(label = basis_label, 
                color= `basis function`,
                y = value),
            parse=T, hjust = 0,nudge_x = 0.01)+
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,1.15))+
  scale_color_viridis_d()+
  guides(color = "none")


top_row_plot = plot_grid(basis_func_plot,basis_penalty_plot,
                         ncol=2,
                         rel_widths = c(1.,1),labels= c("a","b"))
  
full_plot = plot_grid(top_row_plot,basis_sample_plot,
                      nrow=2,labels= c("", "c"),axis = "t",rel_heights = c(1,0.9))

ggsave("figures/example spline basis and penalties.png", full_plot, 
       width=8, height= 7, units = "in",dpi = 400)

