####Setup ####
#all the packages needed for this tutorial are listed here
library(mgcv)
library(MASS)
library(stringr)
library(gamm4)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(viridis)
library(cowplot)
library(kableExtra)
library(docxtools)
library(knitr)
library(tibble)
library(dplyr)
library(gratia)
library(latex2exp)

#Set the default theme for ggplot objects to theme_bw()
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

#Check for the R version number, and if greater than 3.6, switch to using the
#old random number generator, to ensure replicability
#Thanks to Github user @bastistician for pointing this out!

if(getRversion()>= 3.6) RNGversion("3.5.0")





#### Code for part I: Introduction ####


#### Code for II: A review of Generalized Additive Models ####

# example of varying lambda

set.seed(12)

# generate some data
dat <- gamSim(1, n=100, dist="normal", scale=2)
dat$y <- dat$y - (dat$f1 + dat$f0 + dat$f3)
dat$x <- dat$x2
true <- data.frame(x = sort(dat$x),
                   y = dat$f2[order(dat$x)])

## REML fit
b <- gam(y~s(x, k=100), data=dat, method = "REML")

# lambda=0
b.0 <- gam(y~s(x, k=100), data=dat, sp=0)

# lambda=infinity
b.inf <- gam(y~s(x, k=100), data=dat, sp=1e10)

# merging predictions from the models together to make plotting easier
pdat <- with(dat, data.frame(x = seq(min(x), max(x), length = 200)))
p <- cbind(pdat, fit = predict(b, newdata = pdat))
p.0 <- cbind(pdat, fit = predict(b.0, newdata = pdat))
p.inf <- cbind(pdat, fit = predict(b.inf, newdata = pdat))
ylims <- range(p, p.0, p.inf)

lab.l <- labs(x = "x", y = "y")
dat.l <- geom_point(data = dat, aes(x = x, y = y), colour = "darkgrey")
true.l <- geom_line(data = true, aes(x = x, y = y), colour = "blue")
coord.l <- coord_cartesian(ylim = ylims)

#plotting models
p1 <- ggplot(p, aes(x = x, y = fit)) +
  dat.l + true.l +
  geom_line(colour = "darkred") + lab.l + coord.l

p2 <- ggplot(p.0, aes(x = x, y = fit)) +
  dat.l + true.l +
  geom_line(colour = "darkred") + lab.l + coord.l

p3 <- ggplot(p.inf, aes(x = x, y = fit)) +
  dat.l + true.l +
  geom_line(colour = "darkred") + lab.l + coord.l

plot_grid(p1, p2, p3, align = "hv", axis = "lrtb", ncol = 3, labels = "AUTO")


#Code for generating figure 3: examples of basis functions and splines ####
k = 6
plotting_data = data.frame(x = seq(0,1,length=100))

#This creates the basis functions for a thin-plate spline The absorb.cons=FALSE
#setting makes sure that smoothCon does not remove basis functions that have a
#non-zero sum (in this case, the intercept). Absorbing constraints would result
#in having less than k basis functions, which is why fitted terms in mgcv often
#have less than k maximum EDF.
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
#ensure we can reproduce this example
set.seed(6) 
coef_sample = rmvn(n = 1, mu = rep(0, times = k),V = 4*ginv(tp_basis$S[[1]]))

#randomly draw the null space terms from a normal distribution
coef_sample[(k-1):k] = c(-0.1,0.12)  
random_basis = as.matrix(tp_basis_funcr[,1:6])%*%diag(coef_sample) 
random_basis = as.data.frame(random_basis)
names(random_basis) = names(tp_basis_funcr)[1:6]

tp_example_curve = tp_basis_funcr %>%
  #remove the old basis functions
  select(-matches("F[1-9]")) %>% 
  #append the new basis functions multiplied by a random sample of coefficients
  bind_cols(random_basis) %>%
  mutate(Ftotal = rowSums(select(., matches("F[1-9]"))))%>%
  gather(key = `basis function`, value = value, matches("F[1-9]"))

  
bs_func_labels = tp_example_curve %>%
  group_by(`basis function`)%>%
  summarise(value = value[x==max(x)],
            x     = max(x))%>%
  mutate(coef = round(coef_sample,2),
         basis_label = paste("paste(",`basis function`,"%*%", coef, ")", 
                             sep=""))




#### Creating plots for smoothness and null-space penalties
#getting a color-blind palette for the 6 levels, avoiding the yellow and black 
#levels
basis_func_palette = colorblind_pal()(8)
basis_func_palette = basis_func_palette[-c(1,5)]

#Plotting the basis functions
basis_func_plot = ggplot(aes(x=x,y=value,color=func),data=spline_basis_funcr)+
  geom_line()+
  scale_x_continuous(breaks=seq(0,1,length=3),
                     labels=c("0","0.5","1"))+
  facet_wrap(~func)+
  scale_color_manual(values = basis_func_palette)+
  guides(color = "none")

basis_penalty_plot = ggplot(aes(x=basis_x,y=basis_y,fill=value),
                            data=spline_basis_penalties)+
  geom_tile(color="black")+
  scale_fill_gradient2("penalty",
                       high = "#b2182b",
                       low="#2166ac",
                       midpoint = 0,
                       breaks = c(0,0.5,1),
                       labels = c("0", "0.5", "1"))+
  labs(x="", y="")+
  coord_fixed()+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme(axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

basis_sample_plot = ggplot(data= tp_example_curve, aes(x, value))+
  geom_line(aes(y = Ftotal),size=2)+
  geom_line(aes(group = `basis function`,color= `basis function`),size=0.5)+ 
  geom_text(data= bs_func_labels,
            aes(label = basis_label, 
                color= `basis function`,
                y = value),
            parse=TRUE, hjust = 0,nudge_x = 0.01)+
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,1.15))+
  scale_color_manual(values = basis_func_palette)+
  guides(color = "none")

#This is a nessecary step to make sure the top left plot is aligned with the 
#bottom plot. See
# https://cran.r-project.org/web/packages/cowplot/vignettes/plot_grid.html
#for details
aligned_plots = align_plots(basis_func_plot, 
                            basis_sample_plot, 
                            align = 'v', 
                            axis = 'l')

top_row_plot = plot_grid(aligned_plots[[1]],basis_penalty_plot,
                         ncol=2,
                         rel_widths = c(1,1),labels= c("","B"))
  
full_plot = plot_grid(top_row_plot,aligned_plots[[2]],
                      nrow=2,
                      labels= c("A", "C"),
                      rel_heights = c(1,0.9))

full_plot


#### Code for III: What are hierarchical GAMs? ####


#The default CO2 plant variable is ordered;
#This recodes it to an unordered factor (see main text for why).
CO2 <- transform(CO2, Plant_uo=factor(Plant, ordered=FALSE))

#Loading simulated bird movement data
bird_move <- read.csv("data/bird_move.csv", stringsAsFactors = TRUE)

CO2_vis_plot <- ggplot(CO2, aes(x=conc, 
                                y=uptake, 
                                group=Plant,
                                color=Plant, 
                                lty=Plant)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = rep(c("red","blue","black"), times =4))+
  scale_linetype_manual(values = rep(1:4, each=3))+
  guides(color="none",linetype="none")+
  labs(x=expression(CO[2] ~ concentration ~ (mL ~ L^{-1})), 
       y=expression(CO[2] ~ uptake ~ (mu*mol ~ m^{-2})))

bird_vis_plot <- ggplot(dplyr::filter(bird_move, count > 0),
                        aes(x=week, y=latitude, size=count))+
  facet_wrap(~ species) +
  geom_point() +
  scale_size(name = "Count", range = c(0.2, 3)) +
  labs(x = "Week", y = "Latitude") +
  theme(legend.position = "bottom")

plot_grid(CO2_vis_plot, bird_vis_plot, nrow=1, labels=c("A","B"),
          align = "hv", axis = "lrtb")




CO2_modG <- gam(log(uptake) ~ s(log(conc), k=5, bs="tp") +
                  s(Plant_uo, k=12, bs="re"),
                data=CO2, method="REML", family="gaussian")
  
#plot the default gratia plot for the CO2 model
draw(CO2_modG)


# setup prediction data
CO2_modG_pred <- with(CO2,
                      expand.grid(conc=seq(min(conc), max(conc), length=100),
                                  Plant_uo=levels(Plant_uo)))

# make the prediction, add this and a column of standard errors to the prediction
# data.frame. Predictions are on the log scale.
CO2_modG_pred <- cbind(CO2_modG_pred,
                       predict(CO2_modG, 
                               CO2_modG_pred, 
                               se.fit=TRUE, 
                               type="response"))

# make the plot. Note here the use of the exp() function to back-transform the
# predictions (which are for log-uptake) to the original scale
ggplot(data=CO2, aes(x=conc, y=uptake, group=Plant_uo)) +
  facet_wrap(~Plant_uo) +
  geom_ribbon(aes(ymin=exp(fit - 2*se.fit), ymax=exp(fit + 2*se.fit), x=conc),
              data=CO2_modG_pred, 
              alpha=0.3, 
              inherit.aes=FALSE) +
  geom_line(aes(y=exp(fit)), data=CO2_modG_pred) +
  geom_point() +
  labs(x=expression(CO[2] ~ concentration ~ (mL ~ L^{-1})),
       y=expression(CO[2] ~ uptake ~ (mu*mol ~ m^{-2})))


# Note for specifying tensor products: you can either specify bs (basis) and k
# (number of basis functions) as single values, which would assign the same
# basis and k to each marginal value, or pass them as vectors, one value for
# each distinct marginal smoother (see ?mgcv::te for details)



bird_modG <- gam(count ~ te(week, latitude, bs=c("cc", "tp"), k=c(10, 10)),
                 data=bird_move, method="REML", family="poisson",
                 knots=list(week=c(0, 52)))


#gratia draw plot for the two-dimensional tensor product smoother for bird_modG.
draw(bird_modG)


#add the predicted values from the model to bird_move
bird_move <- transform(bird_move, modG = predict(bird_modG, type="response"))

ggplot(bird_move, aes(x=modG, y=count)) +
  facet_wrap(~species) +
  geom_point(alpha=0.1) +
  geom_abline() +
  labs(x="Predicted count", y="Observed count")


CO2_modGS <- gam(log(uptake) ~ s(log(conc), k=5, m=2) + 
                   s(log(conc), Plant_uo, k=5,  bs="fs", m=2),
                 data=CO2, method="REML")


#gratia draw() plot for CO2_modGS
draw(CO2_modGS)


CO2_modGS_pred <- predict(CO2_modGS, se.fit=TRUE)
CO2 <- transform(CO2, 
                 modGS = CO2_modGS_pred$fit, 
                 modGS_se = CO2_modGS_pred$se.fit)

ggplot(data=CO2, aes(x=conc, y=uptake, group=Plant_uo)) +
  facet_wrap(~Plant_uo) +
  geom_ribbon(aes(ymin=exp(modGS-2*modGS_se),
                  ymax=exp(modGS+2*modGS_se)), alpha=0.25) +
  geom_line(aes(y=exp(modGS))) +
  geom_point() +
  labs(x=expression(CO[2] ~ concentration ~ (mL ~ L^{-1})),
       y=expression(CO[2] ~ uptake ~ (mu*mol ~ m^{-2})))


bird_modGS <- gam(count ~ te(week, latitude, bs=c("cc", "tp"),
                             k=c(10, 10), m=2) +
                    t2(week, latitude, species, bs=c("cc", "tp", "re"),
                       k=c(10, 10, 6), m=2, full=TRUE),
                  data=bird_move, method="REML", family="poisson", 
                  knots=list(week=c(0, 52)))


bird_move <- transform(bird_move, modGS = predict(bird_modGS, type="response"))

bird_modGS_indiv <- ggplot(data=bird_move, 
                          aes(x=week, y=latitude, fill=modGS,color=modGS)) +
  geom_tile(size=0.25) +
  facet_wrap(~ species, ncol=6) +
  scale_fill_viridis("Count") +
  scale_color_viridis("Count") +
  scale_x_continuous(expand=c(0, 0), breaks=c(1, 26, 52)) +
  scale_y_continuous(expand=c(0, 0), breaks=c(0, 30, 60)) +
  labs(x = "Week", y = "Latitude") +
  theme(legend.position="right")

bird_modGS_indiv_fit <- ggplot(data=bird_move, aes(x=modGS, y=count)) +
  facet_wrap(~ species, ncol=6) +
  geom_point(alpha=0.1) +
  geom_abline() +
  labs(x="Predicted count (model *GS*)", y= "Observed count")

plot_grid(bird_modGS_indiv, bird_modGS_indiv_fit, 
          ncol=1, 
          align="vh", 
          axis = "lrtb",
          labels=c("A","B"), 
          rel_heights= c(1,1))

## 
## CO2_modGI <- gam(log(uptake) ~ s(log(conc), k=5, m=2, bs="tp") +
##                   s(log(conc), by=Plant_uo, k=5, m=1, bs="tp") +
##                   s(Plant_uo, bs="re", k=12),
##                 data=CO2, method="REML")


#Fitting CO2_modGI 
CO2_modGI <- gam(log(uptake) ~ s(log(conc), k=5, m=2, bs="tp") +
                   s(log(conc), by= Plant_uo, k=5, m=1, bs="tp") +
                   s(Plant_uo, bs="re", k=12),
                 data=CO2, method="REML")

#plotting CO2_modGI 
draw(CO2_modGI, select = c(1,14,8,2,11,5), scales = "fixed")


bird_modGI <- gam(count ~ species +
                    te(week, latitude, bs=c("cc", "tp"), k=c(10, 10), m=2) +
                    te(week, latitude, by=species, bs= c("cc", "tp"),
                       k=c(10, 10), m=1),
                 data=bird_move, method="REML", family="poisson",
                 knots=list(week=c(0, 52)))


CO2_modS <- gam(log(uptake) ~ s(log(conc), Plant_uo, k=5, bs="fs", m=2),
                data=CO2, method="REML")

bird_modS <- gam(count ~ t2(week, latitude, species, bs=c("cc", "tp", "re"),
                            k=c(10, 10, 6), m=2, full=TRUE),
                 data=bird_move, method="REML", family="poisson",
                 knots=list(week=c(0, 52)))


CO2_modI <- gam(log(uptake) ~ s(log(conc), by=Plant_uo, k=5, bs="tp", m=2) +
                  s(Plant_uo, bs="re", k=12),
                data=CO2, method="REML")


bird_modI <- gam(count ~ species + te(week, latitude, by=species,
                                      bs=c("cc", "tp"), k=c(10, 10), m=2),
                 data=bird_move, method="REML", family="poisson",
                 knots=list(week=c(0, 52)))


AIC_table <- AIC(CO2_modG,CO2_modGS, CO2_modGI, CO2_modS, CO2_modI,
             bird_modG, bird_modGS, bird_modGI, bird_modS, bird_modI)%>%
  rownames_to_column(var= "Model")%>%
  mutate(data_source = rep(c("CO2","bird_data"), each =5))%>%
  group_by(data_source)%>%
  mutate(deltaAIC = AIC - min(AIC))%>%
  ungroup()%>%
  dplyr::select(-data_source)%>%
  mutate_at(.vars = vars(df,AIC, deltaAIC), 
            .funs = funs(round,.args = list(digits=0)))


#### Code for IV: Examples ####


zooplankton <- read.csv("data/zooplankton_example.csv",stringsAsFactors = TRUE)%>%
  mutate(year_f = factor(year))

#This is what the data looks like:
str(zooplankton)
levels(zooplankton$taxon)
levels(zooplankton$lake)

# We'll now break it into testing and training data. The training data will be
# used to fit the model, and the testing data will be used to evaluate model
# fit.

#the first training and testing data set will be used to compare dynamics of
#plankton communities in Lake Mendota
zoo_train <- subset(zooplankton, year%%2==0 & lake=="Mendota")
zoo_test  <- subset(zooplankton, year%%2==1 & lake=="Mendota") 

#The second training and testing set will compare Daphnia mendotae dynamics
#among four lakes
daphnia_train <- subset(zooplankton, year%%2==0 & taxon=="D. mendotae")
daphnia_test  <- subset(zooplankton, year%%2==1 & taxon=="D. mendotae")

#This function calculates the deviance of out-of-sample data,
#conditional on their mean predicted value from the model
get_deviance <- function(model, y_pred, y_obs, weights = NULL){
  stopifnot(length(y_obs)==length(y_pred))
  #We don't use the weights term in this paper, but it can be useful if
  #how well the model matters more for some sample points than others
  if(is.null(weights)) weights = rep(1, times= length(y_obs))
  #this uses the deviance residual function from the model family to
  #calculate deviances for individual points
  dev_residuals = model$family$dev.resids(y_obs, y_pred, weights)
  return(sum(dev_residuals))
}


zoo_comm_modS <- gam(density_adj ~ s(taxon, year_f, bs="re") +
                       s(day, taxon, bs="fs", k=10, xt=list(bs="cc")),
                     data=zoo_train, knots=list(day=c(0, 365)),
                     family=Gamma(link="log"), method="REML",
                     drop.unused.levels=FALSE)


# Note that s(taxon, bs="re") has to be explicitly included here, as the 
# day by taxon smoother does not include an intercept
zoo_comm_modI <- gam(density_adj ~ s(day, by=taxon, k=10, bs="cc") + 
                       s(taxon, bs="re") + s(taxon, year_f, bs="re"),
                     data=zoo_train, knots=list(day=c(0, 365)),
                     family=Gamma(link="log"), method="REML",
                     drop.unused.levels=FALSE)

## 
## gam.check(zoo_comm_modS)
## gam.check(zoo_comm_modI)


#Checking residuals and qqplots for GAM fits

#QQ-plot, using gratia's qq_plot function, with simulated confidence intervals.
#We are removing the title and subtitle to simplify the figure
plt1 <- qq_plot(zoo_comm_modI, method = "simulate") +
  labs(title =NULL, subtitle =NULL)
df <- data.frame(log_fitted = log(fitted(zoo_comm_modI)),
                 residuals  = resid(zoo_comm_modI, type = "deviance"))

#fitted versus deviance plot
plt2 <- ggplot(df, aes(x = log_fitted, y = residuals)) +
    geom_point() +
    labs(x = "Linear predictor", y = "Deviance residual")
plot_grid(plt1, plt2, ncol = 2, align = "hv", axis = "lrtb",labels=c("A","B"))

## 
## #individual components of gam.check: the results for k.check
## round(k.check(zoo_comm_modI),2)




#Create synthetic data to use to compare predictions
zoo_plot_data <- expand.grid(day = 1:365, 
                             taxon = factor(levels(zoo_train$taxon)), 
                             year_f = 1980)

#extract predicted values and standard errors for both models. 
#the exclude = "s(taxon,year_f)" term indicates that predictions should be made
#excluding the effect of the taxon by year random effect (effectively setting
#making predictions averaging over year-taxon effects).
zoo_modS_fit <- predict(zoo_comm_modS, 
                        zoo_plot_data, 
                        se.fit = TRUE, 
                        exclude = "s(taxon,year_f)")
zoo_modI_fit <- predict(zoo_comm_modI, 
                        zoo_plot_data, 
                        se.fit = TRUE, 
                        exclude = "s(taxon,year_f)")

zoo_plot_data$modS_fit <- as.numeric(zoo_modS_fit$fit)
zoo_plot_data$modI_fit <- as.numeric(zoo_modI_fit$fit)

zoo_plot_data <- gather(zoo_plot_data, model, fit, modS_fit, modI_fit)
zoo_plot_data <- mutate(zoo_plot_data, se= c(as.numeric(zoo_modS_fit$se.fit),
                                             as.numeric(zoo_modI_fit$se.fit)),
                        upper = exp(fit + (2 * se)),
                        lower = exp(fit - (2 * se)),
                        fit   = exp(fit))

#Plot the model output, with means plus standard deviations for each model.
zoo_plot_model_labels = paste("Model", c("S","I"))
zoo_plot_model_labels = factor(zoo_plot_model_labels, 
                               levels = zoo_plot_model_labels)

zoo_plot <- ggplot(zoo_plot_data) +
  facet_wrap(~taxon, nrow = 4,scales = "free_y")+
  geom_ribbon(aes(x=day,
                  ymin = lower,
                  ymax = upper,
                  fill = model),
              alpha=0.2)+
  geom_point(data= zoo_train, aes(x = day, y = density_adj),size=0.06)+
  geom_point(data= zoo_test, aes(x = day, y = density_adj),
             size=0.06,col="grey")+
  geom_line(aes(x = day, y = fit, color = model))+
  labs(y = expression(atop(Population~density,("10 000"~individuals~m^{-2}))), 
       x = "Day of Year") +
  scale_y_log10(breaks = c(0.1,1,10,100, 1000), 
                labels = c("0.1","1","10","100","1000"))+
  scale_fill_brewer(name = "", palette = "Dark2",
                      labels = zoo_plot_model_labels) +
  scale_colour_brewer(name = "",
                        palette = "Dark2", labels = zoo_plot_model_labels)+
  theme(legend.position = "top")

zoo_plot


#Getting the out-of-sample predictions for both models:

# we need to compare how well this model fits with a null model. here we'll use an
# intercept-only model
zoo_comm_mod0 <- gam(density_adj ~ s(taxon,bs="re"),
                     data=zoo_train,
                     knots = list(day =c(0, 365)),
                     family = Gamma(link ="log"), 
                     method = "REML",
                     drop.unused.levels = FALSE)

#Correlations between fitted and observed values for all species:
#\n is in variable titles to add a line break in the printed table.
zoo_test_summary = zoo_test %>%
  mutate(
    mod0 = predict(zoo_comm_mod0, ., type="response"),
    modS = predict(zoo_comm_modS, ., type="response"),
    modI = predict(zoo_comm_modI, ., type="response"))%>%
  group_by(taxon)%>%
  summarise(
    `Intercept only` = format(get_deviance(zoo_comm_mod0, mod0, density_adj), 
                              scientific = FALSE, 
                              digits=3),
    `Model S` = format(get_deviance(zoo_comm_modS, modS, density_adj), 
                       scientific = FALSE, 
                       digits=3),
    `Model I` = format(get_deviance(zoo_comm_modI, modI, density_adj), 
                       scientific = FALSE, 
                       digits=3))





zoo_daph_modG <- gam(density_adj ~ s(day, bs="cc", k=10) + s(lake, bs="re") +
                       s(lake, year_f, bs="re"),
                     data=daphnia_train, knots=list(day=c(0, 365)),
                     family=Gamma(link="log"), method="REML",
                     drop.unused.levels=FALSE)


zoo_daph_modGS <- gam(density_adj ~ s(day, bs="cc", k=10) +
                        s(day, lake, k=10, bs="fs", xt=list(bs="cc")) +
                        s(lake, year_f, bs="re"),
                      data=daphnia_train, knots=list(day=c(0, 365)),
                      family=Gamma(link="log"), method="REML",
                      drop.unused.levels=FALSE)


zoo_daph_modGI <- gam(density_adj~s(day, bs="cc", k=10) + s(lake, bs="re") +
                        s(day, by=lake, k=10, bs="cc") +
                        s(lake, year_f, bs="re"),
                      data=daphnia_train, knots=list(day=c(0, 365)),
                      family=Gamma(link ="log"), method="REML",
                      drop.unused.levels=FALSE)


#Checking residuals and qqplots for GAM fits

#qqplot, using gratia's qq_plot function, with simulated confidence intervals
pltG <- qq_plot(zoo_daph_modG, method = "simulate")+
  labs(subtitle = NULL, title=NULL)
pltGS <- qq_plot(zoo_daph_modGS, method = "simulate")+
  labs(subtitle = NULL, title=NULL, y=NULL)
pltGI <- qq_plot(zoo_daph_modGI, method = "simulate")+
  labs(subtitle = NULL, title=NULL, y=NULL)

plot_grid(pltG, pltGS,pltGI, 
          ncol = 3, 
          align = "hv", 
          axis = "lrtb",labels=c("A","B","C"))


#Create synthetic data to use to compare predictions
daph_plot_data <- expand.grid(day = 1:365, 
                              lake = factor(levels(zoo_train$lake)),
                              year_f = 1980)

#extract predicted values and standard errors for both models. the 
#exclude ="s(taxon,year_f)" term indicates that predictions should be made 
#excluding the effect of the taxon-by-year random effect (effectively making
#predictions averaging over year-taxon effects).
daph_modG_fit <- predict(zoo_daph_modG, 
                         newdata = daph_plot_data, 
                         se.fit = TRUE, 
                         exclude = "s(lake,year_f)")
daph_modGS_fit <- predict(zoo_daph_modGS, 
                         newdata = daph_plot_data, 
                         se.fit = TRUE, 
                         exclude = "s(lake,year_f)")
daph_modGI_fit <- predict(zoo_daph_modGI, 
                         newdata = daph_plot_data, 
                         se.fit = TRUE, 
                         exclude = "s(lake,year_f)")

daph_plot_data$modG_fit <- as.numeric(daph_modG_fit$fit)
daph_plot_data$modGS_fit <- as.numeric(daph_modGS_fit$fit)
daph_plot_data$modGI_fit <- as.numeric(daph_modGI_fit$fit)

daph_plot_data <- gather(daph_plot_data, 
                         key = model, 
                         value = fit, 
                         modG_fit, 
                         modGS_fit, 
                         modGI_fit)

daph_plot_data <- mutate(daph_plot_data, 
                         se = c(as.numeric(daph_modG_fit$se.fit),
                                as.numeric(daph_modGS_fit$se.fit),
                                as.numeric(daph_modGI_fit$se.fit)),
                         upper = exp(fit + (2 * se)),
                         lower = exp(fit - (2 * se)),
                         fit   = exp(fit))

daph_plot_model_labels = paste("Model", c("G","GS","GI"))
daph_plot_model_labels = factor(daph_plot_model_labels, 
                                levels= daph_plot_model_labels)

daph_plot <- ggplot(daph_plot_data, aes(x=day))+
  facet_wrap(~lake, nrow = 2)+
  geom_ribbon(aes(x = day, ymin = lower, ymax = upper, fill = model), 
                alpha = 0.2) +
  geom_point(data= daphnia_train, 
             aes(x = day, 
                 y = density_adj),
             size=0.06)+
  geom_point(data= daphnia_test, 
             aes(x = day, 
                 y = density_adj),
             size=0.06,
             col="grey")+
  geom_line(aes(x = day, y = fit, colour = model)) +

  labs(y = expression(atop(Population~density,
                           ("10 000"~individuals~m^{-2}))), 
       x = "Day of Year") +
  scale_x_continuous(expand = c(0,0))+
  scale_y_log10()+
    scale_fill_brewer(name = "", 
                      palette = "Dark2",
                      labels = daph_plot_model_labels) +
    scale_colour_brewer(name = "",
                        palette = "Dark2", 
                        labels = daph_plot_model_labels)+
  theme(legend.position = "top")


daph_plot


# we need to compare how well this model fits with a null model. here we'll use
# an intercept-only model
zoo_daph_mod0 <- gam(density_adj~s(lake, bs="re"),
                     data=daphnia_train,
                     knots=list(day =c(0, 365)),
                     family=Gamma(link ="log"),
                     method="REML",
                     drop.unused.levels = FALSE)



# We'll look at the correlation between fitted and observed values for all species:

daph_test_summary <- daphnia_test %>%
  mutate(
    #get out-of-sample predicted fits
    mod0 = as.numeric(predict(zoo_daph_mod0,.,type="response")),
    modG = as.numeric(predict(zoo_daph_modG,.,type="response")),
    modGS = as.numeric(predict(zoo_daph_modGS,.,type="response")),
    modGI = as.numeric(predict(zoo_daph_modGI,.,type="response")))%>%
  group_by(lake)%>%
  summarise(`Intercept only` = format(get_deviance(zoo_daph_mod0, 
                                                   mod0, 
                                                   density_adj), 
                                      scientific = FALSE, 
                                      digits=2),
            `Model G` = format(get_deviance(zoo_daph_modG, 
                                            modG, 
                                            density_adj), 
                               scientific = FALSE, 
                               digits=2),
            `Model GS` = format(get_deviance(zoo_daph_modGS, 
                                             modGS, 
                                             density_adj), 
                               scientific = FALSE, 
                               digits=2),
            `Model GI` = format(get_deviance(zoo_daph_modGI, 
                                             modGI, 
                                             density_adj), 
                               scientific = FALSE, 
                               digits=2))%>%
  rename(Lake = lake)




#### Code for V: Computational and statistical issues when fitting HGAMs ####


#This code will generate the bias-variance tradeoff plot
set.seed(1)

#Calculate the numerical approximation of the 2nd derivative for a a function
#given x and y values. Assumes a constant step size (delta) for x along its path.
calc_2nd_deriv = function(x,y){
  deriv_val = (lag(y) + lead(y) - 2*y)/(x-lag(x))^2
  deriv_val
}

#Generate true regression functions that differ in their frequencies. 
#Higher frequencies correspond to more variable functions. 
freq_vals = c(1/2,1,2,4)
n_reps = 25
noise_levels = c(0.5,1,2)
biasvar_data = crossing(noise = noise_levels,
                        rep = 1:n_reps,
                        x = seq(0,2*pi,length=150),
                        freq = freq_vals
               )%>%
  mutate(y = cos(freq*x) +rnorm(n(), 0, noise),
         grp = paste("frequency = ",freq,sep= ""),
         grp = factor(grp,  
                      levels = paste("frequency = ",freq_vals,sep= "")))

biasvar_fit = biasvar_data %>%
  group_by(noise,rep)%>%
  do(
    #Fit model S (shared smoothness) for the test data
    modS = bam(y~s(x,k=30,grp, bs="fs"), data=.),
    #Fit a model I function (differing smoothness) for the test data
    modI = bam(y~s(x,k=30,by=grp)+s(grp,bs="re"), data=.)
    )

#Extract fitted values for each model for all the test data
biasvar_predict_data = crossing(x = seq(0,2*pi,length=500), 
                                freq = freq_vals)%>%
  mutate(grp = paste("frequency = ",freq,sep= ""),
         grp = factor(grp,  
                      levels = paste("frequency = ",freq_vals,sep= "")),
         y = cos(freq*x))

biasvar_predict_fit = biasvar_fit %>%
  group_by(noise,rep)%>%
  do(fitS = as.numeric(predict(.$modS[[1]],
                               newdata = biasvar_predict_data,
                               type = "response")),
     fitI = as.numeric(predict(.$modI[[1]],
                               newdata = biasvar_predict_data,
                               type="response")))%>%
  unnest(fitS, fitI) %>%
  bind_cols(crossing(noise=noise_levels, 
                     rep= 1:n_reps,
                     biasvar_predict_data))

#turn this into long-format data for plotting, and to make it easier to
#calculate derivatives
biasvar_predict_fit_long = biasvar_predict_fit %>%
  gather(model, value, y, fitS, fitI)%>%
  mutate(model = recode(model, y = "true value",fitS = "model S fit",
                        fitI = "model I fit"),
         model = factor(model, levels=  c("true value","model S fit",
                                          "model I fit")))

biasvar_predict_fit_summary = biasvar_predict_fit_long %>%
  group_by(noise,grp,model, x)%>%
  summarize(lower = min(value),
            upper = max(value),
            value = mean(value)
            )
#estimate 2nd derivatives of each curve, then for each curve calculate the sum
#of squared 2nd derivatives of the true curve and the predictions for both
#models.
deriv_est_data = biasvar_predict_fit %>%
  group_by(grp, rep,noise)%>%
  arrange(grp, x)%>%
  mutate(fitS_deriv = calc_2nd_deriv(x,fitS),
         fitI_deriv = calc_2nd_deriv(x,fitI))%>%
  summarize(freq = freq[1], fitS_int = sum(fitS_deriv^2*(x-lag(x)),
                                           na.rm = TRUE),
            fitI_int = sum(fitI_deriv^2*(x-lag(x)),na.rm = TRUE))%>%
  ungroup()%>%
  mutate(sqr_2nd_deriv = freq^3*(sin(4*pi*freq)+4*pi*freq)/4)%>%
  gather(key=model,value = obs_sqr_deriv,fitS_int,fitI_int)%>%
  mutate(model = factor(ifelse(model=="fitS_int", "model S fit",
                               "model I fit"),
                        levels = c("model S fit",
                                   "model I fit")))



# Uses the Tex function from the latex2exp package to create a math label for 
# the facets. Based off code from
# https://sahirbhatnagar.com/blog/2016/facet_wrap_labels/
noise_labeller <- function(string) {
  signal_to_noise = as.numeric(string)
  signal_to_noise = 0.5/signal_to_noise^2
  signal_to_noise = as.character(signal_to_noise)
  TeX(paste("signal:noise =", signal_to_noise)) 
}

#The derivative plots

deriv_min = -3
deriv_plot =  ggplot(data=deriv_est_data, 
                     aes(x=sqr_2nd_deriv,     
                         y= pmax(obs_sqr_deriv,10^deriv_min),
                         fill= model,
                         group=paste(sqr_2nd_deriv,model)))+
  facet_grid(.~noise, 
             labeller = as_labeller(noise_labeller,
                                    default = label_parsed))+
  geom_dotplot(binaxis = "y", 
               stackdir = "center",
               binwidth = 0.125,
               color=NA)+
  scale_y_log10("Fitted wiggliness",
                limits = c(10^deriv_min,5e+3),expand=c(0,0.1),
                breaks = c(10^deriv_min,  1e+0, 1e+3), 
                labels = list(bquote(""<=10^.(deriv_min)),
                              bquote(10^0), 
                              bquote(10^3))
                )+
  scale_x_log10("True wiggliness", 
                breaks= c(1e-3, 1e+0, 1e+3), 
                labels = list(bquote(10^-3),bquote(10^0), bquote(10^3))
                )+
  scale_fill_brewer(name=NULL,palette= "Set1")+
  geom_abline(color="black")+
  theme(legend.position = "top",
        strip.text.x = element_blank(),
        plot.margin = unit(c(3, 5.5, 5,5, 5.5), "pt"))

fit_colors = c("black",RColorBrewer::brewer.pal(3, "Set1")[1:2])

overfit_vis_plot = ggplot(data=biasvar_predict_fit_summary,
                          aes(x=x,y= value,color=model))+
  facet_grid(grp~noise, labeller = as_labeller(noise_labeller,
                                               default = label_parsed))+
  geom_line(data=filter(biasvar_predict_fit_long,
                        rep==1,
                        model=="true value"),
            color= "black",
            size=0.9)+
  geom_line(data=filter(biasvar_predict_fit_long,rep%in%1:3),
            aes(group=paste(rep,model)))+
  scale_color_manual(name = "",values=fit_colors)+
  scale_y_continuous("Estimated curve", breaks = c(-2,0,2))+
  scale_x_continuous(name = "x")+
  coord_cartesian(ylim=c(-2,2))+
  theme(legend.position = "top",
        strip.text.y = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 2.5, 5.5), "pt"))

#Plot the overfit graphs together.
cowplot::plot_grid(overfit_vis_plot, deriv_plot, ncol=1, labels="AUTO",
                   align="hv", axis="lr",
                   rel_heights = c(1,0.5))


#Note: this code takes quite a long time to run, as it's fitting all 10 models.

#This function extracts the number of penalties used in a model. For Gamma and
#Gaussian families, you need to remove the penalty (i.e. variance) associated
#with the scale term
get_n_pen  = function(model) {
  family = model$family[[1]]
  if(family %in% c("Gamma","gaussian")){
    capture.output({out_val = nrow(gam.vcomp(model))-1})
  }else{
    capture.output({out_val = nrow(gam.vcomp(model))})
  }
  return(out_val)
}

#Extract the number of coefficients from a fitted GAM
get_n_coef = function(model) length(coef(model))

#Get the number of inner and outer iterations needed to fit the final model
get_n_iter = function(model) model$outer.info$iter
get_n_out_iter = function(model) model$iter

#combine results into a single table
comp_resources = crossing(model_number = factor(c("G","GS","GI","S","I"),
                                                levels = c("G","GS","GI","S","I")),     
                          data_source = factor(c("CO2","bird_move"),
                                               levels = c("CO2","bird_move")),
                          time = 0, n_smooths = 0,
                          n_coef = 0)


#Fit each model to the example data sets, and calculate run time for them
comp_resources[1,"time"] = system.time(
  CO2_modG <- gam(log(uptake) ~ s(log(conc),k=5,m=2, bs="tp")+
                                s(Plant_uo, k =12,  bs="re"),
                  data= CO2,
                  method="REML",
                  control = list(keepData=TRUE))
  )[3]

comp_resources[2,"time"] = system.time(
  bird_modG <- gam(count ~ te(week,latitude, bs= c("cc", "tp"), k=c(10,10)),
                   data= bird_move, 
                   method="REML", 
                   family= poisson,
                   knots=list(week=c(0, 52)),
                   control = list(keepData=TRUE))
  )[3]

comp_resources[3,"time"] = system.time(
  CO2_modGS <- gam(log(uptake) ~ s(log(conc),k=5,m=2, bs="tp")+
                                 s(log(conc), Plant_uo, k=5, bs="fs",m=1),
                  data= CO2,
                  method="REML",
                  control = list(keepData=TRUE))
  )[3]


comp_resources[4,"time"] = system.time(
  bird_modGS <- gam(count ~ te(week,latitude, bs= c("cc", "tp"),
                              k=c(10,10),m=2)+
                    t2(week, latitude, species, bs=c("cc", "tp", "re"),
                       k=c(10, 10, 6), m=2, full=TRUE),
                   data= bird_move, 
                   method="REML", 
                   family= poisson,
                   control = list(keepData=TRUE))
  )[3]


comp_resources[5,"time"] = system.time(
  CO2_modGI <- gam(log(uptake) ~ s(log(conc),k=5,m=2, bs="tp")+
                                 s(log(conc),by= Plant_uo, k =5,  bs="ts",m=1)+
                                 s(Plant_uo,bs="re",k=12),
                  data= CO2,
                  method="REML",
                  control = list(keepData=TRUE)))[3]



comp_resources[6,"time"] = system.time(
  bird_modGI <- gam(count ~ te(week,latitude, bs= c("cc", "tp"),
                               k=c(10,10),m=2) +
                            te(week,latitude, bs= c("cc", "tp"),
                               k=c(10,10),m=1,by= species),
                   data= bird_move, 
                   method="REML", 
                   family= poisson,
                   control = list(keepData=TRUE)))[3]


comp_resources[7,"time"] = system.time(
  CO2_modS <- gam(log(uptake) ~ s(log(conc), Plant_uo, k=5,  bs="fs",m=2),
                  data= CO2,
                  method="REML",
                  control = list(keepData=TRUE))
  )[3]


comp_resources[8,"time"] = system.time(
  bird_modS <- gam(count ~ t2(week, latitude, species, bs=c("cc", "tp", "re"),
                            k=c(10, 10, 6), m=2, full=TRUE),
                   data= bird_move, 
                   method="REML", 
                   family= poisson,
                   control = list(keepData=TRUE))
)[3]


comp_resources[9,"time"] = system.time(
  CO2_modI <- gam(log(uptake) ~ s(log(conc),by= Plant_uo, k =5,  bs="tp",m=2) +
                                s(Plant_uo,bs="re",k=12),
                  data= CO2,
                  method="REML",
                  control = list(keepData=TRUE))
)[3]

comp_resources[10,"time"] = system.time(
  bird_modI <- gam(count ~ te(week,latitude,by=species, bs= c("cc", "tp"),
                              k=c(10,10),m = 2),
                   data= bird_move, 
                   method="REML", 
                   family= poisson,
                   control = list(keepData=TRUE))
)[3]

#combine all fitted models into a list
comp_resources$model = list(CO2_modG, bird_modG, CO2_modGS, bird_modGS,
                            CO2_modGI, bird_modGI,CO2_modS, bird_modS,
                            CO2_modI, bird_modI)

#Extract all of the information on computer time and resources needed for each
#model
comp_resources = comp_resources %>%
  group_by(model_number, data_source)%>%
  mutate(n_smooths = get_n_pen(model[[1]]),
         n_coef = get_n_coef(model[[1]]),
         n_iter = get_n_iter(model[[1]]),
         n_iter_out = get_n_out_iter(model[[1]]))



comp_resources_table =comp_resources %>%
  ungroup()%>%
  arrange(data_source,model_number)%>%
  transmute(data_source =data_source, Model=model_number,
            `Relative Time` = time,`Coefficients` = n_coef,
            `Penalties` = n_smooths
            )%>%
  group_by(data_source) %>%
  mutate(#scales processing time relative to model *G*
         `Relative Time` = `Relative Time`/`Relative Time`[1],
         #rounds to illustrate differences in timing.
         `Relative Time` = ifelse(`Relative Time`<10, 
                                  signif(`Relative Time`,1), 
                                  signif(`Relative Time`, 2)) 
         )%>%
  ungroup() %>%
  dplyr::select( - data_source)




#This code calculates the timing it takes to fit the same model (with varying
#amounts of data) for gam, bam, gamm, and gamm4.

# ensures that each new model parameter set is an extension of the old one
set.seed = 1

n_x = 20
x = seq(-2,2, length=n_x)
n_steps = 7

#setting up blank data frame to put results in.
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

#for each number of observations, fit all the models and calculate run times for
#them
for(i in 1:n_steps){

  n_g =  fit_timing_data$n_groups[i]
  
  fac_current = fac_all[1:n_g]
  fac_current = factor(fac_current, levels=  unique(fac_current))

  model_coefs = model_coefs_all%>%
    filter(fac %in% fac_current)%>%
    mutate(fac = factor(fac, levels= unique(fac)))
  
  #ensures that each new data set is an extension of the old one
  set.seed = 1 
  
  model_data = crossing(fac=fac_current, x=x)%>%
    left_join(model_coefs)%>%
    mutate(base_func  = dnorm(x)*10,
           indiv_func = int + x^2*x2 +
                        2*(exp(x*logit_slope)/(1+exp(x*logit_slope))-0.5),
           y = base_func + indiv_func + rnorm(n()))
  
  fit_timing_data$gam[i] = system.time(gam(y~s(x,k=10, bs="cp") + 
                                             s(x,fac, k=10, bs="fs", 
                                               xt=list(bs="cp"), m=2),
                             data= model_data, method="REML")
                             )[3]
  
  fit_timing_data$`bam (discrete = FALSE)`[i] = system.time(
    bam(y~s(x,k=10, bs="cp") + 
          s(x,fac, k=10, bs="fs", xt=list(bs="cp"),m=2),
        data= model_data, 
        discrete = FALSE)
    )[3]
  
  fit_timing_data$`bam (discrete = TRUE)`[i] = system.time(
    bam(y~s(x,k=10, bs="cp") + 
          s(x,fac, k=10, bs="fs", xt=list(bs="cp"),m=2),
        data= model_data,
        discrete=TRUE)
    )[3]
  
  
  fit_timing_data$gamm[i] = system.time(
    gamm(y~s(x,k=10, bs="cp") + 
           s(x,fac, k=10, bs="fs", xt=list(bs="cp"),m=2),
         data= model_data)
    )[3]
  
  
  fit_timing_data$gamm4[i] = system.time(
    gamm4(y~s(x,k=10, bs="cp") + 
            s(x,fac, k=10, bs="fs", xt=list(bs="cp"),m=2),
          data= model_data)
    )[3]
}

#Combine all of the timing data into long format, ready for plotting
fit_timing_long = fit_timing_data %>% 
  gather(model, 
         timing, 
         gam,
         `bam (discrete = FALSE)`, 
         `bam (discrete = TRUE)`, gamm, gamm4)%>%
  mutate(model = factor(model, levels = c("gam",
                                         "bam (discrete = FALSE)",
                                         "bam (discrete = TRUE)",
                                         "gamm", 
                                         "gamm4")))


timing_plot = ggplot(aes(n_groups, timing, color=model, linetype= model), 
                     data=fit_timing_long)+
  geom_line()+
  geom_point(show.legend = FALSE)+
  scale_color_manual(name = "Model",
                     values = c("black", 
                                "#1b9e77",
                                "#1b9e77", 
                                "#d95f02", 
                                "#7570b3"))+
  scale_linetype_manual(name = "Model", values =c(1,1,2,1,1)) +
  scale_y_log10("Run time (seconds)", 
                breaks = c(0.1,1,10,100), 
                labels = c("0.1", "1","10", "100"))+
  scale_x_log10("Number of groups", 
                breaks = c(2,8,32,128))+
  guides(color = guide_legend(nrow = 2, byrow = TRUE))+
  theme(legend.position = "top")
timing_plot


#Load the global function for the bird_move dataset
bird_move_global <- read.csv("data/bird_move_global.csv")


#Simple te() model 
bird_modGS_te <- gam(count ~ te(week, latitude, bs=c("cc", "tp"),
                                k=c(10, 10), m=2) +
                             te(week, latitude, species, bs=c("cc", "tp", "re"),
                                k=c(10, 10, 6), m=2),
                     data=bird_move, method="REML", family="poisson", 
                     knots = list(week = c(0, 52)))

bird_mod_predict = bird_move_global %>%
  mutate(species = "sp1")%>%
  mutate(`te` = predict(bird_modGS_te,
                                               newdata = ., 
                                               type="terms")[,1],
         `t2` =  predict(bird_modGS,
                                             newdata = ., 
                                             type="terms")[,1])%>%
  mutate(`true global function` = `global_scaled_function`)%>%
  gather(key = model, value =`fitted value`, `te`:`true global function`)%>%
  mutate(model = factor(model, 
                        levels = c("true global function",
                                   "te",
                                   "t2")))

bird_global_labels <- data_frame(week = 5, latitude = 55, 
                                 label = c("A", "B","C"),
                                 model = c("true global function",
                                   "te",
                                   "t2"))%>%
  mutate(model = factor(model, 
                        levels = c("true global function",
                                   "te",
                                   "t2")))

bird_global_fitted_plot <- ggplot(bird_mod_predict, 
                                aes(x=week, 
                                    y=latitude))+
  facet_grid(~model)+
  geom_raster(aes(fill = `fitted value`))+
  geom_text(data= bird_global_labels, aes(label =label),size=9)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_viridis_c(name="linear predictor")+
  theme(legend.position = "bottom")

print(bird_global_fitted_plot)
