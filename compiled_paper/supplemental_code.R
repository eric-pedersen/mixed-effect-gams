####Setup ####
#all the packages needed for this tutorial are listed here
library(mgcv)
library(MASS)
library(stringr)
library(gamm4)
library(tidyr)
library(ggplot2)
library(viridis)
library(cowplot)
library(kableExtra)
library(docxtools)
library(knitr)
library(tibble)
library(dplyr)

#Set the default theme for ggplot objects to theme_bw()
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

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

plot_grid(p1, p2, p3, align = "hv", axis = "lrtb", ncol = 3, labels = "auto")

#Code for generating figure 3: examples of basis functions and splines ####
k = 6
plotting_data = data.frame(x = seq(0,1,length=100))

#This creates the basis functions for a thin-plate spline
#The absorb.cons=FALSE setting makes sure that smoothCon does not remove basis 
#functions that have a non-zero sum (in this case, the intercept). Absorbing 
#constraints would result in having less than k basis functions, which is why fitted 
#terms in mgcv often have less than k maximum EDF. 
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
set.seed(6) #ensure we can reproduce this example
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
            parse=TRUE, hjust = 0,nudge_x = 0.01)+
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,1.15))+
  scale_color_viridis_d()+
  guides(color = "none")


top_row_plot = plot_grid(basis_func_plot,basis_penalty_plot,
                         ncol=2,
                         rel_widths = c(1.,1),labels= c("","b"))
  
full_plot = plot_grid(top_row_plot,basis_sample_plot,
                      nrow=2,
                      labels= c("a", "c"),
                      axis = "t",
                      rel_heights = c(1,0.9))

full_plot

#### Code for III: What are hierarchical GAMs? ####

#The default CO2 plant variable is ordered;
#This recodes it to an unordered factor (see main text for why).
CO2 <- transform(CO2, Plant_uo=factor(Plant, ordered=FALSE))

#Loading simulated bird movement data
bird_move <- read.csv("../data/bird_move.csv")

CO2_vis_plot <- ggplot(CO2, aes(x=conc, y=uptake, group=Plant,color=Plant, lty=Plant)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = rep(c("red","blue","black"), times =4))+
  scale_linetype_manual(values = rep(1:4, each=3))+
  guides(color="none",linetype="none")+
  labs(x=expression(CO[2] ~ concentration ~ (mL ~ L^{-1})), y=expression(CO[2] ~ uptake ~ (mu*mol ~ m^{-2})))

bird_vis_plot <- ggplot(dplyr::filter(bird_move, count > 0),
                        aes(x=week, y=latitude, size=count))+
  facet_wrap(~ species) +
  geom_point() +
  scale_size(name = "Count", range = c(0.2, 3)) +
  labs(x = "Week", y = "Latitude") +
  theme(legend.position = "bottom")

plot_grid(CO2_vis_plot, bird_vis_plot, nrow=1, labels=c("a","b"),
          align = "hv", axis = "lrtb")
CO2_mod1 <- gam(log(uptake) ~ s(log(conc), k=5, bs="tp") +
                  s(Plant_uo, k=12, bs="re"),
                data=CO2, method="REML", family="gaussian")

plot(CO2_mod1, pages=1, seWithMean=TRUE)

# setup prediction data
CO2_mod1_pred <- with(CO2,
                      expand.grid(conc=seq(min(conc), max(conc), length=100),
                                  Plant_uo=levels(Plant_uo)))

# make the prediction, add this and a column of standard errors to the prediction
# data.frame. Predictions are on the log scale.
CO2_mod1_pred <- cbind(CO2_mod1_pred,
                       predict(CO2_mod1, CO2_mod1_pred, se.fit=TRUE, type="response"))

# make the plot. Note here the use of the exp() function to back-transform the
# predictions (which are for log-uptake) to the original scale
ggplot(data=CO2, aes(x=conc, y=uptake, group=Plant_uo)) +
  facet_wrap(~Plant_uo) +
  geom_ribbon(aes(ymin=exp(fit - 2*se.fit), ymax=exp(fit + 2*se.fit), x=conc),
              data=CO2_mod1_pred, alpha=0.3, inherit.aes=FALSE) +
  geom_line(aes(y=exp(fit)), data=CO2_mod1_pred) +
  geom_point() +
  labs(x=expression(CO[2] ~ concentration ~ (mL ~ L^{-1})),
       y=expression(CO[2] ~ uptake ~ (mu*mol ~ m^{-2})))
## # Note for specifying tensor products: you can either specify bs (basis) and
## # k (number of basis functions) as single values, which would assign the same
## # basis and k to each marginal value, or pass them as vectors, one value for each
## # distinct marginal smoother (see ?mgcv::te for details)
bird_mod1 <- gam(count ~ te(week, latitude, bs=c("cc", "tp"), k=c(10, 10)),
                 data=bird_move, method="REML", family="poisson",
                 knots = list(week = c(0, 52)))
#mgcv gam plot for the two-dimensional tensor product smoother for bird_mod1.
#scheme=2 displays the color scheme (rather than mgcv's default, which only
#shows contour lines)
plot(bird_mod1, pages=1, scheme=2, rug=FALSE)
box()
#add the predicted values from the model to bird_move
bird_move <- transform(bird_move, mod1 = predict(bird_mod1, type="response"))

ggplot(bird_move, aes(x=mod1, y=count)) +
  facet_wrap(~species) +
  geom_point(alpha=0.1) +
  geom_abline() +
  labs(x="Predicted count", y= "Observed count")
CO2_mod2 <- gam(log(uptake) ~ s(log(conc), k=5, m=2) + 
                  s(log(conc), Plant_uo, k=5,  bs="fs", m=2),
                data=CO2, method="REML")
plot(CO2_mod2, page=1, seWithMean=TRUE)
CO2_mod2_pred <- predict(CO2_mod2, se.fit=TRUE)
CO2 <- transform(CO2, mod2 = CO2_mod2_pred$fit, mod2_se = CO2_mod2_pred$se.fit)

ggplot(data=CO2, aes(x=conc, y=uptake, group=Plant_uo)) +
  facet_wrap(~Plant_uo) +
  geom_ribbon(aes(ymin=exp(mod2-2*mod2_se),
                  ymax=exp(mod2+2*mod2_se)), alpha=0.25) +
  geom_line(aes(y=exp(mod2))) +
  geom_point() +
  labs(x=expression(CO[2] ~ concentration ~ (mL ~ L^{-1})),
       y=expression(CO[2] ~ uptake ~ (mu*mol ~ m^{-2})))
bird_mod2 <- gam(count ~ te(week, latitude, bs=c("cc", "tp"),
                            k=c(10, 10), m=2) +
                   t2(week, latitude, species, bs=c("cc", "tp", "re"),
                      k=c(10, 10, 6), m=2, full=TRUE),
                 data=bird_move, method="REML", family="poisson", 
                 knots = list(week = c(0, 52)))
bird_move <- transform(bird_move, mod2 = predict(bird_mod2, type="response"))

bird_mod2_indiv <- ggplot(data=bird_move, 
                          aes(x=week, y=latitude, fill=mod2,color=mod2)) +
  geom_tile(size=0.25) +
  facet_wrap(~ species, ncol=6) +
  scale_fill_viridis("Count") +
  scale_color_viridis("Count") +
  scale_x_continuous(expand=c(0, 0), breaks=c(1, 26, 52)) +
  scale_y_continuous(expand=c(0, 0), breaks=c(0, 30, 60)) +
  labs(x = "Week", y = "Latitude") +
  theme(legend.position="right")

bird_mod2_indiv_fit <- ggplot(data=bird_move, aes(x=mod2, y=count)) +
  facet_wrap(~ species, ncol=6) +
  geom_point(alpha=0.1) +
  geom_abline() +
  labs(x="Predicted count (model 2)", y= "Observed count")

plot_grid(bird_mod2_indiv, bird_mod2_indiv_fit, ncol=1, align="vh", axis = "lrtb",
          labels=c("a","b"), rel_heights= c(1,1))
## CO2_mod3 <- gam(log(uptake) ~ s(log(conc), k=5, m=2, bs="tp") +
##                   s(log(conc), by=Plant_uo, k=5, m=1, bs="tp") +
##                   s(Plant_uo, bs="re", k=12),
##                 data=CO2, method="REML")
CO2_mod3 <- gam(log(uptake) ~ s(log(conc), k=5, m=2, bs="tp") +
                  s(log(conc), by= Plant_uo, k=5, m=1, bs="tp") +
                  s(Plant_uo, bs="re", k=12),
                data=CO2, method="REML")

op <- par(mfrow=c(2, 3), mar =c(4, 4, 1, 1))
CO_mod3_ylim = c(-0.3,0.3)
plot(CO2_mod3, scale=0, 
     select=1, 
     ylab="Global smoother", 
     seWithMean=TRUE)
plot(CO2_mod3, 
     scale=0, 
     select=14, 
     ylab="Intercept",
     main=NA)
plot(CO2_mod3, 
     scale=0, 
     select=3,  
     ylab="Plant Qn1", 
     seWithMean=TRUE,
     ylim = CO_mod3_ylim)
plot(CO2_mod3, 
     scale=0, 
     select=5,  
     ylab="Plant Qc1",     
     seWithMean=TRUE,
     ylim = CO_mod3_ylim)
plot(CO2_mod3, 
     scale=0, 
     select=10, 
     ylab="Plant Mn1",     
     seWithMean=TRUE,
     ylim = CO_mod3_ylim)
plot(CO2_mod3, 
     scale=0, 
     select=13, 
     ylab="Plant Mc1",
     seWithMean=TRUE,
     ylim = CO_mod3_ylim)
par(op)
bird_mod3 <- gam(count ~ species +
                   te(week, latitude, bs=c("cc", "tp"),
                      k=c(10, 10), m=2) +
                   te(week, latitude, by=species, bs= c("cc", "tp"),
                      k=c(10, 10), m=1),
                 data=bird_move, method="REML", family="poisson",
                 knots = list(week = c(0, 52)))
CO2_mod4 <- gam(log(uptake) ~ s(log(conc), Plant_uo, k=5, bs="fs", m=2),
                data=CO2, method="REML")

bird_mod4 <- gam(count ~ t2(week, latitude, species, bs=c("cc", "tp", "re"),
                            k=c(10, 10, 6), m=c(2, 2, 2)),
                 data=bird_move, method="REML", family="poisson",
                 knots = list(week = c(0, 52)))
CO2_mod5 <- gam(log(uptake) ~ s(log(conc), by=Plant_uo, k=5, bs="tp", m=2) +
                              s(Plant_uo, bs="re", k=12),
                data= CO2, method="REML")


bird_mod5 <- gam(count ~ species + 
                   te(week, latitude, by=species, bs= c("cc", "tp"), 
                      k=c(10, 10), m=2),
                 data=bird_move, method="REML", family="poisson",
                 knots = list(week = c(0, 52)))

AIC_table = AIC(CO2_mod1,CO2_mod2, CO2_mod3, CO2_mod4, CO2_mod5,
             bird_mod1, bird_mod2, bird_mod3, bird_mod4, bird_mod5)%>%
  rownames_to_column(var= "Model")%>%
  mutate(data_source = rep(c("CO2","bird_data"), each =5))%>%
  group_by(data_source)%>%
  mutate(deltaAIC = AIC - min(AIC))%>%
  ungroup()%>%
  dplyr::select(-data_source)%>%
  mutate_at(.vars = vars(df,AIC, deltaAIC), 
            .funs = funs(round,.args = list(digits=0)))

#### Code for IV: Examples ####
zooplankton <- read.csv("../data/zooplankton_example.csv")%>%
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

#This function calculates the root-mean-squared-error for out-of-sample data
get_RMSE <- function(fit, obs) sqrt(mean((fit-obs)^2))

zoo_comm_mod4 <- gam(density_adj ~ s(day, taxon,
                                     bs="fs",
                                     k=10,
                                     xt=list(bs="cc"))+
                                   s(taxon, year_f, bs="re"),
                     data=zoo_train,
                     knots = list(day =c(0, 365)),
                     family = Gamma(link ="log"), 
                     method = "REML",
                     drop.unused.levels = FALSE)
# Note that  s(taxon, bs="re") has to be explicitly included here, as the 
# day  by taxon smoother does not include an intercept
zoo_comm_mod5 <- gam(density_adj ~ s(day, by=taxon,
                                     k=10, bs="cc") + 
                                   s(taxon, bs="re") +
                                   s(taxon, year_f, bs="re"),
                     data=zoo_train,
                     knots = list(day =c(0, 365)),
                     family = Gamma(link ="log"), 
                     method = "REML",
                     drop.unused.levels = FALSE)
## gam.check(zoo_comm_mod4)
## gam.check(zoo_comm_mod5)
## #individual components of gam.check: the results for k.check
## round(k.check(zoo_comm_mod5),2)
#individual components of gam.check: residual plots
par(mfrow= c(1,2))
qq.gam(zoo_comm_mod4)
plot(log(fitted(zoo_comm_mod5)), 
     residuals.gam(zoo_comm_mod5,type = "deviance"), 
     xlab = "linear predictor",
     ylab = "residuals")
#Create synthetic data to use to compare predictions
zoo_plot_data <- expand.grid(day = 1:365, 
                             taxon = factor(levels(zoo_train$taxon)), 
                             year_f = 1980)

#extract predicted values and standard errors for both models. the exclude = "s(taxon,year_f)" 
#term indicates that predictions should be made excluding the effect of the
#taxon by year random effect (effectively setting making predictions averaging
#over year-taxon effects).
zoo_mod4_fit <- predict(zoo_comm_mod4, 
                        zoo_plot_data, 
                        se.fit = TRUE, 
                        exclude = "s(taxon,year_f)")
zoo_mod5_fit <- predict(zoo_comm_mod5, 
                        zoo_plot_data, 
                        se.fit = TRUE, 
                        exclude = "s(taxon,year_f)")

zoo_plot_data$mod4_fit <- as.numeric(zoo_mod4_fit$fit)
zoo_plot_data$mod5_fit <- as.numeric(zoo_mod5_fit$fit)

zoo_plot_data <- gather(zoo_plot_data, model, fit, mod4_fit, mod5_fit)
zoo_plot_data <- mutate(zoo_plot_data, se= c(as.numeric(zoo_mod4_fit$se.fit),
                                             as.numeric(zoo_mod5_fit$se.fit)),
                        upper = exp(fit + (2 * se)),
                        lower = exp(fit - (2 * se)),
                        fit   = exp(fit))

#Plot the model output, with means plus standard deviations for each model.
zoo_plot <- ggplot(zoo_plot_data) +
  facet_wrap(~taxon, nrow = 4,scales = "free_y")+
  geom_ribbon(aes(x=day,
                  ymin = lower,
                  ymax = upper,
                  fill = model),
              alpha=0.2)+
  geom_point(data= zoo_train, aes(x = day, y = density_adj),size=0.1)+
  geom_line(aes(x = day, y = fit, color = model))+
  labs(y = expression(atop(Population~density,("10 000"~individuals~m^{-2}))), 
       x = "Day of Year") +
    scale_fill_brewer(name = "", palette = "Dark2",
                      labels = paste("Model", 4:5)) +
    scale_colour_brewer(name = "",
                        palette = "Dark2", labels = paste("Model", 4:5))+
  theme(legend.position = "top")

zoo_plot

#Getting the out of sample predictions for both models:

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
    mod4 = predict(zoo_comm_mod4, ., type="response"),
    mod5 = predict(zoo_comm_mod5, ., type="response"))%>%
  group_by(taxon)%>%
  summarise(
    `Intercept only` = format(get_RMSE(mod0, density_adj), 
                              scientific = FALSE, 
                              digits=3),
    `Model 4` = format(get_RMSE(mod4, density_adj), 
                       scientific = FALSE, 
                       digits=3),
    `Model 5` = format(get_RMSE(mod5, density_adj), 
                       scientific = FALSE, 
                       digits=3))%>%
  #need to specify this to ensure that species names are italized in the table
  mutate(taxon = cell_spec(taxon, 
                           italic = c(TRUE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE))) 

zoo_daph_mod1 <- gam(density_adj ~ s(day, bs="cc", k=10)+
                       s(lake, bs="re") + 
                       s(lake, year_f,bs="re"),
                     data=daphnia_train,
                     knots=list(day =c(0, 365)),
                     family=Gamma(link ="log"),
                     method="REML",
                     drop.unused.levels = FALSE)
zoo_daph_mod2 <-
  gam(density_adj ~ s(day, bs="cc", k=10) + 
        s(day, lake, k=10, bs="fs", xt=list(bs="cc")) + 
        s(lake, year_f,bs="re"),
      data=daphnia_train, 
      knots=list(day=c(0, 365)), 
      family=Gamma(link ="log"),
      drop.unused.levels = FALSE,
      method="REML")
zoo_daph_mod3 <- gam(density_adj~s(day, bs="cc", k=10) +
                             s(day, by=lake, k=10, bs="cc")+
                             s(lake, bs="re") + 
                             s(lake, year_f,bs="re"),
                     data=daphnia_train,
                     knots=list(day =c(0, 365)),
                     family=Gamma(link ="log"),
                     method="REML",
                     drop.unused.levels = FALSE)
#Create synthetic data to use to compare predictions
daph_plot_data <- expand.grid(day = 1:365, 
                              lake = factor(levels(zoo_train$lake)),
                              year_f = 1980)

#extract predicted values and standard errors for both models. the exclude = "s(taxon,year_f)" 
#term indicates that predictions should be made excluding the effect of the
#taxon by year random effect (effectively setting making predictions averaging
#over year-taxon effects).
daph_mod1_fit <- predict(zoo_daph_mod1, 
                         newdata = daph_plot_data, 
                         se.fit = TRUE, 
                         exclude = "s(lake,year_f)")
daph_mod2_fit <- predict(zoo_daph_mod2, 
                         newdata = daph_plot_data, 
                         se.fit = TRUE, 
                         exclude = "s(lake,year_f)")
daph_mod3_fit <- predict(zoo_daph_mod3, 
                         newdata = daph_plot_data, 
                         se.fit = TRUE, 
                         exclude = "s(lake,year_f)")

daph_plot_data$mod1_fit <- as.numeric(daph_mod1_fit$fit)
daph_plot_data$mod2_fit <- as.numeric(daph_mod2_fit$fit)
daph_plot_data$mod3_fit <- as.numeric(daph_mod3_fit$fit)

daph_plot_data <- gather(daph_plot_data, 
                         key = model, 
                         value = fit, 
                         mod1_fit, 
                         mod2_fit, 
                         mod3_fit)

daph_plot_data <- mutate(daph_plot_data, 
                         se = c(as.numeric(daph_mod1_fit$se.fit),
                                as.numeric(daph_mod2_fit$se.fit),
                                as.numeric(daph_mod3_fit$se.fit)),
                         upper = exp(fit + (2 * se)),
                         lower = exp(fit - (2 * se)),
                         fit   = exp(fit))


daph_plot <- ggplot(daph_plot_data, aes(x=day))+
  facet_wrap(~lake, nrow = 2)+
  geom_ribbon(aes(x = day, ymin = lower, ymax = upper, fill = model), 
                alpha = 0.2) +
  geom_point(data= daphnia_train, aes(x = day, y = density_adj),size=0.1)+
  geom_line(aes(x = day, y = fit, colour = model)) +

  labs(y = expression(atop(Population~density,
                           ("10 000"~individuals~m^{-2}))), 
       x = "Day of Year") +
  scale_x_continuous(expand = c(0,0))+
    scale_fill_brewer(name = "", palette = "Dark2",
                      labels = paste("Model", 1:3)) +
    scale_colour_brewer(name = "",
                        palette = "Dark2", labels = paste("Model", 1:3))


daph_plot
# we need to compare how well this model fits with a null model. here we'll use an
# intercept-only model
zoo_daph_mod0 <- gam(density_adj~s(lake, bs="re"),
                     data=daphnia_train,
                     knots=list(day =c(0, 365)),
                     family=Gamma(link ="log"),
                     method="REML",
                     drop.unused.levels = FALSE)



# We'll look at the correlation between fitted and observed values for all species:

daph_test_summary <- daphnia_test %>%
  mutate(#get out-of-sample predicted fits
    mod0 = as.numeric(predict(zoo_daph_mod0,.,type="response")),
    mod1 = as.numeric(predict(zoo_daph_mod1,.,type="response")),
    mod2 = as.numeric(predict(zoo_daph_mod2,.,type="response")),
    mod3 = as.numeric(predict(zoo_daph_mod3,.,type="response")))%>%
  group_by(lake)%>%
  summarise(`Intercept only` = format(get_RMSE(mod0, density_adj), 
                                      scientific = FALSE, digits=2),
            `Model 1` = format(get_RMSE(mod1, density_adj), 
                               scientific = FALSE, digits=2),
            `Model 2` = format(get_RMSE(mod2, density_adj), 
                               scientific = FALSE, digits=2),
            `Model 3` = format(get_RMSE(mod3, density_adj), 
                               scientific = FALSE, digits=2))%>%
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
freq_vals = c(1/4,1/2,1,2,4)
dat = crossing(x = seq(0,2*pi,length=150),freq = freq_vals)%>%
  mutate(y = sin(freq*x) +rnorm(n(), 0, 0.2),
         grp = paste("frequency = ",freq,sep= ""),
         grp = factor(grp,  levels = paste("frequency = ",freq_vals,sep= "")))

#Fit a type-4 function (shared smoothness) for the test data
mod1 = bam(y~s(x,k=30,grp, bs="fs"), data=dat)

#Fit a type-5 function (differing smoothness) for the test data
mod2 = bam(y~s(x,k=30,by=grp)+s(grp,bs="re"), data=dat)

#Extract fitted values for each model for all the test data
overfit_predict_data = crossing(x = seq(0,2*pi,length=500), 
                                freq = freq_vals)%>%
  mutate(grp = paste("frequency = ",freq,sep= ""),
         grp = factor(grp,  levels = paste("frequency = ",freq_vals,sep= "")),
         y = sin(freq*x))%>%
  mutate(fit1 = as.numeric(predict(mod1,newdata = .,type = "response")),
         fit2 = as.numeric(predict(mod2,newdata = .,type="response")))

#turn this into long-format data for plotting, and to make it easier to
#calculate derivatives
overfit_predict_data_long = overfit_predict_data %>%
  gather(model, value, y, fit1, fit2)%>%
  mutate(model = recode(model, y = "true value",fit1 = "model 4 fit",
                        fit2 = "model 5 fit"),
         model = factor(model, levels=  c("true value","model 4 fit",
                                          "model 5 fit")))

#estimate 2nd derivatives of each curve, then for each curve calculate the sum
#of squared 2nd derivatives of the true curve and the predictions for both
#models.
deriv_est_data = overfit_predict_data%>%
  group_by(grp)%>%
  arrange(grp, x)%>%
  mutate(fit1_deriv = calc_2nd_deriv(x,fit1),
         fit2_deriv = calc_2nd_deriv(x,fit2))%>%
  summarize(freq = freq[1], fit1_int = sum(fit1_deriv^2*(x-lag(x)),
                                           na.rm = TRUE),
            fit2_int = sum(fit2_deriv^2*(x-lag(x)),na.rm = TRUE))%>%
  ungroup()%>%
  mutate(sqr_2nd_deriv = -freq^3*(sin(4*pi*freq)-4*pi*freq)/4)%>%
  gather(key=model,value = obs_sqr_deriv,fit1_int,fit2_int)%>%
  mutate(model = factor(ifelse(model=="fit1_int", "model 4 fit",
                               "model 5 fit"),
                        levels = c("model 4 fit",
                                   "model 5 fit")))

#The derivative plots
deriv_plot =  ggplot(data=deriv_est_data, aes(x=sqr_2nd_deriv, 
                                              y= obs_sqr_deriv,
                                              color= model))+
  geom_point()+
  scale_y_log10("Estimated wiggliness\nof fitted curves")+
  scale_x_log10("Wiggliness of true curve")+
  scale_color_brewer(name=NULL,palette= "Set1")+
  geom_abline(color="black")+
  theme(legend.position = "top")

fit_colors = c("black",RColorBrewer::brewer.pal(3, "Set1")[1:2])

overfit_vis_plot = ggplot(data=overfit_predict_data_long,
                          aes(x=x,y= value,color=model))+
  geom_line()+
  scale_color_manual(values=fit_colors)+
  facet_grid(.~grp)+
  theme(legend.position = "top")
#Plot the overfit graphs together.
cowplot::plot_grid(overfit_vis_plot, deriv_plot, ncol=1, labels="auto",
                   align="hv", axis="lrtb")
#Note: this code takes quite a long time to run! It's fitting all 10 models. Run
#once if possible, then rely on the cached code. There's a reason it's split off
#from the rest of the chunks of code.

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
comp_resources = crossing(model_number = c("1","2","3","4","5"),
                          data_source = factor(c("CO2","bird_move"),
                                               levels = c("CO2","bird_move")),
                          time = 0, n_smooths = 0,
                          n_coef = 0)


#Fit each model to the example data sets, and calculate run time for them
comp_resources[1,"time"] = system.time(
  CO2_mod1 <- gam(log(uptake) ~ s(log(conc),k=5,m=2, bs="tp")+
                                s(Plant_uo, k =12,  bs="re"),
                  data= CO2,
                  method="REML",
                  control = list(keepData=TRUE))
  )[3]

comp_resources[2,"time"] = system.time(
  bird_mod1 <- gam(count ~ te(week,latitude, bs= c("cc", "tp"), k=c(10,10)),
                   data= bird_move, 
                   method="REML", 
                   family= poisson,
                   control = list(keepData=TRUE))
  )[3]

comp_resources[3,"time"] = system.time(
  CO2_mod2 <- gam(log(uptake) ~ s(log(conc),k=5,m=2, bs="tp")+
                                s(log(conc), Plant_uo, k=5, bs="fs",m=1),
                  data= CO2,
                  method="REML",
                  control = list(keepData=TRUE))
  )[3]


comp_resources[4,"time"] = system.time(
  bird_mod2 <- gam(count ~ te(week,latitude, bs= c("cc", "tp"),
                              k=c(10,10),m=c(2,2))+
                           te(week,latitude,species, bs= c("cc", "tp","re"),
                              k=c(10,10,6),m = c(1,1,1)),
                   data= bird_move, 
                   method="REML", 
                   family= poisson,
                   control = list(keepData=TRUE))
  )[3]


comp_resources[5,"time"] = system.time(
  CO2_mod3 <- gam(log(uptake) ~ s(log(conc),k=5,m=2, bs="tp")+
                                s(log(conc),by= Plant_uo, k =5,  bs="ts",m=1)+
                                s(Plant_uo,bs="re",k=12),
                  data= CO2,
                  method="REML",
                  control = list(keepData=TRUE)))[3]



comp_resources[6,"time"] = system.time(
  bird_mod3 <- gam(count ~ te(week,latitude, bs= c("cc", "tp"),
                              k=c(10,10),m=c(2,2)) +
                           te(week,latitude, bs= c("cc", "tp"),
                              k=c(10,10),m=c(1,1),by= species),
                   data= bird_move, 
                   method="REML", 
                   family= poisson,
                   control = list(keepData=TRUE)))[3]


comp_resources[7,"time"] = system.time(
  CO2_mod4 <- gam(log(uptake) ~ s(log(conc), Plant_uo, k=5,  bs="fs",m=2),
                  data= CO2,
                  method="REML",
                  control = list(keepData=TRUE))
  )[3]


comp_resources[8,"time"] = system.time(
  bird_mod4 <- gam(count ~ te(week,latitude,species, bs= c("cc", "tp","re"),
                              k=c(10,10,6),m = 2),
                   data= bird_move, 
                   method="REML", 
                   family= poisson,
                   control = list(keepData=TRUE))
)[3]


comp_resources[9,"time"] = system.time(
  CO2_mod5 <- gam(log(uptake) ~ s(log(conc),by= Plant_uo, k =5,  bs="tp",m=2) +
                                s(Plant_uo,bs="re",k=12),
                  data= CO2,
                  method="REML",
                  control = list(keepData=TRUE))
)[3]

comp_resources[10,"time"] = system.time(
  bird_mod5 <- gam(count ~ te(week,latitude,by=species, bs= c("cc", "tp"),
                              k=c(10,10),m = 2),
                   data= bird_move, 
                   method="REML", 
                   family= poisson,
                   control = list(keepData=TRUE))
)[3]

#combine all fitted models into a list
comp_resources$model = list(CO2_mod1, bird_mod1, CO2_mod2, bird_mod2,
                            CO2_mod3, bird_mod3,CO2_mod4, bird_mod4,
                            CO2_mod5, bird_mod5)

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
  transmute(data_source =data_source, model=model_number,
            `relative time` = time,`coefficients` = n_coef,
            `penalties` = n_smooths
            )%>%
  group_by(data_source) %>%
  mutate(#scales processing time relative to model 1
         `relative time` = `relative time`/`relative time`[1],
         #rounds to illustrate differences in timing.
         `relative time` = ifelse(`relative time`<10, 
                                  signif(`relative time`,1), 
                                  signif(`relative time`, 2)) 
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
           indiv_func = int + x^2*x2 + 2*(exp(x*logit_slope)/(1+exp(x*logit_slope))-0.5),
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
  gather(model, timing, gam,`bam (discrete = FALSE)`, 
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
  scale_color_manual(values = c("black", "#1b9e77","#1b9e77", "#d95f02", "#7570b3"))+
  scale_linetype_manual(values =c(1,1,2,1,1))+
  scale_y_log10("run time (seconds)", 
                breaks = c(0.1,1,10,100), 
                labels = c("0.1", "1","10", "100"))+
  scale_x_log10("number of groups", 
                breaks = c(2,8,32,128))+
  guides(color = guide_legend(nrow = 2, byrow = TRUE))+
  theme(legend.position = "top")
timing_plot
