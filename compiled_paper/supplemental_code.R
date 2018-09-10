knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE, cache=TRUE)
library(mgcv)
library(gamm4)
library(tidyr)
library(ggplot2)
library(viridis)
library(cowplot)
library(kableExtra)
library(knitr)
library(tibble)
library(dplyr)

#Set the default theme for ggplot objects to theme_bw()
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

table_out_format <- ifelse("pdf_document" %in% rmarkdown::all_output_formats("full_document.Rmd"),
                    "latex",
                     ifelse("html_document" %in% rmarkdown::all_output_formats("full_document.Rmd"),
                            "html",
                            NA_character_)
                      )




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

pdat <- with(dat, data.frame(x = seq(min(x), max(x), length = 200)))
p <- cbind(pdat, fit = predict(b, newdata = pdat))
p.0 <- cbind(pdat, fit = predict(b.0, newdata = pdat))
p.inf <- cbind(pdat, fit = predict(b.inf, newdata = pdat))
ylims <- range(p, p.0, p.inf)

lab.l <- labs(x = "x", y = "y")
dat.l <- geom_point(data = dat, aes(x = x, y = y), colour = "darkgrey")
true.l <- geom_line(data = true, aes(x = x, y = y), colour = "blue")
coord.l <- coord_cartesian(ylim = ylims)

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
                data=CO2, method="REML", family = "gaussian")

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
  geom_point() +
  geom_line(aes(y=exp(fit)), data=CO2_mod1_pred) +
  geom_ribbon(aes(ymin=exp(fit - 2*se.fit), ymax=exp(fit + 2*se.fit), x=conc),
              data=CO2_mod1_pred, alpha=0.3, inherit.aes=FALSE) +
  labs(x=expression(CO[2] ~ concentration ~ (mL ~ L^{-1})),
       y=expression(CO[2] ~ uptake ~ (mu*mol ~ m^{-2})))
## # Note for specifying tensor products: you can either specify bs (basis) and
## # k (number of basis functions) as single values, which would assign the same
## # basis and k to each marginal value, or pass them as vectors, one value for each
## # distinct marginal smooth (see ?mgcv::te for details)
bird_mod1 <- gam(count ~ te(week, latitude, bs=c("cc", "tp"), k=c(10, 10)),
                 data=bird_move, method="REML", family=poisson,
                 knots = list(week = c(0.5, 52.5)))
#mgcv gam plot for the two-dimensional tensor product smoother for bird_mod1.
#scheme=2 displays the color scheme (rather than mgcv's default, which only
#shows contour lines)
plot(bird_mod1, pages=1, scheme=2, rug=FALSE)
box()
bird_move <- transform(bird_move, mod1 = predict(bird_mod1, type="response"))

ggplot(bird_move, aes(x=mod1, y=count)) +
  facet_wrap(~species) +
  geom_point() +
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
  geom_point() +
  geom_line(aes(y=exp(mod2))) +
  geom_ribbon(aes(ymin=exp(mod2-2*mod2_se),
                  ymax=exp(mod2+2*mod2_se)), alpha=0.25) +
  labs(x=expression(CO[2] ~ concentration ~ (mL ~ L^{-1})),
       y=expression(CO[2] ~ uptake ~ (mu*mol ~ m^{-2})))
bird_mod2 <- gam(count ~ te(week, latitude, bs=c("cc", "tp"),
                            k=c(10, 10), m=c(2, 2)) +
                         t2(week, latitude, species, bs=c("cc", "tp", "re"),
                            k=c(10, 10, 6), m=c(2, 2, 2), full=TRUE),
                 data=bird_move, method="REML", family=poisson)
bird_move <- transform(bird_move, mod2 = predict(bird_mod2, type="response"))

bird_mod2_indiv <- ggplot(data=bird_move, aes(x=week, y=latitude, fill=mod2,color=mod2)) +
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
  geom_point() +
  geom_abline() +
  labs(x="Predicted count (model 2)", y= "Observed count")

plot_grid(bird_mod2_indiv, bird_mod2_indiv_fit, ncol=1, align="vh", axis = "lrtb",
          labels=c("a","b"), rel_heights= c(1,1))
## CO2_mod3 <- gam(log(uptake) ~ s(log(conc), k=5, m=2, bs="tp") +
##                               s(log(conc), by=Plant_uo, k=5, m=1, bs="ts") +
##                               s(Plant_uo, bs="re", k=12),
##                 data=CO2, method="REML")
CO2_mod3 <- gam(log(uptake) ~ s(log(conc), k=5, m=2, bs="tp") +
                              s(log(conc), by= Plant_uo, k=5, m=1, bs="ts") +
                              s(Plant_uo, bs="re", k=12),
                data=CO2, method="REML")

op <- par(mfrow=c(2, 3), mar =c(4, 4, 1, 1))
plot(CO2_mod3, scale=0, select=1,  ylab="Global smooth", seWithMean=TRUE)
plot(CO2_mod3, scale=0, select=14, ylab="Intercept",     main=NA)
plot(CO2_mod3, scale=0, select=3,  ylab="Plant Qn1",     seWithMean=TRUE)
plot(CO2_mod3, scale=0, select=5,  ylab="Plant Qc1",     seWithMean=TRUE)
plot(CO2_mod3, scale=0, select=10, ylab="Plant Mn1",     seWithMean=TRUE)
plot(CO2_mod3, scale=0, select=13, ylab="Plant Mc1",     seWithMean=TRUE)
par(op)
bird_mod3 <- gam(count ~ te(week, latitude, bs=c("cc", "tp"),
                            k=c(10, 10), m=c(2, 2)) +
                         te(week, latitude, by=species, bs= c("cc", "ts"),
                            k=c(10, 10), m=c(1, 1)),
                 data=bird_move, method="REML", family=poisson)
CO2_mod4 <- gam(log(uptake) ~ s(log(conc), Plant_uo, k=5,  bs="fs", m=2),
                data=CO2, method="REML")

bird_mod4 <- gam(count ~ t2(week, latitude, species, bs=c("cc", "tp", "re"),
                            k=c(10, 10, 6), m = c(2,2,2)),
                 data=bird_move, method="REML", family=poisson)
CO2_mod5 <- gam(log(uptake) ~ s(log(conc), by=Plant_uo, k=5, bs="tp", m=2) +
                              s(Plant_uo, bs="re", k=12),
                data= CO2, method="REML")


bird_mod5 <- gam(count ~ te(week, latitude, by=species, bs= c("cc", "ts"),
                            k=c(10, 10), m=c(2,2)),
                 data=bird_move, method="REML", family=poisson)

AIC_table = AIC(CO2_mod1,CO2_mod2, CO2_mod3, CO2_mod4, CO2_mod5,
             bird_mod1, bird_mod2, bird_mod3, bird_mod4, bird_mod5)%>%
  rownames_to_column(var= "Model")%>%
  mutate_at(.vars = vars(df,AIC), .funs = funs(round,.args = list(digits=0)))


kable(AIC_table, format =table_out_format, caption="AIC table comparing model fits for example datasets", booktabs = T)%>% 
  kable_styling(full_width = F)%>%
  group_rows("A. CO2 models", 1,5)%>%
  group_rows("B. bird_move models", 6,10)


zooplankton <- read.csv("../data/zooplankton_example.csv")

#This is what the data looks like:
str(zooplankton)
levels(zooplankton$taxon)
levels(zooplankton$lake)
zoo_train <- subset(zooplankton, year%%2==0)
zoo_test <- subset(zooplankton, year%%2==1) 
zoo_comm_mod4 <-
  gam(density_scaled ~ s(day, taxon, bs="fs", k=10, xt=list(bs="cc")),
      data=zoo_train, knots=list(day=c(1, 365)), method="ML")
zoo_comm_mod5 <-
  gam(density_scaled ~ taxon + s(day, by=taxon, k=10, bs="cc"),
      data=zoo_train, knots=list(day=c(1, 365)), method="ML")
#Create synthetic data to use to compare predictions
zoo_plot_data <- expand.grid(day = 1:365, taxon = factor(levels(zoo_train$taxon)))

#extract predicted values and standard errors for both models
zoo_mod4_fit <- predict(zoo_comm_mod4, zoo_plot_data, se.fit = T)
zoo_mod5_fit <- predict(zoo_comm_mod5, zoo_plot_data, se.fit = T)

zoo_plot_data$mod4_fit <- as.numeric(zoo_mod4_fit$fit)
zoo_plot_data$mod5_fit <- as.numeric(zoo_mod5_fit$fit)

zoo_plot_data$mod4_se <- as.numeric(zoo_mod4_fit$se.fit)
zoo_plot_data$mod5_se <- as.numeric(zoo_mod5_fit$se.fit)

#Plot the model output, with means plus standard deviations for each model.
zoo_plot = ggplot(zoo_plot_data, aes(x=day))+
  facet_wrap(~taxon, nrow = 2)+
  geom_point(data= zoo_train, aes(y=density_scaled),size=0.1)+
  geom_line(aes(y=mod4_fit, color = "Model 4"))+
  geom_line(aes(y=mod5_fit, color = "Model 5"))+
  geom_ribbon(aes(ymin = mod4_fit - 2*mod4_se,
                  ymax = mod4_fit + 2*mod4_se,
                  fill="Model 4"),
              alpha=0.25)+
  geom_ribbon(aes(ymin = mod5_fit - 2*mod5_se,
                  ymax = mod5_fit + 2*mod5_se,
                  fill="Model 5"),
              alpha=0.25)+
  scale_y_continuous("Scaled log-transformed density")+
  scale_x_continuous(expand = c(0,0))+
  scale_color_manual("", breaks = c("Model 4", "Model 5"), values = c("black","red"))+
  scale_fill_manual("", breaks = c("Model 4", "Model 5"), values = c("black","red"))+
  theme(legend.position = "bottom")

zoo_plot
get_MSE = function(obs, pred) mean((obs-pred)^2)
#Getting the out of sample predictions for both models:
zoo_test$mod4 = as.numeric(predict(zoo_comm_mod4,zoo_test))
zoo_test$mod5 = as.numeric(predict(zoo_comm_mod5,zoo_test))

#Correlations between fitted and observed values for all species:
#\n is in variable titles to add a line break in the printed table.
zoo_test_summary = zoo_test %>%
  group_by(taxon)%>%
  summarise(`model 4 MSE` = round(get_MSE(density_scaled,mod4),2),
            `model 5 MSE` = round(get_MSE(density_scaled,mod5),2))

kable(zoo_test_summary, format = table_out_format, caption="Out-of-sample predictive ability for model 4 and 5 applied to the zooplankton community dataset. MSE values represent the average squared difference between model predictions and observations for test data.", booktabs = TRUE)%>%
  kable_styling(full_width = FALSE)
daphnia_train <- subset(zoo_train, taxon=="D. mendotae")
daphnia_test <- subset(zoo_test, taxon=="D. mendotae")

zoo_daph_mod1 <-
  gam(density_scaled ~ s(day, bs="cc", k=10),
      data=daphnia_train, knots=list(day=c(1, 365)), method="ML")

printCoefmat(summary(zoo_daph_mod1)$s.table)
zoo_daph_mod2 <-
  gam(density_scaled ~ s(day, bs="cc", k=10) +
        s(day, lake, k=10, bs="fs", xt=list(bs="cc")),
      data=daphnia_train, knots=list(day=c(1, 365)), method="ML")

printCoefmat(summary(zoo_daph_mod2)$s.table)
zoo_daph_mod3 <-
  gam(density_scaled ~ lake + s(day, bs="cc", k=10) + 
        s(day, by=lake, k=10, bs="cc"),
      data=daphnia_train, knots=list(day=c(1, 365)), method="ML")

printCoefmat(summary(zoo_daph_mod3)$s.table)
#Create synthetic data to use to compare predictions
daph_plot_data = expand.grid(day = 1:365, lake = factor(levels(zoo_train$lake)))


daph_mod1_fit = predict(zoo_daph_mod1, daph_plot_data, se.fit = T)
daph_mod2_fit = predict(zoo_daph_mod2, daph_plot_data, se.fit = T)
daph_mod3_fit = predict(zoo_daph_mod3, daph_plot_data, se.fit = T)


daph_plot_data$mod1_fit = as.numeric(daph_mod1_fit$fit)
daph_plot_data$mod2_fit = as.numeric(daph_mod2_fit$fit)
daph_plot_data$mod3_fit = as.numeric(daph_mod3_fit$fit)

daph_plot_data$mod1_se = as.numeric(daph_mod1_fit$se.fit)
daph_plot_data$mod2_se = as.numeric(daph_mod2_fit$se.fit)
daph_plot_data$mod3_se = as.numeric(daph_mod3_fit$se.fit)

daph_plot = ggplot(daph_plot_data, aes(x=day))+
  facet_wrap(~lake, nrow = 2)+
  geom_point(data= daphnia_train, aes(y=density_scaled),size=0.1)+
  geom_line(aes(y=mod1_fit, linetype = "Model 1", size = "Model 1"))+
  geom_line(aes(y=mod2_fit,color = "Model 2"))+
  geom_line(aes(y=mod3_fit,color = "Model 3"))+
  geom_ribbon(aes(ymin = mod2_fit - 2*mod2_se,
                  ymax = mod2_fit + 2*mod2_se,fill = "Model 2"),
              alpha=0.25)+
  geom_ribbon(aes(ymin = mod3_fit - 2*mod3_se,
                  ymax = mod3_fit + 2*mod3_se,fill = "Model 3"),
              alpha=0.25)+
  scale_y_continuous("Scaled log-transformed density")+
  scale_x_continuous(expand = c(0,0))+
  scale_color_manual("", breaks = c("Model 2", "Model 3"), values = c("black","red"))+
  scale_fill_manual("", breaks = c("Model 2", "Model 3"), values = c("black","red"))+
  scale_linetype_manual("", breaks = c("Model 1"), values = 2)+
  scale_size_manual("", breaks = c("Model 1"), values = 2)+
  theme(legend.position = "bottom")

print(daph_plot)
#Getting the out of sample predictions for both models:
daphnia_test$mod1 = as.numeric(predict(zoo_daph_mod1,daphnia_test))
daphnia_test$mod2 = as.numeric(predict(zoo_daph_mod2,daphnia_test))
daphnia_test$mod3 = as.numeric(predict(zoo_daph_mod3,daphnia_test))

# We'll look at the correlation between fitted and observed values for all species:
daph_test_summary = daphnia_test %>%
  group_by(lake)%>%
  summarise(`model 1 MSE` = round(get_MSE(density_scaled,mod1),2),
            `model 2 MSE` = round(get_MSE(density_scaled,mod2),2),
            `model 3 MSE` = round(get_MSE(density_scaled,mod3),2))

kable(daph_test_summary,format = table_out_format, caption="Out-of-sample predictive ability for model 1-3 applied to the \\textit{D. mendotae} dataset. MSE values represent the average squared difference between model predictions and observations for held-out data (zero predictive ability would correspond to a MSE of one).", booktabs = T)%>%
  kable_styling(full_width = F)

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

print(cowplot::plot_grid(overfit_vis_plot,
                         deriv_plot,
                         ncol = 1,
                         labels = "auto",
                         label_y = c(1,1.2)))
#Note: this code takes quite a long time to run! It's fitting all 10 models.
#Run once if possible, then rely on the cached code. There's a reason it's split off from the rest of the chunks of code.
source("../code/functions.R")

get_n_pen  = function(model) {
  family = model$family[[1]]
  if(family %in% c("Gamma","gaussian")){
    capture.output({out_val = nrow(gam.vcomp(model))-1})
  }else{
    capture.output({out_val = nrow(gam.vcomp(model))})
  }
  return(out_val)
}

get_n_coef = function(model) length(coef(model))

get_n_iter = function(model) model$outer.info$iter
get_n_out_iter = function(model) model$iter

bird_move = read.csv("../data/bird_move.csv")

comp_resources = crossing(model_number = c("1","2","3","4","5"),
                       data_source = factor(c("CO2","bird_move"),
                                     levels = c("CO2","bird_move")),
                       time = 0, n_smooths = 0,
                       n_coef = 0)

CO2$Plant_uo = factor(CO2$Plant, levels = levels(CO2$Plant), ordered = F)


comp_resources[1,"time"] = system.time(CO2_mod1 <- gam(log(uptake) ~ s(log(conc),k=5,m=2, bs="tp")+s(Plant_uo, k =12,  bs="re"),
                                                    data= CO2,method="REML",
                                                    control = list(keepData=TRUE)))[3]

comp_resources[2,"time"] = system.time(bird_mod1 <- gam(count ~ te(week,latitude, bs= c("cc", "tp"), k=c(10,10)),
                                                     data= bird_move, method="REML", family= poisson,
                                                     control = list(keepData=TRUE)))[3]

comp_resources[3,"time"] = system.time(CO2_mod2 <- gam(log(uptake) ~ s(log(conc),k=5,m=2, bs="tp")+
                                                       s(log(conc), Plant_uo, k=5, bs="fs",m=1),
                                                     data= CO2,method="REML",
                                                     control = list(keepData=TRUE)))[3]


comp_resources[4,"time"] = system.time(bird_mod2 <- gam(count ~ te(week,latitude, bs= c("cc", "tp"),
                                                                 k=c(10,10),m=c(2,2))+
                                                        te(week,latitude,species, bs= c("cc", "tp","re"),
                                                           k=c(10,10,6),m = c(1,1,1)),
                                                      data= bird_move, method="REML", family= poisson,
                                                      control = list(keepData=TRUE)))[3]


comp_resources[5,"time"] = system.time(
  CO2_mod3 <- gam(log(uptake) ~ s(log(conc),k=5,m=2, bs="tp")+
                    s(log(conc),by= Plant_uo, k =5,  bs="ts",m=1)+
                    s(Plant_uo,bs="re",k=12),
                  data= CO2,method="REML",
                  control = list(keepData=TRUE)))[3]



comp_resources[6,"time"] = system.time(
  bird_mod3 <- gam(count ~ te(week,latitude, bs= c("cc", "tp"),
                              k=c(10,10),m=c(2,2)) +
                     te(week,latitude, bs= c("cc", "tp"),
                        k=c(10,10),m=c(1,1),by= species),
                   data= bird_move, method="REML", family= poisson,
                   control = list(keepData=TRUE)))[3]


comp_resources[7,"time"] = system.time(
  CO2_mod4 <- gam(log(uptake) ~ s(log(conc), Plant_uo, k=5,  bs="fs",m=2),
                  data= CO2,method="REML",
                  control = list(keepData=TRUE)))[3]


comp_resources[8,"time"] = system.time(
  bird_mod4 <- gam(count ~ te(week,latitude,species, bs= c("cc", "tp","re"),
                              k=c(10,10,6),m = 2),
                   data= bird_move, method="REML", family= poisson,
                   control = list(keepData=TRUE))
)[3]


comp_resources[9,"time"] = system.time(
  CO2_mod5 <- gam(log(uptake) ~ s(log(conc),by= Plant_uo, k =5,  bs="tp",m=2)+
                    s(Plant_uo,bs="re",k=12), data= CO2,method="REML",
                  control = list(keepData=TRUE))

)[3]

comp_resources[10,"time"] = system.time(
  bird_mod5 <- gam(count ~ te(week,latitude,by=species, bs= c("cc", "tp"),
                              k=c(10,10),m = 2),
                   data= bird_move, method="REML", family= poisson,
                   control = list(keepData=TRUE))
)[3]


comp_resources$model = list(CO2_mod1, bird_mod1, CO2_mod2, bird_mod2,
                            CO2_mod3, bird_mod3,CO2_mod4, bird_mod4,
                            CO2_mod5, bird_mod5)

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
  mutate(`relative time` = `relative time`/`relative time`[1],#scales processing time relative to model 1
         `relative time` = ifelse(`relative time`<10, signif(`relative time`,1), signif(`relative time`, 2)) #rounds to illustrate differences in timing.
         )%>%
  ungroup()%>%
  select(-data_source)

kable(comp_resources_table,format ="latex", caption="Relative computational time and model complexity for different HGAM formulations of the two example data sets from section III. All times are scaled relative to the length of time model 1 takes to fit to that data set. The number of coefficients measures the total number of model parameters (including intercepts). The number of smooths is the total number of unique penalty values estimated by the model.", booktabs = T)%>% #NOTE: change format to "latex" when compiling to pdf, "html" when compiling html
  kable_styling(full_width = F)%>%
  add_header_above(c(" " = 1," "=1, "# of terms"=2))%>%
  group_rows("A. CO2 data", 1,5)%>%
  group_rows("B. bird movement data", 6,10)

## library(brms)
## 
## CO2_mod2_brms <- brm(
##   bf(log(uptake) ~
##        s(log(conc), k=5, m=2) +
##        t2(log(conc), Plant_uo, k=c(5,12),
##                                  bs=c("tp","re"), m=2, full =TRUE)),
##   data = CO2,
##   family = gaussian(),
##   chains = 2,
##   control = list(adapt_delta = 0.95)
## )
## plot(marginal_effects(CO2_mod2_brms), points=TRUE, ask=FALSE, plot=TRUE)
## stancode(CO2_mod2_brms)
#This is just to make sure that the figures occur before the bibliography.
cat('\\FloatBarrier')
