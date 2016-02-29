library(dplyr)
library(ggplot2)
library(mgcv)
library(tidyr)

set.seed(2)

#### Generating predictor data ####
# Generate data w/ 10 groups, and 100 data points per group, as well as 
# group-specific coefficients for the height and displacement of the group-specific
# function
n = 100
g = 12
full_data = as.data.frame(expand.grid(x=seq(-2,2,length=n), 
                                      group = letters[1:g]))
full_data$indiv = 1:(n*g)

group_width_coefs = rlnorm(g,-0.5,1)
group_shift_coefs = rnorm(g,0,1)

#### Generate y data ####
#assuming each y is a function of an overall function of x and a group-specific one

#overall function: Gaussian curve
full_data$global_function = dnorm(full_data$x)
full_data$group_function = with(full_data,
                                0.4*dnorm(x, group_shift_coefs[as.numeric(group)], 
                                     group_width_coefs[as.numeric(group)]))
full_data$total_function = full_data$global_function + full_data$group_function
full_data$y = rnorm(n*g, full_data$total_function,0.2)

#adding group-specific functions


#### split full data into fitting and test (out of sample) data ####
# by splitting each group at the zero point. 
# For the first half of groups, use data where x<0.
# For the other half use x >0 for the fitting data.
fit_data = full_data %>%
  filter((group%in%letters[1:(g/2)]&x<0)|(group%in%letters[(g/2+1):g]&x>0))
test_data = full_data%>%
  filter(!indiv%in%fit_data$indiv)

base_data_plot = qplot(x,y, data= full_data)+
  facet_wrap(~group)+
  geom_point(data=fit_data, col="blue")+
  geom_line(aes(y= total_function),size=2,col="red")+
  geom_smooth(method=gam,formula=y~s(x, k=20))+
  theme_bw()

print(base_data_plot)


#### Models listed from least pooled (model_1) to most pooled (model_5) ####
# models 2a/2b/2c and 4a/4b/4c represent three ways of parameterizing (respectively)
# a model w/ only group-level smoothers sharing a common penalty term (model 2) and
# a single global smoother plus group-specific smoothers w/ shared smoothing
# parameters between groups (model 4)
model_1 = gam(y~s(x,by=group,k=20),method="REML",data=fit_data)
model_2a = gam(y~s(x,group,bs="fs", k=20),method="REML",data=fit_data)
model_2b = gam(y~s(x,group,bs="fs", k=20,m=1),method="REML",data=fit_data)
model_2c = gam(y~te(x,group,bs=c("tp","re"), k=c(20,g))+
                 s(group,bs="re"),method="REML",data=fit_data)
model_3 = gam(y~s(x,k=20) +s(x,by=group,k=20),method="REML",
              data=fit_data,select=T)
model_4a = gam(y~s(x,k=20)+s(x,group,bs="fs", k=20),method="REML",
               data=fit_data)
model_4b = gam(y~s(x,k=20,m=1)+s(x,group,bs="fs", k=20,m=1),method="REML",
               data=fit_data)
model_4c = gam(y~s(x,k=20)+te(x,group,bs=c("tp","re"), k=c(20,g))+
                 s(group, bs="re"),method="REML",data=fit_data)
model_5 = gam(y~s(x, k=20),method="REML",data=fit_data)



#### Compiling models to summarize output ####
model_list= list(model_1, model_2a,model_2b,model_2c,
                 model_3, model_4a,model_4b, model_4c, model_5)


summary_data = data.frame(model = c("model 1", 
                                    "model 2a","model 2b","model 2c",
                                    "model 3",
                                    "model 4a", "model 4b","model 4c",
                                    "model 5"))

summary_data$rmse_out_of_sample = sapply(model_list, 
                                         function(x)mean((test_data$y - predict.gam(x,test_data))^2))
summary_data$rmse_in_sample = sapply(model_list, 
                                         function(x)mean((fit_data$y - predict.gam(x,fit_data))^2))

summary_long_data  = gather(summary_data,fit_measure, value, -model)


#### Plotting model rmse summary data ####
summary_plot = ggplot(aes(x=model,y=value), data=summary_long_data)+
  geom_point()+
  facet_grid(fit_measure~.,scale="free_y")+
  scale_y_log10()+
  theme_bw()

print(summary_plot)
