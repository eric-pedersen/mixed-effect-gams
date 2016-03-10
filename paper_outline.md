# Outline of tutorial paper

**General concept: Teaching ecologists how to use generalized additive models with group-specific predictor functions, using the mgcv package.**

##The basic idea:

Penalized spline regression is one of the most powerful techniques available for modelling smooth nonlinear relationships between predictor (x) and outcome (y) variables. The *mgcv* package for R is the most popular tool for ecologists using this technique. *mgcv* allows for a high degree of flexibility in modelling different functional relationships, including allowing complex interactions between smooth predictors. However, in many ecological applications data fall into discrete groups. We often want to know both how the functional relationship (f(x)) between x and y vary between groups, and if there is a strong relationship on average between x and y across groups.

The *mgcv* package has multiple tools to allow estimate both global and group-level functions, but many ecologists are not aware of the different options available, what the trade-offs between these smoothers might be, and what different choices of group-wise smoothers assumes about the underlying ecological relationships. This paper will cover each of the approaches to group-level smoothing, the options for each one and why a user might choose it, and demonstrate the different approaches on a few case studies. 

**EJP**: *I think we should choose one or couple different examples of groups here to use throughout, such as different species responding to a common gradient, or age-structures of a fish species varying between lakes.*

### Intro to generalized additive modelling with penalized splines

A brief description of the math behind penalized spline regression (aimed at a 
general ecology audience). Include a brief, non-mathematical description of penalization as a tool for model estimation. This would also be a good place to talk about penalizing spline wiggliness vs. penalizing the null space for a given spline. We can make the connection here to hierarchical modelling and statistical pooling: large smoothness penalties pool basis functions toward their null-space (the shape of which depends on the splines used). Large null-space penalties draw the smooth toward a trend-less function. High values for both penalties will result in a flat function. 

The basic concepts we should teach:

1. Defining what we mean by a "smooth" function, and how one function can be smoother than another
2. Basis functions
3. Penalty matrices
4. Penalty selection
5. Null spaces and penalties on null spaces
6. Estimating uncertainty for splines
7. Tensor splines and interactions

![](https://raw.githubusercontent.com/noamross/mixed-effect-gams/master/figures/fig.1%20-%20example%20spline%20basis%20and%20penalties.png)
*figure 1: A) Basis functions and B) penalty matrices for two types of penalized splines with 6 basis functions. Top: cubic spline. Bottom: thin-plate spline.*


### Alternate models of between-group variation in f(x)
There are two simple ways for estimating functions for different groups: a) fit a single function for all groups (the fully pooled case) or b) fit a separate function (with its own smoothness estimate) for each group (the un-pooled case). The first case will ignore any inter-group variation, and if different groups have radically different response functions, it will tend to estimate a flat line, or a single curve that does not fit well for any given group. The second case allows for intergroup variability, but does not allow any sharing of information between groups. If each group has only a few data points in it, this will estimate a noisy relationship for each group. 

As with hierarchical regression, we can gain some of the benefits of both these approaches by partially pooling our estimated group-level functions toward one another. In classic hierarchical linear regression, this would mean estimating a single global parameter and penalizing parameter estimates for each group that deviated from that. However, the situation is a bit more complicated when partially pooling functions, as we're already performing pooling within each group simply by estimating smooth curves. We now have two types of penalty: penalizing group-level functions that differ from an overall average, or penalizing group-level functions for being less smooth than the average function.


There are five types of model that *mgcv* can fit, from most to least pooling of information between groups (figure 2):

1. Fully pooled. All groups are assumed to share the same functional relationship between x and y (*y<sub>j</sub>~f(x)* for all *j* in *g* groups).
2. Global smooth with group-specific deviations from this smooth, with all deviations being similarly smooth (*y<sub>j</sub>~f(x) + f<sub>j</sub>(x)*, *p<sub>j</sub>=p* for all *j*). 
3. Global smooth with group-specific deviations, but groups can differ in how smooth their deviation is (*y<sub>j</sub>~f(x) + f<sub>j</sub>(x)*, *p<sub>j</sub>* differs for each *j*).
4. Each group has its own group-specific smooth, but all groups share the same smoothness penalty (*y<sub>j</sub>~ f<sub>j</sub>(x)*, *p<sub>j</sub>=p* for all *j*).
5. Fully un-pooled: Each group has its own smooth and penalty term. (*y<sub>j</sub>~ f<sub>j</sub>(x)*, *p<sub>j</sub>* differs for each *j*). 


Here *j* indicates the *jth* level of the grouping variable, *g* is the total number of different groups, f(x) is the mean function across all levels of the grouping variable, f<sub>j</sub>(x) is the group-level smoother for group *j*, and *p<sub>j</sub>* is the penalty term for the smoother for group *j*. *n* refers to the total number of data points, and *n_j* is the number of independent observations in group *j* (*sum(n<sub>j</sub>) = n*).  

Realistically, all ecological models will fall into models 3 or 5: there will always be at least some form of functional differences between different groups. The issue becomes how well does a given model fit observed data when we only have a limited number of data points to estimate each relationship. There will be a bias-variance trade-off between models 1-5. If the relationship between x and y differs strongly between each group, model 1 will fit each group poorly. However, if on average *n<sub>j</sub>* is small, model 5 may be able to estimate intra-group relationships more accurately, but the predictions will be highly variable. 

<img src="https://raw.githubusercontent.com/noamross/mixed-effect-gams/master/figures/fig.2%20-%20alternate%20models%20of%20functional%20variability.png" width="500">

*figure 1: A simple diagram, illustrating the concept of a global function and group-level functions derived from it*


## How to code the model in R
Discuss how to code the model in R, why we use each term, and pitfalls to avoid. Throughout, I will use 'x' to refer to the independent continuous variable, 'y' to the dependent outcome variable, and 'group' to the grouping factor. 

Start with the basics of how to code a smooth term in *mgcv*, with group level intercepts. Basic ideas to cover:

1. s() terms;
2. Specifying bases via `bs=` term;
3. Specifying `k=` terms for smooths;
4. specifying intercepts and "standard" random effects;
5. `method = "REML"` or `method="ML"` instead of default "GCV";
6. Adding null-space penalties via `select=T`; 
7. The centering constraint assumption for smooths (contrast this w/ `bs="fs"` later);
8. Using the `predict.gam` function to extract global and group-specific functions and standard error estimates.

**EJP**: *It probably makes sense to try and put all the 
different code structures into a table, grouped by degree of pooling and whether the null-space is penalized or not*

###Model 1: only global smooth terms

1. 
```r
model = gam(y ~ s(x), method="REML")
```
2. 
```r
model = gam(y ~ s(x)+ group, method="REML")
```
3. 
```r
model = gam(y ~ s(x)+ s(group,bs="re"), method="REML", select=T)
```

The first model fits only the global trend, the second fits a global trend w/ group-specific intercepts, and the third penalizes both the smooth term and the group-specific intercepts toward the global mean.


###Model 2: Global smoother w/ group-specific deviations, with equal smoothness penalties for all groups 
**EJP** *note that I'm leaving out essentially identical equivalents here (ie: `gam(y ~ s(x)+s(x, by=group, id=1),select=T, method="REML")` should give almost identical output to using `bs="fs"`*

1. 
```r
gam(y ~ s(x)+s(x, group, bs = "fs"), method="REML")
```
2. 
```r
gam(y ~ ti(x)+ti(x,group, bs=c("tp","re")+ti(group, bs="re"), 
method="REML",select=F)
```
3. 
```r
gam(y ~ ti(x)+ti(x,group, bs=c("tp","re")+ ti(group,bs="re"), 
select=T, method="REML")
```

### Model 3: Global smooth term w/ group-specific deviations, with different smoothness penalties for each group

1. 
```r
gam(y~s(x) + s(x,by=group)+s(group,bs="re"), select=F, 
method="REML")
``` 
2. 
```r
gam(y~s(x) + s(x,by=group)+s(group,bs="re"), select=T, 
method="REML")
``` 

The second case here will allow for null-space penalties for both the global and group-specific functions. In general, #2 is probably a better choice than #1, as it will allow mgcv to smooth the null spaces of each group-term to zero.

**EJP**: *I still haven't convinced myself of the last point; I'm not sure yet if mgcv will set this model up so that the individual-level smooths created via the by= term will have an empty null-space anyway, by applying centering and other constraints.* 

### Model 4: No global smoother, group-level smooths share a penalty factor

1. 
```r
gam(y ~ s(x, group, bs = "fs"), method="REML")
```
2. 
```r
gam(y ~ te(x,group, bs=c("tp","re")+ s(group,bs="re"),
 method="REML",select=T)
```
3. 
```r
gam(y ~ te(x,group, bs=c("tp","re")+s(group,bs="re"), 
select=T, method="REML")
```

### Model 5: No shared information between groups
1. 
```r
gam(y~ s(x,by=group)+s(group,bs="re"), 
select=F, method="REML")
``` 
2. 
```r
gam(y~ s(x,by=group)+s(group,bs="re"), 
select=T,  method="REML")
``` 

## Why does it work, and how well does it work?

It may make sense here to show usage and performance, with in-sample and out-of-sample error rates, for simulated data. I've added some code to the github page to demonstrate what I'm thinking of here. 

This would be most useful to show a limited number of cases, w/ or w/out a global function added, and w/ or w/out different levels of smoothness between group-level functions, and finally with a few different numbers of groups and numbers of total data points, to try and illustrate where the bias-variance trade-offs occur here. 




## Visualizing the model:

* Using plot and predict functions to visualize curves. 
* Discuss the importance of interpreting group-level plots as a sum of mean and 
group-level functions. 
* New code to visualize smooths and uncertainty in smooths, and how to plot these curves in *ggplot2*. 



## Extending the model: 

* Modelling interactions between smooth terms that differ between groups: `s(x,y,group, bs="fs")` or `te(x,y,group, bs=c("tp","tp","re")`. 
* The use of `by=` terms to allow the global function to vary based on known covariates `(s(x,by=cov))`, or to allow the estimated smoothing parameters to differ between groups (`s(x,group,bs="fs",by=cov)` or `te(x,group, bs=c("tp","re"), by=cov)`. 
* m=1 vs. m=2 terms for the `bs="fs"` approach. 

**EJP** *For the last point, I'm personally still not sure when each one is a better approach though, so I'm not sure how good my advice on this would be... m=1 will give more "wiggly" but flatter smooths in general, as its penalizing the square of the first derivative rather than the second, but it's hard to tell when that would be better or worse when looking at group-level smooths.*


## Where to use each approach, and how to compare different approaches
We should discuss some of the basic testing approaches:

1. Anova (need to note to use `method="ML"` when making model comparisons)
2. check.gam
3. diagnostic plots

Important trade-offs to consider:

1. Bias vs. variance. When you have only a few data points per group, lower models (1&2) will generally perform better than higher models (3-5), but it will depend on the shape of the smooths, and the degree of observational/intrinsic noise in the data.
2. Accuracy vs. computing time: these models can take a long time to calculate.
3. Forecasting in-sample vs. interpolating within groups vs. extrapolating outside the range of any group: models 1-3 will in general help more with prediction if many groups are missing blocks of data, as the global trend can be used to interpolate the missing blocks. Further, if you want to make predictions for groups you haven't observed, you need the global trend. If the focus is on making accurate predictions within the range of each group, however, models 3&5 may perform better, as groups may differ substantially in how smooth their functional relationship is, and forcing all groups to have the same degree of smoothness could result in some groups being heavily over-penalized, and others substantially under-penalized. 



## Worked examples

* EJP: Walleye recruitment: A data set of roughly 100 lakes with walleye recruitment (counts of young of year)
time series sampled irregularly for each lake from 1989 onward. Use the "fs" method
to derive how recruitment has changed over time globally, and how it varies between lakes
* EJP: Walleye age structure: A smaller data set (a subset of the recruitment lakes) where counts of adults at multiple lengths were taken, as part of a population estimate. Use this technique to measure the mean age structure across all lakes, how it's changed over time, and how much inter-lake variation there is. 
* Noam: time-series of bat virus antibody prevalence over time, for multiple cohorts of bats over 
their first year, which reveals the seasonality of when they are exposed to a virus
* DLM: spatial models (more details soon!)

For all cases, we need to be careful to note that both time series and spatial data can have short-time scale auto-correlation that straight gam models won't pick up on. This can be handled w/ autoCor functions in `gamm` or `gamm4`, but is out of the scope of this paper, and should be handled with care.  