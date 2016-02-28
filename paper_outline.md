# Outline of tutorial paper

## General concept: Teaching ecologists how to use generalized additive models with group-specific predictor functions, using the mgcv package.

The basic idea is that the "fs" basis in mgcv allows you to fit different smooths for different factor levels of a variable, but share a certain amount of information between smooths in that grouping level. That is, when fitting a smooth for factor levels that have few data points, using the "fs" basis will draw those levels towards the global mean smooth for that variable. Conceptually very close to a mixed effect model with random effects for slopes.

**GLS:** _we should talk about subject-specific smooths in general, thus including `by` smooths as a fixed effect alternative. Then we can make the point that with larger *j* (sensu eqn below) a random effect smooth might be more useful generally. Then there are different ways to accomplish these things with **mgcv**. Putting everything in this broader context would be more generally useful than the very specific efficient random smooths of `fs` terms. Also, need to be clear that `fs` terms will tend to have similar amoungs of wigglyness, but `by` terms need not._

**DLM**: _Definitely want to set out in general (and plain) terms what these smooths mean in a non-mathemtatical sense to give folks an idea of the **kind** of models that they can fit. We can dive into some specific examples but we also want to ensure that someone with a new problem that fits into this approach can work that out from our paper (is this asking too much?). Assuming you were writing the above quickly, as smooths aren't really fixed effects :) -- this brings up some terminology issues in terms of how to describe random effects (which are smooths) vs. smooths (which are random effects) vs. combinations of both (which are combinations of both). Need to be careful to not write sentences as confusing as my last one._

***EJP*** _Fair point about my misuse "fixed effects" to describe smooths; I was writing quickly (and sloppily) above. I should have said that the standard `s(x,"by=group")` model assumes that there is smoothing/penalization/hierarchical relationships between coefficients **within** each level of a group, whereas alternative models (`s(x,group,bs="fs")` or  `s(x) +s(x,group,bs="fs")`) allow for hierarchical relationships/smoothing **between** different groups in a single regression. After thinking about this for a while, I think we should actually avoid the "fixed" vs. "random" effects terminology. I don't think I've seen any discussion involving those terms that doesn't result in a lot of confusion about what someone means by fixed or not. Maybe instead use Gelman-esque terminology, like "unpooled", "partially pooled" and "fully pooled", and contrast pooling within groups vs. pooling estimates between groups?_

For a one-dimensional regression, this model is in effect fitting: 

*y <sub>i,j </sub> ~ f(x) + f <sub>j </sub>(x) + e<sub>i,j </sub>*

Here i indicates individual i, j indicates level j of the grouping variable, 
f(x) is the mean function across all levels of the grouping variable, f <sub>j </sub>(x) is the group-level smoother, and e<sub>i,j </sub> is the individual-level error term. By using the "fs" basis, mgcv will penalize the group-level functions f <sub>j </sub>(x) towards zero, thus drawing them closer to the global function.

*figure 1: A simple diagram, illustrating the concept of a global function and group-level functions derived from it*

![](https://github.com/noamross/mixed-effect-gams/blob/master/figures/fig.1%20-%20example%20spline%20basis%20and%20penalties.png?raw=true)


## Why does it work, and how well does it work?

A brief description of the math behind penalized spline regression (aimed at a 
general ecology audience), and the penalty term assumed by mgcv for the "fs" smoother. 

*figure 2: a plot showing the spline basis for a 1-D smooth (probably a simple 
one to understand, like a cubic spline basis. Potentially use a heatmap matrix
to show what the penalty matrix would look like.*

**GLS:** _a 1-D thin-plate spline is pretty easy to grok (just not compute easily) and those are the defaults so perhaps showing both in a 2-panel plot?_

**DLM**: _agree but default is `tprs` and visualising eigenbases is a bit more tricky. Wonder how we can get around this maintaining both simplicity and fidelity to what actually happens... I think the usual way is to explain the thin plate then say "we do this clever trick to reduce dimensionality". Is that what you were thinking of?_

**EJP**: _I agree with Gavin here, that we could definitely include the plots of the smooth basis curves that are calculated by mgcv for the `tprs` basis, without having to go into the details of **how** they're calculated. I've added sample code to generate this figure, as well as an example figure in the figures folder. I've tried to link it into this document; Hopefully uploading it to Github doesn't mess it up._


## How to code the model in R
Discuss how to code the model in R, why we use each term, and pitfalls to avoid. Throughout, I will use 'x' to refer to the independent continuous variable, 'y' to the dependent outcome variable, and 'group' to the grouping factor. 

It may make sense here to show usage and performance, with in-sample and out-of-sample error rates, for simulated data. I've added some code to the github page to demonstrate what I'm thinking of here. 

**DLM**: _I wonder how much performance will be context-dependent here. May be more interesting to try to thoroughly delineate the differences between the models we have or show how misspecification effects interpretation?_

**EJP**: _This was basically what I had in mind; come up with a few simulated case studies where data comes from different types of model, then show how each approach to modelling group-specific curves performs, both in terms of in-sample and out-of-sample data. The cases I was thinking of would be: group-specific data generated from: a) the same underlying function, w/ noise b) an overall function + subject-specific functions c) subject-specific functions, all with the same degree of smoothness d) group-specific functions w/ different degrees of smoothness._

###Basic model formula:

```r
model = gam(y ~ s(x) + s(x, group, bs = "fs", m = 1))
```

**GLS:** _no need for `select = TRUE` here; the `fs` smooths already have penalties on the null space. The `m` normally refers to the order of the penalty, but Bayeen et al says "...and `m=1` requests shrinkage to obtain wiggly random effects" (page 4 of Bayeen et al [arXiv](http://arxiv.org/pdf/1601.02043v1.pdf)). Also not sure you need `s(x)` here for `fs` smooths, at least that's my reading of Bayeen et al, where they don't include include it always. Depends on focus; this model gives a global function of `x` and then difference splines estimated as random effect factor-spline smooths. Without `s(x)` you seem to get the same model but you don't see the global smoother. This model is easier to work with, but not if you want the global smooth._

_Also, at what point do we need `method = "REML"` or `method = "ML"` here? I tend to default to fitting via one of these rather than GCV these days..._

**DLM**: _Yup definitely want to cover this. Esp. since model selection via REML can only take place if the other terms in the model (than `"fs"`) have penalised nullspaces. I see in `?smooth.construct.fs.smooth.spec` that Simon uses either `gamm` or `method="ML"`, need to discuss this too. If we start talking about using `gamm` then we're in a more sticky situation maths-wise as we really need to say something about random effect-smooth term equivalence and what that means. (That doesn't sound like a paper than normal ecologists want to read.)_

**EJP** _I think we should probably not talk about gam vs. gamm in this paper, beyond a brief mention that it is possible; it's beyond the scope of what we're teaching, and I think we'd be too likely to get deep into the weeds. The `method=ML` vs. `method=REML` question is an important one to cover though; I've just been using "GCV.Cp" while testing code, but ML or REML seem like they should be better options. This is one area that I'm not as familiar w/ the GAM literature on, as I don't often use GAM for model selection per se; I typically have a model in mind when fitting. What is the general advice on this? I know from hierarchical regression the general suggestion is to use REML when comparing models differing in random effects (here, smooth terms) and ML when comparing models differing in fixed effects?_

_As I had mentioned in issue, the model `s(x) + s(x,group, bs="fs")` will have different out-of-sample performance than `s(x,group,bs="fs")` by itself, as in the first case, including the global smoother allows different groups to share information on the shapes of the curves (rather than just the degree of smoothness). It also seems like good advice to include the global smooth generally; if the average curve across groups is close to flat, mgcv will shrink it to zero anyway, at the cost of slightly longer computation time._

Alternate forms, and their effect on performance: 

1. `gam(y ~ s(x, group, bs = "fs"), select = TRUE)`

    *this will over-smooth all functions fitted, as it's trying to pull all groups toward a flat function, rather than the global function.*

    **_GLS_:will it?; it's mainly redundant as `fs` smooths have penalties on the null space already, which is what you are asking for with `select = TRUE`; `s(x, by = group, id = 1), select = TRUE` would be somewhat equivalent** 

**DLM**: _`smooth.r:1709` shows the "fs" construction is similar to what's in pp160-161 (though I have an annotation on my copy that says the equation at the top of page 161 is not correct and should be a matrix with 1 on the the diagonal on the last few rows only). Fitting with and without `select=TRUE` gives identical penalty matrices for the "fs" term -- so this is exactly equivalent I think._

_BUT! This will make a difference when other terms are included, especially if one wants to use `method="REML"` (see above).

_Probably also worth consulting [Giampiero Marra's thesis (PDF)](http://opus.bath.ac.uk/22536/1/UnivBath_PhD_2010_G_Marra.pdf) on this matter._

_I had misunderstood what `bs=fs` was doing in this case. You're right, that it should have no tendency to pull the smooth terms toward zero, and the select=T is entirely redundant in this case. This model won't have any tendency to draw all smooth terms to a global smoother, though, as there is no penalty pooling the same basis functions in different groups towards a common value._


2. `gam(y ~ s(x) + s(x, group, bs = "fs", m = 1) + group, select = TRUE)`

    *This will raise an error, as it's overparameterizing the model. mgcv automatically gives each level of the group a separate intercept.*

     **_GLS_:in `fs` smooths. For `by` smooths are subject to centring constraints and hence you need the parameteric `group` term.**

**DLM**: _So let's add "centring constraints" to the "to explain" list :)_


3. `gam(y ~ s(x) + s(x, group, bs = "fs"), select = TRUE)`

    *the m=1 term for thin-plate splines means that each level of the group will not have an independent trend term.  This may or may not cause problems when fitting.*

    **_GLS_: I'm not quite sure what you mean here; these models do fit "separate" splines. Normally `m = 1` just means the penalty is on the first derivative. Bayeen et al imply this is something to do with shrinkage to obtain wiggly smooths (see above). A quick test with data from `?factor.smooth.interaction` suggests that `m = 1` gives wigglier smooths than without.**

*DLM*: _As above, I think this warrants a look at the source of `smooth.construct.fs.smooth.spec` to really see precisely what's happening (I'm also a little lost by the explanation above, so I'll let you elaborate more before I wade in here)._

**EJP**: _What I had meant was that when setting `m=1`, you only had an intercept term in the null-space, not a linear term. This also implies that if you are extrapolating outside of the range of a specific group when m=1, the model will predict a flat group-specific function. After digging into the code, the `m=1` term does the same thing in the "fs" case as in the standard "tp" case: it does just penalize first derivatives rather than second (and remove the linear term from the null space, I assume as a consequence of penalizing 1st instead of 2nd derivatives)._

4. `gam(y ~ s(x) + s(x, group, bs= "fs"), select = FALSE)`

    *This removes the slight penalty term that pulls the smooths toward zero (rather than a straight line). Not sure yet what affect this has on model performance overall.*

    **_GLS_:This is redundant; fitting like this (assuming we drop the `s(x)` term is how you are supposed to fit these smooth factor interactions because they already have the null-space penalties added. I don't think this removes the null-space penalties in the `fs` term; it just doesn't apply a null-space penalty to the global term, which ordinarily wouldn't have one anyway, but it would if you didn't include `s(x)`**

**DLM**: _See my comments above._


5. `gam(y ~ group + s(x, by = group))`

    **_GLS_:This formulation fits full fixed effects terms for x within the levels of `group`; we have a smooth effect of `x` for each group. This also implies that each group has a potentially different set of basis functions and smoothing parameters.**

**DLM**: _Again not keen on the use of the term "fixed effects" here :)_

6. `gam(y ~ group + s(x, by = group, id = 1), select = TRUE)`

    **_GLS_:This is almost equivalent to the `gam(y ~ s(x, group, bs = "fs"))` form --- see `?factor.smooth.interaction`. The `id = 1` links all smooths such that they share the same basis functions and have the same smoothness penalties.**

**DLM**: _Will they not share the same basis functions without the `id` argument? I think the key here is to emphasize the shapes of the smooths will be different by their amount of wigglyness will remain the same (and say why this is useful IRL)._

7. `gam(y ~ group + s(x) + s(x, by = group, m = 1))`

    **_GLS_:This fits a global spline for `x` plus difference smooths for `x` by group, where penalty is based on first derivative. Essentially the penalty should tend to zero as the differences between the `group`-specific splines and the global spline go to zero. This again is a fixed-effects version/approach.**

### Visualizing the model:

* Using plot and predict functions to visualize curves. 
* Discuss the importance of interpreting group-level plots as a sum of mean and 
group-level functions. 

**DLM**: _The default `plot` method is nice but I think we can make something a little more visually appealing (there are **so many colours** in the default method), specifically highlighting the "global" function nicely._


**EJP**: _I agree, don't think the default 'fs' plot is the best way of depicting these smooths, especially as it doesn't show confidence intervals. Were you thinking of us writing a new plotting function, or just demonstrating alternate ways of plotting the models?_


### Extending the model: 

* The use of 'by=' terms to include fixed effects or group. Although after playing around
a bit with code for this, this can be very finicky; it often raises errors about 
having more coefficients that data, in cases that it should work.

    **GLS:** _I haven't had problems with this, as long as *j* (number of groups) * *k* (basis dimension) isn't larger than your data set size or size of an individual group._

**DLM**: _I've found this tricky once you get up to a fair size of *j*, even with relatively large data sets, though using `id` may solve this._

**EJP** _

* using this for generalized additive models with alternate error terms.

**DLM**: _what does this mean?_

**EJP** _I meant that we should demonstrate how these alternative models perform when `family` equals something other than Gaussian. But thinking about it further, it probably makes sense to just say that the `s(x,group,bs="fs")` approach works for non-Gaussian data as well, without going into details; the real-world examples we use will be using mostly non-Gaussian data anyway._ 

## Where this approach should work well, where it can fail

As with any penalty based approach, this will not always work. We should discuss cases 
when this model will work well (when a function can be treated as a sum of a overall and 
individual level function).

We can also discuss when it will work less well, such as:
*Most importantly: when the true functions for different groups do not have similar degrees of smoothness. If two groups differ substantially in their degree of smoothness, 
* when functions differ by a multiplicative value (f <sub>j</sub>(x) = f <sub>j </sub>(x)* f <sup>*</sup> (x)). This is a sub-case of the first issue. 
* when functions are differ by a phase shift (f <sub>j</sub>(x) = f <sup>*</sup> (x+a<sub>j</sub>)). In this case, `s(x,group,bs="fs")` should work well (as all functions will have equivalent degree of smoothness), but `s(x)+s(x,group,bs="fs")` shouldn't offer any improvement in predictive accuracy.
* when there are hidden parameters that lead to a multi-modal global function; if there are multiple underlying groups that each have their own mean function (say, one function for males and one for females) the mean function derived from averaging across the hidden variation will not represent any of the hidden groupings well. 


## Worked examples

* Walleye recruitment: A data set of roughly 100 lakes with walleye recruitment (counts of young of year)
time series sampled irregularly for each lake from 1989 onward. Use the "fs" method
to derive how recruitment has changed over time globally, and how it varies between lakes
* Walleye age structure: A smaller data set (a subset of the recruitment lakes) where counts of adults at multiple lengths were taken, as part of a population estimate. Use this technique to measure the mean age structure across all lakes, how it's changed over time, and how much inter-lake variation there is. 
* Noam: time-series of bat virus antibody prevalence over time, for multiple cohorts of bats over 
their first year, which reveals the seasonality of when they are exposed to a virus

**GLS:** _we're going to have be somewhat careful when fitting time series to not allow (or warn users) too complex smooth terms unless we explicitly account for residual temporal autocorrelation. This really depends on how much data is available and how strong any autocorrelation is. Thinking about the examples, unless they are quite hi-res I doubt this will be an issue._

**DLM**: _I'm experimenting with these in some spatial models I'm working on, that seems like a fairly general situation that might be interesting to folks. More on this soon!_

**EJP** _In the walleye case at least, I can attest that there is very little auto-correlation... I honestly wish there was more. I came in to the project I'm working on with it hoping to construct dynamic population models, but it turns out walleye recruitment year to year is pretty darn close to a random variable. All we've really found so far is a large-scale trend of declining recruitment from 1995 on._