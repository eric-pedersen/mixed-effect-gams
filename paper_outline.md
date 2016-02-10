# Outline of tutorial paper

## General concept: Teaching ecologists how to use functional mixed effects
using the mgcv package and the "fs" basis. The basic idea is that the "fs"
basis in mgcv allows you to fit different smooths for different factor levels of
a variable, but share a certain amount of information between smooths in that
grouping level. That is, when fitting a smooth for factor levels that have few
data points, using the "fs" basis will draw those levels towards the global mean
smooth for that variable. Conceptually very close to a mixed effect model with
random effects for slopes.

**GLS:** _we should talk about subject-specific smooths in general, thus including `by` smooths as a fixed effect alternative. Then we can make the point that with larger *j* (sensu eqn below) a random effect smooth might be more useful generally. Then there are different ways to accomplish these things with **mgcv**. Putting everything in this broader context would be more generally useful than the very specific efficient random smooths of `fs` terms. Also, need to be clear that `fs` terms will tend to have similar amoungs of wigglyness, but `by` terms need not._

For a one-dimensional ~~Gaussian~~ regression, this model is in effect fitting: 

*y <sub>i,j </sub> ~ f(x) + f <sub>j </sub>(x) + e<sub>i,j </sub>*

Here i indicates individual i, j indicates level j of the grouping variable, 
f(x) is the mean function across all levels of the grouping variable, f <sub>j </sub>(x) is the group-level smoother, and e<sub>i,j </sub> is the individual-level error term. By using the "fs" basis, mgcv will penalize the group-level functions f <sub>j </sub>(x) towards zero, thus drawing them closer to the global function.

*figure 1: A simple diagram, illustrating the concept of a global function and group-level functions derived from it*




## Why does it work, and how well does it work? 
A brief description of the math behind penalized spline regression (aimed at a 
general ecology audience), and the penalty term assumed by mgcv for the "fs" smoother. 

*figure 2: a plot showing the spline basis for a 1-D smooth (probably a simple 
one to understand, like a cubic spline basis. Potentially use a heatmap matrix
to show what the penalty matrix would look like.*

**GLS:** _a 1-D thin-plate spline is pretty easy to grok (just not compute easily) and those are the defaults so perhaps showing both in a 2-panel plot?_

## How to code the model in R
Discuss how to code the model in R, why we use each term, and pitfalls to avoid. Throughout, I will use 'x' to refer to the independent continuous variable, 'y' to the dependent outcome variable, and 'group' to the grouping factor. 

It may make sense here to show usage and performance, with in-sample and out-of-sample 
error rates, for simulated data. I've added some code to the github page to demonstrate what I'm thinking of here. 

###Basic model formula: 

```r
model = gam(y ~ s(x) + s(x, group, bs = "fs", m = 1))
```

**GLS:** _no need for `select = TRUE` here; the `fs` smooths already have penalties on the null space. The `m` normally refers to the order of the penalty, but Bayeen et al says "...and `m=1` requests shrinkage to obtain wiggly random effects" (page 4 of Bayeen et al [arXiv](http://arxiv.org/pdf/1601.02043v1.pdf)). Also not sure you need `s(x)` here for `fs` smooths, at least that's my reading of Bayeen et al, where they don't include include it always. Depends on focus; this model gives a global function of `x` and then difference splines estimated as random effect factor-spline smooths. Without `s(x)` you seem to get the same model but you don't see the global smoother. This model is easier to work with, but not if you want the global smooth._

Alternate forms, and their effect on performance: 

1. `gam(y ~ s(x, group, bs = "fs"), select = TRUE)`

    *this will over-smooth all functions fitted, as it's trying to pull all groups toward a flat function, rather than the global function.*

    **_GLS_:will it?; it's mainly redundant as `fs` smooths have penalties on the null space already, which is what you are asking for with `select = TRUE`; `s(x, by = group, id = 1), select = TRUE` would be somewhat equivalent** 

2. `gam(y ~ s(x) + s(x, group, bs = "fs", m = 1) + group, select = TRUE)`

    *This will raise an error, as it's overparameterizing the model. mgcv automatically gives each level of the group a separate intercept.*

     **_GLS_:in `fs` smooths. For `by` smooths are subject to centring constraints and hence you need the parameteric `group` term.**

3. `gam(y ~ s(x) + s(x, group, bs = "fs"), select = TRUE)`

    *the m=1 term for thin-plate splines means that each level of the group will not have an independent trend term.  This may or may not cause problems when fitting.*

    **_GLS_: I'm not quite sure what you mean here; these models do fit "separate" splines. Normally `m = 1` just means the penalty is on the first derivative. Bayeen et al imply this is something to do with shrinkage to obtain wiggly smooths (see above). A quick test with data from `?factor.smooth.interaction` suggests that `m = 1` gives wigglier smooths than without.**

4. `gam(y ~ s(x) + s(x, group, bs= "fs"), select = FALSE)`

    *This removes the slight penalty term that pulls the smooths toward zero (rather than a straight line). Not sure yet what affect this has on model performance overall.*

    **_GLS_:This is redundant; fitting like this (assuming we drop the `s(x)` term is how you are supposed to fit these smooth factor interactions because they already have the null-space penalties added.**

5. `gam(y ~ group + s(x, by = group))`

    **_GLS_:This formulation fits full fixed effects terms for x within the levels of `group`; we have a smooth effect of `x` for each group. This also implies that each group has a potentially different set of basis functions and smoothing parameters.**

6. `gam(y ~ group + s(x, by = group, id = 1), select = TRUE)`

    **_GLS_:This is almost equivalent to the `gam(y ~ s(x, group, bs = "fs"))` form --- see `?factor.smooth.interaction`. The `id = 1` links all smooths such that they share the same basis functions and have the same smoothness penalties.**

7. `gam(y ~ group + s(x) + s(x, by = group, m = 1))`

    **_GLS_:This fits a global spline for `x` plus difference smooths for `x` by group, where penalty is based on first derivative. Essentially the penalty should tend to zero as the differences between the `group`-specific splines and the global spline go to zero. This again is a fixed-effects version/approach.**

### Visualizing the model:

* Using plot and predict functions to visualize curves. 
* Discuss the importance of interpreting group-level plots as a sum of mean and 
group-level functions. 

### Extending the model: 

* The use of 'by=' terms to include fixed effects or group. Although after playing around
a bit with code for this, this can be very finicky; it often raises errors about 
having more coefficients that data, in cases that it should work. 
* using this for generalized additive models with alternate error terms.

## Where this approach should work well, where it can fail

As with any penalty based approach, this will not always work. We should discuss cases 
when this model will work well (when a function can be treated as a sum of a overall and 
individual level function).

We can also discuss when it will work less well, such as: 
* when functions differ by a multiplicative value (f <sub>j</sub>(x) = f <sub>j </sub>(x)* f <sup>*</sup> (x))
* when functions are differ by a phase shift (f <sub>j</sub>(x) = f <sup>*</sup> (x+a<sub>j</sub>))
* when there are hidden parameters that lead to a "multi-modal" global function; if
there are multiple underlying groups that each have their own mean function (say, one function for males and one for females) the mean function derived from averaging across tha hidden variation will not represent any of the hidden groupings well. 

## Worked examples

* Walleye recruitment: A data set of roughly 100 lakes with walleye recruitment (counts of young of year)
time series sampled irregularly for each lake from 1989 onward. Use the "fs" method
to derive how recruitment has changed over time globally, and how it varies between lakes
* Walleye age structure: A smaller data set (a subset of the recruitment lakes) where counts of adults at multiple lengths were taken, as part of a population estimate. Use this technique to measure the mean age structure across all lakes, how it's changed over time, and how much inter-lake variation there is. 
* Noah: time-series of bat virus antibody prevalence over time, for multiple cohorts of bats over 
their first year, which reveals the seasonality of when they are exposed to a virus

**GLS:** _we're going to have be somewhat careful when fitting time series to not allow (or warn users) too complex smooth terms unless we explicitly account for residual temporal autocorrelation. This really depends on how much data is available._

