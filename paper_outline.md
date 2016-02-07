# Outline of tutorial paper

## General concept: Teaching ecologists how to use functional mixed effects
using the mgcv package and the "fs" basis. The basic idea is that that "fs"
basis in mgcv allows you to fit different smooths for different factor levels of
a variable, but share a certain amount of information between smooths in that
grouping level. That is, when fitting a smooth for factor levels that have few
data points, using the "fs" basis will draw those levels towards the global mean
smooth for that variable. Conceptually very close to a mixed effect model with
random effects for slopes.

For a one-dimensional Gaussian regression, this model is in effect fitting: 

*y <sub>i,j </sub> ~ f \* (x) + f <sub>j </sub>(x) + e<sub>i,j </sub>*

Here i indicates individual i, j indicates level j of the grouping variable, 
f\*(x) is the mean function across all levels of the grouping variable, f <sub>j </sub>(x) is the group-level smoother, and e<sub>i,j </sub> is the 
individual-level error term. By using the "fs" basis, mgcv will penalize the group-level functions f <sub>j </sub>(x) towards zero, thus drawing them closer
to the global function.

*figure 1: A simple diagram, illustrating the concept of a global function and group-level functions derived from it*




## Why does it work, and how well does it work? 
A brief description of the math behind penalized spline regression (aimed at a 
general ecology audience), and the penalty term assumed by mgcv for the "fs" smoother. 

*figure 2: a plot showing the spline basis for a 1-D smooth (probably a simple 
one to understand, like a cubic spline basis. Potentially use a heatmap matrix
to show what the penalty matrix would look like.*


## How to code the model in R
Discuss how to code the model in R, why we use each term, and pitfalls to avoid. Throughout, I will use 'x' to refer to the independent continuous variable, 'y' to the dependent outcome variable, and 'group' to the grouping factor. 

It may make sense here to show usage and performance, with in-sample and out-of-sample 
error rates, for simulated data. I've added some code to the github page to demonstrate what I'm thinking of here. 

###Basic model formula: 

model = gam(y~s(x)+s(x,group, bs= "fs", m=1), select=TRUE)

Alternate forms, and their effect on performance: 

1. gam(y~s(x,group, bs= "fs"), select=TRUE) *this will over-smooth all functions fitted, 
as it's trying to pull all groups toward a flat function, rather than the global function.*
2. gam(y~s(x)+s(x,group, bs= "fs", m=1)+group, select=TRUE).  *This will raise an error, 
as it's overparameterizing the model. mgcv automatically gives each level of the group a separate intercept.*
3. gam(y~s(x)+s(x,group, bs= "fs"), select=TRUE). *the m=1 term for thin-plate splines means that each level of the group will not have an independent trend term. This may or may not cause problems when fitting.*
4. gam(y~s(x)+s(x,group, bs= "fs"), select=FALSE). *This removes the slight penalty term 
that pulls the smooths toward zero (rather than a straight line). Not sure yet what affect
this has on model performance overall.*

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
* Walleye age structure: A smaller data set (a subset of the recruitment lakes) where counts of adults at multiple lengths were taken, as part of a population estimate. Use this technique to measure the mean age structure across all lakes, how it's changed over time,
and how much inter-lake variation there is. 
* Noah: time-series of bat virus antibody prevalence over time, for multiple cohorts of bats over 
their first year, which reveals the seasonality of when they are exposed to a virus

