#funtions ####
generate_smooth_func = function(x,n_funcs=1,length_scale=0.25,
                                amp = 1){
  # generates data from a Gaussian process with a squared exponential kernel
  # x is the input data, n_funcs is the number of functions to sample
  # length_scale determines how rapidly the function oscillates, and
  # amp determines the function's overall variance (over ranges of x >> length_scale,
  # var(func) will go to amp.).
  
  stopifnot(is.numeric(x))
  stopifnot(is.vector(x))
  x_dist = as.matrix(dist(x))
  x_cov = amp*exp(-x_dist^2/(2*length_scale^2))
  funcs = rmvn(n_funcs,rep(0, length=length(x)), x_cov)
  if(is.null(dim(funcs))){
    funcs = matrix(funcs, nrow=1)
  }
  return(t(funcs))
}

calc_rmse = function(fit,y,digits=2) round(sqrt(sum((y-fit)^2)),2)

