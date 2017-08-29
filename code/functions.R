library(assertthat)

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

get_r2 = function(model) {
  model_summary = summary(model)
  return(model_summary$dev.expl)
}

create_connected_pmat = function(grouping_var){
  assert_that(is.factor(grouping_var),
              msg = "Markov random fields in mgcv require factor variables, and will not work properly on either character or numeric variables")
  grp_levels = levels(grouping_var)
  n_levels = length(grp_levels)
  if(!all(grp_levels%in% grouping_var)){
    warning("not all grouping levels are present in the variable")
  }
  pmat = matrix(-1, nrow= n_levels,ncol=n_levels)
  diag(pmat) = n_levels-1
  rownames(pmat) =colnames(pmat)=  grp_levels
  return(pmat)
}


refit_model = function(model, algorithm = c("bam", "gamm", "gamm4")) {
  if(!algorithm[1] %in% c("bam", "gamm", "gamm4")) stop("can only refit models using bam, gamm, and gamm4")

  family = get(model$family$family)(link = model$family$link)
  data = model$data
  formula = model$formula
  
  if(algorithm[1]=="gamm4"){
    if(!require(gamm4)) stop("refitting the model using gamm4 requires the gamm4 package to be installed")
    formula = as.formula(gsub("te\\(", "t2\\(",x =  deparse(formula)))
  }
  
  timing = system.time(refit_model <- get(algorithm[1])(formula = formula, data=data,family=family))
  return(list(timing = as.numeric(timing[3]), refit_model =refit_model))
}
