
################################################
#                                              #
#  Functions for fitting Bivariate Copula      #
#                                              #
################################################

library(VineCopula)


#' @description Find the best fitted Bivariate copulas for consecutive points 
#'              and/or bivariate data
#' @param uts_1 univariate vector 
#' @param uts_2 univariate vector 
#'        Note: if null, them will separate uts_1 into 2 vectors
#' @param print whether to print the result, default to FALSE
#' @param ... additional arguments pass to BiCopSelect
#' @return all objects returned by BiCopSelect augmented with the following entries:
#'         C_{2|1}, C_{1|2}, and c12 of the fitted bivariate copula
fit_bi_copula = function(uts_1, uts_2 = NULL, print = FALSE, ...){
  
  # if only 1 times series vector, then we assume it is X_{t-1} and X_t
  if(is.null(uts_2)){
    n = length(uts_1)
    uvec_1 = uts_1[-c(n)]; uvec_2 = uts_1[-c(1)]
  }else{
    uvec_1 = uts_1; uvec_2 = uts_2
  }
  
  # fit bivariate copulas
  result = BiCopSelect(u1 = uvec_1, u2 = uvec_2, ...)
  
  # add C2|1, C1|2, and c12
  h <- BiCopHfunc(u1 = uvec_1, u2 = uvec_2, result)
  result$C2g1 = h$hfunc1
  result$C1g2 = h$hfunc2
  result$c12 = BiCopPDF(u1 = uvec_1, u2 = uvec_2, result)
  
  if(print){
    cat("Best fitting copula is", result$familyname, "with tau", result$tau, "and AIC", result$AIC)
  }
  
  return(result)
}


# example
gauss_example = function(){
  
  # set parameters
  sim_par1 = matrix(c(0,0.5,0,0),2,2, byrow=T)
  sim_fam = matrix(c(0,1,0,0),2,2, byrow=T)
  rvine_ma = matrix(c(1,1,0,2),2,2, byrow=T)
  RVM <- RVineMatrix(Matrix = rvine_ma, family = sim_fam, par = sim_par1)
  
  # simulation
  set.seed(123)
  sim_sample = RVineSim(N = 5000, RVM = RVM)
  
  # fit
  fit_example = fit_bi_copula(uts_1 = sim_sample[,1], uts_2 = sim_sample[,2], print = TRUE)
  fit_example$par
  # Best fitting copula is Gaussian with tau 0.3304938 and AIC -1365.584
  # 0.4961323
  
  # note: the above result is same as following
  cop1 <- BiCopSelect(sim_sample[,1], sim_sample[,2], familyset = 1:10, indeptest = FALSE, level = 0.05)
  cop1$familyname
  cop1$tau
}



#' @description Find the best fitted Bivariate copulas for univariate Markov order k
#' @param uts univariate vector 
#' @param k Markov order k >= 1, default is k = 1
#' @param print whether to print the result, default to FALSE
#' @param ... additional arguments pass to BiCopSelect
#' @return list of best fitted bivariate copulas, total negative loglikelihood, 
#'         and overall C_{2|1} (transition probability) 
fit_uni_mk = function(uts, k = 1, print = FALSE, ...){
  
  # initialize the value
  cop_list = list()
  tot_nllk = 0
  
  for(i in 1:k){
    
    # fit current level copula
    if(i == 1){
      cop_new = fit_bi_copula(uts_1 = uts, ...)
    }else{
      n_i = length(cop_prev$C1g2)
      cop_new = fit_bi_copula(uts_1 = cop_prev$C1g2[-n_i], uts_2 = cop_prev$C2g1[-1], ...)
    }
    
    # total nllk
    if(k > i){
      tot_nllk = tot_nllk - sum(log(cop_new$c12)[-(1:(k-i))])
    }else if(i == k){
      tot_nllk = tot_nllk - cop_new$logLik
    }
    
    # add the new fitted copula
    cop_list = append(cop_list, list(cop_new)) 
    
    # update previous copula
    cop_prev = cop_new
  }
  
  # AIC
  AIC = 2*k + 2*tot_nllk
  
  # Values
  best_cop = sapply(cop_list, function(x) x$family)
  best_param = sapply(cop_list, function(x) x$par)
  best_param2 = sapply(cop_list, function(x) x$par2)
  best_tau = sapply(cop_list, function(x) x$tau)
  
  if(print){
    cat("Best fitting copula(s) is/are", best_cop, "\n")
    cat("with param", best_param, "\n")
    cat("with param2", best_param2, "\n")
    cat("with tau", best_tau, "\n")
    cat("and overall negative log likelihood", tot_nllk, "\n")
  }
  
  return(list(fitted_cop = cop_list, tot_nllk = tot_nllk, C2g1 = cop_prev$C2g1,
              family_vec = best_cop, par_vec = best_param, par2_vec = best_param2, 
              tau_vec = best_tau))
}


# example
uni_mk_example = function(){
  
  proj = "C:/Users/Chun Fang/OneDrive - UBC/Research"
  source(file.path(proj, "simulation", "sim-univariate-mkp.R"))
  
  ## first example for d=4, all Gumbel
  apar1 = matrix(c(0, 2,2,2, 0,0, 1.5,1.5, 0,0,0,1.8, 0,0,0,0),4,4, byrow=T)
  apar2 = matrix(0,4,4)
  afam4 = matrix(4,4,4) # all Gumbel
  
  uts = sim_uni_mk(n=5000,fam=afam4,param1=apar1,param2=apar2)
  
  # fit 
  fit_uts = fit_uni_mk(uts, k = 3, familyset = c(1,4,5,14), print = TRUE)
}



#' @description Find the best fitted Bivariate copulas for 
#'              Markov order k innovation model
#' @param uts_x univariate vector 
#' @param uts_y univariate vector 
#' @param k Markov order k >= 1, default is k = 1
#' @param print whether to print the result, default to FALSE
#' @param ... additional arguments pass to BiCopSelect
#' @return best fitted families, parameter estimates, total negative loglikelihood
fit_innov_mk = function(uts_x, uts_y, k = 1, print = FALSE, ...){
  
  # fit
  fit_x = fit_uni_mk(uts = uts_x, k = k, ...)
  fit_y = fit_uni_mk(uts = uts_y, k = k, ...)
  fit_innov = fit_bi_copula(uts_1 = fit_x$C2g1, uts_2 = fit_y$C2g1, ...)
  
  # total nllk
  tot_nllk = fit_x$tot_nllk + fit_y$tot_nllk + fit_innov$nllk
  
  if(print){
    
    x_best_cop = sapply(fit_x$fitted_cop, function(x) x$family)
    x_best_param = sapply(fit_x$fitted_cop, function(x) x$par)
    y_best_cop = sapply(fit_y$fitted_cop, function(x) x$family)
    y_best_param = sapply(fit_y$fitted_cop, function(x) x$par)
    
    cat("Best fitting copulas are", "\n")
    cat("X:", x_best_cop, "with parameter", x_best_param, "and nllk", fit_x$tot_nllk, "\n")
    cat("Y:", y_best_cop, "with parameter", y_best_param, "and nllk", fit_y$tot_nllk, "\n")
    cat("C_innov:", BiCopName(fit_innov$family, short = FALSE), "with parameter", fit_innov$par, 
        "and nllk", -fit_innov$logLik)
  }
  
  return(list(C_X = fit_x, C_Y = fit_y, C_innov = fit_innov, tot_nllk = tot_nllk))
}

# example
# need to import sim_bi_innov function from sim-bivariate-mkp-innov.R
innov_example = function(){
  
  proj = "C:/Users/Chun Fang/OneDrive - UBC/Research"
  source(file.path(proj, "simulation", "sim-bivariate-mkp-innov.R"))
  
  innov_dat = sim_bi_innov_mkp(n = 10000, family.innov = 1, par.innov = 0.7, 
                         family.1 = 1, par.1 = 0.3, 
                         family.2 = 1, par.2 = 0.6, trunc = 1000, 
                         unif.margin = TRUE)
  
  testfit = fit_innov_mk(uts_x = innov_dat$dat[,1], uts_y = innov_dat$dat[,2], print = TRUE)
}


# TODO: rotations = FALSE!
#' @description Find the best fitted Bivariate copulas for d-variate
#'              Markov order k innovation model
#' @param dat matrix with number of columns = d; each column corresponds to a 
#'            univariate vector 
#' @param k Markov order k >= 1, default is k = 1
#' @param print whether to print the result, default to FALSE
#' @param ... additional arguments pass to BiCopSelect
#' @return best fitted families, parameter estimates, total negative loglikelihood
fit_innov = function(dat, k = 1, print = FALSE, ...){
  
  # fit
  d = ncol(dat)
  fit_all_col = lapply(1:d, function(i) fit_uni_mk(uts = dat[, i], k = k, ...))
  innov_dat = do.call(cbind, lapply(fit_all_col, function(x) x$C2g1))
  fit_innov = RVineStructureSelect(data = innov_dat, ...)
  
  # total nllk
  tot_nllk = sum(sapply(fit_all_col, function(x) x$tot_nllk)) - fit_innov$logLik
  
  if(print){
    cat("The Negative loglik of the best fitting Rvine innovation model is", tot_nllk)
  }
  
  return(list(univariate_fit = fit_all_col, innov_fit = fit_innov, tot_nllk = tot_nllk))
}



# example
# need to import sim_bi_innov function from sim-bivariate-mk1.R
innov_2_example = function(){
  
  proj = "C:/Users/Chun Fang/OneDrive - UBC/Research"
  source(file.path(proj, "simulation", "sim-bivariate-mk1.R"))
  
  innov_dat = sim_bi_innov(n = 10000, family.innov = 1, par.innov = 0.7, 
                           family.1 = 1, par.1 = 0.3, 
                           family.2 = 1, par.2 = 0.6, trunc = 1000, 
                           unif.margin = TRUE)
  
  testfit = fit_innov(dat = innov_dat$dat, familyset = c(0,1,4,5,14))
  
  test_uni = testfit$univariate_fit
  test_innov = testfit$innov_fit
  
  test_innov$family
  test_innov$par
  
  test_uni[[1]]$fitted_cop[[1]]$family
  test_uni[[1]]$fitted_cop[[1]]$par
  test_uni[[2]]$fitted_cop[[1]]$family
  test_uni[[2]]$fitted_cop[[1]]$par
}
