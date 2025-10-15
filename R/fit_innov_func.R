
################################################
#                                              #
#  Functions for fitting innovation model      #
#                                              #
################################################

# library(VineCopula)


#' Find the best fitted Bivariate copulas for consecutive points and/or bivariate data
#'
#' @param uts_1 univariate vector
#' @param uts_2 univariate vector
#'        Note: if null, them will separate uts_1 into 2 vectors
#' @param print whether to print the result, default to FALSE
#' @param ... additional arguments pass to BiCopSelect
#' @returns all objects returned by BiCopSelect augmented with the following entries:
#'         C_{2|1}, C_{1|2}, and c12 of the fitted bivariate copula
#' @examples
#' # set parameters
#' sim_par1 = matrix(c(0,0.5,0,0),2,2, byrow=T)
#' sim_fam = matrix(c(0,1,0,0),2,2, byrow=T)
#' rvine_ma = matrix(c(1,1,0,2),2,2, byrow=T)
#' RVM <- RVineMatrix(Matrix = rvine_ma, family = sim_fam, par = sim_par1)
#'
#' # simulation
#' set.seed(123)
#' sim_sample = RVineSim(N = 5000, RVM = RVM)
#'
#' # fit
#' fit_example = fit_bi_copula(uts_1 = sim_sample[,1], uts_2 = sim_sample[,2], print = TRUE)
#' fit_example$par
#' # Best fitting copula is Gaussian with tau 0.3304938 and AIC -1365.584
#' # 0.4961323
#'
#' # note: the above result is same as following
#' cop1 <- BiCopSelect(sim_sample[,1], sim_sample[,2], familyset = 1:10, indeptest = FALSE, level = 0.05)
#' cop1$familyname
#' cop1$tau
#'
#' @export
fit_bi_copula = function(uts_1, uts_2 = NULL, print = FALSE, ...){

  # if only 1 times series vector, then we assume it is X_{t-1} and X_t
  if(is.null(uts_2)){
    n = length(uts_1)
    uvec_1 = uts_1[-c(n)]; uvec_2 = uts_1[-c(1)]
  }else{
    uvec_1 = uts_1; uvec_2 = uts_2
  }

  # fit bivariate copulas
  result = VineCopula::BiCopSelect(u1 = uvec_1, u2 = uvec_2, ...)

  # add C2|1, C1|2, and c12
  h <- VineCopula::BiCopHfunc(u1 = uvec_1, u2 = uvec_2, result)
  result$C2g1 = h$hfunc1
  result$C1g2 = h$hfunc2
  result$c12 = VineCopula::BiCopPDF(u1 = uvec_1, u2 = uvec_2, result)

  if(print){
    cat("Best fitting copula is", result$familyname, "with tau", result$tau, "and AIC", result$AIC)
  }

  return(result)
}



#' Find the best fitted Bivariate copulas for univariate Markov order k
#'
#' @param uts univariate vector
#' @param k Markov order k >= 1, default is k = 1
#' @param print whether to print the result, default to FALSE
#' @param separate Logical; whether log-likelihoods are returned point wisely
#'              (default: separate = FALSE)
#' @param ... additional arguments pass to BiCopSelect
#' @returns list of best fitted bivariate copulas, total transitional loglikelihood,
#'         and overall C_{2|1} (transition probability)
#' @examples
#' #source(file.path(proj, "simulation", "sim-univariate-mkp.R"))
#'
#' ## first example for d=4, all Gumbel
#' apar1 = matrix(c(0, 2,2,2, 0,0, 1.5,1.5, 0,0,0,1.8, 0,0,0,0),4,4, byrow=T)
#' apar2 = matrix(0,4,4)
#' afam4 = matrix(4,4,4) # all Gumbel
#'
#' uts = sim_uni_mk(n=5000,fam=afam4,param1=apar1,param2=apar2)
#'
#' # fit
#' fit_uts = fit_uni_mk(uts, k = 3, familyset = c(1,4,5,14), print = TRUE)
#'
#' @export
fit_uni_mk = function(uts, k = 1, print = FALSE, separate = FALSE, ...){

  # initialize the value
  cop_list = list()
  trans_logLik = rep(0, times = (length(uts)-k))

  for(i in 1:k){

    # fit current level copula
    if(i == 1){
      cop_new = fit_bi_copula(uts_1 = uts, ...)
    }else{
      n_i = length(cop_prev$C1g2)
      cop_new = fit_bi_copula(uts_1 = cop_prev$C1g2[-n_i], uts_2 = cop_prev$C2g1[-1], ...)
    }

    # transitional loglik
    if(k > i){
      trans_logLik = trans_logLik + log(cop_new$c12)[-(1:(k-i))]
    }else if(i == k){
      trans_logLik = trans_logLik + log(cop_new$c12)
    }

    # add the new fitted copula
    cop_list = append(cop_list, list(cop_new))

    # update previous copula
    cop_prev = cop_new
  }

  # AIC
  # AIC = 2*k + 2*tot_nllk

  # Values
  best_cop = sapply(cop_list, function(x) x$family)
  best_param = sapply(cop_list, function(x) x$par)
  best_param2 = sapply(cop_list, function(x) x$par2)
  best_tau = sapply(cop_list, function(x) x$tau)

  # sum
  all_logLik = sum(trans_logLik)

  if(print){
    cat("Best fitting copula(s) is/are", best_cop, "\n")
    cat("with param", best_param, "\n")
    cat("with param2", best_param2, "\n")
    cat("with tau", best_tau, "\n")
    cat("and overall transitional log likelihood", all_logLik, "\n")
  }

  if(!separate){
    trans_logLik = all_logLik
  }

  return(list(fitted_cop = cop_list, trans_logLik = trans_logLik, C2g1 = cop_prev$C2g1,
              family_vec = best_cop, par_vec = best_param, par2_vec = best_param2,
              tau_vec = best_tau))
}


#' Fit d-variate Markov order k innovation model
#'
#' @description
#' Find the best fitted Bivariate copulas for d-variate Markov order k
#' innovation model
#'
#' @param dat matrix with number of columns = d; each column corresponds to a
#'            univariate vector
#' @param k Markov order k >= 1, default is k = 1
#' @param print whether to print the result, default to FALSE
#' @param ... additional arguments pass to BiCopSelect
#' @returns best fitted families, parameter estimates, total transitional loglikelihood
#' @examples
#' #source(file.path(proj, "simulation", "sim-bivariate-mk1.R"))
#'
#' innov_dat = sim_bi_innov(n = 10000, family.innov = 1, par.innov = 0.7,
#'                          family.1 = 1, par.1 = 0.3,
#'                          family.2 = 1, par.2 = 0.6, trunc = 1000,
#'                          unif.margin = TRUE)
#'
#' testfit = fit_innov(dat = innov_dat$dat, familyset = c(0,1,4,5,14))
#'
#' test_uni = testfit$univariate_fit
#' test_innov = testfit$innov_fit
#'
#' test_innov$family
#' test_innov$par
#'
#' test_uni[[1]]$fitted_cop[[1]]$family
#' test_uni[[1]]$fitted_cop[[1]]$par
#' test_uni[[2]]$fitted_cop[[1]]$family
#' test_uni[[2]]$fitted_cop[[1]]$par
#'
#' @export
fit_innov = function(dat, k = 1, print = FALSE, ...){

  # fit
  d = ncol(dat)
  fit_all_col = lapply(1:d, function(i) fit_uni_mk(uts = dat[, i], k = k, ...))
  innov_dat = do.call(cbind, lapply(fit_all_col, function(x) x$C2g1))
  fit_innovation = VineCopula::RVineStructureSelect(data = innov_dat, ...)

  # total transitional loglikelihood
  trans_logLik = (sum(sapply(fit_all_col, function(x) x$trans_logLik)) + fit_innovation$logLik)/nrow(innov_dat)

  if(print){
    cat("The transitional loglik of the best fitting Rvine innovation model is", trans_logLik)
  }

  return(list(univariate_fit = fit_all_col, innov_fit = fit_innovation, trans_logLik = trans_logLik))
}


#' Return the loglik on a given dataset based on the univariate Markov order k model
#'
#' @description Return the out of sample transitional loglikelhood of the test
#'              sample based on the univariate Markov order k model
#' @param uni_model output returned by fit_uni_mk
#' @param uts univariate vector
#' @param separate Logical; whether log-likelihoods are returned point wisely
#'              (default: separate = FALSE)
#' @param return_innov Logical; whether return innovations a = F_{k|1,...,k-1}
#'              (default: return_innov = FALSE)
#' @returns out-of-sample transitional loglikelihood on the test set
#' @examples
#' # Load necessary functions
#' source(file.path(proj, "simulation", "sim-univariate-mkp.R"))
#'
#' ## first example for d=2, all Gumbel
#' apar1 = matrix(c(0, 2, 0,0),2,2, byrow=T)
#' apar2 = matrix(0,2,2)
#' afam4 = matrix(4,2,2) # all Gumbel
#'
#' uts_1 = sim_uni_mk(n=5000,fam=afam4,param1=apar1,param2=apar2)
#'
#' # fit
#' fit_uts = fit_uni_mk(uts_1, k = 4, familyset = c(1,4,5,14), print = TRUE)
#'
#' fit_uts$trans_logLik
#'
#' # test function
#' out_llk = UnivariateLogLik(fit_uts, uts_1, separate = FALSE)
#' out_llk
#'
#' ## first example for d=4, all Gumbel
#' apar1 = matrix(c(0, 2,2,2, 0,0, 1.5,1.5, 0,0,0,1.8, 0,0,0,0),4,4, byrow=T)
#' apar2 = matrix(0,4,4)
#' afam4 = matrix(4,4,4) # all Gumbel
#'
#' uts_3 = sim_uni_mk(n=5000,fam=afam4,param1=apar1,param2=apar2)
#'
#' # fit
#' fit_uts = fit_uni_mk(uts_3, k = 3, familyset = c(1,4,5,14), print = TRUE)
#'
#' fit_uts$trans_logLik
#'
#' # test function
#' out_llk = UnivariateLogLik(fit_uts, uts_3, return_innov = TRUE)
#' out_llk
#'
#' @export
UnivariateLogLik = function(uni_model, uts, separate = FALSE, return_innov = FALSE){

  # Markov order k
  k = length(uni_model$fitted_cop)

  # Initialization
  cdfcond1_2_prev = uts
  cdfcond2_1_prev = uts
  pdfcond2_1_prev = rep(1, times = length(uts))

  for(i in 1:k){

    # sample size
    n_i = length(cdfcond1_2_prev)

    # Conditional Distribution
    cdfcond = VineCopula::BiCopHfunc(cdfcond1_2_prev[-n_i], cdfcond2_1_prev[-1],
                                     family = uni_model$family_vec[i],
                                     par = uni_model$par_vec[i],
                                     par2 = uni_model$par2_vec[i])

    # Density
    pdfcond2_1 = VineCopula::BiCopPDF(cdfcond1_2_prev[-n_i], cdfcond2_1_prev[-1],
                                      family = uni_model$family_vec[i],
                                      par = uni_model$par_vec[i],
                                      par2 = uni_model$par2_vec[i])*(pdfcond2_1_prev[-1])

    # update
    cdfcond1_2_prev = cdfcond$hfunc2
    cdfcond2_1_prev = cdfcond$hfunc1
    pdfcond2_1_prev = pdfcond2_1
  }

  if(separate){
    trans_logLik = log(pdfcond2_1)
  }else{
    trans_logLik = sum(log(pdfcond2_1))
  }

  if(return_innov){
    return(list(trans_logLik = trans_logLik, innov = cdfcond2_1_prev))
  }else{
    return(trans_logLik)
  }
}


#' Return the loglik on a given dataset based on the Markov order k innovation model
#'
#' @description Return the out of sample transitional loglikelhood of the test
#'              sample based on the Markov order k innovation model
#' @param innov_model innovation in the form of the output of fit_innov function;
#'              list of length 3 where first gives univariate fit, second element
#'              gives innovations fit and last gives the transitional loglikelihood
#' @param umat matrix with number of columns = d; each column corresponds to a
#'            univariate vector
#' @returns out-of-sample transitional loglikelihood on umat
#' @examples
#' # example code
#' apar1x = matrix(c(0, 2,2,2, 0,0, 1.5,1.5, 0,0,0,1.2, 0,0,0,0),4,4, byrow=T)
#' apar2x = matrix(0,4,4)
#' afam4x = matrix(4,4,4) # all Gumbel
#'
#' apar1y = matrix(c(0, 1.5,1.5,1.5, 0,0, 1.7,1.7, 0,0,0,2, 0,0,0,0),4,4, byrow=T)
#' apar2y = matrix(0,4,4)
#' afam4y = matrix(4,4,4) # all Gumbel
#'
#' testing = sim_bi_innov_mkp(n = 5000, family.innov = 4, par.innov = 2,
#'                            family.1 = afam4x, par.1 = apar1x, par2.1 = apar2x,
#'                            family.2 = afam4y, par.2 = apar1y, par2.2 = apar2y,
#'                            unif.margin = TRUE)
#'
#' testfit = fit_innov(dat = testing$dat, k=3, familyset = c(0,1,4,5,14))
#' testfit$trans_logLik
#' InnovLogLik(testfit, testing$dat)
#'
#' @export
InnovLogLik = function(innov_model, umat){

  # univariate fit
  d = ncol(umat)
  all_uni = lapply(1:d, function(i) UnivariateLogLik(uni_model = innov_model$univariate_fit[[i]],
                                                          uts = umat[, i], return_innov = TRUE))

  # innovation data
  innov_dat = do.call(cbind, lapply(all_uni, function(x) x$innov))

  # loglik based on innovation
  innov_loglik = VineCopula::RVineLogLik(data = innov_dat, RVM = innov_model$innov_fit)

  # total nllk
  trans_logLik = (sum(sapply(all_uni, function(x) x$trans_logLik)) + innov_loglik$loglik)/nrow(innov_dat)

  return(trans_logLik)
}
