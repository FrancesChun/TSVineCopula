
####################################
#                                  #
#  Functions for fitting S-vine    #
#                                  #
####################################


#' Helper function for fit_svine
#'
#' @param dat N x d matrix; each column corresponds to a
#'            univariate vector
#' @param k Markov order k >= 1, default is k = 1
#' @param var_order a vector that contains the desired order of the variable
#'        (If var_order = 1:d, it means the order of the variable stays the same
#'        as the original order in dat)
#' @returns an extended (N−k)×d(k+1) matrix with variables in the order specified
#'         by var_order and Markov order k
#' @examples
#' test_d3 = matrix(data = runif(30), ncol = 3, nrow = 10)
#' extend_df(test_d3, var_order = 1:3)
#'
extend_df = function(dat, var_order, k=1){

  # number of rows
  N = nrow(dat)

  # Reorder the columns
  dat_reordered <- dat[, var_order]

  # Extend the matrix
  dat_ext = dat_reordered[1:(N-k), ]
  for(i in 1:k){
    dat_ext = cbind(dat_ext, dat_reordered[(1+i):(N-k+i), ])
  }

  return(dat_ext)
}


#' Fit S-vine d-variate Markov order k model
#'
#' @description
#' Find the best fitted bivariate copulas for Markov order k, d-variate model
#'
#' @param dat matrix with number of columns = d; each column corresponds to a
#'            univariate vector
#' @param k Markov order k >= 1, default is k = 1
#' @param ... additional arguments pass to RVineCopSelect
#' @returns best fitted families, parameter estimates, negative loglikelihood
#' @examples
#' # need to import sim_bi_mvine function from sim-bivariate-mk1.R
#' source(file.path(proj, "simulation", "sim-bivariate-mk1.R"))
#'
#' # example for d=4, all Gaussian
#' apar1 = matrix(c(0, 0.4, 0.5, 0.4,
#'                  0, 0, 0.2, 0.3,
#'                  0, 0, 0, 0.2,
#'                  0, 0, 0, 0),
#'                4,4, byrow=T)
#' apar2 = matrix(0,4,4)
#' afam4 = matrix(1,4,4) # all Gaussian
#'
#' # simulation
#' mvine_dat = sim_bi_mvine(5000,afam4,apar1,apar2)
#'
#' testfit = fit_svine(dat = mvine_dat)
#' testfit[[1]]$Matrix
#' testfit[[1]]$family
#' testfit[[1]]$trans_logLik
#'
#' @export
fit_svine = function(dat, k = 1, ...){

  # dimensions
  n=nrow(dat)
  d=ncol(dat)

  # possible d-variate S-vine array
  d_matrix_list = svine_mk1(d)

  # possible permutations
  var_order = combinat::permn(1:d)

  # output
  result = list()

  # looping through each possible S-vine array
  for(mat in d_matrix_list){

    mk_d_mat = mk1_to_mkp(mat, p=k)

    # Replace NA with 0
    mk_d_mat[is.na(mk_d_mat)] <- 0

    for(i in var_order){
      df = extend_df(dat, var_order = i, k = k)

      # fits R-vine copula model
      fitRvine = VineCopula::RVineCopSelect(data = df, Matrix = mk_d_mat, ...)
      fitRvine$var_order = i

      # compute Markov order k-1 model
      df_sub = df[,1:(d*k)]
      mat_sub = mk_d_mat[1:(d*k),1:(d*k)]
      fitRvine_sub = VineCopula::RVineCopSelect(data = df_sub, Matrix = mat_sub, ...)

      # transition probability f_{x_{t+k}|x_{t+k-1},...,x_{t}}
      fitRvine$trans_logLik = fitRvine$logLik - fitRvine_sub$logLik
      result = append(result, list(fitRvine))
    }

  }


  return(result)
}


