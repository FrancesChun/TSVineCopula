
####################################
#                                  #
#  Functions for fitting S-vine    #
#                                  #
####################################

#' Convert lower triangular matrix to upper triangular matrix
#'
#' @description
#' Convert lower triangular matrix to upper triangular matrix and convert all
#' NA entries to 0 entries
#'
#' @param mat triangular matrix
#' @param NA_to_0 Boolean, whether to convert NA to 0; default is TRUE
#' @returns Upper triangular matrix
#' @export
to_upper_tri <- function(mat, NA_to_0 = TRUE){

  outmat = mat
  if(NA_to_0){
    outmat[is.na(outmat)] <- 0
  }

  if(fastmatrix::is.lower.tri(outmat)){
    outmat = rotate_180(outmat)
  }

  return(outmat)
}


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


#' Return averaged par matrix
#'
#' @param par par matrix with inconsistent entries
#' @param d Markov order k >= 1, default is k = 1
#' @returns Averaged par matrix
#'
out_svine_par = function(par, d=2){

  mat = par
  n = nrow(par)

  average_spaced_all_starts <- function(vec, d) {
    m <- length(vec)
    newvec <- vec

    if(m > d){
      for (start in 1:d) {
        indices <- seq(start, m, by = d)
        mean_val <- mean(vec[indices])
        newvec[indices] <- mean_val
      }
    }

    return(newvec)
  }

  if(fastmatrix::is.lower.tri(par)){

    for(i in 2:n){
      mat[i, 1:(i-1)] = average_spaced_all_starts(par[i, 1:(i-1)], d)
    }

  }else if(fastmatrix::is.upper.tri(par)){

    for(i in 1:(n-1)){
      mat[i, (i+1):n] = average_spaced_all_starts(par[i, (i+1):n], d)
    }

  }

  return(mat)
}


#' Return Consistent Model
#'
#' @description Return Consistent Model
#' @param Smodel list of length 4, must contain (Matrix, family, par, par2)
#' @param d number of variables, default is d = 2
#' @returns Consistent Model (list of length 4)
#'
return_consistent_model = function(Smodel, d=2){

  # check consistent family
  avg_family = out_svine_par(Smodel$family, d)
  if(!isTRUE(all.equal(avg_family, Smodel$family))){
    stop(paste("Error: Family matrix is not consistent.",
               "Family matrix is", Smodel$family))
  }

  Smodel_new = Smodel
  Smodel_new$par = out_svine_par(Smodel$par, d)
  Smodel_new$par2 = out_svine_par(Smodel$par2, d)

  return(Smodel_new)
}





#' Calculates the log-likelihood of a S-vine copula model
#'
#' @description Calculates the log-likelihood of a S-vine copula model for a
#'              given copula data set.
#' @param data matrix with number of columns = d; each column corresponds to a
#'            univariate vector with uniform margin
#' @param Smodel list of at least length 5, must contain (Matrix, family, par, par2, var_order)
#' @returns log-likelihood
#' @examples
#' # example for d=4, all Gaussian
#' d2_mk1 = matrix(c(1,1,1,3,
#'                   0,2,2,1,
#'                   0,0,3,2,
#'                   0,0,0,4),
#'                 nrow = 4, byrow = TRUE)
#' apar1 = matrix(c(0, 0.4, 0.5, 0.4,
#'                  0, 0, 0.2, 0.3,
#'                  0, 0, 0, 0.2,
#'                  0, 0, 0, 0),
#'                4,4, byrow=T)
#' apar2 = matrix(0,4,4)
#' afam4 = matrix(1,4,4) # all Gaussian
#'
#' # simulation
#' set.seed(123)
#' mvine_dat = sim_svine(n = 5000, vinematrix = d2_mk1, fam = afam4,
#' #'                       param1 = apar1, param2 = apar2)
#'
#' # fit
#' testfit = SVineCopSelect(data = mvine_dat, Matrix = d2_mk1, var_order = c(1,2),
#'                          k = 1, familyset = c(0,1,4,5,14))
#' summary(testfit)
#' contour(testfit)
#'
#' testfit$Matrix
#' testfit$family
#' testfit$par
#' testfit$trans_logLik #0.2896287
#'
#' @export
SVineLogLik = function(data, Smodel){

  # number of variables
  d = ncol(data)

  # Markov order k
  N = nrow(Smodel$Matrix)
  k = (N/d)-1

  # Prepare the Dataset
  df = extend_df(dat = data, var_order = Smodel$var_order, k = k)

  # compute Markov order k-1 model
  Smodel_sub = VineCopula::RVineMatrix(Matrix = to_upper_tri(Smodel$Matrix)[1:(d*k),1:(d*k)],
                           family = to_upper_tri(Smodel$family)[1:(d*k),1:(d*k)],
                           par = to_upper_tri(Smodel$par)[1:(d*k),1:(d*k)],
                           par2 = to_upper_tri(Smodel$par2)[1:(d*k),1:(d*k)])

  # transition probability f_{x_{t+k}|x_{t+k-1},...,x_{t}}
  df_sub = df[,1:(d*k)]
  logLik_full = VineCopula::RVineLogLik(df, RVM = Smodel)
  logLik_sub = VineCopula::RVineLogLik(df_sub, RVM = Smodel_sub)
  trans_logLik = (logLik_full$loglik - logLik_sub$loglik)/nrow(df)

  return(trans_logLik)
}


#' Fit a S-vine model
#'
#' @description Return the fitted S-vine family and parameters, along with
#'             transitional loglikelihood on the dataset
#' @param data matrix with number of columns = d; each column corresponds to a
#'            univariate vector with uniform margin
#' @param Matrix can be either:
#'              1. lower or upper triangular 2d x 2d matrix that defines the
#'                 Markov order 1 S-vine tree structure.
#'              2. lower or upper triangular d(k+1) x d(k+1) matrix that defines the
#'                 Markov order k S-vine tree structure.
#' @param k Markov order k >= 1, default is k = 1
#' @param var_order A permutation of numbers from 1 to d, re-arranges column in
#'              data; default is 1:d (no re-arrangement needed)
#' @param ... additional arguments pass to RVineCopSelect
#' @returns best fitted families, parameter estimates, loglikelihood
#' @export
SVineCopSelect = function(data, Matrix, k = 1, var_order = NULL, ...){

  # number of variables
  d = ncol(data)

  # make sure Matrix is upper triangular
  Matrix = to_upper_tri(Matrix)

  # Markov order k S-vine array
  if(nrow(Matrix) == (d*(k+1)) && ncol(Matrix) == nrow(Matrix)){
    mk_d_mat = Matrix
  }else{
    mk_d_mat = mk1_to_mkp(Matrix, p=k)
  }
  mk_d_mat[is.na(mk_d_mat)] <- 0

  # Prepare the Dataset
  df = extend_df(dat = data, var_order = var_order, k = k)

  # Fit the model
  fitRvine_raw = VineCopula::RVineCopSelect(data = df, Matrix = mk_d_mat, ...)

  # take mean par
  fitRvine = return_consistent_model(fitRvine_raw, d=d)

  fitRvine = VineCopula::RVineMatrix(Matrix = mk_d_mat,
                         family = fitRvine$family,
                         par = fitRvine$par,
                         par2 = fitRvine$par2)

  if(is.null(var_order)){
    fitRvine$var_order = 1:d
  }else{
    fitRvine$var_order = var_order
  }

  fitRvine$trans_logLik = SVineLogLik(data = data, Smodel = fitRvine)

  return(fitRvine)
}





#' Fit S-vine d-variate Markov order k model
#'
#' @description
#' Find the best fitted bivariate copulas for Markov order k, d-variate model
#'
#' @param dat matrix with number of columns = d; each column corresponds to a
#'            univariate vector
#' @param k Markov order k >= 1, default is k = 1
#' @param train_prop Proportion of datasets in the training set, default is
#'          train = 1 (i.e. no test set)
#' @param ... additional arguments pass to RVineCopSelect
#' @returns best fitted families, parameter estimates, negative loglikelihood
#' @examples
#'
#' # example for d=4, all Gaussian
#' d2_mk1 = matrix(c(1,1,1,3,
#'                   0,2,2,1,
#'                   0,0,3,2,
#'                   0,0,0,4),
#'                   nrow = 4, byrow = TRUE)
#' apar1 = matrix(c(0, 0.4, 0.5, 0.4,
#'                  0, 0, 0.2, 0.3,
#'                  0, 0, 0, 0.2,
#'                  0, 0, 0, 0),
#'                4,4, byrow=T)
#' apar2 = matrix(0,4,4)
#' afam4 = matrix(1,4,4) # all Gaussian
#'
#' # simulation
#' mvine_dat = sim_svine(n = 5000, vinematrix = d2_mk1, fam = afam4,
#'                       param1 = apar1, param2 = apar2)
#'
#' testfit = fit_svine(dat = mvine_dat, train_prop = 0.8, familyset = c(0,1,4,5,14))
#' testfit[[1]]$Matrix
#' testfit[[1]]$family
#' testfit[[1]]$trans_logLik
#'
#' @export
fit_svine = function(dat, k = 1, train_prop = 1, Matrix = NULL, ...){

  # dimensions
  n=nrow(dat)
  d=ncol(dat)

  # train split
  if(train_prop < 1){
    ntrain = floor(train_prop*n)
    train = dat[1:ntrain,]
    test = dat[(ntrain+1-k):n,]
  }

  if(is.null(Matrix)){
    # possible d-variate S-vine array
    d_matrix_list = svine_mk1(d)
  }else{
    d_matrix_list = list(Matrix)
  }

  # possible permutations
  var_order = combinat::permn(1:d)

  # output
  result = list()

  # looping through each possible S-vine array
  for(mat in d_matrix_list){

    for(i in var_order){
      if(train_prop == 1){
        fitRvine = SVineCopSelect(data = dat, Matrix = mat, k = k, var_order = i, ...)
      }else if(train_prop < 1){
        fitRvine = SVineCopSelect(data = train, Matrix = mat, k = k, var_order = i, ...)
        fitRvine$test_trans_logLik = SVineLogLik(data = test, Smodel = fitRvine)
      }
      result = append(result, list(fitRvine))
    }

  }

  return(result)
}

