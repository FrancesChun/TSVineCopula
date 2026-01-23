
# library(VineCopula)


# need functions from Rvine-Rosenblatt-v2.R
# source(file.path(proj, "simulation", "Rvine-Rosenblatt-v2.R"))


##############################################################################


#' Simulate d-variate Markov order k time series based on Innovation model
#' @description
#' simulate d-variate time series based on innovations
#' NOT initializing from a stationary state
#' @param n number of bivariate observations simulated
#' @param innov.RVine  RVineMatrix() object for the innovations
#' @param uni.Rvine a list of d RVineMatrix() objects for each of the d univariate
#'                  time series; d>=2
#' @param trunc numeric; number of sample points (starting from the beginning)
#'              to remove, default 0 (i.e. keep all data points)
#' @param thin numeric; keep every thin(th) points in the dataset;
#'             Ex: if thin = 2, then we keep 1st, 3rd, 5th, 7th ... data points
#'             default 1 (i.e. keep all data points)
#' @param seed random number seed, default no seed
#' @returns simulated time series (nxd matrix)
#' @examples
#' var1apar1 = matrix(c(0, 2,2,2, 0,0, 1.5,1.5, 0,0,0,1.2, 0,0,0,0),4,4, byrow=T)
#' var1apar2 = matrix(0,4,4)
#' var1afam4 = matrix(4,4,4) # all Gumbel
#'
#' var2apar1 = matrix(c(0, 1.5,1.5,1.5, 0,0, 1.7,1.7, 0,0,0,2, 0,0,0,0),4,4, byrow=T)
#' var2apar2 = matrix(0,4,4)
#' var2afam4 = matrix(4,4,4) # all Gumbel
#'
#' var3apar1 = matrix(c(0, 1.6,1.6,1.6, 0,0, 1.8,1.8, 0,0,0,2, 0,0,0,0),4,4, byrow=T)
#' var3apar2 = matrix(0,4,4)
#' var3afam4 = matrix(4,4,4) # all Gumbel
#'
#' innovama   = matrix(c(1, 1, 1,
#'                       0, 2, 2,
#'                       0, 0, 3), 3, 3, byrow = T)
#' innovapar1 = matrix(c(0, 0.5,0.5, 0,0, 0.3, 0,0,0),3,3, byrow=T)
#' innovapar2 = matrix(0,3,3)
#' innovafam3 = matrix(1,3,3) # all Gaussian
#'
#' testing = sim_innov(n = 5000, innov.RVine = RVineMatrix(Matrix = innovama, family = innovafam3,
#'                                                         par = innovapar1, par2 = innovapar2),
#'                     uni.Rvine = list(var1_RVM = RVineMatrix(Matrix = dvine(4), family = var1afam4,
#'                                                             par = var1apar1, par2 = var1apar2),
#'                                      var2_RVM = RVineMatrix(Matrix = dvine(4), family = var2afam4,
#'                                                             par = var2apar1, par2 = var2apar2),
#'                                      var3_RVM = RVineMatrix(Matrix = dvine(4), family = var3afam4,
#'                                                             par = var3apar1, par2 = var3apar2)))
#'
#' testing_fit = fit_innov(dat = testing$dat, k = 3, familyset = c(0,1,4,5,14))
#'
#' test_uni = testing_fit$univariate_fit
#' test_innov = testing_fit$innov_fit
#'
#' test_innov$Matrix
#' test_innov$family
#' test_innov$par
#'
#' test_uni[[1]]$fitted_cop[[1]]$family
#' test_uni[[1]]$fitted_cop[[1]]$par
#' test_uni[[2]]$fitted_cop[[1]]$family
#' test_uni[[2]]$fitted_cop[[1]]$par
#'
#' @export
sim_innov = function(n, innov.RVine, uni.Rvine,
                        trunc = 0, thin = 1, seed=FALSE){
  if(is.numeric(seed)){
    set.seed(seed)
  }

  # simulate the innovations (note: the first k innovations won't be used)
  innov = VineCopula::RVineSim(n, innov.RVine)

  # d-variate
  d = length(uni.Rvine)

  # pair copula family matrix of the first time series
  family_1 = uni.Rvine[[1]]$family

  # Markov order k
  if(is.matrix(family_1)){
    k = nrow(family_1)-1
  }else{
    k=1
  }

  # define the output matrix
  umat = matrix(nrow = n, ncol = d)

  # initialize the matrix
  umat[1:k,] = runif((k*d))

  # simulate the time series based on innovations
  for(i in (k+1):n){

    # Markov order 1
    if(k==1){

      for(j in 1:d){
        currentRvine = uni.Rvine[[j]]
        umat[i,j] = VineCopula::BiCopHinv1(u1 = umat[i-1,j], u2 = innov[i, j],
                                           family = to_upper_tri(currentRvine$family)[1,2],
                                           par = to_upper_tri(currentRvine$par)[1,2],
                                           par2 = to_upper_tri(currentRvine$par2)[1,2])
      }

    # Markov order 2 or higher
    }else{
      # use Rosenblatt transform to compute
      # p_{i-k+1}=C_{2|1}(u_{i-k+1}|u_{i-k})
      # ...
      # p_{i-1}=C_{k|1,...,k-1}(u_{i-1}|u_{i-k},...,u_{i-2})

      ulen = i-1 # length of current simulated points

      for(j in 1:d){
        currentRvine = uni.Rvine[[j]]

        uvec_prev = umat[(ulen-k+1):ulen,j]
        lenma = nrow(currentRvine$Matrix)
        pvec_prev = rvinepcond(uvec_prev,
                               to_upper_tri(currentRvine$Matrix)[-lenma,-lenma],
                               to_upper_tri(currentRvine$family)[-lenma,-lenma],
                               to_upper_tri(currentRvine$par)[-lenma,-lenma],
                               to_upper_tri(currentRvine$par2)[-lenma,-lenma])
        pvec_i = c(pvec_prev, innov[i,j])

        # gets u values based on pvec_i via inverse Rosenblatt transform (rvineqcond)
        uvec_i = rvineqcond(pvec_i, currentRvine$Matrix, currentRvine$family,
                            currentRvine$par, currentRvine$par2)

        umat[i,j] = uvec_i[k+1]
      }
    }

  }

  # truncation and thinning
  trunc.thin.index = seq(from = (trunc+1), to = n, by = thin)
  innov = innov[trunc.thin.index,]
  umat = umat[trunc.thin.index,]

  return(list(dat = umat, innov = innov))
}
