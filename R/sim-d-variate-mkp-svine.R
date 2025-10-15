
#library(VineCopula)

# import functions
# source(file.path(proj, "simulation", "Rvine-Rosenblatt-v2.R"))
# source(file.path(proj, "simulation", "mvine-array.R"))

##############################################################################


#' Simulate S-vine d-variate Markov order k-1 time series
#' @param n length of time series x and y.
#' @param vinematrix d-variate Markov order k S-vine matrix, d(k+1)xd(k+1)
#' @param fam d(k+1)xd(k+1)  family code matrix (VineCopula copula family index) .
#'            where the diagonal is in the order of
#'            x1 - y1 - x2 - y2 - ... - xk, yk
#'            Note: we assume CX1Y1 = CX2Y2 --> fam[1,2] = fam[1,4]
#' @param param1 d(k+1)xd(k+1)  parameter matrix (VineCopula first copula parameter)
#'     param1[l,j] for tree l, variable j.
#'     Note: we assume CX1Y1 = CX2Y2 --> param1[1,2] = param1[1,4]
#' @param param2 d(k+1)xd(k+1)  parameter matrix (VineCopula second copula parameter)
#'     param2[l,j] for tree l, variable j; 0 if not needed.
#' @param d numeric, number of variables, must be >=2, default is 2
#' @param trunc numeric; number of sample points (starting from the beginning)
#'              to remove, default 0 (i.e. keep all data points)
#' @param thin numeric; keep every thin(th) points in the dataset;
#'             Ex: if thin = 2, then we keep 1st, 3rd, 5th, 7th ... data points
#'             default 1 (i.e. keep all data points)
#' @param seed random number seed, default no seed
#' @returns bivariate time series of length n
#' @examples
#' ## Example 1: 4-variate Markov order 1, all Gumbel
#' d2_mk1 = matrix(c(1,1,1,3,
#'                   0,2,2,1,
#'                   0,0,3,2,
#'                   0,0,0,4),
#'                 nrow = 4, byrow = TRUE)
#' apar1 = matrix(c(0,2,2.1,2,
#'                  0,0,1.5,1.4,
#'                  0,0,0,1.2,
#'                  0,0,0,0),4,4, byrow=T)
#' apar2 = matrix(0,4,4)
#' afam4 = matrix(4,4,4) # all Gumbel
#'
#' mvine_dat = sim_svine(5000,d2_mk1,afam4,apar1,apar2, seed = 123)
#'
#' testfit = fit_svine(dat = mvine_dat, family, familyset = c(0,1,4,5,14))
#' rotate_180(testfit[[1]]$Matrix)
#' rotate_180(testfit[[1]]$family)
#' rotate_180(testfit[[1]]$par)
#'
#' ## Example 2: Bivariate Markov order 2, all Gumbel
#' d2_mk2 = bimvine_mk(2)
#' apar1 = matrix(c(0, 2, 2.5, 2,   2.5, 2,
#'                  0, 0, 1.5, 1.3, 1.5, 1.3,
#'                  0, 0, 0,   1.2, 1.6, 1.2,
#'                  0, 0, 0,   0,   1.1, 1.1,
#'                  0, 0, 0,   0,   0,   1.1,
#'                  0, 0, 0,   0,   0,   0),
#'                6,6, byrow=T)
#' apar2 = matrix(0,6,6)
#' afam4 = matrix(4,6,6) # all Gumbel
#'
#' # default d=2
#' mvine_dat = sim_svine(5000,d2_mk2,afam4,apar1,apar2, seed = 123)
#'
#' testfit = fit_svine(dat = mvine_dat, k=2, familyset = c(0,1,4,5,14))
#'
#' @export
sim_svine = function(n, vinematrix, fam, param1, param2, d=2,
                          seed = FALSE, ntrunc=0,varname=numeric(0),iprint=FALSE,
                          trunc = 0, thin = 1){

  # convert matrix to upper triangular matrix
  vinematrix = to_upper_tri(vinematrix)
  fam = to_upper_tri(fam)
  param1 = to_upper_tri(param1)
  param2 = to_upper_tri(param2)

  if(is.numeric(seed)){
    set.seed(seed)
  }

  # Markov order k
  k = nrow(vinematrix)/d-1

  # generate U(0,1)
  p = runif(d*n)

  # initialize the matrix
  umat = matrix(nrow = n, ncol = d)

  # for each iteration
  for(i in (k+1):n){

    # pvec_i = p vector in ith iteration
    # if i == k+1, then initialize pvec_i
    if(i == (k+1)){
      pvec_i = p[1:(d*(k+1))]

      # else, compute pvec_i
    }else{
      uvec_prev = uvec_i[(d+1):(d*(k+1))]

      lastd = (d*k+1):(d*(k+1))
      pvec_prev = rvinepcond(uvec_prev,vinematrix[-lastd,-lastd],
                             fam[-lastd,-lastd],
                             param1[-lastd,-lastd],param2[-lastd,-lastd])

      plen = d*i
      pvec_i = c(pvec_prev, p[(plen-d+1):plen])
    }

    # gets u values based on pvec_i via inverse Rosenblatt transform (rvineqcond)
    uvec_i = rvineqcond(pvec_i,vinematrix,fam,param1,param2,
                        ntrunc,varname,iprint)

    if(i == (k+1)){
      # initialize u vector
      umat[1:(k+1),] = matrix(uvec_i, nrow = (k+1), ncol = d, byrow = TRUE)
    }else{
      umat[i,] = uvec_i[(d*k+1):(d*(k+1))]
    }
  }

  # truncation and thinning
  trunc.thin.index = seq(from = (trunc+1), to = n, by = thin)
  umat = umat[trunc.thin.index,]

  return(umat)
}

