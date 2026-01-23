# This document contains the following functions (excluding helper functions):

# dvine :
#   Generate A D-vine array with 1:d on diagonal

# bimvine_mk :
#   Generate A Bivariate M-vine array with
#   X_1, Y_1, ..., X_{k+1}, Y_{k+1} on diagonal
#
# mk1_to_mkp :
#   Extend the given Markov order 1 M-vine array
#   to Markov order p
#
# SVineMatrix_mk1_Check:
#   Testing whether the given upper triangular matrix is a valid
#   Markov order 1 S-vine matrix.
#
# svine_mk1:
#   Generate all possible d-dimensional Markov order 1 S-vine array
#   with 1=X_1, 2=X_2,...,d=X_d, d+1=X_1',...,2d=X_d' on diagonal
#
# svine_mkp:
#   Generate all possible d-dimensional Markov order p S-vine array
#   with 1=X_1, 2=X_2,...,d=X_d, d+1=X_1',...,2d=X_d' on diagonal
#

#============================================================

# library(VineCopula)
# library(combinat)

###################################################################
#                                                                 #
#         d-dimensional Markov order p M-vine array Matrix        #
#                                                                 #
###################################################################



#' Generate A D-vine array with 1:d on diagonal
#'
#' @param d Number of variables, must be an integer
#' @returns A D-vine array with 1:d on diagonal
#' @export
dvine <- function(d){
  dvinearray = matrix(0, nrow = d, ncol = d)
  for(i in 1:d){
    # if last row
    if(i == d){
      dvinearray[i,i] <- i
    }else{
      # otherwise
      dvinearray[i,i:d] <- c(i,1:(d-i))
    }
  }
  return(dvinearray)
}

#' Rotate a matrix by 180 degree
#'
#' @param mat Matrix
#' @returns Same matrix except that it has been rotated by 180 degree
#' @export
rotate_180 <- function(mat){
  mat[nrow(mat):1, ncol(mat):1]
}

#' Generate A Bivariate M-vine array
#'
#' Generate A Bivariate M-vine array with X_1, Y_1, ..., X_{k+1}, Y_{k+1}
#' on diagonal
#'
#' @param k Markov order k >= 1, default is k = 1, must be an integer
#' @returns A Bivariate M-vine array
#' @examples
#' # return the R-vine array for Bivariate M-vine Markov order 2
#' bimvine_mk(2)
#'
#' @export
bimvine_mk <- function(k){

  # length of X_1, Y_1, ..., X_{k+1}, Y_{k+1}
  d = 2*(k+1)

  # initialize the matrix
  mvinearray = matrix(0, nrow = d, ncol = d)

  # sequentially add values to each diagonal
  for(i in 1:d){
    if(i == 1){
      diag(mvinearray) = 1:d
    }else if(i %% 2 == 0){ #even diagonal
      mvinearray[row(mvinearray) + i - 1 == col(mvinearray)] = i
    }else if(i %% 2 == 1){ #odd diagonal
      mvinearray[row(mvinearray) + i - 1 == col(mvinearray)] = i-2
    }
  }
  # first row and even column:
  mvinearray[row(mvinearray) == 1 & (col(mvinearray) %% 2 == 0)] = seq(from=1, to=d, by = 2)

  return(mvinearray)
}


#' Extend the given Markov order 1 M-vine array to Markov order p
#'
#' @description
#' Extend the given Markov order 1 M-vine array to Markov order p
#'
#' @param mk1_mvine_array matrix, Markov order 1 M-vine array, ncol = nrow and is even
#' @param p integer >= 1, Markov order, default 2
#' @returns A Markov order p M-vine array
#' @examples
#' d2_mk1 = matrix(c(1,1,1,3,
#'                   NA,2,2,1,
#'                   NA,NA,3,2,
#'                   NA,NA,NA,4),
#'                 nrow = 4, byrow = TRUE)
#' mk1_to_mkp(d2_mk1,p=2)  #Markov order 2
#' bimvine_mk(2)
#' mk1_to_mkp(d2_mk1,p=3)  #Markov order 3
#'
#' d3_mk1 = matrix(c(1,1,1,1,4,4,
#'                   NA,2,2,2,1,5,
#'                   NA,NA,3,3,2,1,
#'                   NA,NA,NA,4,3,2,
#'                   NA,NA,NA,NA,5,3,
#'                   NA,NA,NA,NA,NA,6),
#'                nrow = 6, byrow = TRUE)
#' mk1_to_mkp(d3_mk1,p=2) #Markov order 2
#' mk1_to_mkp(d3_mk1,p=3) #Markov order 3
#'
#' @export
mk1_to_mkp <- function(mk1_mvine_array, p=2){

  # dimension
  d = nrow(mk1_mvine_array)/2

  # initialize the vine-array
  mkp_array = matrix(0, nrow = d*(p+1), ncol = d*(p+1))
  mkp_array[1:(2*d), 1:(2*d)] = mk1_mvine_array

  if(p >= 2){

    # add diagonal elements
    diag(mkp_array) = 1:nrow(mkp_array)

    # for each time point
    for(t in 3:(p+1)){

      # First, copy the rows and columns
      copy_col_ind = ((t-2)*d+1):((t-1)*d)
      copy_row_ind = 1:(d*(t-2))
      mkp_array[copy_row_ind, (copy_col_ind + d)] = mkp_array[copy_row_ind, copy_col_ind]+d

      # Second, extend diagonally
      for(i in (copy_col_ind + d)){
        for(j in (d*(t-2)+1):(i-1)){
          mkp_array[j,i] = mkp_array[(j-1),(i-1)]
        }
      }
    }
  }

  # convert NA to 0
  mkp_array[is.na(mkp_array)] <- 0

  return(mkp_array)
}


#' Helper function, return corresponding index
#'
#' @param diag_ele  Diagonal elements of a Markov order 1 vine array. A 2d vector.
#' @param index     The index for which we want to get the corresponding index,
#'                  A vector
#' @param is.next   Boolean. If TRUE, returns the next corresponding index.
#'                  Otherwise return the previous index. Default FALSE.
#' @returns         Corresponding index, a vector with the same length as index.
#' @examples
#' A <- c(3, 2, 1, 6, 5, 4)
#' B <- c(4, 5, 4)
#' corr_index(diag_ele=A, index=B)
#'
corr_index <- function(diag_ele, index, is.next = FALSE) {

  # number of dimensions
  d = length(diag_ele)/2

  # return error if d is not an integer
  if (d %% 1 != 0) {
    # Stop execution and return an error
    stop("Error: Length of diag_ele must be even.")
  }

  # return error if index not in diag_ele
  if (!all(index %in% diag_ele)) {
    # Stop execution and return an error
    stop("Error: Index does not exist in diag_ele.")
  }

  diag_1 = diag_ele[1:d]
  diag_2 = diag_ele[(d+1):(2*d)]

  if(is.next){
    out_index = diag_2[match(index, diag_1)]
  }else{
    out_index = diag_1[match(index, diag_2)]
  }

  # # return error if index not in diag_1/diag_2
  # if (anyNA(out_index)) {
  #   # Stop execution and return an error
  #   stop("Error: No corresponding index exists. Consider changing is.next.")
  # }

  return(out_index)
}


#' Check whether the matrix is a valid Markov order 1 S-vine matrix
#'
#' Testing whether the given upper triangular matrix is a valid
#' Markov order 1 S-vine matrix.
#' Requires function RVineMatrixCheck from package VineCopula.
#'
#' @param M     A 2dx2d vine matrix.
#' @returns code  1 for OK;
#'              -5 matrix is a valid RvineMatrix but not a valid Svine;
#'              -4 matrix is neither lower nor upper triangular;
#'              -3 diagonal can not be put in order d:1;
#'              -2 for not permutation of j:d in column d-j;
#'              -1 if cannot find proper binary array from array in natural order.
#' @examples
#' d2_mk1 = matrix(c(1,1,1,3,
#'                   NA,2,2,1,
#'                   NA,NA,3,2,
#'                   NA,NA,NA,4),
#'                 nrow = 4, byrow = TRUE)
#' SVineMatrix_mk1_Check(d2_mk1) # should be 1
#'
#' d2_mk1_F = matrix(c(1,1,1,2,
#'                     NA,2,2,1,
#'                     NA,NA,3,3,
#'                     NA,NA,NA,4),
#'                   nrow = 4, byrow = TRUE)
#' SVineMatrix_mk1_Check(d2_mk1_F) # should NOT be 1
#'
#' d3_mk1 = matrix(c(1,1,1,1,4,4,
#'                   NA,2,2,2,1,5,
#'                   NA,NA,3,3,2,1,
#'                   NA,NA,NA,4,3,2,
#'                   NA,NA,NA,NA,5,3,
#'                   NA,NA,NA,NA,NA,6),
#'                 nrow = 6, byrow = TRUE)
#' SVineMatrix_mk1_Check(d3_mk1)  # should be 1
#'
#' d3_mk1_F = matrix(c(1,1,1,1,4,4,
#'                     NA,2,2,2,1,5,
#'                     NA,NA,3,3,3,1,
#'                     NA,NA,NA,4,2,2,
#'                     NA,NA,NA,NA,5,3,
#'                     NA,NA,NA,NA,NA,6),
#'                   nrow = 6, byrow = TRUE)
#' SVineMatrix_mk1_Check(d3_mk1_F)  # should NOT be 1
#'
#' @export
SVineMatrix_mk1_Check <- function(M){

  # number of dimensions
  d = nrow(M)/2

  # return error if d is not an integer
  if (d %% 1 != 0) {
    # Stop execution and return an error
    stop("Error: M must be in dimension 2dx2d.")
  }

  # Replace NA values with 0
  M[is.na(M)] <- 0

  # RVineMatrixCheck
  code <- VineCopula::RVineMatrixCheck(M)

  # if it is a valid Rvine Matrix, check if it is a valid Svine:
  if(code == 1){

    # check linking variable
    if(M[1,d+1]!=M[1,1]){
      code = -5
      return(code)
    }

    ### check other variables ###

    # extracting fixed elements
    mat_l = M[1:d,1:d]
    mat_r = M[1:d,(d+1):(2*d)]
    mat_l_ele = mat_l[upper.tri(mat_l, diag = FALSE)]
    mat_r_ele = mat_r[upper.tri(mat_r, diag = FALSE)]

    # check whether they are +d of each other
    corr_ind = corr_index(diag_ele = diag(M), index = mat_r_ele)
    if(anyNA(corr_ind) | any(corr_ind != mat_l_ele)){
      code = -5
    }

  }

  return(code)
}


#' Helper function 1
#'
#' Helper function, attempt to fill in all the cells in a given column
#' requires permn from Library combinat
#'
#' @param mat     A 2dx2d vine matrix.
#' @param k       kth Column, an integer.
#' @returns        List of possible matrices.
#' @examples
#'test <- matrix(c(1, 1, NA, NA,
#'                 NA, 2, NA, NA,
#'                 NA, NA, 3, NA,
#'                 NA, NA, NA, 4), nrow = 4, ncol = 4, byrow = TRUE)
#'all_possible_matrices <- generate_specific_upper_tri_matrices(test, 3)
#'print(all_possible_matrices)
#'
generate_specific_upper_tri_matrices <- function(mat, k) {
  # Validate column index
  if (k < 1 || k > ncol(mat)) {
    stop("Error: Column index k is out of bounds")
  }

  # Find the positions of NA elements in the upper triangular part of the k-th column
  na_positions = which(is.na(mat) & upper.tri(mat) & col(mat) == k, arr.ind = TRUE)

  # Generate all combinations of possible values for the NA positions
  possible_values = setdiff(1:k, mat[,k])
  if(length(possible_values)==1){
    combinations <- matrix(data = possible_values)
  }else{
    combinations <- do.call(rbind, combinat::permn(possible_values))
  }

  # Create a list to store all possible matrices
  all_matrices <- list()

  # Fill the NA positions with all possible combinations
  for (i in 1:nrow(combinations)) {
    new_mat <- mat
    new_mat[na_positions[,1], k] <- combinations[i,]
    all_matrices = append(all_matrices, list(new_mat))
  }

  return(all_matrices)
}


#' Helper function 2
#'
#' Helper function, extend column d+2:2d and corresponding rows
#'
#' @param mat     A 2dx2d vine matrix where the upper 1:d+1 columns have been filled.
#' @returns        A filled 2dx2d vine matrix.
#' @examples
#'d2_mk1 <- matrix(c(1, 1, 1, NA,
#'                   NA, 2, 2, NA,
#'                   NA, NA, 3, NA,
#'                   NA, NA, NA, 4), nrow = 4, ncol = 4, byrow = TRUE)
#'extend_to_2d(d2_mk1)
#'
#'d3_mk1 = matrix(c(1, 1, 1, 1, NA,NA,
#'                  NA,2, 2, 2, NA,NA,
#'                  NA,NA,3, 3, NA,NA,
#'                  NA,NA,NA,4, NA,NA,
#'                  NA,NA,NA,NA,5, NA,
#'                  NA,NA,NA,NA,NA,6),
#'                nrow = 6, byrow = TRUE)
#'extend_to_2d(d3_mk1)
#'
extend_to_2d <- function(mat) {

  # number of dimensions
  d = nrow(mat)/2

  # return error if d is not an integer
  if (d %% 1 != 0) {
    # Stop execution and return an error
    stop("Error: mat must be in dimension 2dx2d.")
  }

  # initialize matrix
  filled_mat = mat

  # fill in the d+2:2d columns
  for(j in (d+2):(2*d)){
    # extend diagonally
    filled_mat[(j-d):(j-1),j] = filled_mat[1:d,(d+1)]

    # extend the restricted upper triangular part
    filled_mat[1:(j-d-1),j] = corr_index(diag_ele = diag(mat),
                                         index = filled_mat[1:(j-d-1),(j-d)],
                                         is.next = TRUE)
  }




  return(filled_mat)
}


#' Helper function 3
#'
#' Helper function, extend column d+2:2d and corresponding rows
#' and then check whether the filled array is a valid svine array.
#' @param mat   A 2dx2d vine matrix where the upper 1:d+1 columns have been filled.
#' @returns     A filled 2dx2d vine matrix if it is a valid svine array. NULL
#'              otherwise.
fill_check <- function(mat){
  filled_mat = extend_to_2d(mat)
  if(SVineMatrix_mk1_Check(filled_mat)==1){
    return(filled_mat)
  }else{
    return(NULL)
  }
}


#' Generate all d-dimensional Markov order 1 S-vine array
#'
#' @description
#' Generate all possible d-dimensional Markov order 1 S-vine array
#' with 1=X_1, 2=X_2,...,d=X_d, d+1=X_1',...,2d=X_d' on diagonal
#' @param d dimension k >= 2, default is k = 2, must be an integer
#' @returns List of 2dx2d S-vine arrays
#' @examples
#' svine_d2_list = svine_mk1(2)
#' svine_d3_list = svine_mk1(3)
#' svine_d4_list = svine_mk1(4)
#'
#' @export
svine_mk1 <- function(d=2){

  # fill in fixed cells
  svine_fixed = matrix(nrow = 2*d, ncol = 2*d)
  diag(svine_fixed) = 1:(2*d)
  svine_fixed[1,c(2,d+1)] = c(1,1)

  # initialize list of possible s-vine arrays
  svine_list_prev = list(svine_fixed)

  for(i in 3:(d+1)){
    svine_list_new = list()
    for(mat in svine_list_prev){
      poss_i = generate_specific_upper_tri_matrices(mat, k=i)
      svine_list_new = append(svine_list_new, poss_i)
    }
    svine_list_prev = svine_list_new
  }

  # Filter the list based on SVineMatrix_mk1_Check
  filtered_list <- lapply(svine_list_new, fill_check)
  filtered_list <- filtered_list[!sapply(filtered_list, is.null)]

  return(filtered_list)
}


#' Generate all d-dimensional Markov order p S-vine array
#'
#' @description Generate all possible d-dimensional Markov order p S-vine array
#'              with 1=X_1, 2=X_2,...,d=X_d, d+1=X_1',...,2d=X_d' on diagonal
#' @param d dimension k >= 2, default is k = 2, must be an integer
#' @return List of d(p+1)xd(p+1) S-vine arrays
#' @examples
#' svine_d3_p2_list = svine_mkp(d=3,p=2)
#'
#'@export
svine_mkp <- function(d=2, p=1){

  svine_mk1_list = svine_mk1(d)
  mkp_list = lapply(svine_mk1_list, mk1_to_mkp, p=p)

  return(mkp_list)
}



