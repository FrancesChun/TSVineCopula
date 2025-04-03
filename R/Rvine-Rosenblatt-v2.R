# Rosenblatt transform (forward and inverse) for R-vine based
# on diagonal order in vine array
# See Section 6.9.1 of Joe (2014).

# Scalar versions and then vectorized versions

library(fastmatrix)

#' @description
#' Function to convert a vine array A to a maximum matrix. Algorithm 5 in Joe (2014).
#' 
#' @param A a dxd vine array of R-vine in standard order with 1:d on diagonal, d>=2
#' @param iprint print flag for intermediate calculations
#' @param str string to describe vine for printing if iprint=T
#'
#' @return Returns a list with two components 
#' \enumerate{
#' \item \code{M} - dxd array with m_{kj}= max a_{k1},..,a_{kj}
#' \item \code{icomp} - dxd indicator array on whether back step [k,j] is needed ; icomp[k-1,m_{kj}=1 if a_{kj}<m_{kj} for k>=2
#' }
#'
varray2M = function(A, iprint=FALSE, str="")
{ d = ncol(A)
  d1 = d-1
  M = A
  icomp = matrix(0,d,d)
  
  # Note: if d == 2, only 1 possible A: [1 1; 0 2]
  # in this case, M = A, and icomp = matrix(0,d,d)
  
  if(d >= 3){
    
    # modify M
    for(k in 2:d1){ 
      for(j in (k+1):d) M[k,j] = max(M[k-1,j],A[k,j])
    }
    
    if(iprint) { cat("\n",str,"\n"); print(A); print(M) }
    
    # modify icomp
    for(k in 2:d1){ 
      for(j in (k+1):d){ 
        if(A[k,j]<M[k,j]) icomp[k-1,M[k,j]] = 1
      }
    }
  }
 
  if(iprint) print(icomp)
  list(mxarray=M, icomp=icomp)
}  


#' @description
#' Function to rotate a matrix by 180 degree
#' 
#' @param mat a matrix to be rotated
#'
#' @return Returns the rotated matrix
#'
rotate_180 <- function(mat){
  mat[nrow(mat):1, ncol(mat):1]
}




#======================================================================

#' @description 
#' Rosenblatt transform (forward) for R-vine based on diagonal order in vine array
#'  This is R-vine simulation if p is a vector of independent U(0,1) 
#'  Algorithms 17 and 18 in Joe (2014).
#'  This function could be vectorized for simulations (random sample size n).
#' @param p vector of length d with values in interval (0,1)
#' @param A dxd vine array with 1:d on diagonal, d>=2
#'     if truncated vine, only rows 1 to ntrunc are used.
#' @param fam dxd code matrix (VineCopula copula family index) .
#' @param param1 dxd parameter matrix (VineCopula first copula parameter) 
#'     param1[l,j] for tree l, variable j.
#' @param param2 dxd parameter matrix (VineCopula second copula parameter) 
#'     param2[l,j] for tree l, variable j; 0 if not needed.
#' @param ntrunc truncation level, integer between 1 and d-1 (d-1 for no truncation)
#' @param varname character vector of variable names, optional
#' @param iprint print flag for intermediate calculations
#' @return d-vector of values in (0,1)
#'   If p is vector of independent U(0,1), then output is a random vector from the vine copula
#'     based on A, fam, param1, param2, ntrunc
#' that is, u1=p1, C_{2|1}^{-1}(p2|u1), C_{3|12}^{-1}(p3|u1,u2), ...  C_{d|1..d-1}^{-1}(p_d|u[1:(d-1)])
#' @example
#' to add
#' @details
#' BiCopHinv1 is from VineCopula for inverse of conditional cdf C_{2|1}^{-1}(v|u)
#' BiCopHfunc2 is from VineCopula for conditional cdf C_{1|2}(u|v)
#'
rvineqcond = function(p, A,fam,param1,param2,ntrunc=0,
   varname=numeric(0),iprint=FALSE)
{ 
  # if lower triangular matrix, convert to upper tri matrix
  if(is.lower.tri(A)){
    A = rotate_180(A)
  }
  if(is.lower.tri(fam)){
    fam = rotate_180(fam)
  }
  if(is.lower.tri(param1)){
    param1 = rotate_180(param1)
  }
  if(is.lower.tri(param2)){
    param2 = rotate_180(param2)
  }
  
  d = ncol(A)
  ntrunc = floor(ntrunc)
  if(ntrunc<1 | ntrunc>=d) ntrunc = d-1
  if(length(varname)!=d) varname = paste0("u",1:d)
  # convert diagonal of A to 1:d to use Algorithm 17 of Joe (2014)
  diagA = diag(A)
  dict = data.frame(Col1=c(0, diagA), Col2=0:d)
  A = matrix(dict$Col2[match(A, dict$Col1)], nrow=d, byrow=FALSE)
  # A now has 1:d on diagonal, dict is a dictionary of the mapping
  out = varray2M(A)
  M = out$mxarray
  icomp = out$icomp
  qq = matrix(0,d,d); v = matrix(0,d,d)
  u = rep(0,d)
  u[1] = p[1]; qq[1,1] = p[1]; qq[2,2] = p[2];
  u[2] = BiCopHinv1(p[1],p[2],fam[1,2],param1[1,2],param2[1,2]); 
  qq[1,2] = u[2]
  if(icomp[1,2]==1) 
  { v[1,2] = BiCopHfunc2(u[1],u[2],fam[1,2],param1[1,2],param2[1,2]) }
  
  # the main loop
  if(d >= 3){
    for(j in 3:d){  # variable index
      tt = min(ntrunc,j-1)
      qq[tt+1,j] = p[j]
      if(tt>1){
        for(ell in seq(tt,2)){ # tree index
          s = ifelse(A[ell,j]==M[ell,j], qq[ell,A[ell,j]], v[ell-1,M[ell,j]] )
          qq[ell,j] = BiCopHinv1(s,qq[ell+1,j], fam[ell,j], param1[ell,j], param2[ell,j]); 
        }
      }
      qq[1,j] = BiCopHinv1(u[A[1,j]],qq[2,j], fam[1,j], param1[1,j], param2[1,j])
      u[j] = qq[1,j] 
    
      # set up for next iteration (not needed for last j=d)
      v[1,j] = BiCopHfunc2(u[A[1,j]],u[j], fam[1,j],param1[1,j], param2[1,j])
      if(tt>1){ 
        for(ell in 2:tt){ 
          s = ifelse(A[ell,j]==M[ell,j], qq[ell,A[ell,j]], v[ell-1,M[ell,j]] )
          if(icomp[ell,j]==1){ 
            v[ell,j] = BiCopHfunc2(s, qq[ell,j], fam[ell,j], param1[ell,j], param2[ell,j]); 
          }
        }
      }
    }
  }
  
  if(iprint) { print(qq); print(v) }
  uvec = u[order(dict$Col1[2:(d+1)])]
  names(uvec) = varname
  uvec
}

#------------------------------------------------------------

# Rosenblatt inverse

#' @description 
#' Rosenblatt transform (inverse) for R-vine based on diagonal order in vine array
#'  If uvec is generated from the vine copula based on A,fam,param1,param1,ntrunc
#'    then this function recovers the vector of U(0,1) used for generation.
#'  The steps are embeded in Algorithm 4 in Joe (2014).
#' @param uvec vector of length d with values in interval (0,1)
#' @param A dxd vine array with 1:d on diagonal, d>=2
#'     if truncated vine, only rows 1 to ntrunc are used.
#' @param fam dxd code matrix (VineCopula copula family index) .
#' @param param1 dxd parameter matrix (VineCopula first copula parameter) 
#'     param1[l,j] for tree l, variable j.
#' @param param2 dxd parameter matrix (VineCopula second copula parameter) 
#'     param2[l,j] for tree l, variable j; 0 if not needed.
#' @param ntrunc truncation level, integer between 1 and d-1 (d-1 for no truncation)
#' @param varname character vector of variable names, optional
#' @param iprint print flag for intermediate calculations
#' @return d-vector of values in (0,1)
#' p1=u1, p2=C_{2|1}(u2|u1), p3=C_{3|12}(u3|u1,u2), ...  p_d=C_{d|1..d-1}(u_d|u[1:(d-1)])
#' @example
#' to add
#' @details
#' BiCopHfunc1 is from VineCopula for conditional cdf C_{2|1}(v|u)
#' BiCopHfunc2 is from VineCopula for conditional cdf C_{1|2}(u|v)
#' 
rvinepcond = function(uvec, A,fam,param1,param2,ntrunc=0,
   varname=numeric(0),iprint=FALSE)
{
  # if lower triangular matrix, convert to upper tri matrix
  if(is.lower.tri(A)){
    A = rotate_180(A)
  }
  if(is.lower.tri(fam)){
    fam = rotate_180(fam)
  }
  if(is.lower.tri(param1)){
    param1 = rotate_180(param1)
  }
  if(is.lower.tri(param2)){
    param2 = rotate_180(param2)
  }
  
  d = ncol(A)
  ntrunc = floor(ntrunc)
  if(ntrunc<1 | ntrunc>=d) ntrunc = d-1
  if(length(varname)!=d) varname = paste0("p",1:d)
  # convert diagonal of A to 1:d to use Algorithm 17 of Joe (2014)
  diagA = diag(A)
  dict = data.frame(Col1=c(0, diagA), Col2=0:d)
  A = matrix(dict$Col2[match(A, dict$Col1)], nrow=d, byrow=FALSE)
  # A now has 1:d on diagonal, dict is a dictionary of the mapping
  out = varray2M(A)
  M = out$mxarray
  icomp = out$icomp
  vinepcond = rep(0,d); vinepcond[1] = uvec[1]
  v = rep(0,d); vp = rep(0,d); s = rep(0,d);
  # tree 2
  if(ntrunc>=2)
  { for(j in 2:d) 
    { 
      if(icomp[1,j]==1) vp[j] = BiCopHfunc2(uvec[A[1,j]],uvec[j],fam[1,j],param1[1,j],param2[1,j])
      v[j] = BiCopHfunc1(uvec[A[1,j]],uvec[j],fam[1,j],param1[1,j],param2[1,j])
    }
    vinepcond[2] = v[2]
    for(j in 3:d) 
    { if(A[2,j]<M[2,j]) s[j] = vp[M[2,j]] else s[j] = v[A[2,j]] } 
    w = v; wp = vp
  }
  # remaining trees
  if(ntrunc>=3)
  { for(ell in 3:ntrunc)
    { for(j in ell:d) 
      { 
        if(icomp[ell-1,j]==1) vp[j] = BiCopHfunc2(s[j],w[j],fam[ell-1,j],param1[ell-1,j],param2[ell-1,j])
        v[j] = BiCopHfunc1( s[j],w[j], fam[ell-1,j],param1[ell-1,j],param2[ell-1,j])
      }
      vinepcond[ell] = v[ell]
      for(j in (ell+1):d) 
      { if(A[ell,j]<M[ell,j]) s[j] = vp[M[ell,j]] else s[j] = v[A[ell,j]] } 
      w = v; 
    }
  }
  if(ntrunc>1 & ntrunc<d)
  { ell = ntrunc+1
    for(j in ell:d) 
    { 
      if(icomp[ell-1,j]==1) vp[j] = BiCopHfunc2(s[j],w[j],fam[ell-1,j],param1[ell-1,j],param2[ell-1,j])
      v[j] = BiCopHfunc1(s[j],w[j],fam[ell-1,j],param1[ell-1,j],param2[ell-1,j])
      vinepcond[j] = v[j] 
    }
  }
  if(ntrunc==1)
  { ell = 2
    for(j in ell:d) 
    { v[j] = BiCopHfunc1(uvec[A[1,j]],uvec[j],fam[1,j],param1[ell-1,j],param2[ell-1,j])
      vinepcond[j] = v[j] 
    }
  }
  vinepcond = vinepcond[order(dict$Col1[2:(d+1)])]
  names(vinepcond) = varname
  vinepcond
}

#======================================================================

# Vectorized versions

#' @description 
#' Rosenblatt transform (forward) for R-vine based on diagonal order in vine array
#'  This is R-vine simulation if p is a vector of independent U(0,1) 
#'  Algorithms 17 and 18 in Joe (2014).
#'  This function could be vectorized for simulations (random sample size n).
#' @param p nxd matrix with values in interval (0,1)
#' @param A dxd vine array with 1:d on diagonal, d>=3 
#'     if truncated vine, only rows 1 to ntrunc are used.
#' @param fam dxd code matrix (VineCopula copula family index) .
#' @param param1 dxd parameter matrix (VineCopula first copula parameter) 
#'     param1[l,j] for tree l, variable j.
#' @param param2 dxd parameter matrix (VineCopula second copula parameter) 
#'     param2[l,j] for tree l, variable j; 0 if not needed.
#' @param ntrunc truncation level, integer between 1 and d-1 (d-1 for no truncation)
#' @param varname character vector of variable names, optional
#' @param iprint print flag for intermediate calculations
#' @return nxd matrix of values in (0,1)
#'   If p is matrix of independent U(0,1), then output is a random sample from the vine copula
#'     based on A, fam, param1, param2, ntrunc
#' that is, by row:
#'  u1=p1, C_{2|1}^{-1}(p2|u1), C_{3|12}^{-1}(p3|u1,u2), ...  C_{d|1..d-1}^{-1}(p_d|u[1:(d-1)])
#' @example
#' to add
#' @details
#' BiCopHinv1 is from VineCopula for inverse of conditional cdf C_{2|1}^{-1}(v|u)
#' BiCopHfunc2 is from VineCopula for conditional cdf C_{1|2}(u|v)
#'
rvineqcond_vec = function(p, A,fam,param1,param2,ntrunc=0,
   varname=numeric(0),iprint=FALSE)
{ d = ncol(A)
  ntrunc = floor(ntrunc)
  if(ntrunc<1 | ntrunc>=d) ntrunc = d-1
  if(length(varname)!=d) varname = paste0("u",1:d)
  # convert diagonal of A to 1:d to use Algorithm 17 of Joe (2014)
  diagA = diag(A)
  dict = data.frame(Col1=c(0, diagA), Col2=0:d)
  A = matrix(dict$Col2[match(A, dict$Col1)], nrow=d, byrow=FALSE)
  # A now has 1:d on diagonal, dict is a dictionary of the mapping
  out = varray2M(A)
  M = out$mxarray
  icomp = out$icomp
  n = nrow(p)
  qq = array(0,c(n,d,d)); v = array(0,c(n,d,d))
  u = matrix(0,n,d)
  u[,1] = p[,1]; qq[,1,1] = p[,1]; qq[,2,2] = p[,2];
  u[,2] = BiCopHinv1(p[,1],p[,2],fam[1,2],param1[1,2],param2[1,2]); 
  qq[,1,2] = u[,2]
  if(icomp[1,2]==1) 
  { v[,1,2] = BiCopHfunc2(u[,1],u[,2],fam[1,2],param1[1,2],param2[1,2]) }
  # the main loop 
  for(j in 3:d)  # variable index
  { tt = min(ntrunc,j-1)
    qq[,tt+1,j] = p[,j]
    if(tt>1)
    { for(ell in seq(tt,2)) # tree index
      { if(A[ell,j]==M[ell,j]) { s = qq[,ell,A[ell,j]] }
        else { s = v[,ell-1,M[ell,j]] }
        qq[,ell,j] = BiCopHinv1(s,qq[,ell+1,j], fam[ell,j], param1[ell,j], param2[ell,j]); 
      }
    } 
    qq[,1,j] = BiCopHinv1(u[,A[1,j]],qq[,2,j], fam[1,j], param1[1,j], param2[1,j])
    u[,j] = qq[,1,j] 
    # set up for next iteration (not needed for last j=d)
    v[,1,j] = BiCopHfunc2(u[,A[1,j]],u[,j], fam[1,j],param1[1,j], param2[1,j])
    if(tt>1)
    { for(ell in 2:tt)
      { if(A[ell,j]==M[ell,j]) { s = qq[,ell,A[ell,j]] }
        else { s = v[,ell-1,M[ell,j]] }
        if(icomp[ell,j]==1) 
        { v[,ell,j] = BiCopHfunc2(s, qq[,ell,j], fam[ell,j], param1[ell,j], param2[ell,j]) 
        }
      }
    }
  }
  if(iprint) { print(qq[1:2,]); print(v[1:2,]) }
  u = u[,order(dict$Col1[2:(d+1)])]
  colnames(u) = varname
  u
}

#------------------------------------------------------------

# vectorized version Rosenblatt inverse

#' @description 
#' Rosenblatt transform (inverse) for R-vine based on diagonal order in vine array
#'  If umat is generated from the vine copula based on A,fam,param1,param1,ntrunc
#'    then this function recovers the vector of U(0,1) used for generation.
#'  The steps are embeded in Algorithm 4 in Joe (2014).
#' @param umat nxd matrix with values in interval (0,1)
#' @param A dxd vine array with 1:d on diagonal, d>=3 
#'     if truncated vine, only rows 1 to ntrunc are used.
#' @param fam dxd code matrix (VineCopula copula family index) .
#' @param param1 dxd parameter matrix (VineCopula first copula parameter) 
#'     param1[l,j] for tree l, variable j.
#' @param param2 dxd parameter matrix (VineCopula second copula parameter) 
#'     param2[l,j] for tree l, variable j; 0 if not needed.
#' @param ntrunc truncation level, integer between 1 and d-1 (d-1 for no truncation)
#' @param varname character vector of variable names, optional
#' @param iprint print flag for intermediate calculations
#' @return nxd matrix of values in (0,1), by row:
#' p1=u1, p2=C_{2|1}(u2|u1), p3=C_{3|12}(u3|u1,u2), ...  p_d=C_{d|1..d-1}(u_d|u[1:(d-1)])
#' @example
#' to add
#' @details
#' BiCopHfunc1 is from VineCopula for conditional cdf C_{2|1}(v|u)
#' BiCopHfunc2 is from VineCopula for conditional cdf C_{1|2}(u|v)
#' 
rvinepcond_vec = function(umat, A,fam,param1,param2,ntrunc=0,
   varname=numeric(0),iprint=FALSE)
{ d = ncol(A)
  ntrunc = floor(ntrunc)
  if(ntrunc<1 | ntrunc>=d) ntrunc = d-1
  if(length(varname)!=d) varname = paste0("p",1:d)
  # convert diagonal of A to 1:d to use Algorithm 17 of Joe (2014)
  diagA = diag(A)
  dict = data.frame(Col1=c(0, diagA), Col2=0:d)
  A = matrix(dict$Col2[match(A, dict$Col1)], nrow=d, byrow=FALSE)
  # A now has 1:d on diagonal, dict is a dictionary of the mapping
  out = varray2M(A)
  M = out$mxarray
  icomp = out$icomp
  n = nrow(umat)
  vinepcond = matrix(0,n,d); vinepcond[,1] = umat[,1]
  v = matrix(0,n,d); vp = matrix(0,n,d); s = matrix(0,n,d);
  # tree 2
  if(ntrunc>=2)
  { for(j in 2:d) 
    { 
      if(icomp[1,j]==1) vp[,j] = BiCopHfunc2(umat[,A[1,j]],umat[,j],fam[1,j],param1[1,j],param2[1,j])
      v[,j] = BiCopHfunc1(umat[,A[1,j]],umat[,j],fam[1,j],param1[1,j],param2[1,j])
    }
    vinepcond[,2] = v[,2]
    for(j in 3:d) 
    { if(A[2,j]<M[2,j]) s[,j] = vp[,M[2,j]] else s[,j] = v[,A[2,j]] } 
    w = v; wp = vp
  }
  # remaining trees
  if(ntrunc>=3)
  { for(ell in 3:ntrunc)
    { for(j in ell:d) 
      { 
        if(icomp[ell-1,j]==1) vp[,j] = BiCopHfunc2(s[,j],w[,j],fam[ell-1,j],param1[ell-1,j],param2[ell-1,j])
        v[,j] = BiCopHfunc1( s[,j],w[,j], fam[ell-1,j],param1[ell-1,j],param2[ell-1,j])
      }
      vinepcond[,ell] = v[,ell]
      for(j in (ell+1):d) 
      { if(A[ell,j]<M[ell,j]) s[,j] = vp[,M[ell,j]] else s[,j] = v[,A[ell,j]] } 
      w = v; 
    }
  }
  if(ntrunc>1 & ntrunc<d)
  { ell = ntrunc+1
    for(j in ell:d) 
    { 
      if(icomp[ell-1,j]==1) vp[,j] = BiCopHfunc2(s[,j],w[,j],fam[ell-1,j],param1[ell-1,j],param2[ell-1,j])
      v[,j] = BiCopHfunc1(s[,j],w[,j],fam[ell-1,j],param1[ell-1,j],param2[ell-1,j])
      vinepcond[,j] = v[,j] 
    }
  }
  if(ntrunc==1)
  { ell = 2
    for(j in ell:d) 
    { v[,j] = BiCopHfunc1(umat[,A[1,j]],umat[,j],fam[1,j],param1[ell-1,j],param2[ell-1,j])
      vinepcond[,j] = v[,j] 
    }
  }
  vinepcond = vinepcond[,order(dict$Col1[2:(d+1)])]
  colnames(vinepcond) = varname
  vinepcond
}


#======================================================================

# testing here

wrap=function() {
library(VineCopula)

## d=2 examples
p2 = c(0.3, 0.2)
A2 = matrix(c(1,1, 0,2), nrow = 2, byrow = T)
apar1_d2 = matrix(c(0,2, 0,0), nrow=2, ncol=2, byrow=T)
apar2_d2 = matrix(0, nrow=2, ncol=2)
afam_d2 = matrix(4, nrow=2, ncol=2) # all Gumbel

u2 = rvineqcond(p2,A2,afam_d2,apar1_d2,apar2_d2,iprint=FALSE)
print(u2)
rvinepcond(u2,A2,afam_d2,apar1_d2,apar2_d2,iprint=FALSE)
  
## d=4 examples
umat4 = matrix(c(.1,.4,.6,.7, .2,.9,.4,.3), 2,4, byrow=T)

## first example for d=4, all Gumbel
apar1 = matrix(c(0, 2,2.1,2.2, 0,0, 1.5,1.4, 0,0,0,1.2, 0,0,0,0),4,4, byrow=T)
apar2 = matrix(0,4,4)
afam4 = matrix(4,4,4) # all Gumbel

A4 = matrix(c(1,1,1,1, 0,2,2,2, 0,0,3,3, 0,0,0,4), 4,4,byrow=T)
A4p = matrix(c(4,4,4,4, 0,3,3,3, 0,0,2,2, 0,0,0,1), 4,4,byrow=T)


# check inverse operation
for(i in 1:nrow(umat4))
{ p = rvinepcond(umat4[i,],A4,afam4,apar1,apar2,ntrunc=3,iprint=T)
  print(p)
  u = rvineqcond(p,A4,afam4,apar1,apar2,ntrunc=3,iprint=FALSE)
  print(u)
}
#       p1        p2        p3        p4 
#0.1000000 0.7794878 0.9351206 0.9545815 
# u1  u2  u3  u4 
#0.1 0.4 0.6 0.7 
#        p1         p2         p3         p4 
#0.20000000 0.99443237 0.08760435 0.13443775 
# u1  u2  u3  u4 
#0.2 0.9 0.4 0.3 


for(itrunc in 3:1)
{ p = rvinepcond(umat4[1,],A4,afam4,apar1,apar2,ntrunc=itrunc,iprint=T)
  print(p)
  u = rvineqcond(p,A4,afam4,apar1,apar2,ntrunc=itrunc,iprint=FALSE)
  print(u)
}

#       p1        p2        p3        p4 
#0.1000000 0.7794878 0.9351206 0.9545815  # ntrunc=3
# u1  u2  u3  u4 
#0.1 0.4 0.6 0.7 
#       p1        p2        p3        p4 
#0.1000000 0.7794878 0.9351206 0.9807228   # ntrunc=2
# u1  u2  u3  u4 
#0.1 0.4 0.6 0.7 
#       p1        p2        p3        p4 
#0.1000000 0.7794878 0.9345996 0.9741819   # ntrunc=1
# u1  u2  u3  u4 
#0.1 0.4 0.6 0.7 

pmat = rvinepcond_vec(umat4,A4,afam4,apar1,apar2,ntrunc=3,iprint=T)
print(pmat)
#      p1        p2         p3        p4
#[1,] 0.1 0.7794878 0.93512060 0.9545815
#[2,] 0.2 0.9944324 0.08760435 0.1344377
u = rvineqcond_vec(pmat,A4,afam4,apar1,apar2,ntrunc=3,iprint=FALSE)
print(u)
#      u1  u2  u3  u4
#[1,] 0.1 0.4 0.6 0.7
#[2,] 0.2 0.9 0.4 0.3

pmat = rvinepcond_vec(umat4,A4p,afam4,apar1,apar2,ntrunc=3,iprint=T)
print(pmat)
#            p1         p2        p3  p4
#[1,] 0.9545815 0.93512060 0.7794878 0.1
#[2,] 0.1344377 0.08760435 0.9944324 0.2
# Be careful in inverse with diag of A matrix is not 1:d
u = rvineqcond_vec(pmat[,4:1],A4p,afam4,apar1,apar2,ntrunc=3,iprint=FALSE)
print(u)
#      u1  u2  u3  u4
#[1,] 0.7 0.6 0.4 0.1
#[2,] 0.3 0.4 0.9 0.2

ii = 1:3 # subvine
p3 = rvinepcond(umat4[1,ii],A4[ii,ii],afam4[ii,ii],apar1[ii,ii],apar2[ii,ii],ntrunc=2,iprint=F)
print(p3)
#  0.1000000 0.7794878 0.9351206
u3 = rvineqcond(p3,A4[ii,ii],afam4[ii,ii],apar1[ii,ii],apar2[ii,ii],ntrunc=2,iprint=FALSE)
# [1] 0.1 0.4 0.6

## second vine example with d=4

bpar1= matrix(c(0,0.73,0.66,0.58, 0,0,4.5,3.5, 0,0,0,1.2, 0,0,0,0),4,4,byrow=T)
bpar2= matrix(c(0,2.64,2.06,1.71, 0,0,0.8,0.5, 0,0,0,0, 0,0,0,0),4,4,byrow=T)
bfam4 = matrix(c(0,7,7,7, 0,0,10,20, 0,0,0,5, 0,0,0,0), 4,4, byrow=T)
B4 = matrix(c(1,1,2,3, 0,2,1,2, 0,0,3,1, 0,0,0,4), 4,4,byrow=T)
B4p = matrix(c(4,4,3,2, 0,3,4,3, 0,0,2,4, 0,0,0,1), 4,4,byrow=T)

for(itrunc in 3:1)
{ p = rvinepcond(umat4[1,],B4,bfam4,bpar1,bpar2,ntrunc=itrunc,iprint=T)
  print(p)
  u = rvineqcond(p,B4,bfam4,bpar1,bpar2,ntrunc=itrunc,iprint=FALSE)
  print(u)
}
#       p1        p2        p3        p4 
#0.1000000 0.9763016 0.9908524 0.8942927  ntrunc=3
# u1  u2  u3  u4 
#0.1 0.4 0.6 0.7 
#       p1        p2        p3        p4 
#0.1000000 0.9763016 0.9908524 0.8173441  ntrunc=2
# u1  u2  u3  u4 
#0.1 0.4 0.6 0.7 
#       p1        p2        p3        p4 
#0.1000000 0.9763016 0.8040850 0.7054554  ntrunc=1
# u1  u2  u3  u4 
#0.1 0.4 0.6 0.7 


pmat = rvinepcond_vec(umat4,B4,bfam4,bpar1,bpar2,ntrunc=3,iprint=T)
print(pmat)
#      p1        p2        p3        p4
#[1,] 0.1 0.9763016 0.9908524 0.8942927
#[2,] 0.2 0.9998121 0.0677722 0.2167492
u = rvineqcond_vec(pmat,B4,bfam4,bpar1,bpar2,ntrunc=3,iprint=FALSE)
print(u)
#      u1  u2  u3  u4
#[1,] 0.1 0.4 0.6 0.7
#[2,] 0.2 0.9 0.4 0.3

pmat = rvinepcond_vec(umat4,B4p,bfam4,bpar1,bpar2,ntrunc=3,iprint=T)
print(pmat)
#            p1        p2        p3  p4
#[1,] 0.8942927 0.9908524 0.9763016 0.1
#[2,] 0.2167492 0.0677722 0.9998121 0.2
# permute pmat in same way as B4 -> B4p
u = rvineqcond_vec(pmat[,4:1],B4p,bfam4,bpar1,bpar2,ntrunc=3,iprint=FALSE)
print(u)
#     u1  u2  u3  u4
#[1,] 0.7 0.6 0.4 0.1
#[2,] 0.3 0.4 0.9 0.2

ii = 1:3 # subvine
p3 = rvinepcond(umat4[1,ii],B4[ii,ii],bfam4[ii,ii],bpar1[ii,ii],bpar2[ii,ii],ntrunc=2,iprint=F)
print(p3)
#       p1        p2        p3 
#0.1000000 0.9763016 0.9908524
u3 = rvineqcond(p3,B4[ii,ii],bfam4[ii,ii],bpar1[ii,ii],bpar2[ii,ii],ntrunc=2,iprint=FALSE)
print(u3)
# u1  u2  u3 
#0.1 0.4 0.6


#----------------------------------------------------------------------

## d=5 examples with BB1 family 7, BB10 family 10, and Frank family 5

umat5 = matrix(c(.1,.4,.6,.7,.4, .2,.9,.4,.3,.1), 2,5, byrow=T)

bpar1= matrix(c(0,0.73,0.66,0.58,0.51, 0,0,4.5,3.5,0.4, 0,0,0,1.2,1.1, 0,0,0,0,0.6, 0,0,0,0,0),5,5,byrow=T)
bpar2= matrix(c(0,2.64,2.06,1.71,1.47, 0,0,0.8,0.5,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0),5,5,byrow=T)
bfam5 = matrix(c(0,7,7,7,7, 0,0,10,20,1, 0,0,0,5,5, 0,0,0,0,5, 0,0,0,0,0),
    5,5, byrow=T)
B5 = matrix(c(1,1,1,2,4, 0,2,2,1,2, 0,0,3,3,1, 0,0,0,4,3, 0,0,0,0,5), 5,5,byrow=T)
# interchange 5,4,2,1,3 diag
B5p = matrix(c(5,5,4,2,1, 0,4,5,4,2, 0,0,2,5,4, 0,0,0,1,5, 0,0,0,0,3), 5,5,byrow=T)

for(itrunc in 4:1)
{ p = rvinepcond(umat5[1,],B5,bfam5,bpar1,bpar2,ntrunc=itrunc,iprint=T)
  print(p)
  u = rvineqcond(p,B5,bfam5,bpar1,bpar2,ntrunc=itrunc,iprint=FALSE)
  print(u)
}
#       p1        p2        p3        p4        p5 
#0.1000000 0.9763016 0.8595233 0.9245306 0.4096694  ntrunc=4
# u1  u2  u3  u4  u5 
#0.1 0.4 0.6 0.7 0.4 
#       p1        p2        p3        p4        p5 
#0.1000000 0.9763016 0.8595233 0.9245306 0.4524849  ntrunc=3
# u1  u2  u3  u4  u5 
#0.1 0.4 0.6 0.7 0.4 
#       p1        p2        p3        p4        p5 
#0.1000000 0.9763016 0.8595233 0.9470452 0.3273888  ntrunc=2
# u1  u2  u3  u4  u5 
#0.1 0.4 0.6 0.7 0.4 
#       p1        p2        p3        p4        p5 
#0.1000000 0.9763016 0.9839769 0.8602944 0.2100976  ntrunc=1
# u1  u2  u3  u4  u5 
#0.1 0.4 0.6 0.7 0.4 

p = rvinepcond_vec(umat5,B5p,bfam5,bpar1,bpar2,ntrunc=4,iprint=T)
print(p)
#            p1        p2         p3        p4  p5
#[1,] 0.8942927 0.9908524 0.35344018 0.9763016 0.1
#[2,] 0.2167492 0.0677722 0.04089308 0.9998121 0.2
# no permutation of above?

u = rvineqcond_vec(p[,c(5,4,2,1,3)],B5p,bfam5,bpar1,bpar2,ntrunc=4,iprint=FALSE)
print(u)
#      u1  u2  u3  u4  u5
#[1,] 0.7 0.6 0.4 0.4 0.1
#[2,] 0.3 0.4 0.1 0.9 0.2

print(umat5)
#     [,1] [,2] [,3] [,4] [,5]
#[1,]  0.1  0.4  0.6  0.7  0.4
#[2,]  0.2  0.9  0.4  0.3  0.1


ii = 1:3 # subvine
p3 = rvinepcond(umat5[1,ii],B5[ii,ii],bfam5[ii,ii],bpar1[ii,ii],bpar2[ii,ii],ntrunc=2,iprint=F)
print(p3)
#       p1        p2        p3 
#0.1000000 0.9763016 0.8595233 
u3 = rvineqcond(p3,B5[ii,ii],bfam5[ii,ii],bpar1[ii,ii],bpar2[ii,ii],ntrunc=2,iprint=FALSE)
print(u3)
# u1  u2  u3 
#0.1 0.4 0.6 

ii = 1:2 # subvine (doesn't work)
ii = 1:4 # subvine
p4 = rvinepcond_vec(umat5[,ii],B5[ii,ii],bfam4[ii,ii],bpar1[ii,ii],bpar2[ii,ii],ntrunc=3,iprint=F)
print(p4)
#      p1        p2        p3         p4
#[1,] 0.1 0.9763016 0.8595233 0.92453058
#[2,] 0.2 0.9998121 0.1270160 0.09538316
u4 = rvineqcond_vec(p4,B5[ii,ii],bfam4[ii,ii],bpar1[ii,ii],bpar2[ii,ii],ntrunc=3,iprint=FALSE)
print(u4)
#      u1  u2  u3  u4
#[1,] 0.1 0.4 0.6 0.7
#[2,] 0.2 0.9 0.4 0.3

ii = 1:3 # subvine
p3 = rvinepcond_vec(umat5[,ii],B5[ii,ii],bfam5[ii,ii],bpar1[ii,ii],bpar2[ii,ii],ntrunc=3,iprint=F)
print(p3)
#      p1        p2        p3
#[1,] 0.1 0.9763016 0.8595233
#[2,] 0.2 0.9998121 0.1270160

u3 = rvineqcond_vec(p3,B5[ii,ii],bfam5[ii,ii],bpar1[ii,ii],bpar2[ii,ii],ntrunc=3,iprint=FALSE)
print(u3)
#      u1  u2  u3
#[1,] 0.1 0.4 0.6
#[2,] 0.2 0.9 0.4
} # end of wrap function
