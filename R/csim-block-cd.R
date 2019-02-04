#' A block coordinate descent function
#'
#' \code{block.cd} implements the block-coordinate descent algorithm of Radchenko (2015) to optimize the (active) single index coefficients.
#'
#' @param w.a  a warm start (initial) value of the current nonzero single index coefficients.
#' @param ind.m   a value of the index used in defining each coordinate-wise block.
#' @param ind.old   a vector of the indices corresponding to the current nonzero coefficients.
#' @param u  a vector of the current singe index variable, u = alpha'X.
#' @param y   treatment outcomes, a n-by-1 vector.
#' @param Tr  treatment indicators, a n-by-1 vector; each element represents one of the K available treatment options.
#' @param X.a   a matrix of the pre-treatment covarates, associated with the current nonzero single index coefficients.
#' @param d.link.fn.obj  an object of the class \code{d.link.fn.obj} passed from the previous iteration.
#' @param it.max  an integer value specifying the maximum number of iterations for each coordinate.
#' @param eps  a value specifying the converge criterion of algorithm.
#' @param trace  if \code{TRUE}, show the trace of the fitting procedure; the default is \code{FALSE}.
#' @param nbasis.t  a length K+1 vector; each element specifies the number of B-spline basis funtions for approximating the treatment-specific link function; the last element is for the "main effect" link function; the default is \code{nbasis.t=NULL}, and will be determined depending on the sample size.
#' @param rho.grid  a grid vector of (ridge-type) smoothing parameters for approximating the link functions.
#' @param linear.link  if \code{TRUE}, restrict the link functions to be linear functions; the default is \code{FALSE}.
#' @param ortho.constr  the constraint that separates the interaction effects from the main effect (without this, the interaction effect can be confounded by the main effect); the default is \code{TRUE}.
#'
#' @return a list of information of the optimized single index coefficients including
#' \item{w.a}{a vector of the optimized single index coefficients associated with the active index set.} \item{u}{a vector of the optimized single index variable.} \item{it}{the number of iterations.}
#'
#' @seealso \code{pred.csim},  \code{fit.csim},  \code{fit.csim.cv}
#'
######################################################################
# a function to implement block-coordinate descent given the current index u = alpha'X;
# return: w.a, the updated (active) coefficients; u, the updated single index; it, the number of iterations
######################################################################
block.cd <- function(w.a, ind.m, ind.old,
                     u, y, Tr, X.a, d.link.fn.obj,
                     it.max =60, eps=10^-4, trace=F,
                     nbasis.t = c(6,6,8),
                     rho.grid = c(0, 0.25, 0.5),
                     linear.link= FALSE,
                     ortho.constr = TRUE)
{

  w.a.old <- w.a + 5
  it <- 0;
  chng.fl <- F;
  w.a.start <- w.a;

  while((max(abs(w.a- w.a.old))>eps)&(max(abs(w.a- w.a.start))>eps/100)|(it==0))
  {
    w.a.old <- w.a;
    it <- it + 1
    for(j in 1:length(w.a))
    {
      gr.m <- crossprod(d.link.fn.obj$d.link.fn1, X.a[,ind.m]);
      gr1 <- crossprod(d.link.fn.obj$d.link.fn1, X.a[,j])
      tmp1 <- apply(as.matrix(X.a[,ind.old]), 2, function(v){crossprod(d.link.fn.obj$d.link.fn1, v)})
      tmp2 <- -sign(w.a[ind.old]);
      thr1 <- tmp1*tmp2
      thr2 <- -sign(w.a[ind.m])*crossprod(d.link.fn.obj$d.link.fn1, X.a[,ind.m]);
      rm(tmp1);
      rm(tmp2)
      thr.entr <- max(c(thr1,thr2))

      # it all begins when the gr1 exceeds the threshold...
      if( ( (abs(gr1)>=thr.entr) | (w.a[j]!=0) ) & (j!=ind.m) )
      {
        if(w.a[j]==0) {sn= -sign(gr1)} else{sn= sign(w.a[j])}
        gr.cmb <- crossprod(d.link.fn.obj$d.link.fn1, (X.a[,j] - sign(w.a[ind.m])*sn*X.a[,ind.m]))
        hes.cmb <- crossprod((X.a[,j] - sign(w.a[ind.m])*sn*X.a[,ind.m])^2, d.link.fn.obj$d.link.fn2)
        # compute how much it needs to be updated..
        dj <- -gr.cmb/hes.cmb
        old.w.aj <- w.a[j];
        old.w.am <- w.a[ind.m]
        w.a[j] <- w.a[j]+dj
        if(sign(w.a[j]*old.w.aj) < 0)
        {
          dj <- -old.w.aj;
          w.a[j] <- 0;
          w.a[ind.m] <- w.a[ind.m]*(abs(w.a[ind.m]) + abs(old.w.aj))/abs(w.a[ind.m]);
          dm <- w.a[ind.m]- old.w.am
        }
        else{
          dm <- -dj*sign(w.a[ind.m])*sign(w.a[j]);
          w.a[ind.m] <- w.a[ind.m] - dj*sign(w.a[ind.m])*sign(w.a[j])
        }
        if(sign(w.a[ind.m]*old.w.am)<=0)
        {
          dm <- -old.w.am;
          w.a[ind.m] <- 0;
          chng.fl <- T;
          if(trace) cat("Ch")
          if(old.w.aj!=0)
          {
            w.a[j] <- old.w.aj*(abs(old.w.aj) + abs(old.w.am))/abs(old.w.aj)
          }
          else{
            w.a[j] <- sign(w.a[j])*abs(old.w.am)
          }
          dj <- w.a[j] - old.w.aj
        }

        u <- u + dj*X.a[,j]+ dm*X.a[,ind.m]

        # update the gt and their 1st derivatives
        link.fn.obj <- fit.link.fn.gcv(y, Tr, u= u, nbasis.t = nbasis.t, rho.grid = rho.grid, linear.link = linear.link, ortho.constr = ortho.constr)
        d.link.fn.obj <- deriv.link.fn(link.fn.obj)

        if(trace) cat("wa1=",w.a,"\n");
        if(trace) cat("grads",gr1, gr.m,"\n")
        if(chng.fl)
        {
          ind.m <- which.max(abs(w.a));
          chng.fl <- F
        }
      }
    }
    if(it > it.max) break; if(trace) cat("\n","DID NOT CONVERGE","\n"); break
  }

  results <- list(w.a=w.a, u=u, it=it)
  return(results)
}


#' A wrapper function for the link function estimator
#'
#' This is a wrapper function that fits the (B-spline approximated) link functions, given the current single index u = alpha'X. An optimal smoothing parameter, \code{rho.opt}, is chosen by minimizing the generalized cross validation (GCV) for prediction errors.
#'
#' @param y  treatment outcomes, n-by-1 vector.
#' @param Tr  treatment indicators, n-by-1 vector; each element represents one of the K available treatment options.
#' @param u  a vector of the current singe index variable, u = alpha'X.
#' @param nbasis.t  a length K+1 vector; each element specifies the number of B-spline basis funtions for approximating the treatment-specific link function; the last element is for the "main effect" link function; the default is \code{nbasis.t=NULL}, and will be determined depending on the sample size.
#' @param rho.grid  a grid vector of (ridge-type) smoothing parameters for approximating the link functions.
#' @param linear.link  if \code{TRUE}, restrict the link functions to be linear funtions; the default is \code{FALSE}.
#' @param ortho.constr  the constraint that separates the interaction effects from the main effect (without this, the interaction effect can be confounded by the main effect); the default is \code{TRUE}.
#' @param ini if \code{TRUE}, calculate the MSE using the interaction-effect component only.
#'
#'@return
#' \item{beta.t.coef}{a vector of the estimated treatment-specific B-spline coefficients vectors.} \item{beta.0.coef}{a vector of the estimated B-spline coefficients fitted regardless of the treatment indicators.} \item{smoother}{the \code{smoother} object obtained from \code{smoother.fn} used in fitting the link functions.} \item{resid}{a vector of the residuals from the current fitted constrained single index model.} \item{working.resid}{a vector of the residuals from the current fitted constrained single index model, which accounts for the fitted main effects implied by  \code{beta.0.coef}.} \item{Q}{the criterion value (a MSE estimate) of the current fitted constrained single index model.} \item{rho.opt}{the optimized smoothing parameter used in fitting the link functions.}
#'
#' @seealso \code{smoother.fn},  \code{fit.link.fn},  \code{fit.csim}, \code{block.cd}
######################################################################
# a wrapper function to fit the (B-spline approximated) link functions, given the current index u = alpha'X;
# an optimal smoothing parameter, rho.opt, is chosen by minimizing GCV.
# return: beta.t.coef, beta.0.coef, smoother, resid, and MSE.
######################################################################
fit.link.fn.gcv <- function(y, Tr, u, nbasis.t = NULL, rho.grid = c(0, 0.25, 0.5), linear.link = FALSE, ortho.constr=TRUE, ini=FALSE)
{

  if(is.null(nbasis.t))
  {
    n <- length(y);
    nt <- summary(as.factor(Tr));
    K <- length(nt);
    for(t in 1:K)
    {
      nbasis.t[t] <- floor(nt[t]^{1/5.5})  + 4
    }
    nbasis.t[K+1] <- floor(n^{1/5.5}) + 6
  }

  if(length(rho.grid) >1)
  {
    smoother <-  smoother.fn(Tr, u, nbasis.t = nbasis.t, rho.grid = rho.grid, linear.link = linear.link)
    rho.opt <- fit.link.fn(y,  smoother, ortho.constr=ortho.constr, ini=ini)$rho.opt
  }else{
    rho.opt <- rho.grid
  }
  smoother <-  smoother.fn(Tr, u,  nbasis.t=nbasis.t, rho.grid = rho.opt, linear.link =linear.link)
  link.fn.obj <- fit.link.fn(y, smoother, ortho.constr=ortho.constr, ini=ini)
  return(link.fn.obj)
}



#' A smoother constructor function
#'
#' Given both a vector of the current single index variable u = alpha'X and a vector of the treatment indicators \code{Tr}, \code{smoother.fn} constructs a set of B-spline  smoother matrices, used in estimating the treatment-specific link functions (see \code{fit.link.fn}).
#'
#' The function returns a set of the QR decomposed design matrices (and the knot sequences used in constructing the B-spline design matrices), over the values of the (ridge-type)  smoothing parameters, \code{rho.grid}. Since the ridge-type smoothing is  equivalent to a regular least squares estimation with added observations, some psedo observations are added to the design matrices for the case of the nonzero values of \code{rho.grid}.
#'
#' @param Tr  treatment indicators, n-by-1 vector; each element represents one of the K available treatment options.
#' @param u  a vector of the current singe index variable, u = alpha'X.
#' @param nbasis.t  a length K+1 vector; each element specifies the number of B-spline basis funtions used in approximating the treatment-specific link function; the last element is for the "main effect" link function.
#' @param rho.grid  a grid vector of (ridge-type) smoothing parameters for approximating the link functions.
#' @param linear.link  if \code{TRUE}, restrict the link functions to be linear funtions; the default is \code{FALSE}.
#'
#' @return a list of information of the smoother matrices given the current single index variable including
#' \item{Bt.qr}{the QR decomposed design matrix constructed based on both the single index variable and the treatment indicator.} \item{B0.qr}{the QR decomposed design matrix  based on the single index variable only.} \item{Bt}{the design matrix based on both the single index variable and the treatment indicator.} \item{B0}{the design matrix  based on the single index variable only.}\item{knots.t}{the knot sequences used in constructing the treatment-specific B-spline design matrices.}
#'
#' @seealso \code{fit.link.fn},  \code{fit.link.fn.gcv},  \code{fit.csim}
######################################################################
# a subfunction to construct (B-spline) smoother matrices, given the current index u = alpha'X.
# return: the QR decomposed design matrices (and the knot sequences used in constructing the B-spline design matrices);
######################################################################
smoother.fn <- function(Tr, u, nbasis.t=c(6,6,8), rho.grid = c(0, 0.25, 0.5), linear.link = FALSE)
{

  K <- length(unique(Tr));
  u.t = design.t <- vector("list", K+1);
  # create a list, dat.list, grouped by the treatment indicator
  dat <- data.frame(Tr=Tr, u =u);
  dat_list <- dlply(dat, .(dat$Tr), function(dat) as.matrix(dat[,-1]));
  for(t in 1:K)
  {
    u.t[[t]] <-  dat_list[[t]][,1];     # data points from the tth treatment group
  }
  u.t[[K+1]] <- u;
  u.min <- max(sapply(u.t, min));
  u.max <- min(sapply(u.t, max));

  # construct treatment group-specific design matrices
  design.t = knots.t  <- vector("list", K+1)
  for(t in 1:(K+1))
  {
    if(linear.link)   # if linear.link==TRUE, construct the linear model design matrix
    {
      nbasis.t[t] <- 2
      design.t[[t]] <- cbind(1, u.t[[t]])
    }else{
      #knots.t[[t]] <- c(rep(u.min, 3),  quantile(u.t[[t]], probs = seq(0, 1, length = nbasis.t[t] -2)),  rep(u.max,3))
      knots.t[[t]] <- seq(u.min, u.max, length.out= nbasis.t[t]+4);
      design.t[[t]] <-  splineDesign(knots.t[[t]], x= u.t[[t]], outer.ok = TRUE)
    }
  }

  # construct the block-diagonal matrix consist of the treatment-specific design matrices, to approximate E(Y| u, T)
  design.t.block <- NULL;
  for(t in 1:K)
  {
    design.t.block <- c( design.t.block, list(design.t[[t]]))
  }
  Bt <- Reduce(adiag,  design.t.block)  # Bt is the block-diagonal design matrix
  # QR decomposition of the design matrix Bt, given each value of the smoothness tuning parameter, rho
  Bt.qr <- vector("list", length(rho.grid))
  ncol.Bt <- ncol(Bt)
  D <- diff(diag(ncol.Bt), differences = 2)
  for(r in seq_along(rho.grid)) # a ridge-type smoothing (equivalent to a regular least squares with added observations)
  {
    #Bt.qr[[r]] <- qr(rbind(Bt, diag(sqrt(rho.grid[r]), ncol.Bt)))
    Bt.qr[[r]] <- qr(rbind(Bt, sqrt(rho.grid[r])*D))
  }

  # compute effective degrees of freedom of smoothers, so that later we use GCV to select an optimal smoothing parameter
  edf <- vector("list", length=length(rho.grid))
  svd.Bt <- svd(Bt)
  for(r in seq_along(rho.grid))
  {
    edf[[r]] <-  sum(svd.Bt$d[svd.Bt$d>0]^2/(svd.Bt$d[svd.Bt$d>0]^2 +rho.grid[r] )) #/K
  }

  # QR decomposition of the design matrix B0
  B0 <- design.t[[K+1]]
  B0.qr  <- qr(B0)

  results <- list(Bt.qr= Bt.qr, B0.qr= B0.qr, Bt=Bt, B0= B0,
                  u.t = u.t, u.min = u.min, u.max = u.max,
                  edf= edf, rho.grid = rho.grid, K=K,
                  knots.t = knots.t, nbasis.t = nbasis.t, ncol.Bt = ncol.Bt,
                  linear.link=linear.link)
  return(results)
}





#' A link function estimator
#'
#' Given a vector of responses \code{y} and a smoother object obtained from \code{smoother.fn}, \code{fit.link.fn} computes the treatment-specific link functions of the constrained single index model.
#'
#' This is a subfunction that fits the (B-spline approximated) link functions gt, given the current single index u = alpha'X. An optimal smoothing parameter, \code{rho.opt}, is chosen by minimizing the generalized cross validation (GCV) for prediction errors.
#'
#' @param y  treatment outcomes, n-by-1 vector.
#' @param smoother  a smoother object obtained from \code{smoother.fn}.
#' @param ortho.constr  the constraint that separates the interaction effects from the main effect (without this, the interaction effect can be confounded by the main effect); the default is \code{TRUE}.
#' @param ini if \code{TRUE}, calculate the MSE using the interaction-effect component only.
#'
#' @return a list of information of the fitted link functions given the current single index variable including
#' \item{beta.t.coef}{a vector of the estimated treatment-specific B-spline coefficients vectors.} \item{beta.0.coef}{a vector of the estimated B-spline coefficients fitted regardless of the treatment indicators.} \item{smoother}{the \code{smoother} object obtained from \code{smoother.fn} used in fitting the link functions.} \item{resid}{a vector of the residuals from the current fitted constrained single index model.} \item{working.resid}{a vector of the residuals from the current fitted constrained single index model, which accounts for the fitted main effects implied by  \code{beta.0.coef}.} \item{MSE}{the criterion value (a MSE estimate) of the current fitted constrained single index model.} \item{rho.opt}{the optimized smoothing parameter used in fitting the link functions.}
#'
#' @seealso \code{smoother.fn},  \code{fit.link.fn.gcv},  \code{fit.csim}
######################################################################
# a subfunction to fit the (B-spline approximated) link functions gt, given the current index u = alpha'X.
# an optimal smoothing parameter, rho.opt, is chosen by minimizing GCV.
# return: beta.t.coef, beta.0.coef, smoother, resid, and MSE.
######################################################################
fit.link.fn <- function(y, smoother, ortho.constr = TRUE, ini=FALSE)
{

  options(warn=-1)
  # a ridge-type regularization (equivalent to an OLS with added 0s)
  y.aug <- c(y, rep(0, smoother$ncol.Bt-2))
  n <- length(y)

  # pick an optimal regularization (smoothing) paramter by GCV
  if(length(smoother$rho.grid) >1)
  {
    GCV <- numeric()
    for(s in seq_along(smoother$rho.grid))
    {
      GCV[s] <- sum((y - qr.fitted(smoother$Bt.qr[[s]], y.aug)[1:n] )^2) /(1 - smoother$edf[[s]] / n )^2
    }
    rho.index.opt <- which.min(GCV)
  }else{
    rho.index.opt <- 1
  }

  proj.Vt <- qr.fitted(smoother$Bt.qr[[rho.index.opt]], y.aug)
  if(ortho.constr)
  {
    beta.0.coef <- qr.coef(smoother$B0.qr, proj.Vt[1:n])
    proj.V0 <- drop(smoother$B0 %*% beta.0.coef)
    y.hat <- proj.Vt  - proj.V0
    beta.t.coef <- qr.coef(smoother$Bt.qr[[1]], y.hat)
  }else{
    beta.0.coef <- rep(0, ncol(smoother$B0))
    proj.V0 <- rep(0, n)
    y.hat <- proj.Vt
    beta.t.coef <- qr.coef(smoother$Bt.qr[[1]], y.hat)
  }

  working.resid <-  y - proj.Vt[1:n]
  resid <- y - y.hat[1:n]
  if(ini)
  {
    MSE  <- mean(resid^2)
  }else{
    MSE  <- mean(working.resid^2)
  }

  results <- list(MSE= MSE, resid = resid, working.resid = working.resid,
                  y.hat = y.hat[1:n],
                  working.y.hat = proj.Vt[1:n],
                  proj.V0 = proj.V0,
                  smoother = smoother,
                  beta.t.coef = beta.t.coef, beta.0.coef = beta.0.coef,
                  rho.opt = smoother$rho.grid[rho.index.opt])

  class(results) <- c("gt", "list")
  return(results)
}



#' A first derivative of the link functions estimator
#'
#' A subfunction to compute the first derivative of the link functions, evaluated at the current single index variable u = alpha'X.
#'
#' @param link.fn.obj   a fitted link function object of class \code{link.fn}, obtained from the functions \code{fit.link.fn} or \code{fit.link.fn.gcv}.
#'
#' @return a list of information of the fitted link functions given the current single index variable including
#' \item{d.link.fn}{a n x 1 vector of the 1st derivatives of the estimated link functions evaluated at the current index u = alpha'X.} \item{d.link.fn1}{an associated n x 1 vector used in computing the gradient with respect to the single index coefficients.} \item{d.link.fn2}{an associated n x 1 vector used in computing the updating rule for the single index coefficients.}
#'
#' @seealso \code{fit.link.fn},  \code{fit.link.fn.gcv},  \code{fit.csim}, \code{block.cd}
######################################################################
# a subfunction to compute the 1st derivatives of the estimated link functions, evaluated at the current index u = alpha'X.
# return: d.link.fn, a n x 1 vector of the 1st derivatives of the estimated link functions evaluated at the current index u = alpha'X; also, some other related quantities, d.gt1 and d.gt2.
######################################################################
deriv.link.fn <- function(link.fn.obj)
{

  smoother <- link.fn.obj$smoother
  K <- smoother$K
  u.t <- smoother$u.t
  knots.t <- smoother$knots.t

  t.ind  <- unlist(lapply(1:K, function(x) rep(x, smoother$nbasis.t[-(K+1)][x])))
  beta.t.coef <- split(link.fn.obj$beta.t.coef, t.ind)

  d.design.t <- vector("list", K)
  d.link.fn <- NULL
  for(t in 1:K)
  {
    if(smoother$linear.link)
    {
      d.link.fn <- c(d.link.fn, rep(beta.t.coef[[t]][2], length(u.t[[t]])) )
    }else{
      d.design.t[[t]] <- splineDesign(knots.t[[t]], x=u.t[[t]], derivs=rep(1, length(u.t[[t]])), outer.ok = TRUE)  # compute the 1st derivative of the design functions
      d.link.fn <- c(d.link.fn, d.design.t[[t]] %*% beta.t.coef[[t]])
    }
  }
  rm(d.design.t)
  d.link.fn1 <- - link.fn.obj$working.resid * d.link.fn
  d.link.fn2 <- d.link.fn1^2

  return(list(d.link.fn=d.link.fn, d.link.fn1=d.link.fn1, d.link.fn2=d.link.fn2))
}



#' A coordinate descent function
#'
#' \code{forw.cd} implements an ordinary (not a block-) coordinate descent to opitmize the single index coefficients, given an initial single index variable u = alpha'X.
#'
#' This function can be used to fit an un-regularized constrained single index model. This function is also used in obtaining the "post-selected" (re-fitted) single index coefficients, \code{cfs.pst}; see \code{fit.csim} and \code{block.cd}; see also Radchenko (2015) for the post-selected  single index coefficient estimator.
#'
#'
#' @param w.a  an initial estimate of the single index coefficient vector.
#' @param u  a vector of the current singe index variable, u = alpha'X.
#' @param y  treatment outcomes, n-by-1 vector.
#' @param Tr  treatment indicators, n-by-1 vector; each element represents one of the K available treatment options.
#' @param X.a   a matrix of the pre-treatment covarates, associated with the single index coefficients.
#' @param nbasis.t  a length K+1 vector; each element specifies the number of B-spline basis funtions used in approximating the treatment-specific link function; the last element is for the "main effect" link function.
#' @param rho.grid  a grid vector of (ridge-type) smoothing parameters for approximating the link functions.
#' @param it.max  an integer value specifying the maximum number of iterations for each coordinate.
#' @param eps  a value specifying the converge criterion of algorithm.
#' @param trace  if \code{TRUE}, show the trace of the fitting procedure; the default is \code{FALSE}.
#' @param linear.link  if \code{TRUE}, restrict the link functions to be linear functions; the default is \code{FALSE}.
#' @param ortho.constr  the constraint that separates the interaction effects from the main effect (without this, the interaction effect can be confounded by the main effect); the default is \code{TRUE}.
#' @param i.fx  a value of the index to be fixed throughout the estimation for model identifiability; the default is \code{NULL}, hence it is estimated.
#'
#' @return a list of information of the optimized single index coefficients including
#' \item{w.a}{a vector of the estimated single index coefficients.} \item{u}{a n x 1 vector of the estimated single index variable.} \item{it}{the number of iterations.}
#'
#' @seealso \code{pred.csim},  \code{fit.csim},  \code{fit.csim.cv}
######################################################################
# a function to implement an ordinary (not a block-) coordinate descent given the current index u = alpha'X;
# this function can be used to fit an un-regularized constrained single index model.
# return: w.a, the updated (active) coefficients; u, the updated single index; it, the number of iterations
######################################################################
forw.cd <- function(w.a, u, y, Tr, X.a, nbasis.t = c(6,6,8), rho.grid = c(0,0.25,0.5),
                    it.max=60, eps=10^-4, trace=F, linear.link = FALSE, ortho.constr= TRUE, i.fx=NULL)
{

  i.m <- i.fx;
  w.a.old <- w.a + 5;
  it <- 0

  while(max(abs(w.a-w.a.old)) > eps)
  {
    w.a.old <- w.a;
    it <- it + 1;
    if(is.null(i.fx))  i.m <- which.max(w.a);

    for(j in 1:length(w.a))
    {
      if(j!=i.m)
      {
        link.fn.obj <- fit.link.fn.gcv(y, Tr, u, nbasis.t = nbasis.t, rho.grid = rho.grid, linear.link = linear.link, ortho.constr=ortho.constr)
        d.link.fn.obj <- deriv.link.fn(link.fn.obj)
        gr1 <- crossprod(d.link.fn.obj$d.link.fn1, X.a[,j]);
        hes1 <- crossprod((X.a[,j])^2, d.link.fn.obj$d.link.fn2)
        dj <- -gr1/hes1;
        w.a[j] <- w.a[j] + dj;
        u <- u + dj*X.a[,j]
        if(trace) cat("w.a=",w.a,"\n")
        ####
        if(trace) cat("predictor=", j, "grad= ", gr1, "\n")
      }
    }
    if(it > it.max)
    {
      if(trace) cat("\n","DID NOT CONVERGE","\n");
      break
    }
  }

  list(w.a=w.a, u=u, it=it)
}
#######################################################################################################



######################################################################
## END OF THE FILE
######################################################################
