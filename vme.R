#### Variational Mixed Effects Models ####
# Version:  0.3                          #
# Modified: 07/27/2023                   #

#### require libraries ####
suppressPackageStartupMessages(require(Matrix))
suppressPackageStartupMessages(require(splines))
suppressPackageStartupMessages(require(MASS))
suppressPackageStartupMessages(require(numbers))

#### penalty matrix ####

##### B-spline penalty matrix #####
formPmat  <- function(D, xi, Theta){
  diff0   <- diag(1, D, D)
  diff2   <- matrix(rep(c(1, -2, 1, rep(0, D - 2)), D - 2)[1:((D - 2) * D)], D - 2, D, 
                    byrow = TRUE)
  P0      <- t(Theta) %*% t(diff0) %*% diff0 %*% Theta
  P2      <- t(Theta) %*% t(diff2) %*% diff2 %*% Theta
  P.mat   <- xi * P0 + (1 - xi) * P2
  
  return(P.mat)
}


#### Gaussian case ####
##### KL lower bound for Gaussian VME #####
logp_gauss  <- function(Sigmaq, sigB, Muq, p, Sigqb, Ae, Be, n, Bqe, Au, Bu, Ku, Bqu, 
                        Af= NULL, Bf= NULL, Kf= NULL, Bqf= NULL, fe = NULL){
  if(is.null(fe)){
    if(n/2 < 150){
      out   <- (1/2)*(p + Ku) - (n/2)*log(2*pi) - (p/2)*log(sigB) +
        (1/2)*log(det(Sigmaq)) - (1/(2*sigB))*(t(Muq[1:p])%*%Muq[1:p] + sum(diag(Sigqb))) -
        Ae*log(Be) - (Ae + n/2)*log(Bqe) + log(gamma(Ae + n/2)) - log(gamma(Ae)) +
        Au*log(Bu) - (Au + Ku/2)*log(Bqu) + log(gamma(Au + Ku/2)) - log(gamma(Au))
    } else{
      out   <- (1/2)*(p + Ku) - (n/2)*log(2*pi) - (p/2)*log(sigB) +
        (1/2)*log(det(Sigmaq)) - (1/(2*sigB))*(t(Muq[1:p])%*%Muq[1:p] + sum(diag(Sigqb))) -
        Ae*log(Be) - (Ae + n/2)*log(Bqe) - log(gamma(Ae)) +
        Au*log(Bu) - (Au + Ku/2)*log(Bqu) + log(gamma(Au + Ku/2)) - log(gamma(Au))
    }
  } else{
    if(n/2 < 150){
      out   <- (1/2)*(p + Kf + Ku) - (n/2)*log(2*pi) - (p/2)*log(sigB) +
        (1/2)*log(det(Sigmaq)) - (1/(2*sigB))*(t(Muq[1:p])%*%Muq[1:p] + sum(diag(Sigqb))) -
        Ae*log(Be) - (Ae + n/2)*log(Bqe) + log(gamma(Ae + n/2)) - log(gamma(Ae)) +
        Af*log(Bf) - (Af + Kf/2)*log(Bqf) + log(gamma(Af + Kf/2)) - log(gamma(Af)) +
        Au*log(Bu) - (Au + Ku/2)*log(Bqu) + log(gamma(Au + Ku/2)) - log(gamma(Au))
    } else{
      out   <- (1/2)*(p + Kf + Ku) - (n/2)*log(2*pi) - (p/2)*log(sigB) +
        (1/2)*log(det(Sigmaq)) - (1/(2*sigB))*(t(Muq[1:p])%*%Muq[1:p] + sum(diag(Sigqb))) -
        Ae*log(Be) - (Ae + n/2)*log(Bqe) - log(gamma(Ae)) +
        Af*log(Bf) - (Af + Kf/2)*log(Bqf) + log(gamma(Af + Kf/2)) - log(gamma(Af)) +
        Au*log(Bu) - (Au + Ku/2)*log(Bqu) + log(gamma(Au + Ku/2)) - log(gamma(Au))
    }
  }
  
  return(as.numeric(out))
}



##### log(p(y|theta*)) for Gaussian VME #####
loglikst_gauss <- function(n, Bqe, Ae, Y, XZ, Muq){
  out <- -n/2 - (n/2)*(log(Bqe) - log(Ae + n/2 - 1)) - (1/2)*((Ae + n/2 - 1)/Bqe)*t(Y - XZ%*%Muq)%*%(Y - XZ%*%Muq)
  
  return(as.numeric(out))
}

##### p*D for Gaussian VME #####
eqloglik_gauss <- function(n, Ae, Bqe, XZ, Sigmaq, Y, Muq){
  XtSX <- XZ%*%Sigmaq%*%t(XZ)
  out <- -n/2*log(2*pi) + (n/2)*(digamma(Ae + n/2) - log(Bqe)) - 
    (1/2)*((Ae + n/2)/Bqe)*(sum(diag(XtSX)) + t(Y - XZ%*%Muq)%*%(Y - XZ%*%Muq))
  
  return(as.numeric(out))
}

##### Gaussian VAIC #####
VAIC_vme <- function(n, Ae, Bqe, XZ, Sigmaq, Y, Muq){
  Pst_D   <- 2*loglikst_gauss(n, Bqe, Ae, Y, XZ, Muq) - 2*eqloglik_gauss(n, Ae, Bqe, XZ, Sigmaq, Y, Muq)
  out     <- -2*loglikst_gauss(n, Bqe, Ae, Y, XZ, Muq) + 2*Pst_D
  
  return(as.numeric(out))
}

##### Gaussian VME #####
vme   <- function(fixed, Z = NULL, rint.only = TRUE, fe = NULL, data, ID,
                  ps.spline = list(df = 5, degree = 3, xi = 0.001),
                  re.spline = list(df = NULL, degree = 3, xi = 0.001),
                  priors = list(Ae = 0.01, Be = 0.01, Af = 0.01, Bf = 0.01, 
                                Au = 0.01, Bu = 0.01, sigB = 1000),
                  tol = 1e-5, maxIter = 500, up = 100, dots = up/10,
                  verbose = TRUE){
  
  ## extract data, build matrices and vectors
  mf      <- model.frame(formula = fixed, data = data)
  Y       <- model.response(mf)
  Xc      <- model.matrix(attr(mf,"terms"), data = mf)
  if(!is.null(fe)){
    mfe     <- as.formula(paste(deparse(fe),"- 1"))
    Zl      <- model.matrix(mfe)
    Ks      <- ncol(Zl)
    x       <- 1:Ks
    dfs     <- ps.spline$df
    degs    <- ps.spline$degree
    xi      <- ps.spline$xi
    Th      <- bs(x, df = dfs, intercept = TRUE, degree = degs)
    O       <- formPmat(Ks, xi, Th)
    Tmat    <- list(Th = Th)
    Zf      <- Zl%*%Th
  }
  n       <- length(Y) # total number of samples, ncol(Zi) is number of subjects
  
  if(rint.only == FALSE){
    if(is.null(Z)){
      stop('User must provide general random effects design matrix')
    }
    Kz      <- ncol(Z)
    
    if(is.null(re.spline$df)){
      redf  <- length(unique(data[,ID])) # one knot per subject
    } else {
      redf  <- re.spline$df
    }
    redeg   <- re.spline$degree
    rexi    <- re.spline$xi
    Tz      <- bs(1:Kz, df = redf, intercept = TRUE, degree = redeg)
    Pz      <- formPmat(Kz, rexi, Tz)
    
    ZT      <- Z%*%Tz
    
    if(is.null(fe)){
      Tmat    <- list(Tz = Tz)
    } else{
      Tmat    <- list(Th = Th, Tz = Tz)
    }
  } else{
    id    <- data[,ID]
    Z     <- model.matrix(~ as.factor(id) - 1)
  }
  
  ### build design matrix ###
  if(rint.only == TRUE){
    if(is.null(fe)){
      XZ      <- cbind(Xc, Z)
    } else {
      XZ      <- cbind(Xc, Zf, Z)
    }
  } else {
    if(is.null(fe)){
      XZ      <- cbind(Xc, ZT)
    } else {
      XZ      <- cbind(Xc, Zf, ZT)
    }
  }
  
  ## prior values ##
  Ae    <- priors$Ae
  Be    <- priors$Be
  Af    <- priors$Af
  Bf    <- priors$Bf
  Au    <- priors$Au
  Bu    <- priors$Bu
  sigB  <- priors$sigB
  
  ## initialize algorithm ##
  Bqe         <- 1
  Bqf         <- 1
  Bqu         <- 1
  logp_delta	<- 0.5
  iter		    <- 0
  logp_prev	  <- 0
  
  ## set matrices ##
  CtC         <- t(XZ)%*%XZ
  CtY         <- t(XZ)%*%Y
  p           <- ncol(Xc)
  if(!is.null(fe)){
    Kf         <- ncol(Zf)
  } 
  if(rint.only == TRUE){
    Ku          <- ncol(Z)
  } else {
    Ku          <- ncol(ZT)
  }
  if(is.null(fe)){
    MuMat       <- matrix(0, nrow = maxIter, ncol = p + Ku)
    SigmaMat    <- array(0, dim = c(p + Ku, p + Ku, maxIter))
  } else{
    MuMat       <- matrix(0, nrow = maxIter, ncol = p + Kf + Ku)
    SigmaMat    <- array(0, dim = c(p + Kf + Ku, p + Kf + Ku, maxIter))
  }
  BqeVec      <- rep(0, maxIter)
  BqfVec      <- rep(0, maxIter)
  BquVec      <- rep(0, maxIter)
  lDelta      <- rep(logp_delta, maxIter)
  logpi       <- rep(0, maxIter)
  logpp       <- rep(0, maxIter)
  
  ## algorithm ##
  while(logp_delta > tol){
    iter	<- iter + 1
    
    ## update \Sigma_q ##
    Sigqb   <- (1/sigB)*diag(p)
    if(!is.null(fe)){
      Sigqf   <- ((Af + (1/2)*Kf)/Bqf)*O
    }
    if(rint.only == TRUE){
      Sigqu   <- ((Au + (1/2)*Ku)/Bqu)*diag(Ku) # random intercept variance
    } else {
      Sigqu   <- ((Au + (1/2)*Ku)/Bqu)*Pz # random function variance
    }
    if(is.null(fe)){
      G       <- bdiag(Sigqb, Sigqu)
    } else{
      G       <- bdiag(Sigqb, Sigqf, Sigqu)
    }
    Precq   <- as.matrix(((Ae + n/2)/Bqe)*CtC + G)
    Sigmaq  <- solve(Precq)
    
    ## update \mu_q ##
    Muq     <- (Ae + n/2)/(Bqe)*Sigmaq%*%CtY
    Ctf     <- (Ae + n/2)/(Bqe)*Sigmaq%*%t(XZ)
    
    ## update B_qe ##
    SSY     <- as.numeric(t(Y - XZ%*%Muq)%*%(Y - XZ%*%Muq))
    Bqe     <- Be + (1/2)*(SSY + sum(diag(CtC%*%Sigmaq)))
    
    ## update B_qf ##
    if(!is.null(fe)){
      MltMf   <- as.numeric(t(Muq[(p + 1):(p + Kf)])%*%Muq[(p + 1):(p + Kf)])
      Bqf    <- Bf + (1/2)*(MltMf + sum(diag(Sigqf)))
    }
    
    ## update B_qu ##
    if(is.null(fe)){
      MutMu   <- as.numeric(t(Muq[-c(1:p)])%*%Muq[-c(1:p)])
      Bqu     <- Bu + (1/2)*(MutMu + sum(diag(Sigqu)))
    } else{
      MutMu   <- as.numeric(t(Muq[-c(1:(p + Kf))])%*%Muq[-c(1:(p + Kf))])
      Bqu     <- Bu + (1/2)*(MutMu + sum(diag(Sigqu)))
    }
    
    ## check criteria ##
    if(is.null(fe)){
      logp_i      <- logp_gauss(Sigmaq = Sigmaq, sigB = sigB, Muq = Muq, p = p, 
                                Sigqb = Sigqb, Ae = Ae, Be = Be, n = n, Bqe = Bqe,
                                Au = Au, Bu = Bu, Ku = Ku, Bqu = Bqu)
    } else{
      logp_i      <- logp_gauss(Sigmaq = Sigmaq, sigB = sigB, Muq = Muq, p = p, 
                                Sigqb = Sigqb, Ae = Ae, Be = Be, n = n, Bqe = Bqe,
                                Au = Au, Bu = Bu, Ku = Ku, Bqu = Bqu, 
                                Af = Af, Bf = Bf, Kf = Kf, Bqf = Bqf, fe = fe)
    }
    logp_delta	<- abs(logp_i - logp_prev)
    
    ## track convergence ##  
    lDelta[iter]  <- logp_delta
    logpi[iter]   <- logp_i
    logpp[iter]   <- logp_prev
    logp_prev	    <- logp_i
    
    ## store estimates ##
    MuMat[iter,]      <- as.numeric(Muq)
    SigmaMat[,,iter]  <- as.matrix(Sigmaq)
    BqeVec[iter]      <- Bqe
    BqfVec[iter]      <- Bqf
    BquVec[iter]      <- Bqu
    
    ## console update ##
    if(verbose){
      if(mod(iter, dots) == 0){
        cat('.')
      }
      
      if(mod(iter, up) == 0){
        cat(paste("\n",iter,"samples completed\n"))
      }
    }
    
    if(iter == maxIter){
      warning('Max number of iterations reached; increase max'); break
    }
    
  }
  
  ## algo outs ##
  algOut  <- list(maxIter = maxIter, tol = tol, totalIter = iter, logpd = list(lDelta = lDelta, logp = logpi, prev_logp = logpp))
  
  ## outputs ##
  if(is.null(fe)){
    ## fixed effects ##
    xcoef     <- matrix(MuMat[iter,1:p], p, 1)
    colnames(xcoef) <- "Estimate"
    rownames(xcoef) <- colnames(Xc)
    
    ## random effects ##
    if(rint.only == TRUE){
      U         <- MuMat[iter,(p + 1):ncol(MuMat)]
    } else {
      U         <- Tz%*%MuMat[iter,(p + 1):ncol(MuMat)]
    }
    
    ## design matrix ##
    if(rint.only == TRUE){
      dmat <- list(XZ = XZ, Xc = Xc, Z = Z)
    } else {
      dmat <- list(XZ = XZ, Xc = Xc, ZT = ZT)
    }
    
    ## VAIC ##
    vaic        <- VAIC_vme(n, Ae, Bqe, XZ, Sigmaq, Y, Muq)
    
    ## outs ##
    out       <- list(VBEstimates = MuMat[iter,],
                      fixedEffects = xcoef,
                      randomEffects = U,
                      varComps = list(Sigmaq = SigmaMat[,,iter], Bqe = BqeVec[iter], 
                                      Bqf = BqfVec[iter], Bqu = BquVec[iter]),
                      VAIC = vaic,
                      model = list(fixed = fixed, Z = Z, rint.only = rint.only, fe = fe, data = data, ID = ID,
                                   p = p, Ku = Ku, family = 'Gaussian',
                                   dmat = dmat, Y = Y),
                      priors = priors,
                      algOut = algOut)
  } else{
    ## functional effect ##
    fee       <- as.numeric(Th%*%MuMat[iter,(p + 1):(p + Kf)])
    feevar    <- as.matrix(Th%*%SigmaMat[(p + 1):(p + Kf),(p + 1):(p + Kf),iter]%*%t(Th))
    
    ## fixed effects ##
    xcoef     <- matrix(MuMat[iter,1:p], p, 1)
    colnames(xcoef) <- "Estimate"
    rownames(xcoef) <- colnames(Xc)
    
    ## random effects ##
    if(rint.only == TRUE){
      U         <- MuMat[iter,(p + Kf + 1):ncol(MuMat)]
    } else {
      U         <- Tz%*%MuMat[iter,(p + Kf + 1):ncol(MuMat)]
    }
    
    ## design matrix ##
    if(rint.only == TRUE){
      dmat <- list(XZ = XZ, Xc = Xc, Zf = Zf, Z = Z)
    } else {
      dmat <- list(XZ = XZ, Xc = Xc, Zf = Zf, ZT = ZT)
    }
    
    ## VAIC ##
    vaic        <- VAIC_vme(n, Ae, Bqe, XZ, Sigmaq, Y, Muq)
    
    ## outs ##
    out       <- list(VBEstimates = MuMat[iter,],
                      fixedEffects = xcoef,
                      funcEffect = fee,
                      randomEffects = U,
                      varComps = list(Sigmaq = SigmaMat[,,iter], feVar = feevar,
                                      Bqe = BqeVec[iter], Bqf = BqfVec[iter], 
                                      Bqu = BquVec[iter]),
                      VAIC = vaic,
                      model = list(fixed = fixed, Z = Z, rint.only = rint.only, fe = fe, data = data, ID = ID,
                                   p = p, Kf = Kf, Ku = Ku, O = O, Tmat = Tmat, family = 'Gaussian',
                                   dmat = dmat, Y = Y),
                      priors = priors,
                      algOut = algOut)
  }
  
  class(out) <- 'vme'
  return(out)
  
}


#### class functions ####

VAIC <- function(object, ...){
  UseMethod("VAIC")
}

VAIC.vme <- function(fit){
  return(fit$VAIC)
}

plot.vme <- function(fit, scale = c('lograte', 'rate'), ...){
  if(is.null(fit$model$fe)){
    stop('No plotting for VMEs with no functional effects')
  }
  EX    <- fit$model$exposure
  family    <- fit$model$family
  if(length(scale) > 1){
    scale <- 'lograte'
  }
  if(is.null(EX)){
    lags <- fit$lagedEffects
    switch(family,
           Gaussian = plot(rev(lags), type = 'l', ...),
           Poisson = switch(scale, 
                            rate = plot(exp(rev(lags)), type = 'l', ...), 
                            lograte = plot(rev(lags), type = 'l', ...))
    )
  } else{
    diffs <- fit$lagedEffects$diff
    switch(family,
           Gaussian = plot(rev(diffs), type = 'l', ...),
           Poisson = switch(scale, 
                            rate = plot(exp(rev(diffs)), type = 'l', ...), 
                            lograte = plot(rev(diffs), type = 'l', ...))
    )
  }
}

summary.vme <- function(fit){
  family    <- fit$model$family
  fixedEf   <- fit$model$fixed
  vaic      <- fit$VAIC
  totalIter <- fit$algOut$totalIter
  maxIter   <- fit$algOut$maxIter
  tol       <- fit$algOut$tol
  cat(paste(family, "Variational Mixed Effects Model\n"))
  if(fit$model$rint.only){
    cat(paste("with Random Intercept\n"))
  } else{
    cat(paste("with General Random Function\n"))
  }
  if(totalIter < maxIter){
    cat(paste("Converged on", totalIter, "iterations, tolerance of", tol, "\n"))
    cat(paste("VAIC =", round(vaic,2), "\n\n"))
    cat(paste("Fixed effects:", deparse(fixedEf), "\n\n"))
    coeft <- t(fit$fixedEffects)
    rownames(coeft) <- ''
    print(round(coeft,4))
  } else{
    cat(paste("FAILED to converge on", totalIter, "iterations, tolerance of", tol, "\n"))
    cat("Increase number of iterations")
  }
}

print.vme <- function(fit){
  family    <- fit$model$family
  fixedEf   <- fit$model$fixed
  totalIter <- fit$algOut$totalIter
  maxIter   <- fit$algOut$maxIter
  cat(paste(family, "Variational Mixed Effects Model\n"))
  if(fit$model$rint.only){
    cat(paste("with Random Intercept\n"))
  } else{
    cat(paste("with General Random Function\n"))
  }
  if(totalIter < maxIter){
    cat(paste("Fixed effects:", deparse(fixedEf), "\n\n"))
    coeft <- t(fit$fixedEffects)
    rownames(coeft) <- ''
    print(round(coeft,4))
  } else{
    cat(paste("FAILED to converge on", totalIter, "iterations, tolerance of", tol, "\n"))
    cat("Increase number of iterations")
  }  
}

coef.vme <- function(fit){
  return(as.numeric(fit$fixedEffects))
}



