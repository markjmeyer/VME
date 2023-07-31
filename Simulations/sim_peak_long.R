#### Simulation Study ####
# Vary:
#     n  = 20, 50, 100, T = 50, 100
#     nc = 2, 3, 4, 5
#     beta = Cyclical, *Peak*, Sigmoidal
#     models = random function, random intercept
#

#### Source Code ####
source('vme_v0-2.R')
library(mvnfast)

#### Sim Specs ####

## change ##
n         <- 100 # 20, 50, 100
tT        <- 100 # 50, 100

## set-up ##
B         <- 500
ncVec     <- 2:5
maxIterS  <- 20000
xi0       <- 0.01
df0       <- 8
dots      <- 50
up        <- 1000

sigX2     <- 1
rho       <- 0.5

sigg      <- 1
rhog      <- 0.75

sigu      <- 1
sige      <- 1

tv        <- seq(0, 1, length.out = tT)

## lag variance ##
ar1Corr   <- diag(tT)
for(i in 1:(tT-1)){
  for(j in (i+1):tT){
    ar1Corr[i,j] <- rho^(j-i)
    ar1Corr[j,i] <- rho^(j-i)
  }
}
SigSim  <- sigX2*ar1Corr

## functional effect ##
g0        <- 0.25*dnorm(tv, mean = 0.75, sd = 0.1)

# intercept #
b0        <- 50

#### Storage Matrices ####
# '_f' = true model is function
# '_i' = true model is intercept only

vaicMat_f   <- array(NA, dim = c(length(ncVec), B, 2))
vaicMat_i   <- array(NA, dim = c(length(ncVec), B, 2))

biasMat_f   <- array(NA, dim = c(length(ncVec), B, 2))
biasMat_i   <- array(NA, dim = c(length(ncVec), B, 2))

miseMat_f   <- array(NA, dim = c(length(ncVec), B, 2))
miseMat_i   <- array(NA, dim = c(length(ncVec), B, 2))

#### Run ####

iter <- 0

for(nc in ncVec){
  
  ### data ###
  ID        <- rep(1:n, each = nc)
  
  ## functions ## 
  set.seed(2019)
  Fe        <- rmvn(n*nc, mu = rep(0, tT), sigma = SigSim)
  
  #### generate datasets and run models ####
  for(b in 1:B){
    set.seed(b)
    
    ## random effects ##
    
    # random function #
    g       <- rnorm(n, mean = 0, sd = sigg)
    
    UF     <- matrix(0, nrow = nc*n, ncol = tT*n)
    rind    <- matrix(1:(nc*n), nrow = nc)
    cind    <- matrix(1:(n*ncol(Fe)), nrow = ncol(Fe))
    
    for(i in 1:n){
      ri          <- rind[,i]
      ci          <- cind[,i]
      lsi         <- Fe[ri,]
      UF[ri,ci]  <- lsi
    }
    Tg      <- bs(1:ncol(UF), df = n, intercept = FALSE, degree = 3)
    
    # random intercept #
    u       <- rnorm(n, sd = sigu)
    U       <- model.matrix(~ as.factor(ID) - 1)
    
    # model error #
    e       <- rnorm(nc*n, sd = sige)
    
    ## outcome ##
    Yf      <- b0 + Fe%*%g0 + UF%*%Tg%*%g + e
    Yi      <- b0 + Fe%*%g0 + U%*%u + e
    
    ## data ##
    datai   <- data.frame(Y = Yi, ID = ID)
    dataf   <- data.frame(Y = Yf, ID = ID)
    
    #### random intercept model ####
    model_f_ri  <- vme(Y ~ 1, Z = UF, rint.only = FALSE, fe = ~ Fe, data = datai, ID = 'ID',
                       ps.spline = list(df = df0, degree = 3, xi = xi0),
                       re.spline = list(df = NULL, degree = 3, xi = xi0),
                       priors = list(Ae = 0.01, Be = 0.01, Af = 0.01, Bf = 0.01,
                                     Au = 0.01, Bu = 0.01, sigB = 1000),
                       tol = 1e-5, maxIter = maxIterS, up = 100,
                       verbose = FALSE)
    
    model_i_ri  <- vme(Y ~ 1, Z = NULL, rint.only = TRUE, fe = ~ Fe, data = datai, ID = 'ID',
                       ps.spline = list(df = df0, degree = 3, xi = xi0),
                       re.spline = list(df = NULL, degree = 3, xi = xi0),
                       priors = list(Ae = 0.01, Be = 0.01, Af = 0.01, Bf = 0.01,
                                     Au = 0.01, Bu = 0.01, sigB = 1000),
                       tol = 1e-5, maxIter = maxIterS, up = 100,
                       verbose = FALSE)
    
    #### random function model ####
    model_f_rf  <- vme(Y ~ 1, Z = UF, rint.only = FALSE, fe = ~ Fe, data = dataf, ID = 'ID',
                       ps.spline = list(df = df0, degree = 3, xi = xi0),
                       re.spline = list(df = NULL, degree = 3, xi = xi0),
                       priors = list(Ae = 0.01, Be = 0.01, Af = 0.01, Bf = 0.01,
                                     Au = 0.01, Bu = 0.01, sigB = 1000),
                       tol = 1e-5, maxIter = maxIterS, up = 100,
                       verbose = FALSE)
    
    
    model_i_rf  <- vme(Y ~ 1, Z = NULL, rint.only = TRUE, fe = ~ Fe, data = dataf, ID = 'ID',
                       ps.spline = list(df = df0, degree = 3, xi = xi0),
                       re.spline = list(df = NULL, degree = 3, xi = xi0),
                       priors = list(Ae = 0.01, Be = 0.01, Af = 0.01, Bf = 0.01,
                                     Au = 0.01, Bu = 0.01, sigB = 1000),
                       tol = 1e-5, maxIter = maxIterS, up = 100,
                       verbose = FALSE)
    
    #### VAIC ####
    ni                    <- which(ncVec == nc)
    vaicMat_i[ni, b, 1]   <- VAIC(model_f_ri)
    vaicMat_i[ni, b, 2]   <- VAIC(model_i_ri) # True Model
    
    vaicMat_f[ni, b, 1]   <- VAIC(model_f_rf) # True Model
    vaicMat_f[ni, b, 2]   <- VAIC(model_i_rf)
    
    #### Bias ####
    biasMat_i[ni, b, 1]   <- mean(c(model_f_ri$fixedEffects - b0, model_f_ri$funcEffect - g0))
    biasMat_i[ni, b, 2]   <- mean(c(model_i_ri$fixedEffects - b0, model_i_ri$funcEffect - g0)) # true model
    
    biasMat_f[ni, b, 1]   <- mean(c(model_f_rf$fixedEffects - b0, model_f_rf$funcEffect - g0)) # true model
    biasMat_f[ni, b, 2]   <- mean(c(model_i_rf$fixedEffects - b0, model_i_rf$funcEffect - g0))
    
    #### MISE ####
    miseMat_i[ni, b, 1]   <- mean(c(model_f_ri$fixedEffects - b0, model_f_ri$funcEffect - g0)^2)
    miseMat_i[ni, b, 2]   <- mean(c(model_i_ri$fixedEffects - b0, model_i_ri$funcEffect - g0)^2) # true model
    
    miseMat_f[ni, b, 1]   <- mean(c(model_f_rf$fixedEffects - b0, model_f_rf$funcEffect - g0)^2) # true model
    miseMat_f[ni, b, 2]   <- mean(c(model_i_rf$fixedEffects - b0, model_i_rf$funcEffect - g0)^2)
    
    #### simulation controls ####
    iter <- iter + 1
    
    if(mod(iter, dots) == 0){
      cat('.')
    }
    
    if(mod(iter, up) == 0){
      cat(paste("\n",iter,"datasets completed\n"))
    }
    
  }
}


#### Save Output ####
# pname   <- paste('~Functional/Peak/N', n, '/T', tT, sep = '')
# setwd(pname)

saveRDS(vaicMat_i, 'vaicMat_i.rds')
saveRDS(vaicMat_f, 'vaicMat_f.rds')

saveRDS(biasMat_i, 'biasMat_i.rds')
saveRDS(biasMat_f, 'biasMat_f.rds')

saveRDS(miseMat_i, 'miseMat_i.rds')
saveRDS(miseMat_f, 'miseMat_f.rds')










