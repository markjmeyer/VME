#### Simulation Study ####
# Vary:
#     n  = 20, 50, 100
#     nc = 2, 3, 4, 5
#     beta0 = 25, beta1 = -5, beta2 = -1
#     models = random intercept vs random intercept + slope
#

#### Source Code ####
source('vme_v0-2.R')
library(mvnfast)

#### Sim Specs ####

## change ##
n         <- 20 # 50, 100
B         <- 500
xi0       <- 0.99
ncVec     <- 2:5
maxIterS  <- 20000
dots      <- 50
up        <- 1000

sigu      <- 1
sige      <- 1
Sigb      <- matrix(c(1, 0.5*1,
                      0.5*1, 1), nrow = 2, byrow = TRUE)

# coefficients #
b0        <- 25
b1        <- -5
b2        <- -1

#### Storage Matrices ####
# '_b' = true model is intercept + slope
# '_i' = true model is intercept only

vaicMat_b   <- array(NA, dim = c(length(ncVec), B, 2))
vaicMat_i   <- array(NA, dim = c(length(ncVec), B, 2))

biasMat_b   <- array(NA, dim = c(length(ncVec), B, 3, 2))
biasMat_i   <- array(NA, dim = c(length(ncVec), B, 3, 2))

miseMat_b   <- array(NA, dim = c(length(ncVec), B, 3, 2))
miseMat_i   <- array(NA, dim = c(length(ncVec), B, 3, 2))


#### Run ####

iter <- 0

for(nc in ncVec){
  
  ##### data #####
  ID        <- rep(1:n, each = nc)
  X         <- c(rep(0,n*nc/2), rep(1, n*nc/2))
  time      <- 2*rep(0:(nc-1), n)
  
  ##### generate datasets and run models #####
  for(b in 1:B){
    set.seed(b)
    
    ##### generate error terms and matrices #####
    e       <- rnorm(nc*n, sd = sige)
    u       <- rnorm(n, sd = sigu)
    # s       <- rnorm(n, sd = sigu)
    us      <- rmvn(n, mu = c(0, 0), sigma = Sigb)
    
    U       <- model.matrix(~ as.factor(ID) - 1)
    S       <- model.matrix(~time*as.factor(ID) - 1 - as.factor(ID) - time)
    USList   <- vector('list', length = n)
    
    for(i in 1:n){
      USi           <- matrix(c(rep(1,nc), 2*(0:(nc-1))), nrow = nc)
      USList[[i]]   <- USi
    }
    
    US  <- as.matrix(bdiag(USList))
    
    ##### construct outcomes #####
    Yu      <- b0 + b1*X + b2*time + U%*%u + e
    Yus     <- b0 + b1*X + b2*time + US%*%c(t(us)) + e
    
    ##### build data frames #####
    datai   <- data.frame(Y = Yu, X = X, time = time, ID = ID)
    datais  <- data.frame(Y = Yus, X = X, time = time, ID = ID)
    
    ##### random intercept only #####
    model_is_ru  <- vme(fixed = Y ~ X + time, Z = US, rint.only = FALSE, fe = NULL,
                        data = datai, ID = 'ID',
                        ps.spline = list(df = 5, degree = 3, xi = 0.001),
                        re.spline = list(df = NULL, degree = 1, xi = xi0),
                        priors = list(Ae = 0.01, Be = 0.01, Af = 0.01, Bf = 0.01, 
                                      Au = 0.01, Bu = 0.01, sigB = 1000),
                        tol = 1e-5, maxIter = maxIterS, up = 100,
                        verbose = FALSE)
    
    model_i_ru    <- vme(fixed = Y ~ X + time, Z = NULL, rint.only = TRUE, fe = NULL,
                         data = datai, ID = 'ID',
                         ps.spline = list(df = 5, degree = 3, xi = 0.001),
                         re.spline = list(df = NULL, degree = 1, xi = xi0),
                         priors = list(Ae = 0.01, Be = 0.01, Af = 0.01, Bf = 0.01, 
                                       Au = 0.01, Bu = 0.01, sigB = 1000),
                         tol = 1e-5, maxIter = maxIterS, up = 100,
                         verbose = FALSE)

    ##### random intercept + slope #####
    model_is_rus  <- vme(fixed = Y ~ X + time, Z = US, rint.only = FALSE, fe = NULL,
                         data = datais, ID = 'ID',
                         ps.spline = list(df = 5, degree = 3, xi = 0.001),
                         re.spline = list(df = NULL, degree = 1, xi = xi0),
                         priors = list(Ae = 0.01, Be = 0.01, Af = 0.01, Bf = 0.01, 
                                       Au = 0.01, Bu = 0.01, sigB = 1000),
                         tol = 1e-5, maxIter = maxIterS, up = 100,
                         verbose = FALSE)
    
    model_i_rus   <- vme(fixed = Y ~ X + time, Z = NULL, rint.only = TRUE, fe = NULL,
                         data = datais, ID = 'ID',
                         ps.spline = list(df = 5, degree = 3, xi = 0.001),
                         re.spline = list(df = NULL, degree = 1, xi = xi0),
                         priors = list(Ae = 0.01, Be = 0.01, Af = 0.01, Bf = 0.01, 
                                       Au = 0.01, Bu = 0.01, sigB = 1000),
                         tol = 1e-5, maxIter = maxIterS, up = 100,
                         verbose = FALSE)
    
    #### VAIC ####
    ni                      <- which(ncVec == nc)
    vaicMat_i[ni, b, 1]     <- VAIC(model_is_ru)
    vaicMat_i[ni, b, 2]     <- VAIC(model_i_ru) # True Model
    
    vaicMat_b[ni, b, 1]     <- VAIC(model_is_rus) # True Model
    vaicMat_b[ni, b, 2]     <- VAIC(model_i_rus)
    
    biasMat_i[ni, b, , 1]   <- coef(model_is_ru) - c(b0, b1, b2)
    biasMat_i[ni, b, , 2]   <- coef(model_i_ru) - c(b0, b1, b2) # True Model
    
    biasMat_b[ni, b, , 1]   <- coef(model_is_rus) - c(b0, b1, b2) # True Model
    biasMat_b[ni, b, , 2]   <- coef(model_i_rus) - c(b0, b1, b2)

    miseMat_i[ni, b, , 1]   <- (coef(model_is_ru) - c(b0, b1, b2))^2
    miseMat_i[ni, b, , 2]   <- (coef(model_i_ru) - c(b0, b1, b2))^2 # True Model
    
    miseMat_b[ni, b, , 1]   <- (coef(model_is_rus) - c(b0, b1, b2))^2 # True Model
    miseMat_b[ni, b, , 2]   <- (coef(model_i_rus) - c(b0, b1, b2))^2

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
# pname   <- paste('~InterSlope', '/N', n, sep = '')
# setwd(pname)

saveRDS(vaicMat_i, 'vaicMat_i.rds')
saveRDS(vaicMat_b, 'vaicMat_b.rds')

saveRDS(biasMat_i, 'biasMat_i.rds')
saveRDS(biasMat_b, 'biasMat_b.rds')

saveRDS(miseMat_i, 'miseMat_i.rds')
saveRDS(miseMat_b, 'miseMat_b.rds')
