#### Data load ####

lead <- read.table('~leadlong.txt', header = TRUE)

#### Treatment of Lead Exposure in Children ####
source('vme.R')

##### random intercept model #####
model_ri  <- vme(fixed = level ~ treatment + time, Z = NULL, rint.only = TRUE, fe = NULL,
                 data = lead, ID = 'ID',
                 ps.spline = list(df = 5, degree = 3, xi = 0.001),
                 re.spline = list(df = NULL, degree = 1, xi = 0.99),
                 priors = list(Ae = 0.01, Be = 0.01, Af = 0.01, Bf = 0.01, 
                               Au = 0.01, Bu = 0.01, sigB = 1000),
                 tol = 1e-5, maxIter = 500, up = 100,
                 verbose = TRUE)

##### random intercept and slope on time #####
Zrs   <- model.matrix(~as.factor(ID) + time*as.factor(ID) - 1 - time, data = lead)

model_ris   <- vme(fixed = level ~ treatment + time, Z = Zrs, rint.only = FALSE, fe = NULL,
                   data = lead, ID = 'ID',
                   ps.spline = list(df = 5, degree = 3, xi = 0.001),
                   re.spline = list(df = NULL, degree = 1, xi = 0.99),
                   priors = list(Ae = 0.01, Be = 0.01, Af = 0.01, Bf = 0.01,
                                 Au = 0.01, Bu = 0.01, sigB = 1000),
                   tol = 1e-5, maxIter = 2000, up = 100,
                   verbose = TRUE)

VAIC(model_ri) - 
VAIC(model_ris)

summary(model_ris)

# pdf('lead_logp.pdf')
plot(1:model_ris$algOut$totalIter, model_ris$algOut$logpd$lDelta[1:model_ris$algOut$totalIter], type = 'n',
     xlab = 'Iteration', main = 'Lead Level Study', ylab = '',
     cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25)
title(ylab = expression(paste(Delta, log(p))), cex.lab = 1.5, line = 2)
abline(v = axTicks(1), h = axTicks(2), lty = 6, col = 'lightgray', lwd = 2)
lines(1:model_ris$algOut$totalIter, model_ris$algOut$logpd$lDelta[1:model_ris$algOut$totalIter], type = 'b',
      lwd = 2)
# dev.off()

summary(model_ri)


#### DTI Data (Goldsmith et al. 2012) ####
# The MRI/DTI data were collected at Johns Hopkins University and the Kennedy-Krieger Institute #
library(refund)

##### perform fpca first on functional predictors #####
cca   <- fpca.sc(DTI2$cca)
Fc    <- cca$Yhat

rcst  <- fpca.sc(DTI2$rcst)
Fr    <- rcst$Yhat

### make binary X ###
DTI2$X    <- 1*(DTI2$visit >= 2)

##### CCA #####
cca_fri   <- vme(pasat ~ X, Z = NULL, rint.only = TRUE, fe = ~ Fc, data = DTI2, ID = 'id',
                 ps.spline = list(df = 5, degree = 3, xi = 0.01),
                 re.spline = list(df = NULL, degree = 3, xi = 0.01),
                 priors = list(Ae = 0.01, Be = 0.01, Af = 0.01, Bf = 0.01,
                               Au = 0.01, Bu = 0.01, sigB = 1000),
                 tol = 1e-5, maxIter = 2000, up = 100,
                 verbose = TRUE)

ids     <- unique(DTI2$id)
n       <- length(ids)

BFcList   <- vector('list', length = n)

for(i in 1:n){
  idi           <- ids[i]
  ri            <- which(DTI2$id == idi)
  Fci           <- Fc[ri,]
  BFcList[[i]]  <- Fci
}

Bfc <- as.matrix(bdiag(BFcList))

cca_frf   <- vme(pasat ~ X, Z = Bfc, rint.only = FALSE, fe = ~ Fc, data = DTI2, ID = 'id',
                 ps.spline = list(df = 5, degree = 3, xi = 0.01),
                 re.spline = list(df = NULL, degree = 3, xi = 0.01),
                 priors = list(Ae = 0.01, Be = 0.01, Af = 0.01, Bf = 0.01,
                               Au = 0.01, Bu = 0.01, sigB = 1000),
                 tol = 1e-5, maxIter = 2000, up = 100,
                 verbose = TRUE)

VAIC(cca_fri) -
VAIC(cca_frf)

plot(cca_fri$funcEffect, type = 'l')

# pdf('cca.pdf')
plot(cca_fri$funcEffect, type = 'n', ylab = '', xlab = 't', main = 'CCA',
     cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25)
title(ylab = expression(gamma(t)), line = 2, cex.lab = 1.5)
abline(v = axTicks(1), h = axTicks(2), lty = 6, col = 'lightgray', lwd = 2)
lines(cca_fri$funcEffect, lty = 1, lwd = 2)
# dev.off()

summary(cca_fri)

# pdf('cca_logp.pdf')
plot(1:cca_fri$algOut$totalIter, cca_fri$algOut$logpd$lDelta[1:cca_fri$algOut$totalIter], type = 'n',
     xlab = 'Iteration', main = 'DTI Study: CCA Model', ylab = '',
     cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25)
title(ylab = expression(paste(Delta, log(p))), cex.lab = 1.5, line = 2)
abline(v = axTicks(1), h = axTicks(2), lty = 6, col = 'lightgray', lwd = 2)
lines(1:cca_fri$algOut$totalIter, cca_fri$algOut$logpd$lDelta[1:cca_fri$algOut$totalIter], type = 'b',
      lwd = 2)
# dev.off()


##### RCST #####
rcst_fri   <- vme(pasat ~ X, Z = NULL, rint.only = TRUE, fe = ~ Fr, data = DTI2, ID = 'id',
                  ps.spline = list(df = 5, degree = 3, xi = 0.01),
                  re.spline = list(df = NULL, degree = 3, xi = 0.01),
                  priors = list(Ae = 0.01, Be = 0.01, Af = 0.01, Bf = 0.01,
                                Au = 0.01, Bu = 0.01, sigB = 1000),
                  tol = 1e-5, maxIter = 2000, up = 100,
                  verbose = TRUE)

ids     <- unique(DTI2$id)
n       <- length(ids)

BFrList   <- vector('list', length = n)

for(i in 1:n){
  idi           <- ids[i]
  ri            <- which(DTI2$id == idi)
  Fri           <- Fc[ri,]
  BFrList[[i]]  <- Fri
}

Bfr <- as.matrix(bdiag(BFrList))

rcst_frf   <- vme(pasat ~ X, Z = Bfr, rint.only = FALSE, fe = ~ Fr, data = DTI2, ID = 'id',
                  ps.spline = list(df = 5, degree = 3, xi = 0.01),
                  re.spline = list(df = NULL, degree = 3, xi = 0.01),
                  priors = list(Ae = 0.01, Be = 0.01, Af = 0.01, Bf = 0.01,
                                Au = 0.01, Bu = 0.01, sigB = 1000),
                  tol = 1e-5, maxIter = 2000, up = 100,
                  verbose = TRUE)

VAIC(rcst_fri) - 
VAIC(rcst_frf)

# pdf('rcst.pdf')
plot(rcst_fri$funcEffect, type = 'n', ylab = '', xlab = 't', main = 'RCST',
     cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25)
title(ylab = expression(gamma(t)), line = 2, cex.lab = 1.5)
abline(v = axTicks(1), h = axTicks(2), lty = 6, col = 'lightgray', lwd = 2)
lines(rcst_fri$funcEffect, lty = 1, lwd = 2)
# dev.off()

summary(rcst_fri)

# pdf('rcst_logp.pdf')
plot(1:rcst_fri$algOut$totalIter, rcst_fri$algOut$logpd$lDelta[1:rcst_fri$algOut$totalIter], type = 'n',
     xlab = 'Iteration', main = 'DTI Study: RCST Model', ylab = '',
     cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25)
title(ylab = expression(paste(Delta, log(p))), cex.lab = 1.5, line = 2)
abline(v = axTicks(1), h = axTicks(2), lty = 6, col = 'lightgray', lwd = 2)
lines(1:rcst_fri$algOut$totalIter, rcst_fri$algOut$logpd$lDelta[1:rcst_fri$algOut$totalIter], type = 'b',
      lwd = 2)
# dev.off()