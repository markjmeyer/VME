#### libraries ####
library(xtable)

#### random intercept+slope simulations ####

##### n = 20 #####
# setwd('~InterSlope/N20')

##### n = 50 #####
# setwd('~InterSlope/N50')

##### n = 100 #####
# setwd('~InterSlope/N100')

vaicMat_i_is   <- readRDS('vaicMat_i.rds')
vaicMat_b_is   <- readRDS('vaicMat_b.rds')

biasMat_i_is   <- readRDS('biasMat_i.rds')
biasMat_b_is   <- readRDS('biasMat_b.rds')

miseMat_i_is   <- readRDS('miseMat_i.rds')
miseMat_b_is   <- readRDS('miseMat_b.rds')

Correct_i_is   <- array(NA, dim = c(dim(vaicMat_i_is)[1], dim(vaicMat_i_is)[2]))
VAICDif_i_is   <- array(NA, dim = c(dim(vaicMat_i_is)[1], dim(vaicMat_i_is)[2]))
Within2_i_is   <- matrix(0, nrow = dim(vaicMat_i_is)[1], ncol = dim(vaicMat_i_is)[2])
Within5_i_is   <- matrix(0, nrow = dim(vaicMat_i_is)[1], ncol = dim(vaicMat_i_is)[2])
Within10_i_is  <- matrix(0, nrow = dim(vaicMat_i_is)[1], ncol = dim(vaicMat_i_is)[2])

Correct_b_is   <- array(NA, dim = c(dim(vaicMat_i_is)[1], dim(vaicMat_i_is)[2]))
VAICDif_b_is   <- array(NA, dim = c(dim(vaicMat_i_is)[1], dim(vaicMat_i_is)[2]))
Within2_b_is   <- matrix(0, nrow = dim(vaicMat_i_is)[1], ncol = dim(vaicMat_i_is)[2])
Within5_b_is   <- matrix(0, nrow = dim(vaicMat_i_is)[1], ncol = dim(vaicMat_i_is)[2])
Within10_b_is  <- matrix(0, nrow = dim(vaicMat_i_is)[1], ncol = dim(vaicMat_i_is)[2])

for(i in 1:dim(vaicMat_i_is)[1]){
  for(j in 1:dim(vaicMat_i_is)[2]){
    Correct_i_is[i, j]   <- 1*(vaicMat_i_is[i, j, 2] <= vaicMat_i_is[i, j,1])
    VAICDif_i_is[i, j]   <- vaicMat_i_is[i, j, 2] - vaicMat_i_is[i, j, 1]
    if(VAICDif_i_is[i, j] > 0){
      Within2_i_is[i, j]   <- 1*(VAICDif_i_is[i, j] <= 2)
      Within5_i_is[i, j]   <- 1*(VAICDif_i_is[i, j] <= 5)
      Within10_i_is[i, j]  <- 1*(VAICDif_i_is[i, j] <= 10)
    }
    
    Correct_b_is[i, j]  <- 1*(vaicMat_b_is[i, j, 1] <= vaicMat_b_is[i, j, 2])
    VAICDif_b_is[i, j]  <- vaicMat_b_is[i, j, 1] - vaicMat_b_is[i, j, 2]
    if(VAICDif_b_is[i, j] > 0){
      Within2_b_is[i,j]    <- 1*(VAICDif_b_is[i, j] <= 2)
      Within5_b_is[i,j]    <- 1*(VAICDif_b_is[i, j] <= 5)
      Within10_b_is[i,j]   <- 1*(VAICDif_b_is[i, j] <= 10)     
    }
    
  }
}



vaici_is <- matrix(c(mean(apply(Correct_i_is, MARGIN = 1, FUN = mean)), 
                         mean(apply(Correct_i_is + Within2_i_is, MARGIN = 1, FUN = mean)), 
                         mean(apply(Correct_i_is + Within5_i_is, MARGIN = 1, FUN = mean)),
                         mean(apply(Correct_i_is + Within10_i_is, MARGIN = 1, FUN = mean))),
                       nrow = 1, byrow = TRUE)

vaicb_is <- matrix(c(mean(apply(Correct_b_is, MARGIN = 1, FUN = mean)), 
                         mean(apply(Correct_b_is - Within2_b_is, MARGIN = 1, FUN = mean)),
                         mean(apply(Correct_b_is - Within5_b_is, MARGIN = 1, FUN = mean)), 
                         mean(apply(Correct_b_is - Within10_b_is, MARGIN = 1, FUN = mean))),
                       nrow = 1, byrow = TRUE)

xtable(100*vaici_is, digits = 2)

xtable(100*vaicb_is, digits = 2)

paste(round(mean(biasMat_i_is[,,,2]), digits = 4), ' & ', round(mean(biasMat_i_is[,,,1]), digits = 6)) # true model first
paste(round(mean(miseMat_i_is[,,,2]), digits = 5), ' & ', round(mean(miseMat_i_is[,,,1]), digits = 4)) # true model first


paste(round(mean(biasMat_b_is[,,,2]), digits = 4), ' & ', round(mean(biasMat_b_is[,,,1]), digits = 4)) # true model second
paste(round(mean(miseMat_b_is[,,,2]), digits = 4), ' & ', round(mean(miseMat_b_is[,,,1]), digits = 4)) # true model second


mean(biasMat_i_is[1,,,2])
mean(biasMat_i_is[2,,,2])
mean(biasMat_i_is[3,,,2])
mean(biasMat_i_is[4,,,2])

mean(biasMat_i_is[1,,,1])
mean(biasMat_i_is[2,,,1])
mean(biasMat_i_is[3,,,1])
mean(biasMat_i_is[4,,,1])
