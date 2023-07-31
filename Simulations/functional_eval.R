#### libraries ####
library(xtable)

#### random function simulations ####

##### settings #####
### set working directory to correct path ###
# type  <- 'Sigmoidal'
# N     <- 100
# Tt    <- 100
# fname <- paste('~Functional/', type, '/N', N, '/T', Tt, sep = '' )
# 
# setwd(fname)


vaicMat_i_is   <- readRDS('vaicMat_i.rds')
vaicMat_f_is   <- readRDS('vaicMat_f.rds')

biasMat_i_is   <- readRDS('biasMat_i.rds')
biasMat_f_is   <- readRDS('biasMat_f.rds')

miseMat_i_is   <- readRDS('miseMat_i.rds')
miseMat_f_is   <- readRDS('miseMat_f.rds')

vaicMat_i_is   <- vaicMat_i
vaicMat_f_is   <- vaicMat_f

biasMat_i_is   <- biasMat_i
biasMat_f_is   <- biasMat_f

miseMat_i_is   <- miseMat_i
miseMat_f_is   <- miseMat_f

Correct_i_is   <- array(NA, dim = c(dim(vaicMat_i_is)[1], dim(vaicMat_i_is)[2]))
VAICDif_i_is   <- array(NA, dim = c(dim(vaicMat_i_is)[1], dim(vaicMat_i_is)[2]))
Within2_i_is   <- matrix(0, nrow = dim(vaicMat_i_is)[1], ncol = dim(vaicMat_i_is)[2])
Within5_i_is   <- matrix(0, nrow = dim(vaicMat_i_is)[1], ncol = dim(vaicMat_i_is)[2])
Within10_i_is  <- matrix(0, nrow = dim(vaicMat_i_is)[1], ncol = dim(vaicMat_i_is)[2])

Correct_f_is   <- array(NA, dim = c(dim(vaicMat_i_is)[1], dim(vaicMat_i_is)[2]))
VAICDif_f_is   <- array(NA, dim = c(dim(vaicMat_i_is)[1], dim(vaicMat_i_is)[2]))
Within2_f_is   <- matrix(0, nrow = dim(vaicMat_i_is)[1], ncol = dim(vaicMat_i_is)[2])
Within5_f_is   <- matrix(0, nrow = dim(vaicMat_i_is)[1], ncol = dim(vaicMat_i_is)[2])
Within10_f_is  <- matrix(0, nrow = dim(vaicMat_i_is)[1], ncol = dim(vaicMat_i_is)[2])

for(i in 1:dim(vaicMat_i_is)[1]){
  for(j in 1:dim(vaicMat_i_is)[2]){
    Correct_i_is[i, j]   <- 1*(vaicMat_i_is[i, j, 2] <= vaicMat_i_is[i, j,1])
    VAICDif_i_is[i, j]   <- vaicMat_i_is[i, j, 2] - vaicMat_i_is[i, j, 1]
    if(VAICDif_i_is[i, j] > 0){
      Within2_i_is[i, j]   <- 1*(VAICDif_i_is[i, j] <= 2)
      Within5_i_is[i, j]   <- 1*(VAICDif_i_is[i, j] <= 5)
      Within10_i_is[i, j]  <- 1*(VAICDif_i_is[i, j] <= 10)
    }
    
    Correct_f_is[i, j]  <- 1*(vaicMat_f_is[i, j, 1] <= vaicMat_f_is[i, j, 2])
    VAICDif_f_is[i, j]  <- vaicMat_f_is[i, j, 1] - vaicMat_f_is[i, j, 2]
    if(VAICDif_f_is[i, j] > 0){
      Within2_f_is[i,j]    <- 1*(VAICDif_f_is[i, j] <= 2)
      Within5_f_is[i,j]    <- 1*(VAICDif_f_is[i, j] <= 5)
      Within10_f_is[i,j]   <- 1*(VAICDif_f_is[i, j] <= 10)     
    }
    
  }
}



vaici_is <- matrix(c(mean(apply(Correct_i_is, MARGIN = 1, FUN = mean)), 
                     mean(apply(Correct_i_is + Within2_i_is, MARGIN = 1, FUN = mean)), 
                     mean(apply(Correct_i_is + Within5_i_is, MARGIN = 1, FUN = mean)),
                     mean(apply(Correct_i_is + Within10_i_is, MARGIN = 1, FUN = mean))),
                   nrow = 1, byrow = TRUE)

vaicb_is <- matrix(c(mean(apply(Correct_f_is, MARGIN = 1, FUN = mean)), 
                     mean(apply(Correct_f_is - Within2_f_is, MARGIN = 1, FUN = mean)),
                     mean(apply(Correct_f_is - Within5_f_is, MARGIN = 1, FUN = mean)), 
                     mean(apply(Correct_f_is - Within10_f_is, MARGIN = 1, FUN = mean))),
                   nrow = 1, byrow = TRUE)

xtable(100*vaici_is, digits = 2)

xtable(100*vaicb_is, digits = 2)

paste(round(mean(biasMat_i_is[,,2]), digits = 4), ' & ', round(mean(biasMat_i_is[,,1]), digits = 4)) # true model first
paste(round(mean(miseMat_i_is[,,2]), digits = 4), ' & ', round(mean(miseMat_i_is[,,1]), digits = 4)) # true model first


paste(round(mean(biasMat_f_is[,,2]), digits = 4), ' & ', round(mean(biasMat_f_is[,,1]), digits = 4)) # true model second
paste(round(mean(miseMat_f_is[,,2]), digits = 4), ' & ', round(mean(miseMat_f_is[,,1]), digits = 4)) # true model second


mean(biasMat_f_is[1,,2])
mean(biasMat_f_is[2,,2])
mean(biasMat_f_is[3,,2])
mean(biasMat_f_is[4,,2])

mean(biasMat_f_is[1,,1])
mean(biasMat_f_is[2,,1])
mean(biasMat_f_is[3,,1])
mean(biasMat_f_is[4,,1])

mean(miseMat_i_is[1,,2])
mean(miseMat_i_is[2,,2])
mean(miseMat_i_is[3,,2])
mean(miseMat_i_is[4,,2])

mean(miseMat_i_is[1,,1])
mean(miseMat_i_is[2,,1])
mean(miseMat_i_is[3,,1])
mean(miseMat_i_is[4,,1])
