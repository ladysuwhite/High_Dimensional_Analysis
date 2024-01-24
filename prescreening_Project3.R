####
### this part is mainly focusing on pre-screening process. The output is briefly checked here, and is served as the
### input for Yichi's part.
####
library(glinternet)
library(dplyr)
library(tibble)
library(VariableScreening)
library(devtools)
install_github('wwrechard/screening')
library(screening)


setwd("")
load("SNP-gene data.RData")
dim1 <- dim(X)[1]  #row 
dim2 <- dim(X)[2]  #column

###### all level matrix ######
X.new <- matrix(, nrow = dim1, ncol = 0)
for (jj in 1:dim2){
  level <- length(unique(X[,jj]))
  seq <- c(1:level)-1
  for (kk in seq){
    tempp <- c(rep(0,dim1))
    tempp[which(X[,jj] == kk)] <- 1
    #X.new <- cbind(X.new, tempp)
    #X.new[paste0()] <- "Value"
    # X.new <- X.new %>%
    #   add_column(Empty_Column = tempp) %>%
    X.new <- cbind(X.new, tempp)
    }
  }
kkk <- as.vector(apply(X, 2, function(X){length(unique(X))}))
colnameshere <- rep(colnames(X), times = kkk)
colnames(X.new) <- colnameshere

#### SIRS
answer.screeniid.sisny <- screenIID(X.new, y, method = "MV-SIS-NY")
#print(which(answer.screeniid.sisny$rank<=100))
order.screeniid.sisny.200 <- order(answer.screeniid.sisny$rank)[1:200]
names.order.screeniid.sisny.100 <- unique(colnames(X.new)[order.screeniid.sisny.200])[1:100]
index.screeniid.sisny.100 <- c()
for (kkk in 1:100){
  temp <- which(colnames(X) == names.order.screeniid.sisny.100[kkk])
  index.screeniid.sisny.100 <- c(index.screeniid.sisny.100, temp)
}
save(index.screeniid.sisny.100, file = "First_100_ScreeniidSIRS.Rdata")

##2000
order.screeniid.sisny.3000 <- order(answer.screeniid.sisny$rank)[1:3000]
names.order.screeniid.sisny.2000 <- unique(colnames(X.new)[order.screeniid.sisny.3000])[1:2000]
index.screeniid.sisny.2000 <- c()
for (kkk in 1:2000){
  temp <- which(colnames(X) == names.order.screeniid.sisny.2000[kkk])
  index.screeniid.sisny.2000 <- c(index.screeniid.sisny.2000, temp)
}
save(index.screeniid.sisny.2000, file = "First_2000_ScreeniidSIRS.Rdata")


## same results
answer.screeniid.sirs <- screenIID(X.new, y, method = "SIRS")
print(which(answer.screeniid.sirs$rank<=10))

answer.screeniid.dcsis <- screenIID(X.new, y, method = "DC-SIS")
print(which(answer.screeniid.dcsis$rank<=10))

##SIS
##100
sis.screening.200 <- screening(X.new, y, method = 'sis', num.select = 200)$screen
names.order.sis.100 <- unique(colnames(X.new)[sis.screening.200])[1:100]
index.sis.100 <- c()
for (kkk in 1:100){
  temp <- which(colnames(X) == names.order.sis.100[kkk])
  index.sis.100 <- c(index.sis.100, temp)
}
save(index.sis.100, file = "First_100_SIS.Rdata")

##2000
sis.screening.3000 <- screening(X.new, y, method = 'sis', num.select = 3000)$screen
names.order.sis.2000 <- unique(colnames(X.new)[sis.screening.3000])[1:2000]
index.sis.2000 <- c()
for (kkk in 1:2000){
  temp <- which(colnames(X) == names.order.sis.2000[kkk])
  index.sis.2000 <- c(index.sis.2000, temp)
}
save(index.sis.2000, file = "First_2000_SIS.Rdata")



###RRCS
#output <- screening(X.new, y, method = 'holp', num.select = 10)
rrcs.screening.200 <- screening(X.new, y, method = 'rrcs', num.select = 200)$screen
names.order.rrcs.100 <- unique(colnames(X.new)[rrcs.screening.200])[1:100]
index.rrcs.100 <- c()
for (kkk in 1:100){
  temp <- which(colnames(X) == names.order.rrcs.100[kkk])
  index.rrcs.100 <- c(index.rrcs.100, temp)
}
save(index.rrcs.100, file = "First_100_RRCS.Rdata")

#### 2000
rrcs.screening.3000 <- screening(X.new, y, method = 'rrcs', num.select = 3000)$screen
names.order.rrcs.2000 <- unique(colnames(X.new)[rrcs.screening.3000])[1:2000]
index.rrcs.2000 <- c()
for (kkk in 1:2000){
  temp <- which(colnames(X) == names.order.rrcs.2000[kkk])
  index.rrcs.2000 <- c(index.rrcs.2000, temp)
}
save(index.rrcs.2000, file = "First_2000_RRCS.Rdata")


index.sis.100
index.rrcs.100
index.screeniid.sisny.100  #sirs



#################  briefly check here, not presented in slides  #################

#split dataset
set.seed(790)
train <- sample(dim1, dim1 * 0.8)  #80-20 train/test
X.train <- X[train, ]
X.test <- X[-train, ]
y.train <- y[train]
y.test <- y[-train]

num.levels <- c(rep(3, dim2))

limits.range <- c(1:100)
mse.test.vec <- c()
lambdaHat.cv.glin <- c()
lambdaHat1Std.cv.glin <- c()
cverr.cv.glin <- c()
cverrstd.cv.glin <- c()

for (i in limits.range){
  cv.glin <- glinternet.cv(X.train, y.train, numLevels = num.levels, 
                           #lambda = c(exp(seq(2,-2,,100))),
                           nLambda = 50, interactionCandidates=NULL, interactionPairs=NULL, screenLimit = i,
                           verbose = TRUE, tol = 1e-05, maxIter=5000)
  #plot(cv.glin)
  print(paste0("glinternet.cv$lambda.min for limit is ", i));cv.glin$lambdaHat
  print(paste0("glinternet.cv$lambda.1se for limit is ", i));cv.glin$lambdaHat1Std
  print(paste0("glinternet.cv$cvErr for limit is ", i));cv.glin$cvErr
  print(paste0("glinternet.cv$cvErrStd for limit is ", i));cv.glin$cvErrStd
  lambdaHat.cv.glin <- c(lambdaHat.cv.glin, cv.glin$lambdaHat)
  lambdaHat1Std.cv.glin <- c(lambdaHat1Std.cv.glin, cv.glin$lambdaHat1Std)
  cverr.cv.glin <- rbind(cverr.cv.glin, cv.glin$cvErr)
  cverrstd.cv.glin <- rbind(cverrstd.cv.glin, cv.glin$cvErrStd)
  save(cv.glin, file = paste0("screenlimit=", i, ".RData"))
}

## single
fit.glin <- glinternet(X.train, y.train, numLevels = num.levels, 
                       # lambda = c(exp(seq(2,-2,,100))),
                       nLambda = 100, interactionCandidates=NULL, interactionPairs=NULL, screenLimit = i,
                       verbose = TRUE, tol = 1e-05, maxIter=50000)
#plot(fit.glin, label = TRUE, xvar = "lambda", main = "glinternet plot with different lambdas")

glin.coef <- coef(fit.glin, lambdaIndex = cv.glin$lambdaHat1Std)  #coefficients

mse.test <- mean((predict(fit.glin, X.test, lambda = cv.glin$lambdaHat1Std) - y.test)^2)
print("The test mse is: "); mse.test
mse.test.vec <- c(mse.test.vec, mse.test)

#the fitting output
glin.coef$mainEffects #A list with components cat and cont
glin.coef$mainEffectsCoef #List of coefficients for the main effects in mainEffects
glin.coef$interactions #List of interactions, with components contcont, catcont and catcat, each 2-column matrices of variable indices
glin.coef$interactionsCoef #List of interaction coefficients for interactions



