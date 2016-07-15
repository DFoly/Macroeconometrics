###############################
## VAR with sign restrictions
###############################
source('R/functions.R')
library(Matrix); library(fBasics); library(MASS)
reps = 1000
burn = 700

usdata1 <- read.csv("C:/Users/dfoley/Dropbox/R/Macroeconometrics/Data/usdata1.csv",
                    header = FALSE)
names(usdata1) <- c('FFRate', 'GDP', 'CPI', 'Cons', 'URate',
                   'Inv','Exp','M2', 'Bond', 'Stocks', 'Exrate')

N = ncol(usdata1)
L = 2 # laglength of VAR
Y = usdata1

X <- Y
# loop to create lags
#for (i in 1:2){
  #temp <- lag0(X,i)
  #X <- cbind(X,temp)
#}

X <- cbind(lag0(Y,1), lag0(Y,2))


# delete first L+1 pbservations and add constant
ones <- rep(1,nrow(Y))
X <- cbind(ones, X)
X <- X[(L+1):nrow(X),]
Y <- Y[(L+1):nrow(Y),]

T1 = nrow(X)

#Priors for VAR coefficients
lambdaP <- 1  # tightness of prior on first lag
tauP <- 10*lambdaP
epsilonP <- 1 # tightness of prior on constant
muP <-t(colMeans(Y))
sigmaP <- matrix()
deltaP <- matrix()

# initial OLS estimates
for (i in 1:N){
  ytemp <- matrix(Y[,i])
  xtemp <- cbind(lag0(ytemp,1),matrix(1,nrow = nrow(Y), ncol = 1))
  ytemp <- ytemp[2:nrow(ytemp),]
  xtemp <- xtemp[2:nrow(xtemp),]
  btemp <- solve(t(xtemp) %*% xtemp) %*% (t(xtemp) %*% ytemp)
  etemp <- ytemp - xtemp %*% btemp
  stemp <- (t(etemp) %*% etemp)/length(ytemp)
  deltaP[i] <-  btemp[1]
  sigmaP[i] <-  stemp
}


# build dummy observations that incorporates priors
# return x and y
# result has two outputs so need to return as list and convert back to matrix
result <- create_dummies(lambdaP,tauP,deltaP,epsilonP,L,muP,sigmaP,N)
x <- result[1][1]
# convert back to matrix
x <- matrix(unlist(x), ncol = 23)
y <- result[2][1]
y <- matrix(unlist(y), ncol = N)
# append data to dummies
Y = as.matrix(Y)
X = as.matrix(X)
Y0 <- rbind(Y,y)
X0 <- rbind(X,x)

#conditional mean of the coefficients
mstar <- solve(t(X0) %*% X0) %*% t(X0) %*% Y0
mstar <- vec(mstar)
xx <- t(X0) %*% X0
ixx <- solve(xx)
sigma <- diag(1,N)
out <- array(0, c(reps-burn,36,N))

jj=1

# main gibbs sampling alogrithm
for (i in 1:reps){
  # posterior mean and variance
  vstar <- kronecker(sigma,ixx)
  # had to create function to calculate chol of non positive definite matrix
  # returns lower triangular so need to trasnspose for upper
  beta <- mstar + t(t(matrix(rnorm(N*(N*L+1),1))) %*% t(my_chol(vstar,10)))

  # draw covariance from inverse wishart
  e = Y0 - X0 %*% matrix(beta, nrow = (N*L+1), ncol = N)
  scale = t(e) %*% e
  sigma = IWPQ(T1+nrow(y), ginv(scale))

  if(i > burn){
    # impose sign restrictions
    chck = -1
    while(chck < 0){

      K = matrix(rnorm(N*N),N)
      Q = getqr(K)
      A0hat = chol(sigma)
      A0hat1 = (Q %*% A0hat) # candidate draw

      # check signs
      # first row corresponds to t = 1 of shock
      e1 = A0hat1[1,1] > 0  # Response of R
      e2 = A0hat1[1,2] < 0  # Response of Y
      e3 = A0hat1[1,3] < 0  # Response of Inflation
      e4 = A0hat1[1,4] < 0  # Response of CConsumption
      e5 = A0hat1[1,5] > 0  # Response of U
      e6 = A0hat1[1,6] < 0  # Response of Investment
      e7 = A0hat1[1,7] < 0  # Response of Money
      if(e1+e2+e3+e4+e5+e6+e7 == 7){
        chck =10
      }else{
        # we swtich signs of shock
        e1 = A0hat1[1,1] > 0  # Response of R
        e2 = A0hat1[1,2] < 0  # Response of Y
        e3 = A0hat1[1,3] < 0  # Response of Inflation
        e4 = A0hat1[1,4] < 0  # Response of CConsumption
        e5 = A0hat1[1,5] > 0  # Response of U
        e6 = A0hat1[1,6] < 0  # Response of Investment
        e7 = A0hat1[1,7] < 0  # Response of Money
        if(e1+e2+e3+e4+e5+e6+e7 == 7){
          A0hat1[1,1:N] = -A0Hat1[1,1:N]
          chck = 10
        }
      }

    }# end of while loop
    yhat = matrix(0,36,N) # s6 periods for IRF
    vhat = matrix(0,36,N)
    vhat[3,1] = 1  # shock to FFrate

    #impulse response functions
    for(j in 3:36){
      yhat[j,] <-    t(matrix(c(yhat[j-1,],
                       yhat[j-2,], 0))) %*% matrix(beta, nrow = (N*L+1), ncol = N) +
                       (vhat[j,] %*% A0hat1)
    }
    out[JJ, , ] <- yhat
    jj = jj+1


  }

} # end of gibbs sampling
