# load in lag function
source('R/functions/lag0.R')
Y.df <- read.csv("J:/DFoley/Matlab code/code/chapter1/data/inflation.csv")
names(Y.df) <- c('Date', 'Inflation')
Y <- data.frame(Y.df[,2])

T1 = nrow(Y)

ones = rep(1,T1)
X = matrix(c(ones,lag0(Y,1), lag0(Y,2)),
           nrow = T1, ncol = 3)
Y =  as.matrix(Y, nrow = T1)

# remove missing values at beginning
Y = as.matrix(Y[3:T1,])
X = as.matrix(X[3:T1,])
T1 = nrow(X)
k = ncol(X)
# set up prior mean and variance
# 1 variable 2 lags and constant so k = 3
B0 = c(0,0,0)
B0 <- as.matrix(B0, nrow = 1, ncol = k)
sigma0 <- diag(1,k)

# priors for sigma
T0 = 1    # degrees of freedom
D0 = 0.1  # scale

# starting values
B = B0
sigma2 = 1
reps = 5000
#burn =  4000

out = matrix(0, nrow = reps, ncol = k+1)
colnames(out) <- c('constant', 'beta1','beta2', 'sigma')
#matrix to store forecasts
out1 <- matrix(0, nrow = reps, ncol = 14)

# gibbs loop
# define 1/sigma2 as.numeric to solve
# NB dont use burn until end
# define out matrix that extract all coefficients at each period and
# append by row inot out matrix

for(i in 1:reps){
  M = solve(solve(sigma0) + as.numeric(1/sigma2) * t(X) %*% X) %*%
      (solve(sigma0) %*% B0 + as.numeric(1/sigma2) * t(X) %*% Y)

  V = solve(solve(sigma0) + as.numeric(1/sigma2) * t(X) %*% X)

  chck = -1
  while(chck < 0){   # check for stability

    B <- M + t(rnorm(3) %*% chol(V))
    b = matrix(c(B[2], 1, B[3], 0), nrow = k-1, ncol = k-1)
    ee <- max(sapply(eigen(b)$values,abs))
    if( ee<=1){
      chck=1
    }
  }
  # compute residuals
  resids <- Y- X%*%B
  T2 = T0 + T1
  D1 = D0 + t(resids) %*% resids

  #draw from Inverse Gamma
  z0 = rnorm(T1,1)
  z0z0 = t(z0) %*% z0
  sigma2 = D1/z0z0

  # keeps samples after burn period
  out[i,] <- t(matrix(c(t(B),sigma2)))

  # compute 2 year forecasts
  yhat = rep(0,14)
  end = as.numeric(length(Y))
  yhat[1:2] = Y[(end-1):end,]
  cfactor = sqrt(sigma2)
  for(m in 3:14){
    yhat[m] = c(1,yhat[m-1],yhat[m-2]) %*% B + rnorm(1) * cfactor
  }
  out1[i,] <- yhat
}
# burn first 4000
out <- out[4001:5000,]
# stores forecasts
out2 <- out1[4001:5000,]

# calculate posterior mean which is just col mean of each variable
# similar to matlab answer
post_means <- colMeans(out)
forecasts <- as.matrix(colMeans(out2))

##############################################################
library(matrixStats); library(ggplot2); library(reshape2)

# quantiles for all data points, makes plotting easier


error_bands <- colQuantiles(out1,prob = c(0.16,0.84))
Y_temp = cbind(Y,Y)
error_bands <- rbind(Y_temp, error_bands[3:dim(error_bands)[1],])
all <- c(Y[1:(length(Y)-2)],forecasts)

forecasts.mat <- cbind.data.frame(error_bands[,1],all, error_bands[,2])
names(forecasts.mat) <- c('lower', 'mean', 'upper')


# create date vector for plotting
Date <- seq(as.Date('1948/03/01'), as.Date('2013/06/01'), 'quarters')
# regular cbind loses data class
data.plot <- cbind.data.frame(Date, forecasts.mat)

# puts into long format for plotting
# easiest way to plot I found was to generate
# intervals for all data and forecasts
# intervals for actual data WOULD JUST BE REPEATED as actual data

data.plot_melt <- melt(data.plot, id = 'Date')
ggplot(data.plot_melt, aes(x = Date, y = value, colour = variable)) +
  geom_line()

ggplot(data.plot_melt, aes(x = Date, y = value)) +
  geom_line(data = data.plot_melt, aes(colour = variable))

# plot with confidence intervals
ggplot(data.plot, aes(x = Date, y = mean))+ geom_line(colour = 'blue') +
    geom_ribbon(data = data.plot ,
    aes(ymin = lower, ymax = upper, alpha = 0.2))

