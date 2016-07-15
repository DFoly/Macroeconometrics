######################
# functions.R
########################
#p = lags, k = no columns
lag0 <- function(x,p){
  R = nrow(x)
  k = ncol(x)
  #k = ncol(x)
  # take first R-p rows
  x1 = as.matrix(x[1:(R-p),])
  # preceed them with p zeros
  out = matrix(0, nrow = p, ncol = k)
  out = as.matrix(rbind(out,x1))
  out
}


cols <- function(x){
  out <- dim(x)[2]
}

rows <- function(x){
  out = dim(x)[1]
}





##############################################
# function to check stability of matrix
##############################################
#   coef   (n*l+1) * n matrix with the coef from the VAR
#   l      number of lags
#   n      number of endog variables
#   FF     matrix with all coef
#   S      dummy var: if equal one->stability
stability <- function(beta,n,l){

  FF <- matrix(0, nrow = n*l, ncol = n*l)
  # companion form, skips first row
  # go across to column of number of variables as
  # only making identities for 1 lag so col is p-1 lags (l-1)
  FF[(n+1):(n*l), 1:n*(l-1)] <- diag(1, nrow = n*(l-1), ncol = n*(l-1))

  dim(beta) <- c(n*l+1,n)
  temp <- t(beta[2:(n*l)+1, 1:n])
  # state space companion form
  FF[1:n,1:(n*l)] <- temp
  ee = max(sapply(eigen(FF)$values,abs))
  S = ee > 1
  return(S)
}


#############################
# Cretea Dummies for priors
#############################
# as in Banbura et al 2007
create_dummies <- function(lambda, tau, delta, epsilon, p, mu, sigma, n){
  library(Matrix)
  # initialise empty matrix
  x <- matrix()
  y <- matrix()
  yd1 <- matrix()
  yd2 <- matrix()
  xd1 <- matrix()
  xd2 <- matrix()

  # eqn 5 banbuar et al 2007
  if(lambda > 0){
    if(epsilon > 0){
      temp1 <- as.numeric((sigma * delta)/lambda)
      temp1 <- diag( x= temp1)
      temp2 <- matrix(0, n*(p-1),n)
      temp3 <- diag(sigma)
      temp4 <- matrix(0,1,n)
      yd1 <- rbind(temp1,temp2,temp3,temp4)
      jp <- diag(1:p)

      temp5 <- cbind(kronecker(jp,diag(sigma)/lambda), matrix(0,(n*p),1))
      temp6 <- matrix(0,n,(n*p)+1)
      temp7 <- matrix(0,1, (n*p)+1)
      # epsilon in last row and column
      temp7[n*p+1] <- epsilon
      xd1 <- rbind(temp5,temp6,temp7)

      }else{

        temp1 <- as.numeric((sigma * delta)/lambda)
        temp1 <- diag( x= temp1)
        temp2 <- matrix(0, n*(p-1),n)
        temp3 <- diag(sigma)
        yd1 <- rbind(temp1,temp2,temp3)

        jp = diag(1:p)
        temp5 <- kronecker(jp, diag(sigma)/lambda)
        temp6 <- matrix(0, n, (n*p))
        xd1 <- rbind(temp5,temp6)


    }
  }


  if(tau > 0){
    if(epsilon > 0){
      temp8 <- as.numeric((delta*mu)/tau)
      yd2 <- diag(temp8)
      temp9 <- kronecker((1:p), yd2)
      xd2 <- cbind(t(temp9), matrix(0,n,1))

    }else{
      temp10 <- as.numeric((delta*mu)/tau)
      yd2 <- diag(temp10)
      xd2 <- t(kronecker((1:p), yd2))
    }

  }

  y = rbind(yd1,yd2)
  x = rbind(xd1,xd2)
  return(list(x = x, y = y))


}


########################################
### Draw from the inverse Wishart Dist
########################################
IWPQ <- function(v,ixpx){
  library(MASS)
  k = nrow(ixpx)
  z = matrix(0, nrow = v, ncol = k)
  mu = matrix(0, nrow = k, 1)
  for( i in 1:v){
    z[i,] <- t((t(my_chol(ixpx,k)) %*% rnorm(k,1)))
  }
  out <- ginv(t(z) %*% z)
  out
}


###################################
## QR decompostion of a matrix
###################################
getqr <- function(a){
  # Returns a modified QR decomposition of matrix a, such that the
  # diagonal elements of the 'R' matrix are all positive

  # [Q,R] = QR(A), where A is m-by-n, produces an m-by-n upper triangular
  #  matrix R and an m-by-m unitary matrix Q (i.e. Q Q'=I) so that A = Q*R.
  temp <- qr(a)
  Q <- qr.Q(temp, complete = TRUE) # keep dimensions of a
  R <- qr.R(temp, complete = TRUE)

  # If diagonal elements of R are negative then multiply the corresponding
  # column of Q by -1; Note: the modified Q matrix is still unitary.
  for (i in 1:nrow(Q)){
     if( R[i,i] < 0){
      Q[,i] = - Q[,i]
    }

  }
  out = Q
  return(out)
}

####################################################################
## performs cholesky decompsition on non positive definite matrix
####################################################################
# performs cholesky decompsition on semi positive definite matrix
my_chol_psd = function(a){
  n = dim(a)[1];
  root = matrix(0,n,n);

  for (i in 1:n){
    sum = 0;
    if (i>1){
      sum = sum(root[i,1:(i-1)]^2);
    }

    x = a[i,i] - sum;

    if (x<0){
      x = 0;
    }

    root[i,i] = sqrt(x);

    if (i < n){
      for (j in (i+1):n){

        if (root[i,i] == 0){
          x=0;
        }
        else{
          sum = 0;
          if (i>1) {
            sum = root[i,1:(i-1)] %*% t(t(root[j,1:(i-1)]))
          }
          x = (a[i,j] - sum)/root[i,i];
        }

        root[j,i] = x;
      }
    }
  }
  return(root);
}

##########################################
# CALCULATES PSUEDO INVERSE OF MATRIX
##########################################
my_chol = function(a,b){ # where a is the matrix and b is the block size
  n = dim(a)[1];
  out = a;
  x = floor(n/b-1)*b;
  i=1
  while(i<x){
    out[i:(i+b-1),i:(i+b-1)] = my_chol_psd(out[i:(i+b-1),i:(i+b-1)]);
    out[(i+b):n,i:(i+b-1)] = out[(i+b):n,i:(i+b-1)] %*% ginv(t(out[i:(i+b-1),i:(i+b-1)]));
    out[(i+b):n,(i+b):n] = out[(i+b):n,(i+b):n] - out[(i+b):n,i:(i+b-1)]%*%t(out[(i+b):n,i:(i+b-1)]);

    i = i + b;
  }
  out[(x+1):n,(x+1):n] =my_chol_psd(out[(x+1):n,(x+1):n]);
  for (i in 1:(n-1)){
    out[i,(i+1):n] = rep(0,n-i);
  }
  return(out);
}
