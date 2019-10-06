# Data were collected from a mine in Cobar, NSW, Australia. At each of 38 sampling points,
# several measurements were taken, one of which is the 'true-width' of an ore-bearing rock
# layer. Also given are the coordinates t1 and t2 of of the data sites.
load("hw_data/ore.RData")
str(ore)
library(mgcv)
ore.gam <- gam(width ~ s(t1) + s(t2), data = ore)
plot(ore.gam, pages = 1)


# my attempt
source("ksmoother.R")

k.smoother.loo <- function(x,y, bandwidths = exp(seq(log(0.05), log(5), length.out=30))) {
  risks <- rep(0, length(bandwidths))
  for (b in 1:length(bandwidths)) {
    risks[b] <- k.smoother(x,y, bandwidths[b])$risk
  }
  
  return(
    k.smoother(x,y, bandwidths[which.min(risks)])
  )
}

backfit <- function (X,Y, incr_tol = 1e-04, maxiter = 1e+03) {
  d <- dim(X)[2]
  n <- dim(X)[1]
  
  # initializing
  alpha.hat <- mean(Y)
  M <- matrix(rep(0, d*n), nrow=n, ncol=d)
  r <- matrix(rep(0, d*n), nrow=n, ncol=d)
  
  # starting error
  rss0 <- sum((Y - alpha.hat - rowSums(M)) ^ 2)
  
  # running
  sm.obj <- list()
  iter <- 0
  incr <- incr_tol + 1
  while ((incr > incr_tol) && (iter < maxiter)) {
    for (j in 1:d) {
      # updating residual
      r[,j] <- Y - alpha.hat - M[,-j]
      # smoothing residual
      sm.obj[[j]] <- k.smoother.loo(X[,j], r[,j])
      # recentering
      M[,j] <- sm.obj[[j]]$y - mean(sm.obj[[j]]$y)
    }
    
    # increment control
    rss <- sum((Y - alpha.hat - rowSums(M)) ^ 2)
    incr <- abs(rss - rss0) / rss
    rss0 <- rss
    iter <- iter + 1
  }
  
  # building additive model
  model <- function(x) {
    d <- dim(x)[2]
    n <- dim(x)[1]
    
    pred <- rep(alpha.hat, n)
    for (j in 1:d) {
      pred <- pred + sm.obj[[j]]$func(x[,j]) - mean(sm.obj[[j]]$y)
    }
    return(pred)
  }
  
  return(list(model=model, iterations=iter, M=M))
}


# applying to ore data
X <- cbind(ore$t1, ore$t2)

# fitting
pred <- backfit(X, ore$width)

# plotting component functions over data
par(mfrow=c(1,2))
t2.sorted <- sort(ore$t2, index.return=T)

plot(ore$t1, ore$width)
lines(ore$t1, mean(ore$width) + pred$M[,1])

plot(ore$t2, ore$width)
lines(t2.sorted$x, mean(ore$width) + pred$M[t2.sorted$ix,2])

# comparing to mgcv.gam predictions
gam.pred <- predict(ore.gam)

# ...through mse...
mse <- mean((gam.pred - pred$model(X))^2)
print(mse)

# ...through overlaid plots
gam.terms <- predict(ore.gam, type='terms')
par(mfrow=c(1,2))

plot(ore$t1, gam.terms[,1], type='l')
lines(ore$t1, pred$M[,1], col="red")

plot(t2.sorted$x, gam.terms[t2.sorted$ix,2], type='l')
lines(t2.sorted$x, pred$M[t2.sorted$ix,2], col="red")