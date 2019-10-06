# Data were collected from a mine in Cobar, NSW, Australia. At each of 38 sampling points,
# several measurements were taken, one of which is the 'true-width' of an ore-bearing rock
# layer. Also given are the coordinates t1 and t2 of of the data sites.
load("hw_data/ore.RData")
#str(ore)
library(mgcv)
ore.gam <- gam(width ~ s(t1) + s(t2), data = ore)
plot(ore.gam, pages = 1)


# my attempt
library(KernSmooth)

cv.locpoly <- function(x,y, k=10) {
  CV.indices <- caret::createFolds(y, k=k)

  # cross-validating on 20 log-spaced bandwidths  
  range.x <- range(x)
  meas.x <- range.x[2] - range.x[1]

  # bw <- exp(seq(log(meas.x / length(x)), log(meas.x / 3), length.out=20))
  #bw <- seq(0.2, 0.95, length.out = 10)
  bw <- c(0.7,0.75)
  test.error <- rep(0, length(bw))
  
  for (i in 1:k) {
    # splitting data
    training.X <- x[-CV.indices[[i]]]
    training.Y <- y[-CV.indices[[i]]]
    
    validation.X <- x[CV.indices[[i]]]
    validation.Y <- y[CV.indices[[i]]]
    
    # estimating risk vs bandwidth
    for (b in 1:length(bw)) {  
      #y.hat <- locpoly(training.X, training.Y, degree=2, bandwidth=bw[b], 
      #                gridsize=200, range.x = range.x)
      #y.hat.validation.X <- approx(y.hat$x, y.hat$y, xout=validation.X, rule=2)$y

      poly.obj <- loess(training.Y ~ training.X)
      y.hat.validation.X <- predict(poly.obj, validation.X)
      test.error[b] <- (test.error[b] 
                        + 1/k * mean((y.hat.validation.X - validation.Y) ^ 2))
    }
  }
  
  # fitting using the optimal bandwidth
  bw.opt <- bw[which.min(test.error)]
  return(
    locpoly(x, y, degree=2, bandwidth=bw.opt, gridsize=200, range.x = range.x)
  )
}

backfit <- function (X,Y, incr_tol = 1e-02, maxiter = 1e+03) {
  d <- dim(X)[1]
  n <- dim(X)[2]
  
  row.range <- apply(X, 1, range)
  
  # initializing
  alpha.hat <- mean(Y)
  m <- list()
  M <- matrix(rep(0, d*n), nrow=d, ncol=n)
  r <- matrix(rep(0, d*n), nrow=d, ncol=n)
  
  # starting error
  rss0 <- sum((Y - alpha.hat - colSums(M)) ^ 2)
  
  # running
  #r.smoothed <- list()
  iter <- 0
  incr <- incr_tol + 1
  while ((incr > incr_tol) && (iter < maxiter)) {
    for (j in 1:d) {
      # updating residual
      r[j,] <- Y - alpha.hat - M[-j,] 
      # ^ should be "... - colSums(M[-j,])" BUT since d=2 M[-j,] has only got one row
      # smoothing residual
      # sm.grid <- cv.locpoly(X[j,], r[j,])
      #bb <- dpill(X[j,],r[j,])
      #sm.grid <- locpoly(X[j,], r[j,], degree=1, bandwidth=bb, range.x=row.range[,j], gridsize=200)
      
      # sm.fun <- approxfun(sm.grid$x, sm.grid$y, rule=2)
      
      sm.obj[[j]] <- cv.locpoly(X[,j], r[,j])
      M[,j] <- sm.obj[[j]]$fitted
      m[[j]] <- function(x) predict(sm.obj[[j]], x) - mean(M[,j])
      
      #r.smoothed[[j]] <- smooth.spline(X[j,], r[j,])
      #m[[j]] <- function(v) predict(r.smoothed[[j]], v)$y
      # M[j,] <- sm.fun(X[j,])
      # recentering
      M[j,] <- M[j,] - mean(M[j,])
     # m[[j]] <- approxfun(X[j,], M[j,], rule=2)
    }
    
    # increment control
    rss <- sum((Y - alpha.hat - colSums(M)) ^ 2)
    incr <- abs(rss - rss0) / rss
    rss0 <- rss
    iter <- iter + 1
  }
  
  # building additive model
  model <- function(x) {
    pred <- rep(alpha.hat, dim(x)[2])
    for (j in 1:(dim(x)[1])) {
      pred <- pred + m[[j]](x[j,])
    }
    return(pred)
  }
  
  return(list(m=m, model=model, iterations=iter, M=M))
}


# applying to ore data
X <- rbind(ore$t1, ore$t2)
pred <- backfit(X, ore$width)

# plotting results
x1 <- seq(range(ore$t1)[1], range(ore$t1)[2], length.out=100)
x2 <- seq(range(ore$t2)[1], range(ore$t2)[2], length.out=100)

par(mfrow=c(1,2))

plot(ore$t1, ore$width)
lines(x1, mean(ore$width) + pred$m[[1]](x1))

plot(ore$t2, ore$width)
lines(x2, mean(ore$width) + pred$m[[2]](x2))


# comparing to mgcv.gam predictions
gam.pred <- predict(ore.gam, type='terms')

# ...through mse...
mse.t1 <- mean((gam.pred[,1] - pred$m[[1]](ore$t1))^2)
mse.t2 <- mean((gam.pred[,2] - pred$m[[2]](ore$t2))^2)

print(c(mse.t1, mse.t2))

# ...through overlaid plots
par(mfrow=c(1,2))

plot(ore$t1, gam.pred[,1])
lines(ore$t1, pred$m[[1]](ore$t1))

plot(ore$t2, gam.pred[,2])
lines(sort(ore$t2), pred$m[[2]](sort(ore$t2)))


# splines plot
#par(mfrow=c(1,2))

#plot(ore$t1, ore$width)
#lines(ore$t1, mean(ore$width) + pred$M[1,])

#plot(ore$t2, ore$width)
#points(ore$t2, mean(ore$width) + pred$M[2,])