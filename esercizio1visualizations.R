# data visualization 1: downsampling using wrassp
library(wrassp)

f1 <- read.AsspDataObj("hw_data/hw_data/f1.au")

plot(
  seq(0,numRecs.AsspDataObj(f1) - 1, 100) / rate.AsspDataObj(f1),
  f1$audio[ c( TRUE,rep(FALSE, 99) ) ],
  type='l',
  xlab ='time (s)',
  ylab ='Audio samples'
  )

# data visualization 2: downsampling using tuneR
library(tuneR)

downsample_rate = 11025
f1wv <- Wave( as.numeric(f1$audio), samp.rate = rate.AsspDataObj(f1), bit = 16)
#play(f1wv) # do not run in Markdown
plot(f1wv)
f1w.dwn = downsample(f1wv, samp.rate = downsample_rate)
plot(f1w.dwn)

# data visualization 3: fundamental frequency contour
f0vals = ksvF0("hw_data/hw_data/f1.au", toFile=F)

plot(seq(0,numRecs.AsspDataObj(f0vals) - 1) / rate.AsspDataObj(f0vals) +
       attr(f0vals, 'startTime'),
     f0vals$F0, 
     type='l', 
     xlab='time (s)', 
     ylab='F0 frequency (Hz)')

# data visualization 4: vanilla Fourier Transform
cos.basis <-function(x = (1:500)/500, j.max = length(x)){
  n = length(x)
  L = x(length(x))
  mat.basis <- matrix(NA, nrow = n, ncol = j.max)
  mat.basis[,1] <- 1
  for(j in 2:j.max){
    mat.basis[,j] = sqrt(2 / L) * cos((j - 1) * pi / L * x)
  }
  
  return(mat.basis)
}

L <- 30.
t <- seq(0, 30, length.out = length(f1w.dwn))
max_basis_index <- 1000
basis.mat.t <- sqrt(2 / L) * cos(pi / L * outer((0:max_basis_index-1),t))

f1.freq <- 1 / downsample_rate * f1w.dwn@left %*% basis.mat.t
plot(1:max_basis_index,f1.freq)
