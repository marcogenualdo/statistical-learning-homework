suppressMessages(require(wrassp,  quietly = T))
suppressMessages(require(signal, quietly = T))
suppressMessages(require(tuneR, quietly = T)) 

library(class)
library(randomForest)
library(kernlab)

suppressMessages(require(dimRed, quietly = T))


# DATA PREPROCESSING

n <- 150

x <- read.AsspDataObj("hw_data/hw_data/f1.au")
fs <- rate.AsspDataObj(x)

# STFT parameters
winsize  <- 2 * fs
hopsize  <- fs
noverlap <- winsize - hopsize

sp <- specgram(x = x$audio, Fs = fs, window = winsize, overlap = noverlap)
nwindows <- dim(sp$S)[2] - 3

# Frequency bands selection
nbands   <- 2^3
lowB <- 100
eps  <- .Machine$double.eps

# Number of seconds of the analyzed window
corrtime <- 15

energy.avg <- matrix(rep(0, n * nbands), nrow = n, ncol = nbands)
energy.sd <- matrix(rep(0, n * nbands), nrow = n, ncol = nbands)

for (j in 1:n) {
  x <- read.AsspDataObj(paste("hw_data/hw_data/f", j, ".au", sep=""))
  
  # Sampling rate
  fs <- rate.AsspDataObj(x)
  # Short-time fourier transform
  sp       <- specgram(x = x$audio, Fs = fs, window = winsize, overlap = noverlap)
  ntm      <- ncol(sp$S)  # number of (overlapping) time segments
  
  # Energy of bands
  fco    <- round( c(0, lowB*(fs/2/lowB)^((0:(nbands-1))/(nbands-1)))/fs*winsize)
  energy <- matrix(0, nbands, ntm)
  for (tm in 1:ntm){
    for (i in 1:nbands){
      lower_bound <- 1 + fco[i]
      upper_bound <- min( c( 1 + fco[i + 1], nrow(sp$S) ) )
      energy[i, tm] <- sum( abs(sp$S[ lower_bound:upper_bound, tm ])^2 )
    }
  }
  energy[energy < eps] <- eps
  energy = 10*log10(energy)
  
  energy.avg[j,] <- rowMeans(energy)
  energy.sd[j,] <- sqrt(rowMeans((energy - energy.avg[j,] %o% rep(1,ntm)) ^ 2))
  
  #energy.mat[j,] <- c(energy)
}

# getting labels
labels <- read.table("hw_data/hw_data/Labels.txt")$x


# Coordinate transformation
suppressMessages(require(dimRed, quietly = T))

embed_methods <- c("PCA", "DiffusionMaps", "Isomap", "LLE", "tSNE")
data.emb <- lapply(embed_methods, function(x) embed(energy.avg.centr, x))

plot(data.emb[[2]]@data@data, col=color.map)


library(plotly)
p <- plot_ly(
  x=data.mat[,1],
  y=data.mat[,2],
  z=data.mat[,3],
  size=1,
  color = ~labels
) %>% add_markers()
p

# classification
train.index <- sample(1:n, floor(0.75 * n), replace = F)
train.set <- data.mat[train.index,]
test.set <- data.mat[-train.index,]

# KNN classification
test.pred <- knn(train.set, test.set, labels[train.index], k=3)
conf.matrix <- caret::confusionMatrix(test.pred, labels[-train.index])
conf.matrix

# svm
basic.classifier <- ksvm(train.set, labels[train.index])
prediction <- predict(basic.classifier, test.set)
conf.table <- caret::confusionMatrix(prediction, labels[-train.index])
conf.table

# random forest
rf.classifier <- randomForest(x=train.set, y=labels[train.index], keep.forest = T,
                              xtest=test.set, ytest=labels[-train.index])
prediction <- predict(rf.classifier, test.set)
conf.table <- caret::confusionMatrix(prediction, labels[-train.index])
conf.table
plot(importance(rf.classifier))
