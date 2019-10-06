library(wrassp)
library(tuneR)

library(class)
library(randomForest)
library(kernlab)

suppressMessages(require(dimRed, quietly = T))


# DATA PROCESSING

n <- 150

# mfcc parameters
ncep = 13
winsize = 2.
hoptime = winsize / 2

# getting data and applying mfcc
mfcc.result <- list()

mfcc.avg <- matrix(rep(0,n * ncep), nrow=n, ncol=ncep)
mfcc.sd <- matrix(rep(0,n * ncep), nrow=n, ncol=ncep)

for (j in 1:n) {
  song <- read.AsspDataObj(paste("hw_data/hw_data/f", j, ".au", sep=""))
  song.wv <- Wave( as.numeric(song$audio), samp.rate = rate.AsspDataObj(song), bit = 16)
  #downsample_rate = 11025
  #song.wv.dwn = downsample(song.wv, samp.rate = downsample_rate)

  mfcc.result[[j]] <- melfcc(song.wv, wintime=winsize, hoptime=hoptime, numcep=ncep)
  
  mfcc.avg[j,] <- colMeans(mfcc.result[[j]])
  mfcc.sd[j,] <- sqrt(colMeans((mfcc.result[[j]] 
                                - rep(1, dim(mfcc.result[[j]])[1]) %o% mfcc.avg[j,]) ^ 2))
}


# getting labels
labels <- read.table("hw_data/hw_data/Labels.txt")$x

# cutting data into a fixed length to form a matrix
dim1 <- function(mat) dim(mat)[1]
mfcc.lengths <- unlist(lapply(mfcc.result, dim1))
min.mfcc.length <- min(mfcc.lengths)

mfcc.mat <- matrix(rep(0, n * ncep * min.mfcc.length), nrow=n, ncol=ncep * min.mfcc.length)
for (j in 1:n) {
  mfcc.mat[j,] <- c(mfcc.result[[j]][1:min.mfcc.length,])
}


# CLASSIFICATION
data.mat <- mfcc.mat

npars <- 15
nfolds <- 5

sigmas <- seq(0.1, 40, length.out = npars)
fold.index <- caret::createFolds(labels, k=nfolds)

knn.inaccuracy <- rep(0,npars)
svm.inaccuracy <- rep(0,npars)

for (k in 1:npars) {
  for (ff in 1:nfolds) {
    # splitting data
    train.set <- data.mat[-fold.index[[ff]],]
    test.set <- data.mat[fold.index[[ff]],]
    
    train.labels <- labels[-fold.index[[ff]]]
    test.labels <- labels[fold.index[[ff]]]
    
    # KNN classifier
    knn.test <- knn(train.set, test=test.set, cl=train.labels, k=k)
    knn.inaccuracy[k] <- knn.inaccuracy[k] + 1 / nfolds * mean(test.labels != knn.test)
  
    # SVM classifier
    svm.class <- ksvm(train.set, train.labels, sigma=sigmas[k])
    svm.test <- predict(svm.class, test.set)
    svm.inaccuracy[k] <- svm.inaccuracy[k] + 1 / nfolds * mean(test.labels != svm.test)
  }
}


plot(1:npars, knn.inaccuracy, main="KNN accuracy")
plot(sigmas, svm.inaccuracy, main="SVM accuracy")

best.k <- which.min(knn.inaccuracy)
best.sigma <- which.min(svm.inaccuracy)


# creating training-testing partition
train.index <- sample(1:n, floor(0.8 * n), replace = F)
train.set <- data.mat[train.index,]
test.set <- data.mat[-train.index,]

# KNN classification
test.pred <- knn(train.set, test.set, labels[train.index], k=3)
conf.matrix <- caret::confusionMatrix(test.pred, labels[-train.index])
conf.matrix

# SVM classification
svm <- ksvm(train.set, labels[train.index])
test.pred <- predict(svm, test.set)
conf.matrix <- caret::confusionMatrix(test.pred, labels[-train.index])
conf.matrix

# random forest classification
rf.classifier <- randomForest(x=train.set, y=labels[train.index], 
                              xtest=test.set, ytest=labels[-train.index], keep.forest = T)
prediction <- predict(rf.classifier, test.set)
conf.table <- caret::confusionMatrix(prediction, labels[-train.index])
conf.table


# VISUALIZATION

# coordinate transformation
embed_methods <- c("PCA", "DiffusionMaps", "Isomap", "LLE", "tSNE")
data.emb <- lapply(embed_methods, function(x) embed(mfcc.sd, x))

plot(data.emb[[5]]@data@data, col=color.map)


# building color map
genre.names <- levels(labels)
color.list <- c("black", "red", "blue", "green", "cyan")
color.fun <- function(label) color.list[which(genre.names == label)]
color.map <- sapply(labels, color.fun)

plot(mfcc.svd.proj, col=color.map)