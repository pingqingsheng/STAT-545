dir <- "C:\\Users\\zsz\\Desktop\\STAT-541 Multi Var\\Final Project"
setwd(dir)

library(data.table)
library(bit64)
library(MASS)
library(glmnet)
library(e1071)
library(utils)
library(ROCR)
library(mclust)
library(ggplot2)

# convert -coalesce tsne.gif tsne.png

a_train <- as.data.frame(fread("a_Train.csv", header=TRUE))
b_train <- as.data.frame(fread("b_Train.csv", header=TRUE))
a_test <- as.data.frame(fread("a_Test.csv", header=TRUE))
b_test <- as.data.frame(fread("b_Test.csv", header=TRUE))

Dat_Train <- rbind(a_train,b_train)[,-1]
Dat_Test  <- rbind(a_test, b_test)[,-1]
# -------------------------- Tell a FEG from a neutral FE ----------------------------
plotdir <- "C:\\Users\\zsz\\Desktop\\STAT-541 Multi Var\\Presentation"
setwd(plotdir)

png("PC1_center.png", width=1000, height=700, res=100, units="px")
M_0 <- prcomp(Dat_Train[,-1], center=TRUE, scale=TRUE)
score_1 <- M_0$x[,1]
score_2 <- M_0$x[,2]
plot(score_1~ score_2, col=(Dat_Train[,1]+1), pch=19, 
     main="Scatter: PC1 ~ PC2", xlab="PC2", ylab="PC1")
legend("topright", legend = c("Invalid FE","GFE"), col=c(1,2), pch=19)
dev.off()

# library(plot3D)
# score_3 <- M_0$x[,3]
# test <- plot3d(score_1, score_2, score_3, type="p", col=(Dat_Train[,1]+1))
# 

#  Projection on mean between target and neutral FE
bar <- aggregate(  .~as.factor(Dat_Train[,1]) , FUN=mean, data=Dat_Train)
bar <- bar[1,-1] - bar[2,-1]
bar <- bar/sqrt(sum(bar^2))
proj_0 <- as.matrix(Dat_Train[Dat_Train$Target==0,]) %*% t(as.matrix(bar))
proj_1 <- as.matrix(Dat_Test[Dat_Test$Target==1,]) %*% t(as.matrix(bar))

png("Proj.png", width=1000, height=700, res=100, units="px")
plot(density(proj_0), col="red", lwd=2, main="Density : Projected Points")
lines(density(proj_1), col="blue", lwd=2)
legend("topright", legend=c("Invalid FE","GFE"), lwd=2, col=c("red","blue"), lty=1)
dev.off()

# LDA and QDA
system.time(M_lda <- lda(Target ~., data=Dat_Train))
pred  <- predict(M_lda, Dat_Test[,-1])
(confusion <- table(pred$class, Dat_Test[,1]))
(accuracy  <- sum(diag(confusion))/sum(confusion))
# 0.8624743
# tuning and cv............
system.time(M_qda <- qda(Target ~., data=Dat_Train))
pred  <- predict(M_qda, Dat_Test[,-1])
(confusion <- table(pred$class, Dat_Test[,1]))
(accuracy  <- sum(diag(confusion))/sum(confusion))
# 0.7252424
# tuning and cv .............
system.time(M_log <- glm(Target ~., data=Dat_Train, family=binomial()))
pred  <- exp(predict(M_log, Dat_Test[,-1]))/(1+exp(predict(M_log, Dat_Test[,-1])))
pred[pred>0.5] <- 1
pred[pred<0.5] <- 0
(confusion <- table(pred, Dat_Test[,1]))
(accuracy  <- sum(diag(confusion))/sum(confusion))
# 0.8624743
# tuning and cv............

system.time(M_lasso <- cv.glmnet(x=as.matrix(Dat_Train[,-1]), y=as.matrix(Dat_Train[,1]),
                                 family="binomial", nfolds=10))
load("LASSO.RData")
pred    <- predict(M_lasso, as.matrix(Dat_Test[,-1]), type="class", s="lambda.min")
(confusion <- table(pred, Dat_Test[,1]))
(accuracy  <- sum(diag(confusion))/sum(confusion))
# 0.8442551
Coef_sig <- M_lasso$glmnet.fit$beta
times    <- apply(Coef_sig, 1, function(x) length(which(!(x==0))))

library(jpeg)
face <- jpeg::readJPEG(".\\Points.jpg")
face <- abs(face-1)
# 852 * 758 * 3
par(mar = rep(0, 4))
# plot.window(c(0,1),c(0,1), xaxs = "i", yaxs = "i")
plot(0:1, 0:1, type="n", ann=FALSE, axes=FALSE)
rasterImage(face, 0, 0, 1, 1)
# write.csv(Coor, "Coor.csv")
Coor <- read.csv("Coor.csv", header=TRUE)[,-1]
Coor <- rbind(Coor[100,],Coor[-100,])
Effect <- times[3*(seq_len(101)-1)]
Coor_x <- rep(Coor$X, times=Effect)
Coor_y <- rep(Coor$Y, times=Effect)
Coor   <- data.frame(X=Coor_x, Y=Coor_y)
# Coor   <- data.frame(x=Coor[,1], y=Coor[,2], Effect = Effect)
# grid_heat <- expand.grid(x=seq(0,1,length.out=100), y=seq(1,0,length.out=100))
# Effect    <- rep(0, 100*100)
# for (i in 1:100){
#   Inx <- rownames(grid_heat[round(grid_heat$x,digits=2) %in% round(Coor$x[i],digits=2) &
#                             round(grid_heat$y,digits=2) %in% round(Coor$y[i],digits=2),])
#   Inx <- as.numeric(Inx)
#   Effect[Inx] <- Coor$Effect[i]
#   cat("\r", paste0(i,"%"))
# }

# grid_heat <- cbind(grid_heat, Effect)

# library(akima)
# spline <- interp(grid_heat$x, grid_heat$y, grid_heat$Effect, linear=FALSE)
# grid_Heat <- expand.grid(spline$x, spline$y)
# grid_Heat <- data.frame(x=grid_heat[,1], y=grid_heat[,2], z=as.vector(spline$z))

back_groud_R <- as.vector(face[,,1])
back_groud_G <- as.vector(face[,,2])
back_groud_B <- as.vector(face[,,3])
grid_RGB     <- expand.grid(seq(0,1,length.out = 852), seq(1,0,length.out = 758))
face_raster  <- data.frame(x=grid_RGB[,2], y=grid_RGB[seq(645816,1),1],
                           R=back_groud_R, G=back_groud_G, B=back_groud_B)
# write.csv(face_raster, "face_raster.csv")
png("GFE vs FE.png", width=1000, height=700, res=100, units="px")
ggplot()+
  geom_point(data=face_raster, aes(x=x,y=y, col=rgb(R,G,B)))+
  scale_color_identity()+
  stat_density2d(data=Coor, aes(x=X, y=Y, fill=..level.., alpha=..level..), 
                 size=2, bins=5, geom="polygon")+
  scale_fill_gradient(low="transparent",high="red", guide=FALSE)+
  scale_alpha(range = c(0, 0.8), guide = FALSE)+
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )
dev.off()

# naiveBayes
M_naive <- naiveBayes(as.factor(Target) ~., data=Dat_Train)
system.time(pred    <- predict(M_naive, as.matrix(Dat_Test[,-1])))
(confusion <- table(pred, Dat_Test[,1]))
(accuracy  <- sum(diag(confusion))/sum(confusion))
pref    <- prediction(as.numeric(pred), as.numeric(Dat_Test[,1]))
plot(performance(pref, "tpr", "fpr"))

# SVM
system.time(M_svm <- svm(Target ~., data=Dat_Train, type="C"))
pred  <- predict(M_svm, Dat_Test[,-1], type="class")
(confusion <- table(pred, Dat_Test[,1]))
(accuracy  <- sum(diag(confusion)))/sum(confusion)
# 0.8862768

# ada
library(ada)
M_ada <- ada(x=Dat_Train[,-1], y=Dat_Train[,1], loss="exponential", type="discrete", iter=100)
pred  <- predict(M_ada, newdata=data.frame(Dat_Test[,-1]), type="probs")
(confusion <- table(pred, Dat_Test[,1]))
(accuracy  <- sum(diag(confusion)))/sum(confusion)

# ---------------------------------  Separate Different People  -------------------------------------
Train <- rbind(cbind(cate="a",a_train[,-c(1,2)]), 
               cbind(cate="b",b_train[,-c(1,2)]))
Test  <- rbind(cbind(cate="a",a_test[,-c(1,2)]), 
               cbind(cate="b",b_test[,-c(1,2)]))

#  Projection on mean between a and b
bar <- apply(a_train[,-c(1,2)], 2, mean) -  apply(b_train[,-c(1,2)], 2, mean)
bar <- bar/sqrt(sum(bar^2))
proj_a <- as.matrix(a_train[,-c(1,2)]) %*% (bar)
proj_b <- as.matrix(b_train[,-c(1,2)]) %*% (bar)
png("Proj2.png", width=1000, height=700, res=150, units="px")
plot(density(proj_a), col="red", lwd=2, xlim=c(min(proj_a, proj_b), max(proj_a, proj_b)),
     main="Density: Projected Points")
lines(density(proj_b), col="blue", lwd=2)
legend("topleft", legend=c("A user", "B user"), lty=1, lwd=2, col=c("red","blue"))
dev.off()

# First Two PC
png("PC_AB.png", width=1000, height=700, res=150, units="px")
M_0 <- prcomp(Train[,-1], center=FALSE, scale=TRUE)
score_1 <- M_0$x[,1]
score_2 <- M_0$x[,2]
plot(score_1~ score_2, col=(as.numeric(Train[,1])), pch=19, 
     main="Scatter: PC1 ~ PC2 (A VS B)", xlab="PC2", ylab="PC1")
legend("topright", legend = c("a","b"), col=c(2,1), pch=19)
dev.off()

png("PC_AB_standardizing.png", width=1000, height=700, res=150, units="px")
M_0 <- prcomp(Train[,-1], center=TRUE, scale=TRUE)
score_1 <- M_0$x[,1]
score_2 <- M_0$x[,2]
plot(score_1~ score_2, col=(as.numeric(Train[,1])), pch=19, 
     main="Scatter: PC1 ~ PC2 (A VS B)", xlab="PC2", ylab="PC1")
legend("topright", legend = c("a","b"), col=c(2,1), pch=19)
dev.off()

# LDA and QDA
library(MASS)
system.time(M_lda <- lda(cate~., data=Train))
pred  <- predict(M_lda,Test[,-1])
(confusion <- table(pred$class, Test[,1]))
(accuracy  <- sum(diag(confusion))/sum(confusion))
# 0.9995705

system.time(M_qda <- qda(cate~., data=Train))
pred  <- predict(M_qda, Test[,-1])
(confusion <- table(pred$class, Test[,1]))
(accuracy  <- sum(diag(confusion))/sum(confusion))
# 0.9756621
# Sometimes when you don't tune or make extra design, you get really good result
# maybe you're lucky, or maybe it's not gonna to get a good score

# Logistic and logistic LASSO
system.time(M_logistic <- glm(cate~., data=Train, family=binomial))
pred <- exp(predict(M_logistic, Test[,-1]))/(1+exp(predict(M_logistic, Test[,-1])))
pred[pred>0.5] <- 1
pred[pred<0.5] <- 0
(confusion <- table(pred, Test[,1]))
(accuracy  <- sum(diag(confusion))/sum(confusion))
# 0.9997137
library(glmnet)
system.time(M_logLasso <- cv.glmnet(x=as.matrix(Train[,-1]), y=as.matrix(Train[,1]),
                                    family="binomial"))
pred <- predict(M_logLasso, as.matrix(Test[,-1]), type="class", s="lmabda.min")
(confusion <- table(pred, Test[,1]))
(accuracy  <- sum(diag(confusion))/sum(confusion))

load("M_Lasso.RData")
Coef_sig <- M_logLasso$glmnet.fit$beta
times    <- apply(Coef_sig, 1, function(x) length(which(!(x==0))))

face <- jpeg::readJPEG(".\\Points.jpg")
face <- abs(face-1)
# 852 * 758 * 3
par(mar = rep(0, 4))
# plot.window(c(0,1),c(0,1), xaxs = "i", yaxs = "i")
plot(0:1, 0:1, type="n", ann=FALSE, axes=FALSE)
rasterImage(face, 0, 0, 1, 1)
# write.csv(Coor, "Coor.csv")
Coor <- read.csv("Coor.csv", header=TRUE)[,-1]
Coor <- rbind(Coor[100,],Coor[-100,])
Effect <- times[3*(seq_len(101)-1)]
Coor_x <- rep(Coor$X, times=Effect)
Coor_y <- rep(Coor$Y, times=Effect)
Coor   <- data.frame(X=Coor_x, Y=Coor_y)

back_groud_R <- as.vector(face[,,1])
back_groud_G <- as.vector(face[,,2])
back_groud_B <- as.vector(face[,,3])
grid_RGB     <- expand.grid(seq(0,1,length.out = 852), seq(1,0,length.out = 758))
face_raster  <- data.frame(x=grid_RGB[,2], y=grid_RGB[seq(645816,1),1],
                           R=back_groud_R, G=back_groud_G, B=back_groud_B)

png("AvsB.png", width=1000, height=700, res=100, units="px")
ggplot()+
  geom_point(data=face_raster, aes(x=x,y=y, col=rgb(R,G,B)))+
  scale_color_identity()+
  stat_density2d(data=Coor, aes(x=X, y=Y, fill=..level.., alpha=..level..), 
                 size=2, bins=5, geom="polygon")+
  scale_fill_gradient(low="transparent",high="red", guide=FALSE)+
  scale_alpha(range = c(0, 0.8), guide = FALSE)+
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )
dev.off()

# SVM
library(e1071)
system.time(M_svm <- svm(cate ~., data=Train, type="C"))
pred  <- predict(M_svm, Test[,-1])
# save(M_svm,file="M_svm.RData")
(confusion <- table(pred, Test[,1]))
(accuracy  <- sum(diag(confusion))/sum(confusion))

# naive bayes
system.time(M_naivebayes <- naiveBayes(x=Train[,-1], y=Train[,1], newdata=Test[,-1]))
system.time(pred <- predict(M_naivebayes, Test[,-1]))
(confusion <- table(pred, Test[,1]))
(accuracy  <- sum(diag(confusion))/sum(confusion))
# 0.9932713
# Even naive bayes can do 0.9932713

# ada
system.time(M_ada <- ada(x=Dat_Train[,-1], y=Dat_Train[,1], max.iter = 100,
                         loss="exponential", type="discrete"))
# save(M_ada, file="M_ada.RData")
pred  <- predict(M_ada, Dat_Test[,-1], type="probs")
(confusion <- table(pred, Test[,1]))
(accuracy  <- sum(diag(confusion))/sum(confusion))


# ----------------------------- Tell Different Facial Expression -----------------------------
dir <- "C:\\Users\\zsz\\Desktop\\STAT-541 Multi Var\\Final Project\\"
setwd(dir)

library(data.table)
library(bit64)
library(rgl)
library(e1071)
library(rpart)
library(caret)
library(randomForest)
library(xgboost)
library(ada)
library(VGAM)

Dat <- fread("Categorical.csv", header=TRUE)
Dat <- as.data.frame(Dat)[,-1]

png("PC_cate.png", width=1000, height=700, res=150, units="px")
pc <- prcomp(Dat[,-1])
score_1 <- pc$x[,1]
score_2 <- pc$x[,2]
plot(score_1 ~ score_2, pch=19, col=(Dat[,1]+1))
dev.off()

png("PC_cate_standardize.png", width=1000, height=700, res=150, units="px")
pc <- prcomp(Dat[,-1], center=TRUE, scale=TRUE)
score_1 <- pc$x[,1]
score_2 <- pc$x[,2]
plot(score_1 ~ score_2, pch=19, col=(Dat[,1]+1))
dev.off()

Train <- as.data.frame(fread("Category_Train.csv", header=TRUE))[,-1]
Test  <- as.data.frame(fread("Category_Test.csv", header=TRUE))[,-1]

system.time(M_svm <- svm(cate ~. , data=Train, type="C"))
pred  <- predict(M_svm, Test[,-1])
(confusion <- table(pred, Test[,1]))
(accracy   <- sum(diag(confusion))/sum(confusion))
# 0.86

# If you try to tune a svm, then basically I submit the task to CLEAR and two days later
# it is still there.
# M_svm <- tune(svm, train.x = Train[c(1:1000),-1],  train.y = Train[c(1:1000),1],
#               kernel="radial", range=list(cost=10^c(-1:2), gamma=c(0.5,1,2)))
# print(M_svm)


system.time(M_rpart <- rpart(as.factor(cate) ~. , data=Train,  method="class", xval=5, cp=0.01))
pred <- apply( predict(M_rpart, Test[,-1]), 1, which.max)
plot(M_rpart)
text(M_rpart, use.n=TRUE, cex=0.5)
(confusion <- table(pred, Test[,1]))
(accuracy <- sum(diag(confusion))/sum(confusion))
#  0.6062428
# tune .................

# What about a ada boosting 
# M_ada <- ada(cate ~., )
# ada boost here is not available since it' a binary classifier

# bagging took a lot of time and memery, just leave it there
# M_bagging <- randomForest(cate ~., mtr= 300, data=Train, 
#                           xtest=Test[,-1], ytest=Test[,1])


# What about a  random forest
library(bigmemory)
library(doParallel)

cl <- makeCluster(20)
registerDoParallel(cl)

system.time(my_randomforest <- foreach(
  ntree=rep(50,20), .combine=combine, .multicombine=TRUE, .packages="randomForest"
) %dopar%{
  randomForest(as.factor(cate) ~., ntree=ntree, data=Train, 
               xtest=Test[,-1], ytest=as.factor(Test[,1]),
               type="classification")
})

registerDoSEQ()
stopCluster(cl)

load("rforest.RData")
(confusion <- table(rforest$test$predicted, Test[,1]))
(accuracy   <- sum(diag(confusion))/sum(confusion))
# accuracy = 0.9752291
# You don't tune, and you get a 0.9752291 accuracy, which means you're working on a 
# simple data set. 


# Let's try a stochastic gradient boosting way
Train <- as.big.matrix(Train)
Test  <- as.big.matrix(Test)

cl <- makeCluster(20)
registerDoParallel(cl)

system.time(my_xgboost <- xgboost(data=Train[,-1], label=Train[,1], max.depth=2, eta=1, 
                           objective="multi:softmax", nrounds = 100, nthread=20,
                           num_class=10))

registerDoSEQ()
stopCluster(cl)

load("my_xgboost.RData")
(confusion <- table(predict(my_xgboost, Test[,-1]), Test[,1]))
(accuracy  <- sum(diag(confusion))/sum(confusion))

# You don't tune and you get 0.9682131 accuracy again, which means you're working on a 
# easy data set

# Let's try multinomial again with glmnet
Tra <- sparse.model.matrix(~., Train[,-1])
Tes <- sparse.model.matrix(~., Test[,-1])

load("M_MLasso.RData")
My_lasso <- cv.glmnet(Train[,1] , Tra, family="multinomial")
pred <- predict(My_lasso, Test, type="class", s="lambda.min")
(confusion <- table(predict(my_xgboost, Test[,-1]), Test[,1]))
(accuracy  <- sum(diag(confusion))/sum(confusion))

# I should try neurla network, but I know nothing about it, maybe next time 

# Even I myself forget what is the code below doing, just leave it there. 
convert <- "F:\\ApplySoftwareBase\\ImageMagick-7.0.3-Q16\\ImageMagick-7.0.3-Q16\\convert.exe"
ani.options("convert"=convert)

op <- options(digits.secs = 6)
options(op)
format(Sys.time(), "%s")


# ----------------------------- Dimension Reduction and Visualization --------------------------
# So let's do unsupervised things, to make things convoluted 

library(tsne)

name <- paste0("img",seq_len(1000/5),".png")
i    <- 11
egb <- function(x){
  par(mar = rep(0, 4))
  plot(x, type="p", pch=19, col=as.numeric(Test$cate[c(1:2000)])+1,
       xlim=c(-40,40), ylim=c(-40,40),
       ann=FALSE, axes=FALSE)
}

# saveGIF({
#   iris_tsne <- tsne(as.matrix(iris[,-5]), k=3, perplexity = 40,
#                     epoch_callback = egb, epoch = 5)
# }, interval = 0.2, ani.width=600, ani.height=600)

pc_Test <- prcomp(Test[,-1])
pc_Test <- pc_Test$x[c(1:2000),c(1:30)]

saveGIF({
  Test_tsne <- tsne(pc_Test, k=9, perplexity = 40,
                    epoch_callback = egb, epoch = 5)
}, interval = 0.2, ani.width=600, ani.height=600)

library(Rtsne)

tsne <- Rtsne(as.matrix(Train[c(1:20000),-1]), check_duplicates = FALSE, pca = TRUE,
              perplexity=30, theta=0.5, dims=3)
tsne$Y
embedding <- as.data.frame(tsne$Y)
embedding$Class <- as.numeric(Train[c(1:20000),1])
plot<-plot3d(x=embedding$V1, y=embedding$V2 ,z=embedding$V3,col=embedding$Class, pch = ".")

library(ggplot2)
library(Rtsne)
library(plot3D)
library(rgl)

features <- iris[,-5]

tsne <- Rtsne(as.matrix(features), check_duplicates = FALSE, pca = TRUE,
              perplexity=30, theta=0.5, dims=3)
tsne$Y
embedding <- as.data.frame(tsne$Y)
embedding$Class <- as.numeric(iris[,5])
plot<-plot3d(x=embedding$V1, y=embedding$V2 ,z=embedding$V3,col=embedding$Class, pch = ".")

# ----------------------------------- 3D tsne ---------------------------------------------
# library(Rtsne)
# library(plot3D)
# library(animation)
# 
# 
# 
# color <- rainbow(length(levels(iris[,5])))
# names(color) <- levels(iris[,5])
# 
# saveGIF({
#   for (i in 1:1000){
#     set.seed(123)
#     iris_tsne <- Rtsne(as.matrix(iris[,-5]), dim=3, perplexity=40, pca=FALSE, max_iter=i,
#                        check_duplicates = FALSE,
#                        theta=round(i/1000,digits=1))
#     if(i%%10==0){
#       scatter3D(x=iris_tsne$Y[,1], y=iris_tsne$Y[,2], z=iris_tsne$Y[,3], 
#                 colvar=as.numeric(iris[,5]), pch=19, size=20, type="p", bty="g")
#     }
#     cat("\r",paste0(i/1000*100," %"))
#   }
# }, interval=0.1)






