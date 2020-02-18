dev.off() # close plots
rm(list=ls()) # wipe environment
library(DataExplorer)
library(ggplot2)
library(class)
library(psych)
library(corrplot)
library(factoextra)
library(gmodels)
library(caret)
library(reshape2)
library(PerformanceAnalytics)
library(GGally)
library(gridExtra)
library(randomForest)
library(gbm)
library(e1071)
library(highcharter)

bcancer <- read.csv("E:/College/Analytics/R/data.csv", stringsAsFactors = FALSE)
dim(bcancer) 

#Data preparation
bcancer<- bcancer[-c(1,33)] #bcancer<-bcancer[,-1] #bcancer$X <- NULL
table(bcancer$diagnosis)
bcancer$diagnosis <- factor(ifelse(bcancer$diagnosis=="B","Benign","Malignant"))
str(bcancer)
head(bcancer)
summary(bcancer)
describe(bcancer)
summary(is.na(bcancer))
plot_missing(bcancer)

ggplot(data = melt(bcancer, id.var = "diagnosis"), mapping = aes(x = value)) + 
  geom_histogram(bins = 10, aes(fill=diagnosis), alpha=0.6) +
  facet_wrap(~variable, scales ='free_x')

# Correlation between variables
plot_correlation(bcancer)
corrplot(cor(bcancer[,2:31]), order = "hclust", tl.cex = 0.6)

# Applied different function to each data (mean, se, worst)
chart.Correlation(bcancer[,c(2:11)],histogram=TRUE, col="grey10", pch=1, main="Cancer Mean")

pairs.panels(bcancer[,c(12:21)], method="pearson",
             hist.col = "#1fbbfa", density=TRUE, ellipses=TRUE, show.points = TRUE,
             pch=1, lm=TRUE, cex.cor=1, smoother=F, stars = T, main="Cancer SE")

ggpairs(bcancer[,c(22:31)],)+ theme_bw()+
  labs(title="Cancer Worst")+
  theme(plot.title=element_text(face='bold',color='black',hjust=0.5,size=13))

ggpairs(bcancer[,c(2:11,1)], aes(color=diagnosis, alpha=0.75), lower=list(continuous="smooth"))+ theme_bw()+
  labs(title="Cancer Mean")+
  theme(plot.title=element_text(face='bold',color='black',hjust=0.5,size=12))

ggcorr(bcancer[,c(2:11)], name = "corr", label = TRUE)+
  theme(legend.position="none")+
  labs(title="Cancer Mean")+
  theme(plot.title=element_text(face='bold',color='black',hjust=0.5,size=12))
ggcorr(bcancer[,c(12:21)], name = "corr", label = TRUE)+
  theme(legend.position="none")+
  labs(title="Cancer Mean")+
  theme(plot.title=element_text(face='bold',color='black',hjust=0.5,size=12))
ggcorr(bcancer[,c(22:31)], name = "corr", label = TRUE)+
  theme(legend.position="none")+
  labs(title="Cancer Mean")+
  theme(plot.title=element_text(face='bold',color='black',hjust=0.5,size=12))

bcan_n <-scale(bcancer[-1])
summary(bcan_n)

bcan_pca <- prcomp(bcan_n, center=TRUE,  scale = TRUE)
summary(bcan_pca)

plot(bcan_pca, type="l", main='',color='red')
grid(nx = 10, ny = 14)
title(main = "Principal components weight", sub = NULL, xlab = "Components")
box()

mean_pca <- prcomp(bcan_n[,c(1:10)], scale = TRUE)
summary(mean_pca)

#Percentage of variability
fviz_eig(bcan_pca, addlabels=TRUE, ylim=c(0,60), geom = c("bar", "line"), barfill = "pink", barcolor="grey",linecolor = "red", ncp=10)+
  labs(title = "Cancer All Variances - PCA",
       x = "Principal Components", y = "% of variances")

#PCA Variables
var<-get_pca_var(bcan_pca)
corrplot(var$cos2, is.corr=FALSE)
corrplot(var$contrib, is.corr=FALSE)

# Top 10 features in PCA components
p1 <- fviz_contrib(bcan_pca, choice="var", axes=1, fill="lightgreen", color="grey", top=10)
p2 <- fviz_contrib(bcan_pca, choice="var", axes=2, fill="skyblue", color="grey", top=10)
grid.arrange(p1,p2,ncol=2)

#Train and Test split
set.seed(218)
index <- sample(1:nrow(bcancer), size=0.70*nrow(bcancer))
train <- bcancer[index,]
test <- bcancer[-index,]

prop.table(table(train$diagnosis))
# Random Forest 
learn_rf <- randomForest(diagnosis~., data=train, ntree=500, proximity=T, importance=T)
pre_rf   <- predict(learn_rf, test[,-1])
cm_rf<-confusionMatrix(pre_rf, test$diagnosis)
fourfoldplot(cm_rf$table, conf.level = 0, margin = 1, main=paste("RandomForest (",round(cm_rf$overall[1]*100),"%)",sep=""))

#GBM
test_gbm <- gbm(diagnosis~., data=train, distribution="gaussian",n.trees = 10000,
                shrinkage = 0.01, interaction.depth = 4, bag.fraction=0.5, train.fraction=0.5,n.minobsinnode=10,cv.folds=3,keep.data=TRUE,verbose=FALSE,n.cores=1)
best.iter <- gbm.perf(test_gbm, method="cv",plot.it=FALSE)
fitControl = trainControl(method="cv", number=5, returnResamp="all")
learn_gbm = train(diagnosis~., data=train, method="gbm", distribution="bernoulli", trControl=fitControl, verbose=F, tuneGrid=data.frame(.n.trees=best.iter, .shrinkage=0.01, .interaction.depth=1, .n.minobsinnode=1))
pre_gbm <- predict(learn_gbm, test[,-1])
confusionMatrix(pre_gbm, test$diagnosis)

#SVM
learn_svm <- svm(diagnosis~., data=train)
pre_svm <- predict(learn_svm, test[,-1])
confusionMatrix(pre_svm, test$diagnosis)

#KNN
acc_test <- numeric() 

for(i in 1:30){
  predict <- knn(train=train[,-1], test=test[,-1], cl=train[,1], k=i, prob=T)
  acc_test <- c(acc_test,mean(predict==test[,1]))
}

acc <- data.frame(k= seq(1,30), cnt = acc_test)

opt_k <- subset(acc, cnt==max(cnt))[1,]
sub <- paste("Optimal number of k is", opt_k$k, "(accuracy :", opt_k$cnt,") in KNN")

hchart(acc, 'line', hcaes(k, cnt)) %>%
  hc_title(text = "Accuracy With Varying K (KNN)") %>%
  hc_subtitle(text = sub) %>%
  hc_add_theme(hc_theme_google()) %>%
  hc_xAxis(title = list(text = "Number of Neighbors(k)")) %>%
  hc_yAxis(title = list(text = "Accuracy"))

pre_knn <- knn(train = train[,-1], test = test[,-1], cl = train[,1], k=opt_k$k, prob=T)
confusionMatrix(pre_knn, test$diagnosis)


