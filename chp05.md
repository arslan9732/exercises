## Exercises and solutions for Chapter 5

### Classification 
For this set of exercises we will be using the gene expression and patient annotation data from the glioblastoma patient. You can read the data as shown below:
```{r,readMLdataEx,eval=TRUE}
library(compGenomRData)
# get file paths
fileLGGexp=system.file("extdata",
                      "LGGrnaseq.rds",
                      package="compGenomRData")
fileLGGann=system.file("extdata",
                      "patient2LGGsubtypes.rds",
                      package="compGenomRData")
# gene expression values
gexp=readRDS(fileLGGexp)
# patient annotation
patient=readRDS(fileLGGann)
```

1. Our first task is to not use any data transformation and do classification. Run the k-NN classifier on the data without any transformation or scaling. What is the effect on classification accuracy for k-NN predicting the CIMP and noCIMP status of the patient? [Difficulty: **Beginner**]

**solution:** If we do not perform data transformation then the classification accuracy is 87% which is low. 
```{r}
library(caret)
tgexp <- t(gexp)
tgexp=merge(patient,tgexp,by="row.names")
# push sample ids back to the row names
rownames(tgexp)=tgexp[,1]
tgexp=tgexp[,-1]

set.seed(3031) # set the random number seed for reproducibility 

# get indices for 70% of the data set
intrain <- createDataPartition(y = tgexp[,1], p= 0.7)[[1]]

# seperate test and training sets
training <- tgexp[intrain,]
testing <- tgexp[-intrain,]

knnFit=knn3(x=training[,-1], # training set
            y=training[,1], # training set class labels
            k=5)
# predictions on the test set
trainPred=predict(knnFit,training[,-1])
testPred=predict(knnFit,testing[,-1])


# --------  model performance

# predictions on the test set, return class labels
testPred=predict(knnFit,testing[,-1],type="class")

# compare the predicted labels to real labels
# get different performance metrics
confusionMatrix(data=testing[,1],reference=testPred)

```


2. Bootstrap resampling can be used to measure the variability of the prediction error. Use bootstrap resampling with k-NN for the prediction accuracy. How different is it from cross-validation for different $k$s? [Difficulty: **Intermediate**]

**solution:** There is 0.9% difference in terms of accuracy between bootstrap resampling and cross validation when we use K-NN method.
```{r}

#log transformation
gexp=log10(gexp+1)


# transpose the matrix
tgexp <- t(gexp)


#----------- data filtering and scaling
library(caret)
# remove near zero variation for the columns at least
# 85% of the values are the same
# this function creates the filter but doesn't apply it yet
nzv=preProcess(tgexp,method="nzv",uniqueCut = 15)

# apply the filter using "predict" function
# return the filtered dataset and assign it to nzv_tgexp
# variable
nzv_tgexp=predict(nzv,tgexp)

# get most variable genes
SDs=apply(tgexp,2,sd )
topPreds=order(SDs,decreasing = TRUE)[1:1000]
tgexp=tgexp[,topPreds]

processCenter=preProcess(tgexp, method = c("center"))
tgexp=predict(processCenter,tgexp)


# Data scaling
# create a filter for removing highly correlated variables
# if two variables are highly correlated only one of them
# is removed
corrFilt=preProcess(tgexp, method = "corr",cutoff = 0.9)
tgexp=predict(corrFilt,tgexp)

tgexp=merge(patient,tgexp,by="row.names")

# push sample ids back to the row names
rownames(tgexp)=tgexp[,1]
tgexp=tgexp[,-1]
dim(tgexp)

set.seed(3031) # set the random number seed for reproducibility 

# get indices for 70% of the data set
intrain <- createDataPartition(y = tgexp[,1], p= 0.7)[[1]]

# seperate test and training sets
training <- tgexp[intrain,]
testing <- tgexp[-intrain,]

set.seed(17)
# this method controls everything about training
# we will just set up 100 bootstrap samples and for each 
# bootstrap OOB samples to test the error
trctrlbt <- trainControl(method = "boot",number=20,
                       returnResamp="all")

# we will now train k-NN model
knn_fit.bt <- train(subtype~., data = training, 
                 method = "knn",
                 trControl=trctrlbt,
                 tuneGrid = data.frame(k=1:12))
knn_fit.bt$results
knn_fit.bt$bestTune

set.seed(17)
# this method controls everything about training
# we will just set up 10-fold cross validation
trctrl <- trainControl(method = "cv",number=10)

# we will now train k-NN model
knn_fit <- train(subtype~., data = training, 
                 method = "knn",
                 trControl=trctrl,
                 tuneGrid = data.frame(k=1:12))

# best k value by cross-validation accuracy
knn_fit$results
knn_fit$bestTune

```

3. There are a number of ways to get variable importance for a classification problem. Run random forests on the classification problem above. Compare the variable importance metrics from random forest and the one obtained from DALEX applied on the random forests model. How many variables are the same in the top 10? [Difficulty: **Advanced**]

**solution:** There is no similar variable found in top 10.
```{r}
set.seed(17)

# we will do no resampling based prediction error
# although it is advised to do so even for random forests
trctrl <- trainControl(method = "none")

# we will now train random forest model
rfFit <- caret::train(subtype~., 
               data = training, 
               method = "ranger",
               trControl=trctrl,
               importance="permutation", # calculate importance
               tuneGrid = data.frame(mtry=100,
                                     min.node.size = 1,
                                     splitrule="gini")
               )
# print OOB error
#rfFit$finalModel$prediction.error

plot(varImp(rfFit),top=10)


library(DALEX)
set.seed(102)
# do permutation drop-out
explainer_knn<- DALEX::explain(knn_fit, 
                               label="knn", 
                               data =training[,-1], 
                               y = as.numeric(training[,1]))

viknn=feature_importance(explainer_knn,n_sample=50,type="difference")

```



4. Come up with a unified importance score by normalizing importance scores from random forests and DALEX, followed by taking the average of those scores. [Difficulty: **Advanced**]

**solution:**

```{r}
coming soon
```

### Regression
For this set of problems we will use the regression data set where we tried to predict the age of the sample from the methylation values. The data can be loaded as shown below: 
```{r, readMethAgeex,eval=TRUE}
# file path for CpG methylation and age
fileMethAge=system.file("extdata",
                      "CpGmeth2Age.rds",
                      package="compGenomRData")
# read methylation-age table
ameth=readRDS(fileMethAge)
```

1. Run random forest regression and plot the importance metrics. [Difficulty: **Beginner**]

**solution:**
```{r}
# Data cleaning remove low variation data and less than 0.1 sd
library(caret)
ameth=ameth[,c(TRUE,matrixStats::colSds(as.matrix(ameth[,-1]))>0.1)]
dim(ameth)

set.seed(17)

par(mfrow=c(1,2))

# we are not going to do any cross-validation
# and rely on OOB error
trctrl <- trainControl(method = "none")

# we will now train random forest model
rfregFit <- caret::train(Age~., 
               data = ameth, 
               method = "ranger",
               trControl=trctrl,
               # calculate importance
               importance="permutation", 
               tuneGrid = data.frame(mtry=50,
                                     min.node.size = 5,
                                     splitrule="variance")
               )
# plot Observed vs OOB predicted values from the model
plot(ameth$Age,rfregFit$finalModel$predictions,
     pch=19,xlab="observed Age",
     ylab="OOB predicted Age")
mtext(paste("R-squared",
            format(rfregFit$finalModel$r.squared,digits=2)))

# plot residuals
plot(ameth$Age,(rfregFit$finalModel$predictions-ameth$Age),
     pch=18,ylab="residuals (predicted-observed)",
     xlab="observed Age",col="blue3")
abline(h=0,col="red4",lty=2)
rfregFit$metric

```

2. Split 20% of the methylation-age data as test data and run elastic net regression on the training portion to tune parameters and test it on the test portion. [Difficulty: **Intermediate**] 

**solution:**
```{r}
set.seed(17)
library(glmnet)
# this method controls everything about training
# we will just set up 10-fold cross validation
intrain <- createDataPartition(y = ameth[,1], p= 0.8)[[1]]

# seperate test and training sets
training <- ameth[intrain,]
testing <- ameth[-intrain,]


trctrl <- trainControl(method = "cv",number=10)
# we will now train elastic net model
# it will try
enetFit <- train(Age~., data = training,
method = "glmnet",
trControl=trctrl,
# alpha and lambda paramters to try
tuneGrid = data.frame(alpha=0.5,
lambda=seq(0.1,0.7,0.05)))
# best alpha and lambda values by cross-validation accuracy
enetFit$bestTune

class.res=predict(enetFit,testing[,-1])

```

3. Run an ensemble model for regression using the **caretEnsemble** or **mlr** package and compare the results with the elastic net and random forest model. Did the test accuracy increase?
**HINT:** You need to install these extra packages and learn how to use them in the context of ensemble models. [Difficulty: **Advanced**] 

**solution:**
```{r}
coming soon
```
