library(rpart)
library(rpart.plot)
library(FNN)
library(neuralnet)

setwd("C:/dsge/r")
#setwd("C:/temp/2017/dsge/fm")
################################################################
# Multifactor regression #######################################
################################################################
rawdata <- read.csv("input/inputmap.csv", header=TRUE, sep=",", dec=".")
#rawdata <- rawdata[36:nrow(rawdata),]
Xnames <- c("R_","pi_c_","dy_","dc_","di_","dE_") #,"pi_i_","pi_d_","dimp_","dex_"
Ynames <- names(rawdata)[!names(rawdata) %in% c(Xnames,"Year","Quarter","Recession","pi_i_","pi_d_","dimp_","dex_")]

modeloutput <- data.frame(y=character(),
				 RMSE=double(),
                 R2=double(),
                 R2Adjust=double(),
				 df=integer(),
				 fvar=double(),
				 resfvar=double(),
				 rvar=double(),
				 resrvar=double(),
				 corr=double(),
				 rescorr=double(),
                 Models=character(),
                 stringsAsFactors=FALSE)


lag <-2
residuals <- matrix(,nrow=nrow(rawdata)-lag,ncol=length(Ynames))
colnames(residuals) <- Ynames
for (y in Ynames) {
	Traindata <- rawdata[,names(rawdata) %in% c(y,Xnames,"Recession")]
	
	Xnames1 <- names(Traindata)[!names(Traindata) %in% c("Recession")]

	if (lag > 0) {
		for (i in c(1:lag)){
			for (varname in Xnames1){
				Traindata[(i+1):nrow(Traindata),paste0(varname,i)] <- Traindata[1:(nrow(Traindata)-i),varname]
				Traindata[1:i,paste0(varname,i)] <- NA
			}
		}
	}
	
	Traindata <- na.omit(Traindata[(lag+1):nrow(Traindata),])
	
	Xnames2 <- names(Traindata)[!names(Traindata) %in% c(y,"Recession")]
	
	f <-as.formula(paste(y,"~",paste(Xnames2,collapse="+")))
	
	print(f)

	#linear regression
	lr<-lm(f, data=Traindata)
	rmse<-sqrt(mean(lr$residuals^2))
	r2<-1-sum(lr$residuals^2)/sum((Traindata[,y]-mean(Traindata[,y]))^2)
	r2adjust <- r2-(1-r2)*length(Xnames2)/(nrow(Traindata)-length(Xnames2)-1)
	#aic<-AIC(lr)
	#bic<-BIC(lr)
	df<-nrow(Traindata)-length(Xnames2)-1
	fvar <- var(lr$fitted.values)
	resfvar <- var(lr$fitted.values[Traindata$Recession==1])
	rvar <- var(lr$residuals)
	resrvar <- var(lr$residuals[Traindata$Recession==1])
	corr<-cor(lr$residuals,lr$fitted.values)
	rescorr<-cor(lr$residuals[Traindata$Recession==1],lr$fitted.values[Traindata$Recession==1])
	modeloutput[nrow(modeloutput)+1,] <- c(y, rmse, r2, r2adjust, df, fvar, resfvar, rvar, resrvar, corr,rescorr, "LM")
	write.csv(summary(lr)$coefficients,paste0("lr_",y,".csv"))
	residuals[,y]<-c(rep(NA,nrow(residuals)-length(lr$residuals)),lr$residuals)
	
	#cart
	cart = rpart(f, data = Traindata, cp = 10^(-3),minsplit = 10)
	cartpredict <- predict(cart, Traindata)
	rmse<-sqrt(mean((Traindata[,y]-cartpredict)^2))
	r2<-1-sum((Traindata[,y]-cartpredict)^2)/sum((Traindata[,y]-mean(Traindata[,y]))^2)
	r2adjust <- r2-(1-r2)*(length(unique(cart$frame$var))-1)/(nrow(Traindata)-(length(unique(cart$frame$var))-1)-1)
	#rpart.plot(cart)
	df<-nrow(Traindata)-(length(unique(cart$frame$var))-1)-1
	fvar <- var(cartpredict)
	resfvar <- var(cartpredict[Traindata$Recession==1])
	rvar <- var(Traindata[,y]-cartpredict)
	resrvar <- var((Traindata[,y]-cartpredict)[Traindata$Recession==1])
	corr<-cor((Traindata[,y]-cartpredict),cartpredict)
	rescorr<-cor((Traindata[,y]-cartpredict)[Traindata$Recession==1],cartpredict[Traindata$Recession==1])

	modeloutput[nrow(modeloutput)+1,] <- c(y, rmse, r2, r2adjust, df, fvar, resfvar, rvar, resrvar, corr,rescorr, "CART")

	#knn
	knnreg <- knn.reg(train = Traindata[,names(Traindata) %in% Xnames2] , test=Traindata[,names(Traindata) %in% Xnames2], y=Traindata[,y], k=5, algorithm = "kd_tree")
	rmse<-sqrt(mean((Traindata[,y]-knnreg$pred)^2))
	r2<-1-sum((Traindata[,y]-knnreg$pred)^2)/sum((Traindata[,y]-mean(Traindata[,y]))^2)
	r2adjust <- r2-(1-r2)*length(Xnames2)/(nrow(Traindata)-length(Xnames2)-1)
	df<-nrow(Traindata)-length(Xnames2)-1
	fvar <- var(knnreg$pred)
	resfvar <- var(knnreg$pred[Traindata$Recession==1])
	rvar <- var(Traindata[,y]-knnreg$pred)
	resrvar <- var((Traindata[,y]-knnreg$pred)[Traindata$Recession==1])
	corr<-cor((Traindata[,y]-knnreg$pred),knnreg$pred)
	rescorr<-cor((Traindata[,y]-knnreg$pred)[Traindata$Recession==1],knnreg$pred[Traindata$Recession==1])

	modeloutput[nrow(modeloutput)+1,] <- c(y, rmse, r2, r2adjust, df, fvar, resfvar, rvar, resrvar, corr,rescorr,"KNN")

	#ann
	set.seed(123)

	ann <- neuralnet(f, data=data.matrix(Traindata), hidden=c(10,5), linear.output=TRUE, stepmax = 100000, threshold=0.01, act.fct = "tanh", likelihood = TRUE, lifesign ="full", lifesign.step = 1000)
	#ann <- neuralnet(f, data=data.matrix(Traindata), hidden=c(10), linear.output=TRUE, stepmax = 200000, threshold=0.5, act.fct = "tanh", likelihood = TRUE, lifesign ="full", lifesign.step = 100)
	#ann <- neuralnet(f, data=data.matrix(Traindata), hidden=c(5), linear.output=TRUE, stepmax = 200000, threshold=0.5, act.fct = "tanh", likelihood = TRUE, lifesign ="full", lifesign.step = 100)
	#ann <- neuralnet(f, data=data.matrix(Traindata), hidden=c(5), linear.output=TRUE, stepmax = 200000, threshold=0.5, act.fct = "logistic", likelihood = TRUE, lifesign ="full", lifesign.step = 100)
	
	rmse<-sqrt(mean((Traindata[,y]-ann$net.result[[1]][,1])^2))
	r2<-1-sum((Traindata[,y]-ann$net.result[[1]][,1])^2)/sum((Traindata[,y]-mean(Traindata[,y]))^2)
	r2adjust <- r2-(1-r2)*length(Xnames2)/(nrow(Traindata)-length(Xnames2)-1)
	df<-nrow(Traindata)-length(Xnames2)-1
	tryCatch(
	{
		fvar <- var(ann$net.result[[1]][,1])
		resfvar <- var(ann$net.result[[1]][,1][Traindata$Recession==1])
		rvar <- var(Traindata[,y]-ann$net.result[[1]][,1])
		resrvar <- var((Traindata[,y]-ann$net.result[[1]][,1])[Traindata$Recession==1])
		corr<-cor((Traindata[,y]-ann$net.result[[1]][,1]),ann$net.result[[1]][,1])
		rescorr<-cor((Traindata[,y]-ann$net.result[[1]][,1])[Traindata$Recession==1],ann$net.result[[1]][,1][Traindata$Recession==1])
	},
		error = function(ex) {
#			print("errors")
			fvar <- NA
			resfvar <- NA
			rvar <- NA
			resrvar <- NA
			corr <- NA
			rescorr <- NA
		}
	)
	modeloutput[nrow(modeloutput)+1,] <- c(y, rmse, r2, r2adjust, df, fvar, resfvar, rvar, resrvar, corr,rescorr,"ANN")
	
}

	write.csv(modeloutput,paste0("modeloutput1.csv"))
	write.csv(residuals,paste0("residuals1.csv"))

################################################################
# predict recession ############################################
################################################################
rawdata <- read.csv("inputmap.csv", header=TRUE, sep=",", dec=".")
Xnames <- c("pi_c_","dy_","dc_","di_","dE_")
Ynames <- "Recession"

modeloutput <- data.frame(y=character(),
                 Precision=double(),
                 Recall=double(),
                 FMeasure=double(),
                 Models=character(),
                 stringsAsFactors=FALSE)

lag <-2

Traindata <- rawdata[,names(rawdata) %in% c(Ynames,Xnames)]

Xnames1 <- names(Traindata)[!names(Traindata) %in% c(Ynames)]

if (lag > 0) {
	for (i in c(1:lag)){
		for (varname in Xnames1){
			Traindata[(i+1):nrow(Traindata),paste0(varname,i)] <- Traindata[1:(nrow(Traindata)-i),varname]
			Traindata[1:i,paste0(varname,i)] <- NA
		}
	}
}

Traindata <- na.omit(Traindata[(lag+1):nrow(Traindata),])

Xnames2 <- names(Traindata)[!names(Traindata) %in% c(Ynames)]

f <-as.formula(paste(Ynames,"~",paste(Xnames2,collapse="+")))

print(f)

#linear regression
lr<-lm(f, data=Traindata)

predict <- ifelse(lr$fitted.values>0.5, 1, 0)
predictyes <- sum(predict==1)
predictno<-sum(predict==0)
actual <- Traindata[,Ynames]
actualwhenpredyes <- actual[predict==1]
actualwhenpredno <- actual[predict==0]
tp <- sum(actualwhenpredyes==1)
tn <- sum(actualwhenpredno==0)
fp <- length(actualwhenpredyes) - tp
fn <- length(actualwhenpredno) - tn
precision<-tp/predictyes
recall<-tp/(tp+fn)
F = 2*precision*recall/(precision+recall)
paste("The precision is",precision)
paste("The recall is",recall)
paste("The F is",F)

modeloutput[1,] <- c(Ynames,precision, recall, F, "linear")
write.csv(summary(lr)$coefficients,paste0("lr_",Ynames,".csv"))

# Generalized Linear Model
glmr <- glm(f, data=Traindata, family=binomial)
glmrpredict <- predict(glmr, Traindata)

predict <- ifelse(glmrpredict>0.5, 1, 0)
predictyes <- sum(predict==1)
predictno<-sum(predict==0)
actual <- Traindata[,Ynames]
actualwhenpredyes <- actual[predict==1]
actualwhenpredno <- actual[predict==0]
tp <- sum(actualwhenpredyes==1)
tn <- sum(actualwhenpredno==0)
fp <- length(actualwhenpredyes) - tp
fn <- length(actualwhenpredno) - tn
precision<-tp/predictyes
recall<-tp/(tp+fn)
F = 2*precision*recall/(precision+recall)
paste("The precision is",precision)
paste("The recall is",recall)
paste("The F is",F)

modeloutput[2,] <- c(Ynames, precision, recall, F, "glm")
write.csv(glmr$coefficients,paste0("glmr_",Ynames,".csv"))

#cart
cart = rpart(f, data = Traindata, cp = 10^(-3),minsplit = 10,, method = "class")
cartpredict <- predict(cart, Traindata)
predict <- ifelse(cartpredict[,2]>0.5, 1, 0)

predictyes <- sum(predict==1)
predictno<-sum(predict==0)
actual <- Traindata[,Ynames]
actualwhenpredyes <- actual[predict==1]
actualwhenpredno <- actual[predict==0]
tp <- sum(actualwhenpredyes==1)
tn <- sum(actualwhenpredno==0)
fp <- length(actualwhenpredyes) - tp
fn <- length(actualwhenpredno) - tn
precision<-tp/predictyes
recall<-tp/(tp+fn)
F = 2*precision*recall/(precision+recall)
paste("The precision is",precision)
paste("The recall is",recall)
paste("The F is",F)

modeloutput[3,] <- c(Ynames, precision, recall, F, "cart")

#knn
knnreg <- knn.reg(train = Traindata[,names(Traindata) %in% Xnames2] , test=Traindata[,names(Traindata) %in% Xnames2], y=Traindata[,Ynames], k=5, algorithm = "kd_tree")
predict <- ifelse(knnreg$pred>0.5, 1, 0)
predictyes <- sum(predict==1)
predictno<-sum(predict==0)
actual <- Traindata[,Ynames]
actualwhenpredyes <- actual[predict==1]
actualwhenpredno <- actual[predict==0]
tp <- sum(actualwhenpredyes==1)
tn <- sum(actualwhenpredno==0)
fp <- length(actualwhenpredyes) - tp
fn <- length(actualwhenpredno) - tn
precision<-tp/predictyes
recall<-tp/(tp+fn)
F = 2*precision*recall/(precision+recall)
paste("The precision is",precision)
paste("The recall is",recall)
paste("The F is",F)

modeloutput[4,] <- c(Ynames, precision, recall, F, "KNN")

#ann
set.seed(123)

ann <- neuralnet(f, data=data.matrix(Traindata), hidden=c(5), linear.output=TRUE, stepmax = 100000, threshold=0.001, act.fct = "tanh", likelihood = TRUE, lifesign ="full", lifesign.step = 1000)
#ann <- neuralnet(f, data=data.matrix(Traindata), hidden=c(10), linear.output=TRUE, stepmax = 200000, threshold=0.5, act.fct = "tanh", likelihood = TRUE, lifesign ="full", lifesign.step = 100)
#ann <- neuralnet(f, data=data.matrix(Traindata), hidden=c(5), linear.output=TRUE, stepmax = 200000, threshold=0.5, act.fct = "tanh", likelihood = TRUE, lifesign ="full", lifesign.step = 100)
#ann <- neuralnet(f, data=data.matrix(Traindata), hidden=c(5), linear.output=TRUE, stepmax = 200000, threshold=0.5, act.fct = "logistic", likelihood = TRUE, lifesign ="full", lifesign.step = 100)

predict <- ifelse(ann$net.result[[1]][,1]>0.5, 1, 0)
predictyes <- sum(predict==1)
predictno<-sum(predict==0)
actual <- Traindata[,Ynames]
actualwhenpredyes <- actual[predict==1]
actualwhenpredno <- actual[predict==0]
tp <- sum(actualwhenpredyes==1)
tn <- sum(actualwhenpredno==0)
fp <- length(actualwhenpredyes) - tp
fn <- length(actualwhenpredno) - tn
precision<-tp/predictyes
recall<-tp/(tp+fn)
F = 2*precision*recall/(precision+recall)
paste("The precision is",precision)
paste("The recall is",recall)
paste("The F is",F)

modeloutput[5,] <- c(Ynames, precision, recall, F, "ANN")


write.csv(modeloutput,paste0("modeloutput2.csv"))


################################################################
# Variable selection for linear regression #####################
################################################################
rawdata <- read.csv("inputmap.csv", header=TRUE, sep=",", dec=".")
Xnames <- c("R_","pi_c_","pi_i_","pi_d_","dy_","dc_","di_","dimp_","dex_","dE_")
Ynames <- names(rawdata)[!names(rawdata) %in% c(Xnames,"Year","Quarter","Recession")]

modeloutput <- data.frame(y=character(),
				 RMSE=double(),
                 R2=double(),
                 R2Adjust=double(),
				 df=integer(),
				 fvar=double(),
				 resfvar=double(),
				 rvar=double(),
				 resrvar=double(),
				 corr=double(),
				 rescorr=double(),
                 Models=character(),
                 stringsAsFactors=FALSE)


lag <-2
residuals <- matrix(,nrow=nrow(rawdata)-lag,ncol=length(Ynames))
colnames(residuals) <- Ynames
for (y in Ynames) {
	Traindata <- rawdata[,names(rawdata) %in% c(y,Xnames,"Recession")]
	
	Xnames1 <- names(Traindata)[!names(Traindata) %in% c("Recession")]

	if (lag > 0) {
		for (i in c(1:lag)){
			for (varname in Xnames1){
				Traindata[(i+1):nrow(Traindata),paste0(varname,i)] <- Traindata[1:(nrow(Traindata)-i),varname]
				Traindata[1:i,paste0(varname,i)] <- NA
			}
		}
	}
	
	Traindata <- na.omit(Traindata[(lag+1):nrow(Traindata),])
	
	Xnames2 <- names(Traindata)[!names(Traindata) %in% c(y,"Recession")]
	
	f <-as.formula(paste(y,"~",paste(Xnames2,collapse="+")))
	
	print(f)

	#linear regression
	tryCatch(	
	{
		lr<-lm(f, data=Traindata)
		for (i in Xnames2) {
			if (summary(lr)$coefficients[i,4]>0.3) {
				Traindata<-Traindata[,!names(Traindata) %in% c(i)]
			}
		}

		Xnames2 <- names(Traindata)[!names(Traindata) %in% c(y,"Recession")]

		f <-as.formula(paste(y,"~",paste(Xnames2,collapse="+")))
		
		print(f)

		#linear regression
		lr<-lm(f, data=Traindata)
		
		rmse<-sqrt(mean(lr$residuals^2))
		r2<-1-sum(lr$residuals^2)/sum((Traindata[,y]-mean(Traindata[,y]))^2)
		r2adjust <- r2-(1-r2)*length(Xnames2)/(nrow(Traindata)-length(Xnames2)-1)
		#aic<-AIC(lr)
		#bic<-BIC(lr)
		df<-nrow(Traindata)-length(Xnames2)-1
		fvar <- var(lr$fitted.values)
		resfvar <- var(lr$fitted.values[Traindata$Recession==1])
		rvar <- var(lr$residuals)
		resrvar <- var(lr$residuals[Traindata$Recession==1])
		corr<-cor(lr$residuals,lr$fitted.values)
		rescorr<-cor(lr$residuals[Traindata$Recession==1],lr$fitted.values[Traindata$Recession==1])
		modeloutput[nrow(modeloutput)+1,] <- c(y, rmse, r2, r2adjust, df, fvar, resfvar, rvar, resrvar, corr,rescorr, "LM")
		write.csv(summary(lr)$coefficients,paste0("lr_fs_",y,".csv"))
		residuals[,y]<-c(rep(NA,nrow(residuals)-length(lr$residuals)),lr$residuals)
	},
		error = function(ex) {
	#			print("errors")
			fvar <- NA
			resfvar <- NA
			rvar <- NA
			resrvar <- NA
			corr <- NA
			rescorr <- NA
		}
	)

}

	write.csv(modeloutput,paste0("modeloutput3.csv"))
	write.csv(residuals,paste0("residuals3.csv"))

#loop to repair correlation matrix for non-positive definite
repairall <- function(C){

	tryCatch(
	{
		chol(C)
		sa<-1
		return (sa)
	},
		error = function(ex) {
		sa<-0
		return (sa)
		}
	)
}	

#repair correlation matrix for non-positive definite
repaircorr<-function(C){
	# compute eigenvectors/-values
	E   <- eigen(C, symmetric = TRUE)   
	V   <- E$vectors
	D   <- E$values

	# replace negative eigenvalues by 0.001
	D   <- pmax(D,0.001)

	# reconstruct correlation matrix
	BB  <- V %*% diag(D) %*% t(V)

	# rescale correlation matrix
	T   <- 1/sqrt(diag(BB))
	TT  <- outer(T,T)
	C   <- BB * TT
	return (C)
}

#correlation matrix for expansion periods
um<-"pairwise.complete.obs" #"pairwise.complete.obs" "complete.obs" "all.obs" "na.or.complete"
alldata <- read.csv("residuals1.csv", header=TRUE, sep=",", dec=".")
alldata$Recession <- rawdata$Recession[3:length(rawdata$Recession)] #lag = 2
cordata <- alldata[alldata$Recession == 0,]
cordata <- alldata[,!names(alldata) %in% c("Recession", "X")]
normalcorr <- cor(cordata, use = um, method = "pearson")
normalcorr[is.na(normalcorr)]<-0
normalcorr["aaadefault","aaadefault"]<-1
icount <-0
while (repairall(normalcorr)==0) {
	normalcorr <- repaircorr(normalcorr)
	icount <- icount+1
	#print(icount)
}
normalchol <- chol(normalcorr)

#Alternatively, we can reduce records with NA before calculating correlation matrix.
#Then there is no need to repair the correlation matrix
cordata <- cordata[complete.cases(cordata),]
normalcorr <- cor(cordata, use = um, method = "pearson")
normalcorr[is.na(normalcorr)]<-0
normalcorr["aaadefault","aaadefault"]<-1
icount <-0
while (repairall(normalcorr)==0) {
	normalcorr <- repaircorr(normalcorr)
	icount <- icount+1
	print(icount)
}
normalchol <- chol(normalcorr)

#correlation matrix for recession periods
cordata <- alldata[alldata$Recession == 1,]
cordata <- cordata[,!names(cordata) %in% c("Recession","X")]
recessioncorr <- cor(cordata, use = um, method = "pearson")
recessioncorr[is.na(recessioncorr)]<-0
recessioncorr["aaadefault","aaadefault"]<-1
#recessioncorr <-ifelse(abs(recessioncorr)>abs(normalcorr),recessioncorr,normalcorr)
icount <-0
while (repairall(recessioncorr)==0) {
	recessioncorr <- repaircorr(recessioncorr)
	icount <- icount+1
	print(icount)
}
recessionchol <- chol(recessioncorr)

colnames(normalcorr) <- Ynames
colnames(normalchol) <- Ynames
colnames(recessioncorr) <- Ynames
colnames(recessionchol) <- Ynames
write.csv(normalcorr,paste0("normalcorr.csv"),row.names = FALSE)
write.csv(normalchol,paste0("normalchol.csv"),row.names = FALSE)
write.csv(recessioncorr,paste0("recessioncorr.csv"),row.names = FALSE)
write.csv(recessionchol,paste0("recessionchol.csv"),row.names = FALSE)


################################################################
# Macroeconomic factor VAR #####################################
################################################################
library(vars)
rawdata <- read.csv("input/varinput.csv", header=TRUE, sep=",", dec=".")
Xnames <- c("R_","pi_c_","pi_i_","pi_d_","dy_","dc_","di_","dimp_","dex_","dE_","dS_","dw_","dy_star_","pi_star_","R_star_")
Traindata <- rawdata[,names(rawdata) %in% c(Xnames)]
Traindata <- Traindata[complete.cases(Traindata),]
	
var1 <- VAR(Traindata, p = 1, type = "const") #both
stab1 <- stability(var1, h = 0.15, dynamic = FALSE, rescale = TRUE) #type = c("OLS-CUSUM", "Rec-CUSUM", "Rec-MOSUM","OLS-MOSUM", "RE", "ME", "Score-CUSUM", "Score-MOSUM", "fluctuation"),
plot(stab1)

################################################################
# Fit data to Normal distribution ##############################
################################################################
rawdata <- read.csv("input/varinput.csv", header=TRUE, sep=",", dec=".")
Xnames <- c("R_","pi_c_","pi_i_","pi_d_","dy_","dc_","di_","dimp_","dex_","dE_","dS_","dw_","dy_star_","pi_star_","R_star_")

LL <- 0
for (i in Xnames) {
	Traindata <- rawdata[,names(rawdata) %in% i]
	Traindata <- as.numeric(Traindata)
	Traindata <- Traindata[!is.na(Traindata)]
	fit <- fitdistr(Traindata, densfun="normal")
	LL<-fit$loglik+LL
}

print(LL)

