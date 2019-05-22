library(rpart)
library(rpart.plot)
library(FNN)
library(neuralnet)

rawdata <- read.csv("C:/dsge/r/inputmap.csv", header=TRUE, sep=",", dec=".")
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


lag <-1
residuals <- matrix(,nrow=nrow(rawdata)-lag,ncol=length(Ynames))
colnames(residuals) <- Ynames
for (y in Ynames) {
	Traindata <- rawdata[,names(rawdata) %in% c(y,Xnames,"Recession")]
	
	Xnames1 <- names(Traindata)[!names(Traindata) %in% c("Recession",y)]
	if ((y == "spginfrad") || (y == "psped")) {
		Xnames1 <- Xnames1[!Xnames1 %in% c(y)]
	}
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
	write.csv(summary(lr)$coefficients,paste0("C:/dsge/fm/lr_",y,".csv"))
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

	ann <- neuralnet(f, data=data.matrix(Traindata), hidden=c(5), linear.output=TRUE, stepmax = 100000, threshold=0.01, act.fct = "tanh", likelihood = TRUE, lifesign ="full", lifesign.step = 1000)
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

# predict recession
rawdata <- read.csv("C:/2015 Research/Pension Scenarios/inputmap.csv", header=TRUE, sep=",", dec=".")
Xnames <- c("Unemploy","gpdinv","pconsump","gdpgr")
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
	write.csv(summary(lr)$coefficients,paste0("C:/Data/data/lr_",Ynames,".csv"))

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
	write.csv(glmr$coefficients,paste0("C:/Data/data/glmr_",Ynames,".csv"))

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

}

	write.csv(modeloutput,paste0("C:/Data/data/modeloutput1.csv"))


#Let's try variable selection	
rawdata <- read.csv("C:/2015 Research/Pension Scenarios/inputmap.csv", header=TRUE, sep=",", dec=".")
Xnames <- c("Unemploy","gpdinv","pconsump","gdpgr","m3tb","aa10y","tb10y","cpi")#,"export","import","niinv"
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
	if ((y == "spginfrad") || (y == "psped")) {
		Xnames1 <- Xnames1[!Xnames1 %in% c(y)]
	}
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
		write.csv(summary(lr)$coefficients,paste0("C:/Data/data/vs/lr_",y,".csv"))
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

	write.csv(modeloutput,paste0("C:/Data/data/vs/modeloutput2.csv"))
	write.csv(residuals,paste0("C:/Data/data/vs/residuals2.csv"))

#correlation matrix
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

	# replace negative eigenvalues by zero
	D   <- pmax(D,0)

	# reconstruct correlation matrix
	BB  <- V %*% diag(D) %*% t(V)

	# rescale correlation matrix
	T   <- 1/sqrt(diag(BB))
	TT  <- outer(T,T)
	C   <- BB * TT
	return (C)
}

#removes <- c("aaaboa","aaboa","aboa","bbbboa","aaadefault","spmidd","spmid","spsmalld","spsmall","spdivd","spdiv","spginfrad","spginfra","psped","pspe")
keeps <- c("sp500","spmid","spsmall","spdiv","ge","spginfra","pspe","commod","oil","gold","recession") #
#removes <- c("aaadefault","spmidd","spmid","spsmalld","spsmall","spdivd","spdiv","spginfrad","spginfra","psped","pspe")
#removes <- c("spginfrad","spginfra","psped","pspe")
um<-"pairwise.complete.obs" #"pairwise.complete.obs" "complete.obs" "all.obs" "na.or.complete"
alldata <- read.csv("C:/Data/data/residuals.csv", header=TRUE, sep=",", dec=".")
#alldata <- alldata[,!names(alldata) %in% removes]
alldata <- alldata[,names(alldata) %in% keeps]
cordata <- alldata[,!names(alldata) %in% c("recession")]
normalcorr <- cor(cordata, use = um, method = "pearson")
icount <-0
while (repairall(normalcorr)==0) {
	normalcorr <- repaircorr(normalcorr)
	icount <- icount+1
	print(icount)
}
normalchol <- chol(normalcorr)
#normalcorr[is.na(normalcorr)]<-0
#normalcorr["aaadefault","aaadefault"]<-1

cordata <- alldata[alldata$recession == 1,]
cordata <- cordata[,!names(cordata) %in% c("recession")]
recessioncorr <- cor(cordata, use = um, method = "pearson")
icount <-0
while (repairall(recessioncorr)==0) {
	recessioncorr <- repaircorr(recessioncorr)
	icount <- icount+1
	print(icount)
}
recessionchol <- chol(recessioncorr)
#recessioncorr[is.na(recessioncorr)]<-0
#recessioncorr["aaadefault","aaadefault"]<-1


#normalcorr<-repaircorr(normalcorr)
#recessioncorr<-repaircorr(recessioncorr)
write.csv(normalcorr,paste0("C:/Data/data/normalcorr.csv"))
write.csv(normalchol,paste0("C:/Data/data/normalchol.csv"))
write.csv(recessioncorr,paste0("C:/Data/data/recessioncorr.csv"))
write.csv(recessionchol,paste0("C:/Data/data/recessionchol.csv"))

#check normality of residuals
ksoutput <- data.frame(y=character(),
				 ksstat=double(),
                 ksp=double(),
                 stringsAsFactors=FALSE)
checkdata <- read.csv("C:/Data/data/residuals.csv", header=TRUE, sep=",", dec=".")
checkdata <- checkdata[,!names(checkdata) %in% c("recession","aaadefault")]
for (i in names(checkdata)){
	x <- checkdata[,i]
	x <- x[!is.na(x)]
	mu <- mean(checkdata[,i],na.rm=TRUE)
	std <- sd(checkdata[,i],na.rm=TRUE)
	y<-rnorm(500,mu,std)
	s <- ks.test(x,y)
	ksoutput[nrow(ksoutput)+1,] <- c(i, s$statistic, s$p.value)
}
write.csv(ksoutput,paste0("C:/Data/data/ksoutput.csv"))

#generate correlated random variable
gencrv <- function(nchol,rchol,lrecession,varnames){
	mreturn <- matrix(NA,nrow=nrow(lrecession),ncol=length(varnames))
	colnames(mreturn) <- varnames
	nv <- length(varnames)
	cvarnames <- colnames(nchol)
	nc <- length(cvarnames)
	mreturn[,!colnames(mreturn) %in% cvarnames] <- rnorm((nv-nc)*nrow(lrecession))
	mreturn <- cbind(mreturn,lrecession)
	for (j in c(1:nrow(mreturn))){
		if(mreturn[j,"Recession"] == 0) {
			mreturn[j,colnames(mreturn) %in% cvarnames] <- (t(nchol) %*% rnorm(nc))
		} else {
			mreturn[j,colnames(mreturn) %in% cvarnames] <- (t(rchol) %*% rnorm(nc))
		}
	}
	return (mreturn)
}

#Fundamental risk factor VAR
library(vars)
rawdata <- read.csv("inputmap.csv", header=TRUE, sep=",", dec=".")
Xnames <- c("R_","pi_c_","pi_i_","pi_d_","dy_","dc_","di_","dimp_","dex_","dE_")

Traindata <- rawdata[,names(rawdata) %in% c(Xnames)]
	
var1 <- VAR(Traindata, p = 1, type = "const") #both
stab1 <- stability(var1, h = 0.15, dynamic = FALSE, rescale = TRUE) #type = c("OLS-CUSUM", "Rec-CUSUM", "Rec-MOSUM","OLS-MOSUM", "RE", "ME", "Score-CUSUM", "Score-MOSUM", "fluctuation"),
plot(stab1)
nr <- length(Xnames) + 1 #2
varoutput <- matrix(NA,nrow=nr, ncol=length(Xnames))
colnames(varoutput) <- Xnames
for (i in c(Xnames)) {
	varoutput[,i] <- var1$varresult[i][[1]]$coefficients
}

rownames(varoutput) <- names(var1$varresult[i][[1]]$coefficients)
write.csv(t(varoutput),"varoutput.csv")
write.csv(summary(var1)$corres,"var1corres.csv")

nchol <- chol(summary(var1)$corres)
mreturn <- matrix(NA,nrow=202*100,ncol=length(Xnames))
for (j in c(1:nrow(mreturn))){
	mreturn[j,] <- (t(nchol) %*% rnorm(length(Xnames)))
}
colnames(mreturn) <- Xnames
write.csv(mreturn,paste0("cvar1output.csv"))

#Solve stable means
tvaroutput <- t(varoutput)
A<-tvaroutput[,1:length(Xnames)]
B<-tvaroutput[,length(Xnames)+1]
A<- -A
for (i in c(1:nrow(A))){
	A[i,i] <- A[i,i]+1
}
stablemeans <- solve(A,B)