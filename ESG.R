setwd("C:/temp/2017/dsge/fm")

################################################################
# Read Simulated Macroeconomic Factors by the DSGE Model #######
################################################################
macrofac <- read.csv("input/oo_endo_simul.csv",header = FALSE)
# Select the factors used for multifactor regression
cols <- c("R_","pi_c_","pi_i_","pi_d_","dy_","dc_","di_","dimp_","dex_","dE_")
macrofac <- macrofac[c(66,67,69:75,79),]
macrofac <- t(macrofac)
colnames(macrofac) <- cols

# multifactor model
mapping <- read.csv("input/mapping.csv")
mappingNames <- colnames(mapping)[2:(ncol(mapping)-6)]
normalchol <- read.csv("input/normalchol.csv")
recessionchol <- read.csv("input/recessionchol.csv")
inputmap <- read.csv("input/inputmap.csv", header=TRUE, sep=",", dec=".")
Xnames <- c("R_","pi_c_","pi_i_","pi_d_","dy_","dc_","di_","dimp_","dex_","dE_")
Ynames <- names(inputmap)[!names(inputmap) %in% c(Xnames,"Year","Quarter","Recession")]
histAR <- inputmap[,names(inputmap) %in% Ynames]
histAR <- histAR[(nrow(histAR)-1):nrow(histAR),]
histMF <- inputmap[,names(inputmap) %in% Xnames]
histMF <- histMF[(nrow(histMF)-1):nrow(histMF),]
histRecession <- inputmap[(nrow(inputmap)-1):nrow(inputmap),names(inputmap) %in% c("Recession")]
termmix <- read.csv("input/termmix.csv")
migration <- read.csv("input/migration.csv")

#Recession Logistic function
recession <- function(paras, vals){
	prob <- 1/(1+exp(-sum(paras * vals)))
	if(prob > 0.5) {
		return(1)
	}else{
		return(0)
	}
}

recession_paras <- c(-78.93143715, 58.41928887, -85.67123867, 1.257377168, -3.811050839,
					-158.0435208, 4.149726594, -0.926564319, 8.786617544, -5.035576855,
					-68.25666285, 50.83903724, 1.584409523, -74.9761853, 17.21169545, 86.52820878)
# Intercept, pi_c_, dy_, dc_, di_, dE_, pi_c_1, dy_1, dc_1, di_1, dE_1, pi_c_2, dy_2, dc_2, di_2, dE_2

set.seed(123)
# ESG: create one single scenario
esg <- function(fundmap,histMF,histAR,cholNormal,cholRecession,macrofac,period){
	sim_mf <- rbind(histMF, macrofac[1:period,])
	mappingX <- sim_mf
	sim_mf$recession <- 0
	sim_mf$recession[1:2] <- histRecession
	recessionX = sim_mf[, names(sim_mf) %in% c("pi_c_", "dy_", "dc_", "di_", "dE_")]
	lag <-2
	recessionNames <- c("pi_c_", "dy_", "dc_", "di_", "dE_")
	if (lag > 0) {
		for (i in c(1:lag)){
			for (varname in recessionNames){
				recessionX[(i+1):nrow(recessionX),paste0(varname,i)] <- recessionX[1:(nrow(recessionX)-i),varname]
				recessionX[1:i,paste0(varname,i)] <- NA
			}
		}
	}
	recessionX <- recessionX[(lag+1):nrow(recessionX),]
	recessionX <- cbind(Intercept = 1, recessionX)
	for (i in c(3:length(sim_mf))) {
		sim_mf$recession[i] <- recession(recession_paras,recessionX[i,])
	}
	
	sim_mf <- cbind(sim_mf,histAR)
	
	sdtnormal <- sqrt(fundmap[,names(fundmap) %in% c("tvar")])
	sdtrecession <- sqrt(fundmap[,names(fundmap) %in% c("trvar")])
	sdinormal <- sqrt(fundmap[,names(fundmap) %in% c("ivar")])
	sdirecession <- sqrt(fundmap[,names(fundmap) %in% c("irvar")])
	corrnormal <- fundmap[,names(fundmap) %in% c("tcorr")]
	corrrecession <- fundmap[,names(fundmap) %in% c("rcorr")]

	lag <-2
	if (lag > 0) {
		for (i in c(1:lag)){
			for (varname in Xnames){
				mappingX[(i+1):nrow(mappingX),paste0(varname,i)] <- mappingX[1:(nrow(mappingX)-i),varname]
				mappingX[1:i,paste0(varname,i)] <- NA
			}
		}
	}
	mappingX <- mappingX[, colnames(fundmap)[5:(ncol(fundmap)-6)]]

	for (i in c(1:period)){
		mappingX_i <- mappingX[i+2,]
		mappingX_i <- data.frame(c(1,0,0,mappingX_i))
		colnames(mappingX_i) <- mappingNames

		mappingX_i_m <- as.data.frame(lapply(mappingX_i, rep, length(sdinormal)))
		mappingX_i_m$x2 <- as.numeric(sim_mf[i,c(12:ncol(sim_mf))])
		mappingX_i_m$x1 <- as.numeric(sim_mf[i+1,c(12:ncol(sim_mf))])
		
		mappingfunc <- fundmap[,2:(ncol(fundmap)-6)]
		det_returns = rowSums(mappingfunc * mappingX_i_m)

		if (sim_mf$recession[i+2]==0){
			rnds <- data.matrix(normalchol) %*% (sdinormal * rnorm(length(sdinormal)))
			rnds <- (corrnormal * det_returns + sqrt(1-corrnormal*corrnormal)*rnds) * sdinormal
			rnds <- as.numeric(rnds)/sqrt(corrnormal * corrnormal * sdtnormal * sdtnormal+(1-corrnormal*corrnormal)*sdinormal*sdinormal)
			rnds[14] <- 0
		}else{
			rnds <- data.matrix(recessionchol) %*% (sdirecession * rnorm(length(sdirecession)))
			rnds <- (corrrecession * det_returns + sqrt(1-corrrecession*corrrecession)*rnds) * sdirecession
			rnds <- as.numeric(rnds)/sqrt(corrrecession * corrrecession * sdtrecession * sdtrecession+(1-corrrecession*corrrecession)*sdirecession*sdirecession)
			rnds[14] <- 0
		}
		sim_mf[i+2,c(12:ncol(sim_mf))] <- det_returns + as.numeric(rnds)
		print(i)
		print(det_returns)
		print(as.numeric(rnds))
	}
	return(sim_mf)
}
one_scn <- esg(mapping,histMF,histAR,cholNormal,cholRecession,macrofac,20)


#Curve fitting
curvefitting <- function(value,term,inter,extro,target){
# value: yields at terms defined in argument term
# term: a vector that contains the terms in the yield curve
# inter: interpolation method: 0: linear; 1: cubic spline
# extro: extropolation method: 0: constant forward rate; 1: follow inter method with target rate at term 100 years
# target: long term target rate for extropolation
	if (length(term)<2) {
		return(-1)
	} else {
		startT <- term[1]
		endT <- term[length(term)]
		secondEndT <- term[length(term)-1]
		fullterm <- seq(startT, endT, by=0.25)
		if (inter == 0){
			fullvalue <- approx(term, value, fullterm)$y
		} else {
			fullvalue <- spline(term, value, xout = fullterm)$y		
		}
		if (endT<100 && extro == 0) {
			fullterm2 <- seq((endT+0.25),100,by=0.25)
			cfr <- (((1+value[length(value)])^endT)/((1+value[length(value)-1])^secondEndT))^(1/(endT-secondEndT))-1
			cumfac <- (1+value[length(value)])^endT
			cumfacvec <- cumprod(c((1+value[length(value)])^endT, rep((1+cfr),length(fullterm2))))
			cumfacvec <- cumfacvec[2:length(cumfacvec)]
			fullvalue2 <- cumfacvec^(1/fullterm2)-1
			fullvalue <- c(fullvalue, fullvalue2)
			fullterm <- c(fullterm, fullterm2)
		}
		
		if (endT<100 && extro == 1) {
			fullterm2 <- seq((endT+0.25),100,by=0.25)
			fullterm <- c(fullterm, fullterm2)
			term <- c(term, 100)
			value <- c(value, target)
			if (inter==1) {
				fullvalue <- spline(term, value, xout = fullterm)$y
			} else {
				fullvalue <- approx(term, value, fullterm)$y		
			}
		}

		return(list(term = fullterm,rate = fullvalue))
	}
}

terms <- c(0.25,1,2,3,5,7,10,20,30)
values <- c(0.7,1,1.2,1.5,1.7,1.86,2,2.35,2.9)/100
plot(curvefitting(values,terms,0,0,0.04)$term,curvefitting(values,terms,0,0,0.04)$rate)
plot(curvefitting(values,terms,0,1,0.04)$term,curvefitting(values,terms,0,1,0.04)$rate)
plot(curvefitting(values,terms,1,1,0.04)$term,curvefitting(values,terms,1,1,0.04)$rate)

#bond pricer
bondPrice <- function(maturity,freq,coupon,redemption,value,term,inter,extro,target){
	xi <- curvefitting(value,term,inter,extro,target)$term
	yi <- curvefitting(value,term,inter,extro,target)$rate
	if(maturity>freq) {
		bondvalue <- redemption*(1+coupon*freq)/(1+approx(term, value, maturity)$y)^maturity
		for(i in seq(maturity-freq,0.001,by=-freq)) {
			bondvalue <- bondvalue + redemption * coupon * freq/(1+approx(term, value, i)$y)^i
		}
	}else if(maturity <=0){
		bondvalue <- 0
	}else{
		bondvalue <- redemption*(1+coupon*freq)/(1+approx(term, value, maturity, rule=2)$y)^maturity
	}
	
	return(bondvalue)
}
terms <- c(0.25,1,2,3,5,7,10,20,30)
values <- c(0.7,1,1.2,1.5,1.7,1.86,2,2.35,2.9)/100
bondPrice(20,0.25,0.00,10000,values,terms,0,1,0.04)

bondReturn <- function(sim_mf,mix,migration,period,rating,rebalance,bondfreq,inter,extro,target){
#rating: 1-Govt Bond; 2-AAA; 3-AA; 4-A; 5-BBB

	term <- c(0.25,1,2,3,5,7,10,20,30)
	initialVal <- 1000000
	migrationM <- migration/100
	recoveryM <- migrationM[,7]*100
	#rating=2
	#bondfreq=0.25
	MVM <- termmix[rating,] * initialVal
	value <- sim_mf[1,c(12:20)]
	if (rating>1){
		value = value + sim_mf[1,20+rating*2]
	}
	default <- as.numeric(c(0,sim_mf[1,c(25:28)]/100))

	xnew <- seq(0.25, 30, by=0.25)
	if (inter == 0){
		CouponM <- approx(term, value, xnew)$y
	} else {
		CouponM <- spline(term, value, xout = xnew)$y		
	}

	FVM <- MVM
	for (i in c(1:120)){
		FVM[i] <- FVM[i]/bondPrice(i*0.25,bondfreq,CouponM[i]/100,1,value/100,term,inter=inter,extro=extro,target=target)
	}
	nBSM <- rep(0.0,120)
	RedemptionM <- rep(0.0,120)
	tMVM <- MVM
	tFVM <- FVM
	tBSM <- nBSM
	tRDM <- RedemptionM
	tCRM <- CouponM
	cashRtn <- vector()
	priceRtn <- vector()

	for (i in c(3:nrow(sim_mf))){
		bcurvebase <- sim_mf[i,c(12:20)]
		bcurve <- rbind(bcurvebase,bcurvebase + sim_mf[i,21])
		bcurve <- rbind(bcurve,bcurvebase + sim_mf[i,22])
		bcurve <- rbind(bcurve,bcurvebase + sim_mf[i,23])
		bcurve <- rbind(bcurve,bcurvebase + sim_mf[i,24])
		bcurve <- rbind(bcurve,bcurvebase + 2.499197)

		nMVM <- rep(0.0,120)
		nFVM <- rep(0.0,120)
		nBSM <- rep(0.0,120)
		RedemptionM <- rep(0.0,120)

		totalBS <- 0.0
		xnew <- seq(0.25, 30, by=0.25)
		if (inter == 0){
			nCouponM <- approx(term, value, xnew)$y
		} else {
			nCouponM <- spline(term, value, xout = xnew)$y		
		}

		for (j in c(1:120)){
			for (k in c(1:6)){
				if (k==rating){
					nMVM[j] <- as.numeric(bondPrice(j/4.0-0.25,bondfreq,CouponM[j]/100,FVM[j],bcurve[rating,]/100,term,inter=inter,extro=extro,target=target)*migrationM[rating,rating]*(1-default[rating]*(1-recoveryM[rating])))
					nFVM[j] <- as.numeric(FVM[j]*migrationM[rating,rating]*(1-default[rating]*(1-recoveryM[rating])))
				}else{
					totalBS <- totalBS + as.numeric(bondPrice(j/4.0-0.25,bondfreq,CouponM[j]/100,FVM[j],bcurve[k,]/100,term,inter=inter,extro=extro,target=target)*migrationM[rating,k]*(1-default[rating]*(1-recoveryM[rating])))
				}
			}
			if (j==1){
				RedemptionM[j] = as.numeric(nMVM[j]*(1-default[rating]*(1-recoveryM[rating])))
			}
			if (((j %% as.integer(bondfreq*4))==1 && (bondfreq>0.25)) || (bondfreq==0.25)){
				RedemptionM[j] = RedemptionM[j] + as.numeric(FVM[j]*CouponM[j]/100*bondfreq*(1-default[rating]*(1-recoveryM[rating])))
			}
			print(nMVM[j])
		}
	
		totalBS = totalBS + sum(RedemptionM)
		if (rebalance==0){
			nBSM <- totalBS*termmix[rating,]
			nBSFM <- rep(0.0,120)		
			for (j in c(1:120)){
				nBSFM[j] <- as.numeric(nBSM[j]/bondPrice(j/4.0,bondfreq,nCouponM[j]/100,1,bcurve[rating,]/100,term,inter, extro,target))
			}
			for (j in c(2:120)){
				MVM[j-1] <- nMVM[j]
				FVM[j-1] <- nFVM[j]
			}
			MVM[120]<-0
			FVM[120]<-0
			for (j in c(2:120)){
				if ((FVM[j-1]+nBSFM[j-1]) != 0){
					CouponM[j-1]<-as.numeric((FVM[j-1]*CouponM[j]+nBSFM[j-1]*nCouponM[j-1])/(FVM[j-1]+nBSFM[j-1]))
				}else{
					CouponM[j-1]<-as.numeric(CouponM[j])
				}
			}
			CouponM[120]<-nCouponM[120]
			MVM = MVM + nBSM
			FVM = FVM + nBSFM
		}else if ((i %% as.integer(rebalance*4))==0){
			totalBS <- totalBS + sum(nMVM)
			nBSM <- totalBS*termmix[rating,]
			for (j in c(1:119)){
				nBSM[j] <- nBSM[j]-nMVM[j+1]
			}
			nBSFM <- rep(0.0,120)		
			for (j in c(1:120)){
				if (nBSM[j]<0){
					nBSFM[j] <- as.numeric(nBSM[j]/bondPrice(j/4.0,bondfreq,CouponM[j]/100,1,bcurve[rating,]/100,term,inter=inter,extro=extro,target=target))
				}else{
					nBSFM[j] <- as.numeric(nBSM[j]/bondPrice(j/4.0,bondfreq,nCouponM[j]/100,1,bcurve[rating,]/100,term,inter=inter,extro=extro,target=target))			
				}
			}
			for (j in c(2:120)){
				MVM[j-1] <- nMVM[j]
				FVM[j-1] <- nFVM[j]
			}
			MVM[120] <- 0
			FVM[120] <- 0
			for (j in c(2:120)){
				if (nBSM[j-1]>0 && (FVM[j-1]+nBSFM[j-1])!=0){
					CouponM[j-1] <- as.numeric((FVM[j-1]*CouponM[j]+nBSFM[j-1]*nCouponM[j-1])/(FVM[j-1]+nBSFM[j-1]))
				}else{
					CouponM[j-1] <- CouponM[j]
				}
			}
			
			CouponM[120] <- nCouponM[120]
			MVM <- MVM + nBSM
			FVM <- FVM + nBSFM
		}else{
			nBSM <- totalBS*termmix[rating,]
			nBSFM <- rep(0.0,120)		
			for (j in c(1:120)){
				nBSFM[j] <- as.numeric(nBSM[j]/bondPrice(j/4,bondfreq,nCouponM[j]/100,1,value/100,term,inter=inter,extro=extro,target=target))
			}
			for (j in c(2:120)){
				MVM[j-1] <- nMVM[j]
				FVM[j-1] <- nFVM[j]
			}
			MVM[120] <- 0
			FVM[120] <- 0
			for (j in c(2:120)){
				if ((FVM[j-1]+nBSFM[j-1]) != 0){
					CouponM[j-1] <- as.numeric((FVM[j-1]*CouponM[j]+nBSFM[j-1]*nCouponM[j-1])/(FVM[j-1]+nBSFM[j-1]))
				}else{
					CouponM[j-1] <- CouponM[j]
				}
			}
			CouponM[120] <- nCouponM[120]
			MVM <- MVM + nBSM
			FVM <- FVM + nBSFM
		}
		tMVM <- rbind(tMVM,MVM)
		tFVM <- rbind(tFVM,FVM)
		tBSM <- rbind(tBSM,nBSM)
		tRDM <- rbind(tRDM,RedemptionM)
		tCRM <- rbind(tCRM,CouponM)
		cashRtn = as.numeric(c(cashRtn,sum(RedemptionM)/initialVal))
		priceRtn = as.numeric(c(priceRtn,(sum(MVM)-sum(RedemptionM))/initialVal-1))
		default <- as.numeric(c(0,sim_mf[i+1,c(25:28)]/100))
		initialVal <- sum(MVM)
	}
	#print tMVM[range(0,5),:]
	#print tFVM[range(0,5),:]
	#print tBSM[range(0,5),:]
	#print tRDM[range(0,5),:]
	#print tCRM[range(0,5),:]
	return cashRtn,priceRtn

}
