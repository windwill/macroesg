setwd("C:/dsge/r5/")

################################################################
# Simulated Macroeconomic Factors based on the DSGE Model #######
################################################################
exo_vars <- c("eps_c",
"eps_i",
"eps_g",
"eps_r",
"eps_a",
"eps_z",
"eps_H",
"eps_pi_cbar",
"eps_x",
"eps_d",
"eps_mc",
"eps_mi",
"eps_z_tildestar",
"eps_mu_z",
"eps_Rstar",
"eps_pistar",
"eps_ystar",
"eps_i_k",
"eps_i_y",
"me_w",
"me_E",
"me_pi_d",
"me_pi_i",
"me_y",
"me_c",
"me_i",
"me_imp",
"me_ex",
"me_ystar")

endo_vars <- c("y",
"c",
"i",
"g",
"imp",
"ex",
"q",
"c_m",
"i_m",
"i_d",
"pi_cbar",
"pi_c",
"pi_i",
"pi_d",
"pi_mc",
"pi_mi",
"pi_x",
"mc_d",
"mc_mc",
"mc_mi",
"mc_x",
"r_k",
"w",
"P_k",
"mu_z",
"kbar",
"k",
"u",
"dS",
"x",
"R",
"a",
"psi_z",
"H",
"E",
"gamma_mcd",
"gamma_mid",
"gamma_xstar",
"gamma_cd",
"gamma_id",
"gamma_f",
"R_star",
"pi_star",
"y_star",
"e_ystar",
"e_pistar",
"e_c",
"e_i",
"e_a",
"e_z",
"e_H",
"lambda_x",
"lambda_d",
"lambda_mc",
"lambda_mi",
"z_tildestar",
"i_k",
"i_c",
"i_y",
"i_w",
"R_",
"pi_c_",
"pi_cbar_",
"pi_i_",
"pi_d_",
"dy_",
"dc_",
"di_",
"dimp_",
"dex_",
"dy_star_",
"pi_star_",
"R_star_",
"dE_",
"dS_",
"dw_",
"AUX_ENDO_LAG_25_1")


dr_orders <- c("imp",
"ex",
"q",
"c_m",
"i_m",
"i_d",
"pi_i",
"mc_d",
"mc_mc",
"mc_mi",
"u",
"H",
"gamma_cd",
"gamma_id",
"gamma_f",
"i_c",
"i_y",
"i_w",
"R_",
"pi_c_",
"pi_cbar_",
"pi_i_",
"pi_d_",
"dy_",
"dc_",
"di_",
"dimp_",
"dex_",
"dy_star_",
"pi_star_",
"R_star_",
"dE_",
"dS_",
"dw_",
"y",
"g",
"pi_cbar",
"pi_c",
"mc_x",
"kbar",
"k",
"x",
"R",
"a",
"gamma_mcd",
"gamma_mid",
"gamma_xstar",
"R_star",
"e_ystar",
"e_pistar",
"e_i",
"e_a",
"e_z",
"e_H",
"lambda_x",
"lambda_d",
"lambda_mc",
"lambda_mi",
"z_tildestar",
"AUX_ENDO_LAG_25_1",
"c",
"i",
"pi_d",
"pi_mc",
"pi_mi",
"pi_x",
"w",
"mu_z",
"E",
"y_star",
"e_c",
"i_k",
"r_k",
"P_k",
"dS",
"psi_z",
"pi_star")


#import the first order approximation of the calibrated DSGE model
dsge <- read.csv("input/dsge_foa.csv",header = TRUE)
#import the shocks in the calibrated DSGE model
dsge_shocks <- read.csv("input/dsge_shocks.csv",header = TRUE)

nperiod <- 120
nsims <- 1000
results <- matrix(0,nrow=nperiod*nsims,ncol=length(endo_vars))

steady_s <- dsge$y_s
A <- dsge[,4:(length(endo_vars)+3)]
B <- dsge[,(length(endo_vars)+4):ncol(dsge)]

set.seed(6)

counter <- 1
for (isim in c(1:nsims)){
	scn <- as.matrix(dsge$last_y)
	steady <- as.matrix(dsge$y_s)
	for (iperiod in c(1:nperiod)){
		shocks <- as.matrix(sapply(dsge_shocks$sigma, function(x) rnorm(1,0,x)))
		scn <- steady + as.matrix(A) %*% (scn-steady) + as.matrix(B) %*%  shocks
		results[counter,] <- scn
		counter <- counter + 1
	}
	if (isim %% 10 == 0){
		print(paste0("scenario ",isim," has been generated"))
	}
}

# Select the factors used for multifactor regression
colnames(results) <- dr_orders
keeps <- c("R_",
"pi_c_",
"pi_cbar_",
"pi_i_",
"pi_d_",
"dy_",
"dc_",
"di_",
"dimp_",
"dex_",
"dy_star_",
"pi_star_",
"R_star_",
"dE_",
"dS_",
"dw_"
)

results <- results[,colnames(results) %in% keeps]
results <- as.data.frame(results)
results$quarter <- rep(c(1:nperiod),nsims)

apply(results,2,mean)
apply(results,2,sd)

# multifactor model
mapping <- read.csv("input/mapping.csv")
mappingNames <- colnames(mapping)[2:(ncol(mapping)-6)]
normalchol <- read.csv("input/normalchol.csv")
recessionchol <- read.csv("input/recessionchol.csv")
inputmap <- read.csv("input/inputmap.csv", header=TRUE, sep=",", dec=".")
Xnames <- c("R_","pi_c_","pi_i_","pi_d_","dy_","dc_","di_","dimp_","dex_","dE_","dS_","dw_","dy_star_","pi_star_","R_star_")
Ynames <- names(inputmap)[!names(inputmap) %in% c(Xnames,"Year","Quarter","Recession","aaadefault","spmid","spmidd","spsmall","spsmalld",
												"ereit","ereitinc","mreitret","mreitinc","wage","aa1y","aa2y","aa3y","aa5y","aa7y","aa10y","aa20y","aa30y")]
histAR <- inputmap[,names(inputmap) %in% Ynames]
histAR <- histAR[(nrow(histAR)-1):nrow(histAR),]
histMF <- inputmap[,names(inputmap) %in% Xnames]
histMF <- histMF[(nrow(histMF)-1):nrow(histMF),]
histRecession <- inputmap[(nrow(inputmap)-1):nrow(inputmap),names(inputmap) %in% c("Recession")]

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
esg <- function(fundmap,histMF,histAR,cholNormal,cholRecession,macrofac,period, sim){
	macrofac <- macrofac[, colnames(macrofac) %in% c("R_", "pi_c_", "pi_i_", "pi_d_", "dy_", "dc_", "di_", "dimp_", "dex_", "dE_","dS_","dw_","dy_star_","pi_star_","R_star_")]
	sim_mf <- rbind(histMF, macrofac[(nperiod*(sim-1)+1):(nperiod*(sim-1)+period),])
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
	for (i in c(3:nrow(sim_mf))) {
		sim_mf$recession[i] <- recession(recession_paras,recessionX[i-2,])
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
		mappingX_i_m$x2 <- as.numeric(sim_mf[i,c(17:ncol(sim_mf))])
		mappingX_i_m$x1 <- as.numeric(sim_mf[i+1,c(17:ncol(sim_mf))])
		
		mappingfunc <- fundmap[,2:(ncol(fundmap)-6)]
		det_returns = rowSums(mappingfunc * mappingX_i_m)

		if (sim_mf$recession[i+2]==0){
			rnds <- t(data.matrix(normalchol)) %*% (sdinormal * rnorm(length(sdinormal)))
			rnds <- as.numeric(rnds)
		}else{
			rnds <- t(data.matrix(recessionchol)) %*% (sdirecession * rnorm(length(sdirecession)))
			rnds <- (corrrecession * det_returns + sqrt(1-corrrecession*corrrecession)*rnds) * sdirecession
			rnds <- as.numeric(rnds)/sqrt(corrrecession * corrrecession * sdtrecession * sdtrecession+(1-corrrecession*corrrecession)*sdirecession*sdirecession)
		}
		sim_mf[i+2,c(17:ncol(sim_mf))] <- det_returns + as.numeric(rnds)*1
	}
	return(sim_mf)
}

set.seed(6)
nscns <- 1000
zero_scn <- esg(mapping,histMF,histAR,normalchol,recessionchol,results,nperiod,1)

for (i in c(1:(nscns-1))){
	isim <- i %% 1000
	if (isim==0) {isim=1000}
	one_scn <- esg(mapping,histMF,histAR,normalchol,recessionchol,results,nperiod,isim)
	zero_scn <- rbind(zero_scn, one_scn)
	if (i %% 10 == 0){print(paste0(i," simulations are done."))}
}

zero_scn$quarter <- rep(c((-1):nperiod),nscns)
write.csv(zero_scn,"esg_dsge_1k.csv",row.names=FALSE)

apply(zero_scn,2,mean)
apply(zero_scn,2,sd)

percs <- c(0.005,0.01,0.05, 0.1,0.25,0.5,0.75,0.9,0.95,0.99,0.995)

perc <- function(scns, ynames, percs, periods){
	modeloutput <- data.frame(y=character(),
				 period = double(),
				 min_=double(),
				 perc_005 = double(), 
				 perc_01 = double(), 
				 perc_05 = double(), 
				 perc_1 = double(), 
				 perc_25 = double(), 
				 perc_50 = double(), 
				 perc_75 = double(), 
				 perc_90 = double(), 
				 perc_95 = double(), 
				 perc_99 = double(), 
				 perc_995 = double(), 
				 max_=double(),
                 stringsAsFactors=FALSE)
		
	counter = 1

	for (yname in ynames){
		for (p in c((-1):periods)){
			idata <- scns[which(scns$quarter == p),]
			idata <- idata[,yname]
			min_ <- min(idata)
			percs_val <- quantile(idata, percs)
			max_ <- max(idata)
			modeloutput[counter,] = c(yname, p, min_, percs_val, max_)
			counter = counter + 1
		}
	
	}
	
	return(modeloutput)

}

perc(zero_scn,c("dy_"),percs,20)