##---------------------------------------------------------------------------------------------##
##---- Simulation code for Optimal Individualized Regimen (ITR) on CE analysis       ----------##
##---- Dissertation topic of Yizhe Xu -- University of Utah, Biostatistics           ----------##
##---- Date: 2/14/20                                                                 ----------##
options(scipen = 999)
library(survival);library(dplyr);library(data.table);library(cubature);library(rpart);library(party);library(stats);library(splines)
#setwd("C:/Users/cryst/Box Sync/A_Crystal's Thesis/2 - OptimalITRPartition")
source("./PartitionAnalysisFunctions.R")
source("./PartitionSimFunctionsTC.R")
params <- read.csv("./parameters.csv")

if(length(args <- commandArgs(T))>0){
  stopifnot(length(args)==1)
  case.id <- as.integer(args[[1]])
  message("running for parameter set ", case.id)
}#case.id=1;i=1

# Basic settings
N    <- 1000     
iter <- 500      
tau  <- 20
gammaF <- c(as.numeric(params[case.id,1:2]))
gammaC <- rep(as.numeric(params[case.id,3]),2)
lambda <- as.numeric(params[case.id,4])
cup    <- as.numeric(params[case.id,5])
ct     <- as.character(params[case.id,6])

out <- NULL  
EMO <- NULL     
CCR <- NULL  
for (i in 1:iter){
  tryCatch({
  #####----------------------------------------------------------------------------###  
  #####------------ Data Simulation: Survival time & Cumulative Cost --------------###
  #####----------------------------------------------------------------------------###
  set.seed(999+i^2)
  
  # baseline covariates
  x1 <- rnorm(N,1,2)
  x2 <- rnorm(N,1,2)
  x3 <- rnorm(N)
  x4 <- rnorm(N)
  x5 <- rnorm(N)
  X0 <- data.frame(x1,x2,x3,x4,x5)
  
  # Treatment 
  u0 <- 0.5*x1+0.5*x2+0.9*x3
  p0 <- 1/(1+exp(-u0))
  A  <- rbinom(N,1,p0)#;table(A)
  
  # Suvival time ~ Exponential (PH) model 
  surv.prob <- runif(N)
  SurvT.trt <- sapply(1:N, trt=rep(1,N), gammaF=gammaF, St=surv.prob, X0=X0, surT.Coxexp)     
  SurvT.ctl <- sapply(1:N, trt=rep(0,N), gammaF=gammaF, St=surv.prob, X0=X0, surT.Coxexp)    
  SurvT.obs <- sapply(1:N, trt=A,        gammaF=gammaF, St=surv.prob, X0=X0, surT.Coxexp)  
  #summary(SurvT.trt);summary(SurvT.ctl);summary(SurvT.obs)
  
  # Censoring time: ~ cox exponential PH:  h0=0.01 for low (20%), h0=0.05 for high (50%) 
  betaCen <- c(0.3, 0.3, 0.1, 0.1, 0.1)
  cenT    <- -log(surv.prob)/(cup*exp(as.matrix(X0)%*%betaCen)) 
  FUT.obs <- pmin(SurvT.obs,cenT)
  cen.obs <- as.numeric(SurvT.obs<=cenT)#;table(cen.obs)   # cen.obs=0:censored
  
  # Restricted observed suvival times 
  RFUT.obs <- pmin(FUT.obs,rep(tau,N))
  
  # Cumulative cost ~ Gamma dist. 
  M.CumC.trt <- floor(pmin(SurvT.trt,rep(tau,N))*2)      # accured up to this many 6-months (bounded by tau)
  M.CumC.ctl <- floor(pmin(SurvT.ctl,rep(tau,N))*2)
  M.FUT.obs  <- floor(RFUT.obs*2)
  M.SurvT.obs <- floor(SurvT.obs*2)
  
  # Indicator for whether or not count dying cost
  dy.trt <- ifelse(SurvT.trt<=tau,1,NA)                    # No loss-to-follow-up in counterfactuals 
  dy.ctl <- ifelse(SurvT.ctl<=tau,1,NA)                    # cannot observe death if it happens after tau
  dy.obs <- ifelse(SurvT.obs<=pmin(cenT,rep(tau,N)),1,NA)  # cannot observe death if it happens after cenT or tau
  
  # Indicator for accounting for the costs from which 6-months & 6-monthly event indicator 
  Mtau <- tau*2; M.cenT <- floor(cenT*2)
  I.Mcost.trt<-I.Mcost.ctl<-I.Mcost.obs<-matrix(NA, N, Mtau)
  for (j in 1:Mtau){
    I.Mcost.trt[,j] <- ifelse(rep(j,N)<=M.CumC.trt,1,NA)
    I.Mcost.ctl[,j] <- ifelse(rep(j,N)<=M.CumC.ctl,1,NA)
    I.Mcost.obs[,j] <- ifelse(rep(j,N)<=pmin(M.SurvT.obs,M.cenT),1,NA)  # cost cannot be counted after death or censoring
  }
  
  # Accumulated cost up to each 6-month
  cumC.trt <- 1000*as.matrix(t(sapply(1:N, trt=rep(1,N), gammaC=gammaC, I.delta=I.Mcost.trt, death=dy.trt, X0=X0, cumC.gam))) 
  cumC.ctl <- 1000*as.matrix(t(sapply(1:N, trt=rep(0,N), gammaC=gammaC, I.delta=I.Mcost.ctl, death=dy.ctl, X0=X0, cumC.gam)))
  cumC.obs <- 1000*as.matrix(t(sapply(1:N, trt=A,        gammaC=gammaC, I.delta=I.Mcost.obs, death=dy.obs, X0=X0, cumC.gam)))

  # Combining indicators of inital cost, 6-monthly cost, and death cost 
  I.cost.trt <- cbind(rep(1,N),I.Mcost.trt,dy.trt) 
  I.cost.ctl <- cbind(rep(1,N),I.Mcost.ctl,dy.ctl)
  I.cost.obs <- cbind(rep(1,N),I.Mcost.obs,dy.obs)
  
  # Total costs
  cost_tot.trt <- cumC.trt[,(Mtau+2)]    # if a subject dies or is censored before 20 years, the very last accumulated cost will be carried forward   
  cost_tot.ctl <- cumC.ctl[,(Mtau+2)]       
  cost_tot.obs <- cumC.obs[,(Mtau+2)]
  #summary(cost_tot.trt);summary(cost_tot.ctl);summary(cost_tot.obs)
  
  # NMBs & Truth
  Y.trt <- lambda*pmin(SurvT.trt,rep(tau,N))-cost_tot.trt 
  Y.ctl <- lambda*pmin(SurvT.ctl,rep(tau,N))-cost_tot.ctl
  Y.obs <- lambda*RFUT.obs-cost_tot.obs
  g.opt <- as.numeric(Y.trt-Y.ctl>0)
  Y.opt <- (g.opt*Y.trt+(1-g.opt)*Y.ctl)/10000
  
  # 6-monthly event indicator: if a subject is censored, the event indicator=0 for the rest follow-up;
  #                            if a subject had an event, the event indicator=1 till the time interval he/she dies, then NA for the rest follow-up.
  MyMdelta=I.Mcost.obs
  for (j in 1:N){
    if(cen.obs[j]==0){
      MyMdelta[j,][which(is.na(MyMdelta[j,]))]<-0 
    }
  }
  
  # Observed cost data & censoring indicator in long format 
  cumC.obs2 = cumC.obs*I.cost.obs   # 1000 242
  cumC.obs.test <- data.frame(1:N, cumC.obs2);colnames(cumC.obs.test) <- c("subjectid", paste(rep("cost",(Mtau+2)),as.character(0:(Mtau+1)),sep = "_"))
  cost_long <- reshape(cumC.obs.test, varying = c(2:(Mtau+3)), idvar = "subjectid",direction = "long", sep = "_", timevar = "month")
  Cost_long.obs <- cost_long[order(cost_long$subjectid, cost_long$month),];names(Cost_long.obs)[3] <- "cost.obs"
  
  cenInd.obs <- data.frame(1:N, rep(1,N),MyMdelta,dy.obs);names(cenInd.obs) <- c("subjectid",paste(rep("delta.obs",(Mtau+2)),as.character(0:(Mtau+1)),sep = "_"))
  cenInd_long.obs <- reshape(cenInd.obs, varying = c(2:(Mtau+3)), idvar = "subjectid",direction="long", sep = "_", timevar = "month")
  cenInd_long.obs <- cenInd_long.obs[order(cenInd_long.obs$subjectid, cenInd_long.obs$month),]
  
  long.obs <- merge(x=Cost_long.obs,y=cenInd_long.obs,by=c("subjectid","month"),all = TRUE) # dim(long.obs)
  long.obs <- long.obs[is.na(long.obs$delta.obs)==F, ]            # delete intervals after subjects had events  
  long.obs <- long.obs[order(long.obs$subjectid, long.obs$month),]
  
  
  ####------------------------------------------- Analysis Code ------------------------------------------------####
  ####------ Several extra steps for analyzing panel data compaing to one-patient-per-row data frame -----------####
  #### 1. compute cost increment for each month interval                                         ---------------####
  #### 2. calculate censoring weight for each interval                                          ----------------####
  #### 3. estimate the monthly cost increment under treated and untreated using a posited model ----------------####
  ####----------------------------------------------------------------------------------------------------------####
  
  # Simulation model 2: correct PS model, misspecified outcome model
  est_ps   <- EstPS(A=A,Xs=cbind(x1,x2,x3))
  mainvars <- c("x1","x2","x3","x4","x5");intvars <- c("x1")
  
  # Combine long data with baseline info.
  ObsData <- data.frame(1:N,X0,A,FUT.obs,RFUT.obs,cost_tot.obs,Y.obs,cen.obs,est_ps)
  names(ObsData)[c(1,13:14)] <- c("subjectid","est_ps0","est_ps1")
  ObsDatt <- merge(x=ObsData,y=long.obs,by="subjectid",all = TRUE)
  ObsDat  <- ObsDatt[order(ObsDatt$subjectid, ObsDatt$month),]
  
  # Monthly cost increments within each subject
  ObsDat  <- data.table(ObsDat)
  ObsDat2 <- ObsDat[,lagCost.obs:=c(NA, cost.obs[-.N]), by=subjectid]
  ObsDat2[is.na(ObsDat2$lagCost.obs)==TRUE & is.na(ObsDat2$cost.obs)==FALSE,]$lagCost.obs <- 0
  ObsDat2$MoCost.obs <- ObsDat2$cost.obs-ObsDat2$lagCost.obs
  
  # Estimated Monthyly censoring weights: pooled logistic regression
  ObsDatNonZero <- ObsDat2[ObsDat2$month!=0,]  # assumed no one is censored in month 0
  require(stats)
  cubsp <- data.frame(ns(ObsDatNonZero$month,knots=quantile(ObsDatNonZero$month,prob=c(0.25,0.5,0.75))))  # create natural cubic splines of month
  names(cubsp) <- c("Month1","Month2","Month3","Month4")
  ObsDatNonZero <- data.frame(ObsDatNonZero,cubsp)
  fitcenM <- glm(delta.obs ~ Month1+Month2+Month3+Month4+x1+x2+x3+x4+x5+A+A:x1+A:x2, family = "binomial", data = ObsDatNonZero)
  ObsDatNonZero$Mo.cen.w <- predict(fitcenM, type = "response");summary(ObsDatNonZero$Mo.cen.w)
  ObsDatNonZero <- ObsDatNonZero[, c("subjectid","month","Mo.cen.w")]
  ObsDat2 <- merge(x=ObsDat2,y=ObsDatNonZero, by=c("subjectid","month"),all=T)
  ObsDat2[is.na(ObsDat2$Mo.cen.w)==T,]$Mo.cen.w <- 1  # everyone gets censoring weights 1 in month 0
  
  # Estimated censoring weights 
  ObsData$complete <- ifelse((ObsData$cen.obs==1 | ObsData$FUT.obs>=tau),1,0); table(ObsData$complete)
  cen.wm <- glm(complete ~ x1+x2+x3+x4+x5+A+A:x1+A:x2, family = "binomial", data = ObsData)
  ObsData$cen.w <- predict(cen.wm, type = "response")
  
  Xs   <- dplyr::select(ObsData,mainvars);Ms <- dplyr::select(ObsData,intvars)
  REGT <- Reg.mu.surv(ST=FUT.obs, event = cen.obs, A=A, Xs=Xs, Ms=Ms, tau=tau, data=ObsData)
  muT  <- REGT$mus.survRT

  # Cost types: inital, ongoing, death
  ObsDat2$ctype <- ifelse(ObsDat2$month==0,1,
                        ifelse(ObsDat2$month==41,3,2))

  # Estimated monthly costs using gamma or lognormal model
  Xss <- dplyr::select(ObsDat2,c(mainvars,"ctype"));Mss <- dplyr::select(ObsDat2,intvars)
  M.RegC <- Reg.mu.gam(CC=ObsDat2$MoCost.obs, A=ObsDat2$A, Xs=Xss, Ms=Mss, data=ObsDat2)  # include A*cost type term when estimating monthly cost 
  Mo.muC <- M.RegC$mus.gam
    
  RegC <- Reg.mu.gam(CC=ObsData$cost_tot.obs, A=ObsData$A, Xs=Xs, Ms=Ms,data=ObsData)
  muC  <- RegC$mus.gam
  
  ObsDat3 <- data.frame(ObsDat2,Mo.muC);names(ObsDat3)[(dim(ObsDat3)[2]-1):dim(ObsDat3)[2]] <- c("Mo.muC.ctl", "Mo.muC.trt")#;head(ObsDat3)
  ObsDatLong <- ObsDat3[order(ObsDat3$subjectid,ObsDat3$month),]#;head(ObsDatLong) 
  
  # Non-Partitioned Reg-naive 
  Reg.naive_deltaY_NP <- lambda*(muT[,2]-muT[,1])-(muC[,2]-muC[,1])
  Reg.naive_ITR_NP    <- as.numeric(Reg.naive_deltaY_NP>0)   # table(Reg.naive_ITR_NP)
  
  # Non-Partitioned AIPW  
  ITR.aipw_NP <- ITR_AIPW_NP(Q = ObsData$RFUT.obs, C = ObsData$cost_tot.obs, A = ObsData$A, 
                             est_ps= cbind(ObsData$est_ps0,ObsData$est_ps1), cen = ObsData$complete,
                             est_cen = ObsData$cen.w, lambda = lambda, muQ.reg=muT, muC.reg=muC)
  aipw.Y_NP <- ITR.aipw_NP$Contrast.Y
  aipw.ITRCE_NP <- ITR.aipw_NP$ITR_regCE
  
  # Partitioned IPW  
  ITR.ipw_P <- ITR_IPW_P(Q = ObsData$RFUT.obs, C = ObsDatLong$MoCost.obs, A = ObsData$A, 
                         est_ps= cbind(ObsData$est_ps0,ObsData$est_ps1), cen = ObsData$complete, est_cen = ObsData$cen.w,
                         lambda = lambda, long_ps=cbind(ObsDatLong$est_ps0,ObsDatLong$est_ps1), Mocen=ObsDatLong$delta.obs, 
                         Moest_cen=ObsDatLong$Mo.cen.w, muQ.reg=muT, datalong = ObsDatLong, data = ObsData)
  ipw.Y_P <- ITR.ipw_P$Contrast.Y
  ipw.ITRCE_P <- ITR.ipw_P$ITR_regCE
  
  # Partitioned AIPW  
  ITR.aipw_P <- ITR_AIPW_P(Q = ObsData$RFUT.obs, C = ObsDatLong$MoCost.obs, A = ObsData$A,
                           est_ps= cbind(ObsData$est_ps0,ObsData$est_ps1), cen = ObsData$complete, est_cen = ObsData$cen.w,
                           lambda = lambda, long_ps = cbind(ObsDatLong$est_ps0,ObsDatLong$est_ps1), Mocen=ObsDatLong$delta.obs, 
                           Moest_cen=ObsDatLong$Mo.cen.w, muQ.reg = muT, muC.reg = Mo.muC, datalong = ObsDatLong, data = ObsData)
  aipw.Y_P <- ITR.aipw_P$Contrast.Y
  aipw.ITRCE_P <- ITR.aipw_P$ITR_regCE
  
  # Weighted Classification using CART/PARTY 
  comp.dat <- data.frame(X0, A, g.opt, Y.opt, Y.trt, Y.ctl, Reg.naive_deltaY_NP, Reg.naive_ITR_NP,
                         aipw.Y_NP, aipw.ITRCE_NP, ipw.ITRCE_P, ipw.Y_P, aipw.ITRCE_P, aipw.Y_P)
  Rdata.NP <- data.frame(X0, g.opt, aipw.ITRCE_NP, aipw.Y_NP)
  Cdata.NP <- Rdata.NP[complete.cases(Rdata.NP)==TRUE,]
  
  comp.dat$g.optt <- as.factor(comp.dat$g.opt)
  comp.dat$Reg_deltaY_NP <- abs(comp.dat$Reg.naive_deltaY_NP)
  comp.dat$ipw.ITRCE_P   <- as.factor(comp.dat$ipw.ITRCE_P)
  comp.dat$aipw.ITRCE_P  <- as.factor(comp.dat$aipw.ITRCE_P)
  Cdata.NP$aipw.ITRCE_NP <- as.factor(Cdata.NP$aipw.ITRCE_NP)
  
  if(ct=="DT"){

    # NP -- AIPW
    NP_aipw <- rpart(aipw.ITRCE_NP ~ x1+x2+x3+x4+x5, weights=Cdata.NP$aipw.Y_NP, data=Cdata.NP, method="class",
                     control = rpart.control(xval = 10))
    comp.dat$g.aipw_NP <- as.numeric(predict(NP_aipw, newdata = comp.dat,type="class"))-1
    
    # P -- IPW
    P_ipw <- rpart(ipw.ITRCE_P ~ x1+x2+x3+x4+x5, weights=comp.dat$ipw.Y_P, data=comp.dat, method="class",
                   control = rpart.control(xval = 10))
    comp.dat$g.ipw_P <- as.numeric(predict(P_ipw, newdata = comp.dat,type="class"))-1

    # P -- AIPW
    P_aipw <- rpart(aipw.ITRCE_P ~ x1+x2+x3+x4+x5, weights=comp.dat$aipw.Y_P, data=comp.dat, method="class",
                    control = rpart.control(xval = 10))
    comp.dat$g.aipw_P <- as.numeric(predict(P_aipw,newdata = comp.dat, type="class"))-1
  }else if(ct=="RF"){
    
    # cforest (party)
    # NP -- AIPW
    NP_aipw <- cforest(aipw.ITRCE_NP ~ x1+x2+x3+x4+x5, weights=Cdata.NP$aipw.Y_NP, data=Cdata.NP,
                       control = cforest_control(ntree = 50, mtry=2, maxdepth=5))
    comp.dat$g.aipw_NP <- as.numeric(as.character(predict(NP_aipw, newdata = comp.dat, type="response")))

    # P -- IPW
    P_ipw <- cforest(ipw.ITRCE_P ~ x1+x2+x3+x4+x5, weights=comp.dat$ipw.Y_P, data=comp.dat,
                     control = cforest_control(ntree = 50, mtry=2, maxdepth=5))
    comp.dat$g.ipw_P <- as.numeric(as.character(predict(P_ipw, newdata = comp.dat, type="response")))

    # P -- AIPW
    P_aipw <- cforest(aipw.ITRCE_P ~ x1+x2+x3+x4+x5, weights=comp.dat$aipw.Y_P, data=comp.dat,
                      control = cforest_control(ntree = 50, mtry=2, maxdepth=5))
    comp.dat$g.aipw_P <- as.numeric(as.character(predict(P_aipw,newdata = comp.dat, type="response")))
  }
  
  # EMO
  GS   <- comp.dat$g.opt
  
  #Under GS
  E0.y <- (comp.dat$Y.trt*comp.dat$g.opt+(1-comp.dat$g.opt)*comp.dat$Y.ctl)/10000

  # NP_Reg-naive
  E1.y <- (comp.dat$Y.trt*comp.dat$Reg.naive_ITR_NP+(1-comp.dat$Reg.naive_ITR_NP)*comp.dat$Y.ctl)/10000

  # NP_aipw
  E2.y <- (comp.dat$Y.trt*comp.dat$g.aipw_NP+(1-comp.dat$g.aipw_NP)*comp.dat$Y.ctl)/10000

  # P_IPW
  E3.y <- (comp.dat$Y.trt*comp.dat$g.ipw_P+(1-comp.dat$g.ipw_P)*comp.dat$Y.ctl)/10000

  # P_AIPW
  E4.y <- (comp.dat$Y.trt*comp.dat$g.aipw_P+(1-comp.dat$g.aipw_P)*comp.dat$Y.ctl)/10000

  # Output
  CCR <- rbind(CCR,c(i,ct,mean(comp.dat$Reg.naive_ITR_NP==GS),mean(comp.dat$g.aipw_NP==GS),
                          mean(comp.dat$g.ipw_P==GS),mean(comp.dat$g.aipw_P==GS)))
  EMO <- rbind(EMO,c(i,ct,mean(E0.y),mean(E1.y),mean(E2.y),mean(E3.y),mean(E4.y)))
  out <- rbind(out,data.frame(cbind(comp.dat,E0.y,E1.y,E2.y,E3.y,E4.y)))
  print(paste("caseID=",case.id,", iter=",i,", HTE=",gammaF[1],", WTP=",lambda,", CUP=",cup,", ct=",ct, ", at ", Sys.time(),sep=""))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}#;CCR
#c(mean(as.numeric(CCR[,3])),mean(as.numeric(CCR[,4])),mean(as.numeric(CCR[,5])),mean(as.numeric(CCR[,6])))
names(CCR)  <- c("iter","classTools","C_RegNaive_NP","C_AIPW_NP","C_IPW_P","C_AIPW_P")
names(EMO)  <- c("iter","classTools","Truth","E_RegNaive_NP","E_AIPW_NP","E_IPW_P","E_AIPW_P")
output.long <- do.call("rbind", out)

## saving results
write.csv(CCR,file=paste("./ResultsTC/HTE_",gammaF[1],"_",gammaF[2],"WTP_",lambda,"cenrate_",cup,"ct",ct,"CCR_F3.csv",sep=""),row.names = F)
write.csv(EMO,file=paste("./ResultsTC/HTE_",gammaF[1],"_",gammaF[2],"WTP_",lambda,"cenrate_",cup,"ct",ct,"EMO_F3.csv",sep=""),row.names = F)
write.csv(output.long,file=paste("./ResultsTC/HTE_",gammaF[1],"_",gammaF[2],"WTP_",lambda,"cenrate_",cup,"ct",ct,"output_F3.csv",sep=""),row.names = F)



