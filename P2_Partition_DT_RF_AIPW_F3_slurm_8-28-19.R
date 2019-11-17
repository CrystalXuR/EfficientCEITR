##---------------------------------------------------------------------------------------------##
##---- Simulation code for Optimal Individualized Regimen (ITR) on CE analysis       ----------##
##---- Dissertation topic of Yizhe Xu -- University of Utah, Biostatistics           ----------##
##---- Date: 7/31/19                                                                ----------##
options(scipen = 999)#setwd("Z:/Crystal/OptimalITRPartition")
library(rpart);library(ggplot2);library(gridExtra);library(survival);
library(data.table);library(party);library(dplyr);library(geepack);library(randomForest)
source("./ITR Partition functions for CE analysis_7-17-19.R")
source("./ITR Partition functions for CE simulation_7-17-19_F3TC.R")
params = read.csv("./parameters.csv")

if(length(args <- commandArgs(T))>0){
  stopifnot(length(args)==1)
  case.id <- as.integer(args[[1]])
  message("running for parameter set ", case.id)
}

# Basic settings
N = 1000    # sample size
iter = 100  # replication
tau = 20

costmodel = c("gamma", "lognormal")      # cost models 
classTools = c("DT","RF")               # Classification tools 

#### failure time params
h0 = 0.1                                # baseline hazard
betaF = c(0, 0.7, 0.6, 0.3, 0.3, 0.3)

#### cost params
kappa = 2.5                             # shape param
gammaC = c(0.1, 0.1)
betaC = c(0, 0.3, 0.3, 0.1, 0.2, 0.1)
#case.id=1;z=1
gammaF = as.numeric(params[case.id,1:2])
lambda = as.numeric(params[case.id,3])
cup = as.numeric(params[case.id,4])


output = NULL  
EMO = NULL     
CCR = NULL  
for (z in 1:iter){
  
  set.seed(999+z^2)
  #####----------------------------------------------------------------------------###  
  #####------------ Data Simulation: Survival time & Cumulative Cost --------------###
  #####----------------------------------------------------------------------------###
  # baseline covariates
  x1 = rnorm(N,1,2)
  x2 = rnorm(N,1,2)
  x3 = rnorm(N)
  x4 = rnorm(N)
  x5 = rnorm(N)
  X0 = cbind(x1,x2,x3,x4,x5)

  # simulate treatment 
  u0 = 0.2+0.2*x3
  p0 = 1/(1+exp(-u0))
  A = rbinom(N,1,p0);table(A)
  
  ###------ Suvival time ~ Exponential (PH) model -------###
  surv.prob = runif(N)
  SurvT.trt = sapply(1:N, trt_i=rep(1,N), betaF=betaF, gammaF=gammaF, h0=h0, surT.exp)   # counterfactual failure time under A=1 
  SurvT.ctl = sapply(1:N, trt_i=rep(0,N), betaF=betaF, gammaF=gammaF, h0=h0, surT.exp)   # counterfactual failure time under A=0
  SurvT.A = sapply(1:N, trt_i=A, betaF=betaF, gammaF=gammaF, h0=h0, surT.exp)  
  
  ###------ Counterfactual Restricted Suvival times -------###
  RSurvT.trt = pmin(SurvT.trt,rep(tau,N))
  RSurvT.ctl = pmin(SurvT.ctl,rep(tau,N))
  RSurvT.A = pmin(SurvT.A,rep(tau,N))
  #summary(RSurvT.trt); summary(RSurvT.ctl);summary(RSurvT.A)

  ###------ Censoring time: ~ exp(0,0.3) for low  ~ exp(0,0.8) for high -------###
  cenT=rexp(N,cup);  
  cen.obs = as.numeric(RSurvT.A<=cenT);table(cen.obs) # cen.obs=0:censored; % censored = 13.5%
  FUT.obs = pmin(RSurvT.A,cenT)
  
  # Counterfactual event indicator
  event.trt = as.numeric(SurvT.trt<=tau)
  event.ctl = as.numeric(SurvT.ctl<=tau)
  event.A = as.numeric(SurvT.A<=tau)
  
  ###---- Create monthly censoring indicator for generating cost data -------###
  Mtau = tau*12; M.cenT = cenT*12
  M.SurvT.trt = floor(RSurvT.trt*12); M.SurvT.ctl = floor(RSurvT.ctl*12); M.SurvT.A = floor(RSurvT.A*12)
  
  # create censoring indicator up to Mtau months
  I.delta.trt = I.delta.ctl = I.delta.A = matrix(NA, N, Mtau)
  for (i in 1:N){
      I.delta.trt[i,] = c(rep(1,M.SurvT.trt[i]), rep(NA,(Mtau-M.SurvT.trt[i])))
      I.delta.ctl[i,] = c(rep(1,M.SurvT.ctl[i]), rep(NA,(Mtau-M.SurvT.ctl[i])))
      I.delta.A[i,] = c(rep(1,M.SurvT.A[i]), rep(NA,(Mtau-M.SurvT.A[i])))
    }
  
  I.delta.obs = matrix(NA, N, Mtau)
  for (i in 1:N){
    if(M.SurvT.A[i]<=M.cenT[i]){
      I.delta.obs[i,] = c(rep(1,M.SurvT.A[i]), rep(NA,(Mtau-M.SurvT.A[i])))
    }
    else if(M.SurvT.A[i]>M.cenT[i]){
      I.delta.obs[i,] = c(rep(1,floor(M.cenT[i])),0, rep(NA,(Mtau-1-floor(M.cenT[i]))))
    }}
  
  ###------ Cumulative cost ~ Gamma dist. ------### 
  # cumulative cost up to each month
  cumC.trt = 100*as.matrix(t(sapply(1:N, trt_i=rep(1,N), k=kappa, betaC=betaC, gammaC=gammaC, I.delta=I.delta.trt, event=event.trt, tau=tau, cumC.gam))) 
  cumC.ctl = 100*as.matrix(t(sapply(1:N, trt_i=rep(0,N), k=kappa, betaC=betaC, gammaC=gammaC, I.delta=I.delta.ctl, event=event.ctl, tau=tau, cumC.gam)))
  cumC.A = 100*as.matrix(t(sapply(1:N, trt_i=A, k=kappa, betaC=betaC, gammaC=gammaC, I.delta=I.delta.A, event=event.A, tau=tau, cumC.gam)))     
  cumC.obs = 100*as.matrix(t(sapply(1:N, trt_i=A, k=kappa, betaC=betaC, gammaC=gammaC, I.delta=I.delta.obs, event=cen.obs, tau=tau, cumC.gam)))

  # total cost
  ind.trt <- !is.na(cumC.trt)
  cost_tot.trt = tapply(cumC.trt[ind.trt], row(cumC.trt)[ind.trt], tail, 1)         
  
  ind.ctl <- !is.na(cumC.ctl)
  cost_tot.ctl = tapply(cumC.ctl[ind.ctl], row(cumC.ctl)[ind.ctl], tail, 1)        
  
  ind.A <- !is.na(cumC.A)
  cost_tot.A = tapply(cumC.A[ind.A], row(cumC.A)[ind.A], tail, 1)    # the last non-missing cost value in cumCost
  
  ind <- !is.na(cumC.obs)
  cost_tot.obs = tapply(cumC.obs[ind], row(cumC.obs)[ind], tail, 1)
  #summary(cost_tot.trt);summary(cost_tot.ctl);summary(cost_tot.A);summary(cost_tot.obs)
  
  ### Observed: change cost data & censoring indicator to long format 
  cumC.obs.test = data.frame(1:N, cumC.obs)
  colnames(cumC.obs.test) = c("subjectid", paste(rep("cost",(Mtau+1)),as.character(0:Mtau),sep = "_"))
  cost_long = reshape(cumC.obs.test, varying = c(2:(Mtau+2)), idvar = "subjectid",direction="long", sep = "_", timevar = "month")
  Cost_long.obs = cost_long[order(cost_long$subjectid, cost_long$month),];names(Cost_long.obs)[3] = "cost.obs"
  
  cenInd.obs = data.frame(1:N, rep(1,N),I.delta.obs);names(cenInd.obs)=c("subjectid",paste(rep("delta.obs",(Mtau+1)),as.character(0:Mtau),sep = "_"))
  cenInd_long.obs = reshape(cenInd.obs, varying = c(2:(Mtau+2)), idvar = "subjectid",direction="long", sep = "_", timevar = "month")
  cenInd_long.obs = cenInd_long.obs[order(cenInd_long.obs$subjectid, cenInd_long.obs$month),]
  
  # merge survival time data 
  long.obs = merge(x=Cost_long.obs,y=cenInd_long.obs,by=c("subjectid","month"),all = TRUE);dim(long.obs)
  long.obs = long.obs[is.na(long.obs$delta.obs)==FALSE, ]
  long.obs = long.obs[order(long.obs$subjectid, long.obs$month),];head(long.obs)
  
  ###------ CE outcome ------###
  Y.trt = lambda*RSurvT.trt-cost_tot.trt 
  Y.ctl = lambda*RSurvT.ctl-cost_tot.ctl
  Y.A = lambda*RSurvT.A-cost_tot.A
  Y.obs = lambda*FUT.obs-cost_tot.obs   
  #summary(Y.trt);summary(Y.ctl);summary(Y.A);summary(Y.obs)
  
  ###------ Numerically compute optimal g.opt -------###
  g.opt = as.numeric((lambda*(RSurvT.trt-RSurvT.ctl)-(cost_tot.trt-cost_tot.ctl))>0)#;table(g.opt,A)
  
  ###------ Compare Mean treatment effect under each g.opt ------###
  dat.outcomes = data.frame(cbind(Y.trt, Y.ctl, g.opt))
  Y.g.opt = (g.opt*Y.trt+(1-g.opt)*Y.ctl)/10000
  
  #-----------------------------------------------------#
  # Different options for analyzing the simulation data #
  est_ps = EstPS(A=A,Xs=x3)
  mainvars = c("x2","x3","x4","x5");intvars = c("x1")
  ####------------------------------------------- Analysis Code ------------------------------------------------####
  ####------ Several extra steps for analyzing panel data compaing to one-patient-per-row data frame -----------####
  #### 1. compute cost increment for each month interval                                         ---------------####
  #### 2. calculate censoring weight for each interval                                          ----------------####
  #### 3. estimate the monthly cost increment under treated and untreated using a posited model ----------------####
  ####----------------------------------------------------------------------------------------------------------####
  
  ObsData = data.frame(1:N, A, x1, x2, x3, x4, x5, cenT, FUT.obs, cost_tot.obs, cen.obs, est_ps)
  names(ObsData) = c("subjectid","A","x1", "x2", "x3", "x4", "x5","cenT","FUT.obs", "cost_tot.obs", "cen.obs", "est_ps0","est_ps1")
  ObsDatt = merge(x=ObsData,y=long.obs,by="subjectid",all = TRUE)
  ObsDat = ObsDatt[order(ObsDatt$subjectid, ObsDatt$month),]
  
  # Monthly cost increments within each subject
  ObsDat = data.table(ObsDat)
  ObsDat2 = ObsDat[,lagCost.obs:=c(NA, cost.obs[-.N]), by=subjectid]
  ObsDat2[is.na(ObsDat2$lagCost.obs)==TRUE & is.na(ObsDat2$cost.obs)==FALSE,]$lagCost.obs=0
  ObsDat2$MoCost.obs = ObsDat2$cost.obs-ObsDat2$lagCost.obs#; summary(ObsDat2$MoCost.obs)
  
  # Estimated Monthyly censoring weights: we do this way instead of GEE as we want the weights to be 1 when no one is censored
  Mo.cen.w = list(NA,(Mtau+1))
  for (i in 0:Mtau){
    DatMo = data.frame(ObsDat2[ObsDat2$month==i,])        # extract subjects that contributed to month i
    if(length(levels(as.factor(DatMo$delta.obs)))==1|dim(DatMo)[1]==0){
      Mo.cen.w[[(i+1)]] = rep(1,dim(DatMo)[1])          # if everyone in month i is uncensored (=1), then the censoring weight is 1 for all subjects (no case with all subjects are censored)
    }else{
      cen.wm = glm(delta.obs ~ x3+x4+x5+A*x1+A*x2, family = "binomial",data = DatMo)
      Mo.cen.w[[(i+1)]] = predict(cen.wm, type = "response")
    }
  }
  Mo.cen.w.long = unlist(Mo.cen.w)#;summary(Mo.cen.w.long)
  ObsDat2b = ObsDat2[order(ObsDat2$month,ObsDat2$subjectid),]
  ObsDat2b$Mo.cen.w = Mo.cen.w.long
  ObsDat2 = ObsDat2b[order(ObsDat2b$subjectid,ObsDat2b$month),]
  
  # Estimated censoring weights 
  cen.wm = glm(cen.obs ~ x3+x4+x5+A*x1+A*x2, family = "binomial", data = ObsData)
  ObsData$cen.w = predict(cen.wm, type = "response")
  
  # Estimated restricted survival time T using coxph model
  Xs=dplyr::select(ObsData,mainvars);Ms=dplyr::select(ObsData,intvars)
  REG.ExpT = Reg.mu.Exp(ST=FUT.obs, A=A, Xs=Xs, Ms=Ms, event = cen.obs, data=ObsData)
  muTpre = REG.ExpT$mus.ExpT
  muT = matrix(NA,N,2)
  for (j in 1:2){
    muT[,j] = pmin(rep(tau,N),muTpre[,j])
  }
  
  # Create a time interval variable for cost to distinguish treatment effects in intital and death interval is different from others
  ObsDat2$time = ifelse(ObsDat2$month==0,1,NA)
  lastdat = do.call("rbind",by(ObsDat2, INDICES=ObsDat2$subjectid, FUN=function(DF) DF[which.max(DF$month), ]))
  lastdat$time = ifelse(lastdat$cen.obs==1 & lastdat$month!=0,3,
                        ifelse(lastdat$month!=0,2,4))
  # table(lastdat$time) there are 151 subjects only had 0 month follow-up, so the cost in that month=inital+death
  #  2   3   4 
  # 164 685 151 
  lastdat2 = dplyr::select(lastdat,c("subjectid","month","time"))
  ObsDat3 = merge(x=ObsDat2,y=lastdat2,by=c("subjectid","month"),all=T)
  ObsDat3$time = ifelse(is.na(ObsDat3$time.y)==F,ObsDat3$time.y,
                        ifelse(is.na(ObsDat3$time.y)==T & is.na(ObsDat3$time.x)==F,ObsDat3$time.x,
                               ifelse(is.na(ObsDat3$time.x)==T & is.na(ObsDat3$time.y)==T,2,NA)))#; ObsDat3$time = as.factor(ObsDat3$time)
  #   1    2    3    4 
  # 849  6655  685  151 
  
  # Estimated monthly costs using gamma or lognormal model
  Xss=dplyr::select(ObsDat3,mainvars);Mss=dplyr::select(ObsDat3,c(intvars,"time"))
  for (CM in costmodel){
    if(CM=="gamma"){
      Mo.REGC.gam = Reg.mu.gam(CC=ObsDat3$MoCost.obs, A=ObsDat3$A, Xs=Xss, Ms=Mss, data=ObsDat3)
      Mo.muC = Mo.REGC.gam$mus.gam
      
      REGC.gam = Reg.mu.gam(CC=ObsData$cost_tot.obs, A=ObsData$A, Xs=Xs, Ms=Ms,data=ObsData)
      muC = REGC.gam$mus.gam
    }else if(CM=="lognormal"){
      Mo.REGC.norm = Reg.mu.norm(CC=ObsDat3$MoCost.obs, A=ObsDat3$A, Xs=Xss, Ms=Mss,data=ObsDat3)
      Mo.muC = Mo.REGC.norm$mus.norm
      
      REGC.norm = Reg.mu.norm(CC=ObsData$cost_tot.obs, A=ObsData$A, Xs=Xs, Ms=Ms,data=ObsData)
      muC = REGC.norm$mus.norm
    }
    
    ObsDat4 = data.frame(ObsDat3,Mo.muC);names(ObsDat4)[(dim(ObsDat4)[2]-1):dim(ObsDat4)[2]] = c("Mo.muC.ctl", "Mo.muC.trt")#;head(ObsDat4)
    ObsDatLong = ObsDat4[order(ObsDat4$subjectid,ObsDat4$month),]#;head(ObsDatLong) # 15513 19

    
    ###------- Non-Partitioned Reg-naive -------###
    Reg.naive_deltaY_NP = lambda*(muT[,2]-muT[,1])-(muC[,2]-muC[,1])
    Reg.naive_ITR_NP = as.numeric(Reg.naive_deltaY_NP>0)#;table(Reg.naive_ITR_NP)
    
    ###------- NON-Partitioned AIPW weigths -------###
    ITR.aipw_NP=ITR_AIPW_NP(Q = ObsData$FUT.obs, C = ObsData$cost_tot.obs, A = ObsData$A, est_ps= cbind(ObsData$est_ps0,ObsData$est_ps1),
                        cen = ObsData$cen.obs, est_cen = ObsData$cen.w, lambda = lambda, muQ.reg=muT, muC.reg=muC)
    aipw.Y_NP=ITR.aipw_NP$Contrast.Y
    aipw.ITRCE_NP=ITR.aipw_NP$ITR_regCE
    
    ###------- Partitioned IPW weigths -------###
    ITR.ipw_P=ITR_IPW_P(Q = ObsData$FUT.obs, C = ObsDatLong$MoCost.obs, A = ObsData$A, est_ps= cbind(ObsData$est_ps0,ObsData$est_ps1),
                         cen = ObsData$cen.obs, est_cen = ObsData$cen.w, lambda = lambda,long_ps=cbind(ObsDatLong$est_ps0,ObsDatLong$est_ps1), 
                         Mocen=ObsDatLong$delta.obs, Moest_cen=ObsDatLong$Mo.cen.w, muQ.reg=muT, datalong = ObsDatLong, data = ObsData)
    ipw.Y_P=ITR.ipw_P$Contrast.Y
    ipw.ITRCE_P=ITR.ipw_P$ITR_regCE
    
    ###------- Partitioned AIPW weigths -------###
    ITR.aipw_P=ITR_AIPW_P(Q = ObsData$FUT.obs, C = ObsDatLong$MoCost.obs, A = ObsData$A, est_ps= cbind(ObsData$est_ps0,ObsData$est_ps1),
                          cen = ObsData$cen.obs, est_cen = ObsData$cen.w, long_ps = cbind(ObsDatLong$est_ps0,ObsDatLong$est_ps1), 
                          Mocen=ObsDatLong$delta.obs, Moest_cen=ObsDatLong$Mo.cen.w, lambda = lambda, muQ.reg = muT, muC.reg = Mo.muC,
                          datalong = ObsDatLong, data = ObsData)
    aipw.Y_P=ITR.aipw_P$Contrast.Y
    aipw.ITRCE_P=ITR.aipw_P$ITR_regCE
    
    
    ###------- Weighted Classification using CART/PARTY ----------###
    full.dat = data.frame(Y.trt, Y.ctl, Y.g.opt, X0, A, g.opt,ipw.ITRCE_P,ipw.Y_P,
                          aipw.ITRCE_P,aipw.Y_P,Reg.naive_ITR_NP,Reg.naive_deltaY_NP)
    comp.dat = full.dat[complete.cases(full.dat)==T,]#;dim(comp.dat)  # 1000 23
    
    full.dat_NP = data.frame(X0,g.opt,aipw.ITRCE_NP,aipw.Y_NP)
    comp.dat_NP = full.dat_NP[complete.cases(full.dat_NP)==TRUE,]#;dim(comp.dat_NP)  # 1000 23
    
    comp.dat_NP$aipw.ITRCE_NP = as.factor(comp.dat_NP$aipw.ITRCE_NP)
    comp.dat$ipw.ITRCE_P = as.factor(comp.dat$ipw.ITRCE_P)
    comp.dat$aipw.ITRCE_P = as.factor(comp.dat$aipw.ITRCE_P)
    comp.dat$true.ww  = abs(comp.dat$Y.trt - comp.dat$Y.ctl)
    comp.dat$Reg_deltaY_NP = abs(comp.dat$Reg.naive_deltaY_NP)
    comp.dat$g.optt = as.factor(comp.dat$g.opt)
    
    
    for (ct in classTools){
      if(ct=="DT"){

        # NP -- TRUE
        NP_true=rpart(g.opt ~ x1+x2+x3+x4+x5, weights=comp.dat$true.ww, data=comp.dat, method="class",
                      control = rpart.control(minsplit = 10, cp = 0.01,xval = 10))
        comp.dat$g.true_NP=as.numeric(predict(NP_true, newdata = comp.dat,type="class"))-1
        
        # NP -- AIPW
        NP_aipw=rpart(aipw.ITRCE_NP ~ x1+x2+x3+x4+x5, weights=comp.dat_NP$aipw.Y_NP, data=comp.dat_NP, method="class",
                    control = rpart.control(minsplit = 10, cp = 0.01,xval = 10))
        comp.dat$g.aipw_NP=as.numeric(predict(NP_aipw, newdata = comp.dat,type="class"))-1
        
        # P -- IPW 
        P_ipw=rpart(ipw.ITRCE_P ~ x1+x2+x3+x4+x5, weights=comp.dat$ipw.Y_P, data=comp.dat, method="class",
                      control = rpart.control(minsplit = 10, cp = 0.01,xval = 10))
        comp.dat$g.ipw_P=as.numeric(predict(P_ipw, newdata = comp.dat,type="class"))-1
        
        # P -- AIPW 
        P_aipw=rpart(aipw.ITRCE_P ~ x1+x2+x3+x4+x5, weights=comp.dat$aipw.Y_P, data=comp.dat, method="class",
                     control = rpart.control(minsplit = 10, cp = 0.01,xval = 10))
        comp.dat$g.aipw_P=as.numeric(predict(P_aipw,newdata = comp.dat, type="class"))-1
        
      }else if(ct=="RF"){
        
        # NP -- TRUE
        NP_true=cforest(g.optt ~ x1+x2+x3+x4+x5, weights=comp.dat$true.ww, data=comp.dat,
                        control = cforest_unbiased(ntree = 200, mtry=sqrt(5)))
        comp.dat$g.true_NP=as.numeric(as.character(predict(NP_true, newdata = comp.dat, type="response")))

        # NP -- AIPW
        NP_aipw=cforest(aipw.ITRCE_NP ~ x1+x2+x3+x4+x5, weights=comp.dat_NP$aipw.Y_NP, data=comp.dat_NP,
                        control = cforest_unbiased(ntree = 200, mtry=sqrt(5)))
        comp.dat$g.aipw_NP=as.numeric(as.character(predict(NP_aipw, newdata = comp.dat, type="response")))
        
        # P -- IPW
        P_ipw=cforest(ipw.ITRCE_P ~ x1+x2+x3+x4+x5, weights=comp.dat$ipw.Y_P, data=comp.dat,
                        control = cforest_unbiased(ntree = 200, mtry=sqrt(5)))
        comp.dat$g.ipw_P=as.numeric(as.character(predict(P_ipw, newdata = comp.dat, type="response")))
        
        # P -- AIPW
        P_aipw=cforest(aipw.ITRCE_P ~ x1+x2+x3+x4+x5, weights=comp.dat$aipw.Y_P, data=comp.dat,
                       control = cforest_unbiased(ntree = 200, mtry=sqrt(5)))
        comp.dat$g.aipw_P=as.numeric(as.character(predict(P_aipw,newdata = comp.dat, type="response")))
      }
      
      # compare predicted optimal treatments
      GS = comp.dat$g.opt;TP = sum(GS==1)
      
      # Estimate counterfactual means under g.hat
      Ntest = dim(comp.dat)[1]
      
      # Gold standard: True optimal outcome for Yi
      E0.y = (comp.dat$Y.trt*comp.dat$g.opt+(1-comp.dat$g.opt)*comp.dat$Y.ctl)/10000
      
      # Predicted g and plugged-in outcome  
      # NP_Reg-naive
      E1.y = (comp.dat$Y.trt*comp.dat$Reg.naive_ITR_NP+(1-comp.dat$Reg.naive_ITR_NP)*comp.dat$Y.ctl)/10000
      
      # # NP_aipw
      E2.y = (comp.dat$Y.trt*comp.dat$g.aipw_NP+(1-comp.dat$g.aipw_NP)*comp.dat$Y.ctl)/10000
      
      # P_IPW 
      E3.y = (comp.dat$Y.trt*comp.dat$g.ipw_P+(1-comp.dat$g.ipw_P)*comp.dat$Y.ctl)/10000
      
      # P_AIPW 
      E4.y = (comp.dat$Y.trt*comp.dat$g.aipw_P+(1-comp.dat$g.aipw_P)*comp.dat$Y.ctl)/10000
      
      # Output
      CCR = rbind(CCR,c(z,CM,ct,mean(comp.dat$g.true_NP==GS),mean(comp.dat$Reg.naive_ITR_NP==GS),
                        mean(comp.dat$g.aipw_NP==GS),
                        mean(comp.dat$g.ipw_P==GS),mean(comp.dat$g.aipw_P==GS)))
      EMO = rbind(EMO,c(z,CM,ct,mean(E0.y),mean(E1.y),mean(E2.y),mean(E3.y),mean(E4.y)))
      output = rbind(output,data.frame(cbind(comp.dat,E0.y,E1.y,E2.y,E3.y,E4.y)))            
      print(paste("caseID=",case.id," iter=",z,", costModel=",CM,", classTools=",ct, " at ", Sys.time())) 
    } 
  }
}
colnames(CCR)=c("iter","costModel","classTools","C_True_NP","C_RegNaive_NP","C_AIPW_NP","C_IPW_P","C_AIPW_P");CCR
colnames(EMO)=c("iter","costModel","classTools","Truth","E_RegNaive_NP","E_AIPW_NP","E_IPW_P","E_AIPW_P")
output.long = do.call("rbind", output)#;dim(output.long)
## saving results
write.csv(CCR,file=paste("./Results_slurm_F3/HTE_",gammaF[1],"_",gammaF[2],"WTP_",lambda,"cenrate_",cup,"CCR_F3.csv",sep=""),row.names = F)
write.csv(EMO,file=paste("./Results_slurm_F3/HTE_",gammaF[1],"_",gammaF[2],"WTP_",lambda,"cenrate_",cup,"EMO_F3.csv",sep=""),row.names = F)
write.csv(output.long,file=paste("./Results_slurm_F3/HTE_",gammaF[1],"_",gammaF[2],"WTP_",lambda,"cenrate_",cup,"output_F3.csv",sep=""),row.names = F)

# org.ipw.Y_P = comp.dat$ipw.Y_P*(2*as.numeric(as.character(comp.dat$ipw.ITRCE_P))-1)
# org.aipw.Y_P = comp.dat$aipw.Y_P*(2*as.numeric(as.character(comp.dat$aipw.ITRCE_P))-1)
# org.aipw.Y_NP = comp.dat_NP$aipw.Y_NP*(2*as.numeric(as.character(comp.dat_NP$aipw.ITRCE_NP))-1)
# table(comp.dat$g.optt);table(comp.dat$Reg.naive_ITR_NP);table(comp.dat$ipw.ITRCE_P);table(comp.dat_NP$aipw.ITRCE_NP);table(comp.dat$aipw.ITRCE_P)
# par(mfrow=c(2,3))
# hist(Y.trt-Y.ctl);hist(Reg.naive_deltaY_NP);hist(org.ipw.Y_P);hist(org.aipw.Y_NP);hist(org.aipw.Y_P)
