#####--------------------------------------------------------------------------------#####
#####----- Functions Created for the Paper Entitled "Estimating the Optimal    ------#####
#####----- Individualized Treatment Rule from A Cost-Effectiveness Perspective ------#####
#####--------------------------------------------------------------------------------#####
library(survival);library(party);library(dplyr);library(cubature);require(stats);library(splines)
#----------------------------------------------------------------------------------------------#
# Inputs: covariate matrix - X; 
#         covariates in main terms of survival outcome models- Xs;
#         covariates in interaction terms of outcome models- Ms; 
#         covariates in main terms of cost outcome models- XsC;
#         observed survival time - FT; 
#         event indicator - event; 
#         observed interval cost - C;
#         indicator of NOT being event within an interval - Mnoncen;
#         study length of interest - tau;
#         willingness-to-pay parameter - lambda;
#         weight function - wf;
#         conditional random forest hyperparameters - ntree, mtry, maxdepth;
#         name of subject ID variable - id;
#         name of treatment variable - MA;
#         name of interval index variable - interval;
#         data set in long format with interval cost data 
#----------------------------------------------------------------------------------------------#

EfficientCEITR <- function(X, Xs, Ms, XsC, FT, event, C, Mnoncen, tau, lambda, wf="ITR_AIPW",
                           ntree, mtry, maxdepth, id, MA, interval, data){
  
  # Function: Estimate Propensity Score 
  .EstPS <- function(A,X){
    
    if(sort(unique(A))[1]!=0 | sort(unique(A))[2]!=1) 
      stop("Treatment levels are not 0 and 1")
    
    dat <- data.frame(A,X)
    fit <- glm(A ~., family="binomial", data=dat)
    e   <- predict(fit,dat,type = "response")
    PS  <- data.frame(1-e,e);names(PS) = c("PrA0","PrA1")
    return(PS) 
  }
  
  # Function: Estimate Restricted Survival Time 
  .Reg.mu.surv <- function(FT, event, A, Xs, Ms, tau, data){
    
    N <- dim(data)[1]; Xs <- as.matrix(Xs); Ms <- as.matrix(Ms)
    
    if(sort(unique(A))[1]!=0 | sort(unique(A))[2]!=1) 
      stop("Treatment levels are not 0 and 1")
    
    fitsurvT <- survreg(Surv(time = FT, event = event) ~ Xs+A+A:Ms, dist="lognormal", data=data)
    
    # estimate location parameter ui (on the log scale) for each individual & under each treatment arm 
    log.survT <- matrix(NA,N,2)
    for(j in 1:2){
      newdat <- data; newdat$A <- rep(sort(unique(A))[j],N)
      log.survT[,j] <- log(predict(fitsurvT,newdata=newdat,type = "response"))
    }
    
    mulog0 <- log.survT[,1]; mulog1 <- log.survT[,2]; sdlog <- summary(fitsurvT)$scale
    mus.survRT <- matrix(NA,N,2)
    for (i in 1:N){
      # survival function for each treatment arm
      S_t0 <- function(t){plnorm(t,meanlog=mulog0[i],sdlog=sdlog,lower.tail = FALSE)}
      S_t1 <- function(t){plnorm(t,meanlog=mulog1[i],sdlog=sdlog,lower.tail = FALSE)}
      # integrate from 0 to tau
      mus.survRT[i,1] <- cubintegrate(S_t0,0,tau,method = "cuhre",relTol=1e-04,absTol=1e-10)$integral
      mus.survRT[i,2] <- cubintegrate(S_t1,0,tau,method = "cuhre",relTol=1e-04,absTol=1e-10)$integral
    }
    
    outsurvRT <- list(mus.survRT, fitsurvT); names(outsurvRT) <- c("mus.survRT","SurvModel")
    return(mus.survRT)
  }
  
  # Function: Estimate cumulative cost (long data)
  .Reg.mu.gam <- function(C, A, XsC, Ms, data){
    
    N <- dim(data)[1]; XsC <- as.matrix(XsC); Ms <- as.matrix(Ms)
    
    if(sort(unique(A))[1]!=0 | sort(unique(A))[2]!=1) 
      stop("Treatment levels are not 0 and 1")
    
    Coefini <- coef(glm(C ~ A, family = "Gamma"(link="log"),data = data))
    RegGamma <- glm(C ~ XsC + A + A:Ms, family = "Gamma"(link="log"), data = data, 
                    start=c(Coefini[1], rep(0,dim(XsC)[2]), Coefini[2], rep(0,dim(Ms)[2])))
    
    mus.gam <- matrix(NA,N,2)
    for(j in 1:2){
      newdat <- data; newdat$A <- rep(sort(unique(A))[j],N)
      mus.gam[,j] <- predict(RegGamma, newdata=newdat, type="response")
    }
    
    out.gam <- list(mus.gam, RegGamma);names(out.gam) = c("mus.gam","Gamma")
    return(out.gam)
  }
  
  # Function: Compute IPW estimates for classification weights (long data) 
  .ITR_IPW_P <- function(FT, C, A, noncen, Mnoncen, id, MA, data, lambda, ps, muTreg, Mps, MnoncenW){
    
    N <- length(A)
    LA <- sort(unique(A))
    if(sort(unique(A))[1]!=0 | sort(unique(A))[2]!=1) 
      stop("Treatment levels are not 0 and 1")
    
    # IPW estimates for restricted survival time (FT) and cumulative cost (C)
    muT <- muC <- muT_cen <- matrix(NA,N,2)
    for(k in 1:2){
      muT[,k]  <- noncen*(A==LA[k])*FT/ps[,k]
      muC_long <- (Mnoncen/MnoncenW)*((data[,MA]==LA[k])*C/Mps[,k])
      newdat   <- data.frame(data, muC_long)
      muC[,k]  <- aggregate(newdat$muC_long, by=list(newdat[,id]), FUN=sum, na.rm=TRUE)$x
    }
    
    # IPW estimates for NMB (impute the censored subjects with regression estimates)
    tempdat <- data.frame(A, noncen, muTreg, muT)
    names(tempdat)[(dim(tempdat)[2]-3):dim(tempdat)[2]] <- c("regT0","regT1","muT0","muT1")
    
    muT_cen[,1] <- ifelse(tempdat$noncen==0, (1-tempdat$A)*tempdat$regT0, (1-tempdat$A)*tempdat$muT0)
    muT_cen[,2] <- ifelse(tempdat$noncen==0, tempdat$A*tempdat$regT1, tempdat$A*tempdat$muT1)
    ipwY <- lambda*muT_cen-muC
    
    # Estimated individual-level treatment effect 
    delta.Y <- abs(ipwY[,2]-ipwY[,1])
    
    # Optimal treatment w.r.t. NMB & effectiveness
    ITR_ipwCE <- as.numeric(I(ipwY[,2]>ipwY[,1]))
    
    output <- data.frame(delta.Y, ITR_ipwCE)
    return(output)
  }
  
 # Function: Compute AIPW estimates for classification weights (long data)
  .ITR_AIPW_P <- function(FT, C, A, noncen, Mnoncen, id, MA, data, lambda, ps, muTreg, MmuCreg, Mps, MnoncenW){
    
    N <- length(A)
    LA <- sort(unique(A))
    if(sort(unique(A))[1]!=0 | sort(unique(A))[2]!=1) 
      stop("Treatment levels are not 0 and 1")
    
    # AIPW estimates for restricted survival time (FT) and cumulative cost (C)
    muT <- muC <- muT_cen <- matrix(NA,N,2)
    for(k in 1:2){
      muT[,k]  <- noncen*((A==LA[k])*FT/ps[,k]+(1-(A==LA[k])/ps[,k])*muTreg[,k])
      muC_long <- (Mnoncen/MnoncenW)*((data[,MA]==LA[k])*C/Mps[,k]+(1-(data[,MA]==LA[k])/Mps[,k])*MmuCreg[,k])
      newdat   <- data.frame(data, muC_long)
      muC[,k]  <- aggregate(newdat$muC_long, by=list(newdat[,id]), FUN=sum, na.rm=TRUE)$x
    }
    
    # AIPW estimates for NMB (impute the censored subjects with regression estimates)
    tempdat <- data.frame(A, noncen, muTreg, muT)
    names(tempdat)[(dim(tempdat)[2]-3):dim(tempdat)[2]] <- c("regT0","regT1","muT0","muT1")
    
    muT_cen[,1] <- ifelse(tempdat$noncen==0, tempdat$regT0, tempdat$muT0)
    muT_cen[,2] <- ifelse(tempdat$noncen==0, tempdat$regT1, tempdat$muT1)
    aipwY <- lambda*muT_cen-muC
    
    # Estimated individual-level treatment effect 
    delta.Y <- abs(aipwY[,2]-aipwY[,1])
    
    # Optimal treatment w.r.t. NMB
    ITR_aipwCE <- as.numeric(I(aipwY[,2]>aipwY[,1]))
    
    output <- data.frame(delta.Y, ITR_aipwCE)
    return(output)
  }
  
  A <- data[,MA]
  
  # One-row-per-patient data 
  unique_rows <- as.numeric(rownames(data[!duplicated(data[ ,id]),]))
  X_orpp <- X[unique_rows,]
  A_orpp <- A[unique_rows]
  Xs_orpp <- Xs[unique_rows,]
  Ms_orpp <- Ms[unique_rows]
  FT_orpp <- FT[unique_rows]
  event_orpp <- event[unique_rows]
  noncen_orpp <- ifelse(FT_orpp>tau, 1, event_orpp)      # whether a subject is uncensored during tau years of follow-up 
  data_orpp <- data[unique_rows,]
  
  # propensity scores
  PS <- .EstPS(A = A_orpp, X = X_orpp)
  PSdata <- data.frame(data_orpp[,id],PS); colnames(PSdata) <- c(id,"PS0","PS1")
  data <- merge(x=PSdata, y=data, by=id, all.y=T)
  Mps <- data[,c("PS0","PS1")]
  
  # conditional mean outcomes
  RmuT <- .Reg.mu.surv(FT = FT_orpp, event = event_orpp, A = A_orpp, Xs = Xs_orpp, Ms = Ms_orpp, tau = tau, data = data_orpp)
  RmuC <- .Reg.mu.gam(C = C, A = A, XsC = XsC, Ms = Ms, data = data)$mus.gam
  
  # censoring weights for every interval (3-month)
  cubsp <- data.frame(ns(data[,interval], knots=quantile(data[,interval],prob=c(0.25,0.5,0.75))))  # create natural cubic splines of month
  names(cubsp) <- c("Month1","Month2","Month3","Month4")
  data <- data.frame(data, cubsp)
  fitcenM <- glm(Mnoncen ~ Month1+Month2+Month3+Month4+as.matrix(X)+A+A:Ms, family = "binomial", data = data)
  MnoncenW <- predict(fitcenM, type = "response")
  
  # IPW & AIPW weights
  if(wf=="ITR_IPW"){
    Out <- .ITR_IPW_P(FT = FT_orpp, C = C, A = A_orpp, noncen = noncen_orpp, Mnoncen = Mnoncen, id = id, MA = MA,  
                      data = data, lambda = lambda, ps = PS, muTreg = RmuT, Mps = Mps, MnoncenW = MnoncenW)
  }else{
    Out <- .ITR_AIPW_P(FT = FT_orpp, C = C, A = A_orpp, noncen = noncen_orpp, Mnoncen = Mnoncen, id = id, MA = MA, 
                       data = data, lambda = lambda, ps = PS, muTreg = RmuT, MmuCreg = RmuC, Mps = Mps, MnoncenW = MnoncenW)
  }
  
  # Plug in classification weight (cw) and class label (cl) to conditional random forest
  cw <- Out[,1]
  cl <- as.factor(as.character(Out[,2]))
  tmpdat   <- data.frame(X_orpp, cl, cw)
  comp.dat <- tmpdat[complete.cases(tmpdat)==T,]
  cov <- comp.dat[,-c(dim(comp.dat)[2])]
  
  fit <- cforest(cl ~ ., weights = comp.dat$cw, data = cov, control = cforest_control(ntree = ntree, mtry = mtry, maxdepth = maxdepth))
  EfficientCEITR <- as.numeric(as.character(predict(fit, type="response")))
  return(EfficientCEITR)
  
  # The 'cforest' function in 'party' package only report varimp when the case weights are 1s, i.e., do not use case weights
  #fit <- cforest(cl ~ ., data=cov, control = cforest_control(ntree = ntree, mtry = mtry, maxdepth = maxdepth))
  #EfficientCEITR <- list(NA,2)
  #EfficientCEITR[[1]] <- as.numeric(as.character(predict(fit, type="response")))
  #EfficientCEITR[[2]] <- varimp(fit, conditional = F)
}


