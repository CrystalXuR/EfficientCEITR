#####---------------------------------------------------------------------------#####
#####--------------------- Functions for simulating CE data --------------------#####
#####------------------- Effect Modification on both T and C -------------------#####

# Simulate survival time with exponential COX PH model
surT.Coxexp = function(i,trt,gammaF,St,X0){
  betaF <- c(0, 0.8, 0.8, 0.3, 0.3, 0.3)
  tmpx  <- c(1,X0[i,1],X0[i,2],X0[i,3],X0[i,4],X0[i,5])                                                     
  ht    <- 0.1*exp(t(betaF)%*%tmpx-(gammaF[1]*X0[i,1]+gammaF[2]*X0[i,2])*trt[i])       # constant hazard rate 
  SurvT <- -log(St[i])/ht                                                              # -log(S(t)) = cumulative hazard
  return(SurvT)
}

# cumulative summation tool
sumfun = function(x,start,end){
  return(sum(x[start:end]))
}

# Simulate cumulative cost with gamma dist. 
cumC.gam = function(i,trt,gammaC,I.delta,death,X0){
  betaC <- c(2, 0.02, 0.02, 0.02, 0.01, 0.01)
  tmpx  <- c(1,X0[i,1],X0[i,2],X0[i,3],X0[i,4],X0[i,5])  
  tmpB  <- exp(t(betaC)%*%tmpx+gammaC*trt[i])
  
  inital.C  <- rgamma(1,  shape=2.5, scale=tmpB)
  ongoing.C <- rgamma(40, shape=2.5, scale=0.6*tmpB)
  ongoing.C <- I.delta[i,]*ongoing.C
  ongoing.C <- ifelse(ongoing.C==0, NA, ongoing.C)              # if one is censored in month j, the cost of that month is NA
  dying.C   <- rgamma(1,  shape=2.5, scale=0.2*tmpB) 
  allcost   <- c(inital.C, ongoing.C, dying.C*death[i])         # if one is censored, we don't add the dying cost(=0)
  
  cumCost <- rep(NA,length(allcost))
  allcost[is.na(allcost)==T] <- 0
  for (j in 1:length(allcost)){
    cumCost[j] <- sumfun(allcost, 1, j)                         # accumulated cost up to month j                   
  }
  
  return(cumCost)
}


# # Simulate T according to Bang and Tsiatis 
# surT.unif = function(i,trt,gammaF){
#   X0=X0;betaF <- c(20, 0.5, 0.5, 0.5, 0.3, 0.3)
#   tmpx = c(1,X0[i,1],X0[i,2],X0[i,3],X0[i,4],X0[i,5])                                                     
#   mu = t(betaF)%*%tmpx+(gammaF[1]*X0[i,1]+gammaF[2]*X0[i,2])*trt[i]        
#   SurvT = runif(1,0,2*mu)                                                        
#   return(SurvT)
# }
# surT.exp = function(i,trt,gammaF){
#   X0=X0;betaF <- c(2, 0.2, 0.2, 0.2, 0.1, 0.1)
#   tmpx = c(1,X0[i,1],X0[i,2],X0[i,3],X0[i,4],X0[i,5])   
#   mu = t(betaF)%*%tmpx+(gammaF[1]*X0[i,1]+gammaF[2]*X0[i,2])*trt[i]             # mu is about 6
#   SurvT = exp(mu)                                        
#   return(SurvT)
# }




