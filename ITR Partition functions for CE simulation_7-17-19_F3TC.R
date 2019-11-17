#####---------------------------------------------------------------------------#####
#####--------------------- Functions for simulating CE data --------------------#####
#####------------------- Effect Modification on both T and C -------------------#####

# Simulate survival time with exponential COX PH model  
surT.exp = function(i,trt_i,h0,betaF,gammaF){
  tmpx = c(1, x1[i],x2[i],x3[i],x4[i],x5[i])                                                     
  haz.rate = h0 * exp(t(betaF) %*% tmpx - (gammaF[1]*x1[i]+gammaF[2]*x2[i])*trt_i[i])       # hazard at time t, does not depend on time t since h0 is constant
  return(-log(surv.prob[i])/(haz.rate))                                                     # -log(S(t)) = cumulative hazard. Survival T = Cumulative hazard/hazard at time t
}

# cumulative summation tool
sumfun = function(x,start,end){
  return(sum(x[start:end]))
}

# simulate cumulative cost with gamma dist. 
cumC.gam = function(i, trt_i, k, betaC,gammaC,I.delta,event,tau){
  tmpx = c(1,x1[i],x2[i],x3[i],x4[i],x5[i]) 
  tmp.beta = exp(t(betaC) %*% tmpx + (gammaC[1]*x1[i]+gammaC[2]*x2[i])*trt_i[i])
  
  inital.C = rgamma(1,shape = k,scale=tmp.beta)
  ongoing.C = rgamma((tau*12), shape = k, scale=0.1*tmp.beta)
  ongoing.C = I.delta[i,]*ongoing.C
  ongoing.C = ifelse(ongoing.C==0, NA,ongoing.C)
  dying.C = rgamma(1, shape = k, scale=0.2*tmp.beta)                    # if one is censored, we don't add the dying cost(=0)
  
  cumC = rep(NA,(tau*12+1))
  for (j in 1:(tau*12+1)){
    cost.series = c(inital.C,ongoing.C)           # combine initial cost and ongoing cost into one vector 
    cumC[j] = sumfun(cost.series, 1,j)            # cost accured up to each month
  }
  
  ni = sum(!is.na(cumC))                         # count the number of month with a cost
  cumC[ni] = cumC[ni]+dying.C*event[i]           # add the cumulative cost in the last month to dying cost
  cumCost = cumC
  return(cumCost)
}

# simulate cumulative cost with log-normal dist. 
# cumC.lognorm = function(i, trt_i, k, betaC,gammaC,I.delta,event){
#   tmpx = c(1,x1[i],x2[i],x3[i],x4[i],x5[i]) 
#   tmp.beta = t(betaC) %*% tmpx + (gammaC[1]*x1[i]+gammaC[2]*x2[i])*trt_i[i]
#   
#   inital.C = exp(rnorm(1,mean=2*tmp.beta,sd = k))
#   ongoing.C = exp(rnorm(48,mean=0.2*tmp.beta,sd = 0.1*k))
#   ongoing.C = I.delta[i,]*ongoing.C
#   ongoing.C = ifelse(ongoing.C==0, NA,ongoing.C)
#   dying.C = exp(rnorm(1,mean=0.4*tmp.beta,sd = 0.1*k))                    # if one is censored, we don't add the dying cost(=0)
#   
#   cumC = rep(NA,49)
#   for (j in 1:49){
#     cost.series = c(inital.C,ongoing.C)           # combine initial cost and ongoing cost into one vector 
#     cumC[j] = sumfun(cost.series, 1,j)            # cost accured up to each month
#   }
#   
#   ni = sum(!is.na(cumC))                         # count the number of month with a cost
#   cumC[ni] = cumC[ni]+dying.C*event[i]           # add the cumulative cost in the last month to dying cost
#   cumCost = cumC
#   return(cumCost)
# }

# use only one legend for plots
g_legend = function(a.gplot){
  tmp  =  ggplot_gtable(ggplot_build(a.gplot))
  leg  =  which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend  =  tmp$grobs[[leg]]
  return(legend)}

