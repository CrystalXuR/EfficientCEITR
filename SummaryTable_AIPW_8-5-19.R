options(scipen = 999)
setwd("Z:/Crystal/OptimalITRPartition")

# params 
HTE = c(1.3,1.6)
WTP = c(50000,100000)
cenRate = c(0.01,0.03)
costModel = c("gamma","lognormal")
classTools = c("DT","RF")
params = expand.grid(HTE,WTP,cenRate,costModel,classTools)
colnames(params) = c("HTE","WTP","cenRate","costModel","classTools")

# import the raw results
parameters = read.csv("./parameters.csv")
CCR_raw=NULL
for (k in 1:8){
  filenames=paste("Results_slurm_F3a/HTE_",parameters$HTE1[k],"_",parameters$HTE2[k],"WTP_",parameters$WTP[k],"cenrate_",parameters$cenRate[k],"CCR_F3.csv",sep="")
  CCR_out = read.csv(filenames)
  CCR_raw = rbind(data.frame(rep(parameters$HTE2[k],dim(CCR_out)[1]),rep(parameters$WTP[k],dim(CCR_out)[1]),rep(parameters$cenRate[k],dim(CCR_out)[1]),CCR_out),CCR_raw)
}
colnames(CCR_raw)[1:3] = c("HTE","WTP","cenRate");head(CCR_raw)

EMO_raw=NULL
for (i in 1:8){
  filenames=paste("Results_slurm_F3a/HTE_",parameters$HTE1[i],"_",parameters$HTE2[i],"WTP_",parameters$WTP[i],"cenrate_",parameters$cenRate[i],"EMO_F3.csv",sep="")
  EMO_out = read.csv(filenames)
  EMO_raw = rbind(data.frame(rep(parameters$HTE2[i],dim(CCR_out)[1]),rep(parameters$WTP[i],dim(CCR_out)[1]),rep(parameters$cenRate[i],dim(CCR_out)[1]),EMO_out),EMO_raw)
}
colnames(EMO_raw)[1:3] = c("HTE","WTP","cenRate");head(EMO_raw)


# Summarize CCR
ncol=5
SumCCR=matrix(NA,dim(params)[1],ncol)
for (i in 1:dim(params)[1]){
  datCCR = CCR_raw[CCR_raw[,1]==params[i,1]&CCR_raw[,2]==params[i,2]&CCR_raw[,3]==params[i,3]
                   &CCR_raw[,5]==params[i,4]&CCR_raw[,6]==params[i,5],7:(7+ncol-1)]
  Out = rep(NA,ncol)
  for (k in 1:ncol){
    Out[k]=paste(round(colMeans(datCCR)[k],digits=4)," (",round(var(datCCR[k]),digits = 4),")",sep="")
  }
  SumCCR[i,]=Out
}
SumCCR = data.frame(params[,2],params[,1],params[,3],params[,5],params[,4],SumCCR)
colnames(SumCCR) = c("WTP","HTE","cenRate","classTools","costModel","Truth", "C_RegNaive_NP","C_AIPW_NP","C_IPW_P","C_AIPW_P")
SumCCR = SumCCR[order(SumCCR$WTP,SumCCR$HTE,SumCCR$cenRate,SumCCR$classTools,SumCCR$costModel),];head(SumCCR)
write.csv(SumCCR, "./Results_slurm_F3a/SumCCR.csv")


# Summarize EMO
ncol=5
SumEMO = matrix(NA,dim(params)[1],ncol)
for (i in 1:dim(params)[1]){
  datEMO = EMO_raw[EMO_raw[,1]==params[i,1]&EMO_raw[,2]==params[i,2]&EMO_raw[,3]==params[i,3]&
                     EMO_raw[,5]==params[i,4]&EMO_raw[,6]==params[i,5],7:(7+ncol-1)]
  Out = rep(NA,ncol)
  for (k in 1:ncol){
    Out[k]=paste(round(colMeans(datEMO)[k],digits = 2)," (",round(var(datEMO[k]),digits = 2),")",sep="")
  }
  SumEMO[i,]=Out
}
SumEMO = data.frame(params[,2],params[,1],params[,3],params[,5],params[,4],SumEMO)
colnames(SumEMO) = c("WTP","HTE","cenRate","classTools","costModel","Truth","E_RegNaive_NP","E_AIPW_NP","E_IPW_P","E_AIPW_P")
SumEMO = SumEMO[order(SumEMO$WTP,SumEMO$HTE,SumEMO$cenRate,SumEMO$classTools,SumEMO$costModel),];head(SumEMO)
write.csv(SumEMO, "./Results_slurm_F3a/SumEMO.csv")


#####----------------------------------------------------------------------------#####
#####------------------- Visualization of Estimated ITR   -----------------------#####
#####----------------------------------------------------------------------------##### 
# we only plot for cup=0.03, CM=LN, 
library(data.table); library(ggplot2);library(gridExtra)
parameters = read.csv("parameters.csv");parameters = parameters[parameters$cenRate==0.03,]
for (k in 1:4){
  filenames1=paste("Results_slurm_F3/HTE_",parameters$HTE1[k],"_",parameters$HTE2[k],"WTP_",parameters$WTP[k],"cenrate_",parameters$cenRate[k],"output_F3.csv",sep="")
  output_raw = fread(filenames1)
  output_rawt = data.frame(t(output_raw))
  filenames2 = paste("Results_slurm_F3/HTE_",parameters$HTE1[k],"_",parameters$HTE2[k],"WTP_",parameters$WTP[k],"cenrate_",parameters$cenRate[k],"CCR_F3.csv",sep="")
  CCR_raw = read.csv(filenames2)
  
  dat = data.frame(rep(0.01,dim(CCR_raw)[1]),rep(1.6,dim(CCR_raw)[1]),rep(100000,dim(CCR_raw)[1]),CCR_raw[,2:3])
  names(dat)[1:3] = c("cenRate","HTE","WTP")
  datt = dat[rep(seq_len(nrow(dat)), each=1000),]
  output = data.frame(output_rawt,datt);dim(output) # 14400000  30
  colnames(output)[1:28] = c("Y.trt", "Y.ctl", "Y.g.opt", "X1", "X2", "X3", "X4", "X5", "A", "g.opt",
                             "ipw.ITRCE_P","ipw.Y_P", "aipw.ITRCE_P", "aipw.Y_P", "Reg.naive_ITR_NP","Reg.naive_deltaY_NP",
                             "true.ww","Reg_deltaY_NP","g.optt","g.true_NP","g.aipw_NP","g.ipw_P","g.aipw_P",
                             "E0.y","E1.y","E2.y","E3.y","E4.y")
  
  rawout = output[output$costModel=="gamma",]
  datV.est = rawout[,c("X1","X2","classTools","g.opt","Reg.naive_ITR_NP","g.ipw_P","g.aipw_P")];head(datV.est) #100000 5
  datV.estDT = datV.est[datV.est$classTools=="DT",c(1:2,4:7)];names(datV.estDT)[5:6] = c("g.P-ipw_DT","g.P-aipw_DT")
  datV.estRF = datV.est[datV.est$classTools=="RF",6:7];names(datV.estRF) = c("g.P-ipw_RF","g.P-aipw_RF")
  datV.est2 = data.frame(datV.estDT,datV.estRF);head(datV.est2) #100000 8
  
  ## 2-level graphs
  datV.est2$g.opt = as.factor(as.numeric(as.character(datV.est2$g.opt)));table(datV.est2$g.opt)
  datV.est2$Reg.naive_ITR_NP = as.factor(as.numeric(as.character(datV.est2$Reg.naive_ITR_NP)))
  datV.est2$g.P.ipw_DT = as.factor(as.numeric(as.character(datV.est2$g.P.ipw_DT)))
  datV.est2$g.P.aipw_DT = as.factor(as.numeric(as.character(datV.est2$g.P.aipw_DT)))
  datV.est2$g.P.ipw_RF = as.factor(as.numeric(as.character(datV.est2$g.P.ipw_RF)))
  datV.est2$g.P.aipw_RF = as.factor(as.numeric(as.character(datV.est2$g.P.aipw_RF)))
  datV.est2$X1 = as.numeric(as.character(datV.est2$X1))
  datV.est2$X2 = as.numeric(as.character(datV.est2$X2))
  
  #g.opt
  g.opt = qplot(X1, X2, colour = g.opt, data = datV.est2)+
    scale_colour_manual(values=c("#0000CC", "#FF3300"))+
    ggtitle("g.opt")+
    theme(plot.title = element_text(size=12),legend.position="bottom",
          axis.text = element_text(size=12),axis.title = element_text(size=12),
          legend.text = element_text(size=12),legend.title = element_text(size=12))
 
  #g.reg
  g.reg = qplot(X1, X2, colour = Reg.naive_ITR_NP, data = datV.est2)+
    scale_colour_manual(values=c("#0000CC", "#FF3300"))+
    ggtitle("Reg")+
    theme(plot.title = element_text(size=12),legend.position="none",
          axis.text = element_text(size=12),axis.title = element_text(size=12))
  
  # g.P-ipw_DT
  g.ipw_P.DT = qplot(X1, X2, colour = g.P.ipw_DT, data = datV.est2)+
    scale_colour_manual(values=c("#0000CC", "#FF3300"))+
    ggtitle("P-IPW-DT")+
    theme(plot.title = element_text(size=12),legend.position="none",
          axis.text = element_text(size=12),axis.title = element_text(size=12))
  
  # g.P-ipw_RF
  g.ipw_P.RF = qplot(X1, X2, colour = g.P.ipw_RF, data = datV.est2)+
    scale_colour_manual(values=c("#0000CC", "#FF3300"))+
    ggtitle("P-IPW-RF")+
    theme(plot.title = element_text(size=12),legend.position="none",
          axis.text = element_text(size=12),axis.title = element_text(size=12))
  
  # g.P-aipw_DT
  g.aipw_P.DT = qplot(X1, X2, colour = g.P.aipw_DT, data = datV.est2)+
    scale_colour_manual(values=c("#0000CC", "#FF3300"))+
    ggtitle("P-AIPW-DT")+
    theme(plot.title = element_text(size=12),legend.position="none",
          axis.text = element_text(size=12),axis.title = element_text(size=12))
  
  # g.P-aipw_RF
  g.aipw_P.RF = qplot(X1, X2, colour = g.P.aipw_RF, data = datV.est2)+
    scale_colour_manual(values=c("#0000CC", "#FF3300"))+
    ggtitle("P-AIPW-RF")+
    theme(plot.title = element_text(size=12),legend.position="none",
          axis.text = element_text(size=12),axis.title = element_text(size=12))
  
  
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(g.opt))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  mylegend<-g_legend(g.opt)
  filename = paste("Results_slurm_F3/WTP",parameters$WTP[k],"_HTE",parameters$HTE1[k],"_cenRate",parameters$cenRate[k],output$costModel,".png",sep="")
  png(filename, width = 7, height = 7*3/4, units = 'in', res = 300)
  print(grid.arrange(arrangeGrob(g.opt + theme(legend.position="none"), 
                                 g.ipw_P.DT,g.aipw_P.DT,g.reg,g.ipw_P.RF,g.aipw_P.RF,
                                 nrow=2, ncol=3),mylegend,nrow=2,heights=c(10, 1)))
  dev.off()
}

