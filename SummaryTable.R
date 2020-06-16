
##------------------------------ Our simulation with Nsim=500 ----------------------------------------------------------------------###
#setwd("Z:/Crystal/OptimalITRPartition")
setwd("C:/Users/cryst/Box Sync/A_Crystal's Thesis/2 - OptimalITRPartition")
params = read.csv("./parameters.csv")
CCRsum = NULL
EMOsum = NULL
for (i in 1:16){
  tempparm = as.character(as.matrix(params[i,]))
  CCR=read.csv(paste("./ResultsTC/HTE_",tempparm[1],"_",tempparm[2],"WTP_",tempparm[4],"cenrate_",tempparm[5],"ct",tempparm[6],"CCR_F3.csv",sep=""))
  EMO=read.csv(paste("./ResultsTC/HTE_",tempparm[1],"_",tempparm[2],"WTP_",tempparm[4],"cenrate_",tempparm[5],"ct",tempparm[6],"EMO_F3.csv",sep=""))
  CCRsum = rbind(CCRsum,unlist(c(tempparm,round(colMeans(CCR[,3:dim(CCR)[2]],dims=1),digits=3),round(apply(CCR[,3:dim(CCR)[2]],2,sd),digits=3))))
  EMOsum = rbind(EMOsum,unlist(c(tempparm,round(colMeans(EMO[,3:dim(EMO)[2]],dims=1),digits=3),round(apply(EMO[,3:dim(EMO)[2]],2,sd),digits=3))))
}
write.csv(CCRsum,"ResultsTC/F3CCRsum.csv")
write.csv(EMOsum,"ResultsTC/F3EMOsum.csv")

CCRsum = NULL
EMOsum = NULL
for (i in 1:16){
  tempparm = as.character(as.matrix(params[i,]))
  CCR=read.csv(paste("./ResultsT/HTE_",tempparm[1],"_",tempparm[2],"WTP_",tempparm[4],"cenrate_",tempparm[5],"ct",tempparm[6],"CCR_F3a.csv",sep=""))
  EMO=read.csv(paste("./ResultsT/HTE_",tempparm[1],"_",tempparm[2],"WTP_",tempparm[4],"cenrate_",tempparm[5],"ct",tempparm[6],"EMO_F3a.csv",sep=""))
  CCRsum = rbind(CCRsum,unlist(c(tempparm,round(colMeans(CCR[,3:dim(CCR)[2]],dims=1),digits=3),round(apply(CCR[,3:dim(CCR)[2]],2,sd),digits=3))))
  EMOsum = rbind(EMOsum,unlist(c(tempparm,round(colMeans(EMO[,3:dim(EMO)[2]],dims=1),digits=3),round(apply(EMO[,3:dim(EMO)[2]],2,sd),digits=3))))
}
write.csv(CCRsum,"ResultsT/F3aCCRsum.csv")
write.csv(EMOsum,"ResultsT/F3aEMOsum.csv")





#####----------------------------------------------------------------------------#####
#####------------------- Visualization of Estimated ITR   -----------------------#####
#####----------------------------------------------------------------------------##### 
library(data.table)
output_raw = fread('PartitionParty Output F3_6_5_19/outputALL_F3.csv', header = T, sep = ','); dim(output_raw) # 23 14400001
output_raw = output_raw[,2:dim(output_raw)[2]]
output_rawt = t(output_raw);dim(output_rawt)  # 14400000 23

dat = CCR_raw[,1:7]
datt = dat[rep(seq_len(nrow(dat)), each=1000),]
output = data.frame(output_rawt,datt);dim(output) # 14400000  30
colnames(output) = c("Y.trt", "Y.ctl", "Y.g.opt", "X1", "X2", "X3", "X4", "X5", "A",
                     "g.opt", "ipw.ITRCE_P", "aipw.ITRCE_P", "ipw.Y_P", "aipw.Y_P",
                     "g.ipw_NP","g.aipw_NP","g.ipw_P","g.aipw_P","E0.y","E1.y","E2.y","E3.y","E4.y")

for (i in 1:dim(params)[1]){
  rawout = output[output$HTE==params[i,1]& output$WTP==params[i,2]& output$cenRate==params[i,3]&
                  output$SM==params[i,4]& output$costModel==params[i,5]& output$classTools==params[i,6],]
  datV.est = rawout[,c(10,15:18)];head(datV.est) #100000 5
  
  ## 2-level graphs
  datV.est$g.opt = as.factor(datV.est$g.opt)
  datV.est$g.ipw_NP = as.factor(datV.est$g.ipw_NP)
  datV.est$g.aipw_NP = as.factor(datV.est$g.aipw_NP)
  datV.est$g.ipw_P = as.factor(datV.est$g.ipw_P)
  datV.est$g.aipw_P = as.factor(datV.est$g.aipw_P)
  
  #g.opt
  g.opt.comp = qplot(x1, x2, colour = g.opt, data = datV.est)+
    scale_colour_manual(values=c("#0000CC", "#FF3300"))+
    ggtitle("g.opt")+
    theme(plot.title = element_text(size=12),legend.position="bottom",
          axis.text = element_text(size=12),axis.title = element_text(size=12),
          legend.text = element_text(size=12),legend.title = element_text(size=12))
 
   # g.ipw_NP
  g.ipw_NP.comp = qplot(x1, x2, colour = g.ipw_NP, data = datV.est)+
    scale_colour_manual(values=c("#0000CC", "#FF3300"))+
    ggtitle("g.IPW_NP")+
    theme(plot.title = element_text(size=12),legend.position="none",
          axis.text = element_text(size=12),axis.title = element_text(size=12))
  
  # g.aipw_NP
  g.aipw_NP.comp = qplot(x1, x2, colour = g.aipw_NP, data = datV.est)+
    scale_colour_manual(values=c("#0000CC", "#FF3300"))+
    ggtitle("g.AIPW_NP")+
    theme(plot.title = element_text(size=12),legend.position="none",
          axis.text = element_text(size=12),axis.title = element_text(size=12))
  
  # g.ipw_P
  g.ipw_P.comp = qplot(x1, x2, colour = g.ipw_P, data = datV.est)+
    scale_colour_manual(values=c("#0000CC", "#FF3300"))+
    ggtitle("g.IPW_P")+
    theme(plot.title = element_text(size=12),legend.position="none",
          axis.text = element_text(size=12),axis.title = element_text(size=12))
  
  # g.aipw_P
  g.aipw_P.comp = qplot(x1, x2, colour = g.aipw_P, data = datV.est)+
    scale_colour_manual(values=c("#0000CC", "#FF3300"))+
    ggtitle("g.AIPW_P")+
    theme(plot.title = element_text(size=12),legend.position="none",
          axis.text = element_text(size=12),axis.title = element_text(size=12))
  
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(g.opt.comp))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  mylegend<-g_legend(g.opt.comp)
  filename = paste("WTP",rawout$WTP[1],"_HTE",rawout$HTE[1],"_cenRate",rawout$cenRate[1],"_",
                   rawout$classTools[1],"_",rawout$SM[1],"_",rawout$costModel[1],".png",sep="")
  png(paste("PartitionParty Output F3_6_5_19/",filename,sep=""), width = 7, height = 7*3/4, units = 'in',res = 300)
  print(grid.arrange(arrangeGrob(g.opt.comp + theme(legend.position="none"), 
                     g.ipw_NP.comp,g.aipw_NP.comp, g.ipw_P.comp,g.aipw_P.comp,
                     nrow=2, ncol=3), mylegend,nrow=2,heights=c(10, 1)))
  dev.off()
}




