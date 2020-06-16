########### This file is for merging figures from all sumulation results ##########
#---------- Date: 10/28/19
#---------- Author: Crystal Xu 
#####----------------------------------------------------------------------------#####
#####------------------- Visualization of Estimated ITR   -----------------------#####
#####----------------------------------------------------------------------------##### 
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
library(data.table);library(png);library(ggplot2);library(gridExtra);library(grid);library(lattice)
setwd("Z:/Crystal/OptimalITRPartition")
params = data.frame(read.csv("./parameters.csv"))
params$CT = as.character(params$CT)
foldernames <- c("F3_CR0.05_M","F3_CR0.05_L","F3_CR0.01_M","F3_CR0.01_L")
foldernames2 <- c("F3a_CR0.05_M","F3a_CR0.05_L","F3a_CR0.01_M","F3a_CR0.01_L")

# single figure with 8 panels 
for (k in 1:8){
    tempparm = params[k,]
    filenameDT = paste("./ResultsTC/HTE_",tempparm[1],"_",tempparm[2],"WTP_",tempparm[4],"cenrate_",tempparm[5],"ct",tempparm[6],"output_F3.csv",sep="")
    datV.DT = data.frame(t(fread(filenameDT,header=T)));dim(datV.DT)
    colnames(datV.DT) = c("x1","x2","x3","x4","x5","A", "g.opt", "Y.opt", "Y.trt", "Y.ctl",
                          "Reg.naive_deltaY_NP", "Reg.naive_ITR_NP", "aipw.Y_NP", "aipw.ITRCE_NP",
                          "ipw.ITRCE_P", "ipw.Y_P", "aipw.ITRCE_P", "aipw.Y_P","g.optt","Reg_deltaY_NP",
                          "g.aipw_NP","g.ipw_P","g.aipw_P","E0.y","E1.y","E2.y","E3.y","E4.y")
    
    ## 2-level graphs
    datV.DT$g.opt = as.factor(datV.DT$g.opt)
    datV.DT$Reg.naive_ITR_NP = as.factor(datV.DT$Reg.naive_ITR_NP)
    datV.DT$g.aipw_NP = as.factor(datV.DT$g.aipw_NP)
    datV.DT$g.ipw_P = as.factor(datV.DT$g.ipw_P)
    datV.DT$g.aipw_P = as.factor(datV.DT$g.aipw_P)
    
    #g.opt
    g.opt = qplot(x1, x2, colour = g.opt, data = datV.DT)+
      #geom_point(alpha = 1/10)+
      scale_colour_manual(values=c("#0000CC", "#FF3300"))+
      ggtitle("g.opt")+
      theme(plot.title = element_text(size=12),legend.position="bottom",
            axis.text = element_text(size=12),axis.title = element_text(size=12),
            legend.text = element_text(size=12),legend.title = element_text(size=12))
    # g.reg
    g.reg = qplot(x1, x2, colour = Reg.naive_ITR_NP, data = datV.DT)+
      scale_colour_manual(values=c("#0000CC", "#FF3300"))+
      ggtitle("g.Reg-naive")+
      theme(plot.title = element_text(size=12),legend.position="none",
            axis.text = element_text(size=12),axis.title = element_text(size=12))
    
    # g.aipw-np
    g.aipwNP.DT = qplot(x1, x2, colour = g.aipw_NP, data = datV.DT)+
      scale_colour_manual(values=c("#0000CC", "#FF3300"))+
      ggtitle("g.DT-AIPW-NP")+
      theme(plot.title = element_text(size=12),legend.position="none",
            axis.text = element_text(size=12),axis.title = element_text(size=12))
    
    # g.ow
    g.ipwP.DT = qplot(x1, x2, colour = g.ipw_P, data = datV.DT)+
      #geom_point(alpha = 1/10)+
      scale_colour_manual(values=c("#0000CC", "#FF3300"))+
      ggtitle("g.DT-IPW-P")+
      theme(plot.title = element_text(size=12),legend.position="none",
            axis.text = element_text(size=12),axis.title = element_text(size=12))
    
    # g.aipwP
    g.aipwP.DT = qplot(x1, x2, colour = g.aipw_P, data = datV.DT)+
      scale_colour_manual(values=c("#0000CC", "#FF3300"))+
      ggtitle("g.DT-AIPW-P")+
      theme(plot.title = element_text(size=12),legend.position="none",
            axis.text = element_text(size=12),axis.title = element_text(size=12))
    
    
    ###------------------------------------------RF---------------------------------------------###
    tempparm = params[k+8,]
    filenameRF = paste("./ResultsTC/HTE_",tempparm[1],"_",tempparm[2],"WTP_",tempparm[4],"cenrate_",tempparm[5],"ct",tempparm[6],"output_F3.csv",sep="")
    datV.RF = data.frame(t(fread(filenameRF,header=T)));dim(datV.RF)
    colnames(datV.RF) = c("x1","x2","x3","x4","x5","A", "g.opt", "Y.opt", "Y.trt", "Y.ctl",
                          "Reg.naive_deltaY_NP", "Reg.naive_ITR_NP", "aipw.Y_NP", "aipw.ITRCE_NP",
                          "ipw.ITRCE_P", "ipw.Y_P", "aipw.ITRCE_P", "aipw.Y_P","g.optt","Reg_deltaY_NP",
                          "g.aipw_NP","g.ipw_P","g.aipw_P","E0.y","E1.y","E2.y","E3.y","E4.y")
    
    ## 2-level graphs
    datV.RF$g.opt = as.factor(datV.RF$g.opt)
    datV.RF$Reg.naive_ITR_NP = as.factor(datV.RF$Reg.naive_ITR_NP)
    datV.RF$g.aipw_NP = as.factor(datV.RF$g.aipw_NP)
    datV.RF$g.ipw_P = as.factor(datV.RF$g.ipw_P)
    datV.RF$g.aipw_P = as.factor(datV.RF$g.aipw_P)
    
    # g.aipw-np
    g.aipwNP.RF = qplot(x1, x2, colour = g.aipw_NP, data = datV.RF)+
      scale_colour_manual(values=c("#0000CC", "#FF3300"))+
      ggtitle("g.RF-AIPW-NP")+
      theme(plot.title = element_text(size=12),legend.position="none",
            axis.text = element_text(size=12),axis.title = element_text(size=12))
    
    # g.ow
    g.ipwP.RF = qplot(x1, x2, colour = g.ipw_P, data = datV.RF)+
      #geom_point(alpha = 1/10)+
      scale_colour_manual(values=c("#0000CC", "#FF3300"))+
      ggtitle("g.RF-IPW-P")+
      theme(plot.title = element_text(size=12),legend.position="none",
            axis.text = element_text(size=12),axis.title = element_text(size=12))
    
    # g.aipwP
    g.aipwP.RF = qplot(x1, x2, colour = g.aipw_P, data = datV.RF)+
      scale_colour_manual(values=c("#0000CC", "#FF3300"))+
      ggtitle("g.RF-AIPW-P")+
      theme(plot.title = element_text(size=12),legend.position="none",
            axis.text = element_text(size=12),axis.title = element_text(size=12))
    
    # all
    mylegend<-g_legend(g.opt)
    pic3name = paste("./ResultsTC/HTE_",tempparm[1],"WTP_",tempparm[4],"cenrate_",tempparm[5],".png",sep="")
    png(pic3name, width = 7/3*4, height = 7*3/4, units = 'in', res = 300)
    print(grid.arrange(arrangeGrob(g.opt + theme(legend.position="none"), g.aipwNP.DT, g.ipwP.DT, g.aipwP.DT,
                                   g.reg, g.aipwNP.RF, g.ipwP.RF, g.aipwP.RF, nrow=2, ncol=4), mylegend, nrow=2, heights=c(10,1)))
    dev.off()
  }

# Combine 2 figures with differnt WTP but same CR & HTE
for (i in 1:4){
  setwd(paste("Z:/Crystal/OptimalITRPartition/ResultsTC/",foldernames[i],sep=""))
  files <- list.files(path=".", pattern="*.png", all.files=T, full.names=T)
  filelist <- lapply(files, readPNG)
  names(filelist) <- paste0(basename((files)))
  list2env(filelist, envir=.GlobalEnv)
  
  nameout=paste("Z:/Crystal/OptimalITRPartition/ResultsTC/Combo",foldernames[i],".png",sep="")
  png(nameout, width = 650, height = 700, res = 300)
  par(mar=rep(0,4))
  layout(matrix(1:2, ncol=1, byrow=TRUE))
  
  for(j in 1:2) {
    img <- readPNG(names(filelist[j]))
    plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
    rasterImage(img,0,0,1,1)
  }
  dev.off()
}



  



##--------------------------------- Wang (2016) DB plots ---------------------------------------------##
setwd("Z:/Crystal/P1RevisionOptimallTR")
OutW_P1 = data.frame(t(fread("P1_WangX1X2_F3/HTE_1.5Part1WTP_100000SM_PsTOutFoutput_F3_8085.csv",header=T)))
OutW_P2 = data.frame(t(fread("P1_WangX1X2_F3/HTE_1.5Part2WTP_100000SM_PsTOutFoutput_F3_8085.csv",header=T)))
OutW_P3 = data.frame(t(fread("P1_WangX1X2_F3/HTE_1.5Part3WTP_100000SM_PsTOutFoutput_F3_8085.csv",header=T)))
OutW_P4 = data.frame(t(fread("P1_WangX1X2_F3/HTE_1.5Part4WTP_100000SM_PsTOutFoutput_F3_8085.csv",header=T)))
OutW_P5 = data.frame(t(fread("P1_WangX1X2_F3/HTE_1.5Part5WTP_100000SM_PsTOutFoutput_F3_8085.csv",header=T)))
OutW_P6 = data.frame(t(fread("P1_WangX1X2_F3/HTE_1.5Part6WTP_100000SM_PsTOutFoutput_F3_8085.csv",header=T)))
OutW_P7 = data.frame(t(fread("P1_WangX1X2_F3/HTE_1.5Part7WTP_100000SM_PsTOutFoutput_F3_8085.csv",header=T)))
OutW_P8 = data.frame(t(fread("P1_WangX1X2_F3/HTE_1.5Part8WTP_100000SM_PsTOutFoutput_F3_8085.csv",header=T)))
OutW_P9 = data.frame(t(fread("P1_WangX1X2_F3/HTE_1.5Part9WTP_100000SM_PsTOutFoutput_F3_8085.csv",header=T)))
OutW_P10 = data.frame(t(fread("P1_WangX1X2_F3/HTE_1.5Part10WTP_100000SM_PsTOutFoutput_F3_8085.csv",header=T)))
datV.est = data.frame(rbind(OutW_P1,OutW_P2,OutW_P3,OutW_P4,OutW_P5,OutW_P6,OutW_P7,OutW_P8,OutW_P9,OutW_P10));dim(datV.est)
colnames(datV.est) = c("FUT.obs","RFUT.obs","cost_tot.obs","Y.obs","Y.trt","Y.ctl","Y.g.opt","x1","x2","x3","x4","x5",
                       "est_ps0","est_ps1","g.opt","g.reg","A","ipw.ITRCE","aipw.ITRCE","regB_w","owl_w",
                       "ipw.Y","aipw.Y","Y.g.opt.T","Y.g.opt.rand","Y.trt.rescale","Y.ctl.rescale","surv.prob",
                       "g.regbase","g.ow","g.ipw","g.aipw","Trt","g.BRM80","g.BRM85","g.BROL80","g.BROL85",
                       "g.BROG80","g.BROG85","E0.y","E1.y","E2.y","E3.y","E4.y","E5.y","E6.y",
                       "E7.y","E8.y","E9.y","E10.y","E11.y")

## 2-level graphs
datV.est$g.opt = as.factor(datV.est$g.opt)
datV.est$g.reg = as.factor(datV.est$g.reg)
datV.est$g.regbase = as.factor(datV.est$g.regbase)
datV.est$g.owl = as.factor(datV.est$g.ow)
datV.est$g.ipw = as.factor(datV.est$g.ipw)
datV.est$g.aipw = as.factor(datV.est$g.aipw)
datV.est$g.BRM80 = as.factor(datV.est$g.BRM80)
datV.est$g.BRM85 = as.factor(datV.est$g.BRM85)
datV.est$g.BROL80 = as.factor(datV.est$g.BROL80)
datV.est$g.BROL85 = as.factor(datV.est$g.BROL85)
datV.est$g.BROG80 = as.factor(datV.est$g.BROG80)
datV.est$g.BROG85 = as.factor(datV.est$g.BROG85)

#g.opt
g.opt.comp = qplot(x1, x2, colour = g.opt, data = datV.est)+
  #geom_point(alpha = 1/10)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  ggtitle("g.opt")+
  theme(plot.title = element_text(size=12),legend.position="top",plot.margin = rep(unit(0,"null"),4),
        axis.text = element_text(size=12),axis.title = element_text(size=12),
        legend.text = element_text(size=12),legend.title = element_text(size=12))

# g.reg
g.reg.comp = qplot(x1, x2, colour = g.reg, data = datV.est)+
  #geom_point(alpha = 1/10)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  ggtitle("g.reg-naive")+
  theme(plot.title = element_text(size=12),legend.position="none",
        axis.text = element_text(size=12),axis.title = element_text(size=12))

# g.regbase
g.regbase.comp = qplot(x1, x2, colour = g.regbase, data = datV.est)+
  #geom_point(alpha = 1/10)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  ggtitle("g.reg-based")+
  theme(plot.title = element_text(size=12),legend.position="none",
        axis.text = element_text(size=12),axis.title = element_text(size=12))

# g.ow
g.ow.comp = qplot(x1, x2, colour = g.owl, data = datV.est)+
  #geom_point(alpha = 1/10)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  ggtitle("g.owl")+
  theme(plot.title = element_text(size=12),legend.position="none",
        axis.text = element_text(size=12),axis.title = element_text(size=12))

# g.ipw
g.ipw.comp = qplot(x1, x2, colour = g.ipw, data = datV.est)+
  #geom_point(alpha = 1/10)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  ggtitle("g.ipw")+
  theme(plot.title = element_text(size=12),legend.position="none",
        axis.text = element_text(size=12),axis.title = element_text(size=12))

# g.aipw
g.aipw.comp = qplot(x1, x2, colour = g.aipw, data = datV.est)+
  #geom_point(alpha = 1/10)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  ggtitle("g.aipw")+
  theme(plot.title = element_text(size=12),legend.position="none",
        axis.text = element_text(size=12),axis.title = element_text(size=12))

# g.BRM80
g.BRM80.comp = qplot(x1, x2, colour = g.BRM80, data = datV.est)+
  #geom_point(alpha = 1/10)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  ggtitle("g.BRM80")+
  theme(plot.title = element_text(size=12),legend.position="none",
        axis.text = element_text(size=12),axis.title = element_text(size=12))

# g.BRM85
g.BRM85.comp = qplot(x1, x2, colour = g.BRM85, data = datV.est)+
  #geom_point(alpha = 1/10)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  ggtitle("g.BRM85")+
  theme(plot.title = element_text(size=12),legend.position="none",
        axis.text = element_text(size=12),axis.title = element_text(size=12))

# g.BROL80
g.BROL80.comp = qplot(x1, x2, colour = g.BROL80, data = datV.est)+
  #geom_point(alpha = 1/10)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  ggtitle("g.BROL80")+
  theme(plot.title = element_text(size=12),legend.position="none",
        axis.text = element_text(size=12),axis.title = element_text(size=12))

# g.BROL85
g.BROL85.comp = qplot(x1, x2, colour = g.BROL85, data = datV.est)+
  #geom_point(alpha = 1/10)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  ggtitle("g.BROL85")+
  theme(plot.title = element_text(size=12),legend.position="none",
        axis.text = element_text(size=12),axis.title = element_text(size=12))

# g.BROG80
g.BROG80.comp = qplot(x1, x2, colour = g.BROG80, data = datV.est)+
  #geom_point(alpha = 1/10)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  ggtitle("g.BROG80")+
  theme(plot.title = element_text(size=12),legend.position="none",
        axis.text = element_text(size=12),axis.title = element_text(size=12))

# g.BROG85
g.BROG85.comp = qplot(x1, x2, colour = g.BROG85, data = datV.est)+
  #geom_point(alpha = 1/10)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  ggtitle("g.BROG85")+
  theme(plot.title = element_text(size=12),legend.position="none",
        axis.text = element_text(size=12),axis.title = element_text(size=12))

# DB includes 80 & 85
mylegend<-g_legend(g.opt.comp)
picname = paste("P1_WangX1X2_F3/WangDB8085.png",sep="")
png(picname, width = 7, height = 11, units = 'in', res = 300)
print(grid.arrange(arrangeGrob(g.opt.comp+ theme(legend.position="none"),g.reg.comp, g.regbase.comp, g.ow.comp,g.ipw.comp,g.aipw.comp,
                               g.BRM80.comp,g.BROL80.comp,g.BROG80.comp,g.BRM85.comp,g.BROL85.comp,g.BROG85.comp,
                               nrow=4, ncol=3), mylegend, heights=c(10, 1)))
dev.off()



#####------------------- Visualization of WTP DB   ------------------------------#####
setwd("Z:/Crystal/P1RevisionOptimallTR")
library(data.table);library(tidyverse);library(egg);library(grid);library(gridExtra)
datV.est2 = data.frame(fread("./P1_ResultsX1X2_WTP_F3/HTE_0.9Part2Output_F3_UBLB.csv"))
datV.est7 = data.frame(fread("./P1_ResultsX1X2_WTP_F3/HTE_0.9Part7Output_F3_UBLB.csv"))
datV.est10 = data.frame(fread("./P1_ResultsX1X2_WTP_F3/HTE_0.9Part10Output_F3_UBLB.csv"))
datV.est = rbind(datV.est2,datV.est7,datV.est10);dim(datV.est)
datV.est$g.opt.var = as.factor(datV.est$g.opt.var)
datV.est$g.opt.UB = as.factor(datV.est$g.opt.UB)
datV.est$g.opt.LB = as.factor(datV.est$g.opt.LB)
datV.est$g.aipw = as.factor(datV.est$g.aipw)
datV.est$g.aipw_UB = as.factor(datV.est$g.aipw_UB)
datV.est$g.BRO.LinearK = as.factor(datV.est$g.BRO.LinearK)

# 2-half spheres plots 
fig.opt.var = qplot(x1, x2, colour = g.opt.var, data = datV.est)+
  geom_point(alpha = 1/10)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  ggtitle("g.opt.var")+
  theme(plot.title = element_text(size=12),legend.position="bottom",
        axis.text = element_text(size=12),axis.title = element_text(size=12),
        legend.text = element_text(size=12),legend.title = element_text(size=12))

#datV.estlam_UB = datV.est[datV.est$lambda==500000000,]
#datV.estlam_LB = datV.est[datV.est$lambda==40000000,]
datV.estlam_UB = datV.est[datV.est$x1<=median(datV.est$x1) & datV.est$x2>=-16 & datV.est$x2<=20 & datV.est$x1>=-16&datV.est$x1<=23,]
datV.estlam_LB = datV.est[datV.est$x1>median(datV.est$x1) & datV.est$x2>=-16 & datV.est$x2<=20 & datV.est$x1>=-16&datV.est$x1<=23,]

#fig.opt.var
P1 = qplot(x1, x2, colour = g.opt.var, data = datV.estlam_UB)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="none",
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm"))
P2 = qplot(x1, x2, colour = g.opt.var, data = datV.estlam_LB)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  scale_y_continuous(position = "right")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="none",
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm"))
fig.opt.var2 = grid.arrange(P1,P2, ncol=2,left=textGrob("x2"));fig.opt.var2

#fig.est.aipw
P7 = qplot(x1, x2, colour = g.aipw, data = datV.estlam_UB)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="none",panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm"))
P8 = qplot(x1, x2, colour = g.aipw, data = datV.estlam_LB)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  scale_y_continuous(position = "right")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="none",panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm"))
fig.est.aipw2 = grid.arrange(P7,P8, ncol=2,left=textGrob("x2"))

#fig.est.BRO1
P9 = qplot(x1, x2, colour = g.aipw_UB, data = datV.estlam_UB)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="none",panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm"))
P10 = qplot(x1, x2, colour = g.aipw_UB, data = datV.estlam_LB)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  scale_y_continuous(position = "right")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="none",panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm"))
fig.est.aipwfixed2 = grid.arrange(P9,P10, ncol=2,left=textGrob("x2"))

#fig.est.BRO2
P11 = qplot(x1, x2, colour = g.BRO.LinearK, data = datV.estlam_UB)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="none",panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm"))
P12 = qplot(x1, x2, colour = g.BRO.LinearK, data = datV.estlam_LB)+
  scale_colour_manual(values=c("#0000CC", "#FF3300"))+
  scale_y_continuous(position = "right")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="none",panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm"))
fig.est.BRO2 = grid.arrange(P11,P12, ncol=2,left=textGrob("x2"))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(fig.opt.var)
pic3name = paste("./P1_ResultsX1X2_WTP_F3/halfsWTPSolid_UB_031220.png",sep="")
png(pic3name, width = 6, height = 7*3/4, units = 'in', res = 300)
print(grid.arrange(arrangeGrob(fig.opt.var2,fig.est.aipw2, fig.est.aipwfixed2, fig.est.BRO2,
                               nrow=2, ncol=2), mylegend, heights=c(10, 1)))
dev.off()