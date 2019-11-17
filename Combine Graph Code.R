########### This file is for merging figures from all sumulation results ##########
#---------- Date: 1/26/19
#---------- Author: Crystal Xu 

library("png") # for reading in PNGs
setwd("Z:/Crystal/OptimalITRPartition/Results_slurm_F3/F3_gamma_0.03")
files <- list.files(path=".", pattern="*.png", all.files=T, full.names=T)
filelist <- lapply(files, readPNG)
names(filelist) <- paste0(basename((files)))
list2env(filelist, envir=.GlobalEnv)

png("F3CombinedFigOutFPSTCen0.03gamma.png", width = 1100, height = 700, res = 300)
par(mar=rep(0,4))
layout(matrix(1:4, ncol=2, byrow=TRUE))

for(i in 1:4) {
  img <- readPNG(names(filelist[i]))
  plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
  rasterImage(img,0,0,1,1)
}
dev.off()

# write to PDF
#filename = paste(names(filelist),sep = "")
#dev.print(png, "F3aCombinedFig.png")


## Merge two triangle plots 
# setwd("C:/Users/cryst/Box Sync/A_Crystal's Thesis/1 - Optimization-Classification/2 - Simulations & implementation/Simulation Codes and Results/1-CEanalysisITRcode/PubFigures/Triplots")
# files <- list.files(path=".", pattern="*.png", all.files=T, full.names=T)
# filelist <- lapply(files, readPNG)
# names(filelist) <- paste0(basename((files)))
# list2env(filelist, envir=.GlobalEnv)
# 
# png("F3aCombinedFigTrue2L.png", width = 1100, height = 700, res = 300)
# par(mar=rep(0,4))
# layout(matrix(1:2, ncol=2, byrow=TRUE))
# 
# for(i in 1:2) {
#   img <- readPNG(names(filelist[i]))
#   plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
#   rasterImage(img,0,0,1,1)
# }
# dev.off()
