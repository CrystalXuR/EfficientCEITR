##---------------------------------------------- Analysis --------------------------------------------------##
# Load the "EfficientCEITR" function 
source("./EfficientCEITR.R")

# Load the example data
exampledata <- read.csv("example_data.csv")

# Run the analysis
OutCEITR <- EfficientCEITR(X = exampledata[,c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10")],
                           Xs = exampledata[,c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10")],
                           Ms = exampledata$x2,
                           XsC = exampledata[,c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","interval")],
                           FT = exampledata$FT,
                           event = exampledata$event,
                           C = exampledata$cost_mo,
                           Mnoncen = 1-exampledata$cen_mo,
                           tau = exampledata$tau[1],
                           lambda = exampledata$lambda[1], 
                           wf = "ITR_AIPW", 
                           ntree = 50, 
                           mtry = 5, 
                           maxdepth = 10, 
                           id = "Subjectid", 
                           MA = "A", 
                           interval = "interval",
                           data = exampledata)

# Estimated optimal CE ITR
est_g <- OutCEITR

# # Compare with the true optimal regime in the example data 
# unique_rows <- as.numeric(rownames(exampledata[!duplicated(exampledata[ ,"Subjectid"]),]))
# g.opt <- exampledata$g_opt[unique_rows]
# table(g.opt,est_g)