source("src/R/Simulation Study/simCommon.R")

set.seed(1)
testIdOfInterest = testDs.id$amctx[sample(1:nrow(testDs.id), size=60, replace = F)]
testIdOfInterest = sort(testIdOfInterest, decreasing = F)

simTestDs.id = testDs.id[testDs.id$amctx %in% testIdOfInterest,]
trainingSize = nrow(trainingDs.id)

#Plot their longitudinal and survival times
#multiplot(plotlist = lapply(simTestDs.id$amctx, plotTrueLongitudinal), cols = 3)
#multiplot(plotlist = lapply(simTestDs.id$amctx, plotTrueSurvival), cols = 3)

testTimes = seq(0, 3.5, 0.01)

#First method to get dynamic max risk time
cl = makeCluster(detectCores())
registerDoParallel(cl)

dynamicTrueSurvProb = foreach(id=testIdOfInterest, .packages = c("splines", "JMbayes")) %dopar%{
  
  ds_i = simDs.id[id, ]
  wGamma_i = wGamma[id]
  b_i = b_creatinine[id, ]
  
  sapply(testTimes + maxRiskDt, survivalFunc, i=id, ds_i, wGamma_i, b_i) / sapply(testTimes, survivalFunc, i=id, ds_i, wGamma_i, b_i)
}

stopCluster(cl)

simTestDs.id$stoptime_True = sapply(dynamicTrueSurvProb, FUN = function(x){testTimes[which(x <= minSurv)[1]]})
simTestDs.id = simTestDs.id[!is.na(simTestDs.id$stoptime_True),]
simTestDs.id = simTestDs.id[order(simTestDs.id$stoptime_True, decreasing = F)[1:50],]
simTestDs.id = simTestDs.id[order(simTestDs.id$amctx, decreasing = F),]
testIdOfInterest = simTestDs.id$amctx

save(simTestDs.id, file="Rdata/creatinine_sim_6month_2pt5_percentrisk.Rdata")

