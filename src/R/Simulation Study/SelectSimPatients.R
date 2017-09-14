source("src/R/Simulation Study/simCommon.R")

testIdOfInterest = testDs.id$amctx

trainingSize = nrow(trainingDs.id)

#Plot their longitudinal and survival times
multiplot(plotlist = lapply(simTestDs.id$amctx, plotTrueLongitudinal), cols = 3)
multiplot(plotlist = lapply(simTestDs.id$amctx, plotTrueSurvival), cols = 3)

#The maximum risk for now is 10%. i.e. we want survival prob as 90% at 12 months
maxRisk = 0.1
minSurv = 1-maxRisk
maxRiskDt = 1

testTimes = seq(0, 17, 0.01)

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

testDs.id$stoptime_True = sapply(dynamicTrueSurvProb, FUN = function(x){testTimes[which(x <= minSurv)[1]]})


#Create two data sets, simTestDs, simTestDs.id

# #Second method to get dynamic max risk time
# 
# #We can split the entire time range into very small times and 
# #check where exactly is the maxRisk with maxRiskDt not possible
# testTimes = seq(0, 17, 0.01)
# fineSimTestDs = data.frame(amctx = rep(simTestDs.id$amctx, each=length(testTimes)), 
#                    rec_age_fwp1 = rep(simTestDs.id$rec_age, each=length(testTimes)),
#                    rec_age = rep(simTestDs.id$rec_age, each=length(testTimes)),
#                    d_age = rep(simTestDs.id$d_age, each=length(testTimes)),
#                    rec_bmi = rep(simTestDs.id$rec_bmi, each=length(testTimes)),
#                    tx_dial_days = rep(simTestDs.id$tx_dial_days, each=length(testTimes)),
#                    tx_previoustx = rep(simTestDs.id$tx_previoustx, each=length(testTimes)),
#                    d_gender = rep(simTestDs.id$d_gender, each=length(testTimes)),
#                    rec_gender = rep(simTestDs.id$rec_gender, each=length(testTimes)),
#                    d_cadaveric = rep(simTestDs.id$d_cadaveric, each=length(testTimes)),
#                    tx_dgf = rep(simTestDs.id$tx_dgf, each=length(testTimes)),
#                    tx_dm = rep(simTestDs.id$tx_dm, each=length(testTimes)),
#                    ah_nr = rep(simTestDs.id$ah_nr, each=length(testTimes)),
#                    tx_pra = rep(simTestDs.id$tx_pra, each=length(testTimes)),
#                    visitNumber = rep(1:length(testTimes), nrow(simTestDs.id)),
#                    tx_s_years = rep(testTimes, nrow(simTestDs.id)))
# 
# fineSimTestDs$logCreatinine = c(do.call(cbind, lapply(simTestDs.id$amctx, rLogCreatinine, testTimes, T)))
# 
# cl = makeCluster(8)
# registerDoParallel(cl)
# dynamicTrueSurvProb = foreach(time=testTimes, .packages = c("splines", "JMbayes")) %dopar%{
#   patientDs = fineSimTestDs[fineSimTestDs$tx_s_years <= time,]                        
#   survfitJM(simJointModel_replaced, patientDs, idVar="amctx", 
#             survTimes = max(patientDs$tx_s_years)+maxRiskDt)
# }
# stopCluster(cl)
# 
# simTestDs.id$dynamicMaxRiskTime = NA
# for(i in 1:length(testTimes)){
#   riskProbs = 1 - sapply(1:nrow(simTestDs.id), function(index){
#     dynamicTrueSurvProb[[i]]$summaries[[index]][1, "Mean"]
#   })
#   
#   indexOfInterest = is.na(simTestDs.id$dynamicMaxRiskTime) & riskProbs >= maxRisk
#   simTestDs.id$dynamicMaxRiskTime[indexOfInterest] = testTimes[i]
# }
