source("src/R/Simulation Study/simCommon.R")

set.seed(1)
testIdOfInterest = testDs.id$amctx[sample(1:nrow(testDs.id), size=60, replace = F)]
testIdOfInterest = sort(testIdOfInterest, decreasing = F)

simTestDs.id = testDs.id[testDs.id$amctx %in% testIdOfInterest,]
trainingSize = nrow(trainingDs.id)

#Plot their longitudinal and survival times
#multiplot(plotlist = lapply(simTestDs.id$amctx, plotTrueLongitudinal), cols = 3)
#multiplot(plotlist = lapply(simTestDs.id$amctx, plotTrueSurvival), cols = 3)

#The maximum risk for now is 10%. i.e. we want survival prob as 90% at 12 months
#maxRisk = 0.1
maxRisk = 0.05
minSurv = 1-maxRisk
#maxRiskDt = 1
maxRiskDt = 0.5

testTimes = seq(0, 20, 0.01)

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

save(simTestDs.id, file="Rdata/creatinine_sim_6month_7.5_percentrisk.Rdata")

#############################################
#Load the results from Rdata files, 6 months 5% risk
#############################################
i = 1
# minObsArr = c(6,6,10,10,15,15) #in the same order as below
# for(rdataname in c("u1_dynInfoPar_6mo_nFix_6_k100.Rdata", "u1_dynInfo_mod_6mo_nFix_6_k50.Rdata",
#                    "u1_dynInfoPar_6mo_nFix_10_k100.Rdata","u1_dynInfo_mod_6mo_nFix_10_k50.Rdata",
#                    "u1_dynInfoPar_6mo_nFix_15_k50.Rdata", "u1_dynInfo_mod_6mo_nFix_15_k50.Rdata")){
#   load(paste("Rdata/", rdataname, sep=""))
minObsArr = c(6,10,10,10,10) #in the same order as below
for(rdataname in c("u1_dynInfoPar_6mo_nFix_8_k50", "u1_dynInfo_mod_6mo_nFix_8_k50")){
  load(paste("Rdata/", rdataname, sep=""))

  simTestDs.id[,paste("nObs_",i, sep="")] = sapply(patientDsList, nrow)
  simTestDs.id[,paste("stopTime_",i, sep="")] = sapply(patientDsList, function(x){max(x$tx_s_years)})
  
  filter = simTestDs.id$nObs_fixed <= minObsArr[i]
  simTestDs.id[filter, paste("nObs_",i, sep="")] = simTestDs.id$nObs_fixed[filter]
  simTestDs.id[filter, paste("stopTime_",i, sep="")] = simTestDs.id$stopTime_fixed[filter]
  
  i = i + 1
}

simTestDs.id_long=reshape(simTestDs.id, direction='long', idvar='amctx', timevar = "methodNumber",
                   varying=list(seq(17, ncol(simTestDs.id), 2), seq(18, ncol(simTestDs.id), 2)),
                   v.names=c('nObs', 'stopTime'))
simTestDs.id_long = simTestDs.id_long[order(simTestDs.id_long$amctx, simTestDs.id_long$methodNumber, na.last = T), ]
simTestDs.id_long$offset = simTestDs.id_long$stopTime - simTestDs.id_long$stoptime_True
simTestDs.id_long$offset_fail = simTestDs.id_long$stopTime - simTestDs.id_long$years_tx_gl

temp = simTestDs.id_long[simTestDs.id_long$methodNumber %in% c(1,3,4),]
#temp = simTestDs.id_long[simTestDs.id_long$methodNumber <= 2,]
temp$methodNumber = factor(temp$methodNumber, labels =  c("Fixed", "Personalized_1", "Personalized_2"))
#temp$methodNumber = factor(temp$methodNumber, labels =  c("Fixed", "Personalized"))
a = ggplot(data=temp) + 
  geom_boxplot(aes(reorder(methodNumber, offset, FUN=median), offset_fail)) + 
  xlab("Method") + ylab("Stop time - True fail time; [years]") + coord_flip()

b = ggplot(data=temp) + 
  geom_boxplot(aes(reorder(methodNumber, offset, FUN=median), offset)) + 
  xlab("Method") + ylab("Stop time - True threshold time; [years]") + coord_flip()

c = ggplot(data=temp) + 
  geom_boxplot(aes(reorder(methodNumber, nObs, FUN=median), nObs)) + 
  xlab("Method") + ylab("Number of observations") + coord_flip()

multiplot(a,b,c)

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
