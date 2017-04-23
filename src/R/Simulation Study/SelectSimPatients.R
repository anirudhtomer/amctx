source("src/R/Simulation Study/simCommon.R")

testIdOfInterest = testDs.id$amctx[seq(1,15, by = 1)]

#Step 1: Select a bunch of patients
simTestDs.id = testDs.id[testDs.id$amctx %in% testIdOfInterest,]
simTestDs = testDs[testDs$amctx %in% simTestDs.id$amctx,]

#Plot their longitudinal and survival times
multiplot(plotlist = lapply(simTestDs.id$amctx, plotTrueLongitudinal), cols = 3)
multiplot(plotlist = lapply(simTestDs.id$amctx, plotTrueSurvival), cols = 3)

#The maximum risk for now is 20%. i.e. we want survival prob as 87.5% at 6 months
maxRisk = 0.125
maxRiskDt = 0.5

simTestDs.id$marginalMaxRiskTime = foreach(i=1:nrow(simTestDs.id),.combine='c', 
                                    .packages = c("splines", "JMbayes")) %do%{
  pSurvTime(1-maxRisk,  simTestDs.id[i, "amctx"])                                      
}

#We can split the entire time range into very small times and 
#check where exactly is the maxRisk with maxRiskDt not possible
testTimes = seq(0, 15, 0.01)
fineSimTestDs = data.frame(amctx = rep(simTestDs.id$amctx, each=length(testTimes)), 
                   rec_age_fwp1 = rep(simTestDs.id$rec_age, each=length(testTimes)),
                   rec_age = rep(simTestDs.id$rec_age, each=length(testTimes)),
                   d_age = rep(simTestDs.id$d_age, each=length(testTimes)),
                   rec_bmi = rep(simTestDs.id$rec_bmi, each=length(testTimes)),
                   tx_dial_days = rep(simTestDs.id$tx_dial_days, each=length(testTimes)),
                   tx_previoustx = rep(simTestDs.id$tx_previoustx, each=length(testTimes)),
                   d_gender = rep(simTestDs.id$d_gender, each=length(testTimes)),
                   rec_gender = rep(simTestDs.id$rec_gender, each=length(testTimes)),
                   d_cadaveric = rep(simTestDs.id$d_cadaveric, each=length(testTimes)),
                   tx_dgf = rep(simTestDs.id$tx_dgf, each=length(testTimes)),
                   tx_dm = rep(simTestDs.id$tx_dm, each=length(testTimes)),
                   ah_nr = rep(simTestDs.id$ah_nr, each=length(testTimes)),
                   tx_pra = rep(simTestDs.id$tx_pra, each=length(testTimes)),
                   visitNumber = rep(1:length(testTimes), nrow(simTestDs.id)),
                   tx_s_years = rep(testTimes, nrow(simTestDs.id)))

fineSimTestDs$logCreatinine = c(do.call(cbind, lapply(simTestDs.id$amctx, rLogCreatinine, testTimes, T)))

cl = makeCluster(8)
registerDoParallel(cl)
dynamicTrueSurvProb = foreach(time=testTimes, .packages = c("splines", "JMbayes")) %dopar%{
  patientDs = fineSimTestDs[fineSimTestDs$tx_s_years <= time,]                        
  survfitJM(simJointModel_replaced, patientDs, idVar="amctx", 
            survTimes = max(patientDs$tx_s_years)+maxRiskDt)
}
stopCluster(cl)

simTestDs.id$dynamicMaxRiskTime = NA
for(i in 1:length(testTimes)){
  riskProbs = 1 - sapply(1:nrow(simTestDs.id), function(index){
    dynamicTrueSurvProb[[i]]$summaries[[index]][1, "Mean"]
  })
  
  indexOfInterest = is.na(simTestDs.id$dynamicMaxRiskTime) & riskProbs >= maxRisk
  simTestDs.id$dynamicMaxRiskTime[indexOfInterest] = testTimes[i]
}

save.image("Rdata/simCreatinine.Rdata")
