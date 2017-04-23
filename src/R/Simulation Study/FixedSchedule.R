simTestDs$fixed_pt5yr_survprob = NA
for(j in 1:timesPerSubject){
  ND = simTestDs[simTestDs$visitNumber <= j, ]
  futureTimes = max(ND$tx_s_years) + maxRiskDt
  survivalProb = survfitJM(simJointModel_replaced, ND, idVar="amctx", survTimes = futureTimes)
  
  simTestDs$fixed_pt5yr_survprob[simTestDs$visitNumber==j] = sapply(1:nrow(simTestDs.id), function(i){
                                    survivalProb$summaries[[i]][1, "Mean"]
              }, simplify = T)
}

simTestDs.id$fixedScheduleStopTime = sapply(simTestDs.id$amctx, function(patientId){
  patientDs_i = simTestDs[simTestDs$amctx %in% patientId, ]
  
  stopVisitIndex = which((1 - patientDs_i$fixed_pt5yr_survprob) >= maxRisk)[1]
  if(!is.na(stopVisitIndex)){
    patientDs_i$tx_s_years[stopVisitIndex] 
  }else{
    NA
  }
})

simTestDs.id$fixedScheduleObsCount = sapply(simTestDs.id$amctx, function(patientId){
  patientDs_i = simTestDs[simTestDs$amctx %in% patientId, ]
  
  stopVisitIndex = which((1 - patientDs_i$fixed_pt5yr_survprob) >= maxRisk)[1]
  if(!is.na(stopVisitIndex)){
    stopVisitIndex - minFixedMeasurements
  }else{
    NA
  }
})

save.image("Rdata/simCreatinine.Rdata")

patientId = 583
multiplot(plotTrueLongitudinal(patientId), plotTrueSurvival(patientId), plotDynamicSurvival(patientId), cols=3)
