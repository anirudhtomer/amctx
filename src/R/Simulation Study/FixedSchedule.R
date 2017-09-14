cl = makeCluster(detectCores())
registerDoParallel(cl)

testDs.id$nObs_fixed = NA
testDs.id$stopTime_fixed = NA

fixed_results = foreach(id=testIdOfInterest[1:30], .packages = c("splines", "JMbayes"), .combine="rbind") %dopar%{
  
  ds_i = simDs[simDs$amctx == id, ]
  wGamma_i = wGamma[id]
  b_i = b_creatinine[id, ]
  
  nObs_fixed = 0
  stopTime_fixed = NA
  for(j in 1:nrow(ds_i)){
    ds_i_j = ds_i[1:j, ]
    
    temp = survfitJM(simJointModel_replaced, ds_i_j, idVar="amctx", survTimes = ds_i$tx_s_years[j] + maxRiskDt)
    
    dynSurvProb = temp$summaries[[1]][1, "Median"]
    
    nObs_fixed = j
    if(dynSurvProb <= minSurv){
      stopTime_fixed = ds_i$tx_s_years[j]
      break
    }
  }
  
  return(c(nObs_fixed = nObs_fixed, stopTime_fixed = stopTime_fixed))
}

testDs.id$nObs_fixed[testIdOfInterest[1:30] - trainingSize] = fixed_results[,"nObs_fixed"]
testDs.id$stopTime_fixed[testIdOfInterest[1:30] - trainingSize] = fixed_results[,"stopTime_fixed"]

stopCluster(cl)
