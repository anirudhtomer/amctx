cl = makeCluster(detectCores())
registerDoParallel(cl)

simTestDs.id$nObs_fixed = NA
simTestDs.id$stopTime_fixed = NA

fixed_results = foreach(id=testIdOfInterest, .packages = c("splines", "JMbayes"), .combine="rbind") %dopar%{
  set.seed(sample(1:10000, size=1))
  ds_i = simDs[simDs$amctx == id, ]
  wGamma_i = wGamma[id]
  b_i = b_creatinine[id, ]
  
  nObs_fixed = 0
  stopTime_fixed = NA
  for(j in 1:nrow(ds_i)){
    ds_i_j = ds_i[1:j, ]
    
    temp = survfitJM(simJointModel_replaced, ds_i_j, idVar="amctx", survTimes = ds_i$tx_s_years[j] + maxRiskDt)
    
    dynSurvProb = temp$summaries[[1]][1, "Median"]
    print(dynSurvProb)
    
    nObs_fixed = j
    if(dynSurvProb <= minSurv){
      stopTime_fixed = ds_i$tx_s_years[j]
      break
    }
  }
  
  return(c(nObs_fixed = nObs_fixed, stopTime_fixed = stopTime_fixed))
}

simTestDs.id$nObs_fixed = fixed_results[,"nObs_fixed"]
simTestDs.id$stopTime_fixed = fixed_results[,"stopTime_fixed"]

stopCluster(cl)
save(simTestDs.id, file="Rdata/creatinine_sim_6month_7.5_percentrisk.Rdata")
