cl = makeCluster(detectCores())
registerDoParallel(cl)

maxRiskDt = 0.5
maxRisk = 0.025
minSurv = 1 - maxRisk

fixed_results_pt025 = foreach(id=testDs.id$amctx, .packages = c("splines", "JMbayes"), 
                        .combine="rbind") %dopar%{
  set.seed(2018 + id)                        
  ds_i = simDs[simDs$amctx == id, ]
  wGamma_i = wGamma[id]
  b_i = b_creatinine[id, ]
   
  ds_i$creatinine = exp(ds_i$logCreatinine)
  
  latestRowPointer = 6
  while(ds_i$tx_s_years[latestRowPointer] + maxRiskDt <= 10){
    print(paste(id, "---", latestRowPointer))
    dynSurvProb = survfitJM(mvJoint_creatinine_tdboth_training, ds_i[1:latestRowPointer,], 
                     idVar="amctx", 
                     survTimes = ds_i$tx_s_years[latestRowPointer] + maxRiskDt)$summaries[[1]][1, "Mean"]
    if(dynSurvProb <= minSurv){
      #print(dynSurvProb)
      checkFP = c(F,F,F)
      for(ttt in 1:3){
        newRow = ds_i[1, ]
        newRow$tx_s_years = ds_i$tx_s_years[latestRowPointer] + 1/365
        newRow$logCreatinine = rLogCreatinine(newRow$amctx, newRow$tx_s_years)
        newRow$creatinine = exp(newRow$logCreatinine)
        ds_i = rbind(ds_i[1:latestRowPointer,], newRow, ds_i[latestRowPointer+1:nrow(ds_i),])
        
        latestRowPointer = latestRowPointer + 1
        
        dynSurvProbDt = survfitJM(mvJoint_creatinine_tdboth_training, ds_i[1:latestRowPointer,], idVar="amctx", 
                                  survTimes = ds_i$tx_s_years[latestRowPointer] + maxRiskDt)$summaries[[1]][1, "Mean"]
        
        checkFP[ttt] = dynSurvProbDt <= minSurv
        if(checkFP[ttt] == F){
          print("False positive: cutoff broken")
          break
        } 
      }
      
      if(all(checkFP)){
        break
      }
    }
    
    latestRowPointer = latestRowPointer + 1
  }
  
  return(c(nObs_fixed = latestRowPointer, stopTime_fixed = ds_i$tx_s_years[latestRowPointer]))
}

#testDs.id$nObs_fixed_new = fixed_results[,"nObs_fixed"]
#testDs.id$stopTime_fixed_new = fixed_results[,"stopTime_fixed"]

stopCluster(cl)
#save(testDs.id, file="Rdata/creatinine_sim_6month_2pt5_percentrisk.Rdata")
