source("src/R/Simulation Study/utilFunctions/dynInfoPar.R")

#Technique 1:
#Use dynamic prediction to choose upper limit
#Use EKL without denominator for choosing the time point

minFixedMeasurements = 10

ct1 = makeCluster(detectCores())
registerDoParallel(ct1)

persTestDs = simDs[simDs$visitNumber <= minFixedMeasurements & 
                     simDs$amctx %in% testIdOfInterest[1:1],]
patientDsList = split(persTestDs, persTestDs$amctx)

for(i in 1:length(patientDsList)){
  patientId = patientDsList[[i]]$amctx[1]
  trueStopTime = testDs.id$stoptime_True[testDs.id$amctx == patientId]
  
  print(paste(patientId, "---", trueStopTime))
  
  if(!is.na(patientId)){
    repeat{
      lastVisitTime = max(patientDsList[[i]]$tx_s_years)
      
      dynSurvProbDt = survfitJM(simJointModel_replaced, patientDsList[[i]], idVar="amctx", 
                                survTimes = lastVisitTime + maxRiskDt)$summaries[[1]][1, "Median"]
      if(dynSurvProbDt <= minSurv | lastVisitTime > 10){
        break
      }
      
      maxInfoTime = pDynSurvTime(minSurv, patientDsList[[i]])
      maxInfoDt = maxInfoTime - lastVisitTime
      
      dynInfoRes = dynInfoPar(simJointModel_replaced, newdata = patientDsList[[i]], Dt = maxInfoDt, K = 100, seed = 2017, idVar="amctx")
      info = dynInfoRes$summary$Info
      newTime = dynInfoRes$summary$times[which.max(info)]
      
      #add new row to the patient DS
      newRow = patientDsList[[i]]
      newRow$tx_s_years = newTime
      newRow$logCreatinine = rLogCreatinine(patientId, newTime)
      
      patientDsList[[i]] = rbind(patientDsList[[i]], newRow)
      print(paste("Step", newTime))
    }
  }
  print("Next Patient")
}
stopCluster(ct1)

save(patientDsList, file = "Rdata/u1_k100.Rdata")

nObs_u1 = sapply(patientDsList, nrow) - minFixedMeasurements
stopTime_u1 = sapply(patientDsList, function(x){max(x$tx_s_years)})
