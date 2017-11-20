source("src/R/Simulation Study/utilFunctions/dynInfoPar.R")
source("src/R/Simulation Study/utilFunctions/dynInfo_mod.R")

invDynSurvLastTime <- function (lasttime, u, patientDs, maxRiskDt) {
  u - round(survfitJM(simJointModel_replaced, patientDs, idVar="amctx", last.time = lasttime, survTimes = lasttime + maxRiskDt)$summaries[[1]][1, "Median"], 4)
}

pDynSurvMaxRisk = function(survProb, patientDs, maxRiskDt){
  #Return the time at which the dynamic survival probability is say 90%
  
  Low = max(patientDs$tx_s_years) + 1e-05
  Up <- 20
  tries  = 0
  
  repeat{
    tries = tries + 1
    Root <- try(uniroot(invDynSurvLastTime, interval = c(Low, Up), 
                        u = survProb, patientDs = patientDs, maxRiskDt = maxRiskDt)$root, TRUE)
    
    if(inherits(Root, "try-error")){
      if(tries >= 40){
        return(NA)
      }else{
        Up = Up + 0.25    
      }
    }else{
      return(Root)
    }
  }
}


#Technique 3:
#Use own dynamic prediction to choose upper limit
#Use EKL without denominator for choosing the time point


ct = makeCluster(detectCores())
registerDoParallel(ct)

for(minFixedMeasurements in c(8)){
  
  for(methodName in c("dynInfoPar")){
    #for(methodName in c("dynInfo_mod")){
    
    dynInfoMethod = get(methodName)
    
    persTestDs = simDs[simDs$visitNumber <= minFixedMeasurements & 
                         simDs$amctx %in% simTestDs.id$amctx,]
    persTestDs$creatinine = exp(persTestDs$logCreatinine)
    patientDsList = split(persTestDs, persTestDs$amctx)
    
    for(i in 1:length(patientDsList)){
      patientId = patientDsList[[i]]$amctx[1]
      trueStopTime = simTestDs.id$stoptime_True[simTestDs.id$amctx == patientId]
      
      print(paste(patientId, "---", trueStopTime))
      
      if(!is.na(trueStopTime)){
        repeat{
          lastVisitTime = max(patientDsList[[i]]$tx_s_years)
          
          dynSurvProbDt = survfitJM(simJointModel_replaced, patientDsList[[i]], idVar="amctx", 
                                    survTimes = lastVisitTime + maxRiskDt)$summaries[[1]][1, "Median"]
          if(dynSurvProbDt <= minSurv | lastVisitTime > 15){
            break
          }
          
          maxInfoTime = pDynSurvMaxRisk(minSurv, patientDsList[[i]], maxRiskDt = maxRiskDt)
          maxInfoDt = maxInfoTime - max(patientDsList[[i]]$tx_s_years)
          
          dynInfoRes = dynInfoMethod(simJointModel_replaced, newdata = patientDsList[[i]], Dt = maxInfoDt, K = 25, seed = 4001, idVar="amctx")
          #info = exp(dynInfoRes$summary$Info)/apply(dynInfoRes$full.results,2, function(x){mad(exp(x))})
          info = dynInfoRes$summary$Info
          newTime = dynInfoRes$summary$times[which.max(info)]
          
          #add new row to the patient DS
          newRow = patientDsList[[i]][1, ]
          newRow$tx_s_years = newTime
          newRow$logCreatinine = rLogCreatinine(patientId, newTime)
          newRow$creatinine = exp(newRow$logCreatinine)
          
          patientDsList[[i]] = rbind(patientDsList[[i]], newRow)
          print(paste("Step", newTime))
        }
        print("Next Patient")
      }
    }
    save(patientDsList, file = paste("Rdata/u2_", methodName, "_k25.Rdata", sep=""))
  }
}
stopCluster(ct)
