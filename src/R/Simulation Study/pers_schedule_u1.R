source("src/R/Simulation Study/utilFunctions/dynInfoPar.R")
source("src/R/Simulation Study/utilFunctions/dynInfo_mod.R")

invDynSurvival <- function (t, u, patientDs, lasttime) {
  u - round(survfitJM(simJointModel_replaced, patientDs, idVar="amctx", last.time = lasttime, survTimes = t)$summaries[[1]][1, "Median"],4)
}

pDynSurvTime = function(survProb, patientDs, lasttime=NULL){
  
  Low = max(patientDs$tx_s_years) + 1e-05
  if(!is.null(lasttime)){
    Low = lasttime + 1e-05
  }
  Up <- 20
  tries  = 0
  
  repeat{
    tries = tries + 1
    Root <- try(uniroot(invDynSurvival, interval = c(Low, Up), 
                        u = survProb, patientDs = patientDs, lasttime = lasttime)$root, TRUE)
    
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

#Technique 1:
#Use dynamic prediction to choose upper limit
#Use EKL without denominator for choosing the time point

#Technique 1:
#Use dynamic prediction to choose upper limit

ct = makeCluster(detectCores())
registerDoParallel(ct)

for(minFixedMeasurements in c(8)){
  
  for(methodName in c("dynInfo_mod")){
    #for(methodName in c("dynInfo_mod")){
    
    dynInfoMethod = get(methodName)
    
    persTestDs = simDs[simDs$visitNumber <= minFixedMeasurements & 
                         simDs$amctx %in% testIdOfInterest,]
    persTestDs$creatinine = exp(persTestDs$logCreatinine)
    patientDsList = split(persTestDs, persTestDs$amctx)
    
    for(i in 1:length(patientDsList)){
      patientId = patientDsList[[i]]$amctx[1]
      trueStopTime = simTestDs.id$stoptime_True[testDs.id$amctx == patientId]
      
      print(paste(patientId, "---", trueStopTime))
      
      repeat{
        lastVisitTime = max(patientDsList[[i]]$tx_s_years)
        
        dynSurvProbDt = survfitJM(simJointModel_replaced, patientDsList[[i]], idVar="amctx", 
                                  survTimes = lastVisitTime + maxRiskDt)$summaries[[1]][1, "Median"]
        if(dynSurvProbDt <= minSurv | lastVisitTime > 17){
          break
        }
        
        maxInfoTime = pDynSurvTime(minSurv, patientDsList[[i]])
        maxInfoDt = maxInfoTime - lastVisitTime
        
        dynInfoRes = dynInfoMethod(simJointModel_replaced, newdata = patientDsList[[i]], Dt = maxInfoDt, K = 50, seed = 4001, idVar="amctx")
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
    
    save(patientDsList, file = paste("Rdata/u1_", methodName, "_6mo_nFix_", minFixedMeasurements, "_risk_", maxRisk*100,"_k50.Rdata", sep=""))
  }
}
stopCluster(ct)
#nObs_u1 = sapply(patientDsList, nrow) - minFixedMeasurements
#stopTime_u1 = sapply(patientDsList, function(x){max(x$tx_s_years)})
