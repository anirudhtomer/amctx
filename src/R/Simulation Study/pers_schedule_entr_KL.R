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

for(minFixedMeasurements in c(6)){
  
  for(methodName in c("dynInfo_mod")){
    #for(methodName in c("dynInfoPar")){
    
    dynInfoMethod = get(methodName)
    
    persTestDs = simDs[simDs$visitNumber <= minFixedMeasurements & 
                         simDs$amctx %in% simTestDs.id$amctx,]
    persTestDs$creatinine = exp(persTestDs$logCreatinine)
    patientDsList = split(persTestDs, persTestDs$amctx)
    
    for(i in 1:length(patientDsList)){
      patientId = patientDsList[[i]]$amctx[1]
      trueStopTime = simTestDs.id$stoptime_True[simTestDs.id$amctx == patientId]
      
      print(paste(patientId, "---", trueStopTime))
      
      repeat{
        lastVisitTime = max(patientDsList[[i]]$tx_s_years)
        
        if(lastVisitTime > 10){
          print("Last visit time > 10; breaking out")
        }
        
        dynSurvProbDt = survfitJM(simJointModel_replaced, patientDsList[[i]], idVar="amctx", 
                                  survTimes = lastVisitTime + maxRiskDt)$summaries[[1]][1, "Median"]
        
        if(dynSurvProbDt <= minSurv){
          checkFP = c(F,F,F)
          for(ttt in 1:3){
            newRow = patientDsList[[i]][1, ]
            newRow$tx_s_years = max(patientDsList[[i]]$tx_s_years) + 1/365
            newRow$logCreatinine = rLogCreatinine(newRow$amctx, newRow$tx_s_years)
            newRow$creatinine = exp(newRow$logCreatinine)
            patientDsList[[i]] = rbind(patientDsList[[i]], newRow)
            dynSurvProbDt = survfitJM(simJointModel_replaced, patientDsList[[i]], idVar="amctx", 
                                      survTimes = max(patientDsList[[i]]$tx_s_years) + maxRiskDt)$summaries[[1]][1, "Median"] 
            
            checkFP[ttt] = dynSurvProbDt <= minSurv
            if(checkFP[ttt] == F){
              print("False positive: cutoff broken")
              break
            }  
          }
          
          if(all(checkFP)){
            print("Cutoff broken 3 times")
            break
          }
        }
        
        maxInfoTime = pDynSurvTime(minSurv, patientDsList[[i]])
        #maxInfoTime = pDynSurvMaxRisk(minSurv, patientDsList[[i]], maxRiskDt = maxRiskDt)
        maxInfoDt = maxInfoTime - lastVisitTime
        
        lengthout = 16
        timesToPred = c()
        while(lengthout>2){
          timesToPred = seq(lastVisitTime, maxInfoTime, length.out = lengthout)[-1]
          if((timesToPred[2]-timesToPred[1])>=1/365){
            break
          }
          lengthout = lengthout/2
        }
        
        if(lengthout==2){
          print("Another layer of checking FP says truly below cutoff right now")
          break
        }
        
        dynInfoRes = dynInfoMethod(simJointModel_replaced, newdata = patientDsList[[i]], Dt = maxInfoDt, K =lengthout, seed = 4001, idVar="amctx")
        info = dynInfoRes$summary$Info
        
        ###########Technique 1a
        #info = exp(-dynInfoRes$summary$Info)/apply(dynInfoRes$full.results,2, function(x){mad(exp(-x))})
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
      save(patientDsList, file = paste("Rdata/u1_", methodName, "_6mo_nFix_", minFixedMeasurements, "_risk_", maxRisk*100,"_k15.Rdata", sep=""))
    }
  }
}
stopCluster(ct)
#nObs_u1 = sapply(patientDsList, nrow) - minFixedMeasurements
#stopTime_u1 = sapply(patientDsList, function(x){max(x$tx_s_years)})
