ct = makeCluster(detectCores())
registerDoParallel(ct)

for(minFixedMeasurements in c(2)){
  
  persTestDs = simDs[simDs$visitNumber <= minFixedMeasurements & 
                       simDs$amctx %in% simTestDs.id$amctx,]
  persTestDs$creatinine = exp(persTestDs$logCreatinine)
  patientDsList = split(persTestDs, persTestDs$amctx)
  
  for(i in 1:length(patientDsList)){
    patientId = patientDsList[[i]]$amctx[1]
    trueStopTime = simTestDs.id$stoptime_True[simTestDs.id$amctx == patientId]
    
    print(paste(patientId, "---", trueStopTime))
    
    #if(trueStopTime > generateLongtiudinalTimeBySchedule()[minFixedMeasurements]){
    if(T){
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
        
        lengthout = 16
        timesToPred = c()
        while(lengthout>2){
          timesToPred = seq(max(patientDsList[[i]]$tx_s_years), 
                            pDynSurvTime(minSurv, patientDsList[[i]]), length.out = lengthout)[-1]
          if((timesToPred[2]-timesToPred[1])>=1/365){
            break
          }
          lengthout = lengthout/2
        }
        
        if(lengthout==2){
          print("Another layer of checking FP says truly below cutoff right now")
          break
        }
        
        futureY = predict(simJointModel_replaced, newdata = patientDsList[[i]], 
                          FtTimes = timesToPred,
                          type = "Subject", interval = "confidence", idVar="amctx")
        
        probArr = foreach(k=1:length(timesToPred),.combine='c', 
                          .packages = c("splines", "JMbayes")) %dopar%{
                            
                            possibleY = futureY$all.vals[[1]][,k]
                            
                            dynSurvProbDtArr = sapply(possibleY, function(newy){
                              newRow = patientDsList[[i]][1, ]
                              newRow$tx_s_years = timesToPred[k]
                              newRow$logCreatinine = newy
                              newRow$creatinine = exp(newRow$logCreatinine)
                              patientDsList[[i]] = rbind(patientDsList[[i]], newRow)
                              
                              dynSurvProbDt = survfitJM(simJointModel_replaced, patientDsList[[i]], idVar="amctx", 
                                                        survTimes = timesToPred[k] + maxRiskDt)$summaries[[1]][1, "Median"]  
                            })
                            sum(dynSurvProbDtArr <= minSurv) / length(dynSurvProbDtArr)
                          }
        
        infoArr = probArr / apply(futureY$all.vals[[1]], 2, mad)
        
        newRow = patientDsList[[i]][1, ]
        newRow$tx_s_years = timesToPred[which.max(infoArr)]
        newRow$logCreatinine = rLogCreatinine(newRow$amctx, newRow$tx_s_years)
        newRow$creatinine = exp(newRow$logCreatinine)
        patientDsList[[i]] = rbind(patientDsList[[i]], newRow)
        
        print(paste("Step", newRow$tx_s_years))
      }
      print("Next Patient")
      save(patientDsList, file = paste("Rdata/u3_lessminObs.Rdata", sep=""))
    }else{
      print("Too early true stop time")
    }
  }
}
stopCluster(ct)
