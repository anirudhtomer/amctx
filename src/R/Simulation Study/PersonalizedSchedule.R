source("src/R/Simulation Study/dynInfo_mod.R")
source("src/R/Simulation Study/dynInfoPar.R")
source("src/R/Simulation Study/pDynStopTime.R")

avoid = c(10, 13)
minFixedMeasurements = 15

#Technique 2:
#Use dynamic prediction to choose upper limit
#Use EKL without denominator for choosing the time point
persTestDs = simTestDs[simTestDs$visitNumber <= minFixedMeasurements,]
patientDsListt2 = split(persTestDs, persTestDs$amctx)
for(i in 1:length(patientDsListt2)){
  patientDs_i = patientDsListt2[[i]]
  patientId = patientDs_i$amctx[1]
  print(paste(patientId, "---", simTestDs.id$dynamicMaxRiskTime[i]))
  
  if(!is.na(simTestDs.id$dynamicMaxRiskTime[simTestDs.id$amctx==patientId])){
    repeat{
      dynSurvProbDt = survfitJM(simJointModel_replaced, patientDsListt2[[i]], idVar="amctx", 
                                survTimes = max(patientDsListt2[[i]]$tx_s_years)+maxRiskDt)$summaries[[1]][1, "Mean"]
      if((1-dynSurvProbDt)>=maxRisk |  max(patientDsListt2[[i]]$tx_s_years)>10){
        break
      }
      
      maxInfoTime = pDynSurvTime(1-maxRisk, patientDsListt2[[i]])
      maxInfoDt = maxInfoTime - max(patientDsListt2[[i]]$tx_s_years)
      
      dynInfoRes = dynInfo_mod(simJointModel_replaced, newdata = patientDsListt2[[i]], Dt = maxInfoDt, K = 50, seed = 2017, idVar="amctx")
      info = dynInfoRes$summary$Info
      newTime = dynInfoRes$summary$times[which.max(info)]
      
      #add new row to the patient DS
      newRow = patientDs_i[1, ]
      newRow$tx_s_years = newTime
      newRow$logCreatinine = rLogCreatinine(patientId = patientId, newRow$tx_s_years)
      
      patientDsListt2[[i]] = rbind(patientDsListt2[[i]], newRow)
      print(paste("Step", newTime))
    }
    print("Next Patient")
  }
}

stopCluster(ct2)
#
save(patientDsListt2, file = "Rdata/technique2.Rdata")

simTestDs.id$persScheduleObsCountT2 = sapply(patientDsListt2, nrow) - minFixedMeasurements
simTestDs.id$persScheduleStopTimeT2 = sapply(patientDsListt2, function(x){max(x$tx_s_years)})

result = data.frame(offset=numeric(), nObs=numeric())
result = rbind(result, cbind(offset=simTestDs.id$fixedScheduleStopTime[-avoid]-simTestDs.id$dynamicMaxRiskTime[-avoid],
                             nObs=simTestDs.id$fixedScheduleObsCount[-avoid]))
result = rbind(result, cbind(offset=simTestDs.id$persScheduleStopTimeT2[-avoid]-simTestDs.id$dynamicMaxRiskTime[-avoid],
                             nObs=simTestDs.id$persScheduleObsCountT2[-avoid]))
result = rbind(result, cbind(offset=simTestDs.id$persScheduleStopTimeT1[-avoid]-simTestDs.id$dynamicMaxRiskTime[-avoid],
                             nObs=simTestDs.id$persScheduleObsCountT1[-avoid]))
result$method = rep(c("Fixed", "InfoY", "InfoYandT"), each = nrow(simTestDs.id[-avoid,]))
result$offset = result$offset * 365
p1= ggplot(data=result) + geom_boxplot(aes(method, offset)) + scale_y_continuous(breaks = seq(-100, 400, 20)) + ylab("Offset(days)")
p2 = ggplot(data=result) + geom_boxplot(aes(method, nObs)) + scale_y_continuous(breaks = seq(0, 30, 1)) + ylab("Number of observations")
multiplot(plotlist = list(p1, p2), cols=2)


#Technique 3:
#Use own dynamic prediction to choose upper limit
#Use EKL without denominator for choosing the time point
ct3 = makeCluster(8)
registerDoParallel(ct3)

persTestDs = simTestDs[simTestDs$visitNumber <= minFixedMeasurements,]
patientDsListt3 = split(persTestDs, persTestDs$amctx)
for(i in 1:length(patientDsListt3)){
  patientDs_i = patientDsListt3[[i]]
  patientId = patientDs_i$amctx[1]
  print(paste(patientId, "---", simTestDs.id$dynamicMaxRiskTime[i]))
  
  if(!is.na(simTestDs.id$dynamicMaxRiskTime[simTestDs.id$amctx==patientId])){
    repeat{
      dynSurvProbDt = survfitJM(simJointModel_replaced, patientDsListt3[[i]], idVar="amctx", 
                                survTimes = max(patientDsListt3[[i]]$tx_s_years)+maxRiskDt)$summaries[[1]][1, "Mean"]
      if((1-dynSurvProbDt)>=maxRisk |  max(patientDsListt3[[i]]$tx_s_years)>10){
        break
      }
      
      maxInfoTime = pDynSurvMaxRisk(1-maxRisk, patientDsListt3[[i]], maxRiskDt = maxRiskDt)
      maxInfoDt = maxInfoTime - max(patientDsListt3[[i]]$tx_s_years)
      
      dynInfoRes = dynInfoPar(simJointModel_replaced, newdata = patientDsListt3[[i]], Dt = maxInfoDt, K = 50, seed = 2017, idVar="amctx")
      #info = dynInfoRes$summary$Info
      info = exp(dynInfoRes$summary$Info)/apply(dynInfoRes$full.results,2, function(x){mad(exp(x))})
      newTime = dynInfoRes$summary$times[which.max(info)]
      
      #add new row to the patient DS
      newRow = patientDs_i[1, ]
      newRow$tx_s_years = newTime
      newRow$logCreatinine = rLogCreatinine(patientId = patientId, newRow$tx_s_years)
      
      patientDsListt3[[i]] = rbind(patientDsListt3[[i]], newRow)
      print(paste("Step", newTime))
    }
    print("Next Patient")
  }
}

stopCluster(ct3)
save(patientDsListt3, file = "Rdata/technique3.Rdata")

simTestDs.id$persScheduleObsCountT3 = sapply(patientDsListt3, nrow) - minFixedMeasurements
simTestDs.id$persScheduleStopTimeT3 = sapply(patientDsListt3, function(x){max(x$tx_s_years)})

result = data.frame(offset=numeric(), nObs=numeric())
result = rbind(result, cbind(offset=simTestDs.id$fixedScheduleStopTime[-avoid]-simTestDs.id$dynamicMaxRiskTime[-avoid],
                             nObs=simTestDs.id$fixedScheduleObsCount[-avoid]))
result = rbind(result, cbind(offset=simTestDs.id$persScheduleStopTimeT3[-avoid]-simTestDs.id$dynamicMaxRiskTime[-avoid],
                             nObs=simTestDs.id$persScheduleObsCountT3[-avoid]))
result = rbind(result, cbind(offset=simTestDs.id$persScheduleStopTimeT2[-avoid]-simTestDs.id$dynamicMaxRiskTime[-avoid],
                             nObs=simTestDs.id$persScheduleObsCountT2[-avoid]))
result = rbind(result, cbind(offset=simTestDs.id$persScheduleStopTimeT1[-avoid]-simTestDs.id$dynamicMaxRiskTime[-avoid],
                             nObs=simTestDs.id$persScheduleObsCountT1[-avoid]))
result$method = rep(c("Fixed", "InfoYandTM2" , "InfoY", "InfoYandT"), each = nrow(simTestDs.id[-avoid,]))
result$offset = result$offset * 365
ggplot(data=result) + geom_boxplot(aes(method, offset)) + scale_y_continuous(breaks = seq(-100, 400, 30))
ggplot(data=result) + geom_boxplot(aes(method, nObs))


#Technique 4:
#Use own dynamic prediction to choose upper limit
#Use EKL without denominator for choosing the time point
ct4 = makeCluster(8)
registerDoParallel(ct4)

persTestDs = simTestDs[simTestDs$visitNumber <= minFixedMeasurements,]
patientDsListt4 = split(persTestDs, persTestDs$amctx)
for(i in 1:length(patientDsListt4)){
  patientDs_i = patientDsListt4[[i]]
  patientId = patientDs_i$amctx[1]
  print(paste(patientId, "---", simTestDs.id$dynamicMaxRiskTime[i]))
  
  if(!is.na(simTestDs.id$dynamicMaxRiskTime[simTestDs.id$amctx==patientId])){
    repeat{
      dynSurvProbDt = survfitJM(simJointModel_replaced, patientDsListt4[[i]], idVar="amctx", 
                                survTimes = max(patientDsListt4[[i]]$tx_s_years)+maxRiskDt)$summaries[[1]][1, "Mean"]
      if((1-dynSurvProbDt)>=maxRisk |  max(patientDsListt4[[i]]$tx_s_years)>10){
        break
      }
      
      maxInfoTime = pDynSurvMaxRisk(1-maxRisk, patientDsListt4[[i]], maxRiskDt = maxRiskDt)
      maxInfoDt = maxInfoTime - max(patientDsListt4[[i]]$tx_s_years)
      
      dynInfoRes = dynInfo_mod(simJointModel_replaced, newdata = patientDsListt4[[i]], Dt = maxInfoDt, K = 50, seed = 2017, idVar="amctx")
      #info = dynInfoRes$summary$Info
      info = exp(dynInfoRes$summary$Info)/apply(dynInfoRes$full.results,2, function(x){mad(exp(x))})
      newTime = dynInfoRes$summary$times[which.max(info)]
      
      #add new row to the patient DS
      newRow = patientDs_i[1, ]
      newRow$tx_s_years = newTime
      newRow$logCreatinine = rLogCreatinine(patientId = patientId, newRow$tx_s_years)
      
      patientDsListt4[[i]] = rbind(patientDsListt4[[i]], newRow)
      print(paste("Step", newTime))
    }
    print("Next Patient")
  }
}

stopCluster(ct4)
save(patientDsListt4, file = "Rdata/technique4.Rdata")

simTestDs.id$persScheduleObsCountT4 = sapply(patientDsListt4, nrow) - minFixedMeasurements
simTestDs.id$persScheduleStopTimeT4 = sapply(patientDsListt4, function(x){max(x$tx_s_years)})

result = data.frame(offset=numeric(), nObs=numeric())
result = rbind(result, cbind(offset=simTestDs.id$fixedScheduleStopTime[-avoid]-simTestDs.id$dynamicMaxRiskTime[-avoid],
                             nObs=simTestDs.id$fixedScheduleObsCount[-avoid]))
result = rbind(result, cbind(offset=simTestDs.id$persScheduleStopTimeT4[-avoid]-simTestDs.id$dynamicMaxRiskTime[-avoid],
                             nObs=simTestDs.id$persScheduleObsCountT4[-avoid]))
result = rbind(result, cbind(offset=simTestDs.id$persScheduleStopTimeT3[-avoid]-simTestDs.id$dynamicMaxRiskTime[-avoid],
                             nObs=simTestDs.id$persScheduleObsCountT3[-avoid]))
result = rbind(result, cbind(offset=simTestDs.id$persScheduleStopTimeT2[-avoid]-simTestDs.id$dynamicMaxRiskTime[-avoid],
                             nObs=simTestDs.id$persScheduleObsCountT2[-avoid]))
result$method = rep(c("Fixed", "U2", "U3" , "U1"), each = nrow(simTestDs.id[-avoid,]))
result$offset = result$offset * 365
ggplot(data=result) + geom_boxplot(aes(method, offset)) + scale_y_continuous(breaks = seq(-100, 600, 30))
ggplot(data=result) + geom_boxplot(aes(method, nObs)) + scale_y_continuous(breaks = seq(0, 30, 1))


#################################################
# OLD CODE
#################################################

# for(i in 1:length(patientDsList)){
#   patientDs_i = patientDsList[[i]]
#   patientId = patientDs_i$amctx[1]
#   print(paste(patientId, "---", simTestDs.id$dynamicMaxRiskTime[i]))
#   
#   if(!is.na(simTestDs.id$dynamicMaxRiskTime[simTestDs.id$amctx==patientId])){
#     repeat{
#       dynSurvProbDt = survfitJM(simJointModel_replaced, patientDsList[[i]], idVar="amctx", 
#                                 survTimes = max(patientDsList[[i]]$tx_s_years)+maxRiskDt)$summaries[[1]][1, "Mean"]
#       if((1-dynSurvProbDt)>=maxRisk |  max(patientDsList[[i]]$tx_s_years)>10){
#         break
#       }
#       
#       #Subtracting 2, because we are interested in 2 years survival time
#       #temp <- pDynStopTime(simJointModel_replaced, newdata = patientDsList[[i]], Dt = 2, K = 50, seed = 2015, idVar="amctx", maxRiskTime = maxRiskDt)
#       #maxInfoTime = temp$summary$times[which(exp(temp$summary$Info) <= (1-maxRisk))[1]]
#       #apply(temp$full.results, 2, function(x){exp(HPDinterval(as.mcmc(x)))})
#       
#       #maxInfoTime = temp$summary$times[apply(temp$full.results, 2, function(x){interval = exp(HPDinterval(as.mcmc(x))); (1-maxRisk) > interval[1] & (1-maxRisk) <= interval[2]})][1]
#       
#       #if(is.na(maxInfoTime)){
#       #maxInfoTime = temp$summary$times[1]
#       #}
#       
#       #maxInfoDt = maxInfoTime - max(patientDsList[[i]]$tx_s_years)
#       
#       temp2 = dynInfoPar(simJointModel_replaced, newdata = patientDsList[[i]], Dt = 5, K = 50, seed = 2015, idVar="amctx")
#       #info = temp2$summary$Info
#       #apply(exp(temp2$full.results), 2, median)/apply(exp(temp2$full.results),2, sd)
#       #info = exp(temp$summary$Info)/apply(temp2$full.results,2, function(x){tt = HPDinterval(as.mcmc(exp(x))); tt[2]-tt[1]})
#       info = sapply(1:ncol(temp2$full.results), function(colnum){
#         info_t = exp(temp2$full.results[, colnum])
#         interval_90 = HPDinterval(as.mcmc(info_t), prob=0.9)
#         info_filtered = info_t[info_t >= interval_90[1] & info_t <= interval_90[2]]
#         mean(info_filtered)/sd(info_filtered)
#       })
#       
#       newTime = temp2$summary$times[which.max(info)]
#       
#       #add new row to the patient DS
#       newRow = patientDs_i[1, ]
#       newRow$tx_s_years = newTime
#       newRow$logCreatinine = rLogCreatinine(patientId = patientId, newRow$tx_s_years)
#       
#       patientDsList[[i]] = rbind(patientDsList[[i]], newRow)
#       print(paste("Step", newTime))
#     }
#     print("Next Patient")
#   }
# }

