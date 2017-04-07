#Step 1: Calculate after first 5 (say) meaurements, the time point at which we have maxRisk probability

#Should this be after 2 months?
minFixedMeasurements = 5

persTestDs = simTestDs[simTestDs$visitNumber <= minFixedMeasurements,]
patientDsList = split(persTestDs, persTestDs$amctx)

ND = patientDsList[[1]]
longprof = predict(simJointModel_replaced,ND , type = "Subject", interval = "confidence", 
                   return = TRUE, idVar="amctx", FtTimes = seq(0.5, 6, 0.1))
last.time <- with(longprof, tx_s_years[!is.na(low)][1])
longprof[longprof$tx_s_years>last.time,]$logCreatinine=NA
ggplot(data = longprof, aes(x = tx_s_years, y=pred)) + geom_line() +
  geom_ribbon(aes(ymin=low, ymax=upp), fill="grey", alpha=0.5) +
  geom_point(aes(y=logCreatinine), colour="red", alpha=0.4) +
  geom_vline(xintercept = last.time, linetype="dotted") +
  xlab("Time (years)") + ylab("Predicted log(serum creatinine)")
  

for(i in 1:length(patientDsList)){
  patientDs_i = patientDsList[[i]]
  patientId = patientDs_i$amctx[1]
  
  if(!is.na(simTestDs.id$dynamicMaxRiskTime[simTestDs.id$amctx==patientId])){
    repeat{
      dynSurvProbDt = survfitJM(simJointModel_replaced, patientDsList[[i]], idVar="amctx", 
                   survTimes = max(patientDsList[[i]]$tx_s_years)+maxRiskDt)$summaries[[1]][1, "Mean"]
      if((1-dynSurvProbDt)>=maxRisk |  max(patientDsList[[i]]$tx_s_years)>10){
        break
      }
      
      #Subtracting 2, because we are interested in 2 years survival time
      temp <- pDynStopTime(simJointModel_replaced, newdata = patientDsList[[i]], Dt = 2, K = 50, seed = 2015, idVar="amctx", maxRiskTime = maxRiskDt)
      #maxInfoTime = temp$summary$times[which(exp(temp$summary$Info) <= (1-maxRisk))[1]]
      #apply(temp$full.results, 2, function(x){exp(HPDinterval(as.mcmc(x)))})
      
      maxInfoTime = temp$summary$times[apply(temp$full.results, 2, function(x){interval = exp(HPDinterval(as.mcmc(x))); (1-maxRisk) > interval[1] & (1-maxRisk) <= interval[2]})][1]
      
      if(is.na(maxInfoTime)){
        maxInfoTime = temp$summary$times[1]
      }
      
      maxInfoDt = maxInfoTime - max(patientDsList[[i]]$tx_s_years)
      
      temp2 = dynInfoPar(simJointModel_replaced, newdata = patientDsList[[i]], Dt = 5, K = 50, seed = 2015, idVar="amctx")
      newTime = temp2$summary$times[which.max(temp2$summary$Info)]
      
      #add new row to the patient DS
      newRow = patientDs_i[1, ]
      newRow$tx_s_years = newTime
      newRow$logCreatinine = rLogCreatinine(patientId = patientId, newRow$tx_s_years)
      
      patientDsList[[i]] = rbind(patientDsList[[i]], newRow)
      print("Step")
    }
  print(paste("Patient", i))
 }
}

save.image("Rdata/simCreatinine1.Rdata")

simTestDs.id$persScheduleStopTime = sapply(patientDsList, function(x){max(x$tx_s_years)})


#Should this be after 2 months?
minFixedMeasurements = 5

persTestDs = simTestDs[simTestDs$visitNumber <= minFixedMeasurements,]
patientDsList = split(persTestDs, persTestDs$amctx)

for(i in 1:length(patientDsList)){
  patientDs_i = patientDsList[[i]]
  patientId = patientDs_i$amctx[1]
  
  if(!is.na(simTestDs.id$dynamicMaxRiskTime[simTestDs.id$amctx==patientId])){
    repeat{
      dynSurvProbDt = survfitJM(simJointModel_replaced, patientDsList[[i]], idVar="amctx", 
                                survTimes = max(patientDsList[[i]]$tx_s_years)+maxRiskDt)$summaries[[1]][1, "Mean"]
      if((1-dynSurvProbDt)>=maxRisk |  max(patientDsList[[i]]$tx_s_years)>10){
        break
      }
      
      #Subtracting 2, because we are interested in 2 years survival time
      temp <- pDynStopTime(simJointModel_replaced, newdata = patientDsList[[i]], Dt = 8, K = 50, seed = 2015, idVar="amctx", maxRiskTime = maxRiskDt)
      maxInfoTime = temp$summary$times[which(exp(temp$summary$Info) <= (1-maxRisk))[1]]
      
      if(is.na(maxInfoTime)){
        maxInfoTime = temp$summary$times[1]
      }
      
      maxInfoDt = maxInfoTime - max(patientDsList[[i]]$tx_s_years)
      
      temp2 = dynInfo_mod(simJointModel_replaced, newdata = patientDsList[[i]], Dt = maxInfoDt, K = 50, seed = 2015, idVar="amctx")
      newTime = temp2$summary$times[which.max(temp2$summary$Info)]
      
      #add new row to the patient DS
      newRow = patientDs_i[1, ]
      newRow$tx_s_years = newTime
      newRow$logCreatinine = rLogCreatinine(patientId = patientId, newRow$tx_s_years)
      
      patientDsList[[i]] = rbind(patientDsList[[i]], newRow)
    }
    print(paste("Patient", i))
  }
}
save.image("Rdata/simCreatinine2.Rdata")

######################################################################
#### Taking more measurements in between
######################################################################
#Step 1: sample b|T>t, Y(t), theta1


#Step 2: sample T|T>t