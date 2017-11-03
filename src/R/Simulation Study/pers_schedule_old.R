source("src/R/Simulation Study/dynInfo_mod.R")
source("src/R/Simulation Study/dynInfoPar.R")
source("src/R/Simulation Study/pDynStopTime.R")

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

