tx_s_years = seq(0, 6, by=0.005)

id = 606
temp_606 = cbind(simDs.id[id,], tx_s_years)

ds_i = simDs.id[id, ]
wGamma_i = wGamma[id]
b_i = b_creatinine[id, ]


temp_606$trueSurv = sapply(tx_s_years + maxRiskDt, survivalFunc, i=id, ds_i, wGamma_i, b_i) / sapply(tx_s_years, survivalFunc, i=id, ds_i, wGamma_i, b_i)

temp_606$creatinine = exp(rLogCreatinine(id, tx_s_years, mean = F))
temp_606$dynamicSurv = sapply(1:nrow(temp_606), function(timeIndex){
  survfitJM(simJointModel_replaced, temp_606[1:timeIndex,], idVar="amctx", 
                                           survTimes = temp_606$tx_s_years[timeIndex] + maxRiskDt)$summaries[[1]][1, "Median"]})

coarse1TimeIndices = seq(1, nrow(temp_606), by=3)
temp_606$dynamicSurvCoarse1 = NA
temp_606$dynamicSurvCoarse1[coarse1TimeIndices] = sapply(1:length(coarse1TimeIndices), function(timeIndex){
  ds = temp_606[coarse1TimeIndices,]
  survfitJM(simJointModel_replaced, ds[1:timeIndex,], idVar="amctx", 
            survTimes = ds$tx_s_years[timeIndex] + maxRiskDt)$summaries[[1]][1, "Median"]})

coarse2TimeIndices = seq(1, nrow(temp_606), by=6)
temp_606$dynamicSurvCoarse2 = NA
temp_606$dynamicSurvCoarse2[coarse2TimeIndices] = sapply(1:length(coarse2TimeIndices), function(timeIndex){
  ds = temp_606[coarse2TimeIndices,]
  survfitJM(simJointModel_replaced, ds[1:timeIndex,], idVar="amctx", 
            survTimes = ds$tx_s_years[timeIndex] + maxRiskDt)$summaries[[1]][1, "Median"]})

coarse3TimeIndices = seq(1, nrow(temp_606), by=10)
temp_606$dynamicSurvCoarse3 = NA
temp_606$dynamicSurvCoarse3[coarse3TimeIndices] = sapply(1:length(coarse3TimeIndices), function(timeIndex){
  ds = temp_606[coarse3TimeIndices,]
  survfitJM(simJointModel_replaced, ds[1:timeIndex,], idVar="amctx", seed = 10,
            survTimes = ds$tx_s_years[timeIndex] + maxRiskDt)$summaries[[1]][1, "Median"]})

coarse4TimeIndices = seq(1, nrow(temp_606), by=15)
temp_606$dynamicSurvCoarse4 = NA
temp_606$dynamicSurvCoarse4[coarse4TimeIndices] = sapply(1:length(coarse4TimeIndices), function(timeIndex){
  ds = temp_606[coarse4TimeIndices,]
  survfitJM(simJointModel_replaced, ds[1:timeIndex,], idVar="amctx", seed = 10,
            survTimes = ds$tx_s_years[timeIndex] + maxRiskDt)$summaries[[1]][1, "Median"]})

#3 months = 0.25 year gap
betweenBurstGap = round(0.5 * nrow(temp_606) / max(tx_s_years))
burstLength = 10
withinBurstGap = 1

burst1TimeIndices = c(1)
while(tail(burst1TimeIndices,1) < nrow(temp_606)){
  burst1TimeIndices = c(burst1TimeIndices, seq(tail(burst1TimeIndices,1)+1, 
                                               tail(burst1TimeIndices,1) + burstLength * withinBurstGap, 
                                               by=withinBurstGap))
  burst1TimeIndices = c(burst1TimeIndices, tail(burst1TimeIndices,1) + betweenBurstGap)
}

burst1TimeIndices = burst1TimeIndices[burst1TimeIndices <= nrow(temp_606)]

temp_606$dynamicSurvBurst1 = NA
temp_606$dynamicSurvBurst1[burst1TimeIndices] = sapply(1:length(burst1TimeIndices), function(timeIndex){
  ds = temp_606[burst1TimeIndices,]
  survfitJM(simJointModel_replaced, ds[1:timeIndex,], idVar="amctx", seed = 10,
            survTimes = ds$tx_s_years[timeIndex] + maxRiskDt)$summaries[[1]][1, "Median"]})


fixedIndices = sapply(generateLongtiudinalTimeBySchedule()[1:40], function(x){max(which(tx_s_years<=x))})
temp_606$dynamicSurvFixed = NA
temp_606$dynamicSurvFixed[fixedIndices] = sapply(1:length(fixedIndices), function(timeIndex){
  ds = temp_606[fixedIndices,]
  survfitJM(simJointModel_replaced, ds[1:timeIndex,], idVar="amctx", seed = 10,
            survTimes = ds$tx_s_years[timeIndex] + maxRiskDt)$summaries[[1]][1, "Median"]})

patientDsList[[paste(id)]]$dynamicSurv = sapply(1:nrow(patientDsList[[paste(id)]]), function(timeIndex){
  survfitJM(simJointModel_replaced, patientDsList[[paste(id)]][1:timeIndex,], idVar="amctx",seed=95,
            survTimes = patientDsList[[paste(id)]]$tx_s_years[timeIndex] + maxRiskDt)$summaries[[1]][1, "Median"]})

trueStopTime = simTestDs.id$stoptime_True[simTestDs.id$amctx==id]
xlim = trueStopTime + 0.5

plot_all = ggplot(NULL) + geom_line(data=temp_606, aes(x=tx_s_years, y=trueSurv, colour="True")) + 
  geom_line(data=temp_606,aes(x=tx_s_years, y=dynamicSurv, colour="Dynamic")) +
  geom_line(data=temp_606[!is.na(temp_606$dynamicSurvCoarse1),], aes(x=tx_s_years, y=dynamicSurvCoarse1, colour="Coarse1")) +
  geom_line(data=temp_606[!is.na(temp_606$dynamicSurvCoarse2),], aes(x=tx_s_years, y=dynamicSurvCoarse2, colour="Coarse2")) + 
  geom_line(data=temp_606[!is.na(temp_606$dynamicSurvCoarse3),], aes(x=tx_s_years, y=dynamicSurvCoarse3, colour="Coarse3")) + 
  geom_line(data=temp_606[!is.na(temp_606$dynamicSurvCoarse4),], aes(x=tx_s_years, y=dynamicSurvCoarse4, colour="Coarse4")) + 
  geom_line(data=temp_606[!is.na(temp_606$dynamicSurvBurst1),], aes(x=tx_s_years, y=dynamicSurvBurst1, colour="Burst1")) + 
  geom_line(data=temp_606[!is.na(temp_606$dynamicSurvFixed),], aes(x=tx_s_years, y=dynamicSurvFixed, colour="Fixed")) +
  geom_line(data=patientDsList[[paste(id)]], aes(x=tx_s_years, y=dynamicSurv, colour="Personalized")) + 
  geom_hline(yintercept = 0.95) + xlim(0, 6) + ylim(0.925, 1) + 
  xlab("Time (years)") + ylab("Conditional Survival Probability")

ggplotly(plot_all)
