#Step 1: Calculate after first 5 (say) meaurements, the time point at which 
# we have maxRisk probability

#The following approach is wrong and doesn't work

#Should this be after 2 months?
minFixedMeasurements = 10

persTestDs = simTestDs[simTestDs$visitNumber <= minFixedMeasurements,]
patientDsList = split(persTestDs, persTestDs$amctx)
rm(persTestDs)

for(i in 1:1){
  patientId = patientDsList[[i]]$amctx[1]
  
  repeat{
    newTime <- pDynSurvTime(1-maxRisk, patientDsList[[i]])
    
    #The idea is that P(T>t + maxRiskDt|T>t) = 1-maxRisk
    if(is.na(newTime) | (newTime - maxRiskDt)>10 | ((newTime - maxRiskDt) <= max(patientDsList[[i]]$tx_s_years))){
      break
    }
  
    newRow = patientDsList[[i]][1, ]
    newRow$tx_s_years = newTime - maxRiskDt
    newRow$logCreatinine = rLogCreatinine(patientId, newRow$tx_s_years)  
    patientDsList[[i]] = rbind(patientDsList[[i]], newRow)
  }
}