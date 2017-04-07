expectedCondFailureTime = function(patientDs){
  
  dynamicPredProb = function(futureTimes, patientDs){
    survfitJM(simJointModel_replaced, patientDs, idVar="amctx", survTimes = futureTimes)$summaries[[1]][, "Mean"]
  }
  
  lastVisitTime = max(patientDs$tx_s_years)
  lastVisitTime + integrate(dynamicPredProb, lastVisitTime, Inf, patientDs)$value
}