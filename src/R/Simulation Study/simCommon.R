plotTrueSurvival = function(patientId){
  
  time = 1:15
  survProb = sapply(time, function(t){
    survivalFunc(t, patientId)
  })
  
  byYAxis = (max(survProb) - min(survProb))/10
  
  qplot(x=time, y = survProb, geom = "line", xlab = "Time (years)", ylab = "Probability") + 
    ticksX(from=0, max = 15, by=1) + ticksY(min(survProb), max(survProb), by=byYAxis) + 
    ggtitle(patientId)
}

plotTrueLongitudinal = function(patientId){
  ggplot(data=simDs[simDs$amctx == patientId, ], aes(x=tx_s_years, y=logCreatinine)) + 
    geom_line() + geom_point(color="red") + ggtitle(patientId) + xlab("Time (years)") + 
    ylab("log(Creatinine)")
}

plotDynamicSurvival = function(patientId){
  ggplot(data=simTestDs[simTestDs$amctx==patientId,]) + 
    geom_line(aes(x=tx_s_years, y=fixed_pt5yr_survprob)) + xlab("Time (years)") + 
    ylab("Probability") + ggtitle(patientId)
}

simJointModel_replaced = replaceMCMCContents(mvJoint_creatinine_tdboth_training, jmbayes_creatinine_tdboth_training)

invDynSurvival <- function (t, u, patientDs, lasttime) {
  u - round(survfitJM(simJointModel_replaced, patientDs, idVar="amctx", last.time = lasttime, survTimes = t)$summaries[[1]][1, "Median"],3)
}

invDynSurvLastTime <- function (lasttime, u, patientDs, maxRiskDt) {
  u - round(survfitJM(simJointModel_replaced, patientDs, idVar="amctx", last.time = lasttime, survTimes = lasttime + maxRiskDt)$summaries[[1]][1, "Median"],3)
}

pDynSurvMaxRisk = function(survProb, patientDs, maxRiskDt){
  #Return the time at which the dynamic survival probability is say 90%
  
  Low = max(patientDs$tx_s_years) + 1e-05
  Up <- 10
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

pDynSurvTime = function(survProb, patientDs, lasttime=NULL){
  #Return the time at which the dynamic survival probability is say 90%
  
  Low = max(patientDs$tx_s_years) + 1e-05
  if(!is.null(lasttime)){
    Low = Low + lasttime
  }
  Up <- 10
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

rLogCreatinine =  function(patientId, time, mean=F){
  df_s = data.frame(tx_s_years = time, rec_age_fwp1 = simDs.id[patientId, "rec_age"], 
                    simDs.id[patientId, ])
  
  xi_s_val_creatinine = model.matrix(fixedValueFormula_creatinine, df_s)
  zi_s_val_creatinine = model.matrix(randomValueFormula_creatinine, df_s)
  zib_val_creatinine = zi_s_val_creatinine %*% b_creatinine[patientId, ]
  xBetaZb_s_value_creatinine = xi_s_val_creatinine %*% betas_creatinine + zib_val_creatinine
  
  if(mean==T){
    return(c(xBetaZb_s_value_creatinine))
  }else{
    return(sapply(xBetaZb_s_value_creatinine, rnorm, n=1, sigma.y_creatinine)) 
  }
}