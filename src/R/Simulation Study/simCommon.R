plotTrueSurvival = function(patientId){
  
  time = 1:15
  survProb = sapply(time, function(t){
    survivalFunc(t, patientId, simDs.id[patientId, ], 
                 wGamma[patientId], b_creatinine[patientId, ])
  })
  
  byYAxis = (max(survProb) - min(survProb))/10
  
  qplot(x=time, y = survProb, geom = "line", xlab = "Time (years)", ylab = "Probability") + 
    ticksX(from=0, max = 15, by=1) + ticksY(min(survProb), max(survProb), by=byYAxis) + 
    ggtitle(patientId)
}

plotTrueLongitudinal = function(patientId){
  ds = simDs[simDs$amctx == patientId, ]
  ds$logCreatinine = rLogCreatinine(patientId, ds$tx_s_years, T)
  ggplot(data=ds, aes(x=tx_s_years, y=logCreatinine)) + 
    geom_line() + geom_point(color="red") + ggtitle(patientId) + xlab("Time (years)") + 
    ylab("log(Creatinine)")
}

plotObsLongitudinal = function(patientId){
  ggplot(data=simDs[simDs$amctx == patientId, ], aes(x=tx_s_years, y=logCreatinine)) + 
    geom_line() + geom_point(color="red") + ggtitle(patientId) + xlab("Time (years)") + 
    ylab("log(Creatinine)")
}

plotDynamicSurvival = function(patientId){
  ggplot(data=simTestDs[simTestDs$amctx==patientId,]) + 
    geom_line(aes(x=tx_s_years, y=fixed_pt5yr_survprob)) + xlab("Time (years)") + 
    ylab("Probability") + ggtitle(patientId)
}

generateLongtiudinalTimeBySchedule = function(){
  #20 times in the first year and then 4 times per year after that
  years = c(seq(0, 1, length.out = 20), seq(1.25, 10, by=1/4))
  
  return(years)
}


getBetasCreatinine = function(jointModel, weightedOnes = T){
  if(weightedOnes){
    jointModel$statistics$postwMeans$betas1
  }else{
    jointModel$statistics$postMeans$betas1
  }
}

getSigmaCreatinine = function(jointModel, weightedOnes = T){
  if(weightedOnes){
    jointModel$statistics$postwMeans$sigma1
  }else{
    jointModel$statistics$postMeans$sigma1
  }
}

getD = function(jointModel, weightedOnes = T){
  if(weightedOnes){
    jointModel$statistics$postwMeans$D
  }else{
    jointModel$statistics$postMeans$D
  }
}

getGamma = function(jointModel, weightedOnes = T){
  if(weightedOnes){
    jointModel$statistics$postwMeans$gammas
  }else{
    jointModel$statistics$postMeans$gammas
  }
}

getAlphaCreatinine = function(jointModel, weightedOnes = T){
  if(weightedOnes){
    jointModel$statistics$postwMeans$alphas
  }else{
    jointModel$statistics$postMeans$alphas
  }
}

getBsGammas = function(jointModel, weightedOnes = T){
  if(weightedOnes){
    jointModel$statistics$postwMeans$Bs_gammas
  }else{
    jointModel$statistics$postMeans$Bs_gammas
  }
}

hazardFunc = function (s, patientId, simDs.id_i, wGamma_i, b_i) {
  weibullShape_i = weibullShapes[patientId]
  weibullScale_i = weibullScales[patientId]
  
  baselinehazard_s = (weibullShape_i/weibullScale_i)*(s/weibullScale_i)^(weibullShape_i-1)
  
  df_s = data.frame(tx_s_years = s, simDs.id_i)
  
  xi_s_val_creatinine = model.matrix(fixedValueFormula_creatinine, df_s)
  xi_s_slope_creatinine = model.matrix(fixedSlopeFormula_creatinine, df_s)
  
  zi_s_val_creatinine = model.matrix(randomValueFormula_creatinine, df_s)
  zi_s_slope_creatinine = model.matrix(randomSlopeFormula_creatinine, df_s)
  
  zib_val_creatinine = zi_s_val_creatinine %*% b_i
  zib_slope_creatinine = zi_s_slope_creatinine %*% b_i[-1] #-1 to ignore intercept
  
  xBetaZb_s_value_creatinine = xi_s_val_creatinine %*% betas_creatinine + zib_val_creatinine
  xBetaZb_s_slope_creatinine = xi_s_slope_creatinine %*% betas_creatinine[18:21] + zib_slope_creatinine
  
  y_Alpha_creatinine = cbind(xBetaZb_s_value_creatinine, xBetaZb_s_slope_creatinine) %*% getAlphaCreatinine(fittedJointModel, weightedOnes = F)
    #c(1.9,0.6)
    #
  #y_Alpha_creatinine = cbind(xBetaZb_s_value_creatinine) %*% getAlphaCreatinine(fittedJointModel, weightedOnes = T)
  
  baselinehazard_s * exp(wGamma_i + y_Alpha_creatinine)
}


survivalFunc <- function (t, i, simDs.id_i, wGamma_i, b_i) {
  exp(-integrate(hazardFunc, lower = 0, upper = t, i, simDs.id_i, wGamma_i, b_i)$value)
}

invSurvival <- function (t, u, i, simDs.id_i, wGamma_i, b_i) {
  integrate(hazardFunc, lower = 0, upper = t, i, simDs.id_i, wGamma_i, b_i)$value + log(u)
}

pSurvTime = function(survProb, patientId, simDs.id_i, wGamma_i, b_i){
  Low = 1e-05
  Up <- 35
  tries  = 0
  
  #uniroot(invSurvival, interval = c(Low, Up), 
  #        u = survProb, i = patientId, simDs.id_i, wGamma_i, b_i)
  repeat{
    tries = tries + 1
    Root <- try(uniroot(invSurvival, interval = c(Low, Up), 
                        u = survProb, i = patientId, simDs.id_i, wGamma_i, b_i)$root, TRUE)
    
    if(inherits(Root, "try-error")){
      if(tries >= 5){
        return(NA)
      }else{
        Up = Up + 0.5    
      }
    }else{
      return(Root)
    }
  }
}


generateSimLongitudinalData = function(nSub){
  ###############################################
  # Design matrices for the longitudinal measurement model
  ###############################################
  #For the age variable check the simulated and observed data distribution
  longTimes = do.call(what = c, lapply(1:nSub, function(x){generateLongtiudinalTimeBySchedule()}))
  
  # qplot(x = c(amctx.id$rec_age, rgamma(nrow(amctx.id), shape=60, scale = 0.9)), 
  #       geom="density", color=c(rep("Obs", nrow(amctx.id)), rep("Sim", nrow(amctx.id))))
  # qplot(x = c(amctx.id$d_age, rgamma(nrow(amctx.id), shape=50, scale=1)), 
  #       geom="density", color=c(rep("Obs", nrow(amctx.id)), rep("Sim", nrow(amctx.id))))
  #qplot(x = c(amctx.id$tx_cit, rgamma(nrow(amctx.id), shape=4, scale=250)), 
  #       geom="density", color=c(rep("Obs", nrow(amctx.id)), rep("Sim", nrow(amctx.id))))
  # qplot(x = c(amctx.id$tx_dial_days, rgamma(nrow(amctx.id), shape=2.5, scale=500)), 
  #       geom="density", color=c(rep("Obs", nrow(amctx.id)), rep("Sim", nrow(amctx.id))))
  # qplot(x = c(amctx.id$tx_pra[amctx.id$tx_pra>0], 
  #             rgamma(nrow(amctx.id[amctx.id$tx_pra>0,]), shape=1, scale=12) + 2), 
  #       geom="density", color=c(rep("Obs", nrow(amctx.id[amctx.id$tx_pra>0,])), 
  #                               rep("Sim", nrow(amctx.id[amctx.id$tx_pra>0,]))))
  
  #make rec_age_fwp1 for the longitudinal data set
  rec_age_fwp1 = rgamma(n = nSub, shape=60, scale = 0.9)
  d_age = rgamma(n = nSub, shape=50, scale = 1)
  rec_bmi = rgamma(n = nSub, shape=50, scale = 0.5)
  d_bmi = rgamma(n = nSub, shape=50, scale = 0.5)
  tx_dial_days =  rgamma(nSub, shape=2.5, scale=500)
  tx_hla = rep(0:6, rmultinom(1, nSub, prob=table(amctx.id$tx_hla)/nrow(amctx.id)))
  tx_cit = rgamma(nSub, shape=4, scale=250)
  tx_previoustx = ifelse(rbinom(nSub, size = 1, table(amctx.id$tx_previoustx)["yes"]/nrow(amctx.id))==1, "yes", "no")
  d_gender = ifelse(rbinom(nSub, size = 1, table(amctx.id$d_gender)["M"]/nrow(amctx.id))==1, yes = "M", "F")
  tx_hvdis = ifelse(rbinom(nSub, size = 1, table(amctx.id$tx_hvdis)["no"]/nrow(amctx.id))==1, yes = "no", "yes")
  rec_gender = ifelse(rbinom(nSub, size = 1, table(amctx.id$rec_gender)["M"]/nrow(amctx.id))==1, yes = "M", "F")
  d_cadaveric = ifelse(rbinom(nSub, size = 1, table(amctx.id$d_cadaveric)["yes"]/nrow(amctx.id))==1, yes = "yes", "no")
  tx_dgf = ifelse(rbinom(nSub, size = 1, table(amctx.id$tx_dgf)["yes"]/nrow(amctx.id))==1, yes = "yes", "no")
  tx_dm = ifelse(rbinom(nSub, size = 1, table(amctx.id$tx_dm)["yes"]/nrow(amctx.id))==1, yes = "yes", "no")
  
  prob_ahnr = table(amctx.id$ah_nr)/nrow(amctx.id)
  ahnr_counts = rmultinom(1, size = nSub, prob_ahnr)
  ah_nr = sample(unlist(lapply(1:5, function(x){
    count = ahnr_counts[x, 1]
    return(rep(x-1, count))
  })), 
  replace = F)
  
  #tx_pra-first sample the zeros
  tx_pra_zeros = rep(0, sum(rbinom(nSub, size=1, table(amctx.id$tx_pra==0)["TRUE"]/nrow(amctx.id))))
  tx_pra = sample(c(tx_pra_zeros, round(rgamma(nSub-length(tx_pra_zeros), shape=1, scale=12))+2), replace = F)
  
  subId <- rep(1:nSub)
  simDs.id = data.frame(amctx = subId, rec_age_fwp1, d_age, rec_bmi, 
                        tx_dial_days, tx_previoustx, d_gender, rec_gender, 
                        d_cadaveric, tx_dgf, tx_dm, tx_pra, ah_nr, d_bmi, tx_hla, tx_hvdis, tx_cit)
  
  simDs = data.frame(simDs.id, tx_s_years = sort(longTimes))
  simDs = simDs[order(simDs$amctx, decreasing = F),]
  simDs$visitNumber = rep(1:timesPerSubject, nSub)
  
  X_creatinine = model.matrix(fixedValueFormula_creatinine, data = simDs)
  Z_creatinine = model.matrix(randomValueFormula_creatinine, data = simDs)
  
  b_creatinine <- mvrnorm(nSub, mu = rep(0, nrow(D)), D)
  
  Zb_creatinine = unlist(lapply(1:nSub,function(i){
    Z_creatinine[((i-1)*timesPerSubject + 1):(i*timesPerSubject),] %*% b_creatinine[i, ]
  }))
  
  xBetaZb_creatinine = X_creatinine %*% betas_creatinine + Zb_creatinine
  simDs$logCreatinine = rnorm(length(xBetaZb_creatinine), xBetaZb_creatinine, sigma.y_creatinine)
  
  W <- model.matrix(survivalFormula, data = simDs.id)[,-1] #drop the intercept
  wGamma <- as.vector(W %*% gammas)
  
  list(simDs = simDs, simDs.id = simDs.id, b=b_creatinine, wGamma=wGamma)
}


generateSimJointData = function(nSub, simDs, simDs.id, b, wGamma){
  u <- runif(nSub)
  
  simDs.id$years_tx_gl <- NA
  
  ct = makeCluster(detectCores())
  registerDoParallel(ct)
  
  #simDs.id$years_tx_gl = sapply(1:nSub, function(i){
  #  pSurvTime(u[i], i, simDs.id[i, ], wGamma[i], b[i, ])
  #})
  
  simDs.id$years_tx_gl = foreach(i=1:nSub,.combine='c', 
                                 .export=c("pSurvTime", "invSurvival", "hazardFunc", 
                                           "weibullShapes", "weibullScales", 
                                           "fixedValueFormula_creatinine", "fixedSlopeFormula_creatinine",
                                           "randomValueFormula_creatinine", "randomSlopeFormula_creatinine",
                                           "betas_creatinine", "getAlphaCreatinine", "fittedJointModel"),
                                 .packages = c("splines", "JMbayes")) %dopar%{
    pSurvTime(u[i], i, simDs.id[i, ], wGamma[i], b[i, ])
  }
  
  percentageRejected = sum(is.na(simDs.id$years_tx_gl[1:nSub]))/nSub
  
  stopCluster(ct)
  
  print(paste("Percent reject", percentageRejected*100))
  if(percentageRejected > 0.5){
     stop(paste("Too many NA's sampled ", percentageRejected*100, "%", sep=""))
  }
  
  pid_to_keep = simDs.id[!is.na(simDs.id$years_tx_gl),]$amctx
  
  simDs = simDs[simDs$amctx %in% pid_to_keep,]
  simDs.id = simDs.id[simDs.id$amctx %in% pid_to_keep,]
  b = b[pid_to_keep,]
  wGamma = wGamma[pid_to_keep]
  weibullShapes = weibullShapes[pid_to_keep]
  weibullScales = weibullScales[pid_to_keep]
  
  simDs.id$amctx = 1:nrow(simDs.id)
  simDs$amctx = rep(simDs.id$amctx, each=timesPerSubject)
  
  #Divide into training and test
  trainingSize = round(nrow(simDs.id)*0.75)
  trainingDs.id = simDs.id[1:trainingSize, ]
  testDs.id = simDs.id[(trainingSize+1):nrow(simDs.id),]
  trainingDs = simDs[simDs$amctx %in% trainingDs.id$amctx, ]
  testDs = simDs[simDs$amctx %in% testDs.id$amctx, ]
  
  # simulate censoring times from an exponential distribution for TRAINING DATA SET ONLY
  # and calculate the observed event times, i.e., min(true event times, censoring times)
  
  #Ctimes <- rexp(trainingSize, 1/mean(prias.id[prias.id$progressed==0,]$progression_time))
  Ctimes = runif(trainingSize, 0, 15)
  
  trainingDs.id$gl_failure = trainingDs.id$years_tx_gl <= Ctimes
  trainingDs.id$years_tx_gl = pmin(trainingDs.id$years_tx_gl, Ctimes)
  trainingDs$years_tx_gl = rep(trainingDs.id$years_tx_gl, each=timesPerSubject)
  trainingDs$gl_failure = rep(trainingDs.id$gl_failure, each=timesPerSubject)
  testDs$years_tx_gl = rep(testDs.id$years_tx_gl, each=timesPerSubject)
  
  # drop the longitudinal measurementsthat were taken after the observed event time for each subject.
  trainingDs = trainingDs[trainingDs$tx_s_years <= trainingDs$years_tx_gl, ]
  
  ggplot(data=trainingDs, aes(y=logCreatinine, x=tx_s_years)) + geom_line(aes(group=amctx))
  
  list(simDs = simDs, simDs.id = simDs.id, 
       trainingDs = trainingDs, trainingDs.id = trainingDs.id, 
       testDs = testDs, testDs.id = testDs.id, 
       b=b, wGamma=wGamma, 
       weibullShapes = weibullShapes, weibullScales = weibullScales,
       percentageRejected = percentageRejected)
}


rLogCreatinine =  function(patientId, time, mean=F){
  df_s = data.frame(tx_s_years = time, simDs.id[patientId, ])
  
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

timesPerSubject = length(generateLongtiudinalTimeBySchedule())

#Formulae for simulation
fittedJointModel = mvJoint_creatinine_tdboth_complex

fixedValueFormula_creatinine = ~ 1 + rec_age_fwp1 + 
  rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi + 
  tx_hla + tx_pra + tx_dgf + ah_nr + tx_dm + tx_hvdis + 
  tx_cit + tx_dial_days + d_cadaveric +
  ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6))

randomValueFormula_creatinine = ~ 1 + ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6))

fixedSlopeFormula_creatinine = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6))

randomSlopeFormula_creatinine = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6))
#survivalFormula = ~ 1 + tx_hla + tx_previoustx + d_cadaveric + rec_bmi + tx_cit + tx_dial_days
survivalFormula = ~ 1 + tx_hla + tx_previoustx + tx_cit + tx_dial_days

gammas = getGamma(fittedJointModel, weightedOnes = F)
betas_creatinine <- getBetasCreatinine(fittedJointModel, weightedOnes = F)
sigma.y_creatinine <- getSigmaCreatinine(fittedJointModel, weightedOnes = F)
D <- getD(fittedJointModel, weightedOnes = F)