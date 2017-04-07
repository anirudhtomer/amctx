source("src/R/common.R")

library(MASS)
library(splines)

set.seed(1000)
nSub <- 1000

fittedJointModel = mvJoint_creatinine_tdboth

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

hazardFunc = function (s, i) {
  pdf_s = dweibull(s, shape = weibullShape, scale = weibullScale)
  survival_s = (1-pweibull(q = s,shape= weibullShape, scale = weibullScale))
  baselinehazard_s = pdf_s/survival_s
  
  df_s = data.frame(tx_s_years = s, rec_age_fwp1 = simDs.id[i, "rec_age"], simDs.id[i, ])
  
  xi_s_val_creatinine = model.matrix(fixedValueFormula_creatinine, df_s)
  xi_s_slope_creatinine = model.matrix(fixedSlopeFormula_creatinine, df_s)
  
  zi_s_val_creatinine = model.matrix(randomValueFormula_creatinine, df_s)
  zi_s_slope_creatinine = model.matrix(randomSlopeFormula_creatinine, df_s)

  zib_val_creatinine = zi_s_val_creatinine %*% b_creatinine[i, ]
  zib_slope_creatinine = zi_s_slope_creatinine %*% b_creatinine[i, -1] #-1 to ignore intercept
  
  xBetaZb_s_value_creatinine = xi_s_val_creatinine %*% betas_creatinine + zib_val_creatinine
  xBetaZb_s_slope_creatinine = xi_s_slope_creatinine %*% betas_creatinine[c(8:11, 14:17, 18:21)] + zib_slope_creatinine
  
  y_Alpha_creatinine = cbind(xBetaZb_s_value_creatinine, xBetaZb_s_slope_creatinine) %*% getAlphaCreatinine(fittedJointModel)
  #y_Alpha_creatinine = cbind(xBetaZb_s_value_creatinine) %*% getAlphaCreatinine(fittedJointModel, weightedOnes = T)
  
  baselinehazard_s * exp(wGamma[i] + y_Alpha_creatinine)
}

survivalFunc <- function (t, i) {
  exp(-integrate(hazardFunc, lower = 0, upper = t, i)$value)
}

invSurvival <- function (t, u, i) {
  integrate(hazardFunc, lower = 0, upper = t, i)$value + log(u)
}

pSurvTime = function(survProb, patientId){
  Low = 1e-05
  Up <- 25
  tries  = 0
  
  repeat{
    tries = tries + 1
    Root <- try(uniroot(invSurvival, interval = c(Low, Up), 
                        u = survProb, i = patientId)$root, TRUE)
    
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

###############################################
# parameters for the linear mixed effects model
###############################################
fixedValueFormula_creatinine = ~ 1 + rec_age_fwp1 + rec_gender + d_age + 
  tx_pra + ah_nr + tx_dm + 
  ns(tx_s_years,knots=c(0.082, 0.192, 2.740), Boundary.knots = c(0, 6)) * d_cadaveric + 
  ns(tx_s_years,knots=c(0.082, 0.192, 2.740), Boundary.knots = c(0, 6)) * tx_dgf

randomValueFormula_creatinine = ~ 1 + ns(tx_s_years, knots = c(0.082, 0.192), Boundary.knots = c(0, 6))

fixedSlopeFormula_creatinine = ~ 0 + dns(tx_s_years,knots=c(0.082, 0.192, 2.740), Boundary.knots = c(0, 6)) + 
  I(dns(tx_s_years,knots=c(0.082, 0.192, 2.740), Boundary.knots = c(0, 6))*(as.numeric(d_cadaveric)-1)) + 
  I(dns(tx_s_years,knots=c(0.082, 0.192, 2.740), Boundary.knots = c(0, 6))*(as.numeric(tx_dgf)-1))

randomSlopeFormula_creatinine = ~ 0 + dns(tx_s_years, knots = c(0.082, 0.192), Boundary.knots = c(0, 6))

###############################################
# Design matrices for the longitudinal measurement model
###############################################
#For the age variable check the simulated and observed data distribution
longTimes = do.call(what = c, lapply(1:nSub, function(x){generateLongtiudinalTimeBySchedule()}))

timesPerSubject = length(longTimes)/nSub

# qplot(x = c(amctx.id$rec_age, rgamma(nrow(amctx.id), shape=60, scale = 0.9)), 
#       geom="density", color=c(rep("Obs", nrow(amctx.id)), rep("Sim", nrow(amctx.id))))
# qplot(x = c(amctx.id$d_age, rgamma(nrow(amctx.id), shape=50, scale=1)), 
#       geom="density", color=c(rep("Obs", nrow(amctx.id)), rep("Sim", nrow(amctx.id))))
# qplot(x = c(amctx.id$rec_bmi, rgamma(nrow(amctx.id), shape=50, scale=0.5)), 
#       geom="density", color=c(rep("Obs", nrow(amctx.id)), rep("Sim", nrow(amctx.id))))
# qplot(x = c(amctx.id$tx_dial_days, rgamma(nrow(amctx.id), shape=2.5, scale=500)), 
#       geom="density", color=c(rep("Obs", nrow(amctx.id)), rep("Sim", nrow(amctx.id))))
# qplot(x = c(amctx.id$tx_pra[amctx.id$tx_pra>0], 
#             rgamma(nrow(amctx.id[amctx.id$tx_pra>0,]), shape=1, scale=12) + 2), 
#       geom="density", color=c(rep("Obs", nrow(amctx.id[amctx.id$tx_pra>0,])), 
#                               rep("Sim", nrow(amctx.id[amctx.id$tx_pra>0,]))))

#make rec_age_fwp1 for the longitudinal data set
rec_age = rgamma(n = nSub, shape=60, scale = 0.9)
d_age = rgamma(n = nSub, shape=50, scale = 1)
rec_bmi = rgamma(n = nSub, shape=50, scale = 0.5)
tx_dial_days =  rgamma(nSub, shape=2.5, scale=500)
tx_previoustx = ifelse(rbinom(nSub, size = 1, table(amctx.id$tx_previoustx)["yes"]/nrow(amctx.id))==1, "yes", "no")
d_gender = ifelse(rbinom(nSub, size = 1, table(amctx.id$d_gender)["M"]/nrow(amctx.id))==1, yes = "M", "F")
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
simDs.id = data.frame(amctx = subId, rec_age_fwp1=rec_age, rec_age, d_age, rec_bmi, 
                      tx_dial_days, tx_previoustx, d_gender, rec_gender, 
                      d_cadaveric, tx_dgf, tx_dm, tx_pra, ah_nr)

simDs = data.frame(amctx = rep(subId, each=timesPerSubject), 
                   rec_age_fwp1 = rep(rec_age, each=timesPerSubject),
                   rec_age = rep(rec_age, each=timesPerSubject),
                   d_age = rep(d_age, each=timesPerSubject),
                   rec_bmi = rep(rec_bmi, each=timesPerSubject),
                   tx_dial_days = rep(tx_dial_days, each=timesPerSubject),
                   tx_previoustx = rep(tx_previoustx, each=timesPerSubject),
                   d_gender = rep(d_gender, each=timesPerSubject),
                   rec_gender = rep(rec_gender, each=timesPerSubject),
                   d_cadaveric = rep(d_cadaveric, each=timesPerSubject),
                   tx_dgf = rep(tx_dgf, each=timesPerSubject),
                   tx_dm = rep(tx_dm, each=timesPerSubject),
                   ah_nr = rep(ah_nr, each=timesPerSubject),
                   tx_pra = rep(tx_pra, each=timesPerSubject),
                   visitNumber = rep(1:timesPerSubject, nSub),
                   tx_s_years = longTimes, logCreatinine=NA)

X_creatinine = model.matrix(fixedValueFormula_creatinine, data = simDs)
Z_creatinine = model.matrix(randomValueFormula_creatinine, data = simDs)

betas_creatinine <- getBetasCreatinine(fittedJointModel, weightedOnes = T)
sigma.y_creatinine <- getSigmaCreatinine(fittedJointModel, weightedOnes = T)

D <- getD(fittedJointModel, weightedOnes = T)
b_creatinine <- mvrnorm(nSub, mu = rep(0, nrow(D)), D)

Zb_creatinine = unlist(lapply(1:nSub,function(i){
  Z_creatinine[((i-1)*timesPerSubject + 1):(i*timesPerSubject),] %*% b_creatinine[i, ]
}))

xBetaZb_creatinine = X_creatinine %*% betas_creatinine + Zb_creatinine
simDs$logCreatinine = rnorm(length(xBetaZb_creatinine), xBetaZb_creatinine, sigma.y_creatinine)

ggplot(data=simDs, aes(x=tx_s_years, y=logCreatinine)) + 
  geom_line(aes(group=amctx)) + stat_smooth()

# design matrix for the survival model
survivalFormula = ~ 1 + rec_age + d_age + tx_previoustx + d_gender + rec_bmi + tx_pra + I(tx_dial_days/365)
W <- model.matrix(survivalFormula, data = simDs.id)[,-1] #drop the intercept
gammas = getGamma(fittedJointModel, weightedOnes = T)
wGamma <- as.vector(W %*% gammas)

#weibullShape = 1.70
#weibullScale = 8000
#6250

# bsGammas = getBsGammas(fittedJointModel)
# time = seq(1, 10, by=0.1)
# baselineHazard_fitted = exp(splineDesign(fittedJointModel$control$knots, x = time) %*% bsGammas)
# 
# pdf_sim = dweibull(time, shape = weibullShape, scale = weibullScale)
# survival_sim = (1-pweibull(q = time,shape= weibullShape, scale = weibullScale))
# baselinehazard_sim = pdf_sim/survival_sim
# 
# p1 = qplot(y=c(baselinehazard_sim, baselineHazard_fitted), x = c(time,time), geom="line", 
#            color=c(rep("sim", length(time)), rep("Fitted", length(time)))) + theme(legend.position="none")
# p2 = qplot(x = rweibull(10000, shape = weibullShape, scale = weibullScale), geom="density")
# multiplot(p1, p2, cols=2)

#weibullShape = 1.40
#weibullScale = 80000

weibullShape = 1.40
weibullScale = 80000

#set.seed(432)
u <- runif(nSub)
simDs.id$years_tx_gl <- NA
simDs.id$years_tx_gl = foreach(i=1:nSub,.combine='c', 
                               .packages = c("splines", "JMbayes")) %do%{
  pSurvTime(u[i], i)
}
sum(is.na(simDs.id$years_tx_gl[1:nSub]))/nSub

failureTimeCompDs = data.frame(failuretime = c(simDs.id$years_tx_gl, amctx.id$years_tx_gl[amctx.id$gl_failure==1]),
                               type = c(rep("Sim",nrow(simDs.id)), rep("Obs", length(amctx.id$years_tx_gl[amctx.id$gl_failure==1]))))

ggplot(data = failureTimeCompDs) + geom_density(aes(x=failuretime, color=type, fill=type), alpha=0.3)

#Remove all patients with NA time, and keep a copy of the original DS
origSimDs = simDs
origSimDs.id = simDs.id
orig_b_creatinine = b_creatinine
orig_wGamma = wGamma

pid_to_keep = simDs.id[!is.na(simDs.id$years_tx_gl),]$amctx
simDs = simDs[simDs$amctx %in% pid_to_keep,]
simDs.id = simDs.id[simDs.id$amctx %in% pid_to_keep,]
b_creatinine = b_creatinine[pid_to_keep,]
wGamma = wGamma[pid_to_keep]
simDs.id$amctx = 1:nrow(simDs.id)
simDs$amctx = rep(simDs.id$amctx, each=timesPerSubject)

#Divide into training and test
trainingSize = max(round(nrow(simDs.id) * 0.75), nSub/2)
trainingDs.id = simDs.id[1:trainingSize, ]
testDs.id = simDs.id[(trainingSize + 1):nrow(simDs.id),]
trainingDs = simDs[simDs$amctx %in% trainingDs.id$amctx, ]
testDs = simDs[simDs$amctx %in% testDs.id$amctx, ]

save.image("Rdata/simCreatinine.Rdata")
# simulate censoring times from an exponential distribution for TRAINING DATA SET ONLY
# and calculate the observed event times, i.e., min(true event times, censoring times)
Ctimes <- rexp(trainingSize, 1/mean(amctx.id[amctx.id$gl_failure==0,]$years_tx_gl))
mean(Ctimes)
boxplot(Ctimes)
trainingDs.id$gl_failure = trainingDs.id$years_tx_gl <= Ctimes
trainingDs.id$years_tx_gl = pmin(trainingDs.id$years_tx_gl, Ctimes)
trainingDs$years_tx_gl = rep(trainingDs.id$years_tx_gl, each = timesPerSubject)

# drop the longitudinal measurementsthat were taken after the observed event time for each subject.
trainingDs = trainingDs[trainingDs$tx_s_years <= trainingDs$years_tx_gl, ]

ggplot(data=trainingDs, aes(x=tx_s_years, y=logCreatinine)) + 
  geom_line(aes(group=amctx)) + stat_smooth()

save.image("Rdata/simCreatinine.Rdata")

#########################################
# Fit joint model for the simulated data set
########################################
cox_Model_training = coxph(Surv(years_tx_gl, gl_failure) ~ rec_age + 
                             d_age + tx_previoustx + d_gender + 
                             rec_bmi + tx_pra + I(tx_dial_days/365),
                           data = trainingDs.id, x=T, model=T)

mvglmer_creatinine_training=mvglmer(list(logCreatinine ~ rec_age_fwp1 + 
                                           rec_gender + d_age + 
                                           tx_pra + ah_nr + tx_dm + 
                                           ns(tx_s_years,knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6)) * d_cadaveric + 
                                           ns(tx_s_years,knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6)) * tx_dgf + 
                                           (ns(tx_s_years, knots=c(30, 70)/365, Boundary.knots = c(0, 6))|amctx)),
                                    data = trainingDs, families = list(gaussian))

forms_creatinine_training <- list("logCreatinine" = "value",
                                  "logCreatinine" = list(fixed = ~ 0 + dns(tx_s_years, knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6)) + 
                                                           I(dns(tx_s_years, knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6))*(as.numeric(d_cadaveric)-1)) + 
                                                           I(dns(tx_s_years, knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6))*(as.numeric(tx_dgf)-1)),
                                                         random = ~ 0 + dns(tx_s_years, knots = c(30, 70)/365, Boundary.knots = c(0, 6)), 
                                                         indFixed = c(8:11, 14:17, 18:21), indRandom = 2:4, 
                                                         name = "slope"))

mvJoint_creatinine_tdval_training=mvJointModelBayes(mvglmer_creatinine_training, cox_Model_training, timeVar = "tx_s_years",
                                                    priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

mvJoint_creatinine_tdboth_training <- update(mvJoint_creatinine_tdval_training, Formulas = forms_creatinine_training,
                                             priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

###########################################
#Also fit joint model using the jointModelAPI
###########################################
lme_creatinine_training = lme(data=trainingDs,
                              fixed=logCreatinine ~ rec_age_fwp1 + 
                                rec_gender + d_age + 
                                tx_pra + ah_nr + tx_dm + 
                                ns(tx_s_years,knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6)) * d_cadaveric + 
                                ns(tx_s_years,knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6)) * tx_dgf,
                              random = ~ns(tx_s_years, knots=c(30, 70)/365, Boundary.knots = c(0, 6))|amctx,
                              control = lmeControl(opt = "optim"), method="ML")

jmbayes_creatinine_tdval_training = jointModelBayes(lmeObject = lme_creatinine_training, survObject = cox_Model_training, timeVar = "tx_s_years", control = list(n.iter=1000))
jmbayes_creatinine_tdboth_training = update(jmbayes_creatinine_tdval_training, param = "td-both", 
                                            extraForm = list(fixed = ~ 0 + dns(tx_s_years, knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6)) + 
                                                               I(dns(tx_s_years, knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6))*(as.numeric(d_cadaveric)-1)) + 
                                                               I(dns(tx_s_years, knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6))*(as.numeric(tx_dgf)-1)),
                                                             random = ~ 0 + dns(tx_s_years, knots = c(30, 70)/365, Boundary.knots = c(0, 6)), 
                                                             indFixed = c(8:11, 14:17, 18:21), indRandom = 2:4))

save.image("Rdata/simCreatinine.Rdata")

