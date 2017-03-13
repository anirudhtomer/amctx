library(MASS)
library(splines)

nSub <- 1000 # number of subjects

fittedJointModel = mvJoint_pcr_creatinine_tdboth

generateLongtiudinalTimeBySchedule = function(){
  #20 times in the first year and then 4 times per year after that
  years = c(seq(0, 1, length.out = 20), seq(1.25, 10, by=1/4))
  
  df = data.frame(tx_s_years = years, isPcr = T, isCreatinine = T)
  
  return(df)
}

getBetasPCR = function(jointModel, weightedOnes = T){
  if(weightedOnes){
    jointModel$statistics$postwMeans$betas1
  }else{
    jointModel$statistics$postMeans$betas1
  }
}

getBetasCreatinine = function(jointModel, weightedOnes = T){
  if(weightedOnes){
    jointModel$statistics$postwMeans$betas2
  }else{
    jointModel$statistics$postMeans$betas2 
  }
}

getSigmaPCR = function(jointModel, weightedOnes = T){
  if(weightedOnes){
    jointModel$statistics$postwMeans$sigma1
  }else{
    jointModel$statistics$postMeans$sigma1
  }
}

getSigmaCreatinine = function(jointModel, weightedOnes = T){
  if(weightedOnes){
    jointModel$statistics$postwMeans$sigma2
  }else{
    jointModel$statistics$postMeans$sigma2
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

getAlphaPcr = function(jointModel, weightedOnes = T){
  if(weightedOnes){
    jointModel$statistics$postwMeans$alphas[1:2]
  }else{
    jointModel$statistics$postMeans$alphas[1:2]
  }
}

getAlphaCreatinine = function(jointModel, weightedOnes = T){
  if(weightedOnes){
    jointModel$statistics$postwMeans$alphas[3:4]
  }else{
    jointModel$statistics$postMeans$alphas[3:4]
  }
}

hazardFunc = function (s, i) {
  pdf_s = dweibull(s, shape = weibullShape, scale = weibullScale)
  survival_s = (1-pweibull(q = s,shape= weibullShape, scale = weibullScale))
  baselinehazard_s = pdf_s/survival_s
  
  df_s = data.frame(tx_s_years = s, Age = Age[i])
  
  xi_s_val_pcr = model.matrix(fixedValueFormula_pcr, df_s)
  xi_s_val_creatinine = model.matrix(fixedValueFormula_creatinine, df_s)
  
  xi_s_slope_pcr = model.matrix(fixedSlopeFormula_pcr, df_s)
  xi_s_slope_creatinine = model.matrix(fixedSlopeFormula_creatinine, df_s)
  
  zi_s_val_pcr = model.matrix(randomValueFormula_pcr, df_s)
  zi_s_val_creatinine = model.matrix(randomValueFormula_creatinine, df_s)
  
  zi_s_slope_pcr = model.matrix(randomSlopeFormula_pcr, df_s)
  zi_s_slope_creatinine = model.matrix(randomSlopeFormula_creatinine, df_s)
  
  zib_val_pcr = zi_s_val_pcr %*% b_pcr[i, ]
  zib_val_creatinine = zi_s_val_creatinine %*% b_creatinine[i, ]
  
  zib_slope_pcr = zi_s_slope %*% b_pcr[i, -1] #-1 to ignore intercept
  zib_slope_creatinine = zi_s_slope %*% b_creatinine[i, -1] #-1 to ignore intercept
  
  xBetaZb_s_value_pcr = xi_s_val_pcr %*% betas_pcr + zib_val_pcr
  xBetaZb_s_value_creatinine = xi_s_val_creatinine %*% betas_creatinine + zib_val_creatinine
  
  #-c(1:3) to ignore intercept, age and age^2 in the
  xBetaZb_s_slope_pcr = xi_s_slope_pcr %*% betas_pcr[c(4:7, 9:12, 13:16)] + zib_slope_pcr
  xBetaZb_s_slope_creatinine = xi_s_slope_creatinine %*% betas_creatinine[c(8:11, 14:17, 18:21)] + zib_slope_creatinine
   
  y_Alpha_pcr = cbind(xBetaZb_s_value_pcr, xBetaZb_s_slope_pcr) %*% getAlphaPcr(fittedJointModel)
  y_Alpha_creatinine = cbind(xBetaZb_s_value_creatinine, xBetaZb_s_slope_creatinine) %*% getAlphaCreatinine(fittedJointModel)
  
  baselinehazard_s * exp(wGamma[i] + y_Alpha_pcr + y_Alpha_creatinine)
}

survivalFunc <- function (t, i) {
  exp(-integrate(hazardFunc, lower = 0, upper = t, i)$value)
}

invSurvival <- function (t, u, i) {
  integrate(hazardFunc, lower = 0, upper = t, i)$value + log(u)
}

###############################################
# parameters for the linear mixed effects model
###############################################
boundaryKnots_pcr <- c(0, 5.5)
boundaryKnots_creatinine <- c(0, 6)

fixedKnots_pcr <- c(0.082, 0.192, 2.740)
fixedKnots_creatinine <- c(0.082, 0.192, 2.740)

randomKnots_pcr <- fixedKnots_pcr[1:2]
randomKnots_creatinine <- fixedKnots_creatinine[1:2]

fixedValueFormula_pcr = ~ 1 + rec_gender + d_age + 
  ns(tx_s_years, knots=fixedKnots_pcr, Boundary.knots = boundaryKnots_pcr) * rec_gender + 
  ns(tx_s_years,knots=fixedKnots_pcr, Boundary.knots = boundaryKnots_pcr) * d_cadaveric
  
fixedValueFormula_creatinine = ~ 1 + rec_age_fwp1 + rec_gender + d_age + 
  tx_pra + ah_nr + tx_dm + 
  ns(tx_s_years,knots=fixedKnots_creatinine, Boundary.knots = boundaryKnots_creatinine) * d_cadaveric + 
  ns(tx_s_years,knots=fixedKnots_creatinine, Boundary.knots = boundaryKnots_creatinine) * tx_dgf

randomValueFormula_pcr = ~ 1 + ns(tx_s_years, knots = randomKnots_pcr, Boundary.knots = boundaryKnots_pcr)
randomValueFormula_creatinine = ~ 1 + ns(tx_s_years, knots = randomKnots_creatinine, Boundary.knots = boundaryKnots_creatinine)

fixedSlopeFormula_pcr = ~0 + dns(tx_s_years, knots=fixedKnots_pcr, Boundary.knots = boundaryKnots_pcr) + 
  I(dns(tx_s_years, knots=fixedKnots_pcr, Boundary.knots = boundaryKnots_pcr)*(as.numeric(rec_gender)-1)) + 
  I(dns(tx_s_years, knots=fixedKnots_pcr, Boundary.knots = boundaryKnots_pcr)*(as.numeric(d_cadaveric)-1))

fixedSlopeFormula_creatinine = ~ 0 + dns(tx_s_years,knots=fixedKnots_creatinine, Boundary.knots = boundaryKnots_creatinine) + 
  I(dns(tx_s_years,knots=fixedKnots_creatinine, Boundary.knots = boundaryKnots_creatinine)*(as.numeric(d_cadaveric)-1)) + 
  I(dns(tx_s_years,knots=fixedKnots_creatinine, Boundary.knots = boundaryKnots_creatinine)*(as.numeric(tx_dgf)-1))

randomSlopeFormula_pcr = ~ 0 + dns(tx_s_years, knots = randomKnots_pcr, Boundary.knots = boundaryKnots_pcr)
randomSlopeFormula_creatinine = ~ 0 + dns(tx_s_years, knots = randomKnots_creatinine, Boundary.knots = boundaryKnots_creatinine)

###############################################
# Design matrices for the longitudinal measurement model
###############################################
#For the age variable check the simulated and observed data distribution
longTimes = do.call(what = rbind, lapply(1:nSub, function(x){generateLongtiudinalTimeBySchedule()}))

timesPerSubject = rep(nrow(longTimes) / nSub, nSub)

qplot(x = c(amctx.id$rec_age, rnorm(nrow(amctx.id), mean=mean(amctx.id$rec_age), sqrt(var(amctx.id$rec_age)))), 
      geom="density", color=c(rep("Obs", nrow(amctx.id)), rep("Sim", nrow(amctx.id))))
qplot(x = c(amctx.id$d_age, rnorm(nrow(amctx.id), mean=mean(amctx.id$d_age), sqrt(var(amctx.id$d_age)))), 
      geom="density", color=c(rep("Obs", nrow(amctx.id)), rep("Sim", nrow(amctx.id))))
qplot(x = c(amctx.id$rec_bmi, rnorm(nrow(amctx.id), mean=mean(amctx.id$rec_bmi), sqrt(var(amctx.id$rec_bmi)))), 
      geom="density", color=c(rep("Obs", nrow(amctx.id)), rep("Sim", nrow(amctx.id))))

#tx_dial_days has an outlier at positiion: which.max(amctx.id$tx_dial_days)
qplot(x = c(amctx.id$tx_dial_days[-which.max(amctx.id$tx_dial_days)], 
            rnorm(nrow(amctx.id[-which.max(amctx.id$tx_dial_days),]), 
                  mean=mean(amctx.id$tx_dial_days[-which.max(amctx.id$tx_dial_days)]), 
                  sqrt(var(amctx.id$tx_dial_days[-which.max(amctx.id$tx_dial_days)])))), 
      geom="density", color=c(rep("Obs", nrow(amctx.id[-which.max(amctx.id$tx_dial_days),])), 
                              rep("Sim", nrow(amctx.id[-which.max(amctx.id$tx_dial_days),]))))


#make rec_age_fwp1 for the longitudinal data set
rec_age = rnorm(n = nSub, mean=mean(amctx.id$rec_age), sd=sqrt(var(amctx.id$rec_age)))
d_age = rnorm(n = nSub, mean=mean(amctx.id$d_age), sd=sqrt(var(amctx.id$d_age)))
rec_bmi = rnorm(n = nSub, mean=mean(amctx.id$rec_bmi), sd=sqrt(var(amctx.id$rec_bmi)))
tx_dial_days = rnorm(nSub, mean=mean(amctx.id$tx_dial_days[-which.max(amctx.id$tx_dial_days)]), 
                       sqrt(var(amctx.id$tx_dial_days[-which.max(amctx.id$tx_dial_days)])))
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
tx_pra = sample(c(tx_pra_zeros, round(rgamma(nSub-length(tx_pra_zeros), shape=0.5, scale=20))+2), replace = F)

subId <- rep(1:nSub)
simDs.id = data.frame(amctx = subId, rec_age, d_age, rec_bmi, tx_dial_days, tx_previoustx,
                      d_gender, rec_gender, d_cadaveric, tx_dgf, tx_dm, tx_pra, ah_nr)

simDs = do.call(rbind, lapply(1:nSub, function(i){
  data.frame(amctx = rep(subId[i], each=timesPerSubject[i]), 
                   rec_age_fwp1 = rep(rec_age[i], each=timesPerSubject[i]),
                   d_age = rep(d_age[i], each=timesPerSubject[i]),
                   rec_bmi = rep(rec_bmi[i], each=timesPerSubject[i]),
                   tx_dial_days = rep(tx_dial_days[i], each=timesPerSubject[i]),
                   tx_previoustx = rep(tx_previoustx[i], each=timesPerSubject[i]),
                   d_gender = rep(d_gender[i], each=timesPerSubject[i]),
                   rec_gender = rep(rec_gender[i], each=timesPerSubject[i]),
                   d_cadaveric = rep(d_cadaveric[i], each=timesPerSubject[i]),
                   tx_dgf = rep(tx_dgf[i], each=timesPerSubject[i]),
                   tx_dm = rep(tx_dm[i], each=timesPerSubject[i]),
                   ah_nr = rep(ah_nr[i], each=timesPerSubject[i]),
                   tx_pra = rep(tx_pra[i], each=timesPerSubject[i]))
  }))
simDs = cbind(simDs, longTimes)
simDs$logPcr = NA
simDs$logCreatinine = NA

X_pcr = model.matrix(fixedValueFormula_pcr, data = simDs[simDs$isPcr==T,])
X_creatinine = model.matrix(fixedValueFormula_creatinine, data = simDs[simDs$isCreatinine==T,])

Z_pcr = model.matrix(randomValueFormula_pcr, data = simDs[simDs$isPcr==T,])
Z_creatinine = model.matrix(randomValueFormula_creatinine, data = simDs[simDs$isCreatinine==T,])

betas_pcr <- getBetasPCR(fittedJointModel)
betas_creatinine <- getBetasCreatinine(fittedJointModel)
sigma.y_pcr <- getSigmaPCR(fittedJointModel)
sigma.y_creatinine <- getSigmaCreatinine(fittedJointModel)

D <- getD(fittedJointModel)
b <- mvrnorm(nSub, mu = rep(0, nrow(D)), D)
b_pcr <- b[,1:4]
b_creatinine <- b[,5:8]

visitCumsum = c(0, cumsum(timesPerSubject))
Zb_pcr = unlist(lapply(1:nSub,function(i){
  Z_pcr[(visitCumsum[i] + 1): visitCumsum[i+1],] %*% b_pcr[i, ]
}))

Zb_creatinine = unlist(lapply(1:nSub,function(i){
  Z_creatinine[(visitCumsum[i] + 1): visitCumsum[i+1],] %*% b_creatinine[i, ]
}))

xBetaZb_pcr = X_pcr %*% betas_pcr + Zb_pcr
xBetaZb_creatinine = X_creatinine %*% betas_creatinine + Zb_creatinine

simDs$logPcr[simDs$isPcr==T] = rnorm(length(xBetaZb_pcr), xBetaZb_pcr, sigma.y_pcr)
simDs$logCreatinine[simDs$isCreatinine==T] = rnorm(length(xBetaZb_creatinine), xBetaZb_creatinine, sigma.y_creatinine)

ggplot(data=simDs[simDs$isPcr==T,], aes(x=tx_s_years, y=logPcr)) + 
  geom_line(aes(group=amctx)) + stat_smooth()

ggplot(data=simDs[simDs$isCreatinine==T,], aes(x=tx_s_years, y=logCreatinine)) + 
  geom_line(aes(group=amctx)) + stat_smooth()

# design matrix for the survival model
survivalFormula = ~ 1 + rec_age + d_age + tx_previoustx + d_gender + rec_bmi + tx_pra + I(tx_dial_days/365)
W <- model.matrix(survivalFormula, data = simDs.id)[,-1] #drop the intercept
gammas = getGamma(fittedJointModel)
wGamma <- as.vector(W %*% gammas)

weibullShape = 1
weibullScale = 7

u <- runif(nSub)
simDs.id$years_tx_gl <- NA
for (i in 1:nSub) {
    Up <- 10
    tries <- 5
    Root <- try(uniroot(invSurvival, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    while(inherits(Root, "try-error") && tries > 0) {
        tries <- tries - 1
        Up <- Up + 200
        Root <- try(uniroot(invSurvival, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    }
    simDs.id$years_tx_gl[i] <- if (!inherits(Root, "try-error")) Root else NA
}

pid_to_keep = simDs.id[!is.na(simDs.id$years_tx_gl),]$amctx

simDs = simDs[simDs$amctx %in% pid_to_keep,]
simDs.id = simDs.id[simDs.id$amctx %in% pid_to_keep,]
simDs$amctx = droplevels(simDs$amctx)
simDs.id$amctx = droplevels(simDs.id$amctx)

#Divide into training and test
trainingSize = round(nrow(simDs.id)/2)
trainingDs.id = simDs.id[1:trainingSize, ]
testDs.id = simDs.id[(trainingSize + 1):nSub,]
trainingDs.id$amctx = droplevels(trainingDs.id$amctx)
testDs.id$amctx = droplevels(testDs.id$amctx)

trainingDs = simDs[simDs$amctx %in% trainingDs.id$amctx, ]
trainingDs$amctx = droplevels(trainingDs$amctx)
testDs = simDs[simDs$amctx %in% testDs.id$amctx, ]
testDs$amctx = droplevels(testDs$amctx)

# simulate censoring times from an exponential distribution for TRAINING DATA SET ONLY
# and calculate the observed event times, i.e., min(true event times, censoring times)
Ctimes <- rexp(trainingSize, 1/mean(prias.id[prias.id$gl_failure==0,]$years_tx_gl))
trainingDs.id$gl_failure = trainingDs.id$years_tx_gl > Ctimes
trainingDs.id$years_tx_gl = min(trainingDs.id$years_tx_gl, Ctimes)
trainingDs$years_tx_gl = rep(trainingDs.id$years_tx_gl, each=timesPerSubject)

# drop the longitudinal measurementsthat were taken after the observed event time for each subject.
trainingDs = trainingDs[trainingDs$tx_s_years <= trainingDs$years_tx_gl, ]

# delete all unused objects
rm(y, X, Z, id, n, na.ind, long.na.ind, ind, Ctimes, Time, event, W, 
    betas, sigma.y, gammas, alpha, eta.t, eta.y, phi, mean.Cens, t.max,
    trueTimes, u, Root, invSurvival, D, b, K, 
    times, group, i, tries, Up, boundaryKnots, fixedKnots, DF)
