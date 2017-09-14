library(MASS)
library(splines)

source("src/R/common.R")
source("src/R/Simulation Study/simCommon.R")

set.seed(1002)
nSub <- 1000

# design matrix for the survival model

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
# weibullShape = c(5.5,1.8)
# weibullScale = c(70, 5000)
#weibullShape = c(5.5,1.4)
#weibullScale = c(70, 50000)

weibullShapes = rep(5.5, nSub)
weibullScales = rep(2000, nSub)

simulatedDs = generateSimLongitudinalData(nSub)
simulatedDs = generateSimJointData(nSub, simulatedDs$simDs, simulatedDs$simDs.id,
                                   simulatedDs$b, simulatedDs$wGamma)
weibullScales = simulatedDs$weibullScales
weibullShapes = simulatedDs$weibullShapes

trainingDs = simulatedDs$trainingDs
trainingDs.id = simulatedDs$trainingDs.id
testDs = simulatedDs$testDs
testDs.id = simulatedDs$testDs.id
testDs.id$gl_gailure = 1

simDs = simulatedDs$simDs
simDs.id = simulatedDs$simDs.id
b_creatinine = simulatedDs$b
wGamma = simulatedDs$wGamma

rm(simulatedDs)

#########################################
# Fit joint model for the simulated data set
########################################
cox_Model_training = coxph(Surv(years_tx_gl, gl_failure) ~ rec_age_fwp1 + 
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

mvJoint_creatinine_tdboth_training = mvJointModelBayes(mvglmer_creatinine_training, cox_Model_training, timeVar = "tx_s_years",
                                                       Formulas = forms_creatinine_training, 
                                                       priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

mvJoint_creatinine_tdboth_training <- update(mvJoint_creatinine_tdval_training, Formulas = forms_creatinine_training,
                                             priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

###########################################
#Also fit joint model using the jointModelAPI
###########################################
lme_creatinine_training = lme(data= trainingDs,
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

simJointModel_replaced = replaceMCMCContents(mvJoint_creatinine_tdboth_training, jmbayes_creatinine_tdboth_training)
save.image("Rdata/simCreatinine.Rdata")

