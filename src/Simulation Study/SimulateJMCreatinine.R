library(MASS)
library(splines)

source("src/R/common.R")
source("src/R/Simulation Study/simCommon.R")

set.seed(1202)
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

weibullShapes = rep(1.35, nSub)
weibullScales = rep(11000, nSub)

simulatedDs = generateSimLongitudinalData(nSub)
plotTrueSurvival(sample(1:10, size = 1))

simulatedDs = generateSimJointData(nSub, simulatedDs$simDs, simulatedDs$simDs.id,
                                   simulatedDs$b, simulatedDs$wGamma)

graphdf = data.frame(failtime = c(amctx.id$years_tx_gl, simulatedDs$trainingDs.id$years_tx_gl),
                     indicator = c(amctx.id$gl_failure, simulatedDs$trainingDs.id$gl_failure),
                     type = rep(c("amctx", "sim"), c(nrow(amctx.id), nrow(simulatedDs$trainingDs.id))))

ggplot(data=graphdf) + geom_density(aes(x=failtime, fill = type), alpha=0.5) + facet_grid(.~indicator)

weibullScales = simulatedDs$weibullScales
weibullShapes = simulatedDs$weibullShapes

trainingDs = simulatedDs$trainingDs
trainingDs.id = simulatedDs$trainingDs.id
testDs = simulatedDs$testDs
testDs.id = simulatedDs$testDs.id
testDs.id$gl_failure = 1

simDs = simulatedDs$simDs
simDs.id = simulatedDs$simDs.id
b_creatinine = simulatedDs$b
wGamma = simulatedDs$wGamma

rm(simulatedDs)
trainingDs$creatinine = exp(trainingDs$logCreatinine)
testDs$creatinine = exp(testDs$logCreatinine)

#########################################
# Fit joint model for the simulated data set
########################################
cox_Model_training = coxph(Surv(years_tx_gl, gl_failure) ~ tx_hla +
                             tx_previoustx + tx_cit + tx_dial_days,
                           data = trainingDs.id, x=T, model=T)

mvglmer_creatinine_training=mvglmer(list(log(creatinine) ~ rec_age_fwp1 +
                                           rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi +
                                           tx_hla + tx_pra + tx_dgf + ah_nr + tx_dm + tx_hvdis +
                                           tx_cit + tx_dial_days + d_cadaveric +
                                           ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)) +
                                           (ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6))|amctx)),
                                    data = trainingDs, families = list(gaussian))

forms_creatinine_training <- list("log(creatinine)" = "value",
                                  "log(creatinine)" = list(fixed = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                                                           random = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                                                           indFixed = c(18:21), indRandom = 2:5,
                                                           name = "slope"))

#mvJoint_creatinine_tdval_training=mvJointModelBayes(mvglmer_creatinine_training, cox_Model_training, timeVar = "tx_s_years",
#                                                    priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

mvJoint_creatinine_tdboth_training = mvJointModelBayes(mvglmer_creatinine_training, cox_Model_training, timeVar = "tx_s_years",
                                                       Formulas = forms_creatinine_training,
                                                       priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

###########################################
#Also fit joint model using the jointModelAPI
###########################################
lme_creatinine_training = lme(data= trainingDs,
                              fixed=log(creatinine) ~ rec_age_fwp1 + 
                                rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi + 
                                tx_hla + tx_pra + tx_dgf + ah_nr + tx_dm + tx_hvdis + 
                                tx_cit + tx_dial_days + d_cadaveric +
                                ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                              random = ~ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6))|amctx,
                              control = lmeControl(opt = "optim"), method="ML")

#jmbayes_creatinine_tdval_training = jointModelBayes(lmeObject = lme_creatinine_training, survObject = cox_Model_training, timeVar = "tx_s_years", control = list(n.iter=1000))
jmbayes_creatinine_tdboth_training = jointModelBayes(lme_creatinine_training, param = "td-both", 
                                                     survObject = cox_Model_training, timeVar = "tx_s_years", 
                                                     control = list(n.iter=1000),
                                            extraForm = list(fixed = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                                                             random = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)), 
                                                             indFixed = c(18:21), indRandom = 2:5))

simJointModel_replaced = replaceMCMCContents(mvJoint_creatinine_tdboth_training, jmbayes_creatinine_tdboth_training)
save.image("Rdata/simCreatinine.Rdata")


