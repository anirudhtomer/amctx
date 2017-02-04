library(JMbayes)

######################################
#Without the multivariate functionality
######################################
jointfit_creatinine_nomv = jointModelBayes(model_creatinine, coxModel, timeVar = "tx_s_years", n.iter = 1000)
jointfit_creatinine_tdboth_nomv = update(jointfit_creatinine_nomv, param = "td-both", 
                                         extraForm = list(fixed = ~ 0 + dns(tx_s_years, knots=c(100, 300, 1000)/365, Boundary.knots = c(0, 11)), 
                                                          random = ~ 0 + dns(tx_s_years, knots = c(100, 300)/365, Boundary.knots = c(0, 11)),
                                                          indFixed = 20:23, indRandom = 2:4))

####################################################
# Multivariate functionality
####################################################

#remove subject 346 because
#More specifically, according to the data the number of subjects per response is log(pcr) = 238, log(creatinine) = 239.

model_mv_creatinine=mvglmer(list(log(creatinine) ~ rec_age_fwp1 + rec_gender +
                                   d_age +  tx_dgf + is_cni + d_bmi + tx_hla + tx_previoustx + 
                                   tx_pra + tx_cit + tx_dial_days + tx_dm + rec_bmi + 
                                   ns(tx_s_years, knots=(c(100, 300, 1000)/365), Boundary.knots = c(0, 11)) + 
                                   (ns(tx_s_years, knots=(c(100, 300)/365), Boundary.knots = c(0, 11))|amctx)),
                            data = amctx_merged, families = list(gaussian))

jointFit_creatinine_tdval=mvJointModelBayes(model_mv_creatinine, coxModel, timeVar = "tx_s_years",
                                            priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

forms_creatinine <- list("log(creatinine)" = "value",
                         "log(creatinine)" = list(fixed = ~ 0 + dns(tx_s_years, knots=c(100, 300, 1000)/365, Boundary.knots = c(0, 11)), 
                                                  random = ~ 0 + dns(tx_s_years, knots = c(100, 300)/365, Boundary.knots = c(0, 11)), indFixed = 20:23, indRandom = 2:4, 
                                                  name = "slope"))

jointFit_creatinine_tdboth <- update(jointFit_creatinine_tdval, Formulas = forms_creatinine,
                                     priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))
# So both the current value of creatinine and the slope are useful, now we add pcr and its slope
model_pcr_creatinine=mvglmer(list(log(pcr) ~ d_age + ah_ace + ah_arb + 
                                    ah_raasi  + rec_bmi + d_type + d_bmi + tx_cit+ 
                                    tx_hla+ rec_age_fwp1 + tx_previoustx + tx_dial_days + 
                                    tx_dm + tx_pra + ah_nr + 
                                    ns(tx_s_years,knots=c(100, 200, 350)/365, Boundary.knots = c(0, 11)) +
                                    (ns(tx_s_years, knots=c(100, 200)/365, Boundary.knots = c(0, 11))|amctx),
                                  
                                  log(creatinine) ~ rec_age_fwp1 + rec_gender +
                                    d_age +  tx_dgf + is_cni + d_bmi + tx_hla + tx_previoustx + 
                                    tx_pra + tx_cit + tx_dial_days + tx_dm + rec_bmi + 
                                    ns(tx_s_years, knots=(c(100, 300, 1000)/365), Boundary.knots = c(0, 11)) + 
                                    (ns(tx_s_years, knots=(c(100, 300)/365), Boundary.knots = c(0, 11))|amctx)),
                             data = amctx_merged, families = list(gaussian, gaussian))

jointFit_creatinine_pcr_tdval=mvJointModelBayes(model_pcr_creatinine, coxModel, timeVar = "tx_s_years",
                                                priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

forms_creatinine_pcr <- list("log(creatinine)" = "value",
                             "log(creatinine)" = list(fixed = ~ 0 + dns(tx_s_years, knots=c(100, 300, 1000)/365, Boundary.knots = c(0, 11)), 
                                                      random = ~ 0 + dns(tx_s_years, knots = c(100, 300)/365, Boundary.knots = c(0, 11)), indFixed = 20:23, indRandom = 2:4, 
                                                      name = "slope"),
                             "log(pcr)" = "value",
                             "log(pcr)" = list(fixed = ~ 0 + dns(tx_s_years, knots=c(100, 200, 350)/365, Boundary.knots = c(0, 11)), 
                                               random = ~ 0 + dns(tx_s_years, knots = c(100, 200)/365, Boundary.knots = c(0, 11)), indFixed = 20:23, indRandom = 2:4, 
                                               name = "slope"))
jointFit_creatinine_pcr_tdboth <- update(jointFit_creatinine_pcr_tdval, Formulas = forms_creatinine_pcr,
                                         priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))


jointFit_creatinine_pcr_tdval_noridge=mvJointModelBayes(model_pcr_creatinine, coxModel, timeVar = "tx_s_years")

forms_creatinine_pcr <- list("log(creatinine)" = "value",
                             "log(creatinine)" = list(fixed = ~ 0 + dns(tx_s_years, knots=c(100, 300, 1000)/365, Boundary.knots = c(0, 11)), 
                                                      random = ~ 0 + dns(tx_s_years, knots = c(100, 300)/365, Boundary.knots = c(0, 11)), indFixed = 20:23, indRandom = 2:4, 
                                                      name = "slope"),
                             "log(pcr)" = "value",
                             "log(pcr)" = list(fixed = ~ 0 + dns(tx_s_years, knots=c(100, 200, 350)/365, Boundary.knots = c(0, 11)), 
                                               random = ~ 0 + dns(tx_s_years, knots = c(100, 200)/365, Boundary.knots = c(0, 11)), indFixed = 20:23, indRandom = 2:4, 
                                               name = "slope"))
jointFit_creatinine_pcr_tdboth_noridge <- update(jointFit_creatinine_pcr_tdval_noridge, Formulas = forms_creatinine_pcr)
                                         



forms_creatinine_pcr_2 <- list("log(creatinine)" = "value",
                             "log(creatinine)" = list(fixed = ~ 0 + dns(tx_s_years, knots=c(100, 300, 1000)/365, Boundary.knots = c(0, 11)), 
                                                      random = ~ 0 + dns(tx_s_years, knots = c(100, 300)/365, Boundary.knots = c(0, 11)), indFixed = 20:23, indRandom = 2:4, 
                                                      name = "slope"),
                             "log(pcr)" = "value")
jointFit_creatinine_pcr_tdboth_2 <- update(jointFit_creatinine_pcr_tdval, Formulas = forms_creatinine_pcr_2,
                                         priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))


qplot(y=residuals(longModel, type = "response"), x=fitted(longModel), geom=c("point", "smooth"))
qplot(y=residuals(jointFit, process="Longitudinal", type = "marginal"), x=fitted(jointFit, process="Longitudinal", type = "marginal"), geom=c("point", "smooth"))

plot(density(jointFit_ns1random$mcmc$betas[,2]))


ggsList = as.ggsList(jointFit_creatinine_pcr_tdboth)

#The following shows that not everything has converged,D matrix and random effects b
ggs_density(ggsList[[5]][[1]])
ggs_density(ggsList[[7]][[238]])
