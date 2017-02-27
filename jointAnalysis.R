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
mvglmer_creatinine=mvglmer(list(log(creatinine) ~ rec_age_fwp1 + 
                                   rec_gender + d_age + 
                                   tx_pra + ah_nr + tx_dm + 
                                   ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0, 10.5)) * d_cadaveric + 
                                   ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0, 10.5)) * tx_dgf + 
                                   (ns(tx_s_years, knots=c(50, 100)/365, Boundary.knots = c(0, 10.5))|amctx)),
                            data = amctx_merged, families = list(gaussian))
save.image("feedbackmeeting.Rdata")

mvJoint_creatinine_tdval=mvJointModelBayes(mvglmer_creatinine, coxModel, timeVar = "tx_s_years",
                                            priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

save.image("feedbackmeeting.Rdata")

forms_creatinine <- list("log(creatinine)" = "value",
                         "log(creatinine)" = list(fixed = ~ 0 + dns(tx_s_years, knots=c(50, 100, 900)/365, Boundary.knots = c(0, 10.5)) + 
                                                    I(dns(tx_s_years, knots=c(50, 100, 900)/365, Boundary.knots = c(0, 10.5))*(as.numeric(d_cadaveric)-1)) + 
                                                    I(dns(tx_s_years, knots=c(50, 100, 900)/365, Boundary.knots = c(0, 10.5))*(as.numeric(tx_dgf)-1)),
                                                  random = ~ 0 + dns(tx_s_years, knots = c(50, 100)/365, Boundary.knots = c(0, 10.5)), 
                                                  indFixed = c(8:11, 14:17, 18:21), indRandom = 2:4, 
                                                  name = "slope"))

mvJoint_creatinine_tdslope <- update(mvJoint_creatinine_tdval, Formulas = forms_creatinine,
                                     priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

save.image("feedbackmeeting.Rdata")
##############################################
mvglmer_pcr=mvglmer(list(log(pcr) ~ rec_gender + d_age + 
                           ns(tx_s_years, knots=c(50, 200, 365)/365, Boundary.knots = c(0, 10.5)) * rec_gender + 
                           ns(tx_s_years, knots=c(50, 200, 365)/365, Boundary.knots = c(0, 10.5)) * d_cadaveric + 
                            (ns(tx_s_years, knots=c(50, 200)/365, Boundary.knots = c(0, 10.5))|amctx)),
                            data = amctx_merged, families = list(gaussian))

save.image("feedbackmeeting.Rdata")

mvJoint_pcr_tdval=mvJointModelBayes(mvglmer_pcr, coxModel, timeVar = "tx_s_years",
                                            priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))
save.image("feedbackmeeting.Rdata")

forms_pcr <- list("log(pcr)" = "value",
                  "log(pcr)" = list(fixed = ~ 0 + dns(tx_s_years, knots=c(50, 200, 365)/365, Boundary.knots = c(0, 10.5)) + 
                                      I(dns(tx_s_years, knots=c(50, 200, 365)/365, Boundary.knots = c(0, 10.5))*(as.numeric(rec_gender)-1)) + 
                                      I(dns(tx_s_years, knots=c(50, 200, 365)/365, Boundary.knots = c(0, 10.5))*(as.numeric(d_cadaveric)-1)), 
                                                  random = ~ 0 + dns(tx_s_years, knots = c(50, 200)/365, Boundary.knots = c(0, 10.5)), 
                                                  indFixed = c(4:7, 9:12, 13:16), indRandom = 2:4, 
                                                  name = "slope"))

mvJoint_pcr_tdslope <- update(mvJoint_pcr_tdval, Formulas = forms_pcr,
                                     priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))
save.image("feedbackmeeting.Rdata")

# Both PCR and creatinine together
mvglmer_pcr_creatinine=mvglmer(list(log(pcr) ~ rec_gender + d_age + 
                                    ns(tx_s_years, knots=c(50, 200, 365)/365, Boundary.knots = c(0, 10.5)) * rec_gender + 
                                    ns(tx_s_years, knots=c(50, 200, 365)/365, Boundary.knots = c(0, 10.5)) * d_cadaveric + 
                                    (ns(tx_s_years, knots=c(50, 200)/365, Boundary.knots = c(0, 10.5))|amctx),
                                  
                                  log(creatinine) ~ rec_age_fwp1 + 
                                    rec_gender + d_age + 
                                    tx_pra + ah_nr + tx_dm + 
                                    ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0, 10.5)) * d_cadaveric + 
                                    ns(tx_s_years,knots=c(50, 100, 900)/365, Boundary.knots = c(0, 10.5)) * tx_dgf + 
                                    (ns(tx_s_years, knots=c(50, 100)/365, Boundary.knots = c(0, 10.5))|amctx)),
                             data = amctx_merged, families = list(gaussian, gaussian))

save.image("feedbackmeeting.Rdata")

mvJoint_pcr_creatinine_tdval=mvJointModelBayes(mvglmer_pcr_creatinine, coxModel, timeVar = "tx_s_years",
                                                priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

save.image("feedbackmeeting.Rdata")

forms_pcr_creatinine <- list("log(creatinine)" = "value",
                             "log(creatinine)" = list(fixed = ~ 0 + dns(tx_s_years, knots=c(50, 100, 900)/365, Boundary.knots = c(0, 10.5)) + 
                                                        I(dns(tx_s_years, knots=c(50, 100, 900)/365, Boundary.knots = c(0, 10.5))*(as.numeric(d_cadaveric)-1)) + 
                                                        I(dns(tx_s_years, knots=c(50, 100, 900)/365, Boundary.knots = c(0, 10.5))*(as.numeric(tx_dgf)-1)),
                                                      random = ~ 0 + dns(tx_s_years, knots = c(50, 100)/365, Boundary.knots = c(0, 10.5)), 
                                                      indFixed = c(8:11, 14:17, 18:21), indRandom = 2:4, 
                                                      name = "slope"),
                             "log(pcr)" = "value",
                             "log(pcr)" = list(fixed = ~ 0 + dns(tx_s_years, knots=c(50, 200, 365)/365, Boundary.knots = c(0, 10.5)) + 
                                                 I(dns(tx_s_years, knots=c(50, 200, 365)/365, Boundary.knots = c(0, 10.5))*(as.numeric(rec_gender)-1)) + 
                                                 I(dns(tx_s_years, knots=c(50, 200, 365)/365, Boundary.knots = c(0, 10.5))*(as.numeric(d_cadaveric)-1)), 
                                               random = ~ 0 + dns(tx_s_years, knots = c(50, 200)/365, Boundary.knots = c(0, 10.5)), 
                                               indFixed = c(4:7, 9:12, 13:16), indRandom = 2:4, 
                                               name = "slope"))
save.image("feedbackmeeting.Rdata")


mvJoint_pcr_creatinine_tdboth <- update(mvJoint_pcr_creatinine_tdval, Formulas = forms_pcr_creatinine,
                                         priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

mvJoint_pcr_creatinine_tdboth <- update(mvJoint_pcr_creatinine_tdval, 
                                        coxphObject = coxModel,
                                        Formulas = forms_pcr_creatinine,
                                        priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))


save.image("feedbackmeeting.Rdata")


qplot(y=residuals(longModel, type = "response"), x=fitted(longModel), geom=c("point", "smooth"))
qplot(y=residuals(jointFit, process="Longitudinal", type = "marginal"), x=fitted(jointFit, process="Longitudinal", type = "marginal"), geom=c("point", "smooth"))

plot(density(jointFit_ns1random$mcmc$betas[,2]))


ggsList = as.ggsList(jointFit_creatinine_pcr_tdboth)

#The following shows that not everything has converged,D matrix and random effects b
ggs_density(ggsList[[5]][[1]])
ggs_density(ggsList[[7]][[238]])
