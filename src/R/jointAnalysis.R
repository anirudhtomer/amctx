library(JMbayes)

######################################
#Without the multivariate functionality
######################################
jmbayes_creatinine_tdval = jointModelBayes(lmeObject = main_effect_complex_model5, survObject = coxModel_clinical, timeVar = "tx_s_years", control = list(n.iter=1000))
jmbayes_creatinine_tdboth = update(jmbayes_creatinine_tdval, param = "td-both", 
                                   extraForm = list(fixed = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                                                    random = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)), 
                                                    indFixed = c(18:21), indRandom = 2:5))

jmbayes_creatinine_tdboth_replaced = replaceMCMCContents(fromObj = mvJoint_creatinine_tdboth_complex, toObj = jmbayes_creatinine_tdboth)

save.image("Rdata/feedbackmeeting.Rdata")
####################################################
# Multivariate functionality
####################################################
amctx_creatinine$creatinine = amctx_creatinine$value
mvglmer_creatinine=mvglmer(list(log(creatinine) ~ rec_age_fwp1 + 
                                   rec_gender + d_age + 
                                   tx_pra + ah_nr + tx_dm + 
                                   ns(tx_s_years,knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6)) * d_cadaveric + 
                                   ns(tx_s_years,knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6)) * tx_dgf + 
                                   (ns(tx_s_years, knots=c(30, 70)/365, Boundary.knots = c(0, 6))|amctx)),
                            data = amctx_creatinine, families = list(gaussian))
save.image("Rdata/feedbackmeeting.Rdata")

mvJoint_creatinine_tdval=mvJointModelBayes(mvglmer_creatinine, coxModel, timeVar = "tx_s_years",
                                            priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

save.image("Rdata/feedbackmeeting.Rdata")

forms_creatinine <- list("log(creatinine)" = "value",
                         "log(creatinine)" = list(fixed = ~ 0 + dns(tx_s_years, knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6)) + 
                                                    I(dns(tx_s_years, knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6))*(as.numeric(d_cadaveric)-1)) + 
                                                    I(dns(tx_s_years, knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6))*(as.numeric(tx_dgf)-1)),
                                                  random = ~ 0 + dns(tx_s_years, knots = c(30, 70)/365, Boundary.knots = c(0, 6)), 
                                                  indFixed = c(8:11, 14:17, 18:21), indRandom = 2:4, 
                                                  name = "slope"))

mvJoint_creatinine_tdboth <- update(mvJoint_creatinine_tdval, Formulas = forms_creatinine,
                                     priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

save.image("Rdata/feedbackmeeting.Rdata")

###############################################

mvglmer_creatinine_complex=mvglmer(list(log(creatinine) ~ rec_age_fwp1 + 
                                          rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi + 
                                          tx_hla + tx_pra + tx_dgf + ah_nr + tx_dm + tx_hvdis + 
                                          tx_cit + tx_dial_days + d_cadaveric +
                                          ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)) + 
                                          (ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6))|amctx)),
                                   data = amctx_creatinine, families = list(gaussian))

save(mvglmer_creatinine_complex, file="temp_ani.Rdata")
mvJoint_creatinine_tdval=mvJointModelBayes(mvglmer_creatinine_complex, coxModel, timeVar = "tx_s_years",
                                           priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

forms_creatinine <- list("log(creatinine)" = "value",
                         "log(creatinine)" = list(fixed = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                                                  random = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)), 
                                                  indFixed = c(18:21), indRandom = 2:5, 
                                                  name = "slope"))

mvJoint_creatinine_tdboth_complex <- mvJointModelBayes(mvglmer_creatinine_complex, coxModel_clinical,  timeVar = "tx_s_years",
                                    Formulas = forms_creatinine,
                                    priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

##############################################
mvglmer_pcr=mvglmer(list(log(pcr) ~ rec_gender + d_age + 
                           ns(tx_s_years,knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 5.5)) * rec_gender + 
                           ns(tx_s_years,knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 5.5)) * d_cadaveric + 
                            (ns(tx_s_years, knots=c(30, 70)/365, Boundary.knots = c(0, 5.5))|amctx)),
                            data = amctx_merged, families = list(gaussian))

save.image("Rdata/feedbackmeeting.Rdata")

mvJoint_pcr_tdval=mvJointModelBayes(mvglmer_pcr, coxModel, timeVar = "tx_s_years",
                                            priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))
save.image("Rdata/feedbackmeeting.Rdata")

forms_pcr <- list("log(pcr)" = "value",
                  "log(pcr)" = list(fixed = ~ 0 + dns(tx_s_years,knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 5.5)) + 
                                      I(dns(tx_s_years,knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 5.5))*(as.numeric(rec_gender)-1)) + 
                                      I(dns(tx_s_years,knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 5.5))*(as.numeric(d_cadaveric)-1)), 
                                                  random = ~ 0 + dns(tx_s_years, knots = c(30, 70)/365, Boundary.knots = c(0, 5.5)), 
                                                  indFixed = c(4:7, 9:12, 13:16), indRandom = 2:4, 
                                                  name = "slope"))

mvJoint_pcr_tdboth <- update(mvJoint_pcr_tdval, Formulas = forms_pcr,
                                     priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))
save.image("Rdata/feedbackmeeting.Rdata")

# Both PCR and creatinine together
mvglmer_pcr_creatinine=mvglmer(list(log(pcr) ~ 
                                      ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)) + 
                                      (ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6))|amctx),
                                  
                                    log(creatinine) ~ 
                                      ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)) + 
                                    (ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6))|amctx)),
                             data = amctx_merged, families = list(gaussian, gaussian))

save.image("Rdata/feedbackmeeting.Rdata")

mvJoint_pcr_creatinine_tdval=mvJointModelBayes(mvglmer_pcr_creatinine, coxModel, timeVar = "tx_s_years",
                                                priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

save.image("Rdata/feedbackmeeting.Rdata")

forms_pcr_creatinine <- list("log(creatinine)" = "value",
                             "log(creatinine)" = list(fixed = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                                                      random = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                                                      indFixed = c(2:5), indRandom = 2:5, 
                                                      name = "slope"),
                             "log(pcr)" = "value",
                             "log(pcr)" = list(fixed = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                                               random = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                                               indFixed = c(2:5), indRandom = 2:5, 
                                               name = "slope"))
save.image("Rdata/feedbackmeeting.Rdata")


forms_pcr_creatinine <- list("log(creatinine)" = "value",
                             "log(creatinine)" = list(fixed = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)) +  
                                                        I(dns(tx_s_years, knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6))*(as.numeric(d_cadaveric)-1)) + 
                                                        I(dns(tx_s_years, knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 6))*(as.numeric(tx_dgf)-1)),
                                                      random = ~ 0 + dns(tx_s_years, knots = c(30, 70)/365, Boundary.knots = c(0, 6)), 
                                                      indFixed = c(8:11, 14:17, 18:21), indRandom = 2:4, 
                                                      name = "slope"),
                             "log(pcr)" = "value",
                             "log(pcr)" = list(fixed = ~ 0 + dns(tx_s_years,knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 5.5)) + 
                                                 I(dns(tx_s_years,knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 5.5))*(as.numeric(rec_gender)-1)) + 
                                                 I(dns(tx_s_years,knots=c(30, 70, 1000)/365, Boundary.knots = c(0, 5.5))*(as.numeric(d_cadaveric)-1)), 
                                               random = ~ 0 + dns(tx_s_years, knots = c(30, 70)/365, Boundary.knots = c(0, 5.5)), 
                                               indFixed = c(4:7, 9:12, 13:16), indRandom = 2:4, 
                                               name = "slope"))


mvJoint_pcr_creatinine_tdboth <- mvJointModelBayes(mvglmer_pcr_creatinine, coxNull, timeVar = "tx_s_years", Formulas = forms_pcr_creatinine,
                                         priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

save.image("Rdata/feedbackmeeting.Rdata")


qplot(y=residuals(longModel, type = "response"), x=fitted(longModel), geom=c("point", "smooth"))
qplot(y=residuals(jointFit, process="Longitudinal", type = "marginal"), x=fitted(jointFit, process="Longitudinal", type = "marginal"), geom=c("point", "smooth"))

plot(density(jointFit_ns1random$mcmc$betas[,2]))


ggsList = as.ggsList(jointFit_creatinine_pcr_tdboth)

#The following shows that not everything has converged,D matrix and random effects b
ggs_density(ggsList[[5]][[1]])
ggs_density(ggsList[[7]][[238]])
