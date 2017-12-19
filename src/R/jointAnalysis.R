library(JMbayes)

mvglmer_pcr_complex=mvglmer(list(log(pcr) ~ rec_age_fwp1 + d_age + d_bmi + rec_bmi + tx_hla + 
                                              tx_pra + ah_nr + tx_cit + tx_dial_days + 
                                              rec_gender + tx_previoustx + d_gender +
                                              tx_dgf + tx_dm + tx_hvdis + 
                                              d_cadaveric +
                                              ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)) + 
                                              (ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6))|amctx)),
                                       data = amctx_merged_scaled, families = list(gaussian))

save(mvglmer_pcr_complex, file="Rdata/joint models/mvglmer_pcr_complex.Rdata")

mvJoint_pcr_tdval=mvJointModelBayes(mvglmer_pcr_complex, coxModel_clinical, 
                                               timeVar = "tx_s_years",
                                               priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

save(mvJoint_pcr_tdval, file="Rdata/joint models/mvJoint_pcr_tdval.Rdata")

forms_pcr_complex <- list("log(pcr)" = "value",
                                     "log(pcr)" = list(fixed = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                                                       random = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                                                       indFixed = c(18:21), indRandom = 2:5, 
                                                       name = "slope"))

mvJoint_pcr_tdboth_complex <- mvJointModelBayes(mvglmer_pcr_complex, coxModel_clinical, 
                                                           timeVar = "tx_s_years", Formulas = forms_pcr_complex,
                                                           priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

save(mvJoint_pcr_tdboth_complex, file="Rdata/joint models/mvJoint_pcr_tdboth_complex.Rdata")

############ only creatinine
mvglmer_creatinine_complex=mvglmer(list(log(creatinine) ~ rec_age_fwp1 + d_age + d_bmi + rec_bmi + tx_hla + 
                                   tx_pra + ah_nr + tx_cit + tx_dial_days + 
                                   rec_gender + tx_previoustx + d_gender +
                                   tx_dgf + tx_dm + tx_hvdis + 
                                   d_cadaveric +
                                   ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)) + 
                                   (ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6))|amctx)),
                            data = amctx_merged_scaled, families = list(gaussian))

save(mvglmer_creatinine_complex, file="Rdata/joint models/mvglmer_creatinine_complex.Rdata")

mvJoint_creatinine_tdval=mvJointModelBayes(mvglmer_creatinine_complex, coxModel_clinical, 
                                    timeVar = "tx_s_years",
                                    priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

save(mvJoint_creatinine_tdval, file="Rdata/joint models/mvJoint_creatinine_tdval.Rdata")

forms_creatinine_complex <- list("log(creatinine)" = "value",
                          "log(creatinine)" = list(fixed = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                                            random = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                                            indFixed = c(18:21), indRandom = 2:5, 
                                            name = "slope"))

mvJoint_creatinine_tdboth_complex <- mvJointModelBayes(mvglmer_creatinine_complex, coxModel_clinical, 
                                                timeVar = "tx_s_years", Formulas = forms_creatinine_complex,
                                                priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

save(mvJoint_creatinine_tdboth_complex, file="Rdata/joint models/mvJoint_creatinine_tdboth_complex.Rdata")


# Both PCR and creatinine together
mvglmer_pcr_creatinine=mvglmer(list(log(pcr) ~ 
                                      ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)) + 
                                      (ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6))|amctx),
                                  
                                    log(creatinine) ~ 
                                      ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)) + 
                                    (ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6))|amctx)),
                             data = amctx_merged, families = list(gaussian, gaussian))

mvglmer_pcr_creatinine_complex=mvglmer(list(log(pcr) ~ rec_age_fwp1 + d_age + d_bmi + rec_bmi + tx_hla + 
                                              tx_pra + ah_nr + tx_cit + tx_dial_days + 
                                              rec_gender + tx_previoustx + d_gender +
                                               tx_dgf + tx_dm + tx_hvdis + 
                                               d_cadaveric +
                                      ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)) + 
                                      (ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6))|amctx),
                                    
                                    log(creatinine) ~ rec_age_fwp1 + d_age + d_bmi + rec_bmi + tx_hla + 
                                      tx_pra + ah_nr + tx_cit + tx_dial_days + 
                                      rec_gender + tx_previoustx + d_gender +
                                      tx_dgf + tx_dm + tx_hvdis + 
                                      d_cadaveric +
                                      ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)) + 
                                      (ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6))|amctx)),
                               data = amctx_merged_scaled, families = list(gaussian, gaussian))


save(mvglmer_pcr_creatinine_complex, file="Rdata/mvglmer_pcr_creatinine_complex.Rdata")

mvJoint_pcr_creatinine_tdval=mvJointModelBayes(mvglmer_pcr_creatinine_complex, coxModel_clinical, 
                                               timeVar = "tx_s_years",
                                                priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

save(mvJoint_pcr_creatinine_tdval, file="Rdata/mvJoint_pcr_creatinine_tdval.Rdata")

forms_pcr_creatinine_complex <- list("log(creatinine)" = "value",
                             "log(creatinine)" = list(fixed = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                                                      random = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                                                      indFixed = c(18:21), indRandom = 2:5, 
                                                      name = "slope"),
                             "log(pcr)" = "value",
                             "log(pcr)" = list(fixed = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                                               random = ~ 0 + dns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                                               indFixed = c(18:21), indRandom = 2:5, 
                                               name = "slope"))

mvJoint_pcr_creatinine_tdboth_complex <- mvJointModelBayes(mvglmer_pcr_creatinine_complex, coxModel_clinical, 
                                                   timeVar = "tx_s_years", Formulas = forms_pcr_creatinine_complex,
                                         priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))

save(mvJoint_pcr_creatinine_tdboth_complex, file="Rdata/mvJoint_pcr_creatinine_tdboth_complex.Rdata")