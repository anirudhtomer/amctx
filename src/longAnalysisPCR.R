#####################################################
# Longitudinal analysis for PCR
#####################################################

# First check the trend overall
ds = amctx_pcr[amctx_pcr$tx_s_days>0 & amctx_pcr$tx_s_days<=200,]
ds$amctx = droplevels(ds$amctx)
ggplot(data=ds, aes(x=tx_s_years,y=log(value))) + geom_line(aes(group=amctx)) + stat_smooth() + 
  facet_grid(rec_gender~d_cadaveric, labeller = label_both)

#For d_cadaveric or d_type difference evolution
ggplot(data=amctx_pcr, aes(x=tx_s_years,y=log(value))) + 
  geom_line(aes(group=amctx)) + stat_smooth() + ticks(0, 3, 0.05) + ylab("log(pcr)") + 
  xlab("tx_s_years") + facet_grid(.~d_cadaveric) + xlim(0,3)

#For d_cadaveric or d_type difference evolution
ggplot(data=amctx_pcr, aes(x=tx_s_years,y=log(value))) + 
  geom_line(aes(group=amctx)) + stat_smooth() + ticks(0, 5, 0.5) + ylab("log(pcr)") + 
  xlab("tx_s_years") + facet_grid(.~factor(rec_gender)) + xlim(0,3)

#Firsly log transform is better. 
#Secondly I could see some striations on that plot at the bottom, which are marked now
ggplot(data=amctx_pcr, aes(x=tx_s_days,y=log(value), group=amctx, color=factor(log(value)))) + geom_point() + ticks(200) +
  guides(color=F) + 
  geom_abline(intercept = log(2), slope=0) +
  geom_abline(intercept = log(3), slope=0) +
  geom_abline(intercept = log(4), slope=0) +
  geom_abline(intercept = log(5), slope=0) 

# ah_ace + ah_arb = ah_raasi, perhaps it’s better to not include all three in the model. 
# Preferably no medication use is included, since they have not been collected in a dynamic way 
# (just baseline) and without doses. 
# tx_hla can be modelled as a continuous parameter
# (hessel: maybe not statistically correct, but this is often done also to reduce parameters and error)  
model_pcr_feedback0 = lme(data=amctx_merged[!is.na(amctx_merged$pcr),], fixed=log(pcr)~d_age + rec_bmi + d_type + d_bmi + tx_cit+ 
                  tx_hla+ rec_age_fwp1 + tx_previoustx + tx_dial_days + 
                  tx_dm + tx_pra + 
                  ns(tx_s_years,knots=c(100, 200, 350)/365),
                random = ~ns(tx_s_years,knots=c(100, 200)/365)|amctx,
                control = lmeControl(opt = "optim"), method="ML")

model_pcr_feedback1 = lme(data=amctx_merged[!is.na(amctx_merged$pcr),], 
                          fixed=log(pcr)~rec_age_fwp1 + 
                            rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi + 
                            tx_hla + tx_pra + tx_dgf + I(ah_nr>0) + tx_dm + tx_hvdis + 
                            rr_sys + rr_dias +
                            rr_map + tx_dial_days + 
                            ns(tx_s_years,knots=c(50, 200, 365)/365),
                          random = ~ns(tx_s_years,knots=c(50, 200)/365)|amctx,
                          control = lmeControl(opt = "optim"), method="ML")

model_pcr_feedback2 = lme(data=amctx_merged[!is.na(amctx_merged$pcr),], 
                          fixed=log(pcr)~rec_age_fwp1 + 
                            rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi + 
                            tx_hla + tx_pra + tx_dgf + I(ah_nr>0) + tx_dm + tx_hvdis + 
                            rr_sys + rr_dias +
                            rr_map + tx_dial_days + 
                            ns(tx_s_years,knots=c(50, 200, 365)/365) * rec_gender + 
                            ns(tx_s_years,knots=c(50, 200, 365)/365) * d_cadaveric,
                          random = ~ns(tx_s_years,knots=c(50, 200)/365)|amctx,
                          control = lmeControl(opt = "optim"), method="ML")

model_pcr_feedback3 = lme(data=amctx_pcr, 
                          fixed=log(pcr)~ 
                            rec_gender + d_age + 
                            ns(tx_s_years,knots=c(50, 200, 365)/365) * rec_gender + 
                            ns(tx_s_years,knots=c(50, 200, 365)/365) * d_cadaveric,
                          random = ~ns(tx_s_years,knots=c(50, 200)/365)|amctx,
                          control = lmeControl(opt = "optim"), method="ML")

model_pcr_feedback4 = lme(data=amctx_pcr, 
                          fixed=log(pcr)~ 
                            rec_age_fwp1 + 
                            rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi + 
                            tx_hla + tx_pra + tx_dgf + ah_nr + tx_dm + tx_hvdis + 
                            rr_sys + rr_dias + tx_cit + 
                            rr_map + tx_dial_days + d_cadaveric +
                            ns(tx_s_years,knots=c(30, 80, 365)/365),
                          random = ~ns(tx_s_years,knots=c(30, 80, 365)/365)|amctx,
                          control = lmeControl(opt = "optim"), method="ML")

model_pcr_feedback5 = lme(data=amctx_pcr, 
                          fixed=log(pcr)~ 
                            rec_age_fwp1 + 
                            rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi + 
                            tx_hla + tx_pra + tx_dgf + ah_nr + tx_dm + tx_hvdis + 
                            rr_sys + rr_dias + tx_cit + 
                            rr_map + tx_dial_days + d_cadaveric +
                            ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 5.5)),
                          random = ~ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 5.5))|amctx,
                          control = lmeControl(opt = "optim"), method="ML")

model_pcr_feedback6 = lme(data=amctx_pcr, 
                          fixed=log(pcr)~ 
                            rec_age_fwp1 + 
                            rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi + 
                            tx_hla + tx_pra + tx_dgf + ah_nr + tx_dm + tx_hvdis + 
                            rr_sys + rr_dias + tx_cit + 
                            rr_map + tx_dial_days + d_cadaveric +
                            ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6)),
                          random = ~ns(tx_s_years,knots=c(30, 80, 365)/365, Boundary.knots = c(0.03917808, 6))|amctx,
                          control = lmeControl(opt = "optim"), method="ML")



model1 = fitUnivaritePCRModel(fixedSplineKnots = c(30, 70, 1000)/365, 
                              randomSplineKnots = c(30, 70)/365,
                              boundaryKnots = c(0, 5.5))

plotPCRFittedCurve(list(model_pcr_feedback3, model1), individually = F, transform = F)

anova(model_pcr_feedback1, model_pcr_feedback2, model_pcr_feedback3)
anova(model_pcr_feedback1, model_pcr_feedback2)

