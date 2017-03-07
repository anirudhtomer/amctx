source("src/R/common.R")

#################################
# longitudinal analysis for creatinine
#################################

# First check the trend overall
ds = amctx_creatinine[amctx_creatinine$tx_s_days>0 & amctx_creatinine$tx_s_days<=200,]
ds$amctx = droplevels(ds$amctx)
ggplot(data=ds, aes(x=tx_s_years,y=log(value))) + geom_line(aes(group=amctx)) + stat_smooth() + 
  facet_grid(tx_dgf~d_cadaveric, labeller = label_both)

ggplot(data=amctx_creatinine, aes(x=tx_s_years,y=value)) + geom_point() + stat_smooth()

#For d_cadaveric or d_type difference evolution
ggplot(data=amctx_creatinine, aes(x=tx_s_years,y=log(value))) + 
  geom_line(aes(group=amctx)) + stat_smooth() + ticks(0, 10, 0.5) + ylab("log(serum creatinine)") + 
  xlab("tx_s_years") + facet_grid(.~d_cadaveric)

#Tiny bit of difference for males and females in terms of where to put the knots
ggplot(data=amctx_creatinine[amctx_creatinine$tx_s_years<=1,], aes(x=tx_s_years,y=log(value))) + 
  geom_line(aes(group=amctx)) + stat_smooth() + ticks(0, 1, 0.05) + ylab("log(serum creatinine)") + 
  xlab("tx_s_years") + facet_grid(.~rec_gender)

#Tiny bit of difference for males and females in terms of where to put the knots
ggplot(data=amctx_creatinine, aes(x=tx_s_years,y=log(value))) + 
  geom_line(aes(group=amctx)) + stat_smooth() + ticks(0, 2, 0.5) + ylab("log(serum creatinine)") + 
  xlab("tx_s_years") + facet_grid(.~tx_dgf) + xlim(0,2)

idList = unique(amctx_creatinine$amctx)
ggplot(data=amctx_creatinine[amctx$amctx==idList[1],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[2],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[3],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[4],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[5],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[6],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[7],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[8],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[9],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[10],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[11],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[12],], aes(x=tx_s_days,y=value)) + geom_line()
ggplot(data=amctx_creatinine[amctx$amctx==idList[13],], aes(x=tx_s_days,y=value)) + geom_line()

# take residuals and repeat what you did above
linearmodel=lm(log(value)~rec_age_fwp1 + rec_gender +
                 tx_previoustx + d_age + d_gender + d_bmi +
                 rec_bmi + tx_hla + tx_pra + tx_dgf + 
                 tx_cit+ is_nr + is_aza + 
                 is_cni + is_mmf  + ah_nr + 
                 ah_diur + ah_ace + ah_arb + 
                 ah_raasi + ah_bb + ah_ccb  + 
                 dm_oad + dm_insulin + tx_dm + tx_hvdis + 
                 rr_sys + rr_dias + rr_map + 
                 tx_dial_days + d_type, data=amctx_creatinine)

amctx_creatinine$residuals = linearmodel$residuals

ggplot(data=amctx_creatinine, aes(x=tx_s_days,y=residuals)) + geom_point() + stat_smooth() + 
   ticks(200)

# This plot gives much better idea of where to put knots. females have some outliers which may influence knot selection
ggplot(data=amctx_creatinine, aes(x=tx_s_days,y=residuals)) + geom_point() + stat_smooth() + 
  facet_grid(rec_gender~.) + ticks(200)

# ah_ace + ah_arb = ah_raasi, perhaps it’s better to not include all three in the model. 
# Preferably no medication use is included, since they have not been collected in a dynamic way 
# (just baseline) and without doses. 
# tx_hla can be modelled as a continuous parameter
# (hessel: maybe not statistically correct, but this is often done also to reduce parameters and error)  
# 30 to 60 days for people to change trend
model_creatinine_feedback1 = lme(data=amctx_merged[!is.na(amctx_merged$creatinine),],
                                 fixed=log(creatinine) ~ rec_age_fwp1 + 
                                   rec_gender + tx_previoustx + d_age + d_gender + d_bmi + rec_bmi + 
                                   tx_hla + tx_pra + tx_dgf + ah_nr + tx_dm + tx_hvdis + 
                                   rr_sys + rr_dias +
                                   rr_map + tx_dial_days + d_cadaveric +
                                   ns(tx_s_years,knots=c(50, 100, 900)/365) * d_cadaveric + 
                                   ns(tx_s_years,knots=c(50, 100, 900)/365) * tx_dgf,
                                 random = ~ns(tx_s_years,knots=c(50, 100)/365)|amctx,
                                 control = lmeControl(opt = "optim"), method="ML")

#Despite the fact that BIC advises against interaction of time with tx_dgf, the graphs make it clear that such an interaction is possible
lme_creatinine_feedback2 = lme(data=amctx_merged[!is.na(amctx_merged$creatinine),],
                               fixed=log(creatinine) ~ rec_age_fwp1 + 
                                   rec_gender + d_age + tx_dgf + 
                                   tx_pra + ah_nr + tx_dm + d_cadaveric + 
                                   ns(tx_s_years,knots=c(50, 100, 900)/365) * d_cadaveric + 
                                   ns(tx_s_years,knots=c(50, 100, 900)/365) * tx_dgf,
                                 random = ~ns(tx_s_years,knots=c(50, 100)/365)|amctx,
                                 control = lmeControl(opt = "optim"), method="ML")

fixedSplines = lapply(1:50, function(i){
  second = runif(n=1, 50, 365)
  third = runif(n=1, second, 900)
  c(50, second, third)/365
})

model_splines = foreach(i=1:50, .packages = c("splines", "nlme", "ggplot2"), 
                        .export = c("amctx_merged")) %dopar%{
  fitUnivariteCreatinineModel(fixedSplineKnots = fixedSplines[[i]],
                               randomSplineKnots = fixedSplines[[i]][1:2],
                               boundaryKnots = c(0, 6))
                        }

lapply(model_splines, BIC)

model0 = fitUnivariteCreatinineModel(boundaryKnots = c(0,4))
model1 = fitUnivariteCreatinineModel(boundaryKnots = c(0,6), 
                                     fixedSplineKnots = c(50, 100, 1000)/365, 
                                     randomSplineKnots = c(50, 100)/365)

model2 = fitUnivariteCreatinineModel(boundaryKnots = c(0,6), 
                                    fixedSplineKnots = c(50, 80, 1000)/365, 
                                    randomSplineKnots = c(50, 80)/365)
model3 = fitUnivariteCreatinineModel(boundaryKnots = c(0,6), 
                                     fixedSplineKnots = c(30, 70, 1000)/365, 
                                     randomSplineKnots = c(30, 70)/365)

#After trying various models from the above list, the best fit is by model3
save.image("Rdata/feedbackmeeting.Rdata")

plotCreatinineFittedCurve(list(model1, model2, model3), 
                          individually = F, transform = T)

jmbayes_creatinine_orig = jointModelBayes(lmeObject = model3, survObject = coxModel, timeVar = "tx_s_years", control = list(n.iter=1000))
jmbayes_creatinine_replaced = replaceMCMCContents(mvJoint_creatinine_tdval, jmbayes_creatinine_orig)
